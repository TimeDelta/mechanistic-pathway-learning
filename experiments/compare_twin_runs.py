"""Paired comparison of every two slice configurations that differ in one argument (run_main_model_batch.CONFIGURATIONS),
scored three ways, so each recorded "this intervention moved or did not move the score" rests on the same procedure.

For each pair with five finished disease-cluster folds on both sides:
  - the paired per-fold difference of macro AUPRC with a t interval on four degrees of freedom. A percentile bootstrap over
    five fold values, which earlier documents used, gives an interval about a third narrower than this one;
  - the pooled out-of-fold difference of macro AUPRC and AUROC with a paired bootstrap over perturbations
    (paired_bootstrap_macro_difference), on the perturbations both runs scored;
  - the same pooled AUPRC difference after each score is ranked inside one degree stratum of one test fold
    (rank_normalise_within_groups), which gives no credit for ordering perturbations by degree.

A pair is set aside when the two runs' folds hold different perturbations, when either run was trained on a label
selection other than the one given here (--label-selection; none by default), or when either run's per-fold macro AUPRC
cannot be reproduced from its stored predictions and today's labels (the labels changed after it ran). Under a label
selection, the pairs it sets aside are left out of every reading (ranking_and_calibration_metrics, mask). A is the
configuration carrying the extra argument, or the alphabetically later one when both carry the same argument with
different values. Results are cached per pair under runs/twin_comparisons/, keyed on the fold results' modification
times, so a rerun only scores pairs whose runs changed. Idempotent; rewrites its outputs.

Usage:
  OMP_NUM_THREADS=1 python experiments/compare_twin_runs.py
  OMP_NUM_THREADS=1 python experiments/compare_twin_runs.py --pairs b3_linear_response_cofactors_spectral:b3_linear_response_cofactors
  OMP_NUM_THREADS=1 python experiments/compare_twin_runs.py --graph-dir data/processed/graph_full_neuronal \
      --evidence-dir data/processed/evidence_full --label-selection data/processed/label_selection/better_v1_full.parquet \
      --pairs full_linear_response_properties:full_linear_response --cache-dir runs/twin_comparisons_full \
      --markdown-output docs/twin_comparisons_full.md --json-output runs/twin_comparisons_full.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from experiments.run_main_model_batch import CONFIGURATIONS
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    degree_strata,
    micro_auprc,
    paired_bootstrap_macro_difference,
    paired_bootstrap_micro_difference,
    per_symptom_auprc,
    per_symptom_auroc,
    rank_normalise_within_groups,
    scorable_symptom,
)

RUN_ROOTS = (Path("runs/encoder"), Path("runs/full"), Path("runs"))
MINIMUM_POSITIVES_TO_SCORE = 5
NUM_FOLDS = 5
LABEL_REPRODUCTION_TOLERANCE = 1e-6
DEFAULT_CACHE_DIRECTORY = Path("runs/twin_comparisons")
RUN_DIRECTORY_SUFFIX = "_disease_cluster"  # main() sets it from --group-by and --lockbox, as run_main_model_batch.py names the directories


def run_directory(run: str) -> Path | None:
    for root in RUN_ROOTS:
        directory = root / f"{run}{RUN_DIRECTORY_SUFFIX}"
        if directory.exists():
            return directory
    return None


def fold_directories(run: str) -> list[Path]:
    directory = run_directory(run)
    if directory is None:
        return []
    return sorted(path for path in directory.glob("fold*_seed0") if (path / "DONE").exists() and (path / "results.json").exists())


def read_run(run: str) -> tuple[list[dict], pd.DataFrame | None]:
    """(per-fold results, out-of-fold predictions as perturbation x symptom rows with a fold column)."""
    results, frames = [], []
    for fold in fold_directories(run):
        result = json.loads((fold / "results.json").read_text())
        results.append(result)
        predictions = fold / "test_predictions.npy"
        if predictions.exists() and result.get("test_perturbation_ids") and result.get("symptoms"):
            scores = np.load(predictions)
            if scores.shape == (len(result["test_perturbation_ids"]), len(result["symptoms"])):
                frames.append(pd.DataFrame(scores, index=result["test_perturbation_ids"], columns=result["symptoms"])
                              .rename_axis("perturbation_id").reset_index().assign(fold=result.get("fold")))
    return results, (pd.concat(frames, ignore_index=True) if frames else None)


def pooled_matrix(predictions: pd.DataFrame | None, data) -> tuple[np.ndarray, np.ndarray]:
    """(perturbation x symptom matrix in data's row and column order, mask of rows the run scored). Perturbations the
    current labels do not contain are dropped."""
    matrix = np.full(data.outcomes.shape, np.nan)
    if predictions is None:
        return matrix, np.zeros(len(data.perturbation_ids), dtype=bool)
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    known = predictions[predictions.perturbation_id.isin(position_of)]
    matrix[[position_of[perturbation_id] for perturbation_id in known.perturbation_id]] = known[data.symptoms].to_numpy(dtype=float)
    return matrix, ~np.isnan(matrix).any(axis=1)


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Same rule as run_main_model.macro_auprc: symptoms with at least MINIMUM_POSITIVES_TO_SCORE labelled positives and one
    labelled negative."""
    values = [per_symptom_auprc(predictions, outcomes, symptom, mask) for symptom in range(outcomes.shape[1])
              if scorable_symptom(outcomes, symptom, MINIMUM_POSITIVES_TO_SCORE, mask)]
    return float(np.mean(values)) if values else float("nan")


def file_sha256(path: Path | None) -> str | None:
    if path is None:
        return None
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows_of(mask: np.ndarray | None, rows) -> np.ndarray | None:
    return None if mask is None else mask[rows]


def labels_reproduce(results: list[dict], predictions: pd.DataFrame | None, data) -> bool:
    """True when every fold's stored macro AUPRC is recovered from its stored predictions and today's labels."""
    if predictions is None:
        return False
    matrix, _ = pooled_matrix(predictions, data)
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    for result in results:
        if any(perturbation_id not in position_of for perturbation_id in result["test_perturbation_ids"]):
            return False
        rows = [position_of[perturbation_id] for perturbation_id in result["test_perturbation_ids"]]
        if abs(macro_auprc(matrix[rows], data.outcomes[rows], rows_of(data.label_mask, rows)) - result["macro_auprc"]) > LABEL_REPRODUCTION_TOLERANCE:
            return False
    return True


def fold_t_interval(differences: np.ndarray, confidence: float = 0.95) -> tuple[float, float, float]:
    """(mean, lower, upper) of paired per-fold differences from the t distribution with n - 1 degrees of freedom."""
    mean = float(differences.mean())
    half_width = float(stats.t.ppf(0.5 + confidence / 2, len(differences) - 1) * differences.std(ddof=1) / np.sqrt(len(differences)))
    return mean, mean - half_width, mean + half_width


def describe_interval(lower: float, upper: float) -> str:
    return "excludes zero" if lower > 0 or upper < 0 else "includes zero"


def compare(run: str, reference: str, data, strata: np.ndarray, num_bootstrap: int, require_all_folds: bool = True,
            label_selection_sha256: str | None = None) -> dict:
    """Every reading of run minus reference, or a reason the pair is set aside. label_selection_sha256 is the hash of the
    label selection data was loaded with (None for none); a run trained on another selection is set aside."""
    run_results, run_predictions = read_run(run)
    reference_results, reference_predictions = read_run(reference)
    entry = {"a": run, "b": reference, "folds_a": len(run_results), "folds_b": len(reference_results)}
    minimum_folds = NUM_FOLDS if require_all_folds else 2
    if len(run_results) < minimum_folds or len(reference_results) < minimum_folds:
        return {**entry, "set_aside": f"folds finished: {run} {len(run_results)}, {reference} {len(reference_results)}"}
    run_by_fold = {result["fold"]: result for result in run_results}
    reference_by_fold = {result["fold"]: result for result in reference_results}
    shared_folds = sorted(set(run_by_fold) & set(reference_by_fold))
    if any(sorted(run_by_fold[fold]["test_perturbation_ids"]) != sorted(reference_by_fold[fold]["test_perturbation_ids"]) for fold in shared_folds):
        return {**entry, "set_aside": "the two runs' folds hold different perturbations"}
    for name, results in ((run, run_results), (reference, reference_results)):
        if any(result.get("label_selection_sha256") != label_selection_sha256 for result in results):
            return {**entry, "set_aside": f"{name} was trained on a different label selection from the one given"}
    for name, results, predictions in ((run, run_results, run_predictions), (reference, reference_results, reference_predictions)):
        if not labels_reproduce(results, predictions, data):
            return {**entry, "set_aside": f"{name}'s per-fold scores do not reproduce from today's labels"}
    per_fold_differences = np.array([run_by_fold[fold]["macro_auprc"] - reference_by_fold[fold]["macro_auprc"] for fold in shared_folds])
    mean, lower, upper = fold_t_interval(per_fold_differences)
    entry.update(per_fold_a=float(np.mean([run_by_fold[fold]["macro_auprc"] for fold in shared_folds])),
                 per_fold_a_sd=float(np.std([run_by_fold[fold]["macro_auprc"] for fold in shared_folds], ddof=1)),
                 per_fold_b=float(np.mean([reference_by_fold[fold]["macro_auprc"] for fold in shared_folds])),
                 per_fold_b_sd=float(np.std([reference_by_fold[fold]["macro_auprc"] for fold in shared_folds], ddof=1)),
                 per_fold_differences=per_fold_differences.tolist(), per_fold_difference={"difference": mean, "lower": lower, "upper": upper})
    run_matrix, run_rows = pooled_matrix(run_predictions, data)
    reference_matrix, reference_rows = pooled_matrix(reference_predictions, data)
    shared_rows = run_rows & reference_rows
    outcomes = data.outcomes[shared_rows]
    shared_mask = rows_of(data.label_mask, shared_rows)
    entry["pooled_rows"] = int(shared_rows.sum())
    entry["pooled_macro_auprc"] = paired_bootstrap_macro_difference(run_matrix[shared_rows], reference_matrix[shared_rows], outcomes, per_symptom_auprc, num_bootstrap,
                                                                     mask=shared_mask)
    entry["pooled_macro_auroc"] = paired_bootstrap_macro_difference(run_matrix[shared_rows], reference_matrix[shared_rows], outcomes, per_symptom_auroc, num_bootstrap,
                                                                    mask=shared_mask)
    entry["within_degree_strata_macro_auprc"] = paired_bootstrap_macro_difference(
        rank_normalise_within_groups(run_matrix[shared_rows], strata[shared_rows]),
        rank_normalise_within_groups(reference_matrix[shared_rows], strata[shared_rows]), outcomes, per_symptom_auprc, num_bootstrap, mask=shared_mask)
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    per_fold_micro_differences = []
    for fold in shared_folds:
        rows = np.array([position_of[p] for p in run_by_fold[fold]["test_perturbation_ids"] if p in position_of])
        fold_mask = rows_of(data.label_mask, rows)
        per_fold_micro_differences.append(micro_auprc(run_matrix[rows], data.outcomes[rows], fold_mask) - micro_auprc(reference_matrix[rows], data.outcomes[rows], fold_mask))
    micro_mean, micro_lower, micro_upper = fold_t_interval(np.array(per_fold_micro_differences))
    entry["per_fold_micro_difference"] = {"difference": micro_mean, "lower": micro_lower, "upper": micro_upper}
    entry["pooled_micro_auprc"] = paired_bootstrap_micro_difference(run_matrix[shared_rows], reference_matrix[shared_rows], outcomes, num_bootstrap, mask=shared_mask)
    entry["within_degree_strata_micro_auprc"] = paired_bootstrap_micro_difference(
        rank_normalise_within_groups(run_matrix[shared_rows], strata[shared_rows]),
        rank_normalise_within_groups(reference_matrix[shared_rows], strata[shared_rows]), outcomes, num_bootstrap, mask=shared_mask)
    return entry


def argument_pairs(arguments: list[str]) -> set[tuple[str, str | None]]:
    pairs, position = set(), 0
    while position < len(arguments):
        takes_value = position + 1 < len(arguments) and not arguments[position + 1].startswith("--")
        pairs.add((arguments[position], arguments[position + 1] if takes_value else None))
        position += 2 if takes_value else 1
    return pairs


def twin_pairs(configurations: dict[str, list[str]]) -> list[tuple[str, str, str]]:
    """(A, B, the differing argument) for every two configurations whose arguments differ in exactly one flag."""
    twins = []
    names = sorted(configurations)
    for position, first in enumerate(names):
        for second in names[position + 1:]:
            first_arguments, second_arguments = argument_pairs(configurations[first]), argument_pairs(configurations[second])
            only_first, only_second = first_arguments - second_arguments, second_arguments - first_arguments
            if len({flag for flag, _ in only_first | only_second}) != 1:
                continue
            if only_first and not only_second:
                run, reference = first, second
            elif only_second and not only_first:
                run, reference = second, first
            else:  # same flag, different values: alphabetical order, later as A
                run, reference = second, first
            change = " ".join(f"{flag} {value}" if value is not None else flag for flag, value in sorted(argument_pairs(configurations[run]) - argument_pairs(configurations[reference])))
            twins.append((run, reference, change))
    return twins


def fold_stamp(run: str) -> list[float]:
    return [round((fold / "results.json").stat().st_mtime, 3) for fold in fold_directories(run)]


def cached_compare(run: str, reference: str, data, strata: np.ndarray, num_bootstrap: int, cache_directory: Path = DEFAULT_CACHE_DIRECTORY,
                   require_all_folds: bool = True, label_selection_sha256: str | None = None) -> dict:
    """compare(), reusing a stored result while neither run's folds have changed since it was written."""
    cache_directory.mkdir(parents=True, exist_ok=True)
    cache = cache_directory / f"{run}__{reference}.json"
    run_stamp, reference_stamp = fold_stamp(run), fold_stamp(reference)
    all_folds_present = len(run_stamp) == NUM_FOLDS and len(reference_stamp) == NUM_FOLDS  # then require_all_folds changes nothing
    stamp = {"a": run_stamp, "b": reference_stamp, "num_bootstrap": num_bootstrap, "require_all_folds": None if all_folds_present else require_all_folds}
    if label_selection_sha256 is not None:  # absent without a selection, so caches written before the option stay valid
        stamp["label_selection_sha256"] = label_selection_sha256
    if cache.exists():
        stored = json.loads(cache.read_text())
        if stored.get("stamp") == stamp and ("set_aside" in stored["entry"] or "pooled_micro_auprc" in stored["entry"]):  # entries written before micro AUPRC are recomputed
            return stored["entry"]
    entry = compare(run, reference, data, strata, num_bootstrap, require_all_folds, label_selection_sha256)
    cache.write_text(json.dumps({"stamp": stamp, "entry": entry}, indent=1))
    return entry


def format_interval(reading: dict) -> str:
    return f"{reading['difference']:+.3f} [{reading['lower']:+.3f}, {reading['upper']:+.3f}]"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"), help="the metabolic graph, whose degrees define the strata")
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--pairs", nargs="*", default=None, help="A:B pairs to score instead of every one-argument twin")
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="the label selection both runs were trained on (experiments/build_label_selection.py); its set-aside pairs are left out of every reading")
    parser.add_argument("--group-by", choices=["disease_cluster", "disease_cluster_and_targets"], default="disease_cluster",
                        help="the leakage grouping both runs were split by (disease_cluster_and_targets for every full-graph run)")
    parser.add_argument("--lockbox", type=Path, default=None, help="the lockbox removed before the runs' folds were drawn (their directories end in _development)")
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIRECTORY)
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/twin_comparisons.md"))
    parser.add_argument("--json-output", type=Path, default=Path("runs/twin_comparisons.json"))
    arguments = parser.parse_args()

    global RUN_DIRECTORY_SUFFIX
    RUN_DIRECTORY_SUFFIX = f"_{arguments.group_by}" + ("_development" if arguments.lockbox is not None else "")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    if arguments.lockbox is not None:
        data = restrict_to_perturbations(data, ~read_lockbox(arguments.lockbox, data))
    label_selection_sha256 = file_sha256(arguments.label_selection)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, NUM_FOLDS, 0)
    strata = np.array([fold_by_perturbation[p] for p in data.perturbation_ids]) * 1000 + degree_strata(data.perturbation_degrees)
    if arguments.pairs:
        pairs = [(pair.split(":")[0], pair.split(":")[1], "given on the command line") for pair in arguments.pairs]
    else:
        pairs = twin_pairs(CONFIGURATIONS)
    entries = []
    for run, reference, change in pairs:
        entry = {**cached_compare(run, reference, data, strata, arguments.num_bootstrap, arguments.cache_dir,
                                  label_selection_sha256=label_selection_sha256), "change": change}
        entries.append(entry)
        print(f"{run} vs {reference}: {entry.get('set_aside') or format_interval(entry['pooled_macro_auprc'])}", flush=True)
    arguments.json_output.write_text(json.dumps(entries, indent=1))

    scored = [entry for entry in entries if "set_aside" not in entry]
    lines = ["# One-change twin comparisons (generated by experiments/compare_twin_runs.py)", "",
             (f"Every two slice configurations in experiments/run_main_model_batch.py whose arguments differ in one flag, both with "
             f"{NUM_FOLDS} finished disease-cluster folds. A carries the change. Three readings of A minus B in macro AUPRC: paired per fold "
             f"with a t interval ({NUM_FOLDS - 1} degrees of freedom); pooled out-of-fold with a paired bootstrap over perturbations "
             f"({arguments.num_bootstrap} resamples); and pooled after ranking each score inside one degree stratum of one test fold "
             f"(metabolic-graph degrees, {len(np.unique(degree_strata(data.perturbation_degrees)))} strata), which gives no credit for ordering "
             "perturbations by degree. ± is the sample standard deviation over folds (n - 1); the aggregation tables of "
             "aggregate_main_model_runs.py divide by n, which reads about 11 percent smaller at five folds. "
             "Pooled macro AUROC is in the JSON output. Intervals are 95 percent and are not corrected for the "
             f"{len(scored)} comparisons; at that count about {0.05 * len(scored):.1f} intervals would exclude zero by chance alone."), "",
             "| A (carries the change) | B | change | A per fold | B per fold | per fold, t | pooled, bootstrap | within degree strata |",
             "|---|---|---|---|---|---|---|---|"]
    for entry in sorted(scored, key=lambda entry: entry["pooled_macro_auprc"]["difference"]):
        lines.append(f"| {entry['a']} | {entry['b']} | `{entry['change']}` | {entry['per_fold_a']:.3f} ± {entry['per_fold_a_sd']:.3f} | "
                     f"{entry['per_fold_b']:.3f} ± {entry['per_fold_b_sd']:.3f} | {format_interval(entry['per_fold_difference'])} | "
                     f"{format_interval(entry['pooled_macro_auprc'])} | {format_interval(entry['within_degree_strata_macro_auprc'])} |")
    lines += ["", "## Micro AUPRC", "",
              "The same three readings in micro AUPRC, which ranks every labelled perturbation-symptom pair in one list, so symptoms below "
              "the five-positive cutoff of the macro average count too and pairs of frequent symptoms weigh more.", "",
              "| A (carries the change) | B | per fold, t | pooled, bootstrap | within degree strata |", "|---|---|---|---|---|"]
    for entry in sorted(scored, key=lambda entry: entry["pooled_macro_auprc"]["difference"]):
        lines.append(f"| {entry['a']} | {entry['b']} | {format_interval(entry['per_fold_micro_difference'])} | {format_interval(entry['pooled_micro_auprc'])} | "
                     f"{format_interval(entry['within_degree_strata_micro_auprc'])} |")
    lines += ["", "## Set aside", ""]
    for entry in entries:
        if "set_aside" in entry:
            lines.append(f"- {entry['a']} against {entry['b']} (`{entry['change']}`): {entry['set_aside']}.")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.markdown_output} ({len(scored)} scored, {len(entries) - len(scored)} set aside)")


if __name__ == "__main__":
    main()
