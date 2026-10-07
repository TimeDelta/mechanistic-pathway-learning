"""Aggregate run directories of experiments/run_main_model.py into one comparison table (design sections 6 and 7).

A run directory (runs/<configuration>_<group_by>/) holds one split per subdirectory with results.json,
test_predictions.npy and the test perturbation ids. For grouped folds the out-of-fold predictions are
pooled and scored exactly as the baseline harness scores B0 and B1 (pooled per-symptom AUPRC and AUROC
with bootstrap intervals, per-fold mean and standard deviation of the macro metrics, MRR, hits-at-3,
expected calibration error). For hold-out splits (curated modules or reconstruction subsystems) the
held-out sets are pooled the same way. Module supports and link matrices of noisy-OR runs are
summarized (support sizes, number of symptoms per module above the link threshold).

With --label-selection (experiments/build_label_selection.py) the positive pairs the selection sets aside are left out
of every metric (ranking_and_calibration_metrics, mask), and a run trained on another selection is refused.

Usage:
  python experiments/aggregate_main_model_runs.py --run-dirs runs/b6_default_disease_cluster runs/b3_sigmoid_disease_cluster \
      --baseline-results runs/baselines_disease_cluster/results.json --markdown-output docs/phase3_main_model.md
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    degree_strata,
    expected_calibration_error,
    hits_at_k,
    macro_auprc_by_degree_bin,
    mean_reciprocal_rank,
    micro_auprc,
    paired_bootstrap_macro_difference,
    paired_bootstrap_micro_difference,
    per_symptom_auprc,
    per_symptom_auroc,
    rank_normalise_within_groups,
    scorable_symptom,
    scored_rows,
    stratified_auroc,
)

MINIMUM_POSITIVES_TO_SCORE = 5


def rows_of(mask: np.ndarray | None, rows) -> np.ndarray | None:
    return None if mask is None else mask[rows]


def file_sha256(path: Path | None) -> str | None:
    if path is None:
        return None
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def macro_scores(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None = None) -> tuple[float, float]:
    auprcs, aurocs = [], []
    for symptom_index in range(outcomes.shape[1]):
        if not scorable_symptom(outcomes, symptom_index, MINIMUM_POSITIVES_TO_SCORE, mask):
            continue
        auprcs.append(per_symptom_auprc(predictions, outcomes, symptom_index, mask))
        aurocs.append(per_symptom_auroc(predictions, outcomes, symptom_index, mask))
    return (float(np.mean(auprcs)) if auprcs else float("nan"), float(np.mean(aurocs)) if aurocs else float("nan"))


def aggregate_run_directory(run_directory: Path, data, num_bootstrap: int, label_selection_sha256: str | None = None) -> dict | None:
    split_directories = sorted(path for path in run_directory.iterdir() if path.is_dir() and (path / "DONE").exists() and (path / "results.json").exists())
    split_directories = [path for path in split_directories if json.loads((path / "results.json").read_text()).get("time_split") is None]
    if not split_directories:
        return None
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    pooled_predictions = np.zeros_like(data.outcomes)
    scored = np.zeros(len(data.perturbation_ids), dtype=bool)
    per_split = []
    support_sizes, expected_support_sizes, symptoms_per_module, epochs = [], [], [], []
    code_commits: set[str] = set()
    for split_directory in split_directories:
        results = json.loads((split_directory / "results.json").read_text())
        if results.get("label_selection_sha256") != label_selection_sha256:
            raise ValueError(f"{split_directory}: trained on a different label selection from the one given (--label-selection)")
        predictions = np.load(split_directory / "test_predictions.npy")
        rows = np.array([position_of[p] for p in results["test_perturbation_ids"] if p in position_of])
        if len(rows) != len(results["test_perturbation_ids"]):
            raise ValueError(f"{split_directory}: test perturbations not all present in the current evidence table")
        pooled_predictions[rows] = predictions
        scored[rows] = True
        split_mask = rows_of(data.label_mask, rows)
        macro_auprc, macro_auroc = macro_scores(predictions, data.outcomes[rows], split_mask)
        per_split.append({"split": results["split"], "num_test": len(rows), "macro_auprc": macro_auprc, "macro_auroc": macro_auroc,
                          "micro_auprc": micro_auprc(predictions, data.outcomes[rows], split_mask),
                          "mean_reciprocal_rank": mean_reciprocal_rank(predictions, data.outcomes[rows], split_mask),
                          "hits_at_3": hits_at_k(predictions, data.outcomes[rows], 3, split_mask),
                          "epochs_completed": results.get("epochs_completed"), "best_epoch": results.get("best_epoch")})
        epochs.append(results.get("epochs_completed", 0))
        code_commits.update(entry["commit"][:7] for entry in results.get("code_provenance", []) if entry.get("commit"))
        if "module_symptom_links" in results:
            links = np.array(results["module_symptom_links"])
            support_sizes.append(results["module_support_sizes"])
            expected_support_sizes.append(results.get("module_expected_support_sizes", []))
            symptoms_per_module.append((links > 0.5).sum(axis=1).tolist())
    outcomes, predictions, mask = data.outcomes[scored], pooled_predictions[scored], rows_of(data.label_mask, scored)
    np.save(run_directory / "pooled_predictions.npy", pooled_predictions)
    np.save(run_directory / "pooled_scored_rows.npy", scored)
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        if not scorable_symptom(outcomes, symptom_index, MINIMUM_POSITIVES_TO_SCORE, mask):
            continue
        labelled = scored_rows(len(outcomes), symptom_index, mask)
        if mask is None:
            auprc = bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap)
            auroc = bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap)
        else:
            auprc = bootstrap_interval(lambda p, y, m: per_symptom_auprc(p, y, symptom_index, m), predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask)
            auroc = bootstrap_interval(lambda p, y, m: per_symptom_auroc(p, y, symptom_index, m), predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask)
        per_symptom[symptom] = {"positives": int(outcomes[labelled, symptom_index].sum()), "base_rate": float(outcomes[labelled, symptom_index].mean()),
                                "auprc": auprc.__dict__, "auroc": auroc.__dict__}
    fold_auprcs = [entry["macro_auprc"] for entry in per_split if not np.isnan(entry["macro_auprc"])]
    fold_aurocs = [entry["macro_auroc"] for entry in per_split if not np.isnan(entry["macro_auroc"])]
    fold_micro_auprcs = [entry["micro_auprc"] for entry in per_split if not np.isnan(entry["micro_auprc"])]
    first_results = json.loads((split_directories[0] / "results.json").read_text())
    return {
        "run_directory": str(run_directory), "num_splits": len(split_directories), "num_scored_perturbations": int(scored.sum()),
        "code_commits": sorted(code_commits),
        "arguments": {key: first_results["arguments"].get(key) for key in ("head", "field", "pooling", "num_modules", "group_by", "description_length_coefficient", "learning_rate", "node_state_dim", "num_layers", "permute_labels")},
        "per_symptom": per_symptom,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "per_fold_macro_auprc_mean": float(np.mean(fold_auprcs)) if fold_auprcs else float("nan"), "per_fold_macro_auprc_sd": float(np.std(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auroc_mean": float(np.mean(fold_aurocs)) if fold_aurocs else float("nan"), "per_fold_macro_auroc_sd": float(np.std(fold_aurocs)) if fold_aurocs else float("nan"),
        "micro_auprc": bootstrap_interval(micro_auprc, predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask).__dict__,
        "per_fold_micro_auprc_mean": float(np.mean(fold_micro_auprcs)) if fold_micro_auprcs else float("nan"), "per_fold_micro_auprc_sd": float(np.std(fold_micro_auprcs)) if fold_micro_auprcs else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, outcomes, mask), "hits_at_3": hits_at_k(predictions, outcomes, 3, mask),
        "expected_calibration_error": expected_calibration_error(predictions, outcomes, mask=mask),
        "macro_auprc_by_degree_bin": macro_auprc_by_degree_bin(predictions, outcomes, data.perturbation_degrees[scored], mask=mask),
        "label_selection_sha256": label_selection_sha256,
        "per_split": per_split, "mean_epochs": float(np.mean(epochs)) if epochs else float("nan"),
        "module_support_sizes": support_sizes, "module_expected_support_sizes": expected_support_sizes, "symptoms_per_module_above_half": symptoms_per_module,
    }


def collect_time_split_runs(run_directories: list[Path]) -> list[dict]:
    rows = []
    for run_directory in run_directories:
        for split_directory in sorted(path for path in run_directory.iterdir() if path.is_dir() and (path / "DONE").exists()):
            results = json.loads((split_directory / "results.json").read_text())
            if results.get("time_split"):
                rows.append({"run": run_directory.name, "split": results["split"], **results["time_split"]})
    return rows


def load_pooled_predictions(aggregated: dict, run_directories: list[Path], baseline_directory: Path | None, baseline_split: str, data) -> tuple[dict, dict]:
    """Pooled out-of-split predictions and scored-row masks of every aggregated run and every baseline."""
    predictions_by_name: dict[str, np.ndarray] = {}
    rows_by_name: dict[str, np.ndarray] = {}
    for run_directory in run_directories:
        if run_directory.name in aggregated and (run_directory / "pooled_predictions.npy").exists():
            predictions_by_name[run_directory.name] = np.load(run_directory / "pooled_predictions.npy")
            rows_by_name[run_directory.name] = np.load(run_directory / "pooled_scored_rows.npy")
    if baseline_directory is not None and (baseline_directory / "perturbation_ids.json").exists():
        baseline_ids = json.loads((baseline_directory / "perturbation_ids.json").read_text())
        if baseline_ids == data.perturbation_ids:
            rows_path = baseline_directory / f"scored_rows_{baseline_split}.npy"
            for prediction_path in sorted(baseline_directory.glob(f"predictions_{baseline_split}_*.npy")):
                name = prediction_path.stem[len(f"predictions_{baseline_split}_"):]
                predictions_by_name[name] = np.load(prediction_path)
                rows_by_name[name] = np.load(rows_path) if rows_path.exists() else np.ones(len(data.perturbation_ids), dtype=bool)
        else:
            print(f"warning: {baseline_directory} holds other perturbations ({len(baseline_ids)} against {len(data.perturbation_ids)}; a lockbox removed on one side only?); its baselines are left out of the stratified readings")
    return predictions_by_name, rows_by_name


def within_degree_strata(aggregated: dict, predictions_by_name: dict, rows_by_name: dict, data, num_bootstrap: int, num_folds: int = 5, seed: int = 0) -> tuple[dict, list[dict]]:
    """Scores that give no credit for ordering perturbations by degree: macro AUPRC after ranking each score inside its
    group (rank_normalise_within_groups) and macro AUROC from pairs inside a group only (stratified_auroc), plus the
    paired bootstrap of the first between every run and every baseline. A group is a degree stratum inside one test
    fold: pooled out-of-fold scores carry a fold effect (each fold's model learns its own training base rate, which runs
    opposite to its test fold's), and strata that cut across folds would score it (popularity, constant within a
    symptom and fold, reached a stratified AUROC of 0.419 with strata across folds). Runs and baselines share the
    fold assignment (assign_grouped_folds, same seed), and strata come from all perturbations."""
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, num_folds, seed)
    fold_of_row = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
    strata = fold_of_row * 1000 + degree_strata(data.perturbation_degrees)
    scores, normalised_by_name = {}, {}
    for name, predictions in predictions_by_name.items():
        rows = rows_by_name[name]
        normalised = rank_normalise_within_groups(predictions[rows], strata[rows])
        outcomes, mask = data.outcomes[rows], rows_of(data.label_mask, rows)
        scorable = [s for s in range(outcomes.shape[1]) if scorable_symptom(outcomes, s, MINIMUM_POSITIVES_TO_SCORE, mask)]
        normalised_by_name[name] = (normalised, rows)
        scores[name] = {"macro_auprc_within_degree_strata": float(np.mean([per_symptom_auprc(normalised, outcomes, s, mask) for s in scorable])),
                        "macro_auroc_stratified_by_degree": float(np.nanmean([stratified_auroc(predictions[rows], outcomes, strata[rows], s, mask) for s in scorable])),
                        "micro_auprc_within_degree_strata": micro_auprc(normalised, outcomes, mask),
                        "num_strata": int(len(np.unique(degree_strata(data.perturbation_degrees))))}
    comparisons = []
    for run_name in aggregated:
        if run_name not in normalised_by_name:
            continue
        for other in normalised_by_name:
            if other in aggregated:
                continue
            run_normalised, run_rows = normalised_by_name[run_name]
            other_normalised, other_rows = normalised_by_name[other]
            shared = run_rows & other_rows
            run_positions = np.flatnonzero(run_rows)
            other_positions = np.flatnonzero(other_rows)
            run_pick = np.isin(run_positions, np.flatnonzero(shared))
            other_pick = np.isin(other_positions, np.flatnonzero(shared))
            difference = paired_bootstrap_macro_difference(run_normalised[run_pick], other_normalised[other_pick], data.outcomes[shared], per_symptom_auprc, num_bootstrap,
                                                           mask=rows_of(data.label_mask, shared))
            micro_difference = paired_bootstrap_micro_difference(run_normalised[run_pick], other_normalised[other_pick], data.outcomes[shared], num_bootstrap,
                                                                 mask=rows_of(data.label_mask, shared))
            comparisons.append({"a": run_name, "b": other, "rows": int(shared.sum()), "macro_auprc_within_degree_strata": difference,
                                "micro_auprc_within_degree_strata": micro_difference})
    return scores, comparisons


def paired_comparisons(aggregated: dict, run_directories: list[Path], baseline_directory: Path | None, baseline_split: str, data, num_bootstrap: int, include_baseline_pairs: bool = True) -> list[dict]:
    """Paired bootstrap of the pooled macro AUPRC and AUROC difference between every run and every baseline (and between runs) on the rows both scored."""
    predictions_by_name, rows_by_name = load_pooled_predictions(aggregated, run_directories, baseline_directory, baseline_split, data)
    comparisons = []
    names = list(predictions_by_name)
    for position, first in enumerate(names):
        for second in names[position + 1:]:
            if first not in aggregated and second not in aggregated and not include_baseline_pairs:
                continue
            if first not in aggregated and second in aggregated:
                first, second = second, first  # a run-versus-baseline pair is reported once, with the run as A
            rows = rows_by_name[first] & rows_by_name[second]
            auprc = paired_bootstrap_macro_difference(predictions_by_name[first][rows], predictions_by_name[second][rows], data.outcomes[rows], per_symptom_auprc, num_bootstrap,
                                                      mask=rows_of(data.label_mask, rows))
            auroc = paired_bootstrap_macro_difference(predictions_by_name[first][rows], predictions_by_name[second][rows], data.outcomes[rows], per_symptom_auroc, num_bootstrap,
                                                      mask=rows_of(data.label_mask, rows))
            micro = paired_bootstrap_micro_difference(predictions_by_name[first][rows], predictions_by_name[second][rows], data.outcomes[rows], num_bootstrap,
                                                      mask=rows_of(data.label_mask, rows))
            comparisons.append({"a": first, "b": second, "rows": int(rows.sum()), "macro_auprc": auprc, "macro_auroc": auroc, "micro_auprc": micro})
    return comparisons


def summary_row(name: str, entry: dict) -> str:
    micro = (entry.get("micro_auprc") or {}).get("point", float("nan"))  # baseline results written before micro AUPRC have none
    return (f"| {name} | {entry['macro_auprc']:.3f} | {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f} | {micro:.3f} | "
            f"{entry.get('per_fold_micro_auprc_mean', float('nan')):.3f} ± {entry.get('per_fold_micro_auprc_sd', float('nan')):.3f} | {entry['macro_auroc']:.3f} | "
            f"{entry['per_fold_macro_auroc_mean']:.3f} ± {entry['per_fold_macro_auroc_sd']:.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} | "
            f"{entry.get('expected_calibration_error', float('nan')):.3f} | {entry['num_scored_perturbations']} |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--group-by", choices=["gene", "disease_cluster", "disease_cluster_and_targets"], default="disease_cluster")
    parser.add_argument("--baseline-results", type=Path, default=None, help="results.json of run_baselines.py for the same split, to put B0 and B1 in the same table")
    parser.add_argument("--baseline-split", default="grouped", help="key under 'splits' in the baseline results to show")
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="the label selection the runs were trained on (experiments/build_label_selection.py); its set-aside pairs are left out of every metric")
    parser.add_argument("--lockbox", type=Path, default=None,
                        help="the lockbox removed before the runs' folds were drawn; pass the baseline results of run_baselines.py --lockbox with it")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--title", default="Phase 3: proposed model and sigmoid-head baseline")
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase3_main_model.md"))
    parser.add_argument("--json-output", type=Path, default=Path("runs/phase3_aggregate.json"))
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    if arguments.lockbox is not None:
        data = restrict_to_perturbations(data, ~read_lockbox(arguments.lockbox, data))
    label_selection_sha256 = file_sha256(arguments.label_selection)
    aggregated = {}
    for run_directory in arguments.run_dirs:
        entry = aggregate_run_directory(run_directory, data, arguments.num_bootstrap, label_selection_sha256)
        if entry is None:
            print(f"{run_directory}: no finished splits")
            continue
        aggregated[run_directory.name] = entry
        print(f"{run_directory.name:40s} splits {entry['num_splits']}  macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  "
              f"micro AUPRC {entry['micro_auprc']['point']:.3f} (per fold {entry['per_fold_micro_auprc_mean']:.3f})  "
              f"macro AUROC {entry['macro_auroc']:.3f}  MRR {entry['mean_reciprocal_rank']:.3f}  hits@3 {entry['hits_at_3']:.3f}  ECE {entry['expected_calibration_error']:.3f}")
    baseline_entries = {}
    baseline_directory = None
    if arguments.baseline_results and arguments.baseline_results.exists():
        baseline_results = json.loads(arguments.baseline_results.read_text())
        if (baseline_results.get("label_selection") or {}).get("sha256") != label_selection_sha256:
            raise ValueError(f"{arguments.baseline_results}: baselines fitted on a different label selection from the one given (--label-selection)")
        baseline_entries = baseline_results["splits"].get(arguments.baseline_split, {})
        baseline_directory = arguments.baseline_results.parent
    comparisons = paired_comparisons(aggregated, arguments.run_dirs, baseline_directory, arguments.baseline_split, data, arguments.num_bootstrap)
    strata_scores, strata_comparisons = within_degree_strata(aggregated, *load_pooled_predictions(aggregated, arguments.run_dirs, baseline_directory, arguments.baseline_split, data),
                                                             data, arguments.num_bootstrap)
    for comparison in comparisons:
        print(f"{comparison['a']:34s} vs {comparison['b']:28s} macro AUPRC diff {comparison['macro_auprc']['difference']:+.3f} [{comparison['macro_auprc']['lower']:+.3f}, {comparison['macro_auprc']['upper']:+.3f}]  "
              f"macro AUROC diff {comparison['macro_auroc']['difference']:+.3f} [{comparison['macro_auroc']['lower']:+.3f}, {comparison['macro_auroc']['upper']:+.3f}]")
    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps({"runs": aggregated, "paired_comparisons": comparisons, "within_degree_strata": strata_scores,
                                                 "within_degree_strata_comparisons": strata_comparisons}, indent=1))

    lines = [f"# {arguments.title} (generated by experiments/aggregate_main_model_runs.py)", "",
             f"{len(data.perturbation_ids)} perturbations, {len(data.symptoms)} symptoms, leakage groups by {arguments.group_by}. Out-of-split predictions pooled; per-fold values are mean ± standard deviation over splits. "
             "Baseline rows come from experiments/run_baselines.py on the same split. ECE: expected calibration error over all perturbation-symptom pairs. "
             "Micro AUPRC ranks every labelled perturbation-symptom pair in one list, so symptoms below the five-positive cutoff of the macro average count too (n/a for baseline results written before it existed).", "",
             "| model | macro AUPRC (pooled) | macro AUPRC (per fold) | micro AUPRC (pooled) | micro AUPRC (per fold) | macro AUROC (pooled) | macro AUROC (per fold) | MRR | hits@3 | ECE | scored perturbations |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, entry in baseline_entries.items():
        lines.append(summary_row(name, {**entry, "expected_calibration_error": float("nan")}))
    for name, entry in aggregated.items():
        lines.append(summary_row(name, entry))
    if aggregated:
        degree_bins = list(next(iter(aggregated.values()))["macro_auprc_by_degree_bin"])
        lines += ["", "Macro AUPRC by perturbation degree tercile (pooled predictions):", "", "| model | " + " | ".join(degree_bins) + " |", "|---|" + "---|" * len(degree_bins)]
        for name, entry in list(baseline_entries.items()) + list(aggregated.items()):
            bins = entry.get("macro_auprc_by_degree_bin", {})
            lines.append(f"| {name} | " + " | ".join(f"{bins[b]:.3f}" if b in bins else "n/a" for b in degree_bins) + " |")
    if comparisons:
        lines += ["", "## Paired bootstrap comparisons (design section 7)", "", "Difference in pooled macro AUPRC and macro AUROC between two models on the rows both scored; 95 percent percentile interval of the paired bootstrap over perturbations. The pre-registered primary endpoint asks for a difference of at least 0.05 with an interval excluding zero.", "",
                  "| A | B | rows | macro AUPRC A - B [95% CI] | resamples favoring A | micro AUPRC A - B [95% CI] | macro AUROC A - B [95% CI] |", "|---|---|---|---|---|---|---|"]
        for comparison in comparisons:
            a, m, r = comparison["macro_auprc"], comparison["micro_auprc"], comparison["macro_auroc"]
            lines.append(f"| {comparison['a']} | {comparison['b']} | {comparison['rows']} | {a['difference']:+.3f} [{a['lower']:+.3f}, {a['upper']:+.3f}] | {a['fraction_resamples_favoring_a']:.2f} | "
                         f"{m['difference']:+.3f} [{m['lower']:+.3f}, {m['upper']:+.3f}] | {r['difference']:+.3f} [{r['lower']:+.3f}, {r['upper']:+.3f}] |")
    if strata_scores:
        lines += ["", "## Within degree strata", "",
                  f"Rankings scored inside groups of one degree stratum and one test fold ({next(iter(strata_scores.values()))['num_strata']} degree strata from quantile edges, equal degrees always together), so ordering perturbations by degree earns nothing and the fold effect of pooled out-of-fold scores cannot enter: "
                  "macro AUPRC after each score is replaced by its rank inside its group, and macro AUROC from positive-negative pairs inside a group only. Degree-scaled popularity can still order perturbations of different degree inside a wide stratum.", "",
                  "| model | macro AUPRC within degree strata | micro AUPRC within degree strata | macro AUROC stratified by degree |", "|---|---|---|---|"]
        for name, score in strata_scores.items():
            lines.append(f"| {name} | {score['macro_auprc_within_degree_strata']:.3f} | {score['micro_auprc_within_degree_strata']:.3f} | {score['macro_auroc_stratified_by_degree']:.3f} |")
        lines += ["", "Paired bootstrap of the within-strata macro and micro AUPRC, each run against each baseline (95 percent interval):", "",
                  "| A | B | rows | macro difference [95% CI] | resamples favoring A | micro difference [95% CI] |", "|---|---|---|---|---|---|"]
        for comparison in strata_comparisons:
            d, m = comparison["macro_auprc_within_degree_strata"], comparison["micro_auprc_within_degree_strata"]
            lines.append(f"| {comparison['a']} | {comparison['b']} | {comparison['rows']} | {d['difference']:+.3f} [{d['lower']:+.3f}, {d['upper']:+.3f}] | {d['fraction_resamples_favoring_a']:.2f} | "
                         f"{m['difference']:+.3f} [{m['lower']:+.3f}, {m['upper']:+.3f}] |")
    lines += ["", "## Configurations", "", "Commits are those the splits ran under (recorded from 6 October 2026; 'not recorded' marks earlier runs, which evaluated noisy-OR gates with the Louizos test-time estimator rather than the expected training gate).", "",
              "| run | head | field | pooling | modules | description-length coefficient | learning rate | state dim | layers | labels permuted | splits | mean epochs | commits |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, entry in aggregated.items():
        a = entry["arguments"]
        commits = ", ".join(entry.get("code_commits") or []) or "not recorded"
        lines.append(f"| {name} | {a.get('head')} | {a.get('field')} | {a.get('pooling')} | {a.get('num_modules')} | {a.get('description_length_coefficient')} | {a.get('learning_rate')} | {a.get('node_state_dim')} | {a.get('num_layers')} | {a.get('permute_labels')} | {entry['num_splits']} | {entry['mean_epochs']:.1f} | {commits} |")
    lines += ["", "## Per-symptom AUPRC (pooled; 95 percent bootstrap interval over perturbations; base rate in parentheses)", "",
              "| symptom | " + " | ".join(list(baseline_entries) + list(aggregated)) + " |", "|---|" + "---|" * (len(baseline_entries) + len(aggregated))]
    for symptom in data.symptoms:
        cells = []
        for entry in list(baseline_entries.values()) + list(aggregated.values()):
            s = entry["per_symptom"].get(symptom)
            cells.append("n/a" if s is None else f"{s['auprc']['point']:.3f} [{s['auprc']['lower']:.3f}, {s['auprc']['upper']:.3f}] ({s['base_rate']:.3f})")
        lines.append(f"| {symptom} | " + " | ".join(cells) + " |")
    for name, entry in aggregated.items():
        if entry["module_support_sizes"]:
            lines += ["", f"## Module structure, {name}", "", "Support size counts gates above 0.5 in the evaluation gate (see the commits note under Configurations); expected support is the hard-concrete expected number of nonzero gates. Symptoms per module counts links above 0.5.", "",
                      "| split | support sizes | expected support sizes | symptoms per module above 0.5 |", "|---|---|---|---|"]
            for split_entry, sizes, expected, per_module in zip(entry["per_split"], entry["module_support_sizes"], entry["module_expected_support_sizes"], entry["symptoms_per_module_above_half"]):
                lines.append(f"| {split_entry['split']} | {sizes} | {[round(x) for x in expected] if expected else 'n/a'} | {per_module} |")
    time_split_rows = collect_time_split_runs(arguments.run_dirs)
    if time_split_rows:
        lines += ["", "## Time split (monogenic pairs by OMIM biocuration date)", "", "Training positives are the pairs dated on or before the cutoff; scored pairs are those that could still become positive; permuted: new-positive labels permuted within each symptom's scored pairs.", "",
                  "| run | split | training positives | new positives | scored pairs | macro AUPRC | macro AUPRC (permuted) | macro AUROC | macro AUROC (permuted) |", "|---|---|---|---|---|---|---|---|---|"]
        for row in time_split_rows:
            lines.append(f"| {row['run']} | {row['split']} | {row['training_positive_pairs']} | {row['new_positive_pairs']} | {row['scored_pairs']} | {row['macro_auprc']:.3f} | {row['macro_auprc_permuted']:.3f} | {row['macro_auroc']:.3f} | {row['macro_auroc_permuted']:.3f} |")
    lines += ["", "## Per-split macro metrics", "", "| run | split | test size | macro AUPRC | macro AUROC | MRR | hits@3 | epochs | best epoch |", "|---|---|---|---|---|---|---|---|---|"]
    for name, entry in aggregated.items():
        for split_entry in entry["per_split"]:
            lines.append(f"| {name} | {split_entry['split']} | {split_entry['num_test']} | {split_entry['macro_auprc']:.3f} | {split_entry['macro_auroc']:.3f} | {split_entry['mean_reciprocal_rank']:.3f} | {split_entry['hits_at_3']:.3f} | {split_entry['epochs_completed']} | {split_entry['best_epoch']} |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
