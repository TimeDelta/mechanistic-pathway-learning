"""Score the lockbox once: the confirmatory test of docs/preregistration.md (specification of 8 October 2026).

Inputs, all trained on the development set and scored on configs/lockbox_v1.json:
- the four confirmatory models (run_main_model_batch.py --lockbox ... --score-lockbox), each with five seeds and, per
  seed, a run on the real labels (lockbox_seed<k>), one on labels permuted within degree strata (lockbox_seed<k>_permuted)
  and one on a degree-preserving rewiring of the graph (lockbox_seed<k>_rewired); a run counts only if it was refitted
  after early stopping with the large leakage groups in training (amendments of 8 October 2026), its predictions are
  finite, one row per lockbox perturbation and one column per symptom, its recorded arguments are those of its
  configuration, its graph, evidence, descriptor and cell-class files are the ones scored, and a permuted run used the
  permutation this script draws for its seed;
- the baselines of run_baselines.py --lockbox ... --score-lockbox --seed <k> for the same seeds (the seed sets the label
  permutation, which the trainer and the baselines draw identically, and the TransE initialisation).

Readings, each the seed mean of a model's score minus the score of the best baseline on the same lockbox rows:
- macro AUPRC over the symptoms with five or more kept positives in the lockbox (fixed in the lockbox file);
- micro AUPRC over every scored symptom (five or more kept positives in the whole data);
- both again after ranking every score inside degree strata (degree_strata over all perturbations), which gives no
  credit for ordering perturbations by degree;
- the permuted-label difference in differences: (model minus baseline on the real labels) minus (model minus the same
  baseline, both trained and scored on the permuted labels);
- the rewiring difference: the model on the real graph minus the same model on its rewired graph (real labels).
The best baseline is chosen per reading by its seed-mean score on the lockbox, without reference to any model.

Inference: a paired bootstrap over the lockbox's leakage groups (--bootstrap-unit group, the user's decision of 8 October
2026: whole groups are drawn with replacement, the same resampled rows for every model, seed, baseline and labelling),
one-sided p = (1 + #{resampled difference <= 0}) / (B + 1). The bootstrap over single perturbations, which understates
the variance when perturbations of one group co-vary, is reported beside it as a sensitivity reading. Two hypotheses per model:
- H1 (prediction): the macro, micro and both within-strata readings exceed zero, the macro difference is at least
  --minimum-macro-difference and the micro difference at least --minimum-micro-difference. Intersection-union test:
  p(H1) is the largest of the four p-values.
- H2 (graph content): H1 and both rewiring readings exceed zero; p(H2) = max(p(H1), the two rewiring p-values).
The two permutation readings are secondary (the user's decision of 8 October 2026): reported with their intervals and
p-values, outside both hypotheses.
Each model is tested on its own at one-sided --alpha (the user's decision of 8 October 2026: every model is a
candidate with the full level): H1 first and, only if H1 is confirmed, H2 at the same level. With a missing or refused
run the scorer stops before writing anything, so the run can be resumed; with --allow-incomplete it scores and that
model is not confirmed. No correction is made across the models, so the chance that at least one of four models is
confirmed by luck is above --alpha (at most 1 - (1 - alpha)^4, about 0.096 at 0.025, if the four were independent; less,
since they share the data and the baselines); the output states this beside the results.

The lockbox is scored once. A SCORED marker records when; a second scoring needs --rescore and is listed in the output.

Usage:
  python experiments/score_confirmatory.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox
from mechanistic_pathway_learning.evaluation.negative_controls import degree_stratified_row_permutation
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import degree_strata, rank_normalise_within_groups
from run_main_model import RESUME_CONTROL_ARGUMENTS, array_sha256, build_argument_parser
from run_main_model_batch import CONFIGURATIONS

CONFIRMATORY_MODELS = ("confirmatory_message_passing_noisy_or", "confirmatory_message_passing_sigmoid",
                       "confirmatory_linear_response_noisy_or", "confirmatory_linear_response_sigmoid")
BASELINES = ("popularity", "degree_popularity", "random_walk_with_restart", "knowledge_graph_embedding_transe")
VARIANT_SUFFIXES = {"real": "", "permuted": "_permuted", "rewired": "_rewired"}
MINIMUM_POSITIVES_TO_SCORE = 5
H1_READINGS = ("macro", "micro", "macro_within_degree_strata", "micro_within_degree_strata")
H2_EXTRA_READINGS = ("macro_rewiring", "micro_rewiring")
SECONDARY_READINGS = ("macro_permutation", "micro_permutation")
ALL_READINGS = H1_READINGS + H2_EXTRA_READINGS + SECONDARY_READINGS


def average_precision(scores: np.ndarray, labels: np.ndarray) -> float:
    """sklearn.metrics.average_precision_score for one ranking (thresholds at the distinct scores), without its overhead."""
    if labels.size == 0:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    sorted_scores, sorted_labels = scores[order], labels[order]
    threshold_positions = np.r_[np.flatnonzero(np.diff(sorted_scores)), sorted_labels.size - 1]
    true_positives = np.cumsum(sorted_labels)[threshold_positions]
    if true_positives[-1] == 0 or true_positives[-1] == sorted_labels.size:
        return float("nan")
    precision = true_positives / (threshold_positions + 1)
    recall = true_positives / true_positives[-1]
    return float(np.sum(np.diff(np.r_[0.0, recall]) * precision))


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray, columns: list[int]) -> float:
    """Mean over columns of the AUPRC over labelled rows; a column with no positive or no negative in the rows is skipped."""
    values = []
    for column in columns:
        labelled = mask[:, column]
        value = average_precision(predictions[labelled, column], outcomes[labelled, column])
        if not np.isnan(value):
            values.append(value)
    return float(np.mean(values)) if values else float("nan")


def micro_auprc(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray, columns: list[int]) -> float:
    labelled = mask[:, columns]
    return average_precision(predictions[:, columns][labelled], outcomes[:, columns][labelled])


def fixed_sequence_decisions(complete: bool, p_h1: float, p_h2: float, meets_minimum_differences: bool, alpha: float) -> dict[str, bool]:
    """H1 then H2 at the same level: H2 is tested only once H1 is confirmed, so the chance of any false confirmation stays at
    alpha (the fixed-sequence test, one of the procedures Bretz et al., Statistics in Medicine 28, 586, 2009, write as a graph)."""
    h1_confirmed = bool(complete and meets_minimum_differences and p_h1 <= alpha)
    return {"H1_confirmed": h1_confirmed, "H2_confirmed": bool(h1_confirmed and p_h2 <= alpha)}


def bootstrap_rows(unit: str, members_of_group: list[np.ndarray], num_rows: int, generator: np.random.Generator) -> np.ndarray:
    """Row indices of one bootstrap resample of the lockbox. unit "perturbation" draws num_rows rows with replacement;
    unit "group" draws as many leakage groups as the lockbox holds, with replacement, and takes every row of each drawn
    group (a cluster bootstrap), so the dependence among perturbations sharing disease annotations or targets stays
    inside a resample (the user's decision of 8 October 2026)."""
    if unit == "perturbation":
        return generator.integers(0, num_rows, num_rows)
    if unit == "group":
        drawn = generator.integers(0, len(members_of_group), len(members_of_group))
        return np.concatenate([members_of_group[group] for group in drawn])
    raise ValueError(f"unknown bootstrap unit {unit!r}")


def one_sided_p(resampled_differences: np.ndarray) -> float:
    finite = resampled_differences[np.isfinite(resampled_differences)]
    return float((1 + np.sum(finite <= 0)) / (finite.size + 1)) if finite.size else 1.0


def file_sha256(path) -> str | None:
    """SHA-256 of a file; None for no path (a model without descriptors or cell-class weights)."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest() if path else None


# set by the variant (checked separately) or checked by content (the lockbox's hash)
VARIANT_AND_CONTENT_ARGUMENTS = {"permute_labels", "rewire_swaps_per_edge", "keep_reciprocated_relations_symmetric", "lockbox"}


def recorded_arguments(namespace: argparse.Namespace) -> dict:
    """Arguments in the form run_main_model.py writes them into results.json."""
    return {key: (value if isinstance(value, (int, float, str, bool, list, type(None))) else str(value)) for key, value in vars(namespace).items()}


def expected_run_arguments(model: str, seed: int, group_by: str, lockbox: Path) -> dict:
    """The trainer arguments of a lockbox run of model as runs/full/confirmatory_lockbox.sh launches it
    (run_main_model_batch.py --configuration model --score-lockbox), without the resume controls and the variant."""
    parsed = build_argument_parser().parse_args(["--group-by", group_by, *CONFIGURATIONS[model], "--lockbox", str(lockbox), "--score-lockbox", "--seed", str(seed)])
    return {key: value for key, value in recorded_arguments(parsed).items() if key not in RESUME_CONTROL_ARGUMENTS | VARIANT_AND_CONTENT_ARGUMENTS}


def read_model_run(split_directory: Path, lockbox_ids: list[str], lockbox_sha256: str, selection_sha256: str, variant: str,
                   require_refit: bool = True, symptoms: list[str] | None = None, rewiring_swaps_per_edge: int | None = None,
                   expected_arguments: dict | None = None, input_hashes: dict | None = None,
                   permutation_sha256: str | None = None, require_symmetric_rewiring: bool = False) -> tuple[np.ndarray | None, str]:
    """Lockbox predictions of one finished run, or None with the reason it cannot be used. With require_refit (the
    amendments of 8 October 2026) a run counts only if it kept the large leakage groups in training and was refitted on
    its training and validation perturbations, so it fitted on the same perturbations as the baselines. A run whose
    predictions are not finite, or not one row per lockbox perturbation and one column per symptom of the data, or
    (rewired variant) whose rewiring used another number of swaps per edge, is refused rather than scored. So is a run
    whose recorded arguments differ from expected_arguments (its configuration), whose recorded input hashes differ
    from input_hashes (the files scored now), or (permuted variant) whose label permutation is not the one the scorer
    draws for its seed."""
    if not (split_directory / "DONE").exists() or not (split_directory / "results.json").exists():
        return None, "not finished"
    results = json.loads((split_directory / "results.json").read_text())
    if results.get("test_perturbation_ids") != lockbox_ids:
        return None, "scored other perturbations than the lockbox"
    if (results.get("lockbox") or {}).get("sha256") != lockbox_sha256:
        return None, "trained against another lockbox file"
    if results.get("label_selection_sha256") != selection_sha256:
        return None, "trained on another label selection"
    if bool(results.get("labels_permuted")) != (variant == "permuted") or (results.get("rewiring") is not None) != (variant == "rewired"):
        return None, "its controls do not match the variant"
    if require_refit and not (results.get("refit") and (results.get("arguments") or {}).get("keep_large_groups_in_training")):
        return None, "trained without the refit or with the largest leakage group as its validation set"
    if symptoms is not None and results.get("symptoms") != symptoms:
        return None, "trained on another symptom list"
    if variant == "rewired" and rewiring_swaps_per_edge is not None and (results.get("rewiring") or {}).get("swaps_per_edge") != rewiring_swaps_per_edge:
        return None, f"its rewiring did not use {rewiring_swaps_per_edge} swaps per edge"
    if variant == "rewired" and require_symmetric_rewiring and not (results.get("rewiring") or {}).get("undirected_relations"):
        return None, "its rewiring did not keep the relations stored in both directions symmetric (--keep-reciprocated-relations-symmetric)"
    if expected_arguments is not None:
        recorded = results.get("arguments") or {}
        differing = sorted(key for key, value in expected_arguments.items() if recorded.get(key, value) != value)
        if differing:
            return None, "trained with other arguments than its configuration: " + ", ".join(f"{key} {recorded[key]!r} (expected {expected_arguments[key]!r})" for key in differing)
    if input_hashes is not None:
        for key, value in input_hashes.items():
            if results.get(key) != value:
                return None, f"trained on other input files than those scored ({key})"
    if variant == "permuted" and permutation_sha256 is not None and results.get("permutation_source_rows_sha256") != permutation_sha256:
        return None, "its label permutation differs from the one the scorer draws for its seed"
    predictions = np.load(split_directory / "test_predictions.npy")
    if predictions.shape != (len(lockbox_ids), len(symptoms) if symptoms is not None else predictions.shape[1]):
        return None, f"predictions of shape {predictions.shape}"
    if not np.isfinite(predictions).all():
        return None, "predictions that are not finite"
    return predictions, "ok"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v1_full_v2.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v1.json"))
    parser.add_argument("--run-root", type=Path, default=Path("runs/full"), help="holds <model>_<group-by>_confirmatory for every model")
    parser.add_argument("--models", nargs="*", default=list(CONFIRMATORY_MODELS))
    parser.add_argument("--baseline-root", type=Path, default=Path("runs/full/lockbox_baselines"), help="holds seed<k> from run_baselines.py --score-lockbox --seed <k>")
    parser.add_argument("--baselines", nargs="*", default=list(BASELINES))
    parser.add_argument("--seeds", type=int, nargs="*", default=[0, 1, 2, 3, 4])
    parser.add_argument("--num-bootstrap", type=int, default=4000)
    parser.add_argument("--bootstrap-seed", type=int, default=20261008)
    parser.add_argument("--bootstrap-unit", choices=["group", "perturbation"], default="group",
                        help="what the paired bootstrap resamples for the p-values and intervals that decide: whole leakage groups of the lockbox (the "
                             "user's decision of 8 October 2026) or single perturbations; the other is reported beside it as a sensitivity reading")
    parser.add_argument("--allow-asymmetric-rewiring", action="store_true",
                        help="accept rewired runs without --keep-reciprocated-relations-symmetric (tests on older runs only; the user adopted it on 8 October 2026)")
    parser.add_argument("--alpha", type=float, default=0.025, help="one-sided level of each model's H1 and then H2 (a two-sided 95 percent interval excluding zero)")
    parser.add_argument("--minimum-macro-difference", type=float, default=0.041,
                        help="smallest macro AUPRC difference that counts: the projected 95 percent half-width of a macro difference on a 20 percent hold-out")
    parser.add_argument("--minimum-micro-difference", type=float, default=0.028,
                        help="smallest micro AUPRC difference that counts: the projected 95 percent half-width of a micro difference on a 20 percent hold-out")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/confirmatory"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/confirmatory_results.md"))
    parser.add_argument("--rescore", action="store_true", help="score again although a SCORED marker exists (the output lists every scoring)")
    parser.add_argument("--allow-runs-without-refit", action="store_true",
                        help="smoke tests on runs trained before 8 October 2026 only: accept runs without --refit-on-validation and --keep-large-groups-in-training")
    parser.add_argument("--rewiring-swaps-per-edge", type=int, default=50, help="the swaps per edge every rewired run must have used (runs/full/confirmatory_lockbox.sh)")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="score although a model has a missing or refused run (that model is then not confirmed); without it the scorer stops before "
                             "writing anything, so a failed run can be resumed and the single scoring is not spent")
    arguments = parser.parse_args()

    marker = arguments.output_dir / "SCORED"
    previous_scorings = json.loads(marker.read_text()) if marker.exists() else []
    if previous_scorings and not arguments.rescore:
        raise SystemExit(f"the lockbox was scored already ({previous_scorings}); pass --rescore to score it again, which the output will list")

    label_selection = None if str(arguments.label_selection) in ("", "none") else arguments.label_selection  # none: the slice, for tests of this script
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=label_selection)
    in_lockbox = read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir)
    lockbox = json.loads(arguments.lockbox.read_text())
    lockbox_sha256, selection_sha256 = file_sha256(arguments.lockbox), (file_sha256(label_selection) if label_selection else None)
    rows = np.flatnonzero(in_lockbox)
    lockbox_ids = [data.perturbation_ids[i] for i in rows]
    label_mask = data.label_mask if data.label_mask is not None else np.ones_like(data.outcomes, dtype=bool)
    kept = data.outcomes * label_mask
    scored_columns = [column for column in range(len(data.symptoms)) if kept[:, column].sum() >= MINIMUM_POSITIVES_TO_SCORE]
    macro_columns = [data.symptoms.index(symptom) for symptom in lockbox["symptoms_in_the_macro_average"]]
    recomputed = [column for column in range(len(data.symptoms)) if kept[rows, column].sum() >= MINIMUM_POSITIVES_TO_SCORE]
    if sorted(macro_columns) != recomputed:
        raise SystemExit("the lockbox file's macro symptoms differ from the counts in the data")
    strata = degree_strata(data.perturbation_degrees)[rows]
    outcomes, mask = data.outcomes[rows], label_mask[rows]

    labellings = {"real": {seed: (outcomes, mask) for seed in arguments.seeds}, "permuted": {}}
    permutation_sha256_by_seed: dict[int, str] = {}
    for seed in arguments.seeds:
        source_row = degree_stratified_row_permutation(data.perturbation_degrees, random_seed=seed, partition=in_lockbox)
        permutation_sha256_by_seed[seed] = array_sha256(source_row)
        baseline_source = arguments.baseline_root / f"seed{seed}" / "permutation_source_rows.npy"
        if not baseline_source.exists() or not np.array_equal(np.load(baseline_source), source_row):
            raise SystemExit(f"{baseline_source} is missing or holds another permutation than seed {seed} draws")
        labellings["permuted"][seed] = (data.outcomes[source_row][rows], label_mask[source_row][rows])

    baseline_predictions: dict[str, dict[str, dict[int, np.ndarray]]] = {"real": {}, "permuted": {}}
    for name in arguments.baselines:
        for labelling, infix in (("real", "lockbox"), ("permuted", "lockbox_label_permutation")):
            baseline_predictions[labelling][name] = {}
            for seed in arguments.seeds:
                directory = arguments.baseline_root / f"seed{seed}"
                if json.loads((directory / "perturbation_ids.json").read_text()) != data.perturbation_ids:
                    raise SystemExit(f"{directory} holds other perturbations")
                baseline_predictions[labelling][name][seed] = np.load(directory / f"predictions_{infix}_{name}.npy")[rows]

    model_predictions: dict[str, dict[str, dict[int, np.ndarray]]] = {}
    missing_runs: dict[str, list[str]] = {}
    graph_and_evidence_hashes = {"graph_files_sha256": {name: file_sha256(arguments.graph_dir / f"{name}.parquet") for name in ("nodes", "edges")},
                                 "evidence_records_sha256": file_sha256(arguments.evidence_dir / "evidence_records.parquet")}
    for model in arguments.models:
        model_directory = arguments.run_root / f"{model}_{arguments.group_by}_confirmatory"
        model_predictions[model] = {variant: {} for variant in VARIANT_SUFFIXES}
        missing_runs[model] = []
        configuration_arguments = expected_run_arguments(model, 0, arguments.group_by, arguments.lockbox)
        input_hashes = {**graph_and_evidence_hashes, "node_descriptors_sha256": file_sha256(configuration_arguments.get("node_descriptors")),
                        "cell_class_weights_sha256": file_sha256(configuration_arguments.get("cell_class_weights"))}
        for variant, suffix in VARIANT_SUFFIXES.items():
            for seed in arguments.seeds:
                split_directory = model_directory / f"lockbox_seed{seed}{suffix}"
                predictions, reason = read_model_run(split_directory, lockbox_ids, lockbox_sha256, selection_sha256, variant,
                                                     require_refit=not arguments.allow_runs_without_refit, symptoms=data.symptoms,
                                                     rewiring_swaps_per_edge=arguments.rewiring_swaps_per_edge,
                                                     expected_arguments=expected_run_arguments(model, seed, arguments.group_by, arguments.lockbox),
                                                     input_hashes=input_hashes, permutation_sha256=permutation_sha256_by_seed.get(seed),
                                                     require_symmetric_rewiring=not arguments.allow_asymmetric_rewiring)
                if predictions is None:
                    missing_runs[model].append(f"{split_directory}: {reason}")
                else:
                    model_predictions[model][variant][seed] = predictions
    complete_models = [model for model in arguments.models if not missing_runs[model]]
    if len(complete_models) < len(arguments.models) and not arguments.allow_incomplete:
        raise SystemExit("not scored (no output, no SCORED marker): runs missing or refused:\n" + "\n".join(
            reason for model in arguments.models for reason in missing_runs[model]) + "\nresume the runs, or pass --allow-incomplete")

    def reading_scores(index: np.ndarray) -> dict:
        """Every score on the resampled rows index: baselines and models, per seed, labelling and reading."""
        strata_of_rows = strata[index]

        def scores(predictions: np.ndarray, labelling: str, seed: int) -> dict[str, float]:
            labels, labelled = labellings[labelling][seed]
            labels, labelled, chosen = labels[index], labelled[index], predictions[index]
            ranked = rank_normalise_within_groups(chosen, strata_of_rows)
            return {"macro": macro_auprc(chosen, labels, labelled, macro_columns), "micro": micro_auprc(chosen, labels, labelled, scored_columns),
                    "macro_within_degree_strata": macro_auprc(ranked, labels, labelled, macro_columns),
                    "micro_within_degree_strata": micro_auprc(ranked, labels, labelled, scored_columns)}

        result = {"baselines": {labelling: {name: {seed: scores(p, labelling, seed) for seed, p in by_seed.items()} for name, by_seed in by_name.items()}
                                for labelling, by_name in baseline_predictions.items()},
                  "models": {}}
        for model in complete_models:
            result["models"][model] = {variant: {seed: scores(p, "permuted" if variant == "permuted" else "real", seed) for seed, p in by_seed.items()}
                                       for variant, by_seed in model_predictions[model].items()}
        return result

    def seed_mean(per_seed: dict[int, dict[str, float]], reading: str) -> float:
        return float(np.mean([per_seed[seed][reading] for seed in arguments.seeds]))

    observed = reading_scores(np.arange(len(rows)))
    base_readings = ("macro", "micro", "macro_within_degree_strata", "micro_within_degree_strata")
    best_baseline = {reading: max(arguments.baselines, key=lambda name: seed_mean(observed["baselines"]["real"][name], reading)) for reading in base_readings}

    def differences(scored: dict, model: str) -> dict[str, float]:
        per_variant = scored["models"][model]
        real_baselines, permuted_baselines = scored["baselines"]["real"], scored["baselines"]["permuted"]
        values = {}
        for reading in base_readings:
            values[reading] = float(np.mean([per_variant["real"][seed][reading] - real_baselines[best_baseline[reading]][seed][reading] for seed in arguments.seeds]))
        for reading in ("macro", "micro"):
            permuted_advantage = np.mean([per_variant["permuted"][seed][reading] - permuted_baselines[best_baseline[reading]][seed][reading] for seed in arguments.seeds])
            values[f"{reading}_permutation"] = float(values[reading] - permuted_advantage)
            values[f"{reading}_rewiring"] = float(np.mean([per_variant["real"][seed][reading] - per_variant["rewired"][seed][reading] for seed in arguments.seeds]))
        return values

    point = {model: differences(observed, model) for model in complete_models}
    group_of_row = np.array([data.group_ids[i] for i in rows])
    members_of_group = [np.flatnonzero(group_of_row == group) for group in sorted(set(group_of_row.tolist()))]
    bootstrap_units = [arguments.bootstrap_unit] + [unit for unit in ("group", "perturbation") if unit != arguments.bootstrap_unit]
    resampled_by_unit = {}
    for unit_index, unit in enumerate(bootstrap_units):  # the deciding unit first, then the sensitivity reading
        generator = np.random.default_rng(arguments.bootstrap_seed + unit_index)
        resampled_by_unit[unit] = {model: {reading: [] for reading in ALL_READINGS} for model in complete_models}
        for _ in range(arguments.num_bootstrap if complete_models else 0):
            scored = reading_scores(bootstrap_rows(unit, members_of_group, len(rows), generator))
            for model in complete_models:
                for reading, value in differences(scored, model).items():
                    resampled_by_unit[unit][model][reading].append(value)
    resampled = resampled_by_unit[arguments.bootstrap_unit]

    def interval_and_p(values: list[float]) -> dict:
        values = np.array(values, dtype=float)
        finite = values[np.isfinite(values)]
        return {"p_one_sided": one_sided_p(values), "lower_95": float(np.quantile(finite, 0.025)) if finite.size else float("nan"),
                "upper_95": float(np.quantile(finite, 0.975)) if finite.size else float("nan"), "num_resamples_defined": int(finite.size)}

    entries = {}
    for model in arguments.models:
        entry = {"complete": model in complete_models, "missing_runs": missing_runs[model]}
        if model in complete_models:
            readings = {reading: {"difference": point[model][reading], **interval_and_p(resampled[model][reading])} for reading in ALL_READINGS}
            entry["readings"] = readings
            sensitivity_unit = bootstrap_units[1]
            entry["sensitivity_bootstrap"] = {"unit": sensitivity_unit, "readings": {reading: interval_and_p(resampled_by_unit[sensitivity_unit][model][reading]) for reading in ALL_READINGS}}
            entry["sensitivity_bootstrap"]["p_h1"] = max(entry["sensitivity_bootstrap"]["readings"][reading]["p_one_sided"] for reading in H1_READINGS)
            entry["meets_minimum_differences"] = bool(readings["macro"]["difference"] >= arguments.minimum_macro_difference
                                                      and readings["micro"]["difference"] >= arguments.minimum_micro_difference)
            entry["p_h1"] = max(readings[reading]["p_one_sided"] for reading in H1_READINGS)
            entry["p_h2"] = max([entry["p_h1"]] + [readings[reading]["p_one_sided"] for reading in H2_EXTRA_READINGS])
            entry["model_scores"] = {variant: {reading: seed_mean(observed["models"][model][variant], reading) for reading in base_readings} for variant in VARIANT_SUFFIXES}
            entry["seed_scores_real"] = {str(seed): observed["models"][model]["real"][seed] for seed in arguments.seeds}
        entry.update(fixed_sequence_decisions(entry["complete"], entry.get("p_h1", 1.0), entry.get("p_h2", 1.0), entry.get("meets_minimum_differences", False), arguments.alpha))
        entries[model] = entry

    scoring = {"scored_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
               "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
               "tracked_changes": subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True).stdout.splitlines(),
               "incomplete_allowed": arguments.allow_incomplete}
    output = {
        "scorings": previous_scorings + [scoring], "lockbox": str(arguments.lockbox), "lockbox_sha256": lockbox_sha256, "label_selection_sha256": selection_sha256,
        "num_lockbox_perturbations": int(len(rows)), "macro_symptoms": [data.symptoms[c] for c in macro_columns], "micro_symptoms": [data.symptoms[c] for c in scored_columns],
        "seeds": arguments.seeds, "num_bootstrap": arguments.num_bootstrap, "bootstrap_unit": arguments.bootstrap_unit, "num_lockbox_groups": len(members_of_group),
        "alpha_one_sided_per_model": arguments.alpha,
        "minimum_macro_difference": arguments.minimum_macro_difference, "minimum_micro_difference": arguments.minimum_micro_difference, "runs_without_refit_allowed": arguments.allow_runs_without_refit,
        "best_baseline": best_baseline,
        "baseline_scores": {labelling: {name: {reading: seed_mean(by_seed, reading) for reading in base_readings} for name, by_seed in by_name.items()}
                            for labelling, by_name in observed["baselines"].items()},
        "models": entries,
    }
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "lockbox_scores.json").write_text(json.dumps(output, indent=1) + "\n")
    marker.write_text(json.dumps(previous_scorings + [scoring]) + "\n")

    lines = ["# Confirmatory test on the lockbox (generated by experiments/score_confirmatory.py)", "",
             f"Lockbox {arguments.lockbox} ({len(rows)} perturbations, SHA-256 {lockbox_sha256[:12]}); label selection {str(selection_sha256)[:12]}; seeds {arguments.seeds}; "
             f"{arguments.num_bootstrap} paired bootstrap resamples over lockbox {'leakage groups (' + str(len(members_of_group)) + ')' if arguments.bootstrap_unit == 'group' else 'perturbations'}; "
             f"each model tested on its own, H1 and then H2 at one-sided {arguments.alpha}; "
             f"minimum macro difference {arguments.minimum_macro_difference}, minimum micro difference {arguments.minimum_micro_difference}. With {len(entries)} models "
             f"and no correction across them, the chance that at least one is confirmed by luck is at most {1 - (1 - arguments.alpha) ** len(entries):.3f} if they were "
             f"independent, and less since they share the data and baselines. Scorings: {len(previous_scorings) + 1}.", "",
             "Best baseline per reading (the highest seed-mean lockbox score on real labels; no model score enters the choice): "
             + "; ".join(f"{reading} {name} ({output['baseline_scores']['real'][name][reading]:.3f})" for reading, name in best_baseline.items()) + ".", "",
             "| model | complete | macro | micro | macro within strata | micro within strata | macro rewiring | micro rewiring | macro permutation (secondary) | micro permutation (secondary) | p(H1) | H1 confirmed | p(H2) | H2 confirmed |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for model, entry in entries.items():
        if not entry["complete"]:
            lines.append(f"| {model} | no ({len(entry['missing_runs'])} runs missing) |" + " n/a |" * 8 + " 1 | no | 1 | no |")
            continue
        cells = [f"{entry['readings'][reading]['difference']:+.3f} [{entry['readings'][reading]['lower_95']:+.3f}, {entry['readings'][reading]['upper_95']:+.3f}]"
                 for reading in ALL_READINGS]
        lines.append(f"| {model} | yes | " + " | ".join(cells) + f" | {entry['p_h1']:.4f} | {'yes' if entry['H1_confirmed'] else 'no'} | {entry['p_h2']:.4f} | {'yes' if entry['H2_confirmed'] else 'no'} |")
    lines += ["", f"Sensitivity reading, resampling {bootstrap_units[1]}s instead: "
              + "; ".join(f"{model} p(H1) {entry['sensitivity_bootstrap']['p_h1']:.4f}" for model, entry in entries.items() if entry["complete"]) + "."]
    lines += ["", "Each cell: seed-mean difference with the 2.5 and 97.5 percentiles of its bootstrap distribution. H1 takes the first four readings, H2 adds "
              "the two rewiring readings. Rewiring: the model on the real graph minus the same model on its rewired graph. Permutation (secondary, in neither "
              "hypothesis): the advantage on the real labels minus the advantage when model and baseline are trained and scored on labels permuted within "
              "degree strata."]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
