"""Score each noisy-OR run only on the held-out genes it gives a distinct prediction, beside the share of genes that is.

Without a degree offset the noisy-OR head gives two perturbations the same prediction only when their module activations
agree, and with sparse support gates many held-out genes reach no module and all receive the bias-level prediction
(docs/b6_module_diagnosis.md, held_out_within_tolerance_of_modal). The model makes no ranking among those genes. The
main tables score every held-out gene; this script gives the reading of a classifier that may abstain: accuracy among
the genes it decides, beside its coverage (the share it decides).

Per fold, the undecided genes are the largest group of held-out genes whose prediction rows lie within
PREDICTION_TIE_TOLERANCE, in every symptom, of one member of the group (the member with the most such neighbours); a fold
with no group of two or more has no undecided genes. Every other held-out gene is decided. Genes are pooled across folds,
since one fold holds too few to score. The decided set is chosen from the run's own predictions, never from labels, so
an abstaining model could use it at prediction time; but it can also select genes that are easier for every method, so
each comparison model is scored on the same rows, and popularity (a constant per symptom within a fold) shows the
stratum's base rate. Within the undecided genes the run's prediction is constant within a fold, so its pooled score
there only measures between-fold differences in the bias.

Usage: python experiments/score_by_module_coverage.py [--runs NAME ...] [--comparators NAME ...]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.compare_twin_runs import read_run
from experiments.score_by_pathway_split import BASELINES, macro_auprc, pathway_split_strata, pooled_run_predictions
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import paired_bootstrap_macro_difference

PREDICTION_TIE_TOLERANCE = 1e-6  # the tolerance of experiments/diagnose_noisy_or_modules.py
MINIMUM_POSITIVES_TO_SCORE = 5
DEFAULT_RUNS = ("b6_linear_response_gate_time_scales_cofactors", "b6_mechanistic_gate_time_scales")
DEFAULT_COMPARATORS = ("b3_linear_response_cofactors", "b3_typed_nodes")


def largest_tie_group(scores: np.ndarray, tolerance: float = PREDICTION_TIE_TOLERANCE) -> np.ndarray:
    """Boolean mask over rows: the rows within tolerance (max over columns) of the row with the most such neighbours,
    or no rows when no two rows are within tolerance of each other."""
    if len(scores) < 2:
        return np.zeros(len(scores), dtype=bool)
    distances = np.abs(scores[:, None, :] - scores[None, :, :]).max(axis=2)
    neighbour_counts = (distances <= tolerance).sum(axis=1)
    center = int(neighbour_counts.argmax())
    group = distances[center] <= tolerance
    return group if group.sum() >= 2 else np.zeros(len(scores), dtype=bool)


def undecided_rows(run: str, data) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, int]]] | None:
    """(undecided mask over data's rows, pooled out-of-fold predictions, per-fold (fold, held-out, undecided) counts),
    or None when the run lacks a fold."""
    _, predictions = read_run(run)
    if predictions is None or predictions.fold.nunique() < 5:
        return None
    position_of = {perturbation: index for index, perturbation in enumerate(data.perturbation_ids)}
    undecided = np.zeros(len(data.perturbation_ids), dtype=bool)
    pooled = np.full(data.outcomes.shape, np.nan)
    fold_counts = []
    for fold, fold_predictions in predictions.groupby("fold"):
        scores = fold_predictions[data.symptoms].to_numpy(dtype=float)
        rows = [position_of[perturbation] for perturbation in fold_predictions.perturbation_id]
        group = largest_tie_group(scores)
        undecided[rows] = group
        pooled[rows] = scores
        fold_counts.append((int(fold), len(rows), int(group.sum())))
    if np.isnan(pooled).any():
        return None
    return undecided, pooled, fold_counts


def score_rows(name_of_run: str, predictions: dict[str, np.ndarray], outcomes: np.ndarray, rows: np.ndarray,
               num_bootstrap: int, seed: int) -> list[str]:
    """Table rows: every model's macro AUPRC on rows, its lift over popularity, and the run minus it, paired."""
    subset_outcomes = outcomes[rows]
    popularity_score, scored_symptoms = macro_auprc(predictions["popularity"][rows], subset_outcomes)
    if scored_symptoms == 0:
        return [f"No symptom has {MINIMUM_POSITIVES_TO_SCORE} positives among these {len(rows)} genes; not scorable.", ""]
    lines = [f"{scored_symptoms} symptoms scored.", "",
             f"| model | macro AUPRC | lift over popularity | {name_of_run} minus this model, paired [95% CI] |", "|---|---|---|---|"]
    for name, model_predictions in predictions.items():
        score, _ = macro_auprc(model_predictions[rows], subset_outcomes)
        if name == name_of_run:
            paired = "this run"
        else:
            comparison = paired_bootstrap_macro_difference(predictions[name_of_run][rows], model_predictions[rows], subset_outcomes,
                                                           num_bootstrap=num_bootstrap, random_seed=seed,
                                                           minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
            paired = f"{comparison['difference']:+.3f} [{comparison['lower']:+.3f}, {comparison['upper']:+.3f}]"
        lines.append(f"| {name} | {score:.3f} | {score - popularity_score:+.3f} | {paired} |")
    return lines + [""]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--baseline-dir", type=Path, default=Path("runs/baselines_disease_cluster"))
    parser.add_argument("--runs", nargs="*", default=list(DEFAULT_RUNS))
    parser.add_argument("--comparators", nargs="*", default=list(DEFAULT_COMPARATORS),
                        help="other runs scored on the same rows as each noisy-OR run")
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("docs/score_by_module_coverage.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    position_of = {perturbation: index for index, perturbation in enumerate(data.perturbation_ids)}
    outcomes = (data.outcomes >= 0.5).astype(float)
    baseline_ids = json.loads((arguments.baseline_dir / "perturbation_ids.json").read_text())
    if baseline_ids != list(data.perturbation_ids):
        raise ValueError("baseline predictions are not in the evidence table's perturbation order")
    shared_predictions = {name: np.load(arguments.baseline_dir / f"predictions_grouped_{name}.npy") for name in BASELINES}
    skipped_comparators = []
    for comparator in arguments.comparators:
        pooled = next((found for found in (pooled_run_predictions(directory, position_of, data.outcomes.shape)
                                           for directory in (Path("runs/encoder") / f"{comparator}_disease_cluster",
                                                             Path("runs") / f"{comparator}_disease_cluster")
                                           if directory.is_dir()) if found is not None), None)
        if pooled is None:
            skipped_comparators.append(comparator)
        else:
            shared_predictions[comparator] = pooled
    pathway_strata = pathway_split_strata(data, arguments.num_folds, arguments.seed)
    degrees = np.asarray(data.perturbation_degrees, dtype=float)
    positives_per_gene = outcomes.sum(axis=1)

    lines = [
        "# Scores among the genes a noisy-OR run decides, beside its coverage (generated by experiments/score_by_module_coverage.py)",
        "",
        f"{len(data.perturbation_ids)} perturbations, {arguments.num_folds} disease-cluster folds (seed {arguments.seed}); out-of-fold "
        "predictions pooled across folds. Per fold, a held-out gene is undecided when it belongs to the largest group of held-out genes "
        f"whose prediction rows lie within {PREDICTION_TIE_TOLERANCE:g} of one member in every symptom; the rest are decided. The run "
        "chooses the decided set from its own predictions, without labels. Macro AUPRC over symptoms with at least "
        f"{MINIMUM_POSITIVES_TO_SCORE} positives among the rows; paired bootstrap over perturbations ({arguments.num_bootstrap} "
        "resamples), with no correction for the number of intervals. Popularity is a constant per symptom within a fold, so it scores "
        "the rows' base rate, and the lift column corrects for subsets with different base rates. Within the undecided genes the "
        "run's prediction is constant within a fold, so its score there only reflects between-fold differences in the bias."
        + (f" Comparators skipped, fewer than five folds: {', '.join(skipped_comparators)}." if skipped_comparators else ""),
        "",
    ]
    for run in arguments.runs:
        found = undecided_rows(run, data)
        if found is None:
            lines += [f"## {run}: fewer than five folds, skipped", ""]
            continue
        undecided, run_predictions, fold_counts = found
        decided = ~undecided
        predictions = {run: run_predictions, **shared_predictions}
        lines += [f"## {run}", "",
                  f"Coverage: {int(decided.sum())} of {len(decided)} genes decided ({decided.mean():.1%}).", "",
                  "| fold | held-out genes | undecided | undecided share |", "|---|---|---|---|"]
        lines += [f"| {fold} | {held_out} | {tied} | {tied / held_out:.3f} |" for fold, held_out, tied in fold_counts]
        lines += ["", "| genes | count | median metabolic-graph degree | mean labelled positives per gene |", "|---|---|---|---|"]
        for label, mask in (("decided", decided), ("undecided", undecided)):
            lines.append(f"| {label} | {int(mask.sum())} | {np.median(degrees[mask]):.0f} | {positives_per_gene[mask].mean():.2f} |")
        lines += ["", "| pathway stratum (experiments/score_by_pathway_split.py) | genes | decided | decided share |", "|---|---|---|---|"]
        for stratum_name, stratum_rows in pathway_strata.items():
            lines.append(f"| {stratum_name} | {len(stratum_rows)} | {int(decided[stratum_rows].sum())} | "
                         f"{decided[stratum_rows].mean():.3f} |")
        lines.append("")
        for label, rows in (("decided genes", np.where(decided)[0]), ("undecided genes", np.where(undecided)[0]),
                            ("all genes", np.arange(len(decided)))):
            lines += [f"### {run}: {label} ({len(rows)})", ""]
            lines += score_rows(run, predictions, outcomes, rows, arguments.num_bootstrap, arguments.seed)
    arguments.output.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
