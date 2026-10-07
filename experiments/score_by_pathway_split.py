"""Score the out-of-fold predictions separately for genes whose pathway the split divides and genes whose pathway it does not.

The disease-cluster split puts most held-out genes in a pathway that also has training genes, and those pairs share
symptom labels above baseline (docs/fold_mechanism_sharing.md). If a model's advantage comes from that, it should be
concentrated in the split-pathway stratum. Each held-out gene goes to one stratum, by its primary Human-GEM subsystem
(the subsystem most of its reactions belong to, catch-all subsystems excluded):

- split pathway: a gene in another fold has the same primary subsystem, so the pathway has labelled training members
- unsplit: the primary subsystem has no gene in another fold, or the gene has no pathway subsystem at all

AUPRC is not comparable across strata with different base rates, so each stratum reports every model's macro AUPRC
next to popularity's (a constant per symptom, which scores the stratum's base rate) and the lift over it, and a paired
bootstrap of each model against the random walk inside the stratum. Predictions are pooled across folds, since one
fold holds too few unsplit genes to score.

Usage: python experiments/score_by_pathway_split.py [--runs NAME ...] [--output docs/score_by_pathway_split.md]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    CATCH_ALL_SUBSYSTEMS,
    assign_grouped_folds,
    primary_subsystem_by_gene_node,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    paired_bootstrap_macro_difference,
    per_symptom_auprc,
)

MINIMUM_POSITIVES_TO_SCORE = 5
BASELINES = ("popularity", "degree_popularity", "random_walk_with_restart")
DEFAULT_RUNS = ("b3_typed_nodes", "b6_mechanistic_gate_time_scales", "b6_linear_response_gate_time_scales_cofactors",
                "b3_linear_response_cofactors")


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray) -> tuple[float, int]:
    """Macro AUPRC over symptoms with at least MINIMUM_POSITIVES_TO_SCORE positives and one negative; returns it and the count."""
    values = [per_symptom_auprc(predictions, outcomes, symptom) for symptom in range(outcomes.shape[1])
              if MINIMUM_POSITIVES_TO_SCORE <= outcomes[:, symptom].sum() < outcomes.shape[0]]
    return (float(np.mean(values)) if values else float("nan")), len(values)


def pathway_split_strata(data, num_folds: int, seed: int) -> dict[str, np.ndarray]:
    """Row indices of each stratum in the module docstring, keyed by the stratum's heading."""
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, num_folds, seed)
    primary_subsystem = primary_subsystem_by_gene_node(data.node_subsystem, data.edge_source, data.edge_target,
                                                       data.edge_relation, data.relation_types.index("catalyzed_by"))

    def pathway_of(row: int) -> str | None:
        seeds = data.perturbation_seeds[row]
        subsystem = primary_subsystem.get(int(seeds[0])) if len(seeds) == 1 else None
        return subsystem if subsystem and subsystem not in CATCH_ALL_SUBSYSTEMS else None

    pathway_by_row = [pathway_of(row) for row in range(len(data.perturbation_ids))]
    folds_by_pathway: dict[str, set[int]] = {}
    for row, pathway in enumerate(pathway_by_row):
        if pathway:
            folds_by_pathway.setdefault(pathway, set()).add(fold_by_perturbation[data.perturbation_ids[row]])
    is_split = np.array([bool(pathway) and len(folds_by_pathway[pathway] - {fold_by_perturbation[data.perturbation_ids[row]]}) > 0
                         for row, pathway in enumerate(pathway_by_row)])
    has_pathway = np.array([pathway is not None for pathway in pathway_by_row])
    return {"split pathway": np.where(is_split)[0],
            "unsplit: pathway absent from training": np.where(has_pathway & ~is_split)[0],
            "unsplit: no pathway subsystem (transporters, catch-all and non-metabolic genes)": np.where(~has_pathway)[0],
            "unsplit, both kinds together": np.where(~is_split)[0]}


def pooled_run_predictions(run_directory: Path, position_of: dict[str, int], shape: tuple[int, int]) -> np.ndarray | None:
    pooled = np.full(shape, np.nan)
    fold_directories = [path for path in sorted(run_directory.glob("fold*_seed0")) if (path / "DONE").exists()]
    if len(fold_directories) < 5:
        return None
    for fold_directory in fold_directories:
        results = json.loads((fold_directory / "results.json").read_text())
        rows = [position_of[perturbation] for perturbation in results["test_perturbation_ids"]]
        pooled[rows] = np.load(fold_directory / "test_predictions.npy")
    return None if np.isnan(pooled).any() else pooled


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--baseline-dir", type=Path, default=Path("runs/baselines_disease_cluster"))
    parser.add_argument("--run-root", type=Path, default=Path("runs/encoder"))
    parser.add_argument("--runs", nargs="*", default=list(DEFAULT_RUNS))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("docs/score_by_pathway_split.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    position_of = {perturbation: index for index, perturbation in enumerate(data.perturbation_ids)}
    strata = pathway_split_strata(data, arguments.num_folds, arguments.seed)

    baseline_ids = json.loads((arguments.baseline_dir / "perturbation_ids.json").read_text())
    if baseline_ids != list(data.perturbation_ids):
        raise ValueError("baseline predictions are not in the evidence table's perturbation order")
    predictions = {name: np.load(arguments.baseline_dir / f"predictions_grouped_{name}.npy") for name in BASELINES}
    skipped_runs = []
    for run in arguments.runs:
        # slice runs live under runs/encoder or, for the earliest configurations, directly under runs/
        candidate_directories = [arguments.run_root / f"{run}_disease_cluster", Path("runs") / f"{run}_disease_cluster"]
        pooled = next((found for found in (pooled_run_predictions(directory, position_of, data.outcomes.shape)
                                           for directory in candidate_directories if directory.is_dir()) if found is not None), None)
        if pooled is None:
            skipped_runs.append(run)
        else:
            predictions[run] = pooled
    outcomes = (data.outcomes >= 0.5).astype(float)

    lines = [
        "# Scores for genes whose pathway the split divides, and genes whose it does not (generated by experiments/score_by_pathway_split.py)",
        "",
        f"{len(data.perturbation_ids)} perturbations, {arguments.num_folds} disease-cluster folds (seed {arguments.seed}); out-of-fold "
        f"predictions pooled across folds; macro AUPRC over symptoms with at least {MINIMUM_POSITIVES_TO_SCORE} positives in the "
        f"stratum; paired bootstrap over perturbations ({arguments.num_bootstrap} resamples). A gene is in the split-pathway stratum when "
        "a gene in another fold has the same primary Human-GEM subsystem (catch-all subsystems excluded)."
        + (f" Skipped, fewer than five folds: {', '.join(skipped_runs)}." if skipped_runs else ""),
        "",
    ]
    for stratum_name, rows in strata.items():
        stratum_outcomes = outcomes[rows]
        popularity_score, scored_symptoms = macro_auprc(predictions["popularity"][rows], stratum_outcomes)
        if scored_symptoms == 0:
            lines += [f"## {stratum_name}: {len(rows)} genes, no symptom has {MINIMUM_POSITIVES_TO_SCORE} positives; not scorable", ""]
            continue
        lines += [f"## {stratum_name}: {len(rows)} genes, {scored_symptoms} symptoms scored", "",
                  "| model | macro AUPRC | lift over popularity | minus random walk, paired [95% CI] |", "|---|---|---|---|"]
        for name, model_predictions in predictions.items():
            score, _ = macro_auprc(model_predictions[rows], stratum_outcomes)
            if name == "random_walk_with_restart":
                paired = "reference"
            else:
                comparison = paired_bootstrap_macro_difference(model_predictions[rows], predictions["random_walk_with_restart"][rows],
                                                               stratum_outcomes, num_bootstrap=arguments.num_bootstrap, random_seed=arguments.seed,
                                                               minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
                paired = f"{comparison['difference']:+.3f} [{comparison['lower']:+.3f}, {comparison['upper']:+.3f}]"
            lines.append(f"| {name} | {score:.3f} | {score - popularity_score:+.3f} | {paired} |")
        lines.append("")
    arguments.output.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
