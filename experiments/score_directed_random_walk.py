"""Random walks that follow edge direction, without signs, against the undirected walk the null result is measured on.

The baseline random walk with restart (B1) symmetrises the graph, so the signed linear-response encoder has only been
compared with a baseline that has neither direction nor sign. Silverbush et al. (2019, Nature Communications,
doi:10.1038/s41467-019-10887-6) report that diffusion over an oriented network ranks drug targets better than diffusion
over an unoriented one. This script asks the same question on the slice, with no signs anywhere.

B1 scores a held-out gene by the walk mass that lands on the training genes labelled with the symptom. A walk that
follows the metabolic graph's direction never returns to a gene (genes have out-edges only), so that readout is
undefined for it. Every arm here therefore uses one readout, fixed before scoring: the walk's stationary distribution
over reaction and metabolite nodes (genes and currency metabolites excluded), and the score of a held-out gene for a
symptom is the cosine similarity between its distribution and the mean distribution of the training genes labelled with
the symptom. The arms differ only in the edges the walk may follow:

- undirected: every edge in both directions (B1's adjacency)
- encoder edges: the directions the linear-response encoder propagates along, catalyzed_by (gene -> reaction),
  product_of (reaction -> metabolite), substrate_of (metabolite -> reaction) and the encoder's derived depletion edge
  (reaction -> substrate), all unsigned
- downstream only: catalyzed_by, substrate_of and product_of in their stored directions, so mass flows from substrates
  through reactions to products and never back

Reversible reactions keep their stored direction in every arm, as they do in the encoder. Mass at a node with no
outgoing edge returns to the restart distribution (the usual convention for personalised PageRank), so a gene whose
products dead-end is not scored lower for it. Restart probability and iteration count are B1's. B1 itself, with its
own readout, is reproduced from its saved predictions and listed for reference.

Usage: python experiments/score_directed_random_walk.py [--output docs/directed_random_walk.md]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sparse

from experiments.compare_twin_runs import describe_interval, fold_t_interval
from experiments.score_by_pathway_split import macro_auprc
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    degree_strata,
    paired_bootstrap_macro_difference,
    rank_normalise_within_groups,
)
from mechanistic_pathway_learning.models.baselines.random_walk_with_restart_baseline import (
    RandomWalkWithRestartBaseline,
    build_normalized_adjacency,
)

RESTART_PROBABILITY = 0.3  # B1's default in experiments/run_baselines.py
NUM_ITERATIONS = 50  # B1's default; (1 - 0.3) ** 50 is below 1e-7
MINIMUM_POSITIVES_TO_SCORE = 5
REFERENCE_ARM = "undirected"


def arm_edges(data, arm: str) -> tuple[np.ndarray, np.ndarray]:
    """(source, target) of the edges the arm's walk may follow, currency metabolites removed."""
    relation_index = {name: index for index, name in enumerate(data.relation_types)}
    keep = ~(data.is_currency[data.edge_source] | data.is_currency[data.edge_target])
    source, target, relation = data.edge_source[keep], data.edge_target[keep], data.edge_relation[keep]
    if arm == "undirected":
        return np.concatenate([source, target]), np.concatenate([target, source])
    stored = np.isin(relation, [relation_index["catalyzed_by"], relation_index["product_of"], relation_index["substrate_of"]])
    if arm == "downstream only":
        return source[stored], target[stored]
    if arm == "encoder edges":
        substrate = relation == relation_index["substrate_of"]
        return np.concatenate([source[stored], target[substrate]]), np.concatenate([target[stored], source[substrate]])
    raise ValueError(arm)


def walk_distributions(num_nodes: int, source: np.ndarray, target: np.ndarray, seeds_by_perturbation: list[np.ndarray]) -> np.ndarray:
    """[perturbations, nodes] stationary distributions of a walk with restart to each perturbation's seeds."""
    adjacency = sparse.coo_matrix((np.ones(len(source)), (target, source)), shape=(num_nodes, num_nodes)).tocsr()
    adjacency.data[:] = 1.0  # repeated edges count once
    out_degree = np.asarray(adjacency.sum(axis=0)).ravel()
    dangling = out_degree == 0
    transition = adjacency @ sparse.diags(1.0 / np.where(dangling, 1.0, out_degree))
    restart = np.zeros((num_nodes, len(seeds_by_perturbation)))
    for column, seeds in enumerate(seeds_by_perturbation):
        if len(seeds):
            restart[seeds, column] = 1.0 / len(seeds)
    state = restart.copy()
    for _ in range(NUM_ITERATIONS):
        dangling_mass = state[dangling].sum(axis=0, keepdims=True)
        state = (1.0 - RESTART_PROBABILITY) * (transition @ state + dangling_mass * restart) + RESTART_PROBABILITY * restart
    return state.T


def profile_scores(profiles: np.ndarray, outcomes: np.ndarray, train: np.ndarray, test: np.ndarray) -> np.ndarray:
    """Cosine similarity of each test profile with the mean profile of the training positives of each symptom."""
    unit_profiles = profiles / np.maximum(np.linalg.norm(profiles, axis=1, keepdims=True), 1e-300)
    scores = np.zeros((int(test.sum()), outcomes.shape[1]))
    for symptom in range(outcomes.shape[1]):
        positives = train & (outcomes[:, symptom] > 0)
        if positives.any():
            symptom_profile = profiles[positives].mean(axis=0)
            scores[:, symptom] = unit_profiles[test] @ (symptom_profile / max(np.linalg.norm(symptom_profile), 1e-300))
    return scores


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--baseline-dir", type=Path, default=Path("runs/baselines_disease_cluster"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("docs/directed_random_walk.md"))
    parser.add_argument("--predictions-dir", type=Path, default=Path("runs/directed_random_walk"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    outcomes = (data.outcomes >= 0.5).astype(float)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_row = np.array([fold_by_perturbation[perturbation] for perturbation in data.perturbation_ids])
    strata = fold_of_row * 1000 + degree_strata(data.perturbation_degrees)  # as in experiments/compare_twin_runs.py
    num_nodes = len(data.node_ids)
    profile_nodes = (np.asarray(data.node_types) != "gene") & ~np.asarray(data.is_currency, dtype=bool)

    predictions: dict[str, np.ndarray] = {}
    reproduced_b1 = np.zeros_like(outcomes)
    b1_adjacency = build_normalized_adjacency(num_nodes, data.edge_source, data.edge_target, np.where(data.is_currency)[0])
    for fold in range(arguments.num_folds):
        test, train = fold_of_row == fold, fold_of_row != fold
        model = RandomWalkWithRestartBaseline(restart_probability=RESTART_PROBABILITY, num_iterations=NUM_ITERATIONS).fit(
            b1_adjacency, [data.perturbation_seeds[row] for row in np.where(train)[0]], outcomes[train])
        reproduced_b1[test] = model.predict([data.perturbation_seeds[row] for row in np.where(test)[0]])
    saved_b1 = np.load(arguments.baseline_dir / "predictions_grouped_random_walk_with_restart.npy")
    b1_largest_difference = float(np.abs(saved_b1 - reproduced_b1).max())
    predictions["B1 (undirected, mass on training genes)"] = reproduced_b1

    reachable_share = {}
    for arm in ("undirected", "encoder edges", "downstream only"):
        source, target = arm_edges(data, arm)
        profiles = walk_distributions(num_nodes, source, target, list(data.perturbation_seeds))[:, profile_nodes]
        reachable_share[arm] = float(np.mean(profiles.sum(axis=1) > 0))
        arm_predictions = np.zeros_like(outcomes)
        for fold in range(arguments.num_folds):
            test, train = fold_of_row == fold, fold_of_row != fold
            arm_predictions[test] = profile_scores(profiles, outcomes, train, test)
        predictions[arm] = arm_predictions
    arguments.predictions_dir.mkdir(parents=True, exist_ok=True)
    for name, arm_predictions in predictions.items():
        np.save(arguments.predictions_dir / f"predictions_grouped_{name.split(' (')[0].replace(' ', '_')}.npy", arm_predictions)
    (arguments.predictions_dir / "perturbation_ids.json").write_text(json.dumps(list(data.perturbation_ids)))

    lines = [
        "# Random walks that follow edge direction, without signs (generated by experiments/score_directed_random_walk.py)",
        "",
        f"{len(data.perturbation_ids)} perturbations, {arguments.num_folds} disease-cluster folds (seed {arguments.seed}), restart "
        f"probability {RESTART_PROBABILITY}, {NUM_ITERATIONS} iterations, currency metabolites removed. The three arms share one readout "
        "(cosine similarity of the held-out gene's walk distribution over reaction and metabolite nodes with the mean distribution of "
        "the symptom's training positives) and differ only in the edges the walk may follow; see the script's docstring. B1, the "
        "committed baseline with its own readout, is reproduced for reference "
        f"(largest difference from its saved predictions {b1_largest_difference:.1e}). Macro AUPRC over symptoms with at least "
        f"{MINIMUM_POSITIVES_TO_SCORE} positives; per-fold standard deviation with n - 1; paired differences against the undirected "
        f"arm, read three ways as in docs/twin_comparisons.md: per fold (t interval, 4 df), pooled (bootstrap over perturbations, "
        f"{arguments.num_bootstrap} resamples) and pooled within fold-by-degree strata ({len(np.unique(strata))} strata of the "
        "metabolic-graph degree, scores rank-normalised inside each stratum, so ordering genes by degree earns nothing).",
        "",
        "| walk | genes with mass beyond their own node | pooled macro AUPRC | per-fold mean (SD) | minus undirected, per fold [95% t] | "
        "minus undirected, pooled [95% bootstrap] | minus undirected, within degree strata [95% bootstrap] |",
        "|---|---|---|---|---|---|---|",
    ]
    reference = predictions[REFERENCE_ARM]
    for name, arm_predictions in predictions.items():
        pooled_score, _ = macro_auprc(arm_predictions, outcomes)
        per_fold = np.array([macro_auprc(arm_predictions[fold_of_row == fold], outcomes[fold_of_row == fold])[0]
                             for fold in range(arguments.num_folds)])
        reach = f"{reachable_share[name]:.3f}" if name in reachable_share else "not applicable"
        if name == REFERENCE_ARM:
            per_fold_cell = pooled_cell = strata_cell = "reference"
        else:
            reference_per_fold = np.array([macro_auprc(reference[fold_of_row == fold], outcomes[fold_of_row == fold])[0]
                                           for fold in range(arguments.num_folds)])
            mean_difference, lower, upper = fold_t_interval(per_fold - reference_per_fold)
            per_fold_cell = f"{mean_difference:+.3f} [{lower:+.3f}, {upper:+.3f}], {describe_interval(lower, upper)}"
            comparison = paired_bootstrap_macro_difference(arm_predictions, reference, outcomes, num_bootstrap=arguments.num_bootstrap,
                                                           random_seed=arguments.seed, minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
            pooled_cell = (f"{comparison['difference']:+.3f} [{comparison['lower']:+.3f}, {comparison['upper']:+.3f}], "
                           f"{describe_interval(comparison['lower'], comparison['upper'])}")
            within_strata = paired_bootstrap_macro_difference(rank_normalise_within_groups(arm_predictions, strata),
                                                              rank_normalise_within_groups(reference, strata), outcomes,
                                                              num_bootstrap=arguments.num_bootstrap, random_seed=arguments.seed,
                                                              minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
            strata_cell = (f"{within_strata['difference']:+.3f} [{within_strata['lower']:+.3f}, {within_strata['upper']:+.3f}], "
                           f"{describe_interval(within_strata['lower'], within_strata['upper'])}")
        lines.append(f"| {name} | {reach} | {pooled_score:.3f} | {per_fold.mean():.3f} ({per_fold.std(ddof=1):.3f}) | {per_fold_cell} | "
                     f"{pooled_cell} | {strata_cell} |")
    lines.append("")
    arguments.output.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
