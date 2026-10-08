"""A trained readout on random-walk distributions, against the untrained random walk.

Every trained configuration on the slice ties or trails the random walk with restart (docs/research_summary.md, row
5), but the walk is untrained while every configuration fits a readout to the labels, so the comparison changes the
encoder and the training together. Huang et al. (2024, Nature Medicine, doi:10.1038/s41591-024-03233-x) has the same
confound the other way round: a trained typed model against untrained diffusion statistics
(docs/literature_appraisal_staging.md, slot 17). This control holds the diffusion fixed and adds only supervision: if a
readout trained on the walk's own distributions cannot beat the walk on 451 genes, the labels may be too few for any
trained readout to beat it, and the encoders' deficit is not evidence against their mechanism; if it can, the deficit
lies in the encoders.

Features, fixed before scoring: each gene's walk-with-restart distribution over reaction and metabolite nodes (genes
and currency metabolites excluded), from experiments/score_directed_random_walk.py, square-rooted (the Hellinger
transform, which tempers the mass near the seed), reduced to NUM_COMPONENTS by a truncated SVD fitted on the training
genes of each fold only, and standardised. Readout: one L2-regularised logistic regression per symptom, its inverse
regularisation strength chosen from REGULARISATION_GRID by an inner INNER_FOLDS-fold cross-validation on the training
genes scored by average precision. A symptom with fewer than MINIMUM_TRAINING_POSITIVES training positives gets its
training base rate. The labels never enter the features, so nothing leaks across folds; the inner folds are not
grouped by disease cluster, which can only make the chosen strength optimistic, not the test score.

Usage: python experiments/score_trained_diffusion_readout.py [--output docs/trained_diffusion_readout.md]
"""
from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler

from experiments.compare_twin_runs import describe_interval, fold_t_interval
from experiments.score_by_pathway_split import macro_auprc
from experiments.score_directed_random_walk import arm_edges, walk_distributions
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    degree_strata,
    paired_bootstrap_macro_difference,
    rank_normalise_within_groups,
)

NUM_COMPONENTS = 32
REGULARISATION_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)
INNER_FOLDS = 3
MINIMUM_TRAINING_POSITIVES = 5
MINIMUM_POSITIVES_TO_SCORE = 5
WALKS = ("undirected", "downstream only")


def trained_readout_predictions(features: np.ndarray, outcomes: np.ndarray, train: np.ndarray, test: np.ndarray, seed: int) -> np.ndarray:
    """[test rows, symptoms] probabilities from per-symptom logistic regressions fitted on the training rows."""
    reducer = TruncatedSVD(n_components=NUM_COMPONENTS, random_state=seed).fit(features[train])
    scaler = StandardScaler().fit(reducer.transform(features[train]))
    train_features = scaler.transform(reducer.transform(features[train]))
    test_features = scaler.transform(reducer.transform(features[test]))
    predictions = np.zeros((int(test.sum()), outcomes.shape[1]))
    for symptom in range(outcomes.shape[1]):
        labels = outcomes[train, symptom]
        if labels.sum() < MINIMUM_TRAINING_POSITIVES or labels.sum() > len(labels) - INNER_FOLDS:
            predictions[:, symptom] = labels.mean()
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            model = LogisticRegressionCV(Cs=list(REGULARISATION_GRID), cv=INNER_FOLDS, scoring="average_precision",
                                         max_iter=5000, random_state=seed).fit(train_features, labels)
        predictions[:, symptom] = model.predict_proba(test_features)[:, 1]
    return predictions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--baseline-dir", type=Path, default=Path("runs/baselines_disease_cluster"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("docs/trained_diffusion_readout.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster")
    outcomes = (data.outcomes >= 0.5).astype(float)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_row = np.array([fold_by_perturbation[perturbation] for perturbation in data.perturbation_ids])
    strata = fold_of_row * 1000 + degree_strata(data.perturbation_degrees_for_strata)
    profile_nodes = (np.asarray(data.node_types) != "gene") & ~np.asarray(data.is_currency, dtype=bool)

    walk_predictions = np.load(arguments.baseline_dir / "predictions_grouped_random_walk_with_restart.npy")
    predictions = {"random walk with restart (B1, untrained)": walk_predictions,
                   "popularity": np.load(arguments.baseline_dir / "predictions_grouped_popularity.npy")}
    for walk in WALKS:
        source, target = arm_edges(data, walk)
        features = np.sqrt(walk_distributions(len(data.node_ids), source, target, list(data.perturbation_seeds))[:, profile_nodes])
        walk_readout = np.zeros_like(outcomes)
        for fold in range(arguments.num_folds):
            test, train = fold_of_row == fold, fold_of_row != fold
            walk_readout[test] = trained_readout_predictions(features, outcomes, train, test, arguments.seed)
        predictions[f"trained readout on the {walk} walk"] = walk_readout

    def per_fold_scores(model_predictions: np.ndarray) -> np.ndarray:
        return np.array([macro_auprc(model_predictions[fold_of_row == fold], outcomes[fold_of_row == fold])[0]
                         for fold in range(arguments.num_folds)])

    reference_name = "random walk with restart (B1, untrained)"
    reference_per_fold = per_fold_scores(walk_predictions)
    lines = [
        "# A trained readout on random-walk distributions (generated by experiments/score_trained_diffusion_readout.py)",
        "",
        f"{len(data.perturbation_ids)} perturbations, {arguments.num_folds} disease-cluster folds (seed {arguments.seed}). Features: the walk "
        f"distribution over reaction and metabolite nodes, square-rooted, reduced to {NUM_COMPONENTS} components by a truncated SVD fitted "
        f"on each fold's training genes, standardised. Readout: per-symptom L2 logistic regression, strength chosen from "
        f"{list(REGULARISATION_GRID)} by {INNER_FOLDS}-fold inner cross-validation on average precision; symptoms with fewer than "
        f"{MINIMUM_TRAINING_POSITIVES} training positives get the training base rate. Differences are against the untrained walk "
        "(B1) and are read three ways as in docs/twin_comparisons.md: per fold (t interval, 4 df), pooled (paired bootstrap over "
        f"perturbations, {arguments.num_bootstrap} resamples) and pooled within fold-by-degree strata; 95 percent intervals, not corrected.",
        "",
        "| model | pooled macro AUPRC | per-fold mean (SD) | minus B1, per fold [95% t] | minus B1, pooled | minus B1, within degree strata |",
        "|---|---|---|---|---|---|",
    ]
    for name, model_predictions in predictions.items():
        pooled_score, _ = macro_auprc(model_predictions, outcomes)
        per_fold = per_fold_scores(model_predictions)
        if name == reference_name:
            cells = ["reference"] * 3
        else:
            mean_difference, lower, upper = fold_t_interval(per_fold - reference_per_fold)
            pooled = paired_bootstrap_macro_difference(model_predictions, walk_predictions, outcomes, num_bootstrap=arguments.num_bootstrap,
                                                       random_seed=arguments.seed, minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
            within_strata = paired_bootstrap_macro_difference(rank_normalise_within_groups(model_predictions, strata),
                                                              rank_normalise_within_groups(walk_predictions, strata), outcomes,
                                                              num_bootstrap=arguments.num_bootstrap, random_seed=arguments.seed,
                                                              minimum_positives=MINIMUM_POSITIVES_TO_SCORE)
            cells = [f"{mean_difference:+.3f} [{lower:+.3f}, {upper:+.3f}], {describe_interval(lower, upper)}"]
            cells += [f"{reading['difference']:+.3f} [{reading['lower']:+.3f}, {reading['upper']:+.3f}], "
                      f"{describe_interval(reading['lower'], reading['upper'])}" for reading in (pooled, within_strata)]
        lines.append(f"| {name} | {pooled_score:.3f} | {per_fold.mean():.3f} ({per_fold.std(ddof=1):.3f}) | {' | '.join(cells)} |")
    lines.append("")
    arguments.output.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
