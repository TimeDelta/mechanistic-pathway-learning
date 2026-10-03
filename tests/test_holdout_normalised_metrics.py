"""Within-hold-out rank normalisation and stratified AUROC remove the base-rate artifact of pooled hold-out scores (review v0.4, finding 1)."""
from __future__ import annotations

import numpy as np

from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import per_symptom_auroc, rank_normalise_within_groups, stratified_auroc


def test_constant_predictor_scores_the_same_in_every_hold_out_and_stratified_auroc_is_half() -> None:
    group = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 1])
    predictions = np.array([[0.9]] * 3 + [[0.1]] * 7)  # one constant per hold-out: a fold-level base-rate signal, no within-fold ranking
    outcomes = np.array([[1.0], [0.0], [0.0]] + [[1.0]] * 4 + [[0.0]] * 3)
    normalised = rank_normalise_within_groups(predictions, group)
    assert np.allclose(normalised[:3, 0], 0.5) and np.allclose(normalised[3:, 0], 0.5)  # (n + 1) / 2 / (n + 1) in each hold-out
    assert per_symptom_auroc(predictions, outcomes, 0) < 0.5  # raw pooling rewards the fold with the lower base rate
    assert stratified_auroc(predictions, outcomes, group, 0) == 0.5


def test_ranks_are_tie_averaged_and_rows_outside_groups_are_untouched() -> None:
    group = np.array([0, 0, 0, -1])
    predictions = np.array([[0.2, 5.0], [0.2, 1.0], [0.7, 3.0], [0.5, 0.5]])
    normalised = rank_normalise_within_groups(predictions, group)
    assert np.allclose(normalised[:3, 0], [1.5 / 4, 1.5 / 4, 3 / 4]) and np.allclose(normalised[:3, 1], [3 / 4, 1 / 4, 2 / 4])
    assert np.allclose(normalised[3], [0.5, 0.5])
    outcomes = np.array([[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [0.0, 0.0]])
    assert stratified_auroc(predictions, outcomes, group, 0) == 1.0 and stratified_auroc(predictions, outcomes, group, 1) == 1.0
