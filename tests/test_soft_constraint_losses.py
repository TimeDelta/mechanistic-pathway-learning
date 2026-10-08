"""Tests for the soft constraint losses."""
import pytest
import torch

from mechanistic_pathway_learning.models.soft_constraint_losses import (
    anchored_link_prior_penalty,
    bernoulli_kl_divergence,
    evidence_weighted_binary_cross_entropy,
    prior_strength_to_pseudo_observation_count,
)


def test_smoothed_positive_target_is_minimized_at_point_nine_nine() -> None:
    probabilities = torch.tensor([0.90, 0.99, 0.999999])
    outcome = torch.ones(3)
    weight = torch.ones(3)
    losses = [
        evidence_weighted_binary_cross_entropy(probabilities[i : i + 1], outcome[i : i + 1], weight[i : i + 1]).item()
        for i in range(3)
    ]
    assert losses[1] < losses[0]
    assert losses[1] < losses[2]  # certainty is penalized relative to the 0.99 target


def test_zero_weight_excludes_observation() -> None:
    probabilities = torch.tensor([0.01, 0.99])
    outcome = torch.tensor([1.0, 1.0])
    weight = torch.tensor([0.0, 1.0])
    loss = evidence_weighted_binary_cross_entropy(probabilities, outcome, weight)
    only_second = evidence_weighted_binary_cross_entropy(probabilities[1:], outcome[1:], weight[1:])
    assert torch.allclose(loss, only_second)


def test_prior_strength_maps_to_pseudo_observations() -> None:
    assert abs(prior_strength_to_pseudo_observation_count(0.99) - 99.0) < 1e-9
    assert prior_strength_to_pseudo_observation_count(0.0) == 0.0
    with pytest.raises(ValueError):
        prior_strength_to_pseudo_observation_count(1.0)


def test_anchor_penalty_is_zero_at_prior_and_grows_away_from_it() -> None:
    prior = torch.full((2, 3), 0.99)
    pseudo_count = torch.full((2, 3), 99.0)
    assert anchored_link_prior_penalty(prior.clone(), prior, pseudo_count).item() < 1e-6
    moved = torch.full((2, 3), 0.5)
    assert anchored_link_prior_penalty(moved, prior, pseudo_count).item() > 1.0
    assert torch.all(bernoulli_kl_divergence(prior, moved) >= 0.0)


def test_log_space_loss_matches_the_clamped_loss_inside_the_clamp_and_keeps_the_gradient_beyond_it() -> None:
    logits = torch.tensor([-2.0, 0.5, 3.0, 15.0, -15.0], requires_grad=True)
    outcome = torch.tensor([1.0, 0.0, 1.0, 0.0, 1.0])
    weight = torch.ones(5)
    inside = slice(0, 3)
    clamped = evidence_weighted_binary_cross_entropy(torch.sigmoid(logits[inside]), outcome[inside], weight[inside])
    in_log_space = evidence_weighted_binary_cross_entropy(torch.sigmoid(logits[inside]), outcome[inside], weight[inside],
                                                          log_probability=torch.nn.functional.logsigmoid(logits[inside]),
                                                          log_complement=torch.nn.functional.logsigmoid(-logits[inside]))
    assert torch.allclose(clamped, in_log_space, atol=1e-6)
    # a negative at logit 15 and a positive at logit -15: no gradient through the clamp, the full one in log space
    clamped_gradient, = torch.autograd.grad(evidence_weighted_binary_cross_entropy(torch.sigmoid(logits), outcome, weight), logits)
    log_space_gradient, = torch.autograd.grad(evidence_weighted_binary_cross_entropy(
        torch.sigmoid(logits), outcome, weight, log_probability=torch.nn.functional.logsigmoid(logits), log_complement=torch.nn.functional.logsigmoid(-logits)), logits)
    assert clamped_gradient[3] == 0.0 and clamped_gradient[4] == 0.0
    assert log_space_gradient[3] > 0.15 and log_space_gradient[4] < -0.15  # about (sigmoid(z) - target) / 5


def test_equal_drug_loss_shares_gives_every_drug_the_same_row_total_and_leaves_genes_alone() -> None:
    import numpy as np

    from mechanistic_pathway_learning.models.soft_constraint_losses import equal_drug_loss_shares

    pair_weights = np.array([
        [1.0, 1.0, 1.0, 1.0],  # drug, total 4
        [0.5, 0.0, 0.0, 0.5],  # drug, total 1
        [2.0, 0.0, 1.0, 0.0],  # gene, total 3
        [0.0, 0.0, 0.0, 0.0],  # drug with every pair masked
        [3.0, 3.0, 0.0, 0.0],  # drug outside the fitted set, total 6
    ], dtype=np.float32)
    is_drug = np.array([True, True, False, True, True])
    fit_indices = np.array([0, 1, 2, 3])
    rescaled = equal_drug_loss_shares(pair_weights, is_drug, fit_indices)

    common_total = (4.0 + 1.0) / 2
    assert rescaled.dtype == pair_weights.dtype
    assert np.allclose(rescaled.sum(axis=1)[[0, 1, 4]], common_total)
    assert np.allclose(rescaled[2], pair_weights[2])
    assert np.allclose(rescaled[3], 0.0)
    assert np.isclose(rescaled[[0, 1]].sum(), pair_weights[[0, 1]].sum())
    # the pattern inside a row is kept: only its scale changes
    assert np.allclose(rescaled[1] / rescaled[1].sum(), pair_weights[1] / pair_weights[1].sum())


def test_equal_drug_loss_shares_without_fitted_drugs_returns_the_weights_unchanged() -> None:
    import numpy as np

    from mechanistic_pathway_learning.models.soft_constraint_losses import equal_drug_loss_shares

    pair_weights = np.array([[1.0, 2.0], [3.0, 0.0]])
    rescaled = equal_drug_loss_shares(pair_weights, np.array([False, True]), np.array([0]))
    assert np.array_equal(rescaled, pair_weights)
