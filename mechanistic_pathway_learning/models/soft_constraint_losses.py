"""Soft constraint losses (experiment design, section 5.4).

Nothing symptom-level is a hard constraint. Strength is expressed in two
places, both reversible under enough contradicting evidence:

1. Outcome level. Every observation (perturbation, symptom, relation, outcome)
   enters a weighted binary cross-entropy with its evidence-grade weight. Strong
   positives are smoothed to a target of 0.99 rather than 1.0, which is the
   "0.99 no-change" idea from the project discussion: the model is asked to keep
   the link near certain, not to make it certain.

2. Parameter level. Designated module-to-symptom links can be anchored to a prior
   probability with a Bernoulli KL penalty whose coefficient is a pseudo-observation
   count. A prior strength of 0.99 maps to strength / (1 - strength) = 99
   pseudo-observations, so moving the link away from its prior costs as much as
   contradicting 99 observations. Set the count to zero to disable the anchor.

Unobserved (perturbation, symptom) pairs are not true negatives. The training loop
is expected to pass them as degree-matched sampled negatives with a lower weight
(positive-unlabeled treatment); this module does not decide that policy.
"""
from __future__ import annotations

import numpy as np
import torch
from torch import Tensor

PROBABILITY_EPSILON = 1e-6


def equal_drug_loss_shares(pair_weights: np.ndarray, is_drug: np.ndarray, fit_indices: np.ndarray) -> np.ndarray:
    """Loss weights [perturbations, symptoms] with every drug's row rescaled to one common total: the mean row total
    of the drugs among fit_indices. Gene rows are unchanged, and so is the summed weight of the fitted drugs.

    The loss divides by the batch's summed weight, so a perturbation's share of the gradient follows its row total. A
    drug with many weighted positives (often one with many targets) would otherwise pull harder than a drug with one;
    with equal shares each drug counts the same. A drug row with no weight (every pair masked) stays at zero.
    """
    row_totals = pair_weights.sum(axis=1)
    weighted_drug = np.asarray(is_drug, dtype=bool) & (row_totals > 0)
    fitted_drugs = np.asarray(fit_indices, dtype=int)[weighted_drug[np.asarray(fit_indices, dtype=int)]]
    if not len(fitted_drugs):
        return pair_weights
    common_total = row_totals[fitted_drugs].mean()
    scale = np.ones(len(pair_weights), dtype=np.float64)
    scale[weighted_drug] = common_total / row_totals[weighted_drug]
    return (pair_weights * scale[:, None]).astype(pair_weights.dtype)


def evidence_weighted_binary_cross_entropy(
    symptom_probability: Tensor,
    observed_outcome: Tensor,
    evidence_weight: Tensor,
    positive_target: float = 0.99,
    negative_target: float = 0.0,
    log_probability: Tensor | None = None,
    log_complement: Tensor | None = None,
) -> Tensor:
    """Weighted BCE with label smoothing on positives; all tensors broadcast to the same shape.

    symptom_probability: model output in (0, 1).
    observed_outcome: 1.0 for an observed link, 0.0 for a sampled negative.
    evidence_weight: per-observation weight from the evidence grade (0 excludes the observation).
    log_probability, log_complement: log P and log(1 - P) computed by the head in log space. Given both, the loss uses
    them; otherwise it takes the logarithms of P clamped to [1e-6, 1 - 1e-6], where a pair beyond the clamp (a sigmoid
    logit below -13.8 for a positive, above 13.8 for a negative) has a constant loss and no gradient. PyTorch's
    BCEWithLogitsLoss documents the same point: combining the sigmoid and the loss "is more numerically stable than using
    a plain Sigmoid followed by a BCELoss".
    """
    smoothed_target = observed_outcome * positive_target + (1.0 - observed_outcome) * negative_target
    if log_probability is not None and log_complement is not None:
        per_observation_loss = -(smoothed_target * log_probability + (1.0 - smoothed_target) * log_complement)
        total_weight = evidence_weight.sum().clamp_min(PROBABILITY_EPSILON)
        return (per_observation_loss * evidence_weight).sum() / total_weight
    clamped_probability = symptom_probability.clamp(PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)
    per_observation_loss = -(
        smoothed_target * torch.log(clamped_probability) + (1.0 - smoothed_target) * torch.log1p(-clamped_probability)
    )
    total_weight = evidence_weight.sum().clamp_min(PROBABILITY_EPSILON)
    return (per_observation_loss * evidence_weight).sum() / total_weight


def bernoulli_kl_divergence(prior_probability: Tensor, model_probability: Tensor) -> Tensor:
    """KL(Bernoulli(prior) || Bernoulli(model)), elementwise."""
    prior = prior_probability.clamp(PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)
    model = model_probability.clamp(PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)
    return prior * (torch.log(prior) - torch.log(model)) + (1.0 - prior) * (torch.log1p(-prior) - torch.log1p(-model))


def prior_strength_to_pseudo_observation_count(prior_strength: float) -> float:
    """Map a strength in [0, 1) to the equivalent number of pseudo-observations: 0.99 -> 99."""
    if not 0.0 <= prior_strength < 1.0:
        raise ValueError("prior_strength must be in [0, 1); 1.0 would be a hard constraint, which this project does not use")
    return prior_strength / (1.0 - prior_strength)


def anchored_link_prior_penalty(
    link_probability: Tensor,
    prior_link_probability: Tensor,
    prior_pseudo_observation_count: Tensor,
) -> Tensor:
    """Sum over module-symptom links of pseudo_count * KL(prior || model).

    All three tensors have shape [num_pathway_modules, num_symptoms]. A zero pseudo-count
    disables the anchor for that link. Because the penalty is finite, enough contradicting
    data can still move the link; that is the difference from a hard constraint.
    """
    return (prior_pseudo_observation_count * bernoulli_kl_divergence(prior_link_probability, link_probability)).sum()


def attribution_hinge_penalty(attribution_entity_to_symptom: Tensor, literature_weight: Tensor, attribution_threshold_scale: float) -> Tensor:
    """Version 2 soft constraint (design section 5.4): max(0, scale * weight - attribution), summed.

    attribution_entity_to_symptom: attribution mass from a literature-named entity to a symptom,
    one value per literature triple. literature_weight: the triple's evidence weight. The loop
    that computes attributions (for example integrated gradients on the node-state field) is
    not implemented in version 1; the penalty is kept here so the config can reference it.
    """
    return torch.relu(attribution_threshold_scale * literature_weight - attribution_entity_to_symptom).sum()
