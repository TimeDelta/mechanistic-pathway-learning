"""Learned evidence reliability (experiment design, section 4.3, learned weighting).

Fixed grade weights and count-based weights share a flaw: they are set by hand.
This module learns how much to trust each evidence source from the pattern of
agreement across sources, treating every source as a noisy sensor of the true
perturbation-symptom link (Dawid and Skene 1979; Raykar et al. 2010 for the
feature-dependent form).

Model. Item i is a (perturbation, symptom, relation) triple with a latent truth
z_i in {0, 1} and prior P(z_i = 1) = prevalence. Source s reports r_{i,s} in {0, 1}
for every item it covers (a source covers an item when it has any record for that
item's perturbation, so silence inside coverage is a 0 report, and items outside
coverage carry no report). Each source has a sensitivity P(r = 1 | z = 1) and a
specificity P(r = 0 | z = 0). With rubric features (for example from the LLM
appraisal), sensitivity and specificity become logistic functions of the features,
so a well-designed study in a source counts for more than a poorly designed one in
the same source.

Why this is learnable when a bare loss weight is not: the parameters are inside the
likelihood of the reports, not multipliers on a penalty, so the EM objective rewards
explaining agreement rather than discarding evidence. Beta priors on sensitivity and
specificity keep sources from degenerating when they are few.

Output. The posterior P(z_i = 1 | reports) is the soft label and loss weight for the
training observation; a prior strength of 0.99 for a link is no longer hand-set but
read off the posterior and scaled by a single tuned hyperparameter.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

PROBABILITY_EPSILON = 1e-6


@dataclass
class ReliabilityParameters:
    source_names: list[str]
    prevalence: float
    sensitivity: np.ndarray  # [num_sources]
    specificity: np.ndarray  # [num_sources]
    log_likelihood_trace: list[float] = field(default_factory=list)


def _clip(values: np.ndarray) -> np.ndarray:
    return np.clip(values, PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)


def posterior_truth_probability(
    reports: np.ndarray,
    coverage: np.ndarray,
    prevalence: float,
    sensitivity: np.ndarray,
    specificity: np.ndarray,
) -> np.ndarray:
    """E-step: P(z_i = 1 | reports), shape [num_items].

    reports: [num_items, num_sources] in {0, 1}; coverage: same shape, 1 where the source covers the item.
    """
    sensitivity = _clip(sensitivity)[None, :]
    specificity = _clip(specificity)[None, :]
    log_likelihood_true = coverage * (reports * np.log(sensitivity) + (1 - reports) * np.log(1 - sensitivity))
    log_likelihood_false = coverage * (reports * np.log(1 - specificity) + (1 - reports) * np.log(specificity))
    log_odds = np.log(prevalence) - np.log(1 - prevalence) + (log_likelihood_true - log_likelihood_false).sum(axis=1)
    return 1.0 / (1.0 + np.exp(-log_odds))


def fit_dawid_skene(
    reports: np.ndarray,
    coverage: np.ndarray,
    source_names: list[str],
    num_iterations: int = 100,
    sensitivity_prior_counts: tuple[float, float] = (8.0, 2.0),
    specificity_prior_counts: tuple[float, float] = (8.0, 2.0),
    prevalence_prior_counts: tuple[float, float] = (1.0, 1.0),
    convergence_tolerance: float = 1e-6,
) -> tuple[ReliabilityParameters, np.ndarray]:
    """EM for per-source sensitivity and specificity with Beta priors (MAP M-step).

    Returns the fitted parameters and the posterior truth probability per item.
    Initialization is the majority vote over covering sources (the usual Dawid-Skene
    start), which avoids the symmetric local optimum that a flat start can reach.
    A two-class latent model needs at least three conditionally independent sources
    to be identified; with fewer, the priors carry the estimate.
    """
    reports = np.asarray(reports, dtype=float)
    coverage = np.asarray(coverage, dtype=float)
    num_items, num_sources = reports.shape
    sensitivity = np.full(num_sources, 0.8)
    specificity = np.full(num_sources, 0.8)
    prevalence = 0.5
    trace: list[float] = []
    covering_counts = coverage.sum(axis=1)
    posterior = np.where(covering_counts > 0, (reports * coverage).sum(axis=1) / np.maximum(covering_counts, 1.0), 0.5)
    posterior = _clip(posterior)
    for _ in range(num_iterations):
        # M-step with Beta priors
        prevalence = (posterior.sum() + prevalence_prior_counts[0]) / (num_items + sum(prevalence_prior_counts))
        covered_true = (posterior[:, None] * coverage)
        covered_false = ((1.0 - posterior)[:, None] * coverage)
        sensitivity = (covered_true * reports).sum(axis=0) + sensitivity_prior_counts[0] - 1.0
        sensitivity = sensitivity / (covered_true.sum(axis=0) + sum(sensitivity_prior_counts) - 2.0)
        specificity = (covered_false * (1 - reports)).sum(axis=0) + specificity_prior_counts[0] - 1.0
        specificity = specificity / (covered_false.sum(axis=0) + sum(specificity_prior_counts) - 2.0)
        sensitivity, specificity = _clip(sensitivity), _clip(specificity)
        # E-step
        new_posterior = posterior_truth_probability(reports, coverage, prevalence, sensitivity, specificity)
        log_likelihood = marginal_log_likelihood(reports, coverage, prevalence, sensitivity, specificity)
        trace.append(float(log_likelihood))
        converged = np.max(np.abs(new_posterior - posterior)) < convergence_tolerance
        posterior = new_posterior
        if converged:
            break
    return ReliabilityParameters(source_names, float(prevalence), sensitivity, specificity, trace), posterior


def marginal_log_likelihood(reports: np.ndarray, coverage: np.ndarray, prevalence: float, sensitivity: np.ndarray, specificity: np.ndarray) -> float:
    sensitivity = _clip(sensitivity)[None, :]
    specificity = _clip(specificity)[None, :]
    log_true = np.log(prevalence) + (coverage * (reports * np.log(sensitivity) + (1 - reports) * np.log(1 - sensitivity))).sum(axis=1)
    log_false = np.log(1 - prevalence) + (coverage * (reports * np.log(1 - specificity) + (1 - reports) * np.log(specificity))).sum(axis=1)
    stacked = np.stack([log_true, log_false], axis=1)
    maximum = stacked.max(axis=1, keepdims=True)
    return float((maximum[:, 0] + np.log(np.exp(stacked - maximum).sum(axis=1))).sum())


def fit_feature_dependent_sensitivity(
    report_features: np.ndarray,
    reports: np.ndarray,
    posterior: np.ndarray,
    regularization_strength: float = 1.0,
):
    """Raykar-style M-step for one source: logistic sensitivity and specificity as functions of rubric features.

    report_features: [num_covered_items, num_features] for the items this source covers;
    reports: [num_covered_items]; posterior: [num_covered_items] truth probabilities.
    Returns two fitted scikit-learn classifiers (sensitivity model fit on items weighted by
    posterior, specificity model fit on items weighted by 1 - posterior) or None for a side
    with no effective weight. The caller alternates this with posterior_truth_probability
    using per-item sensitivity and specificity instead of per-source constants.
    """
    from sklearn.linear_model import LogisticRegression

    def fit_weighted(targets: np.ndarray, weights: np.ndarray):
        if weights.sum() < PROBABILITY_EPSILON or len(np.unique(targets[weights > PROBABILITY_EPSILON])) < 2:
            return None
        classifier = LogisticRegression(C=1.0 / regularization_strength, max_iter=1000)
        classifier.fit(report_features, targets, sample_weight=weights)
        return classifier

    sensitivity_model = fit_weighted(reports.astype(int), posterior)
    specificity_model = fit_weighted((1 - reports).astype(int), 1.0 - posterior)
    return sensitivity_model, specificity_model


def per_item_reliability(sensitivity_model, specificity_model, report_features: np.ndarray, fallback_sensitivity: float, fallback_specificity: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-item sensitivity and specificity from the fitted feature models, with constant fallbacks."""
    num_items = report_features.shape[0]
    sensitivity = np.full(num_items, fallback_sensitivity) if sensitivity_model is None else sensitivity_model.predict_proba(report_features)[:, 1]
    specificity = np.full(num_items, fallback_specificity) if specificity_model is None else specificity_model.predict_proba(report_features)[:, 1]
    return _clip(sensitivity), _clip(specificity)


def observation_weights_from_posterior(posterior: np.ndarray, global_scale: float = 1.0) -> np.ndarray:
    """Loss weight per observation: the posterior truth probability times one tuned global scale.

    The global scale is the only hyperparameter left for validation-set tuning; everything
    else about evidence strength is learned from the reports.
    """
    return global_scale * posterior


def anchor_pseudo_counts_from_posterior(posterior: np.ndarray, maximum_pseudo_count: float = 99.0) -> np.ndarray:
    """Pseudo-observation count for a link anchor derived from the posterior.

    Posterior 0.99 maps to 99 pseudo-observations (0.99 / 0.01), capped at maximum_pseudo_count,
    so the near-hard anchor discussed for strong evidence is read off the data rather than set.
    """
    clipped = np.clip(posterior, 0.0, 1.0 - PROBABILITY_EPSILON)
    return np.minimum(clipped / (1.0 - clipped), maximum_pseudo_count)
