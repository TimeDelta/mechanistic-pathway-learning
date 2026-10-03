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

Weighted reports. The report table (evidence_reports.py) gives several reports per
(item, source) cell, each with a rubric weight in (0, 1] and a value of 1 (presence) or
0 (absence). The cell then enters the likelihood through two weighted counts,
positive_weights[i, s] and negative_weights[i, s]: sensitivity^pos (1 - sensitivity)^neg
under z = 1 and (1 - specificity)^pos specificity^neg under z = 0. A source's silence on
an item whose perturbation it covers is one implicit negative report of weight
implicit_negative_weight. The {0, 1} model above is the special case pos = reports x
coverage, neg = (1 - reports) x coverage, and fit_dawid_skene stays as that wrapper.

One report's worth per source record. Several reports of one source record on one item
(a disease entry annotated to three descendant terms of the same symptom, several labels
of one drug mentioning the same term) are not independent observations, so inside a cell
the weights are first reduced per source_record_id to the largest rubric weight among the
record's positive reports and, separately, among its negative reports, and only then summed
over distinct records. A cell therefore carries at most one report's weight per record and
side, and a pair listed under three cognitive terms of one synopsis weighs the same as a
pair listed once.

Specificity needs explicit negatives. The M-step for specificity divides the posterior-
weighted negative count by the posterior-weighted total; for a source with no explicit
negative report the only negatives are silence, and the estimate becomes a decreasing
function of how much the source reports rather than a measurement of anything. Such a
source keeps its specificity at the prior mean (estimate_specificity False for it), which
the fit records in specificity_held_at_prior.

Prevalence per item group. No source bridges gene items and drug items, so one pooled
prevalence would let the agreement of the two HPO provenances set the prior of every drug
pair. item_groups (the evidence class in the assembler) gives each group its own
prevalence; prevalence_by_group holds them and prevalence the pooled value.

Degeneracy. Clipping keeps sensitivity and specificity inside (0, 1) but not their sum
above 1, and EM can settle on the label-switched optimum for one source, where a positive
report lowers the posterior (likelihood ratio sensitivity / (1 - specificity) below 1) and
silence raises it. A source whose fitted sensitivity + specificity is at most 1 is listed
in degenerate_sources; the assembler refuses --weighting reliability on such a fit.

Identifiability. A two-class latent model needs three conditionally independent
sources covering the same items. The report table has two HPO provenances on gene
items and one SIDER source per relation on drug items, so no item is covered by three
sources: ReportReliabilityFit.weakly_identified is True on every such table and the
posterior is written as a column but is not the default weight (design section 12,
open question 12).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

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


@dataclass
class ReportReliabilityFit:
    """Fitted sensor model over weighted reports, with the posterior per item in table order."""

    source_names: list[str]
    prevalence: float
    sensitivity: np.ndarray  # [num_sources]
    specificity: np.ndarray  # [num_sources]
    item_keys: list[tuple[str, str, str]]  # (perturbation_id, symptom, relation)
    posterior: np.ndarray  # [num_items] P(z_i = 1 | reports)
    num_iterations: int
    converged: bool
    weakly_identified: bool  # True when no item is covered by three or more sources
    log_likelihood_trace: list[float] = field(default_factory=list)
    per_source_summary: dict[str, dict] = field(default_factory=dict)
    prevalence_by_group: dict[str, float] = field(default_factory=dict)  # one prevalence per item group (evidence class); prevalence is the pooled value
    specificity_held_at_prior: np.ndarray | None = None  # [num_sources] True where specificity was not estimated (no explicit negative report)
    degenerate_sources: list[str] = field(default_factory=list)  # sources with fitted sensitivity + specificity <= 1


def report_count_matrices(
    reports: pd.DataFrame,
    item_keys: list[tuple[str, str, str]],
    source_names: list[str],
    implicit_negative_weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """positive_weights, negative_weights, coverage, each [num_items, num_sources], from a report table with rubric_weight filled.

    positive_weights[i, s] holds, summed over the distinct source records (source_record_id) of source s on
    item i, the largest rubric weight among that record's reports with report_value 1; negative_weights[i, s]
    the same over reports with report_value 0. A record annotated to several descendant terms of one symptom,
    or a drug with several labels mentioning one term, therefore contributes at most one report's weight per
    side (module docstring). A table without a source_record_id column treats every report as its own
    record. coverage[i, s] is 1 when source s has any report for the item's perturbation_id and relation (for
    any symptom); a covered cell with no explicit report gets implicit_negative_weight in negative_weights.
    An item is never created by this function: items are the keys passed in, and a source outside its
    coverage contributes nothing.
    """
    item_index = {key: index for index, key in enumerate(item_keys)}
    source_index = {name: index for index, name in enumerate(source_names)}
    positive_weights = np.zeros((len(item_keys), len(source_names)))
    negative_weights = np.zeros_like(positive_weights)
    explicit = np.zeros_like(positive_weights, dtype=bool)
    covered_perturbations: dict[tuple[str, str], set[str]] = {}
    if len(reports):
        weights = reports.rubric_weight.to_numpy(dtype=float)
        values = reports.report_value.to_numpy(dtype=int)
        record_ids = reports.source_record_id.astype(str).tolist() if "source_record_id" in reports.columns else [str(row_number) for row_number in range(len(reports))]
        largest_weight_per_record: dict[tuple[int, int, str, int], float] = {}
        for row_number, (perturbation_id, symptom, relation, source) in enumerate(zip(reports.perturbation_id, reports.symptom, reports.relation, reports.source)):
            column = source_index.get(source)
            if column is None:
                continue
            covered_perturbations.setdefault((perturbation_id, relation), set()).add(source)
            row = item_index.get((perturbation_id, symptom, relation))
            if row is None:
                continue
            explicit[row, column] = True
            record_key = (row, column, record_ids[row_number], int(values[row_number] == 1))
            largest_weight_per_record[record_key] = max(largest_weight_per_record.get(record_key, 0.0), weights[row_number])
        for (row, column, _, is_positive), weight in largest_weight_per_record.items():
            if is_positive:
                positive_weights[row, column] += weight
            else:
                negative_weights[row, column] += weight
    coverage = np.zeros_like(positive_weights)
    for row, (perturbation_id, _, relation) in enumerate(item_keys):
        for source in covered_perturbations.get((perturbation_id, relation), ()):
            coverage[row, source_index[source]] = 1.0
    negative_weights += implicit_negative_weight * coverage * (~explicit)
    return positive_weights, negative_weights, coverage


def weighted_log_likelihood_terms(positive_weights: np.ndarray, negative_weights: np.ndarray, sensitivity: np.ndarray, specificity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-item log likelihood of the weighted counts under z = 1 and under z = 0, summed over sources."""
    sensitivity = _clip(np.asarray(sensitivity, dtype=float))[None, :]
    specificity = _clip(np.asarray(specificity, dtype=float))[None, :]
    log_true = (positive_weights * np.log(sensitivity) + negative_weights * np.log(1 - sensitivity)).sum(axis=1)
    log_false = (positive_weights * np.log(1 - specificity) + negative_weights * np.log(specificity)).sum(axis=1)
    return log_true, log_false


def weighted_posterior_truth_probability(positive_weights: np.ndarray, negative_weights: np.ndarray, prevalence: float | np.ndarray, sensitivity: np.ndarray, specificity: np.ndarray) -> np.ndarray:
    """E-step over weighted counts: P(z_i = 1 | reports), shape [num_items]; prevalence is one value or one per item."""
    log_true, log_false = weighted_log_likelihood_terms(positive_weights, negative_weights, sensitivity, specificity)
    prevalence = np.asarray(prevalence, dtype=float)
    log_odds = np.log(prevalence) - np.log(1 - prevalence) + log_true - log_false
    return 1.0 / (1.0 + np.exp(-log_odds))


def weighted_marginal_log_likelihood(positive_weights: np.ndarray, negative_weights: np.ndarray, prevalence: float | np.ndarray, sensitivity: np.ndarray, specificity: np.ndarray) -> float:
    log_true, log_false = weighted_log_likelihood_terms(positive_weights, negative_weights, sensitivity, specificity)
    prevalence = np.asarray(prevalence, dtype=float)
    stacked = np.stack([np.log(prevalence) + log_true, np.log(1 - prevalence) + log_false], axis=1)
    maximum = stacked.max(axis=1, keepdims=True)
    return float((maximum[:, 0] + np.log(np.exp(stacked - maximum).sum(axis=1))).sum())


def fit_weighted_dawid_skene(
    positive_weights: np.ndarray,
    negative_weights: np.ndarray,
    source_names: list[str],
    item_keys: list[tuple[str, str, str]] | None = None,
    num_iterations: int = 100,
    sensitivity_prior_counts: tuple[float, float] = (8.0, 2.0),
    specificity_prior_counts: tuple[float, float] = (8.0, 2.0),
    prevalence_prior_counts: tuple[float, float] = (1.0, 1.0),
    convergence_tolerance: float = 1e-6,
    item_groups: list[str] | None = None,
    estimate_specificity: np.ndarray | None = None,
) -> ReportReliabilityFit:
    """EM for per-source sensitivity and specificity over weighted positive and negative counts, Beta priors (MAP M-step).

    Initialization is the weighted majority vote (positive weight over total weight per item, 0.5 for an
    item with no weight); the M-step divides posterior-weighted counts, sensitivity_s = (sum_i post_i pos[i, s]
    + a - 1) / (sum_i post_i (pos + neg)[i, s] + a + b - 2) and its mirror for specificity; the E-step is
    weighted_posterior_truth_probability; convergence is a maximum absolute posterior change below the
    tolerance. weakly_identified is True when no item is covered (any weight) by three or more sources.

    item_groups, one label per item, gives each group its own prevalence (the E-step reads the item's group);
    None pools every item. estimate_specificity, one boolean per source, holds the specificity of a source at
    the prior mean a / (a + b) where False (a source without explicit negative reports, module docstring);
    None estimates every source. A source whose fitted sensitivity + specificity is at most 1 is listed in
    degenerate_sources.
    """
    positive_weights = np.asarray(positive_weights, dtype=float)
    negative_weights = np.asarray(negative_weights, dtype=float)
    num_items, num_sources = positive_weights.shape
    keys = list(item_keys) if item_keys is not None else [(str(index), "", "") for index in range(num_items)]
    group_labels = [str(label) for label in item_groups] if item_groups is not None else ["all"] * num_items
    if len(group_labels) != num_items:
        raise ValueError(f"item_groups has {len(group_labels)} labels for {num_items} items")
    distinct_groups = sorted(set(group_labels))
    group_index = np.array([distinct_groups.index(label) for label in group_labels], dtype=int)
    specificity_estimated = np.ones(num_sources, dtype=bool) if estimate_specificity is None else np.asarray(estimate_specificity, dtype=bool)
    if specificity_estimated.shape != (num_sources,):
        raise ValueError(f"estimate_specificity has shape {specificity_estimated.shape} for {num_sources} sources")
    specificity_prior_mean = specificity_prior_counts[0] / sum(specificity_prior_counts)
    total_weights = positive_weights + negative_weights
    covering_sources_per_item = (total_weights > 0).sum(axis=1)
    weakly_identified = bool(num_items == 0 or covering_sources_per_item.max() < 3)
    sensitivity = np.full(num_sources, 0.8)
    specificity = np.full(num_sources, 0.8)
    prevalence = 0.5
    prevalence_per_group = np.full(len(distinct_groups), 0.5)
    trace: list[float] = []
    item_total = total_weights.sum(axis=1)
    posterior = np.where(item_total > 0, positive_weights.sum(axis=1) / np.maximum(item_total, PROBABILITY_EPSILON), 0.5)
    posterior = _clip(posterior)
    converged = False
    iterations_run = 0
    for iterations_run in range(1, num_iterations + 1):
        prevalence = (posterior.sum() + prevalence_prior_counts[0]) / (num_items + sum(prevalence_prior_counts))
        group_posterior_sums = np.bincount(group_index, weights=posterior, minlength=len(distinct_groups))
        group_sizes = np.bincount(group_index, minlength=len(distinct_groups))
        prevalence_per_group = (group_posterior_sums + prevalence_prior_counts[0]) / (group_sizes + sum(prevalence_prior_counts))
        sensitivity = (posterior[:, None] * positive_weights).sum(axis=0) + sensitivity_prior_counts[0] - 1.0
        sensitivity = sensitivity / ((posterior[:, None] * total_weights).sum(axis=0) + sum(sensitivity_prior_counts) - 2.0)
        specificity = ((1.0 - posterior)[:, None] * negative_weights).sum(axis=0) + specificity_prior_counts[0] - 1.0
        specificity = specificity / (((1.0 - posterior)[:, None] * total_weights).sum(axis=0) + sum(specificity_prior_counts) - 2.0)
        specificity = np.where(specificity_estimated, specificity, specificity_prior_mean)
        sensitivity, specificity = _clip(sensitivity), _clip(specificity)
        prevalence_per_item = prevalence_per_group[group_index]
        new_posterior = weighted_posterior_truth_probability(positive_weights, negative_weights, prevalence_per_item, sensitivity, specificity)
        trace.append(weighted_marginal_log_likelihood(positive_weights, negative_weights, prevalence_per_item, sensitivity, specificity))
        converged = bool(np.max(np.abs(new_posterior - posterior)) < convergence_tolerance) if num_items else True
        posterior = new_posterior
        if converged:
            break
    degenerate_sources = [name for name, source_sensitivity, source_specificity in zip(source_names, sensitivity, specificity) if source_sensitivity + source_specificity <= 1.0]
    return ReportReliabilityFit(
        source_names=list(source_names), prevalence=float(prevalence), sensitivity=sensitivity, specificity=specificity, item_keys=keys, posterior=posterior,
        num_iterations=iterations_run, converged=converged, weakly_identified=weakly_identified, log_likelihood_trace=trace,
        prevalence_by_group={label: float(value) for label, value in zip(distinct_groups, prevalence_per_group)},
        specificity_held_at_prior=~specificity_estimated, degenerate_sources=degenerate_sources,
    )


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
    """EM for per-source sensitivity and specificity with Beta priors (MAP M-step) over {0, 1} reports.

    Returns the fitted parameters and the posterior truth probability per item. This is
    fit_weighted_dawid_skene with pos = reports x coverage and neg = (1 - reports) x coverage.
    Initialization is the majority vote over covering sources (the usual Dawid-Skene
    start), which avoids the symmetric local optimum that a flat start can reach.
    A two-class latent model needs at least three conditionally independent sources
    to be identified; with fewer, the priors carry the estimate.
    """
    reports = np.asarray(reports, dtype=float)
    coverage = np.asarray(coverage, dtype=float)
    fit = fit_weighted_dawid_skene(reports * coverage, (1.0 - reports) * coverage, source_names, None, num_iterations,
                                   sensitivity_prior_counts, specificity_prior_counts, prevalence_prior_counts, convergence_tolerance)
    return ReliabilityParameters(fit.source_names, fit.prevalence, fit.sensitivity, fit.specificity, fit.log_likelihood_trace), fit.posterior


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
