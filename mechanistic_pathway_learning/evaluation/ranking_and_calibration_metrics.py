"""Metrics (design section 6.2): per-symptom AUROC and AUPRC with bootstrap intervals
over perturbations, mean reciprocal rank and hits-at-k for symptom ranking per
perturbation, and expected calibration error. Degree-bin and grade stratification
is done by the caller passing the relevant row subset.

Every metric takes an optional label mask (the same shape as outcomes, True where a pair is labelled). A pair masked
out by a label selection (experiment_data.load_experiment_data, label_selection) is neither a positive nor a negative:
it is left out of the per-symptom rankings, of the per-perturbation symptom rankings and of calibration. With no mask
every pair is scored, as before.
"""
from __future__ import annotations

from dataclasses import dataclass

import warnings

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


@dataclass
class IntervalEstimate:
    point: float
    lower: float
    upper: float


def bootstrap_interval(statistic, predictions: np.ndarray, outcomes: np.ndarray, num_bootstrap: int = 1000, random_seed: int = 0, confidence: float = 0.95,
                       mask: np.ndarray | None = None) -> IntervalEstimate:
    """Percentile bootstrap over rows (perturbations) of a statistic(predictions, outcomes), or of
    statistic(predictions, outcomes, mask) when a label mask is given (resampled with the same rows)."""
    generator = np.random.default_rng(random_seed)
    num_rows = predictions.shape[0]

    def evaluate(rows):
        if mask is None:
            return statistic(predictions[rows], outcomes[rows])
        return statistic(predictions[rows], outcomes[rows], mask[rows])

    point = evaluate(np.arange(num_rows))
    resampled_values = []
    for _ in range(num_bootstrap):
        row_indices = generator.integers(0, num_rows, size=num_rows)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                resampled_value = evaluate(row_indices)
        except ValueError:  # a resample with a single class, older scikit-learn
            continue
        if np.isnan(resampled_value):  # a resample with a single class, newer scikit-learn
            continue
        resampled_values.append(resampled_value)
    if not resampled_values:
        return IntervalEstimate(point, float("nan"), float("nan"))
    lower_quantile = (1.0 - confidence) / 2.0
    return IntervalEstimate(point, float(np.quantile(resampled_values, lower_quantile)), float(np.quantile(resampled_values, 1.0 - lower_quantile)))


def scored_rows(num_rows: int, symptom_index: int, mask: np.ndarray | None) -> np.ndarray:
    """Boolean rows that carry a label for one symptom: every row without a mask."""
    return np.ones(num_rows, dtype=bool) if mask is None else np.asarray(mask[:, symptom_index], dtype=bool)


def per_symptom_auroc(predictions: np.ndarray, outcomes: np.ndarray, symptom_index: int, mask: np.ndarray | None = None) -> float:
    rows = scored_rows(len(outcomes), symptom_index, mask)
    return float(roc_auc_score(outcomes[rows, symptom_index], predictions[rows, symptom_index]))


def per_symptom_auprc(predictions: np.ndarray, outcomes: np.ndarray, symptom_index: int, mask: np.ndarray | None = None) -> float:
    rows = scored_rows(len(outcomes), symptom_index, mask)
    return float(average_precision_score(outcomes[rows, symptom_index], predictions[rows, symptom_index]))


def scorable_symptom(outcomes: np.ndarray, symptom_index: int, minimum_positives: int, mask: np.ndarray | None = None) -> bool:
    """At least minimum_positives positives and at least one negative among the labelled rows of a symptom."""
    rows = scored_rows(len(outcomes), symptom_index, mask)
    positives = outcomes[rows, symptom_index].sum()
    return bool(positives >= minimum_positives and positives < rows.sum())


def mean_reciprocal_rank(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None = None) -> float:
    """For each perturbation with at least one positive symptom, 1 / rank of the best-ranked positive (masked
    symptoms of that perturbation are left out of its ranking)."""
    reciprocal_ranks = []
    for row_index, (row_predictions, row_outcomes) in enumerate(zip(predictions, outcomes)):
        if mask is not None:
            row_predictions, row_outcomes = row_predictions[mask[row_index]], row_outcomes[mask[row_index]]
        if row_outcomes.sum() == 0:
            continue
        ranking = np.argsort(-row_predictions)
        positive_positions = np.where(row_outcomes[ranking] > 0)[0]
        reciprocal_ranks.append(1.0 / (positive_positions[0] + 1))
    return float(np.mean(reciprocal_ranks)) if reciprocal_ranks else float("nan")


def hits_at_k(predictions: np.ndarray, outcomes: np.ndarray, k: int = 3, mask: np.ndarray | None = None) -> float:
    hits = []
    for row_index, (row_predictions, row_outcomes) in enumerate(zip(predictions, outcomes)):
        if mask is not None:
            row_predictions, row_outcomes = row_predictions[mask[row_index]], row_outcomes[mask[row_index]]
        if row_outcomes.sum() == 0:
            continue
        top_k = np.argsort(-row_predictions)[:k]
        hits.append(float(row_outcomes[top_k].sum() > 0))
    return float(np.mean(hits)) if hits else float("nan")


def expected_calibration_error(predictions: np.ndarray, outcomes: np.ndarray, num_bins: int = 10, mask: np.ndarray | None = None) -> float:
    flat_predictions = predictions.ravel() if mask is None else predictions[np.asarray(mask, dtype=bool)]
    flat_outcomes = outcomes.ravel() if mask is None else outcomes[np.asarray(mask, dtype=bool)]
    bin_edges = np.linspace(0.0, 1.0, num_bins + 1)
    error = 0.0
    for lower, upper in zip(bin_edges[:-1], bin_edges[1:]):
        in_bin = (flat_predictions >= lower) & (flat_predictions < upper if upper < 1.0 else flat_predictions <= upper)
        if in_bin.sum() == 0:
            continue
        error += in_bin.mean() * abs(flat_predictions[in_bin].mean() - flat_outcomes[in_bin].mean())
    return float(error)


def macro_auprc_by_degree_bin(predictions: np.ndarray, outcomes: np.ndarray, perturbation_degrees: np.ndarray, num_bins: int = 3, minimum_positives: int = 5,
                              mask: np.ndarray | None = None) -> dict[str, float]:
    """Macro AUPRC over symptoms inside each degree bin of the perturbations (design section 6.2, assumption A9).

    Bins are degree quantiles over the rows given; a symptom is scored inside a bin when it has at least
    minimum_positives positives and one negative there. Hub-driven predictors score well in the top bin only.
    """
    order = np.argsort(perturbation_degrees, kind="stable")  # equal-count bins; ties broken by position so no bin is empty
    bin_of_row = np.empty(len(perturbation_degrees), dtype=int)
    bin_of_row[order] = np.minimum(np.arange(len(perturbation_degrees)) * num_bins // max(1, len(perturbation_degrees)), num_bins - 1)
    result: dict[str, float] = {}
    for bin_index in range(num_bins):
        rows = np.where(bin_of_row == bin_index)[0]
        edges = [perturbation_degrees[rows].min() if len(rows) else float("nan"), perturbation_degrees[rows].max() if len(rows) else float("nan")]
        values = []
        for symptom_index in range(outcomes.shape[1]):
            labelled = rows[scored_rows(len(outcomes), symptom_index, mask)[rows]]
            positives = outcomes[labelled, symptom_index].sum()
            if positives < minimum_positives or positives == len(labelled):
                continue
            values.append(float(average_precision_score(outcomes[labelled, symptom_index], predictions[labelled, symptom_index])))
        label = f"degree_bin_{bin_index}_[{edges[0]:.0f},{edges[1]:.0f}]_n{len(rows)}"
        result[label] = float(np.mean(values)) if values else float("nan")
    return result


def micro_auprc(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Average precision over every labelled (perturbation, symptom) pair ranked together in one list.

    Every symptom enters, however few positives it has, so the symptoms the macro average leaves out (fewer than five
    positives in a test set) still count. Pairs of frequent symptoms weigh more, and a predictor earns credit for ranking
    symptoms by their base rates, which is why it is read beside the macro average and against the same baselines
    (popularity among them) rather than alone. NaN when the labelled pairs hold no positive or no negative."""
    labelled = np.ones(outcomes.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    labels = outcomes[labelled]
    if labels.size == 0 or labels.sum() == 0 or labels.sum() == labels.size:
        return float("nan")
    return float(average_precision_score(labels, predictions[labelled]))


def paired_bootstrap_difference(predictions_a: np.ndarray, predictions_b: np.ndarray, outcomes: np.ndarray, aggregate, num_bootstrap: int = 1000,
                                random_seed: int = 0, confidence: float = 0.95, mask: np.ndarray | None = None) -> dict[str, float]:
    """Paired bootstrap over perturbations of aggregate(A) - aggregate(B) on the same rows (design section 7), for any
    aggregate(predictions, outcomes, mask) such as the macro average of paired_bootstrap_macro_difference or micro_auprc.

    Both prediction matrices are resampled with the same row indices, so the interval is for the paired difference;
    resamples where the difference is undefined are skipped.
    """
    generator = np.random.default_rng(random_seed)
    point = aggregate(predictions_a, outcomes, mask) - aggregate(predictions_b, outcomes, mask)
    differences = []
    for _ in range(num_bootstrap):
        rows = generator.integers(0, outcomes.shape[0], size=outcomes.shape[0])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resampled_mask = None if mask is None else mask[rows]
            difference = aggregate(predictions_a[rows], outcomes[rows], resampled_mask) - aggregate(predictions_b[rows], outcomes[rows], resampled_mask)
        if not np.isnan(difference):
            differences.append(difference)
    lower_quantile = (1.0 - confidence) / 2.0
    return {"difference": float(point), "lower": float(np.quantile(differences, lower_quantile)) if differences else float("nan"),
            "upper": float(np.quantile(differences, 1.0 - lower_quantile)) if differences else float("nan"),
            "fraction_resamples_favoring_a": float(np.mean(np.array(differences) > 0)) if differences else float("nan"), "num_resamples": len(differences)}


def macro_statistic(statistic=per_symptom_auprc, minimum_positives: int = 5):
    """aggregate(predictions, outcomes, mask): the mean of statistic over the symptoms with at least minimum_positives
    positives and one negative among the labelled rows; NaN when none qualifies."""
    def macro(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None) -> float:
        values = []
        for symptom_index in range(outcomes.shape[1]):
            if not scorable_symptom(outcomes, symptom_index, minimum_positives, mask):
                continue
            if mask is None:
                values.append(statistic(predictions, outcomes, symptom_index))
            else:
                values.append(statistic(predictions, outcomes, symptom_index, mask))
        return float(np.mean(values)) if values else float("nan")
    return macro


def paired_bootstrap_macro_difference(predictions_a: np.ndarray, predictions_b: np.ndarray, outcomes: np.ndarray, statistic=per_symptom_auprc, num_bootstrap: int = 1000,
                                      random_seed: int = 0, confidence: float = 0.95, minimum_positives: int = 5, mask: np.ndarray | None = None) -> dict[str, float]:
    """Paired bootstrap over perturbations of macro(statistic of A) - macro(statistic of B) on the same rows (design section 7);
    symptoms with fewer than minimum_positives positives or no negatives in a resample are skipped."""
    return paired_bootstrap_difference(predictions_a, predictions_b, outcomes, macro_statistic(statistic, minimum_positives), num_bootstrap, random_seed, confidence, mask)


def paired_bootstrap_micro_difference(predictions_a: np.ndarray, predictions_b: np.ndarray, outcomes: np.ndarray, num_bootstrap: int = 1000,
                                      random_seed: int = 0, confidence: float = 0.95, mask: np.ndarray | None = None) -> dict[str, float]:
    """Paired bootstrap over perturbations of micro_auprc(A) - micro_auprc(B) on the same rows."""
    return paired_bootstrap_difference(predictions_a, predictions_b, outcomes, micro_auprc, num_bootstrap, random_seed, confidence, mask)


def degree_strata(perturbation_degrees: np.ndarray, num_strata: int = 5) -> np.ndarray:
    """Stratum index per perturbation from quantile edges of its degree, with equal degrees always in one stratum.

    Degree is coarse and heavily tied (on the metabolic slice 30 percent of perturbations have degree 1), so the
    edges are deduplicated and fewer than num_strata strata can result. Used to score rankings within degree strata,
    which gives a model no credit for ordering perturbations by degree (design section 6.2).
    """
    edges = np.unique(np.quantile(perturbation_degrees, np.linspace(0.0, 1.0, num_strata + 1)[1:-1]))
    return np.searchsorted(edges, perturbation_degrees, side="right")


def rank_normalise_within_groups(predictions: np.ndarray, group_of_row: np.ndarray) -> np.ndarray:
    """Replace each score by its tie-averaged rank divided by (group size + 1), per symptom column, inside each group.

    Pooling raw scores across hold-outs whose base rates differ lets a predictor whose score scale tracks the
    training base rate look worse or better than it is (design section 6.2; review v0.4, finding 1). Ranks
    inside each hold-out remove the scale; the (n + 1) denominator keeps a constant predictor at the same
    value in every hold-out, which rank / n does not. Rows with group -1 are left unchanged.
    """
    from scipy.stats import rankdata

    normalised = np.array(predictions, dtype=float, copy=True)
    for group in np.unique(group_of_row):
        if group < 0:
            continue
        rows = np.flatnonzero(group_of_row == group)
        if len(rows) == 0:
            continue
        for column in range(predictions.shape[1]):
            normalised[rows, column] = rankdata(predictions[rows, column], method="average") / (len(rows) + 1.0)
    return normalised


def stratified_auroc(predictions: np.ndarray, outcomes: np.ndarray, group_of_row: np.ndarray, symptom_index: int, mask: np.ndarray | None = None) -> float:
    """AUROC from positive-negative pairs formed inside each group only (a within-hold-out Mann-Whitney statistic).

    Pairs that cross hold-outs never enter, so differences of base rate or score scale between hold-outs cannot
    move it. Returns NaN when no group has both a positive and a negative.
    """
    concordant = 0.0
    pairs = 0.0
    for group in np.unique(group_of_row):
        if group < 0:
            continue
        rows = np.flatnonzero((group_of_row == group) & scored_rows(len(outcomes), symptom_index, mask))
        scores = predictions[rows, symptom_index]
        labels = outcomes[rows, symptom_index] > 0.5
        positive_scores, negative_scores = scores[labels], scores[~labels]
        if len(positive_scores) == 0 or len(negative_scores) == 0:
            continue
        comparison = positive_scores[:, None] - negative_scores[None, :]
        concordant += float((comparison > 0).sum() + 0.5 * (comparison == 0).sum())
        pairs += float(comparison.size)
    return concordant / pairs if pairs else float("nan")
