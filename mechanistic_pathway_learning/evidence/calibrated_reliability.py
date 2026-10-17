"""A feature-dependent reliability model calibrated on held-out labels (decision 9.59, the owner, 16 October 2026).

What it estimates. For a soft-evidence pair (evidence/soft_priors.py: one row per perturbation, symptom, class and
direction, with its sources and its deduplicated record count), the probability that the pair is a label-grade
statement, read off the pairs the fit trains on: a logistic regression of the label outcome (1 where the training
fold's label table holds the pair as a positive, 0 where it does not) on the row's features. Its prediction replaces
the fixed prior probability of configs/soft_prior_tiers.csv; the tier's unit strength and the tempered count still
set how hard the prior pulls, so the regression says what to believe and the tiers say how much.

Features. The evidence class, the direction, each source present, the logarithm of the record count, the species the
evidence was measured in (a property of the source, SPECIES_OF_SOURCE, mouse for IMPC and human otherwise; a primate
or other source takes its own level when one enters), the perturbation type, and the perturbation's literature
popularity (log of its distinct PubMed references in the literature table, 0 without one). Species is a feature and not
a separate fit, as the owner asked, so a non-human absence is weighted by what the regression learns about mouse rows
from the pairs that human sources also cover.

The comparison's own statistic (decision 9.48, option 3: "Carrying stat to calibrated proper is the way"). Where a
pair's evidence carries one (evidence/soft_priors.py keeps the most informative comparison's), four more features:
statistic:logit_excess, the log-odds of its posterior probability that the drug raises the event rate, clipped to
+-LOGIT_CLIP and divided by it, so it lies in [-1, 1]; its square, statistic:logit_excess_squared, because the label
rate need not rise monotonically with the statistic (on the development pairs of v5, 17 October 2026, FAERS pairs
whose posterior put the excess below 0.025 were label positives more often than pairs with an inconclusive posterior:
docs/open_decisions.md, 9.48); statistic:log_information, the logarithm of one plus its information; and
statistic:present. A pair without a statistic takes 0 in all four, so the indicator carries the difference. The
regression then says what each level of the statistic is worth on the training-fold labels, which is the calibration
the decision asked for. When no pair carries a statistic the columns are not made, so a run on tables built before
the statistics fits the model it fitted before.

Popularity. It is a feature so that the volume of reports on a well-studied perturbation is attributed to reporting
propensity rather than to truth (the SAR reading of positive-only sources, docs/soft_priors.md); the prediction is
read with popularity set to its median over the training rows, so the prior carries the report-type effect and not
the popularity effect. The raw prediction is kept beside it in the summary for the reader.

Where it fits. On the perturbations of the rows a fit trains on, never on the lockbox and never on a held-out fold;
the trainer passes those rows (--soft-prior-probability calibrated). Below MINIMUM_CALIBRATION_CASES positives or
negatives among them the tiers are kept and the summary says why.

What it is not. Not Dawid-Skene: there is no latent truth and no conditional-independence assumption; the target is
the label table, with its own ascertainment, which is why the popularity term exists and why the arm is read against
the permuted prior (docs/soft_priors.md, section 3).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

SPECIES_OF_SOURCE = {"IMPC-knockout": "mouse"}
LOGIT_CLIP = 10.0
DEFAULT_SPECIES = "human"
MINIMUM_CALIBRATION_CASES = 20
REGULARISATION_C = 1.0


@dataclass
class Calibration:
    features: list[str]
    coefficients: dict[str, float]
    intercept: float
    training_rows: int
    training_positives: int
    median_log_popularity: float
    predicted_by_class_and_direction: dict[str, float]
    raw_predicted_by_class_and_direction: dict[str, float]
    kept_tiers: bool
    reason: str = ""


def species_of_sources(sources: str) -> str:
    """The species of a pair's evidence: the first non-human species among its sources, else human."""
    for source in str(sources).split(";"):
        species = SPECIES_OF_SOURCE.get(source.strip())
        if species:
            return species
    return DEFAULT_SPECIES


def feature_frame(pairs: pd.DataFrame, perturbation_types: dict[str, str], log_popularity: dict[str, float]) -> pd.DataFrame:
    """One row per resolved pair, with the columns named in the module docstring; categorical levels are one-hot
    columns named level:value, every source its own column source:name."""
    frame = pd.DataFrame(index=pairs.index)
    for level in sorted(set(pairs.evidence_class)):
        frame[f"class:{level}"] = (pairs.evidence_class == level).astype(float)
    frame["direction:positive"] = (pairs.direction == "positive").astype(float)
    sources = sorted({source for cell in pairs.sources for source in str(cell).split(";") if source})
    for source in sources:
        frame[f"source:{source}"] = pairs.sources.map(lambda cell: float(source in str(cell).split(";")))
    frame["log_records"] = np.log1p(pairs.n_distinct.astype(float))
    species = pairs.sources.map(species_of_sources)
    for level in sorted(set(species) - {DEFAULT_SPECIES}):
        frame[f"species:{level}"] = (species == level).astype(float)
    kinds = pairs.perturbation_id.map(lambda identifier: perturbation_types.get(identifier, ""))
    frame["perturbation:drug"] = (kinds == "drug").astype(float)
    frame["log_popularity"] = pairs.perturbation_id.map(lambda identifier: float(log_popularity.get(identifier, 0.0)))
    if "statistic_probability_of_excess" in pairs.columns and pd.to_numeric(pairs.statistic_probability_of_excess, errors="coerce").notna().any():
        probability = pd.to_numeric(pairs.statistic_probability_of_excess, errors="coerce")
        information = pd.to_numeric(pairs.get("statistic_information", pd.Series(np.nan, index=pairs.index)), errors="coerce")
        present = probability.notna() & information.notna()
        clipped = probability.clip(1e-12, 1.0 - 1e-12)
        scaled = np.clip(np.log(clipped / (1.0 - clipped)), -LOGIT_CLIP, LOGIT_CLIP) / LOGIT_CLIP
        frame["statistic:logit_excess"] = np.where(present, scaled, 0.0)
        frame["statistic:logit_excess_squared"] = np.where(present, scaled ** 2, 0.0)
        frame["statistic:log_information"] = np.where(present, np.log1p(information.clip(lower=0.0)), 0.0)
        frame["statistic:present"] = present.astype(float)
    return frame


def cell_table(pairs: pd.DataFrame, features: pd.DataFrame, perturbation_ids: list[str], symptoms: list[str]) -> pd.DataFrame:
    """The resolved pairs that fall inside a data set's perturbation and symptom universe, with their row and column and
    their features; the rows are what the calibration fits on and what it predicts for."""
    row_of = {identifier: index for index, identifier in enumerate(perturbation_ids)}
    column_of = {symptom: index for index, symptom in enumerate(symptoms)}
    rows = pairs.perturbation_id.map(row_of)
    columns = pairs.symptom.map(column_of)
    inside = rows.notna() & columns.notna()
    table = features[inside].copy()
    table.insert(0, "column", columns[inside].astype(int).to_numpy())
    table.insert(0, "row", rows[inside].astype(int).to_numpy())
    table["evidence_class"] = pairs.evidence_class[inside].to_numpy()
    table["direction"] = pairs.direction[inside].to_numpy()
    return table.reset_index(drop=True)


def permute_cell_rows(cells: pd.DataFrame, source_row: np.ndarray) -> pd.DataFrame:
    """After --permute-labels moved row source_row[i] to position i, a prior cell on the old row r sits at the position
    whose source is r."""
    position_of_source = np.argsort(np.asarray(source_row))
    moved = cells.copy()
    moved["row"] = position_of_source[moved.row.to_numpy()]
    return moved


FEATURE_COLUMNS_EXCLUDED = ("row", "column", "evidence_class", "direction")


def fit_calibration(cells: pd.DataFrame, outcomes: np.ndarray, label_mask: np.ndarray | None, fit_rows: np.ndarray,
                    minimum_cases: int = MINIMUM_CALIBRATION_CASES, regularisation_c: float = REGULARISATION_C) -> Calibration:
    """The regression on the prior cells whose row the fit trains on and whose label the table states: the target is 1
    where the cell is a label positive, 0 where it is labelled and not positive (a pair the selection set aside is
    not a statement and is left out). Returns a Calibration that keeps the tiers when either side has fewer than
    minimum_cases cells."""
    feature_names = [column for column in cells.columns if column not in FEATURE_COLUMNS_EXCLUDED]
    in_fit = np.zeros(outcomes.shape[0], dtype=bool)
    in_fit[np.asarray(fit_rows, dtype=int)] = True
    labelled = np.ones(outcomes.shape, dtype=bool) if label_mask is None else np.asarray(label_mask, dtype=bool)
    rows, columns = cells.row.to_numpy(), cells.column.to_numpy()
    usable = in_fit[rows] & labelled[rows, columns]
    targets = (outcomes[rows, columns] > 0).astype(int)
    training = cells[usable]
    training_targets = targets[usable]
    positives, negatives = int(training_targets.sum()), int((1 - training_targets).sum())
    median_log_popularity = float(np.median(training.log_popularity)) if len(training) else 0.0
    if positives < minimum_cases or negatives < minimum_cases:
        return Calibration(feature_names, {}, 0.0, int(len(training)), positives, median_log_popularity, {}, {}, True,
                           f"{positives} positive and {negatives} negative calibration cells on the fit rows; at least {minimum_cases} of each are needed")
    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(C=regularisation_c, max_iter=2000)
    model.fit(training[feature_names].to_numpy(dtype=float), training_targets)
    coefficients = {name: float(value) for name, value in zip(feature_names, model.coef_[0])}
    calibration = Calibration(feature_names, coefficients, float(model.intercept_[0]), int(len(training)), positives, median_log_popularity, {}, {}, False)
    predicted = predict(calibration, cells, at_median_popularity=True)
    raw = predict(calibration, cells, at_median_popularity=False)
    groups = cells.evidence_class.astype(str) + ":" + cells.direction.astype(str)
    calibration.predicted_by_class_and_direction = {str(group): float(np.mean(predicted[groups == group])) for group in sorted(set(groups))}
    calibration.raw_predicted_by_class_and_direction = {str(group): float(np.mean(raw[groups == group])) for group in sorted(set(groups))}
    return calibration


def predict(calibration: Calibration, cells: pd.DataFrame, at_median_popularity: bool = True) -> np.ndarray:
    """The regression's probability per cell; with at_median_popularity the popularity feature is held at the
    training median so the prior carries the report-type effect and not the popularity effect."""
    features = cells[calibration.features].to_numpy(dtype=float)
    if at_median_popularity and "log_popularity" in calibration.features:
        features = features.copy()
        features[:, calibration.features.index("log_popularity")] = calibration.median_log_popularity
    logits = calibration.intercept + features @ np.array([calibration.coefficients[name] for name in calibration.features])
    return 1.0 / (1.0 + np.exp(-np.clip(logits, -30.0, 30.0)))


def calibrated_probability_matrix(calibration: Calibration, cells: pd.DataFrame, tier_probability: np.ndarray) -> np.ndarray:
    """The prior-probability matrix with the regression's prediction written on every prior cell; the tiers' matrix is
    returned unchanged when the calibration kept the tiers."""
    if calibration.kept_tiers:
        return tier_probability
    probability = np.array(tier_probability, copy=True)
    probability[cells.row.to_numpy(), cells.column.to_numpy()] = predict(calibration, cells, at_median_popularity=True)
    return probability


def summarise_calibration(calibration: Calibration) -> dict:
    return {"kept_tiers": calibration.kept_tiers, "reason": calibration.reason, "training_cells": calibration.training_rows,
            "training_positives": calibration.training_positives, "median_log_popularity": calibration.median_log_popularity,
            "intercept": calibration.intercept, "coefficients": calibration.coefficients,
            "mean_prediction_at_median_popularity": calibration.predicted_by_class_and_direction,
            "mean_prediction_raw": calibration.raw_predicted_by_class_and_direction}

