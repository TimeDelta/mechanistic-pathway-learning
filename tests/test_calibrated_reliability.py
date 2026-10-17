"""The feature-dependent reliability regression calibrated on held-out labels (evidence/calibrated_reliability.py, decision 9.59)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.calibrated_reliability import (
    MINIMUM_CALIBRATION_CASES,
    calibrated_probability_matrix,
    cell_table,
    feature_frame,
    fit_calibration,
    permute_cell_rows,
    predict,
    species_of_sources,
    summarise_calibration,
)


def pairs_frame(rows):
    return pd.DataFrame(rows, columns=["perturbation_id", "symptom", "evidence_class", "direction", "n_distinct", "prior_probability", "strength", "sources"])


def test_species_is_read_from_the_sources():
    assert species_of_sources("IMPC-knockout") == "mouse" and species_of_sources("ClinicalTrials-arms;FAERS-disproportionality") == "human"
    assert species_of_sources("HPO-OMIM;IMPC-knockout") == "mouse" and species_of_sources("") == "human"


def test_the_feature_frame_one_hots_classes_sources_species_and_kind_and_logs_the_counts():
    pairs = pairs_frame([("G1", "anxiety", "animal", "negative", 1, 0.0, 0.3, "IMPC-knockout"),
                         ("D1", "anxiety", "human_association", "positive", 3, 0.7, 0.5, "ClinicalTrials-arms;FAERS-disproportionality")])
    frame = feature_frame(pairs, {"G1": "gene", "D1": "drug"}, {"D1": 2.5})
    assert list(frame.columns) == ["class:animal", "class:human_association", "direction:positive", "source:ClinicalTrials-arms", "source:FAERS-disproportionality",
                                   "source:IMPC-knockout", "log_records", "species:mouse", "perturbation:drug", "log_popularity"]
    assert frame.loc[0].tolist() == pytest.approx([1, 0, 0, 0, 0, 1, np.log1p(1), 1, 0, 0.0])
    assert frame.loc[1].tolist() == pytest.approx([0, 1, 1, 1, 1, 0, np.log1p(3), 0, 1, 2.5])


def test_cells_keep_only_pairs_inside_the_universe_and_rows_follow_a_label_permutation():
    pairs = pairs_frame([("G1", "anxiety", "animal", "negative", 1, 0.0, 0.3, "IMPC-knockout"),
                         ("G2", "insomnia", "animal", "positive", 2, 0.6, 0.3, "IMPC-knockout"),
                         ("G9", "anxiety", "animal", "positive", 1, 0.6, 0.3, "IMPC-knockout"),
                         ("G1", "unknown", "animal", "positive", 1, 0.6, 0.3, "IMPC-knockout")])
    frame = feature_frame(pairs, {}, {})
    cells = cell_table(pairs, frame, ["G1", "G2", "G3"], ["anxiety", "insomnia"])
    assert cells.row.tolist() == [0, 1] and cells.column.tolist() == [0, 1] and cells.direction.tolist() == ["negative", "positive"]
    moved = permute_cell_rows(cells, np.array([2, 0, 1]))  # position 0 now holds old row 2, position 1 old row 0, position 2 old row 1
    assert moved.row.tolist() == [1, 2]


def toy_calibration(n: int = 40, seed: int = 0):
    """Two sources on two symptoms: source A's cells are label positives three times in four, source B's one in four."""
    generator = np.random.default_rng(seed)
    rows = []
    perturbation_ids = [f"P{i}" for i in range(n)]
    for index, perturbation in enumerate(perturbation_ids):
        source = "A-source" if index % 2 == 0 else "B-source"
        rows.append((perturbation, "anxiety", "human_association", "positive", 1 + index % 3, 0.7, 0.5, source))
    pairs = pairs_frame(rows)
    frame = feature_frame(pairs, {p: ("drug" if i % 3 == 0 else "gene") for i, p in enumerate(perturbation_ids)}, {p: float(i % 5) for i, p in enumerate(perturbation_ids)})
    cells = cell_table(pairs, frame, perturbation_ids, ["anxiety", "insomnia"])
    outcomes = np.zeros((n, 2))
    for index in range(n):
        probability = 0.75 if index % 2 == 0 else 0.25
        outcomes[index, 0] = float(generator.random() < probability)
    return cells, outcomes, perturbation_ids


def test_the_calibration_keeps_the_tiers_below_the_minimum_and_fits_above_it():
    cells, outcomes, _ = toy_calibration()
    few = fit_calibration(cells, outcomes, None, np.arange(6))
    assert few.kept_tiers and "calibration cells" in few.reason and few.training_rows == 6
    tier_matrix = np.full(outcomes.shape, np.nan)
    tier_matrix[cells.row, cells.column] = 0.7
    assert calibrated_probability_matrix(few, cells, tier_matrix) is tier_matrix
    fitted = fit_calibration(cells, outcomes, None, np.arange(len(outcomes)), minimum_cases=5)
    assert not fitted.kept_tiers and fitted.training_rows == len(outcomes)
    assert fitted.coefficients["source:A-source"] > fitted.coefficients["source:B-source"], "the source whose cells are labels three times in four scores higher"
    predicted = predict(fitted, cells)
    a_cells, b_cells = cells["source:A-source"].to_numpy() == 1, cells["source:B-source"].to_numpy() == 1
    assert predicted[a_cells].mean() > 0.5 > predicted[b_cells].mean()
    matrix = calibrated_probability_matrix(fitted, cells, tier_matrix)
    assert np.isnan(matrix[0, 1]) and matrix[0, 0] == pytest.approx(predicted[0]) and not np.isnan(tier_matrix[0, 0])
    summary = summarise_calibration(fitted)
    assert set(summary["mean_prediction_at_median_popularity"]) == {"human_association:positive"} and summary["training_positives"] == int(outcomes[:, 0].sum())


def test_only_fit_rows_and_labelled_cells_train_the_regression():
    cells, outcomes, _ = toy_calibration()
    label_mask = np.ones(outcomes.shape, dtype=bool)
    label_mask[:10, 0] = False  # the first ten pairs were set aside by the selection
    fitted = fit_calibration(cells, outcomes, label_mask, np.arange(30), minimum_cases=3)
    assert fitted.training_rows == 20, "rows 10 to 29: in the fit and labelled"
    assert MINIMUM_CALIBRATION_CASES == 20


def test_the_prediction_is_read_at_the_median_popularity():
    cells, outcomes, _ = toy_calibration()
    fitted = fit_calibration(cells, outcomes, None, np.arange(len(outcomes)), minimum_cases=5)
    fitted.coefficients["log_popularity"] = 2.0  # force a popularity effect so the two readings differ
    at_median, raw = predict(fitted, cells, at_median_popularity=True), predict(fitted, cells, at_median_popularity=False)
    assert not np.allclose(at_median, raw)
    median_cells = cells.log_popularity.to_numpy() == fitted.median_log_popularity
    assert np.allclose(at_median[median_cells], raw[median_cells])


def test_the_statistic_enters_as_three_features_only_when_a_pair_carries_one():
    pairs = pairs_frame([("D1", "anxiety", "human_association", "positive", 2, 0.7, 0.5, "FAERS-disproportionality"),
                         ("D2", "anxiety", "human_association", "negative", 1, 0.0, 0.3, "ClinicalTrials-arms"),
                         ("G1", "anxiety", "animal", "negative", 1, 0.0, 0.3, "IMPC-knockout")])
    plain = feature_frame(pairs, {"D1": "drug", "D2": "drug", "G1": "gene"}, {})
    assert not [column for column in plain.columns if column.startswith("statistic:")], "tables without statistics fit the model they fitted before"
    with_statistic = pairs.assign(statistic_probability_of_excess=[0.9, 1.0, np.nan], statistic_information=[20.0, 3.0, np.nan])
    frame = feature_frame(with_statistic, {"D1": "drug", "D2": "drug", "G1": "gene"}, {})
    assert frame["statistic:logit_excess"].tolist() == pytest.approx([np.log(9.0) / 10.0, 1.0, 0.0]), "clipped at the logit bound and scaled to [-1, 1]; absent is 0"
    assert frame["statistic:logit_excess_squared"].tolist() == pytest.approx([(np.log(9.0) / 10.0) ** 2, 1.0, 0.0]), "a label rate need not rise monotonically"
    assert frame["statistic:log_information"].tolist() == pytest.approx([np.log1p(20.0), np.log1p(3.0), 0.0])
    assert frame["statistic:present"].tolist() == [1.0, 1.0, 0.0]
