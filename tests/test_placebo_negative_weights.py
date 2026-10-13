"""Amendment item 3: the never-above-placebo negatives take a per-pair weight where the drug's label states the two
arm sizes, and the registered 0.5 where it does not."""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from run_main_model import (
    REGISTERED_PLACEBO_NEGATIVE_WEIGHTS,
    loss_pair_weights,
    placebo_negative_weight_matrix,
)

REGISTERED_DEFAULT_PLACEBO_WEIGHT = 0.5
REGISTERED_MEDIAN_EFFECTIVE_SIZE_OVER_DRUGS = 214.39  # amendment item 3, corrected 13 October 2026 to the constant the file was built with


def weights_file(directory: Path, rows: str) -> Path:
    path = directory / "placebo_negative_weights.csv"
    path.write_text("perturbation_id,symptom,weight\n" + rows)
    return path


def test_the_file_fills_only_the_pairs_it_names(tmp_path: Path) -> None:
    matrix, summary = placebo_negative_weight_matrix(
        weights_file(tmp_path, "D1,insomnia,0.62\nD2,fatigue,0.5\n"), ["D1", "D2"], ["insomnia", "fatigue"])
    assert matrix[0, 0] == pytest.approx(0.62) and matrix[1, 1] == pytest.approx(0.5)
    assert np.isnan(matrix[0, 1]) and np.isnan(matrix[1, 0])
    assert summary["pairs_matched"] == 2 and summary["pairs_not_in_this_evidence_table"] == []


def test_a_pair_the_evidence_table_no_longer_holds_is_counted_rather_than_applied(tmp_path: Path) -> None:
    matrix, summary = placebo_negative_weight_matrix(
        weights_file(tmp_path, "D1,insomnia,0.62\nDROPPED,fatigue,0.7\nD1,NOT_A_SYMPTOM,0.7\n"), ["D1"], ["insomnia"])
    assert summary["pairs_matched"] == 1
    assert summary["pairs_not_in_this_evidence_table"] == ["DROPPED/fatigue", "D1/NOT_A_SYMPTOM"]
    assert matrix.shape == (1, 1)


def test_the_wrong_file_is_refused_and_so_is_a_weight_below_zero(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="wrong file"):
        placebo_negative_weight_matrix(weights_file(tmp_path, "OTHER,insomnia,0.62\n"), ["D1"], ["insomnia"])
    with pytest.raises(ValueError, match="below zero"):
        placebo_negative_weight_matrix(weights_file(tmp_path, "D1,insomnia,-0.5\n"), ["D1"], ["insomnia"])
    no_weight_column = tmp_path / "other.csv"
    no_weight_column.write_text("perturbation_id,symptom\nD1,insomnia\n")
    with pytest.raises(ValueError, match="not a placebo negative weight file"):
        placebo_negative_weight_matrix(no_weight_column, ["D1"], ["insomnia"])


def experiment_data_with_one_evidenced_negative_per_row(per_pair: np.ndarray | None) -> SimpleNamespace:
    outcomes = np.zeros((2, 2))
    outcomes[0, 0] = 1.0  # a positive, which the per-pair weight must not touch
    evidenced = np.array([[True, True], [False, True]])
    return SimpleNamespace(outcomes=outcomes, weights=np.full_like(outcomes, 0.8), evidenced_negative=evidenced,
                           perturbation_types=["drug", "drug"], placebo_negative_weight_matrix=per_pair)


def test_the_per_pair_weight_overrides_the_default_for_an_evidenced_negative_only() -> None:
    arguments = SimpleNamespace(negative_weight=0.2, placebo_negative_weight=REGISTERED_DEFAULT_PLACEBO_WEIGHT, equal_drug_shares=False)
    per_pair = np.array([[0.9, 0.62], [0.9, np.nan]])  # [0,0] is a positive and [1,0] is not an evidenced negative
    weights = loss_pair_weights(experiment_data_with_one_evidenced_negative_per_row(per_pair), None, arguments, np.arange(2))
    assert weights[0, 0] == pytest.approx(0.8)                                 # positive, at its evidence weight
    assert weights[0, 1] == pytest.approx(0.62)                                # evidenced negative the file names
    assert weights[1, 0] == pytest.approx(0.2)                                 # an ordinary negative
    assert weights[1, 1] == pytest.approx(REGISTERED_DEFAULT_PLACEBO_WEIGHT)   # evidenced negative the file leaves out


def test_without_a_file_every_evidenced_negative_keeps_the_registered_default() -> None:
    arguments = SimpleNamespace(negative_weight=0.2, placebo_negative_weight=REGISTERED_DEFAULT_PLACEBO_WEIGHT, equal_drug_shares=False)
    weights = loss_pair_weights(experiment_data_with_one_evidenced_negative_per_row(None), None, arguments, np.arange(2))
    assert weights[0, 1] == pytest.approx(REGISTERED_DEFAULT_PLACEBO_WEIGHT) and weights[1, 1] == pytest.approx(REGISTERED_DEFAULT_PLACEBO_WEIGHT)


def test_the_registered_file_holds_the_counts_and_the_formula_the_amendment_states() -> None:
    """The file and item 3 have to agree, which they did not until 13 October 2026: the item named the median over the
    pairs (196.56) and the file was built with the median over the drugs (214.39). The item now states the drug
    median, and this test holds the two together."""
    import pandas as pd

    table = pd.read_csv(REGISTERED_PLACEBO_NEGATIVE_WEIGHTS)
    measured = table[table.effective_size.notna()]
    assert len(table) == 61 and len(measured) == 41
    assert (table[table.effective_size.isna()].weight == REGISTERED_DEFAULT_PLACEBO_WEIGHT).all()
    effective_size = 1.0 / (1.0 / measured.drug_arm + 1.0 / measured.placebo_arm)
    np.testing.assert_allclose(measured.effective_size, effective_size)

    median_over_drugs = measured.drop_duplicates("drug").effective_size.median()
    assert measured.drug.nunique() == 25
    assert median_over_drugs == pytest.approx(REGISTERED_MEDIAN_EFFECTIVE_SIZE_OVER_DRUGS, abs=0.005)
    np.testing.assert_allclose(measured.weight, effective_size / (effective_size + median_over_drugs))
    assert measured.weight.min() == pytest.approx(0.057, abs=0.0005) and measured.weight.max() == pytest.approx(0.772, abs=0.0005)
