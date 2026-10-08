"""Strata of binding measurement coverage, and the macro AUPRC read inside them."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.measurement_coverage_strata import (
    GENE_PERTURBATION_STRATUM,
    UNMEASURED_DRUG_STRATUM,
    measurement_coverage_strata,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import macro_auprc_by_degree_bin, macro_auprc_by_stratum


def coverage_table(rows: list[tuple[str, str, int]]) -> pd.DataFrame:
    table = pd.DataFrame(rows, columns=["perturbation_id", "perturbation_type", "genes_with_an_affinity"])
    return table.set_index("perturbation_id")


def test_genes_and_unmeasured_drugs_are_their_own_strata() -> None:
    coverage = coverage_table([("gene:HMBS", "gene", -1), ("drug:1", "drug", -1), ("drug:2", "drug", 3), ("drug:3", "drug", 40)])
    stratum_of_row, labels = measurement_coverage_strata(["gene:HMBS", "drug:1", "drug:2", "drug:3"], coverage, num_bins=2)
    assert stratum_of_row[0] == GENE_PERTURBATION_STRATUM
    assert stratum_of_row[1] == UNMEASURED_DRUG_STRATUM
    assert stratum_of_row[2] != stratum_of_row[3]  # the thinly and the heavily measured drug are separate
    assert labels[:2] == [GENE_PERTURBATION_STRATUM, UNMEASURED_DRUG_STRATUM]
    assert labels[2:] == [stratum_of_row[2], stratum_of_row[3]]  # bins follow, in ascending order of coverage
    assert "n1" in stratum_of_row[2] and "[3,3]" in stratum_of_row[2]


def test_bins_hold_equal_counts_and_carry_their_range() -> None:
    measured = [(f"drug:{index}", "drug", index) for index in range(9)]
    stratum_of_row, labels = measurement_coverage_strata([row[0] for row in measured], coverage_table(measured), num_bins=3)
    assert len(labels) == 3
    assert [stratum_of_row.count(label) for label in labels] == [3, 3, 3]
    assert labels[0].endswith("_n3") and "[0,2]" in labels[0] and "[6,8]" in labels[2]


def test_a_perturbation_missing_from_the_table_is_not_a_measured_drug() -> None:
    coverage = coverage_table([("drug:1", "drug", 5)])
    stratum_of_row, labels = measurement_coverage_strata(["drug:1", "drug:absent"], coverage, num_bins=2)
    assert stratum_of_row[1] == GENE_PERTURBATION_STRATUM  # no row means no drug measurement, so never a measured stratum
    assert sum(label.startswith("measured_genes_bin") for label in labels) == 1


def test_macro_auprc_by_stratum_scores_each_stratum_apart() -> None:
    # eight perturbations, one symptom: the model ranks the first four perfectly and the last four backwards
    outcomes = np.array([[1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0], [1.0]])
    outcomes[1::2] = 0.0
    predictions = np.array([[0.9], [0.1], [0.8], [0.2], [0.1], [0.9], [0.2], [0.8]])
    strata = ["well_measured"] * 4 + ["thinly_measured"] * 4
    scores = macro_auprc_by_stratum(predictions, outcomes, strata, minimum_positives=2)
    assert list(scores) == ["well_measured", "thinly_measured"]  # first-seen order
    assert scores["well_measured"] == 1.0
    assert scores["thinly_measured"] < 0.6
    assert np.isnan(macro_auprc_by_stratum(predictions, outcomes, strata, minimum_positives=5)["well_measured"])


def test_degree_bins_still_read_as_before() -> None:
    # each tercile holds one positive and one negative of the one scored symptom, ranked correctly, so each bin reads 1.0
    outcomes = np.array([[1.0], [0.0], [1.0], [0.0], [1.0], [0.0]])
    predictions = np.array([[0.9], [0.1], [0.9], [0.1], [0.9], [0.1]])
    degrees = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    bins = macro_auprc_by_degree_bin(predictions, outcomes, degrees, num_bins=3, minimum_positives=1)
    assert list(bins) == ["degree_bin_0_[1,2]_n2", "degree_bin_1_[3,4]_n2", "degree_bin_2_[5,6]_n2"]
    assert all(value == 1.0 for value in bins.values())
    assert np.isnan(macro_auprc_by_degree_bin(predictions, outcomes, degrees, num_bins=3, minimum_positives=2)["degree_bin_0_[1,2]_n2"])
