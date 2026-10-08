"""Strata of how much binding measurement each drug has, for the overstudied-drug reading (docs/off_target_scoping.md).

A drug that was run through a binding panel has measured affinities at dozens of human targets; a drug measured only at
the target its label names has one or two. That is a property of who screened it, not of its pharmacology, so a reading
that uses measured off-targets, counts them or occupies them, can rise with how well a drug was studied rather than
with any mechanism. The cure is not a correction term but a stratified reading: report the metric inside strata of
measurement coverage and look at whether the ordering holds in the thinly measured stratum too.

experiments/scope_off_target_binding.py writes the coverage table (one row per perturbation, --coverage-table). This
module turns it into one stratum label per scored row.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

GENE_PERTURBATION_STRATUM = "gene_perturbation"
UNMEASURED_DRUG_STRATUM = "drug_without_binding_measurements"
NOT_IN_THE_TABLE_STRATUM = "not_in_the_coverage_table"
COVERAGE_COLUMN = "genes_with_an_affinity"


def read_measurement_coverage(path: Path) -> pd.DataFrame:
    """The coverage table, indexed by perturbation id."""
    table = pd.read_parquet(path)
    missing = {"perturbation_id", "perturbation_type", COVERAGE_COLUMN} - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing {sorted(missing)}; regenerate it with experiments/scope_off_target_binding.py")
    return table.set_index("perturbation_id")


def measurement_coverage_strata(perturbation_ids: list[str], coverage: pd.DataFrame, num_bins: int = 3) -> tuple[list[str], list[str]]:
    """One stratum label per perturbation, and the labels in reading order.

    A perturbation that is not a drug, or a drug with no row in either binding source, is its own stratum: neither has a
    measured off-target, so pooling them with the measured drugs would hide the comparison rather than make it. The
    measured drugs are cut into num_bins equal-count bins by how many human genes carry an affinity for them, ties
    broken by position so no bin is empty. A perturbation the table does not hold at all (a coverage table built on
    another evidence table) is a stratum of its own, never silently a gene.
    """
    rows = coverage.reindex(perturbation_ids)
    in_table = rows.perturbation_type.notna().to_numpy(dtype=bool)
    is_drug = (rows.perturbation_type == "drug").to_numpy(dtype=bool)
    measured_genes = rows[COVERAGE_COLUMN].to_numpy(dtype=float)
    is_measured_drug = is_drug & np.isfinite(measured_genes) & (measured_genes >= 0)
    measured_positions = np.where(is_measured_drug)[0]
    stratum_of_row = [NOT_IN_THE_TABLE_STRATUM if not known else GENE_PERTURBATION_STRATUM if not drug else UNMEASURED_DRUG_STRATUM
                      for known, drug in zip(in_table, is_drug)]
    bin_labels: list[str] = []
    if len(measured_positions):
        order = measured_positions[np.argsort(measured_genes[measured_positions], kind="stable")]
        bin_of_measured = np.minimum(np.arange(len(order)) * num_bins // len(order), num_bins - 1)
        for bin_index in range(num_bins):
            positions = order[bin_of_measured == bin_index]
            if not len(positions):
                continue
            counts = measured_genes[positions]
            label = f"measured_genes_bin_{bin_index}_[{counts.min():.0f},{counts.max():.0f}]_n{len(positions)}"
            bin_labels.append(label)
            for position in positions:
                stratum_of_row[position] = label
    present = set(stratum_of_row)
    labels = [label for label in (GENE_PERTURBATION_STRATUM, UNMEASURED_DRUG_STRATUM, NOT_IN_THE_TABLE_STRATUM) if label in present] + bin_labels
    return stratum_of_row, labels
