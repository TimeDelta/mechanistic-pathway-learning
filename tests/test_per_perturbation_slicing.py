"""restrict_to_perturbations must slice every per-perturbation array.

A field left out of PER_PERTURBATION_ARRAYS keeps its full length while outcomes is cut to the development rows,
and the first loss that multiplies the two raises a broadcast error. evidenced_negative and
placebo_negative_weight_matrix were both missing until 13 October 2026, which broke every development run using a
placebo-negative label selection. This test is written against the dataclass rather than against a list of names,
so a field added later is covered without anyone remembering to add it here.
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from mechanistic_pathway_learning.evaluation.experiment_data import (
    PER_PERTURBATION_ARRAYS,
    ExperimentData,
    restrict_to_perturbations,
)

NUM_PERTURBATIONS = 4
NUM_SYMPTOMS = 3


def experiment_data_with_every_array_filled() -> ExperimentData:
    fields = {field.name: field for field in dataclasses.fields(ExperimentData)}
    values = {
        "node_ids": ["N0", "N1"], "node_index": {"N0": 0, "N1": 1},
        "edge_source": np.array([0]), "edge_target": np.array([1]), "edge_relation": np.array([0]),
        "relation_types": ["r"],
        "perturbation_ids": [f"P{i}" for i in range(NUM_PERTURBATIONS)],
        "perturbation_types": ["gene"] * NUM_PERTURBATIONS,
        "perturbation_labels": [f"L{i}" for i in range(NUM_PERTURBATIONS)],
        "perturbation_seeds": [[0] for _ in range(NUM_PERTURBATIONS)],
        "group_ids": [f"g{i}" for i in range(NUM_PERTURBATIONS)],
        "symptoms": [f"s{i}" for i in range(NUM_SYMPTOMS)],
        "node_types": ["gene", "gene"], "is_currency": np.array([False, False]), "node_degree": np.array([1, 1]),
        "perturbation_signs": [[-1.0] for _ in range(NUM_PERTURBATIONS)],
        "perturbation_magnitudes": [[1.0] for _ in range(NUM_PERTURBATIONS)],
        "outcomes": np.arange(NUM_PERTURBATIONS * NUM_SYMPTOMS, dtype=float).reshape(NUM_PERTURBATIONS, NUM_SYMPTOMS),
        "weights": np.ones((NUM_PERTURBATIONS, NUM_SYMPTOMS)),
    }
    # Filled by name, NOT from PER_PERTURBATION_ARRAYS: a fixture built from the constant under test would leave a
    # missing field unset, the field would be None, and the check below would skip the very omission it is for.
    for name in ("in_metabolic_layer", "evidence_dates", "evidence_date_is_publication"):
        values[name] = np.ones(NUM_PERTURBATIONS, dtype=bool if name != "evidence_dates" else float)
    for name in ("label_mask", "evidenced_negative", "masked_for_soft_evidence"):
        values[name] = np.ones((NUM_PERTURBATIONS, NUM_SYMPTOMS), dtype=bool)
    values["frequencies"] = np.ones((NUM_PERTURBATIONS, NUM_SYMPTOMS))
    values["placebo_negative_weight_matrix"] = np.full((NUM_PERTURBATIONS, NUM_SYMPTOMS), 0.5)
    unknown = set(PER_PERTURBATION_ARRAYS) - set(values)
    assert not unknown, f"PER_PERTURBATION_ARRAYS names {unknown}, which this fixture does not fill; add them here by name"
    missing = [name for name, field in fields.items()
               if name not in values and field.default is dataclasses.MISSING and field.default_factory is dataclasses.MISSING]
    assert not missing, f"the fixture does not supply required fields {missing}"
    return ExperimentData(**{name: value for name, value in values.items() if name in fields})


def test_every_array_whose_first_axis_is_the_perturbations_is_sliced():
    data = experiment_data_with_every_array_filled()
    keep = np.array([True, False, True, False])
    restricted = restrict_to_perturbations(data, keep)
    assert len(restricted.perturbation_ids) == 2
    for field in dataclasses.fields(ExperimentData):
        value = getattr(data, field.name)
        if not isinstance(value, np.ndarray) or value.shape[:1] != (NUM_PERTURBATIONS,):
            continue
        sliced = getattr(restricted, field.name)
        assert sliced is not None and sliced.shape[0] == 2, (
            f"{field.name} is a per-perturbation array that restrict_to_perturbations did not slice; "
            f"add it to PER_PERTURBATION_ARRAYS")


def test_the_loss_weights_still_broadcast_after_restricting():
    """The failure the omission produced: a full-length evidenced_negative against development-length outcomes."""
    import sys
    from pathlib import Path
    from types import SimpleNamespace

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
    from run_main_model import loss_pair_weights

    data = experiment_data_with_every_array_filled()
    restricted = restrict_to_perturbations(data, np.array([True, False, True, False]))
    arguments = SimpleNamespace(negative_weight=0.2, placebo_negative_weight=0.5, equal_drug_shares=False)
    weights = loss_pair_weights(restricted, restricted.label_mask, arguments, np.arange(2))
    assert weights.shape == (2, NUM_SYMPTOMS) and np.isfinite(weights).all()
