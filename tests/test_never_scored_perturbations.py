"""Preregistration item 12c: a receptor-activation row is never scored, and "training only" is not "always trained"."""
from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from run_main_model import never_scored_rows, withhold_never_scored_rows

from mechanistic_pathway_learning.evidence.combined_perturbations import RECEPTOR_ACTIVATION_PERTURBATION_TYPE

ACTIVATION = RECEPTOR_ACTIVATION_PERTURBATION_TYPE


def data_with(types: list[str], groups: list[str]) -> SimpleNamespace:
    return SimpleNamespace(perturbation_types=types, group_ids=groups,
                           perturbation_ids=[f"P{i}" for i in range(len(types))])


def test_a_never_scored_row_leaves_the_test_set_even_when_its_group_is_drawn():
    data = data_with(["gene", ACTIVATION], ["g1", "g1"])
    kept_test, _, taken_out, _ = withhold_never_scored_rows(data, np.array([True, True]), np.array([], dtype=int))
    assert kept_test.tolist() == [True, False] and taken_out == 1


def test_a_never_scored_row_is_dropped_from_training_when_its_group_is_held_out():
    """The leak the rule exists for: training the activation row of a gene whose knockout is in the lockbox shows the
    model that node perturbed, which is what the leakage groups exist to stop."""
    data = data_with(["gene", ACTIVATION, "drug"], ["g1", "g1", "g2"])
    _, train_pool, _, dropped = withhold_never_scored_rows(data, np.array([True, False, False]), np.array([1, 2]))
    assert train_pool.tolist() == [2] and dropped == 1


def test_a_never_scored_row_trains_when_its_group_is_not_held_out():
    data = data_with(["gene", ACTIVATION, "drug"], ["g1", "g2", "g3"])
    _, train_pool, _, dropped = withhold_never_scored_rows(data, np.array([True, False, False]), np.array([1, 2]))
    assert train_pool.tolist() == [1, 2] and dropped == 0


def test_the_held_out_set_over_every_other_perturbation_is_untouched():
    data = data_with(["gene", "drug", ACTIVATION, "gene"], ["g1", "g2", "g3", "g4"])
    kept_test, train_pool, _, _ = withhold_never_scored_rows(data, np.array([True, True, False, False]), np.array([2, 3]))
    assert kept_test.tolist() == [True, True, False, False]  # the ordinary rows are exactly as drawn
    assert train_pool.tolist() == [2, 3]


def test_nothing_changes_when_no_row_is_of_a_never_scored_type():
    data = data_with(["gene", "drug"], ["g1", "g2"])
    kept_test, train_pool, taken_out, dropped = withhold_never_scored_rows(data, np.array([True, False]), np.array([1]))
    assert kept_test.tolist() == [True, False] and train_pool.tolist() == [1] and (taken_out, dropped) == (0, 0)
    assert not never_scored_rows(data).any()
