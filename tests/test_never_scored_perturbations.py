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


def test_a_never_scored_row_does_not_reach_the_validation_set_either():
    """The validation set picks the stopping epoch by macro AUPRC, which is a score. Item 12c keeps these rows out
    of every reading because macro AUPRC is taken per symptom and 43 percent of their labels are one symptom, so a
    row that reaches validation chooses the epoch on the very quantity it is excluded from."""
    import numpy as np
    from argparse import Namespace

    from run_main_model import split_indices

    num_perturbations = 40
    types = ["gene"] * num_perturbations
    types[5] = types[11] = ACTIVATION
    outcomes = np.zeros((num_perturbations, 3))
    outcomes[:, 0] = 1.0
    data = SimpleNamespace(
        perturbation_ids=[f"P{i}" for i in range(num_perturbations)], perturbation_types=types,
        # each never-scored row shares its gene's group, as a real receptor-activation row does; with a group of its
        # own the group would vanish from the test set the moment the row was taken out of it
        group_ids=[f"g{i - 1 if i in (5, 11) else i}" for i in range(num_perturbations)], outcomes=outcomes,
        label_mask=np.ones_like(outcomes, dtype=bool), perturbation_seeds=[[0]] * num_perturbations)
    arguments = Namespace(time_split_cutoff=None, holdout_module="", holdout_subsystem="", fold=0, num_folds=5,
                          validation_fraction=0.2, seed=0, group_by="disease_cluster")
    train, validation, test, _ = split_indices(data, arguments)
    never_scored = {5, 11}
    assert not (set(validation.tolist()) & never_scored), "a never-scored row reached the validation set"
    assert not (set(test.tolist()) & never_scored), "a never-scored row reached the test set"
    # Each such row is either trained on or dropped because its own group is held out; neither is a loss, and a row
    # whose group is held out MUST be dropped rather than trained, which is the leak the rule exists for.
    held_out_groups = {data.group_ids[index] for index in test.tolist()}  # visible because the partner row stays in test
    for index in never_scored:
        in_training = index in set(train.tolist())
        assert in_training != (data.group_ids[index] in held_out_groups), (
            f"P{index} is {'trained' if in_training else 'dropped'} while its group is "
            f"{'held out' if data.group_ids[index] in held_out_groups else 'not held out'}")
