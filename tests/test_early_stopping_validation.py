"""Tests for the early-stopping validation split (experiments/run_main_model.py, early_stopping_validation): without
--keep-large-groups-in-training the largest leakage group goes to the validation set first, with it a group larger
than half the expected validation set stays in training, and when no group is that large the flag changes nothing."""
from argparse import Namespace
from types import SimpleNamespace

import numpy as np

from experiments.run_main_model import early_stopping_validation


def toy_data(group_sizes: list[int]) -> SimpleNamespace:
    perturbation_ids, group_ids = [], []
    for group_index, size in enumerate(group_sizes):
        for member in range(size):
            perturbation_ids.append(f"P{group_index}_{member}")
            group_ids.append(f"group{group_index}")
    return SimpleNamespace(perturbation_ids=perturbation_ids, group_ids=group_ids)


def split(data, keep_large_groups_in_training: bool, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    arguments = Namespace(validation_fraction=0.15, seed=seed, keep_large_groups_in_training=keep_large_groups_in_training)
    return early_stopping_validation(data, np.arange(len(data.perturbation_ids)), arguments)


def test_without_the_flag_one_large_group_is_the_whole_validation_set() -> None:
    data = toy_data([60] + [2] * 120)  # 300 perturbations; the expected validation set is 300 / 7, about 43
    for seed in range(3):
        _, validation = split(data, keep_large_groups_in_training=False, seed=seed)
        assert {data.group_ids[i] for i in validation} == {"group0"}


def test_with_the_flag_the_large_group_stays_in_training() -> None:
    data = toy_data([60] + [2] * 120)
    training, validation = split(data, keep_large_groups_in_training=True)
    assert "group0" not in {data.group_ids[i] for i in validation}
    assert all(f"P0_{member}" in {data.perturbation_ids[i] for i in training} for member in range(60))
    assert len(validation) > 0
    assert not {data.group_ids[i] for i in validation} & {data.group_ids[i] for i in training}
    assert sorted(np.concatenate([training, validation]).tolist()) == list(range(300))


def test_with_no_large_group_the_flag_changes_nothing() -> None:
    data = toy_data([12, 9, 7] + [2] * 130)
    for seed in range(3):
        without_flag = split(data, keep_large_groups_in_training=False, seed=seed)
        with_flag = split(data, keep_large_groups_in_training=True, seed=seed)
        assert all(np.array_equal(a, b) for a, b in zip(without_flag, with_flag))
