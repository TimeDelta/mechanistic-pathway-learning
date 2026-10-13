"""Preregistration item 26: the nosology readout's clustering, its stability reading and its controls."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from nosology_readout import (
    cluster_at,
    cluster_count_by_stability,
    correlation_distances,
    link_concentration,
    partition_agreement,
)


def profiles_with_known_groups(num_coordinates: int = 8, seed: int = 0) -> np.ndarray:
    """Six symptoms in three pairs: within a pair the profiles agree up to noise, across pairs they do not."""
    generator = np.random.default_rng(seed)
    group_profiles = generator.normal(size=(num_coordinates, 3))
    return np.concatenate([group_profiles[:, [group, group]] + 0.01 * generator.normal(size=(num_coordinates, 2))
                           for group in range(3)], axis=1)


def test_the_distance_is_one_minus_correlation_and_a_flat_profile_is_far_from_everything():
    profiles = np.array([[1.0, 2.0, 0.0], [2.0, 4.0, 0.0], [3.0, 6.0, 0.0]])
    distances = correlation_distances(profiles)
    assert distances[0, 1] == pytest.approx(0.0, abs=1e-12)  # the second symptom is twice the first
    assert np.allclose(np.diag(distances), 0.0)
    assert distances[0, 2] == pytest.approx(1.0) and distances[2, 1] == pytest.approx(1.0)  # the flat one
    assert np.allclose(distances, distances.T)


def test_two_flat_profiles_are_far_from_each_other_and_not_a_cluster():
    """Item 26f: symptoms with no profile must not group together by having none. Correlation is undefined between
    two constant columns, and filling it in as 1 would make every thin symptom one cluster."""
    profiles = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    distances = correlation_distances(profiles)
    assert distances[1, 2] == pytest.approx(1.0)
    assert distances[1, 1] == pytest.approx(0.0) and distances[2, 2] == pytest.approx(0.0)


def test_the_clustering_recovers_groups_that_are_there():
    labels = cluster_at(correlation_distances(profiles_with_known_groups()), 3)
    assert [labels[0] == labels[1], labels[2] == labels[3], labels[4] == labels[5]] == [True, True, True]
    assert len(set(labels)) == 3


def test_partition_agreement_is_one_for_a_relabelling_and_falls_when_the_grouping_differs():
    assert partition_agreement(np.array([1, 1, 2, 2]), np.array([7, 7, 3, 3])) == pytest.approx(1.0)
    assert partition_agreement(np.array([1, 1, 2, 2]), np.array([1, 2, 1, 2])) < 0.5


def test_the_number_of_clusters_is_taken_by_stability_and_finds_the_planted_count():
    best, agreement_by_count = cluster_count_by_stability(profiles_with_known_groups(num_coordinates=16), seed=1)
    assert best == 3
    assert agreement_by_count[3] >= agreement_by_count[max(agreement_by_count)]


def test_noise_alone_is_less_stable_than_structure():
    """Ben-Hur et al.'s second claim, that the criterion 'can also detect the lack of structure in data': the planted
    grouping has to score above profiles with nothing in them."""
    generator = np.random.default_rng(3)
    _, structured = cluster_count_by_stability(profiles_with_known_groups(num_coordinates=16), seed=1)
    _, unstructured = cluster_count_by_stability(generator.normal(size=(16, 6)), seed=1)
    assert max(structured.values()) > max(unstructured.values())


def test_link_concentration_is_one_for_a_single_symptom_module_and_even_for_a_flat_one():
    links = np.array([[0.4, 0.0, 0.0, 0.0], [0.1, 0.1, 0.1, 0.1], [0.0, 0.0, 0.0, 0.0]])
    concentration = link_concentration(links)
    assert concentration[0] == pytest.approx(1.0)
    assert concentration[1] == pytest.approx(0.25)
    assert np.isnan(concentration[2])  # a module with no link mass is not a concentrated module
