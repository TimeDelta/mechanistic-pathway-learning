"""Negative controls (design section 6.3).

1. permute_symptom_labels_within_degree_strata: shuffles outcome rows among
   perturbations with similar node degree so hub structure is preserved while the
   biology is destroyed.
2. degree_preserving_rewiring: randomizes edges within each relation type while
   keeping every node's in- and out-degree (double-edge swaps).
3. Grade shuffling is a one-line permutation of the weight column in the training loop.
4. The peripheral-event control is a data selection (events with no central
   mechanism) handled by the evidence loaders.
"""
from __future__ import annotations

import numpy as np


def permute_symptom_labels_within_degree_strata(outcomes: np.ndarray, perturbation_degrees: np.ndarray, num_strata: int = 5, random_seed: int = 0) -> np.ndarray:
    generator = np.random.default_rng(random_seed)
    permuted = outcomes.copy()
    stratum_edges = np.quantile(perturbation_degrees, np.linspace(0.0, 1.0, num_strata + 1))
    stratum_of_row = np.clip(np.searchsorted(stratum_edges, perturbation_degrees, side="right") - 1, 0, num_strata - 1)
    for stratum_index in range(num_strata):
        rows = np.where(stratum_of_row == stratum_index)[0]
        permuted[rows] = outcomes[generator.permutation(rows)]
    return permuted


def degree_preserving_rewiring(edge_index: np.ndarray, edge_relation_type: np.ndarray, num_swaps_per_edge: int = 10, random_seed: int = 0) -> np.ndarray:
    """Return a rewired edge_index [2, num_edges]; each relation type is rewired separately."""
    generator = np.random.default_rng(random_seed)
    rewired = edge_index.copy()
    for relation_index in np.unique(edge_relation_type):
        edge_positions = np.where(edge_relation_type == relation_index)[0]
        num_edges = len(edge_positions)
        if num_edges < 2:
            continue
        for _ in range(num_swaps_per_edge * num_edges):
            first, second = generator.choice(edge_positions, size=2, replace=False)
            source_first, target_first = rewired[:, first]
            source_second, target_second = rewired[:, second]
            if source_first == target_second or source_second == target_first:
                continue
            rewired[1, first], rewired[1, second] = target_second, target_first
    return rewired
