"""Negative controls (design section 6.3).

1. permute_symptom_labels_within_degree_strata: shuffles outcome rows among
   perturbations with similar node degree so hub structure is preserved while the
   biology is destroyed.
   degree_stratified_row_permutation returns the row order behind it, so outcomes, pair
   weights and a label mask can move together, optionally only within a partition (the
   lockbox and the development set).
2. degree_preserving_rewiring: randomizes edges within each relation type while
   keeping every node's in- and out-degree (double-edge swaps).
   fast_degree_preserving_rewiring applies the same swap rule at O(1) per swap, for the
   full graph (267k edges at 50 swaps per edge).
3. Grade shuffling is a one-line permutation of the weight column in the training loop.
4. The peripheral-event control is a data selection (events with no central
   mechanism) handled by the evidence loaders.
"""
from __future__ import annotations

import numpy as np


def degree_stratified_row_permutation(perturbation_degrees: np.ndarray, num_strata: int = 5, random_seed: int = 0, partition: np.ndarray | None = None) -> np.ndarray:
    """Source row for every row: rows are shuffled among perturbations of the same degree quintile (quantiles over all
    rows) and, with a partition, of the same partition value, so no label crosses from the lockbox to the development
    set. Without a partition the order equals the one permute_symptom_labels_within_degree_strata has always used."""
    generator = np.random.default_rng(random_seed)
    stratum_edges = np.quantile(perturbation_degrees, np.linspace(0.0, 1.0, num_strata + 1))
    stratum_of_row = np.clip(np.searchsorted(stratum_edges, perturbation_degrees, side="right") - 1, 0, num_strata - 1)
    source_row = np.arange(len(perturbation_degrees))
    cells = [None] if partition is None else sorted(set(np.asarray(partition).tolist()))
    for cell in cells:
        in_cell = np.ones(len(source_row), dtype=bool) if cell is None else np.asarray(partition) == cell
        for stratum_index in range(num_strata):
            rows = np.where(in_cell & (stratum_of_row == stratum_index))[0]
            source_row[rows] = generator.permutation(rows)
    return source_row


def permute_symptom_labels_within_degree_strata(outcomes: np.ndarray, perturbation_degrees: np.ndarray, num_strata: int = 5, random_seed: int = 0) -> np.ndarray:
    return outcomes[degree_stratified_row_permutation(perturbation_degrees, num_strata, random_seed)]


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


def fast_degree_preserving_rewiring(edge_index: np.ndarray, edge_relation_type: np.ndarray, num_swaps_per_edge: int = 10, random_seed: int = 0,
                                    draws_per_batch: int = 1_000_000) -> np.ndarray:
    """Degree-preserving rewiring within each relation at O(1) per swap, for the full graph. Two edges of one relation
    exchange targets unless that makes a self-loop or repeats an edge already present, the Maslov and Sneppen rule
    (Science 2002, doi:10.1126/science.1065103: "we do not allow ... multiple edges"); degree_preserving_rewiring
    allows the repeat. Edge pairs are drawn as integer positions in batches; num_swaps_per_edge counts attempts,
    rejected ones and a drawn pair of one edge with itself included. The random stream differs from
    degree_preserving_rewiring's, so one seed gives different graphs."""
    generator = np.random.default_rng(random_seed)
    rewired = edge_index.copy()
    for relation_index in np.unique(edge_relation_type):
        edge_positions = np.where(edge_relation_type == relation_index)[0]
        num_edges = len(edge_positions)
        if num_edges < 2:
            continue
        sources = rewired[0, edge_positions].tolist()
        targets = rewired[1, edge_positions].tolist()
        edge_multiplicity: dict[tuple[int, int], int] = {}
        for edge in zip(sources, targets):
            edge_multiplicity[edge] = edge_multiplicity.get(edge, 0) + 1
        attempts_left = num_swaps_per_edge * num_edges
        while attempts_left > 0:
            num_draws = min(attempts_left, draws_per_batch)
            first_edges = generator.integers(0, num_edges, size=num_draws).tolist()
            second_edges = generator.integers(0, num_edges, size=num_draws).tolist()
            for first, second in zip(first_edges, second_edges):
                source_first, source_second = sources[first], sources[second]
                target_first, target_second = targets[first], targets[second]
                if first == second or source_first == target_second or source_second == target_first:
                    continue
                if (source_first, target_second) in edge_multiplicity or (source_second, target_first) in edge_multiplicity:
                    continue
                for removed in ((source_first, target_first), (source_second, target_second)):
                    edge_multiplicity[removed] -= 1
                    if edge_multiplicity[removed] == 0:
                        del edge_multiplicity[removed]
                edge_multiplicity[(source_first, target_second)] = 1
                edge_multiplicity[(source_second, target_first)] = 1
                targets[first], targets[second] = target_second, target_first
            attempts_left -= num_draws
        rewired[1, edge_positions] = targets
    return rewired


def duplicate_edge_count(edge_index: np.ndarray, edge_relation_type: np.ndarray) -> int:
    """Edges repeating an earlier (source, target, relation) triple; rewiring can create them."""
    triples = np.stack([edge_index[0], edge_index[1], edge_relation_type], axis=1)
    return int(len(triples) - len(np.unique(triples, axis=0)))
