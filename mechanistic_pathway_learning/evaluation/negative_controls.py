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


def reciprocated_relations(edge_index: np.ndarray, edge_relation_type: np.ndarray, minimum_share: float = 0.9) -> list[int]:
    """Relations in which at least minimum_share of the edges have their reverse stored in the same relation: undirected
    relations stored in both directions (binds in the full graph, 0.986)."""
    relations = []
    for relation_index in np.unique(edge_relation_type):
        positions = np.where(edge_relation_type == relation_index)[0]
        edges = set(zip(edge_index[0, positions].tolist(), edge_index[1, positions].tolist()))
        if edges and sum((target, source) in edges for source, target in edges) / len(edges) >= minimum_share:
            relations.append(int(relation_index))
    return relations


def _rewire_undirected_units(unit_first: list[int], unit_second: list[int], occupied: dict, num_attempts: int, generator, draws_per_batch: int) -> None:
    """Double-edge swaps on undirected edges {a, b} and {c, d}, which become {a, d} and {c, b}, or {a, c} and {d, b}
    when the second edge is read the other way round (one random bit per draw). A swap that makes a self-loop, or an
    edge already present in either direction (occupied counts directed edges, both directions of every undirected one),
    is rejected. Every node keeps its undirected degree. Lists are changed in place."""
    num_units = len(unit_first)

    def add(edge):
        occupied[edge] = occupied.get(edge, 0) + 1

    def remove(edge):
        occupied[edge] -= 1
        if occupied[edge] == 0:
            del occupied[edge]

    while num_attempts > 0:
        num_draws = min(num_attempts, draws_per_batch)
        first_units = generator.integers(0, num_units, size=num_draws).tolist()
        second_units = generator.integers(0, num_units, size=num_draws).tolist()
        reversed_reading = generator.integers(0, 2, size=num_draws).tolist()
        for first, second, read_reversed in zip(first_units, second_units, reversed_reading):
            if first == second:
                continue
            a, b = unit_first[first], unit_second[first]
            c, d = (unit_second[second], unit_first[second]) if read_reversed else (unit_first[second], unit_second[second])
            if a == d or c == b or (a, d) in occupied or (d, a) in occupied or (c, b) in occupied or (b, c) in occupied:
                continue
            for edge in ((a, b), (b, a), (c, d), (d, c)):
                remove(edge)
            for edge in ((a, d), (d, a), (c, b), (b, c)):
                add(edge)
            unit_second[first], unit_first[second], unit_second[second] = d, c, b
        num_attempts -= num_draws


def fast_degree_preserving_rewiring(edge_index: np.ndarray, edge_relation_type: np.ndarray, num_swaps_per_edge: int = 10, random_seed: int = 0,
                                    draws_per_batch: int = 1_000_000, undirected_relations=(), fixed_relations=()) -> np.ndarray:
    """Degree-preserving rewiring within each relation at O(1) per swap, for the full graph. Two edges of one relation
    exchange targets unless that makes a self-loop or repeats an edge already present, as in the randomisation of
    Maslov and Sneppen (Science 2002, doi:10.1126/science.1065103); degree_preserving_rewiring allows the repeat. Edge pairs are drawn as integer positions in batches; num_swaps_per_edge counts attempts,
    rejected ones and a drawn pair of one edge with itself included. The random stream differs from
    degree_preserving_rewiring's, so one seed gives different graphs.

    undirected_relations (reciprocated_relations finds them): in these relations an edge stored with its reverse is one
    undirected edge, rewired by undirected double-edge swaps and written back in both directions, so the rewired
    relation stays symmetric where the real one is (rewiring the two directions apart leaves 1,664 of the full graph's
    31,118 binds edges with their reverse, from 30,682); its edges without a reverse are rewired as directed edges. Every
    node keeps its in- and out-degree in every relation either way. With none given the output is the one this function
    has always produced for a seed.

    fixed_relations: relations left as they are, such as encodes on a split graph (gene -> its own protein; swapping it
    would join a gene to another gene's protein). The random stream of the other relations is unchanged when the fixed
    relations come after them in index order, as encodes does (appended last by graph/gene_protein_split.py)."""
    generator = np.random.default_rng(random_seed)
    rewired = edge_index.copy()
    undirected_relations = {int(relation) for relation in undirected_relations}
    fixed_relations = {int(relation) for relation in fixed_relations}
    for relation_index in np.unique(edge_relation_type):
        if int(relation_index) in fixed_relations:
            continue
        edge_positions = np.where(edge_relation_type == relation_index)[0]
        num_edges = len(edge_positions)
        if num_edges < 2:
            continue
        if int(relation_index) in undirected_relations:
            _rewire_reciprocated_relation(rewired, edge_positions, num_swaps_per_edge, generator, draws_per_batch)
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


def _rewire_reciprocated_relation(rewired: np.ndarray, edge_positions: np.ndarray, num_swaps_per_edge: int, generator, draws_per_batch: int) -> None:
    """One undirected relation: pairs of reverse edges become undirected units rewired by undirected swaps, then the
    edges without a reverse are rewired as directed ones; both see the same occupied edges, so neither duplicates the
    other. rewired is changed in place at edge_positions."""
    positions_by_edge: dict[tuple[int, int], list[int]] = {}
    for position in edge_positions.tolist():
        positions_by_edge.setdefault((int(rewired[0, position]), int(rewired[1, position])), []).append(position)
    unit_positions, directed_positions = [], []
    for (source, target), positions in sorted(positions_by_edge.items()):
        reverse = positions_by_edge.get((target, source), [])
        if source < target:  # each reverse pair once, from its lower endpoint
            num_pairs = min(len(positions), len(reverse))
            unit_positions += list(zip(positions[:num_pairs], reverse[:num_pairs]))
            directed_positions += positions[num_pairs:]
        elif source > target:
            directed_positions += positions[min(len(positions), len(reverse)):]
        else:
            directed_positions += positions
    occupied: dict[tuple[int, int], int] = {}
    for position in edge_positions.tolist():
        edge = (int(rewired[0, position]), int(rewired[1, position]))
        occupied[edge] = occupied.get(edge, 0) + 1
    unit_first = [int(rewired[0, forward]) for forward, _ in unit_positions]
    unit_second = [int(rewired[1, forward]) for forward, _ in unit_positions]
    if len(unit_positions) >= 2:
        _rewire_undirected_units(unit_first, unit_second, occupied, num_swaps_per_edge * len(unit_positions), generator, draws_per_batch)
    for (forward, backward), first, second in zip(unit_positions, unit_first, unit_second):
        rewired[0, forward], rewired[1, forward] = first, second
        rewired[0, backward], rewired[1, backward] = second, first
    sources = [int(rewired[0, position]) for position in directed_positions]
    targets = [int(rewired[1, position]) for position in directed_positions]
    attempts_left = num_swaps_per_edge * len(directed_positions) if len(directed_positions) >= 2 else 0
    while attempts_left > 0:
        num_draws = min(attempts_left, draws_per_batch)
        first_edges = generator.integers(0, len(directed_positions), size=num_draws).tolist()
        second_edges = generator.integers(0, len(directed_positions), size=num_draws).tolist()
        for first, second in zip(first_edges, second_edges):
            source_first, source_second, target_first, target_second = sources[first], sources[second], targets[first], targets[second]
            if first == second or source_first == target_second or source_second == target_first:
                continue
            if (source_first, target_second) in occupied or (source_second, target_first) in occupied:
                continue
            if (target_second, source_first) in occupied or (target_first, source_second) in occupied:  # would form a reverse pair the real relation lacks
                continue
            for removed in ((source_first, target_first), (source_second, target_second)):
                occupied[removed] -= 1
                if occupied[removed] == 0:
                    del occupied[removed]
            occupied[(source_first, target_second)] = 1
            occupied[(source_second, target_first)] = 1
            targets[first], targets[second] = target_second, target_first
        attempts_left -= num_draws
    rewired[1, directed_positions] = targets


def rewire_walk_graph(num_nodes: int, source_index: np.ndarray, target_index: np.ndarray, excluded_nodes: np.ndarray | None = None,
                      num_swaps_per_edge: int = 50, random_seed: int = 0, draws_per_batch: int = 1_000_000) -> np.ndarray:
    """The random walk's own graph rewired: the undirected simple graph build_normalized_adjacency makes (relations
    merged, both directions, each node pair once, excluded nodes dropped) rewired by undirected double-edge swaps, so
    every node keeps the number of neighbours the walk sees. Returns its edges [2, num_pairs], each pair once.

    Rewiring the stored directed edges per relation and collapsing afterwards does not keep that number: two stored
    edges joining one pair (binds in both directions, a reversible reaction's substrate and product edges) are moved
    apart, so the walk's graph gains neighbours (slice: 112,920 adjacency entries on the real graph against about
    135,000 after rewiring)."""
    keep = np.ones(len(source_index), dtype=bool)
    if excluded_nodes is not None and len(excluded_nodes):
        excluded = np.zeros(num_nodes, dtype=bool)
        excluded[np.asarray(excluded_nodes)] = True
        keep = ~(excluded[source_index] | excluded[target_index])
    first, second = np.minimum(source_index[keep], target_index[keep]), np.maximum(source_index[keep], target_index[keep])
    pairs = np.unique(np.stack([first, second], axis=1), axis=0)
    self_loops = pairs[pairs[:, 0] == pairs[:, 1]]  # kept where they are
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    unit_first, unit_second = pairs[:, 0].tolist(), pairs[:, 1].tolist()
    occupied: dict[tuple[int, int], int] = {}
    for a, b in zip(unit_first, unit_second):
        occupied[(a, b)] = 1
        occupied[(b, a)] = 1
    for loop_node in self_loops[:, 0].tolist():
        occupied[(loop_node, loop_node)] = 1
    if len(unit_first) >= 2:
        _rewire_undirected_units(unit_first, unit_second, occupied, num_swaps_per_edge * len(unit_first), np.random.default_rng(random_seed), draws_per_batch)
    return np.array([unit_first + self_loops[:, 0].tolist(), unit_second + self_loops[:, 1].tolist()], dtype=np.int64)


def duplicate_edge_count(edge_index: np.ndarray, edge_relation_type: np.ndarray) -> int:
    """Edges repeating an earlier (source, target, relation) triple; rewiring can create them."""
    triples = np.stack([edge_index[0], edge_index[1], edge_relation_type], axis=1)
    return int(len(triples) - len(np.unique(triples, axis=0)))
