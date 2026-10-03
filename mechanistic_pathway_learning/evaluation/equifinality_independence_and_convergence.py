"""Equifinality indices (experiment design, section 6.5).

For a symptom s with active modules M_s = {k : link_{k,s} > threshold}:

- independence_index = 1 - max pairwise Jaccard overlap of the module supports.
  Values near 1 mean the modules occupy disjoint regions of the graph.
- convergence_index = mean pairwise Jaccard overlap of the downstream node sets
  reachable from each support within max_hops. Values near 1 mean the modules
  converge on a shared downstream region, the final common pathway pattern.

Both are reported for every symptom; neither outcome is treated as failure.
The sufficiency test (module ablation on held-out perturbations) lives in the
evaluation loop because it needs model predictions.
"""
from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from itertools import combinations


def jaccard_index(set_a: set, set_b: set) -> float:
    union_size = len(set_a | set_b)
    if union_size == 0:
        return 0.0
    return len(set_a & set_b) / union_size


def active_modules_for_symptom(link_probability_by_module: Mapping[str, float], link_threshold: float) -> list[str]:
    return [module_identifier for module_identifier, link_probability in link_probability_by_module.items() if link_probability > link_threshold]


def independence_index(module_supports: Mapping[str, set], active_modules: Sequence[str]) -> float:
    """1 - max pairwise Jaccard of supports; NaN when fewer than two modules are active."""
    if len(active_modules) < 2:
        return math.nan
    max_overlap = max(jaccard_index(module_supports[first], module_supports[second]) for first, second in combinations(active_modules, 2))
    return 1.0 - max_overlap


def downstream_reachable_nodes(adjacency: Mapping[str, Iterable[str]], seed_nodes: Iterable[str], max_hops: int) -> set[str]:
    """Nodes reachable from the seeds within max_hops directed steps, excluding the seeds themselves."""
    frontier = set(seed_nodes)
    visited = set(frontier)
    for _ in range(max_hops):
        next_frontier: set[str] = set()
        for node in frontier:
            for neighbor in adjacency.get(node, ()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    next_frontier.add(neighbor)
        frontier = next_frontier
        if not frontier:
            break
    return visited - set(seed_nodes)


def convergence_index(
    module_supports: Mapping[str, set],
    active_modules: Sequence[str],
    adjacency: Mapping[str, Iterable[str]],
    max_hops: int,
) -> float:
    """Mean pairwise Jaccard of downstream reachable sets; NaN when fewer than two modules are active."""
    if len(active_modules) < 2:
        return math.nan
    downstream_by_module = {
        module_identifier: downstream_reachable_nodes(adjacency, module_supports[module_identifier], max_hops)
        for module_identifier in active_modules
    }
    pairwise_overlaps = [
        jaccard_index(downstream_by_module[first], downstream_by_module[second]) for first, second in combinations(active_modules, 2)
    ]
    return sum(pairwise_overlaps) / len(pairwise_overlaps)
