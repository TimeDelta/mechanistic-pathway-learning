"""Per-seed forward-cone normalisation (decision 9.16, reversed to per seed on 14 October 2026).

Item 28 divides each edge by the destination's in-degree summed over every relation, so an edge that can carry no
signal for this perturbation still takes a share of the divisor. The owner's proposal: "shouldn't the normalization
be the total in degree for the taken routes on that perturbation?", then "oh ... do the per-seed cone divisor
then", after the union version was shown not to be functionally identical.

What this computes, for one seed: the set of nodes reachable forward from that seed, and then a divisor per
destination counting only the edges whose SOURCE lies in that set. An edge outside the cone can never deliver a
message from this seed, so excluding it from the divisor stops a permanent zero from diluting a live route.

The contraction is not at risk. The encoder's bound comes from the learned gains - "every row of sum_r |S_r| sums
to at most one; gains are bounded by contraction < 1" - and a cone divisor leaves the row sums at most one, because
it divides by the sum of exactly the weights that remain in the row.

What it does NOT fix is decision 9.12's promiscuity: a wrongly coupled edge lies inside its own receptor's cone, so
the divisor keeps it. That needed the coupling filter of 9.23 and is a different problem with the same symptom.
"""
from __future__ import annotations

import torch
from torch import Tensor

TINY = torch.finfo(torch.float32).tiny


def forward_reachable_nodes(edge_source: Tensor, edge_target: Tensor, seed_nodes: Tensor,
                            num_graph_nodes: int, maximum_steps: int | None = None) -> Tensor:
    """Boolean mask of the nodes a message from seed_nodes can reach, the seeds included.

    maximum_steps bounds the walk to the propagation depth the encoder actually runs; None walks to closure, which
    is the right default because the fixed point is the limit of that walk.
    """
    reached = torch.zeros(num_graph_nodes, dtype=torch.bool)
    reached[seed_nodes.reshape(-1).long()] = True
    steps_taken = 0
    while True:
        source_is_reached = reached[edge_source]
        newly_reached = torch.zeros_like(reached)
        newly_reached[edge_target[source_is_reached]] = True
        if bool((newly_reached & ~reached).any()) is False:
            return reached
        reached |= newly_reached
        steps_taken += 1
        if maximum_steps is not None and steps_taken >= maximum_steps:
            return reached


def cone_restricted_in_degree(edge_source: Tensor, edge_target: Tensor, edge_weight: Tensor,
                              reachable_nodes: Tensor, num_graph_nodes: int) -> Tensor:
    """Per-destination divisor counting only the edges whose source lies inside the cone.

    A destination with no edge from inside the cone would divide by zero, so its divisor is clamped to the smallest
    positive float; nothing propagates into such a node from this seed anyway, so the value is never used.
    """
    source_is_in_the_cone = reachable_nodes[edge_source]
    weights_inside_the_cone = torch.where(source_is_in_the_cone, edge_weight.abs(), torch.zeros_like(edge_weight))
    divisor = torch.zeros(num_graph_nodes, dtype=edge_weight.dtype)
    divisor.index_add_(0, edge_target, weights_inside_the_cone)
    return divisor.clamp_min(TINY)


def union_of_seed_cones(edge_source: Tensor, edge_target: Tensor, seed_nodes_per_seed: list[Tensor],
                        num_graph_nodes: int, maximum_steps: int | None = None) -> Tensor:
    """The union version the owner first asked for, kept so the two can be compared rather than assumed equal."""
    all_seed_nodes = torch.cat([seed.reshape(-1) for seed in seed_nodes_per_seed]) if seed_nodes_per_seed else torch.empty(0, dtype=torch.long)
    return forward_reachable_nodes(edge_source, edge_target, all_seed_nodes, num_graph_nodes, maximum_steps)


def divisor_inflation_against_the_union(edge_source: Tensor, edge_target: Tensor, edge_weight: Tensor,
                                        seed_nodes_per_seed: list[Tensor], num_graph_nodes: int,
                                        maximum_steps: int | None = None) -> dict:
    """How much larger the union divisor is than each seed's own, which is what choosing the union costs.

    The two coincide only where a perturbation's seeds share a cone. The ratio is reported over the destinations
    that lie inside the seed's own cone, since those are the only ones whose divisor is used.
    """
    union_reachable = union_of_seed_cones(edge_source, edge_target, seed_nodes_per_seed, num_graph_nodes, maximum_steps)
    union_divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight, union_reachable, num_graph_nodes)
    ratios: list[Tensor] = []
    for seed_nodes in seed_nodes_per_seed:
        seed_reachable = forward_reachable_nodes(edge_source, edge_target, seed_nodes, num_graph_nodes, maximum_steps)
        seed_divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight, seed_reachable, num_graph_nodes)
        destinations_in_use = seed_reachable & (seed_divisor > TINY)
        if bool(destinations_in_use.any()):
            ratios.append(union_divisor[destinations_in_use] / seed_divisor[destinations_in_use])
    if not ratios:
        return {"destinations_compared": 0}
    all_ratios = torch.cat(ratios)
    return {"destinations_compared": int(all_ratios.numel()),
            "median_union_over_per_seed": float(all_ratios.median()),
            "mean_union_over_per_seed": float(all_ratios.mean()),
            "ninetieth_percentile": float(all_ratios.quantile(0.9)),
            "maximum": float(all_ratios.max()),
            "share_of_destinations_where_they_agree": float((all_ratios <= 1.0 + 1e-6).float().mean())}
