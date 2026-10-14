"""Preregistration item 20, built 14 October 2026: the per-edge stoichiometric demand.

Until now `bottleneck_edge_demand` was None and every demand was 1, so a reaction consuming two molecules counted
as one - the owner's acetylcholine example. The graph build now writes `stoichiometric_coefficient` on its reaction
edges and the encoder maps it onto the bottleneck adjacency.
"""
from __future__ import annotations

import torch

from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder


def small_encoder() -> tuple[LinearResponseEncoder, dict]:
    """Four nodes: two substrates into one reaction, which makes one product."""
    relation_types = ["substrate_of", "product_of"]
    edges = {
        "edge_source": torch.tensor([0, 1, 2]),
        "edge_target": torch.tensor([2, 2, 3]),
        "edge_relation": torch.tensor([0, 0, 1]),
        "edge_sign": torch.tensor([1.0, 1.0, 1.0]),
    }
    encoder = LinearResponseEncoder(num_graph_nodes=4, relation_types=relation_types,
                                    node_features=torch.zeros((4, 1)), node_state_dim=2, **edges)
    return encoder, edges


def test_the_demand_starts_at_one_everywhere():
    encoder, _ = small_encoder()
    assert encoder.bottleneck_edge_demand is None, "the registered default before item 20 was built"


def test_a_stated_coefficient_reaches_the_bottleneck_entry():
    """Two molecules of node 0 per turnover, one of node 1."""
    encoder, edges = small_encoder()
    counts = encoder.set_bottleneck_edge_demand_from_edges(
        edges["edge_source"], edges["edge_target"], edges["edge_relation"],
        torch.tensor([2.0, 1.0, 3.0]))
    demand = encoder.bottleneck_edge_demand
    assert demand is not None and demand.shape == encoder.bottleneck_edge_source.shape
    stated = {}
    for position in range(len(demand)):
        key = (int(encoder.bottleneck_edge_relation[position]), int(encoder.bottleneck_edge_source[position]),
               int(encoder.bottleneck_edge_target[position]))
        stated[key] = float(demand[position])
    assert stated[(0, 0, 2)] == 2.0, "node 0 is consumed twice per turnover"
    assert stated[(0, 1, 2)] == 1.0
    assert stated[(1, 2, 3)] == 3.0
    assert counts["entries_above_one"] == 2
    assert counts["largest_demand"] == 3.0


def test_a_derived_relation_keeps_demand_one():
    """`depletes_substrate` and `product_opposes` are built by the encoder, not stated by the graph."""
    encoder, edges = small_encoder()
    encoder.set_bottleneck_edge_demand_from_edges(
        edges["edge_source"], edges["edge_target"], edges["edge_relation"], torch.tensor([2.0, 1.0, 3.0]))
    stated_relations = {0, 1}
    for position in range(len(encoder.bottleneck_edge_demand)):
        if int(encoder.bottleneck_edge_relation[position]) not in stated_relations:
            assert float(encoder.bottleneck_edge_demand[position]) == 1.0


def test_a_mismatched_coefficient_count_is_refused():
    encoder, edges = small_encoder()
    try:
        encoder.set_bottleneck_edge_demand_from_edges(edges["edge_source"], edges["edge_target"],
                                                      edges["edge_relation"], torch.tensor([2.0]))
    except ValueError as error:
        assert "one per edge" in str(error)
    else:
        raise AssertionError("a coefficient per edge is required")


def test_a_fractional_coefficient_is_clamped_to_one_and_counted():
    """Human-GEM writes half a molecule of O2. A demand below 1 would raise availability across a step and take the
    bottleneck out of [0, 1], which `set_bottleneck_edge_demand` refuses, so the fraction is clamped and reported -
    clamping being a deviation from the stoichiometry rather than a reading of it."""
    encoder, edges = small_encoder()
    counts = encoder.set_bottleneck_edge_demand_from_edges(edges["edge_source"], edges["edge_target"],
                                                            edges["edge_relation"], torch.tensor([0.5, 0.25, 1.0]))
    assert float(encoder.bottleneck_edge_demand.min()) >= 1.0
    assert counts["entries_clamped_up_from_a_fraction"] == 2
    assert counts["entries_above_one"] == 0
