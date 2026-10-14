"""Decision 4.11, built 14 October 2026: the linear-response encoder reads edge_weight in its propagation adjacency.

Before this, edge_weight scaled a drug's mechanism (entry) edges only, so item 13's six conductance weights were
inert in this encoder while the message-passing encoder read them, and the two registered encoders disagreed about
the physics. The weight enters as a **share**: the edge and its destination's divisor are weighted together, so an
unsigned row still sums to 1 (item 28) and the quantity is a conductance fraction, which is what the Millman
equation the electrical layer linearises actually uses.
"""
from __future__ import annotations

import pytest
import torch

from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

RELATIONS = ["activates", "inhibits"]


def adjacency(edge_weight=None, normalisation="total_in_degree"):
    """Two sources feeding one target, so a weight changes their shares against each other."""
    return LinearResponseEncoder.signed_stacked_adjacency(
        num_graph_nodes=3, relation_types=RELATIONS,
        edge_source=torch.tensor([0, 1]), edge_target=torch.tensor([2, 2]),
        edge_relation=torch.tensor([0, 0]), edge_sign=torch.tensor([1.0, 1.0]),
        normalisation=normalisation, edge_weight=edge_weight)[0].coalesce()


def test_no_weight_reproduces_the_unweighted_adjacency():
    """A graph with no weight column must build exactly the adjacency it always built."""
    for normalisation in ("total_in_degree", "spectral", "per_relation_in_degree"):
        try:
            unweighted = adjacency(None, normalisation)
            ones = adjacency(torch.tensor([1.0, 1.0]), normalisation)
        except ValueError:
            continue  # a normalisation this build does not offer
        assert torch.allclose(unweighted.values(), ones.values()), normalisation


def test_an_unsigned_row_still_sums_to_one_under_any_weighting():
    """Item 28's property: the same contraction and no global shrink."""
    for weight in ([1.0, 1.0], [6.0, 1.0], [0.25, 3.5], [100.0, 0.001]):
        values = adjacency(torch.tensor(weight)).values().abs()
        assert float(values.sum()) == pytest.approx(1.0, abs=1e-6), f"row sum moved under {weight}"


def test_the_weight_acts_as_a_share_between_edges_feeding_one_node():
    """A conductance of 6 against 1 gives 6/7 and 1/7, not 6 and 1."""
    values = adjacency(torch.tensor([6.0, 1.0])).values()
    shares = sorted(float(v) for v in values)
    assert shares[1] == pytest.approx(6.0 / 7.0, abs=1e-6)
    assert shares[0] == pytest.approx(1.0 / 7.0, abs=1e-6)


def test_a_mismatched_weight_count_is_refused():
    try:
        adjacency(torch.tensor([1.0]))
    except ValueError as error:
        assert "one value per edge" in str(error)
    else:
        raise AssertionError("one weight per edge is required")


def test_a_derived_relation_inherits_its_parent_edges_weight():
    """`depletes_substrate` is the substrate edge read backwards, so it carries that edge's conductance.

    One metabolite (node 0) is the substrate of two reactions (nodes 1 and 2) at weights 9 and 1, so both depletion
    edges land on node 0 and their shares must be 0.9 and 0.1. A target with only one inbound edge would normalise
    to 1 whatever its weight, which is the share design working and would hide the inheritance."""
    stacked, names = LinearResponseEncoder.signed_stacked_adjacency(
        num_graph_nodes=3, relation_types=["substrate_of"],
        edge_source=torch.tensor([0, 0]), edge_target=torch.tensor([1, 2]),
        edge_relation=torch.tensor([0, 0]), edge_sign=torch.tensor([1.0, 1.0]),
        normalisation="total_in_degree", edge_weight=torch.tensor([9.0, 1.0]))
    assert "depletes_substrate" in names
    coalesced = stacked.coalesce()
    depletion_row = names.index("depletes_substrate") * 3 + 0  # into node 0
    shares = sorted(abs(float(value)) for row, value in zip(coalesced.indices()[0], coalesced.values())
                    if int(row) == depletion_row)
    assert shares == pytest.approx([0.1, 0.9], abs=1e-6), f"shares must follow the parents' weights, got {shares}"
