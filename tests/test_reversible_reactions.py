"""Preregistration item 27: a perturbation that accumulates a reversible reaction's product drives it backwards.

Linearised mass action gives dv = k_f dS - k_r dP, so the product enters the reaction's rate at sign -1. These tests
run the encoder's own adjacency to a fixed point on a three-node chain and read the signs, which is the claim: a
rise in the product makes the rate deviation negative, the substrate rises and the product is pulled back down.
"""
from __future__ import annotations

import pytest
import torch

from mechanistic_pathway_learning.models.linear_response_encoder import (
    PRODUCT_OPPOSES_RELATION,
    LinearResponseEncoder,
)

SUBSTRATE, REACTION, PRODUCT = 0, 1, 2
RELATIONS = ["substrate_of", "product_of"]
GAIN = 0.5  # below the contraction bound, so the iteration settles


def settled_field(reversible: bool, seed_node: int, seed_value: float = 1.0) -> tuple[torch.Tensor, list[str]]:
    """S -> R -> P, with the reaction flagged reversible or not, seeded at one node and iterated to its fixed point
    under one gain for every relation."""
    edge_source = torch.tensor([SUBSTRATE, REACTION])
    edge_target = torch.tensor([REACTION, PRODUCT])
    edge_relation = torch.tensor([0, 1])
    edge_sign = torch.tensor([1.0, 1.0])
    reversible_reactions = torch.tensor([False, reversible, False])
    stacked, relation_names = LinearResponseEncoder.signed_stacked_adjacency(
        3, RELATIONS, edge_source, edge_target, edge_relation, edge_sign, reversible_reactions=reversible_reactions)
    operator = stacked.to_dense().view(len(relation_names), 3, 3).sum(dim=0) * GAIN
    seed = torch.zeros(3)
    seed[seed_node] = seed_value
    state = torch.zeros(3)
    for _ in range(500):
        state = operator @ state + seed
    return state, relation_names


def test_an_accumulating_product_runs_a_reversible_reaction_backwards():
    state, relation_names = settled_field(reversible=True, seed_node=PRODUCT)
    assert PRODUCT_OPPOSES_RELATION in relation_names
    assert state[REACTION] < 0, "the rate deviation should be negative: the reaction is pushed backwards"
    assert state[SUBSTRATE] > 0, "running backwards regenerates the substrate"


def test_an_irreversible_reaction_is_insensitive_to_its_product():
    """A thermodynamically irreversible step has k_r near 0, so a product rise reaches neither the rate nor the
    substrate. The relation still exists, with no edges, so the model's parameters do not depend on whether a given
    graph or slice happens to contain a reversible reaction."""
    state, relation_names = settled_field(reversible=False, seed_node=PRODUCT)
    assert PRODUCT_OPPOSES_RELATION in relation_names
    assert state[REACTION] == pytest.approx(0.0) and state[SUBSTRATE] == pytest.approx(0.0)


def test_the_forward_direction_is_unchanged():
    """A substrate rise still drives the reaction forwards and raises the product, reversible or not."""
    for reversible in (True, False):
        state, _ = settled_field(reversible=reversible, seed_node=SUBSTRATE)
        assert state[REACTION] > 0 and state[PRODUCT] > 0


def test_nothing_cancels_at_a_participant():
    """The defect this replaces: under the old mirror pair every participant received +1 and -1 from one reaction.
    Now the substrate receives only the depletion edge and the product only the production edge."""
    edge_source, edge_target = torch.tensor([SUBSTRATE, REACTION]), torch.tensor([REACTION, PRODUCT])
    stacked, relation_names = LinearResponseEncoder.signed_stacked_adjacency(
        3, RELATIONS, edge_source, edge_target, torch.tensor([0, 1]), torch.tensor([1.0, 1.0]),
        reversible_reactions=torch.tensor([False, True, False]))
    incoming = stacked.to_dense().view(len(relation_names), 3, 3)
    for participant in (SUBSTRATE, PRODUCT):
        signs = incoming[:, participant, REACTION]
        signs = signs[signs != 0]
        assert len(signs) == 1, f"node {participant} receives {len(signs)} edges from the reaction"


def test_a_graph_still_carrying_the_old_mirror_pair_is_refused():
    """On the old build the nominal direction is lost, so the reverse edge cannot be derived; running anyway would
    put product_opposes on top of the cancelling mirror."""
    edge_source = torch.tensor([SUBSTRATE, REACTION, REACTION, PRODUCT])
    edge_target = torch.tensor([REACTION, PRODUCT, SUBSTRATE, REACTION])
    edge_relation = torch.tensor([0, 1, 1, 0])  # the mirror: R -> S product_of and P -> R substrate_of
    with pytest.raises(ValueError, match="mirror pair"):
        LinearResponseEncoder.signed_stacked_adjacency(
            3, RELATIONS, edge_source, edge_target, edge_relation, torch.ones(4),
            reversible_reactions=torch.tensor([False, True, False]))
