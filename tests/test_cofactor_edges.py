"""Carrier edges: the name list, the recurring-pair rule (per edge) and the encoder's carrier relations."""
import numpy as np
import torch

from mechanistic_pathway_learning.graph.cofactor_edges import cofactor_edge_mask
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

RELATION_TYPES = ["substrate_of", "product_of", "catalyzed_by"]


def test_carrier_list_and_recurring_pairs_mark_edges_but_not_main_pairs() -> None:
    # nodes: 0 phenylalanine, 1 tyrosine, 2 tetrahydrobiopterin, 3 dihydrobiopterin, 4 reaction PAH,
    # 5 glutamate, 6 2-oxoglutarate, 7..9 three transamination reactions with substrates 10..12 and products 13..15
    names = ["phenylalanine", "tyrosine", "tetrahydrobiopterin", "dihydrobiopterin", "PAH", "glutamate", "AKG",
             "TA1", "TA2", "TA3", "x1", "x2", "x3", "y1", "y2", "y3"]
    base = ["phe", "tyr", "bh4", "bh2", "", "glu", "akg", "", "", "", "x1", "x2", "x3", "y1", "y2", "y3"]
    edges = [(0, 4, "substrate_of"), (2, 4, "substrate_of"), (4, 1, "product_of"), (4, 3, "product_of")]
    for reaction, substrate, product in ((7, 10, 13), (8, 11, 14), (9, 12, 15)):
        edges += [(5, reaction, "substrate_of"), (substrate, reaction, "substrate_of"), (reaction, 6, "product_of"), (reaction, product, "product_of")]
    source = np.array([edge[0] for edge in edges])
    target = np.array([edge[1] for edge in edges])
    relation = np.array([RELATION_TYPES.index(edge[2]) for edge in edges])
    mask = cofactor_edge_mask(source, target, relation, RELATION_TYPES, np.array(base), np.array(names), np.zeros(len(names), dtype=bool), minimum_pair_reactions=3)
    marked = {(int(s), int(t)) for s, t, m in zip(source, target, mask) if m}
    assert (2, 4) in marked and (4, 3) in marked  # biopterins are on the carrier list
    assert (0, 4) not in marked and (4, 1) not in marked  # phenylalanine -> tyrosine is the main pair
    assert {(5, 7), (5, 8), (5, 9), (7, 6), (8, 6), (9, 6)} <= marked  # glutamate -> 2-oxoglutarate recurs in three reactions
    assert (10, 7) not in marked and (7, 13) not in marked


def test_encoder_moves_carrier_edges_to_their_own_relations_without_a_depletion_edge() -> None:
    edges = [(2, 1, "substrate_of", 1.0), (0, 1, "catalyzed_by", 1.0), (1, 3, "product_of", 1.0), (4, 1, "substrate_of", 1.0)]
    cofactor = torch.tensor([False, False, False, True])
    encoder = LinearResponseEncoder(
        5, RELATION_TYPES, edge_source=torch.tensor([e[0] for e in edges]), edge_target=torch.tensor([e[1] for e in edges]),
        edge_relation=torch.tensor([RELATION_TYPES.index(e[2]) for e in edges]), edge_sign=torch.tensor([e[3] for e in edges]),
        node_features=torch.eye(5), node_state_dim=2, cofactor_edges=cofactor,
    )
    assert encoder.relation_names == RELATION_TYPES + ["cosubstrate_of", "coproduct_of", "depletes_substrate"]
    dense = encoder.stacked_adjacency.to_dense().view(len(encoder.relation_names), 5, 5)
    names = encoder.relation_names
    assert dense[names.index("cosubstrate_of"), 1, 4] > 0 and dense[names.index("substrate_of"), 1, 4] == 0
    assert dense[names.index("depletes_substrate"), 2, 1] < 0  # the main substrate is depleted
    assert torch.all(dense[:, 4, 1] == 0)  # a recycled carrier is not: no edge from the reaction back to it
