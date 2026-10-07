"""Tests for the per-node cell-class weights of the cell-class propagation channels."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_weights import (ALL_CELLS_CLASS, extracellular_pool_nodes, node_cell_class_weights,
                                                                         relative_expression_weights)


def test_relative_weight_is_expression_over_the_top_class_plus_the_detection_threshold():
    table = pd.DataFrame({"neuron": [388.0, 2.0, 112.0], "glia": [6.0, 2.0, 5606.0]}, index=["TH", "FLAT", "MBP"])
    weights = relative_expression_weights(table)
    assert np.isclose(weights.loc["TH", "neuron"], 388 / 389)
    assert np.isclose(weights.loc["FLAT", "glia"], 2 / 3)
    assert weights.loc["MBP", "neuron"] < 0.03


def test_nodes_without_expression_get_one_and_reactions_follow_their_rule():
    nodes = pd.DataFrame({
        "node_id": ["g1", "g2", "r1", "r2", "m1", "m2"],
        "node_type": ["gene", "gene", "reaction", "reaction", "metabolite", "metabolite"],
        "ensembl_gene_id": ["E1", "E9", None, None, None, None],
        "gene_symbol": ["A", "Z", None, None, None, None],
        "gene_reaction_rule": [None, None, "E1 and E2", "", None, None],
        "compartment": [None, None, None, None, "c", "e"],
    })
    expression = pd.DataFrame({"neuron": [10.0, 1.0], "glia": [0.0, 9.0]}, index=pd.Index(["E1", "E2"], name="ensembl_gene_id"))
    weights = node_cell_class_weights(nodes, expression)
    assert np.isclose(weights.loc["g1", "neuron"], 10 / 11) and weights.loc["g1", "glia"] == 0.0
    assert (weights.loc["g2"] == 1.0).all()  # not in the table
    assert np.isclose(weights.loc["r1", "neuron"], 1 / 2) and weights.loc["r1", "glia"] == 0.0  # min(10, 1) = 1 and min(0, 9) = 0
    assert (weights.loc[["r2", "m1", "m2"]] == 1.0).all().all()  # no rule, metabolites
    assert (weights[ALL_CELLS_CLASS] == 1.0).all()
    assert extracellular_pool_nodes(nodes).tolist() == [False, False, False, False, False, True]
