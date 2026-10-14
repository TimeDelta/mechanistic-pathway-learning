"""Tests for the per-node cell-class weights of the cell-class propagation channels."""
import numpy as np
import pandas as pd

import pytest

from mechanistic_pathway_learning.graph.brain_expression_weights import (ALL_CELLS_CLASS, complex_limiting_member_expression,
                                                                         extracellular_pool_nodes, node_cell_class_weights,
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


# --- Decision 4.2: a protein complex takes its limiting member's expression, per class ---

def complex_fixture() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Two subunits with opposite class specificities, a complex of them, and a complex with no expressed member."""
    nodes = pd.DataFrame({
        "node_id": ["gene:ABUNDANT", "gene:SCARCE", "gene:UNKNOWN", "complex:BOTH", "complex:UNKNOWN_ONLY",
                    "reaction:DOWNSTREAM"],
        "node_type": ["gene", "gene", "gene", "protein_entity", "protein_entity", "reaction"],
        "ensembl_gene_id": ["ABUNDANT", "SCARCE", None, None, None, None],
        "gene_symbol": ["ABUNDANT", "SCARCE", "UNKNOWN", None, None, None],
        "gene_reaction_rule": [None, None, None, None, None, ""],
    })
    edges = pd.DataFrame({
        "source_id": ["gene:ABUNDANT", "gene:SCARCE", "gene:UNKNOWN", "complex:BOTH"],
        "target_id": ["complex:BOTH", "complex:BOTH", "complex:UNKNOWN_ONLY", "reaction:DOWNSTREAM"],
        "relation_type": ["member_of", "member_of", "member_of", "catalyzed_by"],
    })
    class_expression = pd.DataFrame({"neuron": [1000.0, 10.0], "astrocyte": [10.0, 1000.0]},
                                    index=pd.Index(["ABUNDANT", "SCARCE"], name="ensembl_gene_id"))
    return nodes, edges, class_expression


def test_a_complex_takes_the_minimum_raw_expression_of_its_members_per_class():
    nodes, edges, class_expression = complex_fixture()
    limiting = complex_limiting_member_expression(nodes, edges, class_expression)
    assert list(limiting.index) == ["complex:BOTH"], "a complex with no expressed member must be absent, not zero"
    assert limiting.loc["complex:BOTH", "neuron"] == pytest.approx(10.0)
    assert limiting.loc["complex:BOTH", "astrocyte"] == pytest.approx(10.0)


def test_the_limiting_member_is_chosen_within_each_class_and_can_differ_between_them():
    """"per type": the scarce subunit in the neuron column is not the scarce one in the astrocyte column."""
    nodes, edges, class_expression = complex_fixture()
    limiting = complex_limiting_member_expression(nodes, edges, class_expression)
    assert limiting.loc["complex:BOTH", "neuron"] == class_expression.loc["SCARCE", "neuron"]
    assert limiting.loc["complex:BOTH", "astrocyte"] == class_expression.loc["ABUNDANT", "astrocyte"]


def test_the_minimum_is_on_raw_expression_and_not_on_the_members_relative_weights():
    """The registered order. Taking it the other way round drives a complex of differently specific subunits
    towards zero in every class, which is what docs/complex_limiting_member_weight.md measures."""
    nodes, edges, class_expression = complex_fixture()
    registered = relative_expression_weights(
        complex_limiting_member_expression(nodes, edges, class_expression))
    member_weights = relative_expression_weights(class_expression)
    the_other_order = member_weights.min(axis=0)
    assert registered.loc["complex:BOTH"].max() == pytest.approx(10.0 / 11.0)
    assert the_other_order.max() < 0.02
    assert registered.loc["complex:BOTH"].max() > 20 * the_other_order.max()


def test_a_complex_with_no_expressed_member_keeps_weight_one_in_every_class():
    nodes, edges, class_expression = complex_fixture()
    weights = node_cell_class_weights(nodes, class_expression, member_edges=edges)
    assert (weights.loc["complex:UNKNOWN_ONLY", ["neuron", "astrocyte"]] == 1.0).all()


def test_omitting_the_member_edges_leaves_every_complex_flat():
    """The behaviour before decision 4.2, kept available so the change can be measured against it."""
    nodes, edges, class_expression = complex_fixture()
    without = node_cell_class_weights(nodes, class_expression)
    with_members = node_cell_class_weights(nodes, class_expression, member_edges=edges)
    assert (without.loc["complex:BOTH", ["neuron", "astrocyte"]] == 1.0).all()
    assert (with_members.loc["complex:BOTH", ["neuron", "astrocyte"]] < 1.0).all()
    for gene in ("gene:ABUNDANT", "gene:SCARCE", "gene:UNKNOWN"):
        assert without.loc[gene].equals(with_members.loc[gene]), "gene weights must not move"


# --- Decision 4.2b: the rewiring control recomputes the complex gate from the rewired membership ---

def test_the_rewired_gate_follows_the_rewired_membership():
    """The point of 4.2b: a complex's gate must come from whichever genes the rewired graph says belong to it."""
    from mechanistic_pathway_learning.graph.rewired_expression_features import rewired_complex_cell_class_weights
    nodes, edges, class_expression = complex_fixture()
    node_index = {node_id: position for position, node_id in enumerate(nodes.node_id)}
    # the real membership: ABUNDANT and SCARCE both in complex:BOTH
    real = np.array([[node_index["gene:ABUNDANT"], node_index["gene:SCARCE"]],
                     [node_index["complex:BOTH"], node_index["complex:BOTH"]]])
    # a rewiring that moves SCARCE out, leaving the complex with the abundant subunit only
    rewired = np.array([[node_index["gene:ABUNDANT"], node_index["gene:SCARCE"]],
                        [node_index["complex:BOTH"], node_index["complex:UNKNOWN_ONLY"]]])
    membership = np.array([True, True])
    before, _ = rewired_complex_cell_class_weights(nodes, real, membership, class_expression,
                                                   include_all_cells_class=False)
    after, counts = rewired_complex_cell_class_weights(nodes, rewired, membership, class_expression,
                                                       include_all_cells_class=False)
    assert before.loc["complex:BOTH", "neuron"] == pytest.approx(10.0 / 11.0), "min of the two subunits"
    assert after.loc["complex:BOTH", "neuron"] == pytest.approx(1000.0 / 1001.0), "the abundant subunit alone"
    assert counts["member_of_edges_rewired"] == 2


def test_every_complex_gets_a_row_so_a_real_gate_cannot_survive_the_rewiring():
    """A complex whose rewired members carry no expression must fall back to 1, not keep its real gate."""
    from mechanistic_pathway_learning.graph.rewired_expression_features import rewired_complex_cell_class_weights
    nodes, edges, class_expression = complex_fixture()
    node_index = {node_id: position for position, node_id in enumerate(nodes.node_id)}
    # the rewiring gives complex:BOTH only the gene with no expression
    rewired = np.array([[node_index["gene:UNKNOWN"]], [node_index["complex:BOTH"]]])
    weights, counts = rewired_complex_cell_class_weights(nodes, rewired, np.array([True]), class_expression,
                                                         include_all_cells_class=False)
    assert set(weights.index) == {"complex:BOTH", "complex:UNKNOWN_ONLY"}, "every complex needs a row"
    assert (weights.loc["complex:BOTH"] == 1.0).all()
    assert counts["complexes_with_a_rewired_gate"] == 0
    assert counts["complexes_left_flat_by_the_rewiring"] == 2


def test_the_rewired_gate_carries_the_all_cells_column_when_the_stored_table_does():
    from mechanistic_pathway_learning.graph.rewired_expression_features import rewired_complex_cell_class_weights
    nodes, edges, class_expression = complex_fixture()
    node_index = {node_id: position for position, node_id in enumerate(nodes.node_id)}
    rewired = np.array([[node_index["gene:ABUNDANT"]], [node_index["complex:BOTH"]]])
    weights, _ = rewired_complex_cell_class_weights(nodes, rewired, np.array([True]), class_expression,
                                                     include_all_cells_class=True)
    assert ALL_CELLS_CLASS in weights.columns
    assert (weights[ALL_CELLS_CLASS] == 1.0).all()
