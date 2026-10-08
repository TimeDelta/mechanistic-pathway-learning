"""Tests for the reaction expression that follows a rewired graph (mechanistic_pathway_learning/graph/rewired_expression_features.py)."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.brain_expression_weights import node_cell_class_weights
from mechanistic_pathway_learning.graph.rewired_expression_features import (
    rewired_gene_reaction_rules, rewired_reaction_cell_class_weights, substitute_rule_genes)


def toy_nodes() -> pd.DataFrame:
    return pd.DataFrame({"node_id": ["GENE:A", "GENE:B", "GENE:C", "PROTEIN:P", "R1", "R2", "R3"],
                         "node_type": ["gene", "gene", "gene", "protein_entity", "reaction", "reaction", "reaction"],
                         "ensembl_gene_id": ["EA", "EB", "EC", None, None, None, None],
                         "gene_reaction_rule": [None, None, None, None, "EA and EB", "EC", None]})


def test_substitution_keeps_the_rule_structure() -> None:
    assert substitute_rule_genes("(EA and EB) or EC", {"EA": "EX", "EC": "EY"}) == "( EX and EB ) or EY"


def test_rules_follow_the_rewired_catalysts_and_unchanged_reactions_keep_their_rule() -> None:
    nodes = toy_nodes()
    catalysis = np.ones(4, dtype=bool)
    before = np.array([[0, 1, 2, 3], [4, 4, 5, 6]])  # A, B -> R1; C -> R2; P -> R3
    after = np.array([[0, 1, 2, 3], [5, 4, 4, 6]])   # A -> R2; B, C -> R1; P -> R3 (every degree kept)
    rules, counts = rewired_gene_reaction_rules(nodes, before, after, catalysis)
    assert rules["R1"] == "EC and EB"  # B stays, A's place goes to C
    assert rules["R2"] == "EA"
    assert rules["R3"] == ""
    assert counts["rules_changed"] == 2 and counts["rule_genes_replaced"] == 2
    after_with_protein = np.array([[0, 1, 2, 3], [4, 4, 6, 5]])  # C -> R3; P -> R2
    rules, counts = rewired_gene_reaction_rules(nodes, before, after_with_protein, catalysis)
    assert rules["R2"] == "NO_EXPRESSION:PROTEIN:P" and counts["rule_genes_replaced_by_a_catalyst_without_expression"] == 1
    identity_rules, identity_counts = rewired_gene_reaction_rules(nodes, before, before, catalysis)
    assert identity_counts["rules_changed"] == 0 and identity_rules["R1"] == "EA and EB"


def test_weights_from_unchanged_rules_equal_the_built_weights_and_follow_a_rewritten_rule() -> None:
    nodes = toy_nodes()
    class_expression = pd.DataFrame({"neuron": [100.0, 1.0, 5.0], "glia": [1.0, 50.0, 5.0]}, index=["EA", "EB", "EC"])
    built = node_cell_class_weights(nodes, class_expression)
    rules = pd.Series(nodes.gene_reaction_rule.fillna("").to_numpy(), index=nodes.node_id)
    recomputed = rewired_reaction_cell_class_weights(nodes, rules, class_expression)
    pd.testing.assert_frame_equal(recomputed, built.loc[["R1", "R2", "R3"]])
    rewritten = rules.copy()
    rewritten["R2"] = "EA"
    assert np.allclose(rewired_reaction_cell_class_weights(nodes, rewritten, class_expression).loc["R2"].to_numpy(), built.loc["GENE:A"].to_numpy())
