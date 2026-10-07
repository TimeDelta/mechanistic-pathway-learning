import math

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.brain_expression_descriptors import (
    CELL_CLASS_OF_CLUSTER_TYPE, brain_expression_blocks, evaluate_gene_rule, hpa_cell_class_table, hpa_region_table)


def test_gene_rule_takes_minimum_over_and_and_maximum_over_or():
    values = {"A": 1.0, "B": 5.0, "C": 3.0}
    assert evaluate_gene_rule("A or B", values) == 5.0
    assert evaluate_gene_rule("A and B", values) == 1.0
    assert evaluate_gene_rule("A or B and C", values) == 3.0  # and binds tighter: A or (B and C)
    assert evaluate_gene_rule("(A or B) and C", values) == 3.0
    assert evaluate_gene_rule("(A and B) or (B and C)", values) == 3.0


def test_gene_rule_skips_unknown_genes_and_is_missing_when_none_known():
    values = {"A": 2.0}
    assert evaluate_gene_rule("A and UNKNOWN", values) == 2.0
    assert evaluate_gene_rule("UNKNOWN or A", values) == 2.0
    assert math.isnan(evaluate_gene_rule("X or Y", values))
    assert math.isnan(evaluate_gene_rule("", values))


def test_gene_rule_rejects_unbalanced_parentheses():
    with pytest.raises(ValueError):
        evaluate_gene_rule("(A or B", {"A": 1.0})


def test_cell_class_table_takes_the_largest_value_in_each_class_and_rejects_unknown_types():
    clusters = pd.DataFrame({"Gene": ["G1", "G1", "G1"], "Gene name": ["X"] * 3,
                             "Cluster type": ["CGE interneuron", "MGE interneuron", "astrocyte"], "nCPM": [2.0, 7.0, 1.0]})
    table = hpa_cell_class_table(clusters)
    assert table.loc["G1", "class_cortical_interneuron"] == 7.0
    assert table.loc["G1", "class_astrocyte"] == 1.0
    with pytest.raises(ValueError):
        hpa_cell_class_table(clusters.assign(**{"Cluster type": ["no such type"] * 3}))
    assert len(set(CELL_CLASS_OF_CLUSTER_TYPE.values())) == 10


def test_blocks_fill_gene_and_reaction_rows_only():
    nodes = pd.DataFrame({
        "node_id": ["GENE:A", "GENE:B", "R1", "R2", "M1"],
        "node_type": ["gene", "gene", "reaction", "reaction", "metabolite"],
        "ensembl_gene_id": ["EA", "EB", None, None, None],
        "gene_reaction_rule": [None, None, "EA and EB", "", None],
    })
    regions = pd.DataFrame({"Gene": ["EA", "EB"], "Gene name": ["A", "B"], "Brain region": ["amygdala"] * 2,
                            "TPM": [1.0, 1.0], "pTPM": [1.0, 1.0], "nTPM": [10.0, 100.0]})
    expression = hpa_region_table(regions)
    blocks = brain_expression_blocks(nodes, expression)
    assert list(blocks.index) == list(nodes.node_id)
    assert (blocks.loc["M1"] == 0.0).all()
    assert blocks.loc["GENE:A", "gene_brain_has_expression"] == 1.0
    assert blocks.loc["R1", "reaction_brain_has_expression"] == 1.0
    assert blocks.loc["R2", "reaction_brain_has_expression"] == 0.0
    assert blocks.loc["GENE:A", "gene_brain_region_amygdala"] < blocks.loc["GENE:B", "gene_brain_region_amygdala"]
    assert (blocks.filter(like="reaction_brain_").loc[["GENE:A", "GENE:B"]] == 0.0).all().all()
    assert np.isfinite(blocks.to_numpy()).all()


def test_gene_nodes_without_an_ensembl_id_are_found_by_unambiguous_symbol():
    from mechanistic_pathway_learning.graph.brain_expression_descriptors import unambiguous_symbol_to_ensembl
    mapping = unambiguous_symbol_to_ensembl(pd.Series(["A", "B", "B"]), pd.Series(["EA", "EB1", "EB2"]))
    assert mapping == {"A": "EA"}
    nodes = pd.DataFrame({
        "node_id": ["GENE:A", "GENE:B"], "node_type": ["gene", "gene"], "gene_symbol": ["A", "B"],
        "ensembl_gene_id": [None, None], "gene_reaction_rule": [None, None],
    })
    regions = pd.DataFrame({"Gene": ["EA"], "Gene name": ["A"], "Brain region": ["amygdala"], "TPM": [1.0], "pTPM": [1.0], "nTPM": [10.0]})
    blocks = brain_expression_blocks(nodes, hpa_region_table(regions), mapping)
    assert blocks.loc["GENE:A", "gene_brain_has_expression"] == 1.0
    assert blocks.loc["GENE:B", "gene_brain_has_expression"] == 0.0
