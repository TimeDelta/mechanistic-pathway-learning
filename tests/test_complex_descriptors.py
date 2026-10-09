"""The Reactome entities' annotation vectors: union, intersection, the disjunction rule for sets, and the block they
enter (mechanistic_pathway_learning/graph/complex_descriptors.py, the user's decision of 9 October 2026)."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.complex_descriptors import (entity_annotation_rows, member_genes_of_entity,
                                                                    standardise_with_proteome_statistics)
from mechanistic_pathway_learning.graph.node_descriptors import assemble_node_descriptor_table

BLOCK_COLUMNS = {"pfam": ["pfam:PF1", "pfam:PF2"], "go_component": ["go_component:GO:1", "go_component:GO:2"]}
MEMBER_GENES = ["AAA", "BBB", "CCC"]


def annotation_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Three genes. AAA and BBB each carry one Pfam family and share component GO:1; CCC has no component annotation."""
    targets = pd.DataFrame([[1.0, 0.0, 1.0, 1.0],
                            [0.0, 1.0, 1.0, 0.0],
                            [1.0, 1.0, 0.0, 0.0]], index=pd.Index(MEMBER_GENES, name="gene"),
                           columns=BLOCK_COLUMNS["pfam"] + BLOCK_COLUMNS["go_component"])
    observed = pd.DataFrame({"pfam": [True, True, True], "go_component": [True, True, False]}, index=targets.index)
    return targets, observed


def test_a_complex_takes_the_union_of_the_domains_and_the_intersection_of_the_locations():
    targets, observed = annotation_tables()
    values, annotated = entity_annotation_rows({"RCTE_1": ["AAA", "BBB"]}, {"RCTE_1": "Complex"}, targets, observed, BLOCK_COLUMNS)
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["pfam"]]) == [1.0, 1.0], "a complex carries every domain its subunits bring"
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["go_component"]]) == [1.0, 0.0], "it sits only where both subunits sit"
    assert bool(annotated.loc["RCTE_1", "go_component"])


def test_an_unannotated_member_does_not_empty_the_intersection():
    """CCC has no component annotation, so it is absent from the intersection rather than a negative for every term."""
    targets, observed = annotation_tables()
    values, annotated = entity_annotation_rows({"RCTE_1": ["AAA", "CCC"]}, {"RCTE_1": "Complex"}, targets, observed, BLOCK_COLUMNS)
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["go_component"]]) == [1.0, 1.0]
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["pfam"]]) == [1.0, 1.0], "a missing Pfam family is a negative, so CCC counts there"


def test_an_entity_with_no_annotated_member_is_unannotated_in_that_block():
    targets, observed = annotation_tables()
    _, annotated = entity_annotation_rows({"RCTE_1": ["CCC"]}, {"RCTE_1": "Complex"}, targets, observed, BLOCK_COLUMNS)
    assert not bool(annotated.loc["RCTE_1", "go_component"])
    assert bool(annotated.loc["RCTE_1", "pfam"])


def test_a_defined_set_takes_the_mean_because_it_is_a_disjunction():
    targets, observed = annotation_tables()
    values, _ = entity_annotation_rows({"RCTE_1": ["AAA", "BBB"]}, {"RCTE_1": "DefinedSet"}, targets, observed, BLOCK_COLUMNS)
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["pfam"]]) == [0.5, 0.5]
    assert list(values.loc["RCTE_1", BLOCK_COLUMNS["go_component"]]) == [1.0, 0.5], "the mean holds in every block, locations included"


def test_a_one_member_entity_keeps_its_member_row():
    targets, observed = annotation_tables()
    values, _ = entity_annotation_rows({"RCTE_1": ["AAA"]}, {"RCTE_1": "Complex"}, targets, observed, BLOCK_COLUMNS)
    assert list(values.loc["RCTE_1"]) == list(targets.loc["AAA"])


def test_the_standardisation_uses_the_proteome_statistics_not_the_entities():
    targets, observed = annotation_tables()
    values, annotated = entity_annotation_rows({"RCTE_1": ["AAA"], "RCTE_2": ["BBB"]}, {"RCTE_1": "Complex", "RCTE_2": "Complex"},
                                               targets, observed, BLOCK_COLUMNS)
    standardised = standardise_with_proteome_statistics(values, annotated, targets, observed, BLOCK_COLUMNS)
    pfam_positions = [values.columns.get_loc(column) for column in BLOCK_COLUMNS["pfam"]]
    first_column = targets[BLOCK_COLUMNS["pfam"][0]].to_numpy()
    expected = (1.0 - first_column.mean()) / first_column.std() / np.sqrt(len(pfam_positions))
    assert np.isclose(standardised[0, pfam_positions[0]], expected)


def test_member_genes_come_from_the_member_of_edges():
    edges = pd.DataFrame({"source_id": ["GENE:AAA", "GENE:BBB", "GENE:CCC", "RCTE_1"],
                          "target_id": ["RCTE_1", "RCTE_1", "RCTE_9", "reaction_1"],
                          "relation_type": ["member_of", "member_of", "member_of", "substrate_of"]})
    assert member_genes_of_entity(edges, ["RCTE_1"]) == {"RCTE_1": ["AAA", "BBB"]}


def test_entities_share_the_protein_block_when_the_columns_match():
    nodes = pd.DataFrame({"node_id": ["GENE:AAA", "RCTE_1", "RCTE_2"], "node_type": ["gene", "protein_entity", "protein_entity"],
                          "gene_symbol": ["AAA", None, None]})
    protein_table = pd.DataFrame([[0.5, -0.5]], index=pd.Index(["AAA"], name="gene"), columns=["rrr_1", "rrr_2"])
    complex_table = pd.DataFrame([[1.5, 2.5]], index=pd.Index(["RCTE_1"], name="node_id"), columns=["rrr_1", "rrr_2"])
    table = assemble_node_descriptor_table(nodes, protein_table=protein_table, complex_table=complex_table)
    assert list(table.columns) == ["protein_rrr_1", "protein_rrr_2", "protein_has_protein_descriptors", "protein_has_complex_descriptors"]
    assert list(table.loc["RCTE_1", ["protein_rrr_1", "protein_rrr_2"]]) == [1.5, 2.5]
    assert table.loc["RCTE_1", "protein_has_complex_descriptors"] == 1.0
    assert table.loc["RCTE_2", "protein_has_complex_descriptors"] == 0.0, "an entity without a vector keeps a zero block"
    assert table.loc["GENE:AAA", "protein_has_complex_descriptors"] == 0.0
    assert list(table.loc["GENE:AAA", ["protein_rrr_1", "protein_rrr_2"]]) == [0.5, -0.5]


def test_entities_get_their_own_block_when_the_columns_differ():
    nodes = pd.DataFrame({"node_id": ["GENE:AAA", "RCTE_1"], "node_type": ["gene", "protein_entity"], "gene_symbol": ["AAA", None]})
    protein_table = pd.DataFrame([[0.5]], index=pd.Index(["AAA"], name="gene"), columns=["rrr_1"])
    complex_table = pd.DataFrame([[1.5, 2.5]], index=pd.Index(["RCTE_1"], name="node_id"), columns=["annotation_1", "annotation_2"])
    table = assemble_node_descriptor_table(nodes, protein_table=protein_table, complex_table=complex_table)
    assert "complex_annotation_1" in table.columns and "complex_has_complex_descriptors" in table.columns
    assert table.loc["GENE:AAA", "complex_annotation_1"] == 0.0
