"""The plasma substrate binders as carriers: the UniProt ligand parser, the primary-gene filter, and the carriage edges
(mechanistic_pathway_learning/graph/plasma_binding.py, the user's decision of 9 October 2026)."""
import csv
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.plasma_binding import (BINDS_RELATION, CARRIAGE_SIGN, CURATED_EVIDENCE, UNIPROT_EVIDENCE,
                                                               UNIPROT_SOURCE, add_plasma_carriage, curated_carriage_rows,
                                                               ligands_of_binding_features, uniprot_carriage_rows)

CURATED_TABLE = Path("docs/curated_plasma_carriage.csv")
ALBUMIN_FEATURES = ('BINDING 264; /ligand="(4Z,15Z)-bilirubin IXalpha"; /ligand_id="ChEBI:CHEBI:57977"; '
                    '/evidence="ECO:0000269|PubMed:656055"; BINDING 268; /ligand="Ca(2+)"; /ligand_id="ChEBI:CHEBI:29108"; '
                    '/ligand_label="1"; /evidence="ECO:0000250|UniProtKB:P02769"')
BARE_LIGAND_THEN_CHEBI = ('BINDING 100; /ligand="substrate"; /evidence="ECO:0000250"; '
                          'BINDING 200; /ligand="Zn(2+)"; /ligand_id="ChEBI:CHEBI:29105"')


def graph_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Two binder genes, a cargo with an extracellular and a cytosolic copy, and one cargo edge already present."""
    nodes = pd.DataFrame({
        "node_id": ["GENE:ALB", "GENE:SERPINA6", "MAM01615e", "MAM01615c", "MAM01396e", "GENE:OTHER"],
        "node_type": ["gene", "gene", "metabolite", "metabolite", "metabolite", "gene"],
        "display_name": ["ALB", "SERPINA6", "cortisol", "cortisol", "BIL", "OTHER"],
        "compartment": [None, None, "e", "c", "e", None],
        "base_metabolite_id": [None, None, "MAM01615", "MAM01615", "MAM01396", None],
        "degree": [0, 0, 0, 0, 0, 0]})
    edges = pd.DataFrame({"source_id": ["MAM01615e"], "target_id": ["GENE:SERPINA6"], "relation_type": [BINDS_RELATION],
                          "sign": [0.0], "evidence_source": ["OmniPath small_molecule_protein"]})
    return nodes, edges


def carriage_frame(rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["binder_gene", "cargo_name", "cargo_chebi", "source"])


def test_a_ligand_without_a_chebi_identifier_is_not_paired_with_the_next_feature():
    assert ligands_of_binding_features(BARE_LIGAND_THEN_CHEBI) == [("Zn(2+)", "CHEBI:29105")]


def test_every_ligand_of_a_feature_field_is_read_with_its_chebi_identifier():
    assert ligands_of_binding_features(ALBUMIN_FEATURES) == [("(4Z,15Z)-bilirubin IXalpha", "CHEBI:57977"), ("Ca(2+)", "CHEBI:29108")]


def test_the_primary_gene_name_decides_which_entries_are_binders():
    """FBF1 lists ALB among its gene names; only the primary name may match, or a cilium protein gets albumin's cargo."""
    table = pd.DataFrame({"Gene Names (primary)": ["ALB", "FBF1"], "Gene Names": ["ALB", "FBF1 ALB KIAA1863"],
                          "Binding site": [ALBUMIN_FEATURES, ALBUMIN_FEATURES]})
    rows = uniprot_carriage_rows(table, ("ALB", "SERPINA6"))
    assert set(rows.binder_gene) == {"ALB"}
    assert set(rows.source) == {UNIPROT_SOURCE}


def test_the_carriage_edge_runs_from_the_extracellular_cargo_to_the_binder_gene():
    nodes, edges = graph_frames()
    new_nodes, new_edges, relations, summary = add_plasma_carriage(
        nodes, edges, [BINDS_RELATION], carriage_frame([("ALB", "bilirubin", "CHEBI:16990", "pubmed")]),
        lambda chebi: {"MAM01396"})
    added = new_edges.iloc[len(edges):]
    assert list(added.source_id) == ["MAM01396e"] and list(added.target_id) == ["GENE:ALB"], "metabolite -> gene, the graph's own direction"
    assert list(added.relation_type) == [BINDS_RELATION] and list(added.sign) == [CARRIAGE_SIGN]
    assert list(added.evidence_source) == [CURATED_EVIDENCE]
    assert summary["edges_added"] == 1 and relations == [BINDS_RELATION], "binds already exists, so no relation is added"


def test_only_the_extracellular_copy_of_a_cargo_is_joined():
    nodes, edges = graph_frames()
    _, new_edges, _, summary = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                                   carriage_frame([("ALB", "cortisol", "CHEBI:17650", "pubmed")]),
                                                   lambda chebi: {"MAM01615"})
    assert list(new_edges.iloc[len(edges):].source_id) == ["MAM01615e"], "plasma is the extracellular compartment"
    assert summary["edges_added"] == 1


def test_an_edge_the_graph_already_holds_is_not_duplicated():
    nodes, edges = graph_frames()
    _, new_edges, _, summary = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                                   carriage_frame([("SERPINA6", "cortisol", "CHEBI:17650", UNIPROT_SOURCE)]),
                                                   lambda chebi: {"MAM01615"})
    assert len(new_edges) == len(edges) and summary["edges_added"] == 0
    assert summary["edges_already_in_the_graph"] == ["MAM01615e -> GENE:SERPINA6"]


def test_a_row_that_maps_to_nothing_is_reported_rather_than_invented():
    nodes, edges = graph_frames()
    unmapped = carriage_frame([("ALB", "heparin", "CHEBI:28304", UNIPROT_SOURCE),
                               ("MISSING", "cortisol", "CHEBI:17650", "pubmed")])
    new_nodes, new_edges, _, summary = add_plasma_carriage(nodes, edges, [BINDS_RELATION], unmapped, lambda chebi: set())
    assert len(new_edges) == len(edges) and len(new_nodes) == len(nodes), "no node and no edge is invented"
    assert summary["binder_genes_absent_from_graph"] == ["MISSING"]
    assert summary["cargo_not_mapped_to_human_gem"] == ["ALB heparin (CHEBI:28304)"]


def test_a_cargo_without_an_extracellular_copy_is_reported():
    nodes, edges = graph_frames()
    nodes = nodes[nodes.node_id != "MAM01396e"].reset_index(drop=True)
    _, _, _, summary = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                           carriage_frame([("ALB", "bilirubin", "CHEBI:16990", "pubmed")]), lambda chebi: {"MAM01396"})
    assert summary["cargo_without_an_extracellular_node"] == ["ALB bilirubin (CHEBI:16990 -> MAM01396)"]


def test_the_degrees_count_the_new_edges():
    nodes, edges = graph_frames()
    new_nodes, _, _, _ = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                             carriage_frame([("ALB", "bilirubin", "CHEBI:16990", "pubmed")]), lambda chebi: {"MAM01396"})
    degree_of = dict(zip(new_nodes.node_id, new_nodes.degree))
    assert degree_of["GENE:ALB"] == 1 and degree_of["MAM01396e"] == 1
    assert degree_of["GENE:SERPINA6"] == 1 and degree_of["GENE:OTHER"] == 0


def test_the_uniprot_evidence_source_says_where_an_edge_came_from():
    nodes, edges = graph_frames()
    _, new_edges, _, _ = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                             carriage_frame([("ALB", "bilirubin", "CHEBI:16990", UNIPROT_SOURCE)]),
                                             lambda chebi: {"MAM01396"})
    assert list(new_edges.iloc[len(edges):].evidence_source) == [UNIPROT_EVIDENCE]


def test_every_curated_row_carries_its_source_and_a_quote():
    rows = list(csv.DictReader(CURATED_TABLE.open()))
    assert rows, "the curated carriage table must not be empty"
    for row in rows:
        assert row["binder_gene"] and row["cargo_name"] and row["cargo_chebi"].startswith("CHEBI:")
        assert row["source"] in {"pubmed", "uniprot_function"}, row["source"]
        assert len(row["supporting_quote"]) > 40, row
        if row["source"] == "pubmed":
            assert row["pmid"].isdigit(), row
    frame = curated_carriage_rows(pd.read_csv(CURATED_TABLE))
    assert list(frame.columns) == ["binder_gene", "cargo_name", "cargo_chebi", "source"]
    assert "ORM1" in set(frame.binder_gene) and "ORM2" in set(frame.binder_gene), "orosomucoid is modelled (the user's request)"


def test_on_a_split_graph_the_carriage_sits_on_the_protein_node():
    """Binding is the protein's property; on a gene/protein split graph the gene node holds only expression."""
    nodes, edges = graph_frames()
    nodes = pd.concat([nodes, pd.DataFrame([{"node_id": "PROTEIN:ALB", "node_type": "protein", "display_name": "ALB",
                                             "compartment": None, "base_metabolite_id": None, "degree": 0}])], ignore_index=True)
    _, new_edges, _, summary = add_plasma_carriage(nodes, edges, [BINDS_RELATION],
                                                   carriage_frame([("ALB", "bilirubin", "CHEBI:16990", "pubmed")]),
                                                   lambda chebi: {"MAM01396"}, {"ALB": ["PROTEIN:ALB"]})
    assert list(new_edges.iloc[len(edges):].target_id) == ["PROTEIN:ALB"]
    assert summary["binder_nodes"] == ["PROTEIN:ALB"]


def test_the_split_binder_nodes_come_from_the_gene_to_protein_map(tmp_path):
    from experiments.build_plasma_binder_variant import binder_nodes_on_split
    assert binder_nodes_on_split(tmp_path, ("ALB",)) is None, "a merged graph has no map, so the gene node is the binder"
    pd.DataFrame({"gene_node_id": ["GENE:ALB", "GENE:OTHER"], "protein_node_id": ["PROTEIN:ALB", "PROTEIN:OTHER"]}
                 ).to_parquet(tmp_path / "gene_to_protein.parquet")
    assert binder_nodes_on_split(tmp_path, ("ALB",)) == {"ALB": ["PROTEIN:ALB"]}


def test_the_report_names_the_endpoint_the_build_actually_used():
    """The generated document said "gene node" whatever graph it described, which is wrong for a split graph."""
    from experiments.build_plasma_binder_variant import carriage_report
    summary = {"edges_added": 27, "carriage_on": "protein nodes (gene/protein split graph)", "cargo_by_binder": {"ALB": ["bilirubin"]},
               "edges_by_binder": {"ALB": 1}, "edges_by_source": {"pubmed": 1}, "carriage_rows_read": 1,
               "binder_genes_absent_from_graph": [], "cargo_not_mapped_to_human_gem": [],
               "cargo_without_an_extracellular_node": [], "edges_already_in_the_graph": []}
    report = carriage_report(summary, Path("graph_split"), Path("graph_split_binders"), pd.DataFrame([{"binder_gene": "ALB"}]),
                             [Path("graph_other_binders")])
    assert "protein nodes (gene/protein split graph)" in report
    assert "binder's gene node" not in report
    assert "graph_other_binders" in report, "the other variants the same run built are named"
