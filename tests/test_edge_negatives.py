"""IntAct's curated non-interacting pairs against the graph's protein edges (mechanistic_pathway_learning/graph/edge_negatives.py)."""
from __future__ import annotations

import pandas as pd

from mechanistic_pathway_learning.graph.edge_negatives import NegativePair, audit, contradicted_edges, read_human_negative_pairs, symbol_of_node

HEADER = ["#ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A", "Alt. ID(s) interactor B", "Alias(es) interactor A", "Alias(es) interactor B",
          "Interaction detection method(s)", "Publication 1st author(s)", "Publication Identifier(s)", "Taxid interactor A", "Taxid interactor B",
          "Interaction type(s)", "Source database(s)", "Interaction identifier(s)", "Confidence value(s)"]


def mitab_row(acc_a: str, acc_b: str, gene_a: str, gene_b: str, taxid_b: str = "taxid:9606(human)|taxid:9606(Homo sapiens)") -> list[str]:
    return [f"uniprotkb:{acc_a}", f"uniprotkb:{acc_b}", "-", "-", f"psi-mi:x_human(display_long)|uniprotkb:{gene_a}(gene name)",
            f"psi-mi:y_human(display_long)|uniprotkb:{gene_b}(gene name)", 'psi-mi:"MI:0018"(two hybrid)', "Doe J (2020)", "pubmed:12345",
            "taxid:9606(human)|taxid:9606(Homo sapiens)", taxid_b, 'psi-mi:"MI:0915"(physical association)', 'psi-mi:"MI:0469"(IntAct)', "intact:EBI-1", "-"]


def write_mitab(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as mitab:
        mitab.write("\t".join(HEADER) + "\n")
        for row in rows:
            mitab.write("\t".join(row) + "\n")


def test_human_pairs_are_read_with_isoforms_dropped_and_other_species_left_out(tmp_path):
    path = tmp_path / "intact_negative.txt"
    write_mitab(path, [mitab_row("P04637-2", "Q00987", "TP53", "MDM2"), mitab_row("P12345", "P23456", "AAA", "BBB", taxid_b="taxid:10090(mouse)|taxid:10090(Mus musculus)")])
    pairs = read_human_negative_pairs(path)
    assert [(pair.accession_a, pair.accession_b, pair.gene_a, pair.gene_b) for pair in pairs] == [("P04637", "Q00987", "TP53", "MDM2")]
    assert pairs[0].publications == "pubmed:12345" and "physical association" in pairs[0].interaction_type


def test_contradicted_edges_are_found_in_either_direction_on_gene_or_protein_nodes():
    pairs = [NegativePair("P1", "P2", "TP53", "MDM2", "physical association", "two hybrid", "pubmed:1")]
    nodes = pd.DataFrame({"node_id": ["GENE:TP53", "PROTEIN:MDM2", "GENE:BRCA1", "MAM1c"], "node_type": ["gene", "protein", "gene", "metabolite"],
                          "gene_symbol": ["TP53", "MDM2", "BRCA1", None]})
    edges = pd.DataFrame({"source_id": ["PROTEIN:MDM2", "GENE:BRCA1", "GENE:TP53"], "target_id": ["GENE:TP53", "GENE:TP53", "MAM1c"],
                          "relation_type": ["inhibits", "binds", "substrate_of"], "sign": [-1.0, 0.0, 1.0], "evidence_source": ["OmniPath", "OmniPath", "Human-GEM"]})
    contradicted = contradicted_edges(pairs, edges, nodes)
    assert contradicted.source_id.tolist() == ["PROTEIN:MDM2"] and contradicted.intact_publications.iloc[0] == "pubmed:1"
    assert symbol_of_node(nodes) == {"GENE:TP53": "TP53", "PROTEIN:MDM2": "MDM2", "GENE:BRCA1": "BRCA1"}


def test_the_audit_judges_a_source_on_the_pairs_the_anchor_set_names():
    pairs = [NegativePair("P1", "P2", "TP53", "MDM2", "physical association", "two hybrid", "pubmed:1"),
             NegativePair("P3", "P4", "AKT1", "BRCA1", "direct interaction", "pull down", "pubmed:2")]
    nodes = pd.DataFrame({"node_id": ["GENE:TP53", "GENE:MDM2", "GENE:BRCA1", "GENE:AKT1", "GENE:EGFR"], "node_type": ["gene"] * 5,
                          "gene_symbol": ["TP53", "MDM2", "BRCA1", "AKT1", "EGFR"]})
    edges = pd.DataFrame({"source_id": ["GENE:MDM2", "GENE:TP53", "GENE:EGFR", "GENE:AKT1"], "target_id": ["GENE:TP53", "GENE:BRCA1", "GENE:AKT1", "GENE:BRCA1"],
                          "relation_type": ["inhibits", "binds", "activates", "binds"], "sign": [-1.0, 0.0, 1.0, 0.0],
                          "evidence_source": ["OmniPath", "OmniPath", "OmniPath", "CollecTRI"]})
    result = audit(pairs, edges, nodes)
    # OmniPath: three edges, two join proteins the set names (MDM2-TP53, TP53-BRCA1); one is a negative pair (MDM2-TP53); EGFR is unnamed
    assert result["by_source"]["OmniPath"] == {"edges_on_named_proteins": 2, "edges_on_a_negative_pair": 1}
    assert result["by_source"]["CollecTRI"] == {"edges_on_named_proteins": 1, "edges_on_a_negative_pair": 1}
    assert result["human_negative_pairs"] == 2 and result["proteins_named"] == 4 and result["edges_on_a_negative_pair"] == 2
