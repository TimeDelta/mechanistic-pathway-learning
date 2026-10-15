"""Decision 9.5: BRENDA's enzyme regulation enters the graph unweighted, one edge per pair."""
from __future__ import annotations

import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.brenda_enzyme_regulation import (
    SIGN_OF_RELATION,
    add_brenda_enzyme_regulation,
    enzyme_node_of_accession,
    representative_metabolite_node,
)

NODE_COLUMNS = ["node_id", "node_type", "display_name", "compartment", "base_metabolite_id", "degree", "uniprot_entry"]
EDGE_COLUMNS = ["source_id", "target_id", "relation_type", "sign", "evidence_source"]


def graph_with_two_compartments() -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.DataFrame(
        [["MAM00001c", "metabolite", "citrate", "c", "MAM00001", 0, None],
         ["MAM00001m", "metabolite", "citrate", "m", "MAM00001", 0, None],
         ["MAM00002m", "metabolite", "malate", "m", "MAM00002", 0, None],
         ["PROTEIN:PFKL", "protein", "PFKL", None, None, 0, "P17858"],
         ["PROTEIN:IMPA1", "protein", "IMPA1", None, None, 0, "P29218"]],
        columns=NODE_COLUMNS)
    edges = pd.DataFrame([["MAM00002m", "PROTEIN:PFKL", "activates", 1.0, "OmniPath"]], columns=EDGE_COLUMNS)
    return nodes, edges


def test_the_cytosolic_copy_is_the_representative_and_a_base_without_one_falls_back():
    nodes, _ = graph_with_two_compartments()
    representative = representative_metabolite_node(nodes)
    assert representative["MAM00001"] == "MAM00001c"  # cytosol preferred over mitochondrion
    assert representative["MAM00002"] == "MAM00002m"  # no cytosolic copy, so the only one present


def test_one_edge_per_pair_rather_than_one_per_compartment_copy():
    """Proteins carry no compartment, so emitting from both copies of citrate would double PFKL's in-degree."""
    nodes, edges = graph_with_two_compartments()
    _, new_edges, _, summary = add_brenda_enzyme_regulation(
        nodes, edges, ["activates"], [{"base_metabolite_id": "MAM00001", "uniprot_accession": "P17858",
                                       "relation_type": "inhibits", "brenda_reference_count": "3",
                                       "enzyme_routes": "accession stated by BRENDA"}])
    added = new_edges[new_edges.relation_type == "inhibits"]
    assert len(added) == 1
    assert added.source_id.iloc[0] == "MAM00001c"
    assert summary["edges_emitted"] == 1


def test_signs_follow_the_relation_and_cofactor_of_is_positive():
    nodes, edges = graph_with_two_compartments()
    rows = [{"base_metabolite_id": "MAM00001", "uniprot_accession": "P29218", "relation_type": relation,
             "brenda_reference_count": "1", "enzyme_routes": "accession stated by BRENDA"}
            for relation in ("inhibits", "activates", "cofactor_of")]
    _, new_edges, new_relations, summary = add_brenda_enzyme_regulation(nodes, edges, ["activates"], rows)
    emitted = new_edges[new_edges.evidence_source.astype(str).str.startswith("BRENDA")]
    assert dict(zip(emitted.relation_type, emitted.sign)) == {"inhibits": -1.0, "activates": 1.0, "cofactor_of": 1.0}
    assert SIGN_OF_RELATION["cofactor_of"] == 1.0  # an enzyme without its cofactor does not turn over
    assert "cofactor_of" in new_relations and "inhibits" in new_relations


def test_a_pair_the_graph_already_carries_is_not_doubled():
    nodes, edges = graph_with_two_compartments()
    _, new_edges, _, summary = add_brenda_enzyme_regulation(
        nodes, edges, ["activates"], [{"base_metabolite_id": "MAM00002", "uniprot_accession": "P17858",
                                       "relation_type": "activates", "brenda_reference_count": "2",
                                       "enzyme_routes": "accession stated by BRENDA"}])
    assert summary["skipped"]["already_present"] == 1
    assert summary["edges_emitted"] == 0
    assert len(new_edges) == len(edges)


def test_a_repeated_pair_inside_brenda_is_emitted_once():
    nodes, edges = graph_with_two_compartments()
    row = {"base_metabolite_id": "MAM00001", "uniprot_accession": "P29218", "relation_type": "inhibits",
           "brenda_reference_count": "5", "enzyme_routes": "accession stated by BRENDA"}
    _, _, _, summary = add_brenda_enzyme_regulation(nodes, edges, ["activates"], [row, dict(row)])
    assert summary["edges_emitted"] == 1
    assert summary["skipped"]["duplicate_within_brenda"] == 1


def test_a_row_the_graph_cannot_place_is_counted_rather_than_invented():
    nodes, edges = graph_with_two_compartments()
    _, new_edges, _, summary = add_brenda_enzyme_regulation(
        nodes, edges, ["activates"],
        [{"base_metabolite_id": "MAM99999", "uniprot_accession": "P29218", "relation_type": "inhibits",
          "brenda_reference_count": "1", "enzyme_routes": "accession stated by BRENDA"},
         {"base_metabolite_id": "MAM00001", "uniprot_accession": "Q99999", "relation_type": "inhibits",
          "brenda_reference_count": "1", "enzyme_routes": "accession stated by BRENDA"}])
    assert summary["skipped"]["metabolite_absent"] == 1
    assert summary["skipped"]["enzyme_absent"] == 1
    assert summary["edges_emitted"] == 0


def test_the_reference_count_and_route_travel_as_provenance_not_as_a_weight():
    """Decision 9.5: no source weighting, because no edge source asserts a negative, so specificity is unestimable."""
    nodes, edges = graph_with_two_compartments()
    _, new_edges, _, summary = add_brenda_enzyme_regulation(
        nodes, edges, ["activates"], [{"base_metabolite_id": "MAM00001", "uniprot_accession": "P29218",
                                       "relation_type": "inhibits", "brenda_reference_count": "7",
                                       "enzyme_routes": "EC class naming one human accession"}])
    added = new_edges[new_edges.relation_type == "inhibits"].iloc[0]
    assert "7 references" in added.evidence_source
    assert "EC class naming one human accession" in added.evidence_source
    assert abs(added.sign) == 1.0  # the sign alone, never scaled by the reference count
    assert summary["weighting"].startswith("none")


def test_the_split_side_decides_whether_regulation_arrives_at_a_protein_or_a_gene():
    """The module sits on either side of the gene/protein split, which moves gene relations onto the protein."""
    nodes, _ = graph_with_two_compartments()
    assert enzyme_node_of_accession(nodes)["P29218"] == "PROTEIN:IMPA1"
    pre_split = pd.DataFrame(
        [["GENE:IMPA1", "gene", "IMPA1", None, None, 0, "P29218"]], columns=NODE_COLUMNS)
    assert enzyme_node_of_accession(pre_split)["P29218"] == "GENE:IMPA1"


@pytest.mark.parametrize("relation", sorted(SIGN_OF_RELATION))
def test_every_admitted_relation_carries_a_unit_sign(relation):
    assert abs(SIGN_OF_RELATION[relation]) == 1.0
