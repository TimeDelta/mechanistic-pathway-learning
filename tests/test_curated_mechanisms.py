"""graph/curated_mechanisms.py and docs/curated_dead_end_mechanisms.csv: curated edges for dead-end perturbations
(decision 9.68 (a))."""
import re
from pathlib import Path

import pandas as pd
import pytest

from experiments.build_neuronal_graph_variant import CURATED_REACTIONS
from mechanistic_pathway_learning.graph.curated_mechanisms import (DEFAULT_CURATED_MECHANISMS_PATH, EVIDENCE_PREFIX, NODE_TYPES_OF_RELATION,
                                                                   REQUIRED_COLUMNS, SIGN_OF_RELATION, add_curated_mechanisms,
                                                                   problems_of_row, pubmed_ids_of, read_curated_mechanisms)
from mechanistic_pathway_learning.graph.edge_evidence import CURATED_PMID, family_of

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = REPOSITORY_ROOT / DEFAULT_CURATED_MECHANISMS_PATH


def small_graph():
    nodes = pd.DataFrame({"node_id": ["GENE:A", "GENE:B", "GENE:C", "MAR1", "MAM1m"],
                          "node_type": ["gene", "gene", "gene", "reaction", "metabolite"],
                          "degree": [1, 1, 0, 1, 0]})
    edges = pd.DataFrame({"source_id": ["GENE:A"], "target_id": ["MAR1"], "relation_type": ["catalyzed_by"], "sign": [1.0],
                          "evidence_source": ["Human-GEM GPR"], "stoichiometric_coefficient": [float("nan")], "weight": [float("nan")]})
    relations = ["substrate_of", "product_of", "catalyzed_by", "activates", "inhibits", "cofactor_of"]
    return nodes, edges, relations


def table_of(*rows):
    defaults = {"sign": "1", "decision": "9.68", "pmids": "123456", "quotes": "\"a quote\" (PMID 123456)", "reading": "why"}
    return pd.DataFrame([{**defaults, **row} for row in rows], columns=list(REQUIRED_COLUMNS))


def test_rows_become_edges_with_each_pubmed_id_in_the_evidence():
    nodes, edges, relations = small_graph()
    table = table_of({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "activates", "pmids": "111111;222222"},
                     {"source_id": "MAM1m", "target_id": "GENE:C", "relation_type": "cofactor_of"},
                     {"source_id": "GENE:B", "target_id": "MAR1", "relation_type": "catalyzed_by"})
    new_nodes, new_edges, new_relations, summary = add_curated_mechanisms(nodes, edges, relations, table)
    added = new_edges.iloc[len(edges):]
    assert summary["edges_added"] == 3 and summary["already_present"] == []
    assert list(zip(added.source_id, added.target_id, added.relation_type)) == [("GENE:B", "GENE:C", "activates"), ("MAM1m", "GENE:C", "cofactor_of"),
                                                                               ("GENE:B", "MAR1", "catalyzed_by")]
    assert set(added.sign) == {1.0}
    first = added.iloc[0].evidence_source
    assert first.startswith(EVIDENCE_PREFIX) and "decision 9.68" in first
    assert CURATED_PMID.findall(first) == ["111111", "222222"], "the edge evidence table must find every id"
    assert family_of(first) == "curated", "the edge evidence table reads the edge as a curated one"
    assert dict(zip(new_nodes.node_id, new_nodes.degree))["GENE:C"] == 2
    assert new_relations == relations
    assert new_edges.dtypes.to_dict() == edges.dtypes.to_dict()


def test_an_edge_the_graph_already_holds_is_left_out_and_reported():
    nodes, edges, relations = small_graph()
    table = table_of({"source_id": "GENE:A", "target_id": "MAR1", "relation_type": "catalyzed_by"})
    _, new_edges, _, summary = add_curated_mechanisms(nodes, edges, relations, table)
    assert len(new_edges) == len(edges)
    assert summary["already_present"] == ["GENE:A -catalyzed_by-> MAR1"]


def test_an_inhibition_carries_a_negative_sign():
    nodes, edges, relations = small_graph()
    table = table_of({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "inhibits", "sign": "-1"})
    _, new_edges, _, _ = add_curated_mechanisms(nodes, edges, relations, table)
    assert new_edges.iloc[-1].sign == -1.0


@pytest.mark.parametrize("row, expected", [
    ({"source_id": "GENE:Z", "target_id": "GENE:C", "relation_type": "activates"}, "source GENE:Z is not a node"),
    ({"source_id": "GENE:B", "target_id": "MAR1", "relation_type": "activates"}, "target MAR1 is a reaction node"),
    ({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "activates", "sign": "-1"}, "does not match activates"),
    ({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "binds"}, "is not one this table may add"),
    ({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "activates", "pmids": "PMC123"}, "not a PubMed id"),
    ({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "activates", "pmids": ""}, "no PubMed id"),
    ({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "activates", "quotes": " "}, "no quote"),
    ({"source_id": "GENE:B", "target_id": "MAM1m", "relation_type": "cofactor_of"}, "source GENE:B is a gene node"),
])
def test_a_row_that_cannot_be_added_as_written_raises(row, expected):
    nodes, edges, relations = small_graph()
    with pytest.raises(ValueError, match=expected):
        add_curated_mechanisms(nodes, edges, relations, table_of(row))


def test_a_relation_the_graph_lacks_is_a_problem():
    nodes, _, _ = small_graph()
    row = table_of({"source_id": "GENE:B", "target_id": "GENE:C", "relation_type": "inhibits", "sign": "-1"}).iloc[0]
    problems = problems_of_row(row, dict(zip(nodes.node_id, nodes.node_type)), ["activates"])
    assert problems == ["relation 'inhibits' is not among the graph's relation types"]


def test_pubmed_ids_are_split_on_semicolons():
    assert pubmed_ids_of("28132689; 42497073") == ["28132689", "42497073"]
    assert pubmed_ids_of("") == []
    with pytest.raises(ValueError):
        pubmed_ids_of("doi:10.1/x")


def test_the_table_has_every_column_and_each_row_is_consistent():
    table = read_curated_mechanisms(TABLE_PATH)
    assert set(REQUIRED_COLUMNS) <= set(table.columns)
    for row in table.itertuples():
        assert row.relation_type in SIGN_OF_RELATION
        assert float(row.sign) == SIGN_OF_RELATION[row.relation_type]
        assert row.decision == "9.68"
        pmids = pubmed_ids_of(row.pmids)
        assert pmids, f"{row.source_id} -> {row.target_id} cites no paper"
        for identifier in pmids:
            assert f"PMID {identifier}" in row.quotes, f"{row.source_id} -> {row.target_id}: PMID {identifier} has no quote"
        assert len(str(row.dois).split(";")) == len(pmids), "one DOI per PubMed id, in the same order"
        quoted = re.findall(r'"([^"]+)"', row.quotes)
        assert quoted and all(len(text) > 20 for text in quoted), "each paper is quoted, not paraphrased"
        source_type, target_type = NODE_TYPES_OF_RELATION[row.relation_type]
        assert row.source_id.startswith("GENE:") == (source_type == "gene")
        assert row.target_id.startswith("GENE:") == (target_type == "gene")


def test_the_table_holds_the_mechanisms_decision_9_68_names():
    table = read_curated_mechanisms(TABLE_PATH)
    edges = set(zip(table.source_id, table.target_id, table.relation_type))
    for hydroxylase in ("PAH", "TH", "TPH1", "TPH2"):
        assert ("GENE:DNAJC12", f"GENE:{hydroxylase}", "activates") in edges
    assert ("GENE:MMACHC", "MAR08613", "catalyzed_by") in edges, "MMACHC processes cob(III)alamin to cob(II)alamin"
    assert ("MAM01600m", "GENE:MMUT", "cofactor_of") in edges, "adenosylcobalamin is methylmalonyl-CoA mutase's cofactor"


def test_the_lysosomal_cobalamin_release_is_a_curated_reaction_catalysed_by_both_exporters():
    reactions = {reaction.reaction_id: reaction for reaction in CURATED_REACTIONS}
    release = reactions["MAR_LYSOSOMAL_COBALAMIN_RELEASE"]
    assert release.substrates == ("MAM01361e",) and release.products == ("MAM01361c",), "extracellular to cytosolic aquacob(III)alamin"
    assert set(release.catalysts) == {"LMBRD1", "ABCD4"}
    assert "PMID 28572511" in release.evidence
    assert release.release is None and release.membrane_potential_sign is None
