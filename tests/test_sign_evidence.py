"""graph/sign_evidence.py: the PubMed ids behind an edge's sign rather than its pair, and CollecTRI's sign decisions
(decision 9.70; docs/knowledge_graph_practices_review.md, M1 and M2)."""
import json

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.edge_evidence import curated_edge_labels
from mechanistic_pathway_learning.graph.sign_evidence import (BASIS_NOT_FOUND, BASIS_UNSIGNED, read_evidences, references_of, sign_decisions_of,
                                                              sign_evidence_columns, sign_is_evidenced, sign_reference_table)


def entry(resource, references, dataset="omnipath", via=None, decision=None):
    return {"resource": resource, "references": references, "dataset": dataset, "via": via, "attrs": {"sign_decision": decision} if decision else {}}


SIGNALLING_ROWS = [
    {"source": "P1", "target": "P2", "source_genesymbol": "A", "target_genesymbol": "B",
     "evidences": {"positive": [entry("SIGNOR", ["11", "12"]), entry("ProtMapper", ["13"])], "negative": [entry("SPIKE", ["21"])],
                   "directed": [entry("SIGNOR", ["31"])], "undirected": []}},
    {"source": "COMPLEX:P3_P4", "target": "P5", "source_genesymbol": "C_D", "target_genesymbol": "E",
     "evidences": {"positive": [entry("Wang", [])], "negative": [], "directed": [], "undirected": []}},
    {"source": "P6", "target": "P7", "source_genesymbol": "F", "target_genesymbol": "G",
     "evidences": {"positive": [entry("KEGG", ["41"], dataset="kinaseextra")], "negative": [], "directed": [], "undirected": []}},
]
TRANSCRIPTION_ROWS = [
    {"source": "P8", "target": "P9", "source_genesymbol": "TF1", "target_genesymbol": "T1",
     "evidences": {"positive": [entry("CollecTRI", ["51", "52"], dataset="collectri", decision="PMID"), entry("ExTRI", ["53"], dataset="collectri", via="CollecTRI")],
                   "negative": [], "directed": [], "undirected": []}},
    {"source": "P8", "target": "P10", "source_genesymbol": "TF1", "target_genesymbol": "T2",
     "evidences": {"positive": [entry("CollecTRI", ["61"], dataset="collectri", decision="default activation")], "negative": [], "directed": [], "undirected": []}},
    {"source": "P8", "target": "P11", "source_genesymbol": "TF1", "target_genesymbol": "T3",
     "evidences": {"positive": [entry("CollecTRI", ["71", "72"], dataset="collectri", decision=["PMID", "default activation"])], "negative": [],
                   "directed": [], "undirected": []}},
]


def test_references_leave_out_text_mining_and_other_datasets():
    entries = SIGNALLING_ROWS[0]["evidences"]["positive"]
    assert references_of(entries) == {"11", "12"}, "ProtMapper is text mining"
    assert references_of(SIGNALLING_ROWS[2]["evidences"]["positive"], datasets=frozenset({"omnipath"})) == set(), "another dataset's entry"
    transcription = TRANSCRIPTION_ROWS[0]["evidences"]["positive"]
    assert references_of(transcription, resource="CollecTRI") == {"51", "52"}, "CollecTRI's own entry, not ExTRI's through it"


def test_a_sign_decision_can_be_a_list():
    assert sign_decisions_of(TRANSCRIPTION_ROWS[2]["evidences"]["positive"]) == {"PMID", "default activation"}
    assert sign_decisions_of(TRANSCRIPTION_ROWS[1]["evidences"]["positive"]) == {"default activation"}


def test_the_reference_table_keeps_the_signs_apart_and_expands_complexes():
    table = sign_reference_table(SIGNALLING_ROWS, "signalling").set_index(["source_symbol", "target_symbol"])
    assert table.loc[("A", "B"), "positive_pubmed_ids"] == "11;12" and table.loc[("A", "B"), "negative_pubmed_ids"] == "21"
    assert table.loc[("P1", "P2"), "positive_pubmed_ids"] == "11;12", "keyed by accession too"
    assert table.loc[("C", "E"), "positive_pubmed_ids"] == "" and table.loc[("D", "E"), "positive_pubmed_ids"] == "", "a complex's members, no reference"
    regulation = sign_reference_table(TRANSCRIPTION_ROWS, "transcription").set_index(["source_symbol", "target_symbol"])
    assert regulation.loc[("TF1", "T2"), "sign_decisions"] == "default activation"
    assert regulation.loc[("TF1", "T3"), "sign_decisions"] == "PMID+default activation"
    with pytest.raises(ValueError):
        sign_reference_table(SIGNALLING_ROWS, "binding")


def graph_with_evidence():
    edges = pd.DataFrame({"source_id": ["PROTEIN:A", "PROTEIN:A", "PROTEIN:A", "PROTEIN:C", "PROTEIN:X", "PROTEIN:TF1", "PROTEIN:TF1", "PROTEIN:TF1", "MAM1c"],
                          "target_id": ["PROTEIN:B", "PROTEIN:B", "PROTEIN:B", "PROTEIN:E", "PROTEIN:Y", "GENE:T1", "GENE:T2", "GENE:T3", "R1"],
                          "relation_type": ["activates", "inhibits", "binds", "activates", "activates", "regulates_transcription_of", "regulates_transcription_of",
                                            "regulates_transcription_of", "substrate_of"],
                          "sign": [1.0, -1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]})
    families = pd.Series(["signalling"] * 5 + ["transcription"] * 3 + ["metabolic"])
    tables = pd.concat([sign_reference_table(SIGNALLING_ROWS, "signalling"), sign_reference_table(TRANSCRIPTION_ROWS, "transcription")], ignore_index=True)
    return edges, families, sign_evidence_columns(edges, families, tables)


def test_each_edge_gets_the_ids_behind_its_own_sign():
    _, _, columns = graph_with_evidence()
    assert columns.loc[0, "sign_records"] == 2 and columns.loc[0, "sign_pubmed_ids"] == "11;12" and columns.loc[0, "sign_basis"] == "references"
    assert columns.loc[1, "sign_records"] == 1 and columns.loc[1, "sign_pubmed_ids"] == "21", "an inhibition reads the negative evidence"
    assert np.isnan(columns.loc[2, "sign_records"]) and columns.loc[2, "sign_basis"] == BASIS_UNSIGNED
    assert columns.loc[3, "sign_records"] == 0 and columns.loc[3, "sign_basis"] == "no reference", "a pathway resource's sign without a paper"
    assert np.isnan(columns.loc[4, "sign_records"]) and columns.loc[4, "sign_basis"] == BASIS_NOT_FOUND, "absent from the downloads: unknown, not zero"
    assert list(columns.loc[5:7, "sign_basis"]) == ["PMID", "default activation", "PMID+default activation"]
    assert columns.loc[8, "sign_basis"] == "" and np.isnan(columns.loc[8, "sign_records"])


def test_the_sign_rule_concerns_signed_edges_and_needs_a_paper_for_a_regulation():
    edges, families, columns = graph_with_evidence()
    evidence = pd.concat([pd.DataFrame({"family": families}), columns], axis=1)
    assert sign_is_evidenced(evidence, 2).tolist() == [True, False, True, False, False, True, False, True]  + [True]
    assert sign_is_evidenced(evidence, 1).tolist() == [True, True, True, False, False, True, False, True, True]
    with pytest.raises(ValueError, match="sign columns"):
        sign_is_evidenced(evidence.drop(columns="sign_basis"), 1)


def test_the_label_rule_takes_the_sign_minimum_only_when_asked():
    edges, families, columns = graph_with_evidence()
    evidence = pd.concat([pd.DataFrame({"family": families, "kind": ["curated"] * 9, "records": [3] * 9, "contested": [False] * 9,
                                        "human_records": [2] * 9}), columns], axis=1)
    evidence.loc[8, "family"] = "curated"
    plain = curated_edge_labels(evidence)
    assert plain.all(), "without the minimum the sign columns are not read"
    signed = curated_edge_labels(evidence, minimum_sign_records=2)
    assert signed.tolist() == [True, False, True, False, False, True, False, True, True]


def test_read_evidences_needs_a_list(tmp_path):
    good, bad = tmp_path / "good.json", tmp_path / "bad.json"
    good.write_text(json.dumps(SIGNALLING_ROWS))
    bad.write_text(json.dumps({"error": "unknown field"}))
    assert len(read_evidences(good)) == 3
    with pytest.raises(ValueError, match="not an OmniPath interactions list"):
        read_evidences(bad)
