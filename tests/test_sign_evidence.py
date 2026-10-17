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


# decision 9.72: the two CollecTRI releases in one download, and the regulations whose sign they dispute
from mechanistic_pathway_learning.graph.sign_evidence import (DISPUTE_COLUMN, LATER_TRANSCRIPTION_RELEASE, releases_dispute_the_sign,  # noqa: E402
                                                              sign_is_undisputed, signs_stated_by)

RELEASE_ROWS = [
    {"source": "Q1", "target": "Q2", "source_genesymbol": "TP53", "target_genesymbol": "RB1",
     "evidences": {"positive": [entry(LATER_TRANSCRIPTION_RELEASE, ["1608930", "9233768"], dataset="collectri2"),
                                entry("TRRUST", [], dataset="collectri2", via=LATER_TRANSCRIPTION_RELEASE)],
                   "negative": [entry("CollecTRI", ["1608930"], dataset="collectri", decision="PMID")], "directed": [], "undirected": []}},
    {"source": "Q1", "target": "Q3", "source_genesymbol": "TP53", "target_genesymbol": "CDKN1A",
     "evidences": {"positive": [entry("CollecTRI", ["81", "82"], dataset="collectri", decision="PMID"), entry(LATER_TRANSCRIPTION_RELEASE, ["81"], dataset="collectri2")],
                   "negative": [], "directed": [], "undirected": []}},
    {"source": "Q4", "target": "Q5", "source_genesymbol": "GATA4", "target_genesymbol": "GATA6",
     "evidences": {"positive": [entry(LATER_TRANSCRIPTION_RELEASE, ["91"], dataset="collectri2")],
                   "negative": [entry(LATER_TRANSCRIPTION_RELEASE, ["91"], dataset="collectri2")], "directed": [], "undirected": []}},
    {"source": "Q6", "target": "Q7", "source_genesymbol": "SP1", "target_genesymbol": "PRKCA",
     "evidences": {"positive": [entry("CollecTRI", ["95"], dataset="collectri", decision="default activation")], "negative": [],
                   "directed": [], "undirected": []}},
]


def test_a_release_states_a_sign_only_through_its_own_entries():
    evidences = RELEASE_ROWS[0]["evidences"]
    assert signs_stated_by(evidences, "CollecTRI") == {-1}
    assert signs_stated_by(evidences, LATER_TRANSCRIPTION_RELEASE) == {1}
    assert signs_stated_by(evidences, "TRRUST") == set(), "TRRUST's entry came through CollecTRI2 (via), so TRRUST stated nothing itself"
    assert signs_stated_by(RELEASE_ROWS[2]["evidences"], LATER_TRANSCRIPTION_RELEASE) == {1, -1}


def test_the_releases_dispute_a_sign_when_they_oppose_or_one_states_both():
    assert releases_dispute_the_sign({-1}, {1}) is True
    assert releases_dispute_the_sign({1}, {1}) is False
    assert releases_dispute_the_sign({1}, set()) is False, "a silent release disputes nothing"
    assert releases_dispute_the_sign(set(), {1, -1}) is True, "a release stating both signs is a dispute by itself"
    assert releases_dispute_the_sign(set(), set()) is False


def disputed_graph():
    edges = pd.DataFrame({"source_id": ["GENE:TP53", "GENE:TP53", "GENE:GATA4", "GENE:SP1", "PROTEIN:A", "MAM1c"],
                          "target_id": ["GENE:RB1", "GENE:CDKN1A", "GENE:GATA6", "GENE:PRKCA", "PROTEIN:B", "R1"],
                          "relation_type": ["regulates_transcription_of"] * 4 + ["activates", "substrate_of"],
                          "sign": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]})
    families = pd.Series(["transcription"] * 4 + ["signalling", "metabolic"])
    tables = pd.concat([sign_reference_table(RELEASE_ROWS, "transcription"), sign_reference_table(SIGNALLING_ROWS, "signalling")], ignore_index=True)
    return edges, families, tables


def test_the_reference_table_keeps_each_release_s_signs():
    _, _, tables = disputed_graph()
    regulation = tables[tables.family == "transcription"].set_index(["source_symbol", "target_symbol"])
    assert regulation.loc[("TP53", "RB1"), "release_2023_signs"] == "-1" and regulation.loc[("TP53", "RB1"), "later_release_signs"] == "1"
    assert regulation.loc[("Q1", "Q2"), "later_release_signs"] == "1", "keyed by accession too"
    assert regulation.loc[("GATA4", "GATA6"), "release_2023_signs"] == "" and regulation.loc[("GATA4", "GATA6"), "later_release_signs"] == "-1;1"
    signalling = tables[tables.family == "signalling"]
    assert (signalling.release_2023_signs == "").all() and (signalling.later_release_signs == "").all()


def test_each_regulation_is_marked_when_its_releases_dispute_the_sign():
    edges, families, tables = disputed_graph()
    columns = sign_evidence_columns(edges, families, tables)
    assert columns[DISPUTE_COLUMN].tolist() == [True, False, True, False, False, False]
    assert columns.loc[0, "sign_records"] == 0, "the stored +1 has none of the 2023 release's papers behind it, only CollecTRI2's"
    assert columns.loc[1, "sign_records"] == 2 and columns.loc[3, "sign_basis"] == "default activation"


def test_the_label_rule_drops_a_disputed_regulation_only_when_asked():
    edges, families, tables = disputed_graph()
    columns = sign_evidence_columns(edges, families, tables)
    evidence = pd.concat([pd.DataFrame({"family": families, "kind": ["curated"] * 6, "records": [3] * 6, "contested": [False] * 6,
                                        "human_records": [2] * 6}), columns], axis=1)
    evidence.loc[5, "family"] = "curated"
    assert curated_edge_labels(evidence).all()
    assert curated_edge_labels(evidence, exclude_disputed_signs=True).tolist() == [False, True, False, True, True, True]
    both = curated_edge_labels(evidence, minimum_sign_records=1, exclude_disputed_signs=True)
    assert both.tolist() == [False, True, False, False, True, True], "the default-activation regulation fails the sign rule, not the dispute rule"
    assert sign_is_undisputed(evidence).tolist() == [False, True, False, True, True, True]
    with pytest.raises(ValueError, match=DISPUTE_COLUMN):
        sign_is_undisputed(evidence.drop(columns=DISPUTE_COLUMN))
