"""The physiology audit's checks on small fixtures (experiments/audit_physiology.py)."""
from __future__ import annotations

import json

import pandas as pd

from experiments.audit_physiology import (
    EXPECTED_SIGN,
    compartment_consistency,
    coverage,
    graph_integrity,
    hpo_ancestors,
    perturbation_targets,
    read_validation_set,
    term_read_as,
    validation_set_checks,
)


def graph():
    nodes = pd.DataFrame({
        "node_id": ["MAM1c", "MAM1m", "MAM2c", "R1", "R2", "GENE:A", "GENE:B", "GENE:LONE"],
        "node_type": ["metabolite", "metabolite", "metabolite", "reaction", "reaction", "gene", "gene", "gene"],
        "display_name": ["x", "x", "y", "r1", "r2", "A", "B", "LONE"], "compartment": ["c", "m", "c", None, None, None, None, None],
        "is_transport": [False, False, False, False, True, False, False, False], "gene_symbol": [None, None, None, None, None, "A", "B", "LONE"],
        "is_currency": [True, False, False, False, False, False, False, False]})
    edges = pd.DataFrame({
        "source_id": ["MAM1c", "R1", "MAM1c", "R2", "GENE:A", "GENE:A", "GENE:A", "GENE:B"],
        "target_id": ["R1", "MAM1m", "R2", "MAM2c", "GENE:B", "GENE:B", "GENE:B", "GENE:B"],
        "relation_type": ["substrate_of", "product_of", "substrate_of", "product_of", "activates", "activates", "inhibits", "binds"],
        "sign": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 0.0], "evidence_source": ["Human-GEM"] * 4 + ["OmniPath"] * 4})
    return nodes, edges


def test_graph_integrity_counts_duplicates_contradictions_unsigned_edges_self_loops_and_isolates():
    nodes, edges = graph()
    result = graph_integrity(nodes, edges)
    assert result["exact_duplicate_edge_rows"] == 2 and result["exact_duplicate_edge_pairs"] == 1  # GENE:A activates GENE:B twice
    assert result["contradictory_sign_pairs"] == 0, "activates and inhibits are different relations; a contradiction is both signs under one relation"
    assert result["unsigned_fraction_by_relation"]["binds"] == 1.0 and result["self_loops"] == 1 and result["isolated_nodes"] == 1
    assert result["stoichiometry"] == {"column_present": False}
    assert result["top_degree_nodes"][0]["node"] in ("GENE:B", "GENE:A")


def test_a_pair_carrying_both_signs_under_one_relation_is_a_contradiction():
    nodes, edges = graph()
    edges = pd.concat([edges, pd.DataFrame([{"source_id": "GENE:A", "target_id": "GENE:B", "relation_type": "activates", "sign": -1.0, "evidence_source": "CollecTRI"}])])
    result = graph_integrity(nodes, edges)
    assert result["contradictory_sign_pairs"] == 1 and result["contradictory_examples"][0]["source_id"] == "GENE:A"


def test_stoichiometry_is_reported_when_the_column_is_present():
    nodes, edges = graph()
    edges["stoichiometric_coefficient"] = [2.0, 1.0, 1.0, 1.0, None, None, None, None]
    result = graph_integrity(nodes, edges)["stoichiometry"]
    assert result == {"column_present": True, "metabolic_edges": 4, "edges_with_coefficient_other_than_one": 1}


def test_a_currency_tag_on_a_node_that_is_not_a_metabolite_is_reported():
    """The fourth currency defect: GENE:CA2 was tagged because its symbol spells the name list's calcium."""
    nodes, edges = graph()
    assert graph_integrity(nodes, edges)["currency_tags_on_other_node_types"] == []
    nodes.loc[nodes.node_id == "GENE:A", "is_currency"] = True
    assert graph_integrity(nodes, edges)["currency_tags_on_other_node_types"] == ["GENE:A"]


def test_a_reaction_spanning_compartments_without_the_transport_flag_is_found():
    nodes, edges = graph()
    result = compartment_consistency(nodes, edges)
    assert result["reactions_spanning_compartments_without_transport_flag"] == 1 and result["examples"] == ["R1"], "R1 takes x in c and makes x in m without being transport; R2 is flagged"
    assert result["metabolites_by_compartment"] == {"c": 2, "m": 1}


def test_perturbation_targets_read_symbols_and_signs_from_the_node_json():
    assert perturbation_targets(json.dumps([["GENE:DRD2", -1.0, 1.0], ["PROTEIN:HTR2A", 1.0, 0.5]])) == {"DRD2": -1.0, "HTR2A": 1.0}
    assert perturbation_targets("") == {} and perturbation_targets("not json") == {}


def evidence_table():
    return pd.DataFrame([
        {"perturbation_id": "CID1", "perturbation_type": "drug", "perturbation_label": "Haloperidol", "symptom": "parkinsonism", "relation": "induces", "grade": "B",
         "positive_report_count": 2, "negative_report_count": 0, "evidence_class": "pharmacological", "perturbation_nodes": json.dumps([["GENE:DRD2", -1.0, 1.0]])},
        {"perturbation_id": "CID1", "perturbation_type": "drug", "perturbation_label": "Haloperidol", "symptom": "elevated_mood_or_mania", "relation": "induces", "grade": "C",
         "positive_report_count": 1, "negative_report_count": 2, "evidence_class": "human_association", "perturbation_nodes": json.dumps([["GENE:DRD2", -1.0, 1.0]])},
        {"perturbation_id": "CID2", "perturbation_type": "drug", "perturbation_label": "domperidone", "symptom": "parkinsonism", "relation": "induces", "grade": "B",
         "positive_report_count": 1, "negative_report_count": 0, "evidence_class": "pharmacological", "perturbation_nodes": json.dumps([["GENE:DRD2", -1.0, 1.0]])},
        {"perturbation_id": "HTT", "perturbation_type": "gene", "perturbation_label": "HTT", "symptom": "depressed_mood", "relation": "induces", "grade": "A",
         "positive_report_count": 1, "negative_report_count": 0, "evidence_class": "monogenic", "perturbation_nodes": json.dumps([["GENE:HTT", -1.0, 1.0]])},
    ])


def test_validation_verdicts_cover_supported_contradicted_respected_absent_and_sign_disagreement():
    nodes = pd.DataFrame({"node_id": ["GENE:DRD2", "GENE:HTT"], "node_type": ["gene", "gene"], "gene_symbol": ["DRD2", "HTT"]})
    rows = [
        {"kind": "chain", "perturbation_type": "drug", "perturbation": "haloperidol", "target_gene": "DRD2", "direction": "antagonist", "symptom": "parkinsonism", "relation": "induces", "tier": "A"},
        {"kind": "negative", "perturbation_type": "drug", "perturbation": "haloperidol", "target_gene": "DRD2", "direction": "antagonist", "symptom": "elevated_mood_or_mania", "relation": "absent", "tier": "B"},
        {"kind": "negative", "perturbation_type": "drug", "perturbation": "domperidone", "target_gene": "DRD2", "direction": "antagonist", "symptom": "parkinsonism", "relation": "absent", "tier": "A"},
        {"kind": "chain", "perturbation_type": "gene", "perturbation": "HTT", "target_gene": "HTT", "direction": "gain_of_function", "symptom": "depressed_mood", "relation": "induces", "tier": "A"},
        {"kind": "chain", "perturbation_type": "drug", "perturbation": "loperamide", "target_gene": "OPRM1", "direction": "agonist", "symptom": "somnolence_or_hypersomnia", "relation": "induces", "tier": "A"},
        {"kind": "withdrawn", "perturbation_type": "drug", "perturbation": "x", "target_gene": "Y", "direction": "agonist", "symptom": "psychosis", "relation": "absent", "tier": "B", "status": "withdrawn: the label records it"},
    ]
    result = validation_set_checks(rows, nodes, evidence_table())
    verdicts = {(row["perturbation"], row["symptom"]): row["verdict"] for row in result["results"]}
    assert verdicts[("haloperidol", "parkinsonism")] == "supported"
    assert verdicts[("haloperidol", "elevated_mood_or_mania")] == "respected_with_negative", "a grade-C row with negatives is not a label, so the negative is respected"
    assert verdicts[("domperidone", "parkinsonism")] == "contradicted"
    assert verdicts[("HTT", "depressed_mood")] == "supported" and result["sign_disagreements"][0]["perturbation"] == "HTT"
    assert verdicts[("loperamide", "somnolence_or_hypersomnia")] == "perturbation_absent" and verdicts[("x", "psychosis")] == "withdrawn"
    assert result["summary"]["chain"] == {"supported": 2, "perturbation_absent": 1} and result["summary"]["negative"] == {"respected_with_negative": 1, "contradicted": 1}
    assert result["by_tier"]["negative:A"] == {"contradicted": 1}
    assert EXPECTED_SIGN["antagonist"] == -1.0 and EXPECTED_SIGN["gain_of_function"] == 1.0


def test_a_chain_whose_term_the_crosswalk_excludes_is_a_construct_decision_not_a_missing_label():
    parents = {"HP:0001249": {"HP:0011446"}, "HP:0100543": {"HP:0011446"}, "HP:0000742": {"HP:0100716"}}
    expansion = {"cognitive_impairment": ({"HP:0100543"}, set()), "self_injury": ({"HP:0100716"}, {"HP:5200330"})}
    assert term_read_as("HP:0001249", "cognitive_impairment", expansion, parents) is False
    assert term_read_as("HP:0000742", "self_injury", expansion, parents) is True
    assert term_read_as("HP:0000742", "unknown", expansion, parents) is None
    assert hpo_ancestors("HP:0000742", parents) == {"HP:0100716"}
    nodes = pd.DataFrame({"node_id": ["GENE:PAH"], "node_type": ["gene"], "gene_symbol": ["PAH"]})
    evidence = pd.DataFrame([{"perturbation_id": "PAH", "perturbation_type": "gene", "perturbation_label": "PAH", "symptom": "anxiety", "relation": "induces", "grade": "A",
                              "positive_report_count": 1, "negative_report_count": 0, "evidence_class": "monogenic", "perturbation_nodes": json.dumps([["GENE:PAH", -1.0, 1.0]])}])
    rows = [{"kind": "chain", "perturbation_type": "gene", "perturbation": "PAH", "target_gene": "PAH", "direction": "loss_of_function", "symptom": "cognitive_impairment",
             "relation": "induces", "tier": "A", "hpo_term": "HP:0001249"}]
    result = validation_set_checks(rows, nodes, evidence, expansion, parents)
    assert result["results"][0]["verdict"] == "term_outside_crosswalk"


def test_coverage_counts_labels_soft_rows_and_negatives_per_symptom():
    result = coverage(evidence_table())
    by_symptom = {row["symptom"]: row for row in result["by_symptom"]}
    assert by_symptom["parkinsonism"] == {"symptom": "parkinsonism", "label_pairs": 2, "label_genes": 0, "label_drugs": 2, "soft_pairs": 0, "pairs_with_a_negative_report": 0,
                                          "pairs_all_negative": 0, "classes": ["pharmacological"]}
    assert by_symptom["elevated_mood_or_mania"]["soft_pairs"] == 1 and by_symptom["elevated_mood_or_mania"]["pairs_with_a_negative_report"] == 1
    assert result["symptoms_without_gene_labels"] == ["parkinsonism", "elevated_mood_or_mania"] and result["symptoms_with_any_negative"] == 1


def test_the_committed_validation_set_reads_with_every_row_complete():
    rows = read_validation_set()
    assert len(rows) >= 50
    for row in rows:
        assert row["kind"] in ("chain", "negative", "withdrawn") and row["perturbation_type"] in ("drug", "gene") and row["tier"] in ("A", "B", "C")
        assert row["direction"] in EXPECTED_SIGN and row["citation"] and row["quote"]
        assert row["relation"] in ("induces", "relieves", "absent")
        if row["kind"] == "withdrawn":
            assert row["status"].startswith("withdrawn")


def test_the_therapeutic_direction_counts_the_loss_from_disorder_to_symptom(tmp_path):
    import csv as csv_module
    import gzip as gzip_module

    from experiments.audit_physiology import therapeutic_direction

    crosswalk = tmp_path / "crosswalk.csv"
    with crosswalk.open("w", newline="", encoding="utf-8") as handle:
        writer = csv_module.writer(handle)
        writer.writerow(["target_symptom", "meddra_preferred_terms"])
        writer.writerow(["anxiety", "Anxiety;Nervousness"])
        writer.writerow(["insomnia", "Insomnia"])
    indications = tmp_path / "meddra_all_indications.tsv.gz"
    rows = [  # STITCH_flat, UMLS_label_cui, detection_method, label_concept_name, MedDRA_type, UMLS_meddra_cui, meddra_concept_name
        ["CID100000001", "C1", "NLP_indication", "x", "PT", "C1", "Anxiety"],            # a study drug, a mapped symptom
        ["CID100000001", "C2", "NLP_indication", "x", "PT", "C2", "Schizophrenia"],      # a study drug, a disorder the crosswalk does not map
        ["CID100000002", "C3", "NLP_indication", "x", "PT", "C3", "Insomnia"],           # the second study drug
        ["CID100000002", "C4", "NLP_precondition", "x", "PT", "C4", "Anxiety"],          # not an indication
        ["CID100000002", "C5", "NLP_indication", "x", "LLT", "C5", "Anxiety"],           # not preferred-term level
        ["CID100000009", "C6", "NLP_indication", "x", "PT", "C6", "Anxiety"],            # not a study drug
    ]
    with gzip_module.open(indications, "wt", encoding="utf-8") as handle:
        csv_module.writer(handle, delimiter="\t").writerows(rows)
    evidence = pd.DataFrame([
        {"perturbation_id": "CID100000001", "perturbation_type": "drug", "relation": "relieves", "symptom": "anxiety"},
        {"perturbation_id": "CID100000002", "perturbation_type": "drug", "relation": "induces", "symptom": "insomnia"},
        {"perturbation_id": "CID100000003", "perturbation_type": "drug", "relation": "induces", "symptom": "anxiety"},
        {"perturbation_id": "GENE1", "perturbation_type": "gene", "relation": "induces", "symptom": "anxiety"},
    ])
    result = therapeutic_direction(evidence, indications, crosswalk)
    assert result["indications_read"] and result["crosswalk_preferred_terms"] == 3 and result["study_drugs"] == 3
    assert result["relieves_pairs_in_the_evidence"] == 1
    assert result["indication_rows"] == 4, "preferred-term NLP_indication rows only, over every drug"
    assert result["indication_rows_on_a_study_drug"] == 3 and result["indication_rows_naming_a_study_symptom"] == 2
    assert result["distinct_pairs_naming_a_study_symptom"] == 2
    assert result["study_drugs_with_no_indication_row"] == 1, "CID100000003 has none"
    assert result["most_common_indication_terms_outside_the_crosswalk"] == [("Schizophrenia", 1)]
    assert therapeutic_direction(evidence, tmp_path / "absent.tsv.gz", crosswalk) == {"relieves_pairs_in_the_evidence": 1, "indications_read": False}
