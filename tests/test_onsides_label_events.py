"""OnSIDES label statements: deduplication across generic labels, section and LLT handling, combination labels, qualification and the report column set."""
from __future__ import annotations

import json

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.load_onsides_label_events import ONSIDES_EXTRA_COLUMNS, aggregate_label_statements, meddra_terms_to_target_symptoms, onsides_reports
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import IngredientBridge

CROSSWALK = {"anxiety": "anxiety", "insomnia": "insomnia"}
# a toy lower-level term (not a MedDRA name) that SIDER's label-term pairing would place under the crosswalk preferred term
LOWER_LEVEL_TO_PREFERRED = {"toy anxious feeling": {"anxiety"}}


def toy_tables() -> dict[str, pd.DataFrame]:
    """Labels 1 to 3 are US sertraline labels, label 4 a Japanese one (section NA), label 5 an EU two-ingredient product (sertraline with ingredient 99999)."""
    return {
        "product_adverse_effect": pd.DataFrame({
            "product_label_id": [1, 2, 2, 3, 4, 5, 1], "effect_id": [10, 11, 12, 13, 14, 15, 16], "label_section": ["AR", "AR", "BW", "AR", "NA", "NA", "AR"],
            "effect_meddra_id": [100, 100, 100, 200, 100, 200, 400], "match_method": ["PMB", "PMB", "PMB", "PMB", "string", "PMB", "PMB"], "pred0": [0.0] * 7, "pred1": [4.0, 5.0, 5.5, 4.5, float("nan"), 4.2, 4.1],
        }),
        "product_label": pd.DataFrame({"label_id": [1, 2, 3, 4, 5], "source": ["US", "US", "US", "JP", "EU"], "source_product_name": list("abcde"), "source_product_id": list("12345"), "source_label_url": [""] * 5}),
        "product_to_rxnorm": pd.DataFrame({"label_id": [1, 2, 3, 4, 5], "rxnorm_product_id": ["p1", "p2", "p3", "p4", "p5"]}),
        "vocab_rxnorm_ingredient_to_product": pd.DataFrame({"product_id": ["p1", "p2", "p3", "p4", "p5", "p5"], "ingredient_id": ["36437", "36437", "36437", "36437", "36437", "99999"]}),
        "vocab_rxnorm_ingredient": pd.DataFrame({"rxnorm_id": ["36437", "99999"], "rxnorm_name": ["sertraline", "toy partner"], "rxnorm_term_type": ["Ingredient", "Ingredient"]}),
        "vocab_meddra_adverse_effect": pd.DataFrame({"meddra_id": [100, 200, 300, 400], "meddra_name": ["Anxiety", "Insomnia", "Anxiety", "Toy anxious feeling"], "meddra_term_type": ["PT", "LLT", "LLT", "LLT"]}),
    }


def test_lower_level_terms_match_through_their_preferred_term() -> None:
    vocabulary = toy_tables()["vocab_meddra_adverse_effect"]
    by_name_only = meddra_terms_to_target_symptoms(vocabulary, CROSSWALK)
    assert dict(zip(by_name_only.meddra_id, by_name_only.target_symptom)) == {100: "anxiety", 200: "insomnia", 300: "anxiety"}  # without a pairing an LLT matches on its exact name only
    paired = meddra_terms_to_target_symptoms(vocabulary, CROSSWALK, LOWER_LEVEL_TO_PREFERRED)
    assert dict(zip(paired.meddra_id, paired.target_symptom)) == {100: "anxiety", 400: "anxiety"}  # LLT rows go through the pairing; the PT row by name
    assert paired[paired.meddra_id == 400].meddra_term_type.iloc[0] == "LLT"  # the source term type is kept


def test_statements_deduplicate_labels_and_keep_sections_countries_and_combination_labels_apart() -> None:
    statements = aggregate_label_statements(toy_tables(), CROSSWALK, LOWER_LEVEL_TO_PREFERRED)
    keyed = {(row.target_symptom, int(row.meddra_id), row.label_section, row.source_country): row for row in statements.itertuples()}
    assert keyed[("anxiety", 100, "AR", "US")].label_count == 2 and keyed[("anxiety", 100, "AR", "US")].max_pred1 == 5.0  # two US labels, one statement
    assert keyed[("anxiety", 100, "AR", "US")].single_ingredient_label_count == 2
    assert keyed[("anxiety", 100, "BW", "US")].label_count == 1  # the boxed-warning mention is its own statement
    assert keyed[("anxiety", 400, "AR", "US")].meddra_term_type == "LLT"  # the LLT under the crosswalk PT is a statement of its own, at the PT's symptom
    assert ("insomnia", 200, "AR", "US") not in keyed  # an LLT no longer matches by its own name once the pairing is given
    assert keyed[("anxiety", 100, "NA", "JP")].match_methods == "string"
    assert int(keyed[("anxiety", 100, "AR", "US")].labels_of_ingredient_in_country) == 3
    combination = statements[statements.ingredient_id == "99999"]
    assert len(combination) == 0  # the pairing drops meddra 200; the partner has no other event
    sertraline_eu = statements[(statements.ingredient_id == "36437") & (statements.source_country == "EU")]
    assert len(sertraline_eu) == 0
    without_pairing = aggregate_label_statements(toy_tables(), CROSSWALK)
    eu_rows = {(row.ingredient_id, row.target_symptom): row for row in without_pairing[without_pairing.source_country == "EU"].itertuples()}
    assert eu_rows[("36437", "insomnia")].label_count == 1 and eu_rows[("36437", "insomnia")].single_ingredient_label_count == 0  # only the two-ingredient label carries it
    assert eu_rows[("99999", "insomnia")].single_ingredient_label_count == 0


def sertraline_bridge(**overrides) -> IngredientBridge:
    fields = dict(ingredient_id="36437", ingredient_name="sertraline", identifier_system="RXCUI", bridge_method="unii_unichem", chembl_parent="CHEMBL809",
                  perturbation_id="CID100068617", unified_with_sider=True, sider_drug_name="sertraline", atc_codes_rxnav=["N06AB06"])
    fields.update(overrides)
    return IngredientBridge(**fields)


MECHANISMS = {"CHEMBL809": [{"molecule_chembl_id": "CHEMBL809", "target_chembl_id": "CHEMBL228", "action_type": "INHIBITOR"}]}
TARGETS = {"CHEMBL228": {"organism": "Homo sapiens", "gene_symbols": ["SLC6A4"], "target_type": "SINGLE PROTEIN", "pref_name": "Serotonin transporter"}}
NODES = {"SLC6A4": "GENE:SLC6A4"}


def test_reports_qualify_by_single_target_and_atc_and_carry_the_shared_columns() -> None:
    statements = aggregate_label_statements(toy_tables(), CROSSWALK)
    reports, counts = onsides_reports(statements, {"36437": sertraline_bridge()}, MECHANISMS, TARGETS, NODES)
    assert list(reports.columns[: len(REPORT_COLUMNS)]) == list(REPORT_COLUMNS) and set(ONSIDES_EXTRA_COLUMNS) <= set(reports.columns)
    sertraline = reports[reports.ingredient_id == "36437"]
    combination_only = sertraline[sertraline.single_ingredient_label_count == 0]
    assert len(combination_only) == 1 and not combination_only.qualifies.iloc[0] and combination_only.disqualified_reason.iloc[0] == "combination_label_only"
    assert "combination-product labels only" in combination_only.limitations.iloc[0]
    assert sertraline[sertraline.single_ingredient_label_count > 0].qualifies.all()
    assert counts["qualifying_reports"] == int((sertraline.single_ingredient_label_count > 0).sum()) and counts["disqualified"] == {"combination_label_only": 1}
    assert set(reports.evidence_code) == {"adverse_reactions_section", "boxed_warning_section", "non_us_label_section"}
    assert reports.perturbation_id.unique().tolist() == ["CID100068617"] and reports.source.unique().tolist() == ["OnSIDES-label"]
    assert reports[reports.label_section == "BW"].rubric_boxed_warning.iloc[0] == 1.0
    japanese = reports[reports.source_country == "JP"].iloc[0]
    assert japanese.rubric_nlp_score == 0.0  # string matches carry no model score
    assert japanese.label_section == "NA" and japanese.rubric_adverse_reactions_section == 1.0 and japanese.rubric_boxed_warning == 0.0 and japanese.rubric_us_label == 0.0
    assert "non-US label, section not recorded by OnSIDES" in japanese.limitations and "unsectioned" not in japanese.limitations
    assert all(len(json.loads(nodes)) == 1 for nodes in reports.perturbation_nodes)
    assert (reports.rubric_atc_known == 1.0).all()


def test_reports_apply_the_sider_rule_to_atc_codes_and_name_missing_target_records() -> None:
    statements = aggregate_label_statements(toy_tables(), CROSSWALK)
    sertraline_statements = int((statements.ingredient_id == "36437").sum())  # the partner ingredient has no bridge entry and produces no row
    rejected, rejected_counts = onsides_reports(statements, {"36437": sertraline_bridge(atc_codes_rxnav=["C01AA05"])}, MECHANISMS, TARGETS, NODES)
    assert not rejected.qualifies.any() and rejected_counts["disqualified"] == {"atc_not_nervous_system": sertraline_statements} and rejected_counts["ingredients_without_bridge"] == 1
    waived, waived_counts = onsides_reports(statements, {"36437": sertraline_bridge(atc_codes_rxnav=[])}, MECHANISMS, TARGETS, NODES)
    assert not waived.qualifies.any() and waived_counts["disqualified"] == {"no_atc_code": sertraline_statements}  # no waiver: the SIDER rule needs an ATC N code
    assert (waived.rubric_atc_known == 0.0).all() and "no ATC code known" in waived.limitations.iloc[0]
    single_ingredient_statements = sertraline_statements - 1  # the combination-only statement keeps its own reason ahead of the node test
    unrecorded, unrecorded_counts = onsides_reports(statements, {"36437": sertraline_bridge()}, MECHANISMS, {}, NODES)
    assert not unrecorded.qualifies.any() and unrecorded_counts["disqualified"] == {"target_record_missing": single_ingredient_statements, "combination_label_only": 1}  # a missing target record is not a missing graph gene
    unmapped, unmapped_counts = onsides_reports(statements, {"36437": sertraline_bridge()}, MECHANISMS, TARGETS, {})
    assert unmapped_counts["disqualified"] == {"no_graph_node": single_ingredient_statements, "combination_label_only": 1}


def test_a_two_target_ingredient_qualifies_only_when_the_target_cap_is_lifted() -> None:
    statements = aggregate_label_statements(toy_tables(), CROSSWALK)
    two_mechanisms = {"CHEMBL809": MECHANISMS["CHEMBL809"] + [{"molecule_chembl_id": "CHEMBL809", "target_chembl_id": "CHEMBL238", "action_type": "INHIBITOR"}]}
    two_targets = TARGETS | {"CHEMBL238": {"organism": "Homo sapiens", "gene_symbols": ["SLC6A3"], "target_type": "SINGLE PROTEIN", "pref_name": "Dopamine transporter"}}
    nodes = NODES | {"SLC6A3": "GENE:SLC6A3"}
    capped, capped_counts = onsides_reports(statements, {"36437": sertraline_bridge()}, two_mechanisms, two_targets, nodes)
    assert not capped.qualifies.any() and "no_single_mechanism_target" in capped_counts["disqualified"]
    lifted, _ = onsides_reports(statements, {"36437": sertraline_bridge()}, two_mechanisms, two_targets, nodes, max_drug_targets=None)
    qualifying = lifted[lifted.qualifies]
    assert len(qualifying) and all(sorted(node for node, _, _ in json.loads(seeds)) == ["GENE:SLC6A3", "GENE:SLC6A4"] for seeds in qualifying.perturbation_nodes)
