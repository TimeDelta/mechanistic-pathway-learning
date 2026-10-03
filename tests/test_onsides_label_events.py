"""OnSIDES label statements: deduplication across generic labels, section and LLT handling, qualification and the report column set."""
from __future__ import annotations

import pandas as pd

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.load_onsides_label_events import aggregate_label_statements, onsides_reports
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import IngredientBridge


def toy_tables() -> dict[str, pd.DataFrame]:
    return {
        "product_adverse_effect": pd.DataFrame({
            "product_label_id": [1, 2, 2, 3, 4], "effect_id": [10, 11, 12, 13, 14], "label_section": ["AR", "AR", "BW", "AR", "AR"],
            "effect_meddra_id": [100, 100, 100, 200, 100], "match_method": ["PMB", "PMB", "PMB", "PMB", "string"], "pred0": [0.0] * 5, "pred1": [4.0, 5.0, 5.5, 4.5, float("nan")],
        }),
        "product_label": pd.DataFrame({"label_id": [1, 2, 3, 4], "source": ["US", "US", "US", "JP"], "source_product_name": ["a", "b", "c", "d"], "source_product_id": ["1", "2", "3", "4"], "source_label_url": [""] * 4}),
        "product_to_rxnorm": pd.DataFrame({"label_id": [1, 2, 3, 4], "rxnorm_product_id": ["p1", "p2", "p3", "p4"]}),
        "vocab_rxnorm_ingredient_to_product": pd.DataFrame({"product_id": ["p1", "p2", "p3", "p4"], "ingredient_id": ["36437", "36437", "36437", "36437"]}),
        "vocab_rxnorm_ingredient": pd.DataFrame({"rxnorm_id": ["36437"], "rxnorm_name": ["sertraline"], "rxnorm_term_type": ["Ingredient"]}),
        "vocab_meddra_adverse_effect": pd.DataFrame({"meddra_id": [100, 200, 300], "meddra_name": ["Anxiety", "Insomnia", "Anxiety"], "meddra_term_type": ["PT", "LLT", "LLT"]}),
    }


def test_statements_deduplicate_labels_and_keep_sections_and_countries_apart() -> None:
    statements = aggregate_label_statements(toy_tables(), {"anxiety": "anxiety", "insomnia": "insomnia"})
    keyed = {(row.target_symptom, int(row.meddra_id), row.label_section, row.source_country): row for row in statements.itertuples()}
    assert keyed[("anxiety", 100, "AR", "US")].label_count == 2 and keyed[("anxiety", 100, "AR", "US")].max_pred1 == 5.0  # two US labels, one statement
    assert keyed[("anxiety", 100, "BW", "US")].label_count == 1  # the boxed-warning mention is its own statement
    assert keyed[("insomnia", 200, "AR", "US")].meddra_term_type == "LLT"  # an LLT whose name equals a crosswalk term is kept
    assert keyed[("anxiety", 100, "AR", "JP")].match_methods == "string"
    assert int(keyed[("anxiety", 100, "AR", "US")].labels_of_ingredient_in_country) == 3


def test_reports_qualify_by_single_target_and_atc_and_carry_the_shared_columns() -> None:
    statements = aggregate_label_statements(toy_tables(), {"anxiety": "anxiety", "insomnia": "insomnia"})
    bridge = IngredientBridge(ingredient_id="36437", ingredient_name="sertraline", identifier_system="RXCUI", bridge_method="unii_unichem", chembl_parent="CHEMBL809",
                              perturbation_id="CID100068617", unified_with_sider=True, sider_drug_name="sertraline", atc_codes_rxnav=["N06AB06"])
    mechanisms = {"CHEMBL809": [{"molecule_chembl_id": "CHEMBL809", "target_chembl_id": "CHEMBL228", "action_type": "INHIBITOR"}]}
    targets = {"CHEMBL228": {"organism": "Homo sapiens", "gene_symbols": ["SLC6A4"], "target_type": "SINGLE PROTEIN", "pref_name": "Serotonin transporter"}}
    reports, counts = onsides_reports(statements, {"36437": bridge}, mechanisms, targets, {"SLC6A4": "GENE:SLC6A4"})
    assert list(reports.columns[: len(REPORT_COLUMNS)]) == list(REPORT_COLUMNS)
    assert counts["qualifying_reports"] == len(statements) and reports.qualifies.all()
    assert set(reports.evidence_code) == {"adverse_reactions_section", "boxed_warning_section"}
    assert reports.perturbation_id.unique().tolist() == ["CID100068617"] and reports.source.unique().tolist() == ["OnSIDES-label"]
    assert reports[reports.label_section == "BW"].rubric_boxed_warning.iloc[0] == 1.0
    assert reports[reports.source_country == "JP"].rubric_nlp_score.iloc[0] == 0.0  # string matches carry no model score
    assert all(len(__import__("json").loads(nodes)) == 1 for nodes in reports.perturbation_nodes)
    non_cns = IngredientBridge(ingredient_id="36437", ingredient_name="sertraline", identifier_system="RXCUI", bridge_method="unii_unichem", chembl_parent="CHEMBL809",
                               perturbation_id="CID100068617", unified_with_sider=True, sider_drug_name="sertraline", atc_codes_rxnav=["C01AA05"])
    rejected, rejected_counts = onsides_reports(statements, {"36437": non_cns}, mechanisms, targets, {"SLC6A4": "GENE:SLC6A4"})
    assert not rejected.qualifies.any() and rejected_counts["disqualified"] == {"atc_not_nervous_system": len(statements)}
