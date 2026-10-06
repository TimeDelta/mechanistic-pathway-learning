"""External report tables (OnSIDES, literature) join the report table and the reliability fit without changing the grades of the built-in classes."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assemble_evidence_table import append_external_reports, evidence_record_from_positive_reports, predication_type_for_literature_sources
from mechanistic_pathway_learning.evidence.assign_evidence_grades import assign_evidence_grade, loss_weight_for_record
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS, RubricWeightDefaults


def report(**overrides) -> dict:
    row = {column: "" for column in REPORT_COLUMNS}
    row.update({"report_value": 1, "frequency": None, "frequency_denominator": None, "placebo_flag": False, "onset": None, "sex": None, "pubmed_reference_count": 0,
                "evidence_date": None, "rubric_log_sample_size": 0.0, "rubric_frequency_known": 0.0, "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0,
                "rubric_evidence_code_iea": 0.0, "rubric_placebo_controlled": 0.0, "rubric_curated_synopsis": 0.0, "rubric_causal_association": 0.0, "perturbation_nodes": "[]"})
    row.update(overrides)
    return row


def test_append_keeps_qualifying_rows_fills_defaults_and_preserves_weights(tmp_path: Path) -> None:
    internal = pd.DataFrame([report(report_id="a", perturbation_id="CID1", symptom="anxiety", relation="induces", evidence_class="pharmacological", source="SIDER-label", evidence_code="label_text")])
    internal["rubric_weight"] = 0.8
    external = pd.DataFrame([
        {"report_id": "o1", "perturbation_id": "CID1", "symptom": "anxiety", "relation": "induces", "evidence_class": "pharmacological", "source": "OnSIDES-label", "report_value": 1, "evidence_code": "adverse_reactions_section", "model_description": "human; drug label (OnSIDES v3.1.1, US); ChEMBL targets CHEMBL228:INHIBITOR", "qualifies": True, "rubric_nlp_score": 0.7, "label_count": 4},
        {"report_id": "o2", "perturbation_id": "CID9", "symptom": "anxiety", "relation": "induces", "evidence_class": "pharmacological", "source": "OnSIDES-label", "report_value": 1, "evidence_code": "adverse_reactions_section", "model_description": "x", "qualifies": False, "rubric_nlp_score": 0.1, "label_count": 1},
        {"report_id": "l1", "perturbation_id": "GENEA", "symptom": "anxiety", "relation": "induces", "evidence_class": "literature", "source": "PubTator3-cause", "report_value": 1, "evidence_code": "human_case_report", "model_description": "human; case report", "qualifies": True, "rubric_weight": 0.6, "rubric_human_species": 1.0},
    ])
    path = tmp_path / "external.parquet"
    external.to_parquet(path, index=False)
    combined = append_external_reports(internal, path, RubricWeightDefaults())
    assert combined.report_id.tolist() == ["a", "o1", "l1"]  # the disqualified OnSIDES row is dropped
    assert "label_count" not in combined.columns and "rubric_nlp_score" in combined.columns  # only shared and rubric columns are kept
    assert combined.loc[combined.report_id == "l1", "rubric_weight"].iloc[0] == 0.6  # the literature appraisal's own weight survives
    assert combined.loc[combined.report_id == "o1", "frequency"].isna().all() and combined.loc[combined.report_id == "o1", "limitations"].iloc[0] == ""


def test_records_grade_by_primary_class_and_literature_only_pairs_are_soft_priors() -> None:
    drug_pair = pd.DataFrame([
        report(report_id="s", perturbation_id="CID1", symptom="anxiety", relation="induces", evidence_class="pharmacological", source="SIDER-label", evidence_code="label_text", source_term_id="C1", model_description="human; drug label; ChEMBL targets CHEMBL228:INHIBITOR"),
        report(report_id="o", perturbation_id="CID1", symptom="anxiety", relation="induces", evidence_class="pharmacological", source="OnSIDES-label", evidence_code="adverse_reactions_section", source_term_id="10001718", model_description="human; drug label (OnSIDES v3.1.1, US); ChEMBL targets CHEMBL228:INHIBITOR"),
    ])
    record = evidence_record_from_positive_reports(drug_pair)
    assert assign_evidence_grade(record) == "B" and record.source == "OnSIDES-label;SIDER-label"  # mixed sources are named; SIDER-only pairs keep "SIDER 4.1"
    gene_pair = pd.DataFrame([
        report(report_id="h", perturbation_id="GENEA", symptom="anxiety", relation="induces", evidence_class="monogenic", source="HPO-Orphanet", evidence_code="TAS", source_record_id="ORPHA:1", association_type="Major susceptibility factor in", rubric_causal_association=0.0),
        report(report_id="l", perturbation_id="GENEA", symptom="anxiety", relation="induces", evidence_class="literature", source="PubTator3-cause", evidence_code="human_case_report", source_record_id="PMID:1", references="PMID:1"),
    ])
    record = evidence_record_from_positive_reports(gene_pair)
    assert record.evidence_class == "monogenic" and record.independent_case_series_count == 1 and assign_evidence_grade(record) == "C"  # the literature row neither counts as a case series nor lifts the grade
    literature_only = pd.DataFrame([report(report_id="l", perturbation_id="GENEB", symptom="anxiety", relation="induces", evidence_class="literature", source="PubTator3-cause", evidence_code="human_case_report", source_record_id="PMID:2")])
    record = evidence_record_from_positive_reports(literature_only)
    assert record.evidence_class == "literature" and assign_evidence_grade(record) == "E" and record.predication_type == "CAUSES" and loss_weight_for_record(record) == 0.10
    assert predication_type_for_literature_sources({"PubTator3-associate"}, "induces") == "ASSOCIATED_WITH" and predication_type_for_literature_sources({"CTD-curated"}, "relieves") == "AFFECTS"
    assert predication_type_for_literature_sources({"PubTator3-inhibit"}, "induces") == "ASSOCIATED_WITH" and predication_type_for_literature_sources({"PubTator3"}, "induces") == "COMENTION"
