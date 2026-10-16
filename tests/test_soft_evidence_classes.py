"""The soft evidence classes (human_association, grade C; animal, grade D) through the assembler: they enter the fit and
the priors, carry negatives against label pairs, and never make or unmake a label (design section 4.2; the owner,
15 October 2026, on incorporating negatives)."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence import clinical_trial_adverse_events as trials
from mechanistic_pathway_learning.evidence import faers_disproportionality as faers
from mechanistic_pathway_learning.evidence.assemble_evidence_table import (
    EVIDENCE_CLASS_RANK,
    SOFT_EVIDENCE_CLASSES,
    aggregate_reports_to_observations,
    evidence_record_from_positive_reports,
    grade_weighting_rows,
    primary_report,
)
from mechanistic_pathway_learning.evidence.assign_evidence_grades import assign_evidence_grade
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS, RubricWeightDefaults, rubric_weights_for_table
from mechanistic_pathway_learning.evaluation.experiment_data import DEFAULT_LABEL_GRADES, SOFT_PRIOR_ONLY_GRADES

TARGETS = {"CID1": "CHEMBL1:ANTAGONIST"}


def report(**overrides) -> dict:
    row = {column: "" for column in REPORT_COLUMNS}
    row.update({"perturbation_type": "drug", "perturbation_label": "examplone", "relation": "induces", "report_value": 1, "frequency": None,
                "frequency_denominator": None, "placebo_flag": False, "onset": None, "sex": None, "pubmed_reference_count": 0, "evidence_date": None,
                "rubric_log_sample_size": 0.0, "rubric_frequency_known": 0.0, "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0,
                "rubric_evidence_code_iea": 0.0, "rubric_placebo_controlled": 0.0, "rubric_curated_synopsis": 0.0, "rubric_causal_association": 0.0,
                "model_description": "human; drug label; ChEMBL targets CHEMBL1:ANTAGONIST", "perturbation_nodes": json.dumps([["GENE:HTR2A", -1.0, 1.0]])})
    row.update(overrides)
    return row


def observations(rows: list[dict]) -> pd.DataFrame:
    table = pd.DataFrame(rows, columns=list(REPORT_COLUMNS))
    table["rubric_weight"] = rubric_weights_for_table(table, RubricWeightDefaults())
    return aggregate_reports_to_observations(table, "curated_synopsis_provenance", {}, set(), TARGETS)


def label(symptom: str = "insomnia", **overrides) -> dict:
    return report(report_id=f"s|{symptom}", perturbation_id="CID1", symptom=symptom, evidence_class="pharmacological", source="SIDER-label", evidence_code="label_text", **overrides)


def trial(symptom: str = "insomnia", value: int = 0, **overrides) -> dict:
    return report(report_id=f"t|{symptom}|{value}", perturbation_id="CID1", symptom=symptom, evidence_class=trials.EVIDENCE_CLASS, source=trials.SOURCE_NAME,
                  evidence_code=trials.EVIDENCE_CODE, report_value=value, frequency=0.0 if value == 0 else 0.1, frequency_denominator=250, placebo_flag=True,
                  rubric_placebo_controlled=1.0, rubric_frequency_known=1.0, **overrides)


def test_the_two_sources_use_the_design_grade_c_class():
    assert trials.EVIDENCE_CLASS == "human_association" and faers.EVIDENCE_CLASS == "human_association"
    assert set(SOFT_EVIDENCE_CLASSES) == {"human_association", "animal"}
    assert "human_association" not in DEFAULT_LABEL_GRADES and "C" not in DEFAULT_LABEL_GRADES and "D" in SOFT_PRIOR_ONLY_GRADES


def test_a_label_pair_keeps_its_grade_and_gains_a_negative_report():
    table = observations([label(), trial(value=0)])
    row = table.iloc[0]
    assert row.grade == "B" and row.evidence_class == "pharmacological"
    assert row.positive_report_count == 1 and row.negative_report_count == 1 and row.report_count == 2
    assert row.weight > 0


def test_a_trial_only_pair_is_grade_c_and_never_a_label():
    table = observations([trial("anxiety", value=1)])
    row = table.iloc[0]
    assert row.evidence_class == "human_association" and row.grade == "C" and row.grade not in DEFAULT_LABEL_GRADES
    # the design holds grade C at weight 0 ("validation only in version 1", DEFAULT_GRADE_WEIGHTS); the row is in the table
    # and the reliability fit, and a per-source tier strength is what would give it a loss weight
    assert row.positive_report_count == 1 and row.weight == 0.0


def test_a_trial_only_absence_has_no_grade_weighting_row_and_weight_zero():
    """A trial arm with no event has frequency 0 and a denominator, which for a phenotype.hpoa row would be the 0 of N
    absence the grade weighting reads; the rule is the monogenic loader's and does not reach a trial row."""
    rows = pd.DataFrame([trial("anxiety", value=0)], columns=list(REPORT_COLUMNS))
    assert grade_weighting_rows(rows).empty
    row = observations([trial("anxiety", value=0)]).iloc[0]
    assert row.positive_report_count == 0 and row.negative_report_count == 1 and row.weight == 0.0


def test_the_primary_class_is_the_best_present_whatever_the_order():
    """A soft report that arrived first (an --extra-reports file passed before OnSIDES) must not demote a label pair."""
    soft_first = pd.DataFrame([trial(value=0), label()], columns=list(REPORT_COLUMNS))
    assert primary_report(soft_first).evidence_class == "pharmacological"
    assert evidence_record_from_positive_reports(soft_first).evidence_class == "pharmacological"
    assert observations([trial(value=0), label()]).iloc[0].grade == "B"
    assert EVIDENCE_CLASS_RANK["monogenic"] < EVIDENCE_CLASS_RANK["pharmacological"] < EVIDENCE_CLASS_RANK["human_association"] < EVIDENCE_CLASS_RANK["animal"] < EVIDENCE_CLASS_RANK["literature"]


def test_a_knockout_only_gene_pair_is_grade_d_grouped_as_its_own_cluster():
    knockout = report(report_id="k|hyperactivity", perturbation_id="DRD2", perturbation_type="gene", perturbation_label="DRD2", symptom="hyperactivity",
                      evidence_class="animal", source="IMPC-knockout", evidence_code="impc_significant", model_description="mouse knockout Drd2",
                      perturbation_nodes=json.dumps([["GENE:DRD2", -1.0, 1.0]]))
    table = pd.DataFrame([knockout], columns=list(REPORT_COLUMNS))
    table["rubric_weight"] = rubric_weights_for_table(table, RubricWeightDefaults())
    row = aggregate_reports_to_observations(table, "curated_synopsis_provenance", {}, set()).iloc[0]
    assert row.evidence_class == "animal" and row.grade == "D" and row.grade in SOFT_PRIOR_ONLY_GRADES
    assert row.group_id == "DRD2" and row.disease_cluster_id == "cluster:DRD2" and row.omim_entry_count == 0 and row.source == "IMPC-knockout"


def test_a_gene_label_pair_with_a_knockout_report_keeps_its_disease_provenance():
    omim = report(report_id="o|1", perturbation_id="HMBS", perturbation_type="gene", perturbation_label="HMBS", symptom="psychosis", evidence_class="monogenic",
                  source="HPO-OMIM", evidence_code="TAS", source_record_id="OMIM:176000", references="OMIM:176000", rubric_curated_synopsis=1.0,
                  rubric_evidence_code_tas=1.0, model_description="human loss-of-function; OMIM:176000", perturbation_nodes=json.dumps([["GENE:HMBS", -1.0, 1.0]]))
    knockout = report(report_id="k|1", perturbation_id="HMBS", perturbation_type="gene", perturbation_label="HMBS", symptom="psychosis", evidence_class="animal",
                      source="IMPC-knockout", evidence_code="impc_not_significant", report_value=0, source_record_id="MGI:1", perturbation_nodes=json.dumps([["GENE:HMBS", -1.0, 1.0]]))
    table = pd.DataFrame([knockout, omim], columns=list(REPORT_COLUMNS))
    table["rubric_weight"] = rubric_weights_for_table(table, RubricWeightDefaults())
    row = aggregate_reports_to_observations(table, "curated_synopsis_provenance", {"HMBS": "cluster:porphyria"}, set()).iloc[0]
    assert row.evidence_class == "monogenic" and row.grade == "A" and row.disease_identifiers == "OMIM:176000" and row.omim_entry_count == 1
    assert row.disease_cluster_id == "cluster:porphyria" and row.negative_report_count == 1 and "IMPC-knockout" in row.source


@pytest.mark.parametrize("evidence_class, grade", [("human_association", "C"), ("animal", "D")])
def test_the_grade_ladder_places_the_soft_classes_below_the_labels(evidence_class, grade):
    from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord
    assert assign_evidence_grade(EvidenceRecord("x", "y", "induces", evidence_class)) == grade
