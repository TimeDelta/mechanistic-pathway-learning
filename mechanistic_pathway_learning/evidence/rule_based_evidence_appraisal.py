"""Rule-based appraisal of a literature report from document metadata only (evidence class E3, design section 4.3).

The document-grounded language-model appraisal (llm_evidence_appraisal.py) extracts a study-design rubric from
the text of a paper; until a client for it is configured, this module fills the same rubric from what PubMed
and PubTator3 record about the paper, never from outside knowledge: the publication types (randomised trial,
clinical trial, cohort or case-control study, case report, genetic association study, review or meta-analysis)
and the species the paper annotates (human against non-human). The rubric becomes numeric features of the
report and a rubric weight in (0, 1] by a fixed, documented function; the limitations text lists what the
metadata do not say. Sample size, cohort identity and measurement specifics are not available from metadata
and are reported as not_reported.
"""
from __future__ import annotations

from dataclasses import dataclass

HUMAN_TAXONOMY_ID = "9606"
REVIEW_PUBLICATION_TYPES = ("Review", "Systematic Review", "Meta-Analysis")
RANDOMISED_TRIAL_TYPES = ("Randomized Controlled Trial",)
CLINICAL_STUDY_TYPES = ("Clinical Trial", "Controlled Clinical Trial", "Clinical Trial, Phase II", "Clinical Trial, Phase III", "Clinical Trial, Phase IV", "Cohort Studies", "Case-Control Studies", "Observational Study", "Comparative Study", "Multicenter Study", "Pragmatic Clinical Trial")
CASE_REPORT_TYPES = ("Case Reports",)
GENETIC_ASSOCIATION_TYPES = ("Genome-Wide Association Study", "Genetic Association Studies")

STUDY_DESIGN_WEIGHTS: dict[str, float] = {
    "human_randomized_trial": 1.0,
    "human_cohort_or_case_control": 0.9,
    "human_genetic_association": 0.8,
    "human_case_report": 0.6,
    "animal_genetic_perturbation": 0.5,
    "animal_pharmacological": 0.5,
    "review_or_secondary": 0.3,
    "not_reported": 0.5,
}
DIAGNOSIS_LEVEL_DESCRIPTOR_FACTOR = 0.7


@dataclass
class MetadataRubric:
    study_design: str
    species: list[str]
    human_subjects: bool
    animal_only: bool
    measurement_level: str  # "symptom", "diagnosis" or "not_reported"
    publication_year_known: bool
    limitations: str


def study_design_from_metadata(publication_types: list[str], species_taxonomy_ids: list[str], perturbation_type: str) -> str:
    """The design vocabulary of llm_evidence_appraisal.STUDY_DESIGNS from publication types and annotated species."""
    types = set(publication_types or [])
    human = HUMAN_TAXONOMY_ID in species_taxonomy_ids
    non_human = any(taxonomy_id != HUMAN_TAXONOMY_ID for taxonomy_id in species_taxonomy_ids)
    if types & set(REVIEW_PUBLICATION_TYPES):
        return "review_or_secondary"
    if types & set(RANDOMISED_TRIAL_TYPES):
        return "human_randomized_trial"
    if types & set(CLINICAL_STUDY_TYPES):
        return "human_cohort_or_case_control"
    if types & set(GENETIC_ASSOCIATION_TYPES):
        return "human_genetic_association"
    if types & set(CASE_REPORT_TYPES):
        return "human_case_report"
    if non_human and not human:
        return "animal_pharmacological" if perturbation_type == "drug" else "animal_genetic_perturbation"
    return "not_reported"


def appraise_from_metadata(publication_types: list[str], species_taxonomy_ids: list[str], species_names: list[str], perturbation_type: str, descriptor_level: str, publication_year_known: bool) -> MetadataRubric:
    design = study_design_from_metadata(publication_types, species_taxonomy_ids, perturbation_type)
    human = HUMAN_TAXONOMY_ID in species_taxonomy_ids
    animal_only = bool(species_taxonomy_ids) and not human
    measurement_level = "symptom" if descriptor_level == "symptom" else "diagnosis" if descriptor_level == "diagnosis" else "not_reported"
    clauses: list[str] = []
    if not publication_types:
        clauses.append("publication type not reported")
    if not species_taxonomy_ids:
        clauses.append("species not annotated")
    elif animal_only:
        clauses.append("non-human species only (" + ", ".join(sorted(set(species_names))[:3]) + ")")
    if measurement_level == "diagnosis":
        clauses.append("descriptor is a diagnosis, not a symptom")
    if design == "review_or_secondary":
        clauses.append("secondary literature")
    if design == "not_reported":
        clauses.append("study design not inferable from metadata")
    clauses.append("sample size, cohort and measurement not reported in metadata")
    if not publication_year_known:
        clauses.append("publication year unknown")
    return MetadataRubric(design, sorted(set(species_names)), human, animal_only, measurement_level, publication_year_known, "; ".join(clauses))


def rubric_weight_from_metadata(rubric: MetadataRubric) -> float:
    """Fixed weight in (0, 1]: the study-design weight times 0.7 when the descriptor is a diagnosis rather than a symptom."""
    weight = STUDY_DESIGN_WEIGHTS.get(rubric.study_design, STUDY_DESIGN_WEIGHTS["not_reported"])
    if rubric.measurement_level == "diagnosis":
        weight *= DIAGNOSIS_LEVEL_DESCRIPTOR_FACTOR
    return max(0.05, min(1.0, weight))


def rubric_features(rubric: MetadataRubric) -> dict[str, float]:
    return {
        "rubric_human_species": 1.0 if rubric.human_subjects else 0.0,
        "rubric_animal_only": 1.0 if rubric.animal_only else 0.0,
        "rubric_randomized_trial": 1.0 if rubric.study_design == "human_randomized_trial" else 0.0,
        "rubric_case_report": 1.0 if rubric.study_design == "human_case_report" else 0.0,
        "rubric_cohort_or_case_control": 1.0 if rubric.study_design == "human_cohort_or_case_control" else 0.0,
        "rubric_genetic_association": 1.0 if rubric.study_design == "human_genetic_association" else 0.0,
        "rubric_review_or_secondary": 1.0 if rubric.study_design == "review_or_secondary" else 0.0,
        "rubric_symptom_level_descriptor": 1.0 if rubric.measurement_level == "symptom" else 0.0,
        "rubric_publication_year_known": 1.0 if rubric.publication_year_known else 0.0,
    }
