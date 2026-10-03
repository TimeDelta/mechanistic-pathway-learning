"""Evidence grading (experiment design, section 4.3).

Every perturbation-symptom observation carries an evidence record. The grade
decides the role (training label, soft prior or validation only) and the loss
weight. Grades are deliberately coarse; the weights are the initial values from
the design and are overridden by configs/evidence_assembly.yaml.

Grade A  human loss of function in a graph gene with the symptom documented in
         at least two independent case series or an OMIM clinical synopsis
Grade B  CNS-penetrant drug with a dominant target and the event on the label
Grade C  human association (gene-level GWAS statistic, metabolomic association,
         or a monogenic report below the grade A bar); validation only in version 1
Grade D  animal perturbation with a symptom analogue; soft prior
Grade E  literature predication or co-mention; soft prior, weight by predication type
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date

EvidenceGrade = str  # "A", "B", "C", "D" or "E"

DEFAULT_PREDICATION_WEIGHTS: dict[str, float] = {
    "CAUSES": 0.10,
    "DISRUPTS": 0.08,
    "PREDISPOSES": 0.08,
    "AFFECTS": 0.08,
    "ASSOCIATED_WITH": 0.05,
    "COMENTION": 0.02,
}

DEFAULT_GRADE_WEIGHTS: dict[EvidenceGrade, float] = {
    "A": 1.0,
    "B": 0.6,
    "C": 0.0,  # validation only in version 1
    "D": 0.2,
    "E": 0.05,
}


@dataclass
class EvidenceRecord:
    perturbation_identifier: str
    symptom_identifier: str
    relation: str  # "induces" or "relieves"
    evidence_class: str  # "monogenic", "pharmacological", "human_association", "animal", "literature"
    source: str = ""
    evidence_available_date: date | None = None
    independent_case_series_count: int = 0
    has_omim_clinical_synopsis: bool = False
    cns_penetrant: bool = False
    has_dominant_target: bool = False
    label_event_frequency: float | None = None  # fraction in [0, 1] when the label reports it
    predication_type: str | None = None


def assign_evidence_grade(record: EvidenceRecord) -> EvidenceGrade:
    if record.evidence_class == "monogenic":
        if record.independent_case_series_count >= 2 or record.has_omim_clinical_synopsis:
            return "A"
        return "C"
    if record.evidence_class == "pharmacological":
        if record.cns_penetrant and record.has_dominant_target:
            return "B"
        return "E"  # label event without a clean molecular entry point: soft prior at most
    if record.evidence_class == "human_association":
        return "C"
    if record.evidence_class == "animal":
        return "D"
    if record.evidence_class == "literature":
        return "E"
    raise ValueError(f"unknown evidence_class: {record.evidence_class!r}")


def loss_weight_for_record(
    record: EvidenceRecord,
    grade_weights: dict[EvidenceGrade, float] | None = None,
    predication_weights: dict[str, float] | None = None,
) -> float:
    """Loss weight for one observation; zero means the observation is not trained on."""
    grade_weights = grade_weights or DEFAULT_GRADE_WEIGHTS
    predication_weights = predication_weights or DEFAULT_PREDICATION_WEIGHTS
    grade = assign_evidence_grade(record)
    base_weight = grade_weights[grade]
    if grade == "B" and record.label_event_frequency is not None:
        # frequency scaling keeps common label events heavier than rare ones; floor avoids zeroing rare but real events
        return base_weight * max(0.25, min(1.0, record.label_event_frequency * 10.0))
    if grade == "E" and record.predication_type is not None:
        return predication_weights.get(record.predication_type.upper(), predication_weights["COMENTION"])
    return base_weight
