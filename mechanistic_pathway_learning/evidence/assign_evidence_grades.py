"""Evidence grading (experiment design, section 4.3).

Every perturbation-symptom observation carries an evidence record. The grade
decides the role (training label, soft prior or validation only) and the loss
weight. Grades are deliberately coarse; the weights are the initial values from
the design and are overridden by configs/evidence_assembly.yaml.

Grade A  human loss of function in a graph gene with the symptom in a curated clinical
         synopsis (OMIM or Orphanet provenance of the HPO annotation) or documented in at
         least two independent case series
Grade B  CNS-penetrant drug with a dominant target and the event on the label
Grade C  human association (gene-level GWAS statistic, metabolomic association,
         or a monogenic report below the grade A bar); validation only in version 1
Grade D  animal perturbation with a symptom analogue; soft prior
Grade E  literature predication or co-mention; soft prior, weight by predication type

Grade A policies. "curated_synopsis_provenance" (default) accepts any HPO annotation with an
OMIM or Orphanet disease entry, because both are curated clinical synopses, and lets the
frequency qualifier set the weight. "two_distinct_disease_entries" is the version 0.3 proxy
(two distinct disease identifiers or an OMIM entry); it is kept for comparison because it
counted an OMIM entry and its Orphanet counterpart as independent case series and zeroed
single-entry genes such as HMBS (acute intermittent porphyria).

Frequency scaling. Grade A weights are scaled by the largest HPO frequency reported for the
pair, and grade B weights by the label frequency, with the same shape: weight = base *
clip(frequency * scale, floor, 1). Unknown frequency keeps the base weight.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

EvidenceGrade = str  # "A", "B", "C", "D" or "E"

GRADE_A_POLICIES = ("curated_synopsis_provenance", "two_distinct_disease_entries")
DEFAULT_GRADE_A_POLICY = "curated_synopsis_provenance"

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

# HPO frequency 25 percent or more keeps the full weight; Occasional (midpoint 0.17) -> 0.68; Very rare -> floor
MONOGENIC_FREQUENCY_SCALE = 4.0
MONOGENIC_FREQUENCY_FLOOR = 0.25
# label frequency 10 percent or more keeps the full weight; floor avoids zeroing rare but real events
LABEL_FREQUENCY_SCALE = 10.0
LABEL_FREQUENCY_FLOOR = 0.25


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
    disease_identifiers: list[str] = field(default_factory=list)  # OMIM:/ORPHA: entries behind a monogenic record
    omim_entry_count: int = 0
    orpha_entry_count: int = 0
    annotation_row_count: int = 0
    max_annotation_frequency: float | None = None  # largest HPO frequency midpoint across the pair's annotations
    annotation_patient_count: int | None = None  # patients behind n/m fractions, summed, when reported
    distinct_reference_count: int = 0  # distinct reference strings (PMID:, OMIM:, ORPHA:, ...) across the pair's (disease, term) annotations in phenotype.hpoa
    distinct_pubmed_reference_count: int = 0  # the PMID: subset of those; descriptive only, never enters the weight
    association_type_known: bool = False  # True when at least one positive report has a typed gene-disease association (genes_to_disease.txt or Orphadata product 6)
    causal_association_count: int = 0  # positive reports whose association is disease-causing (MENDELIAN or an Orphanet disease-causing type)
    evidence_date_source: str = ""  # source of evidence_available_date: publication, biocuration, release or empty


def frequency_scaled_weight(base_weight: float, frequency: float | None, scale: float, floor: float) -> float:
    if frequency is None:
        return base_weight
    return base_weight * max(floor, min(1.0, frequency * scale))


def assign_evidence_grade(record: EvidenceRecord, grade_a_policy: str = DEFAULT_GRADE_A_POLICY) -> EvidenceGrade:
    if grade_a_policy not in GRADE_A_POLICIES:
        raise ValueError(f"unknown grade_a_policy {grade_a_policy!r}; choose from {GRADE_A_POLICIES}")
    if record.evidence_class == "monogenic":
        if record.association_type_known and record.causal_association_count == 0:
            return "C"  # polygenic locus, susceptibility factor, candidate, modifier or biomarker only: a human association, not a loss of function (review v0.4, finding 4)
        if record.independent_case_series_count >= 2 or record.has_omim_clinical_synopsis:
            return "A"
        if grade_a_policy == "curated_synopsis_provenance" and (record.omim_entry_count > 0 or record.orpha_entry_count > 0):
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
    grade_a_policy: str = DEFAULT_GRADE_A_POLICY,
    scale_by_frequency: bool = True,
) -> float:
    """Loss weight for one observation; zero means the observation is not trained on."""
    grade_weights = grade_weights or DEFAULT_GRADE_WEIGHTS
    predication_weights = predication_weights or DEFAULT_PREDICATION_WEIGHTS
    grade = assign_evidence_grade(record, grade_a_policy)
    base_weight = grade_weights[grade]
    if grade == "A" and record.evidence_class == "monogenic" and scale_by_frequency:
        return frequency_scaled_weight(base_weight, record.max_annotation_frequency, MONOGENIC_FREQUENCY_SCALE, MONOGENIC_FREQUENCY_FLOOR)
    if grade == "B" and scale_by_frequency:
        return frequency_scaled_weight(base_weight, record.label_event_frequency, LABEL_FREQUENCY_SCALE, LABEL_FREQUENCY_FLOOR)
    if grade == "E" and record.predication_type is not None:
        return predication_weights.get(record.predication_type.upper(), predication_weights["COMENTION"])
    return base_weight
