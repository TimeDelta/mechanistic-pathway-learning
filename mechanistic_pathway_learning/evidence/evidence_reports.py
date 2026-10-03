"""The report table: one row per piece of evidence behind a perturbation-symptom pair (design sections 4.2 and 4.3).

A report is one phenotype.hpoa annotation row behind one (gene, target symptom) or one SIDER frequency
statement, label mention or indication behind one (drug, target symptom). Keeping the rows apart, rather
than one aggregated row per pair, is what lets evidence be weighted differently by the model it rests on
(human loss of function in a named disease entry, a drug label), by its evidence code and sample size,
and lets evidence that asserts absence (the NOT qualifier, the Excluded frequency) enter as a report with
value 0 instead of disappearing or, as before, counting as a positive. evidence_records.parquet
aggregates this table to one row per pair (assemble_evidence_table.aggregate_reports_to_observations).

Two things in this module are fixed functions of the columns and never learned:

  rubric_weight     the weight of one report inside the reliability likelihood, in (0, 1]:
                    clip(evidence_code_factor x sample_size_factor, minimum_rubric_weight, 1). It stands in
                    for the feature-dependent sensitivity of Raykar et al. until a third source makes that
                    form estimable; the defaults are overridable from a JSON file.
  limitations_text  the clauses that describe what the report cannot show, for display on mechanism cards;
                    it never enters a weight.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_NAMES: tuple[str, ...] = ("HPO-OMIM", "HPO-Orphanet", "SIDER-label", "SIDER-indication")
HPO_SOURCES: tuple[str, ...] = ("HPO-OMIM", "HPO-Orphanet")
SIDER_SOURCES: tuple[str, ...] = ("SIDER-label", "SIDER-indication")
UNJOINED_EVIDENCE_CODE = "unjoined"
SIDER_FREQUENCY_EVIDENCE_CODES: tuple[str, ...] = ("placebo_controlled_frequency", "frequency")
SMALL_SAMPLE_PATIENT_COUNT = 20


@dataclass
class EvidenceReport:
    """One report; one field per column of evidence_reports.parquet, in column order."""

    report_id: str
    perturbation_id: str
    perturbation_type: str  # "gene" or "drug"
    perturbation_label: str
    symptom: str
    relation: str  # "induces" or "relieves"
    evidence_class: str  # "monogenic" or "pharmacological"
    source: str  # one of SOURCE_NAMES
    report_value: int  # 1 asserts presence, 0 asserts absence
    source_record_id: str  # disease id or STITCH flat id
    source_record_label: str
    source_term_id: str  # HPO term id or MedDRA concept id
    source_term_label: str
    evidence_code: str  # PCS, TAS, IEA, other verbatim codes, unjoined; placebo_controlled_frequency, frequency, label_text
    model_description: str
    frequency: float | None
    frequency_denominator: int | None
    placebo_flag: bool
    onset: str | None
    sex: str | None
    references: str
    pubmed_reference_count: int
    evidence_date: str | None  # ISO date
    rubric_log_sample_size: float
    rubric_frequency_known: float
    rubric_evidence_code_pcs: float
    rubric_evidence_code_tas: float
    rubric_evidence_code_iea: float
    rubric_placebo_controlled: float
    rubric_curated_synopsis: float
    limitations: str = ""
    perturbation_nodes: str = ""  # JSON list of [node_id, sign, magnitude]; filled by the assembler


REPORT_COLUMNS: tuple[str, ...] = tuple(report_field.name for report_field in fields(EvidenceReport))


def reports_to_dataframe(reports: list[EvidenceReport]) -> pd.DataFrame:
    """Reports as a table with the columns of REPORT_COLUMNS in order and stable dtypes, empty-safe."""
    table = pd.DataFrame([asdict(report) for report in reports], columns=list(REPORT_COLUMNS))
    table["report_value"] = table.report_value.astype("int64")
    table["pubmed_reference_count"] = table.pubmed_reference_count.astype("int64")
    table["frequency"] = pd.to_numeric(table.frequency, errors="coerce").astype("float64")
    table["frequency_denominator"] = pd.array([None if value is None or (isinstance(value, float) and math.isnan(value)) else int(value) for value in table.frequency_denominator], dtype="Int64")
    table["placebo_flag"] = table.placebo_flag.astype(bool)
    for column in ("rubric_log_sample_size", "rubric_frequency_known", "rubric_evidence_code_pcs", "rubric_evidence_code_tas", "rubric_evidence_code_iea", "rubric_placebo_controlled", "rubric_curated_synopsis"):
        table[column] = table[column].astype("float64")
    return table


def limitations_text(report: EvidenceReport) -> str:
    """The clauses that apply to one report, joined with "; "; for display, never for a weight."""
    clauses: list[str] = []
    if report.source in HPO_SOURCES:
        if report.frequency is None:
            clauses.append("no frequency reported")
        elif report.frequency_denominator is None:
            clauses.append("frequency from a curated qualifier, no patient count")
        if report.frequency_denominator is not None and report.frequency_denominator < SMALL_SAMPLE_PATIENT_COUNT:
            clauses.append(f"n = {report.frequency_denominator} patients")
        if report.evidence_code == "IEA":
            clauses.append("inferred from electronic annotation (IEA)")
        if report.evidence_code == "TAS":
            clauses.append("traceable author statement, no primary study cited (TAS)")
        if report.pubmed_reference_count == 0:
            clauses.append("no PubMed reference")
        if report.source == "HPO-Orphanet":
            clauses.append("annotation imported from the Orphanet clinical-sign table, undated")
        if report.evidence_code == UNJOINED_EVIDENCE_CODE:
            clauses.append("no phenotype.hpoa row joined")
    else:
        if report.evidence_code == "label_text":
            clauses.append("label text only, no frequency")
        if report.evidence_code == "frequency":
            clauses.append("frequency not placebo-controlled")
        clauses.append("SIDER reports no patient counts")
        if report.evidence_code in SIDER_FREQUENCY_EVIDENCE_CODES:
            clauses.append("label frequencies are not comparable across labels")
    return "; ".join(clauses)


DEFAULT_EVIDENCE_CODE_FACTORS: dict[str, float] = {
    "PCS": 1.0,  # published clinical study
    "TAS": 0.8,  # traceable author statement
    "IEA": 0.6,  # inferred from electronic annotation
    "other": 0.7,  # any other code, and rows with no phenotype.hpoa row joined
    "placebo_controlled_frequency": 1.0,
    "frequency": 0.9,
    "label_text": 0.8,
}


@dataclass(frozen=True)
class RubricWeightDefaults:
    """Parameters of the fixed rubric weight; every field can be overridden from a JSON file.

    implicit_negative_weight is the weight of a source's silence on a covered perturbation inside the
    reliability likelihood; it lives here because silence inside a curated synopsis is weaker evidence of
    absence than a NOT annotation, and it is the parameter to revisit first once a third source exists.
    """

    evidence_code_factors: dict[str, float] = field(default_factory=lambda: dict(DEFAULT_EVIDENCE_CODE_FACTORS))
    sample_size_reference: int = 20
    sample_size_unknown_factor: float = 1.0
    minimum_rubric_weight: float = 0.1
    implicit_negative_weight: float = 1.0

    def as_dict(self) -> dict:
        return asdict(self)


def load_rubric_weight_defaults(overrides_path: Path | None) -> RubricWeightDefaults:
    """Defaults, with any subset of the fields replaced from a JSON file; unknown keys raise ValueError.

    evidence_code_factors in the file is merged into the default factors, so one factor can be changed alone.
    """
    if overrides_path is None:
        return RubricWeightDefaults()
    overrides = json.loads(Path(overrides_path).read_text())
    if not isinstance(overrides, dict):
        raise ValueError(f"{overrides_path}: expected a JSON object with rubric weight fields")
    known_fields = {default_field.name for default_field in fields(RubricWeightDefaults)}
    unknown_keys = sorted(set(overrides) - known_fields)
    if unknown_keys:
        raise ValueError(f"{overrides_path}: unknown rubric weight fields {unknown_keys}; known fields are {sorted(known_fields)}")
    factors = dict(DEFAULT_EVIDENCE_CODE_FACTORS)
    factors.update({str(code): float(value) for code, value in (overrides.get("evidence_code_factors") or {}).items()})
    keyword_arguments = {key: value for key, value in overrides.items() if key != "evidence_code_factors"}
    return RubricWeightDefaults(evidence_code_factors=factors, **keyword_arguments)


def evidence_code_factor(evidence_code: str, defaults: RubricWeightDefaults) -> float:
    return float(defaults.evidence_code_factors.get(evidence_code, defaults.evidence_code_factors["other"]))


def sample_size_factor(frequency_denominator: int | float | None, defaults: RubricWeightDefaults) -> float:
    """1 when no patient count is reported (a curated qualifier or a label frequency is not evidence of a small sample);
    otherwise 0.5 + 0.5 x min(1, log1p(n) / log1p(reference)), so n = 1 gives 0.61, n = 5 gives 0.79 and n >= 20 gives 1."""
    if frequency_denominator is None or (isinstance(frequency_denominator, float) and math.isnan(frequency_denominator)):
        return float(defaults.sample_size_unknown_factor)
    return 0.5 + 0.5 * min(1.0, math.log1p(float(frequency_denominator)) / math.log1p(float(defaults.sample_size_reference)))


def rubric_weight(report: EvidenceReport, defaults: RubricWeightDefaults) -> float:
    """clip(evidence_code_factor x sample_size_factor, minimum_rubric_weight, 1) for one report."""
    product = evidence_code_factor(report.evidence_code, defaults) * sample_size_factor(report.frequency_denominator, defaults)
    return float(min(1.0, max(defaults.minimum_rubric_weight, product)))


def rubric_weights_for_table(reports: pd.DataFrame, defaults: RubricWeightDefaults) -> np.ndarray:
    """rubric_weight for every row of a report table, in row order."""
    code_factors = np.array([evidence_code_factor(str(code), defaults) for code in reports.evidence_code], dtype=float)
    denominators = reports.frequency_denominator.to_numpy(dtype=float, na_value=np.nan) if len(reports) else np.zeros(0)
    size_factors = np.array([sample_size_factor(None if np.isnan(value) else value, defaults) for value in denominators], dtype=float)
    return np.clip(code_factors * size_factors, defaults.minimum_rubric_weight, 1.0)
