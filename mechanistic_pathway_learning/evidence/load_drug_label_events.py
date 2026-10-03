"""Evidence class E2: drug-label events from SIDER 4.1 and OnSIDES (design section 4.2).

SIDER 4.1 (http://sideeffects.embl.de/media/download/), labels current to 2015:
  meddra_all_se.tsv.gz        STITCH_flat, STITCH_stereo, UMLS_label_cui, concept_type (LLT or PT),
                              UMLS_meddra_cui, side_effect_name
  meddra_freq.tsv.gz          STITCH_flat, STITCH_stereo, UMLS_label_cui, placebo_flag, frequency_text,
                              frequency_lower, frequency_upper, concept_type, UMLS_meddra_cui, side_effect_name
  meddra_all_indications.tsv.gz  STITCH_flat, UMLS_label_cui, detection_method, label_concept_name,
                              concept_type, UMLS_meddra_cui, meddra_concept_name
  drug_names.tsv              STITCH_flat, drug_name
  drug_atc.tsv                STITCH_flat, ATC_code (one row per code)
A STITCH flat identifier CID1xxxxxxxx encodes PubChem CID xxxxxxxx (the leading 1 marks the
flat, stereo-merged form); strip "CID1" and leading zeros to recover the PubChem CID.

OnSIDES (https://github.com/tatonetti-lab/onsides), current labels: loader still to be written
against the pinned release schema; keep the two sources in separate columns because the time
split depends on them.

MedDRA preferred-term names are licensed: only the mapping from preferred term to target
symptom (docs/symptom_crosswalk.csv) is committed, never the SIDER tables themselves.

Reports. meddra_freq.tsv.gz has one row per (label, frequency statement), treatment and placebo arms
both; a (drug, preferred term) pair may carry several treatment rows. drug_label_reports turns one
included DrugLabelEvent into one EvidenceReport per treatment-arm frequency row, or one label_text
report when the label gives no frequency, and one label_text report per indication (evidence_reports.py).
The event's label_frequency, the mean of the treatment midpoints, is unchanged and still sets the grade
weight.
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from mechanistic_pathway_learning.evidence.evidence_reports import EvidenceReport, limitations_text


@dataclass
class SiderFrequencyRow:
    """One row of meddra_freq.tsv.gz at the preferred-term level; ordinal is its 1-based position within its
    (drug, preferred term, arm) group in file order, so treatment rows are numbered apart from placebo rows."""

    stitch_flat_id: str
    label_cui: str
    placebo_flag: str  # "" for the treatment arm, "placebo" for the placebo arm
    frequency_text: str
    frequency_lower: float
    frequency_upper: float
    meddra_cui: str
    preferred_term: str
    ordinal: int

    @property
    def midpoint(self) -> float:
        return (self.frequency_lower + self.frequency_upper) / 2.0

    @property
    def is_placebo(self) -> bool:
        return self.placebo_flag == "placebo"


@dataclass
class DrugLabelEvent:
    stitch_flat_id: str
    drug_name: str
    pubchem_cid: int
    meddra_preferred_term: str
    relation: str  # "induces" for adverse events, "relieves" for indications
    target_symptom: str | None = None
    label_frequency: float | None = None  # midpoint of the reported range when SIDER has it
    atc_codes: list[str] = field(default_factory=list)
    source: str = "SIDER 4.1"
    indication_detection_method: str | None = None  # SIDER: NLP_indication, NLP_precondition or text_mention
    meddra_cui: str | None = None  # UMLS concept id of the preferred term
    frequency_rows: list[SiderFrequencyRow] = field(default_factory=list)  # every meddra_freq row of the pair, both arms


def pubchem_cid_from_stitch_flat(stitch_flat_id: str) -> int:
    return int(stitch_flat_id[4:])


def read_two_column_tsv(path: Path) -> dict[str, list[str]]:
    values_by_key: dict[str, list[str]] = defaultdict(list)
    with open(path, encoding="utf-8") as tsv_file:
        for row in csv.reader(tsv_file, delimiter="\t"):
            if len(row) >= 2:
                values_by_key[row[0]].append(row[1])
    return values_by_key


def read_preferred_term_to_target_symptom(crosswalk_path: Path) -> dict[str, str]:
    preferred_term_to_symptom: dict[str, str] = {}
    with open(crosswalk_path, encoding="utf-8") as crosswalk_file:
        for row in csv.DictReader(crosswalk_file):
            for preferred_term in (row.get("meddra_preferred_terms") or "").split(";"):
                if preferred_term.strip():
                    preferred_term_to_symptom[preferred_term.strip().lower()] = row["target_symptom"]
    return preferred_term_to_symptom


def load_sider_frequency_rows(sider_directory: Path) -> dict[tuple[str, str], list[SiderFrequencyRow]]:
    """(STITCH flat id, preferred term lowercase) -> its meddra_freq rows at the preferred-term level, both arms, file order.

    Rows whose bounds do not parse as numbers are skipped, as load_sider_frequencies always did.
    """
    rows_by_pair: dict[tuple[str, str], list[SiderFrequencyRow]] = defaultdict(list)
    ordinal_by_pair_and_arm: dict[tuple[str, str, str], int] = defaultdict(int)
    with gzip.open(sider_directory / "meddra_freq.tsv.gz", "rt", encoding="utf-8") as frequency_file:
        for row in csv.reader(frequency_file, delimiter="\t"):
            stitch_flat_id, label_cui, placebo_flag, frequency_text, lower, upper, concept_type, meddra_cui, term_name = row[0], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]
            if concept_type != "PT":
                continue
            try:
                lower_value, upper_value = float(lower), float(upper)
            except ValueError:
                continue
            pair = (stitch_flat_id, term_name.lower())
            ordinal_by_pair_and_arm[(*pair, placebo_flag)] += 1
            rows_by_pair[pair].append(SiderFrequencyRow(stitch_flat_id, label_cui, placebo_flag, frequency_text, lower_value, upper_value, meddra_cui, term_name, ordinal_by_pair_and_arm[(*pair, placebo_flag)]))
    return dict(rows_by_pair)


def mean_treatment_midpoint(frequency_rows: list[SiderFrequencyRow]) -> float | None:
    """Mean of the treatment-arm range midpoints in file order (the label frequency of the grade table), None without any."""
    midpoints = [frequency_row.midpoint for frequency_row in frequency_rows if not frequency_row.is_placebo]
    return sum(midpoints) / len(midpoints) if midpoints else None


def load_sider_frequencies(sider_directory: Path) -> dict[tuple[str, str], float]:
    """(STITCH flat id, preferred term lowercase) -> mean of reported range midpoints, placebo rows excluded."""
    frequencies: dict[tuple[str, str], float] = {}
    for pair, frequency_rows in load_sider_frequency_rows(sider_directory).items():
        mean_midpoint = mean_treatment_midpoint(frequency_rows)
        if mean_midpoint is not None:
            frequencies[pair] = mean_midpoint
    return frequencies


def load_sider_events(
    sider_directory: Path,
    crosswalk_path: Path,
    keep_unmapped_terms: bool = False,
    indication_methods: tuple[str, ...] = ("NLP_indication",),
) -> list[DrugLabelEvent]:
    """Adverse events (relation induces) and indications (relation relieves) at the preferred-term level.

    Indications are kept only when SIDER's detection method is in indication_methods; the default keeps
    NLP_indication and drops NLP_precondition and text_mention, which include mentions of the condition
    in label text that are not indications.
    """
    drug_names = {key: values[0] for key, values in read_two_column_tsv(sider_directory / "drug_names.tsv").items()}
    atc_codes = read_two_column_tsv(sider_directory / "drug_atc.tsv")
    preferred_term_to_symptom = read_preferred_term_to_target_symptom(crosswalk_path)
    frequency_rows_by_pair = load_sider_frequency_rows(sider_directory)
    events: list[DrugLabelEvent] = []
    seen_pairs: set[tuple[str, str, str]] = set()
    with gzip.open(sider_directory / "meddra_all_se.tsv.gz", "rt", encoding="utf-8") as side_effect_file:
        for row in csv.reader(side_effect_file, delimiter="\t"):
            stitch_flat_id, concept_type, meddra_cui, term_name = row[0], row[3], row[4], row[5]
            if concept_type != "PT":
                continue
            target_symptom = preferred_term_to_symptom.get(term_name.lower())
            if target_symptom is None and not keep_unmapped_terms:
                continue
            pair_key = (stitch_flat_id, term_name.lower(), "induces")
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            frequency_rows = frequency_rows_by_pair.get((stitch_flat_id, term_name.lower()), [])
            events.append(DrugLabelEvent(stitch_flat_id, drug_names.get(stitch_flat_id, stitch_flat_id), pubchem_cid_from_stitch_flat(stitch_flat_id), term_name, "induces", target_symptom, mean_treatment_midpoint(frequency_rows), atc_codes.get(stitch_flat_id, []),
                                         meddra_cui=meddra_cui, frequency_rows=frequency_rows))
    with gzip.open(sider_directory / "meddra_all_indications.tsv.gz", "rt", encoding="utf-8") as indication_file:
        for row in csv.reader(indication_file, delimiter="\t"):
            stitch_flat_id, detection_method, concept_type, meddra_cui, term_name = row[0], row[2], row[4], row[5], row[6]
            if concept_type != "PT" or detection_method not in indication_methods:
                continue
            target_symptom = preferred_term_to_symptom.get(term_name.lower())
            if target_symptom is None and not keep_unmapped_terms:
                continue
            pair_key = (stitch_flat_id, term_name.lower(), "relieves")
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            events.append(DrugLabelEvent(stitch_flat_id, drug_names.get(stitch_flat_id, stitch_flat_id), pubchem_cid_from_stitch_flat(stitch_flat_id), term_name, "relieves", target_symptom, None, atc_codes.get(stitch_flat_id, []), indication_detection_method=detection_method, meddra_cui=meddra_cui))
    return events


def load_onsides_events(onsides_directory: Path, crosswalk_path: Path):
    raise NotImplementedError("pin the OnSIDES release and parse its tables; keep section and confidence columns")


def is_nervous_system_atc(atc_codes: list[str]) -> bool:
    """Crude version 1 proxy for a centrally acting drug: any ATC code in anatomical group N.

    This is not blood-brain barrier penetration; it is replaced by a measured BBB flag when one is added.
    """
    return any(code.startswith("N") for code in atc_codes)


def drug_label_reports(event: DrugLabelEvent, model_description: str, perturbation_nodes_json: str) -> list[EvidenceReport]:
    """Reports of one included event: one per treatment-arm frequency row, else one label_text report; indications are label_text.

    SIDER records no absences, so report_value is 1 throughout. A treatment row whose pair also has a placebo
    row is evidence_code placebo_controlled_frequency (placebo_flag true), otherwise frequency; the placebo-arm
    frequency itself is not a column of the report table.
    """
    if event.target_symptom is None:
        return []
    source = "SIDER-indication" if event.relation == "relieves" else "SIDER-label"
    meddra_cui = event.meddra_cui or ""
    treatment_rows = [frequency_row for frequency_row in event.frequency_rows if not frequency_row.is_placebo]
    has_placebo_row = any(frequency_row.is_placebo for frequency_row in event.frequency_rows)

    def report(evidence_code: str, frequency: float | None, ordinal: int) -> EvidenceReport:
        placebo_controlled = evidence_code == "placebo_controlled_frequency"
        built = EvidenceReport(
            report_id=f"{source}|{event.stitch_flat_id}|{event.stitch_flat_id}|{meddra_cui}|{event.target_symptom}|{ordinal}",
            perturbation_id=event.stitch_flat_id, perturbation_type="drug", perturbation_label=event.drug_name,
            symptom=event.target_symptom, relation=event.relation, evidence_class="pharmacological", source=source,
            report_value=1, source_record_id=event.stitch_flat_id, source_record_label=event.drug_name,
            source_term_id=meddra_cui, source_term_label=event.meddra_preferred_term, evidence_code=evidence_code,
            model_description=model_description, frequency=frequency, frequency_denominator=None, placebo_flag=placebo_controlled,
            onset=None, sex=None, references="", pubmed_reference_count=0, evidence_date=None,
            rubric_log_sample_size=0.0, rubric_frequency_known=1.0 if frequency is not None else 0.0,
            rubric_evidence_code_pcs=0.0, rubric_evidence_code_tas=0.0, rubric_evidence_code_iea=0.0,
            rubric_placebo_controlled=1.0 if placebo_controlled else 0.0, rubric_curated_synopsis=0.0,
            perturbation_nodes=perturbation_nodes_json,
        )
        built.limitations = limitations_text(built)
        return built

    if event.relation == "induces" and treatment_rows:
        return [report("placebo_controlled_frequency" if has_placebo_row else "frequency", frequency_row.midpoint, frequency_row.ordinal) for frequency_row in treatment_rows]
    return [report("label_text", None, 0)]
