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
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


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


def load_sider_frequencies(sider_directory: Path) -> dict[tuple[str, str], float]:
    """(STITCH flat id, preferred term lowercase) -> mean of reported range midpoints, placebo rows excluded."""
    midpoints_by_pair: dict[tuple[str, str], list[float]] = defaultdict(list)
    with gzip.open(sider_directory / "meddra_freq.tsv.gz", "rt", encoding="utf-8") as frequency_file:
        for row in csv.reader(frequency_file, delimiter="\t"):
            stitch_flat_id, placebo_flag, lower, upper, concept_type, term_name = row[0], row[3], row[5], row[6], row[7], row[9]
            if concept_type != "PT" or placebo_flag == "placebo":
                continue
            try:
                midpoints_by_pair[(stitch_flat_id, term_name.lower())].append((float(lower) + float(upper)) / 2.0)
            except ValueError:
                continue
    return {pair: sum(values) / len(values) for pair, values in midpoints_by_pair.items()}


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
    frequencies = load_sider_frequencies(sider_directory)
    events: list[DrugLabelEvent] = []
    seen_pairs: set[tuple[str, str, str]] = set()
    with gzip.open(sider_directory / "meddra_all_se.tsv.gz", "rt", encoding="utf-8") as side_effect_file:
        for row in csv.reader(side_effect_file, delimiter="\t"):
            stitch_flat_id, concept_type, term_name = row[0], row[3], row[5]
            if concept_type != "PT":
                continue
            target_symptom = preferred_term_to_symptom.get(term_name.lower())
            if target_symptom is None and not keep_unmapped_terms:
                continue
            pair_key = (stitch_flat_id, term_name.lower(), "induces")
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            events.append(DrugLabelEvent(stitch_flat_id, drug_names.get(stitch_flat_id, stitch_flat_id), pubchem_cid_from_stitch_flat(stitch_flat_id), term_name, "induces", target_symptom, frequencies.get((stitch_flat_id, term_name.lower())), atc_codes.get(stitch_flat_id, [])))
    with gzip.open(sider_directory / "meddra_all_indications.tsv.gz", "rt", encoding="utf-8") as indication_file:
        for row in csv.reader(indication_file, delimiter="\t"):
            stitch_flat_id, detection_method, concept_type, term_name = row[0], row[2], row[4], row[6]
            if concept_type != "PT" or detection_method not in indication_methods:
                continue
            target_symptom = preferred_term_to_symptom.get(term_name.lower())
            if target_symptom is None and not keep_unmapped_terms:
                continue
            pair_key = (stitch_flat_id, term_name.lower(), "relieves")
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            events.append(DrugLabelEvent(stitch_flat_id, drug_names.get(stitch_flat_id, stitch_flat_id), pubchem_cid_from_stitch_flat(stitch_flat_id), term_name, "relieves", target_symptom, None, atc_codes.get(stitch_flat_id, []), indication_detection_method=detection_method))
    return events


def load_onsides_events(onsides_directory: Path, crosswalk_path: Path):
    raise NotImplementedError("pin the OnSIDES release and parse its tables; keep section and confidence columns")


def is_nervous_system_atc(atc_codes: list[str]) -> bool:
    """Crude version 1 proxy for a centrally acting drug: any ATC code in anatomical group N.

    This is not blood-brain barrier penetration; it is replaced by a measured BBB flag when one is added.
    """
    return any(code.startswith("N") for code in atc_codes)
