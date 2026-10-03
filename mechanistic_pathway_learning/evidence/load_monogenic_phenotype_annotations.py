"""Evidence class E1: monogenic gene-to-phenotype annotations from HPO (design section 4.2).

Input files (https://hpo.jax.org, annotation downloads):
  genes_to_phenotype.txt  tab-separated with a header line; columns include
                          ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id
  hp.obo                  the ontology, used to expand each target symptom to its descendants

Expansion. Each target symptom in docs/symptom_crosswalk.csv names root HPO terms (hpo_ids) and,
optionally, descendant terms to exclude together with their own descendants (excluded_hpo_ids).
Exclusions exist because the is_a closure of a symptom term can contain diagnosis-level terms
(Attention deficit hyperactivity disorder under Cognitive impairment) or a different construct
(Restless legs under Restlessness); experiments/audit_hpo_term_expansion.py measures what every
descendant contributes so each exclusion is auditable (design assumption A7).

Frequency. The frequency column holds an HPO frequency qualifier (HP:0040280 Obligate through
HP:0040285 Excluded), a patient fraction such as 3/5, a percentage, or "-". It is parsed to a
midpoint in [0, 1]; rows qualified as Excluded assert the absence of the feature and are dropped.

Provenance. disease_id is an OMIM or Orphanet entry. HPO annotations to OMIM diseases are
biocurated from the OMIM clinical synopsis and Orphanet annotations from Orphanet's expert
clinical-sign tables, so provenance is recorded per record (omim_entry_count, orpha_entry_count)
and the grade A policy in assign_evidence_grades can be applied without a separate OMIM lookup.
The two-distinct-entries proxy used in version 0.3 counted an OMIM entry and its Orphanet
counterpart as two independent case series; it is kept only as an alternative policy.

Dates. phenotype.hpoa (the disease-level annotation file of the same release) carries a biocuration
field with curator and date per (disease, term) annotation. OMIM-sourced annotations are dated from
2009 onward; every Orphanet annotation carries the import date of the release, which is no date at
all. A (gene, symptom) pair therefore gets evidence_available_date = the earliest biocuration date
among the OMIM (disease, term) annotations behind it, and no date when only Orphanet annotations
support it. The date is when the HPO team recorded the annotation, an upper bound on when the
observation was published (the reference column holds the PMID for most OMIM rows), which the time
split of design section 6.1 states as its caveat.

Output: one EvidenceRecord per (gene, target symptom) with evidence_class "monogenic".
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from datetime import date
from pathlib import Path

from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord

FREQUENCY_QUALIFIER_MIDPOINTS: dict[str, float] = {
    "HP:0040280": 1.0,  # Obligate, 100%
    "HP:0040281": 0.895,  # Very frequent, 80-99%
    "HP:0040282": 0.545,  # Frequent, 30-79%
    "HP:0040283": 0.17,  # Occasional, 5-29%
    "HP:0040284": 0.025,  # Very rare, 1-4%
    "HP:0040285": 0.0,  # Excluded, 0%
}
EXCLUDED_FREQUENCY_QUALIFIER = "HP:0040285"


def load_hpo_is_a_parents_from_obo(obo_path: Path) -> dict[str, set[str]]:
    """Return term id -> set of direct parent ids from an OBO file (is_a lines only)."""
    parents_by_term: dict[str, set[str]] = defaultdict(set)
    current_term: str | None = None
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
            elif line.startswith("id: HP:"):
                current_term = line[len("id: "):]
                parents_by_term.setdefault(current_term, set())
            elif line.startswith("is_a:") and current_term is not None:
                parent_id = line[len("is_a:"):].strip().split(" ")[0]
                parents_by_term[current_term].add(parent_id)
    return dict(parents_by_term)


def load_hpo_term_names_from_obo(obo_path: Path) -> dict[str, str]:
    """Return term id -> name from an OBO file."""
    names_by_term: dict[str, str] = {}
    current_term: str | None = None
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
            elif line.startswith("id: HP:"):
                current_term = line[len("id: "):]
            elif line.startswith("name:") and current_term is not None:
                names_by_term[current_term] = line[len("name:"):].strip()
    return names_by_term


def descendants_of(term_id: str, parents_by_term: Mapping[str, set[str]]) -> set[str]:
    """All terms whose is_a ancestry includes term_id, plus term_id itself."""
    children_by_term: dict[str, set[str]] = defaultdict(set)
    for child, parents in parents_by_term.items():
        for parent in parents:
            children_by_term[parent].add(child)
    collected = {term_id}
    frontier = [term_id]
    while frontier:
        current = frontier.pop()
        for child in children_by_term.get(current, ()):
            if child not in collected:
                collected.add(child)
                frontier.append(child)
    return collected


def read_crosswalk_hpo_terms(crosswalk_path: Path) -> dict[str, tuple[list[str], list[str]]]:
    """Return target symptom -> (root HPO ids, excluded descendant HPO ids) from the crosswalk."""
    terms_by_symptom: dict[str, tuple[list[str], list[str]]] = {}
    with open(crosswalk_path, encoding="utf-8") as crosswalk_file:
        for row in csv.DictReader(crosswalk_file):
            root_ids = [term.strip() for term in (row.get("hpo_ids") or "").split(";") if term.strip()]
            excluded_ids = [term.strip() for term in (row.get("excluded_hpo_ids") or "").split(";") if term.strip()]
            terms_by_symptom[row["target_symptom"]] = (root_ids, excluded_ids)
    return terms_by_symptom


def expand_symptom_terms(
    target_symptom_to_hpo_ids: Mapping[str, Iterable[str]],
    parents_by_term: Mapping[str, set[str]],
    excluded_hpo_ids_by_symptom: Mapping[str, Iterable[str]] | None = None,
) -> dict[str, list[str]]:
    """HPO term id -> target symptoms whose expansion contains it, after removing excluded subtrees."""
    symptoms_by_term: dict[str, list[str]] = defaultdict(list)
    for target_symptom, root_ids in target_symptom_to_hpo_ids.items():
        excluded_terms: set[str] = set()
        for excluded_id in (excluded_hpo_ids_by_symptom or {}).get(target_symptom, ()):
            excluded_terms |= descendants_of(excluded_id, parents_by_term)
        for root_id in root_ids:
            for descendant_id in descendants_of(root_id, parents_by_term):
                if descendant_id not in excluded_terms and target_symptom not in symptoms_by_term[descendant_id]:
                    symptoms_by_term[descendant_id].append(target_symptom)
    return dict(symptoms_by_term)


def parse_frequency_qualifier(raw_frequency: str | None) -> float | None:
    """HPO frequency column -> midpoint fraction in [0, 1], or None when not reported ("-" or empty).

    Accepts the HP:00402xx qualifiers, patient fractions such as "3/5" and percentages such as "25%".
    """
    raw_frequency = (raw_frequency or "").strip()
    if not raw_frequency or raw_frequency == "-":
        return None
    if raw_frequency in FREQUENCY_QUALIFIER_MIDPOINTS:
        return FREQUENCY_QUALIFIER_MIDPOINTS[raw_frequency]
    try:
        if "/" in raw_frequency:
            numerator, denominator = raw_frequency.split("/", 1)
            denominator_value = float(denominator)
            return min(1.0, max(0.0, float(numerator) / denominator_value)) if denominator_value > 0 else None
        if raw_frequency.endswith("%"):
            return min(1.0, max(0.0, float(raw_frequency[:-1]) / 100.0))
    except ValueError:
        return None
    return None


def parse_frequency_denominator(raw_frequency: str | None) -> int | None:
    """Number of patients behind a fraction such as "3/5" (5); None for qualifiers and percentages."""
    raw_frequency = (raw_frequency or "").strip()
    if "/" in raw_frequency:
        try:
            return int(raw_frequency.split("/", 1)[1])
        except ValueError:
            return None
    return None


BIOCURATION_DATE_PATTERN = re.compile(r"\[(\d{4})-(\d{2})-(\d{2})\]")


def load_hpo_annotation_dates(phenotype_hpoa_path: Path, dated_provenance_prefixes: tuple[str, ...] = ("OMIM:",)) -> dict[tuple[str, str], date]:
    """(disease id, HPO term id) -> earliest biocuration date, for annotations from the dated provenances.

    Rows with the NOT qualifier are skipped. Orphanet rows are excluded by default because their
    biocuration date is the release import date.
    """
    dates: dict[tuple[str, str], date] = {}
    with open(phenotype_hpoa_path, encoding="utf-8") as hpoa_file:
        header: list[str] | None = None
        for raw_line in hpoa_file:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                column = {name: index for index, name in enumerate(header)}
                continue
            disease_id = fields[column["database_id"]]
            if not disease_id.startswith(dated_provenance_prefixes) or fields[column["qualifier"]].strip().upper() == "NOT":
                continue
            matches = BIOCURATION_DATE_PATTERN.findall(fields[column["biocuration"]])
            if not matches:
                continue
            earliest = min(date(int(year), int(month), int(day)) for year, month, day in matches)
            key = (disease_id, fields[column["hpo_id"]])
            if key not in dates or earliest < dates[key]:
                dates[key] = earliest
    return dates


def parse_genes_to_phenotype(genes_to_phenotype_path: Path) -> list[dict[str, str]]:
    """Rows of genes_to_phenotype.txt as dictionaries keyed by the header names."""
    with open(genes_to_phenotype_path, encoding="utf-8") as annotation_file:
        header_line = annotation_file.readline().lstrip("#").strip()
        field_names = header_line.split("\t")
        reader = csv.DictReader(annotation_file, fieldnames=field_names, delimiter="\t")
        return [row for row in reader]


def monogenic_evidence_records(
    annotation_rows: Iterable[dict[str, str]],
    target_symptom_to_hpo_ids: Mapping[str, Iterable[str]],
    parents_by_term: Mapping[str, set[str]],
    genes_in_graph: set[str],
    excluded_hpo_ids_by_symptom: Mapping[str, Iterable[str]] | None = None,
    annotation_dates: Mapping[tuple[str, str], date] | None = None,
) -> list[EvidenceRecord]:
    """One record per (gene symbol, target symptom) for genes present in the physiology graph.

    Rows whose frequency qualifier is Excluded (HP:0040285) are dropped: they assert that the
    feature is absent in that disease. A pair whose rows are all Excluded yields no record.
    With annotation_dates (load_hpo_annotation_dates), each record carries the earliest date of
    the dated (disease, term) annotations behind it, or None when none of them is dated.
    """
    symptoms_by_term = expand_symptom_terms(target_symptom_to_hpo_ids, parents_by_term, excluded_hpo_ids_by_symptom)
    disease_ids_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    frequencies_by_pair: dict[tuple[str, str], list[float]] = defaultdict(list)
    denominators_by_pair: dict[tuple[str, str], int] = defaultdict(int)
    rows_by_pair: dict[tuple[str, str], int] = defaultdict(int)
    dates_by_pair: dict[tuple[str, str], date] = {}
    for row in annotation_rows:
        gene_symbol = row.get("gene_symbol", "")
        if gene_symbol not in genes_in_graph:
            continue
        target_symptoms = symptoms_by_term.get(row.get("hpo_id", ""))
        if not target_symptoms:
            continue
        raw_frequency = (row.get("frequency") or "").strip()
        if raw_frequency == EXCLUDED_FREQUENCY_QUALIFIER:
            continue
        frequency = parse_frequency_qualifier(raw_frequency)
        denominator = parse_frequency_denominator(raw_frequency)
        annotation_date = annotation_dates.get((row.get("disease_id", ""), row.get("hpo_id", ""))) if annotation_dates else None
        for target_symptom in target_symptoms:
            pair = (gene_symbol, target_symptom)
            disease_ids_by_pair[pair].add(row.get("disease_id", ""))
            rows_by_pair[pair] += 1
            if annotation_date is not None and (pair not in dates_by_pair or annotation_date < dates_by_pair[pair]):
                dates_by_pair[pair] = annotation_date
            if frequency is not None:
                frequencies_by_pair[pair].append(frequency)
            if denominator is not None:
                denominators_by_pair[pair] += denominator
    records: list[EvidenceRecord] = []
    for (gene_symbol, target_symptom), disease_ids in disease_ids_by_pair.items():
        omim_entries = {disease_id for disease_id in disease_ids if disease_id.startswith("OMIM:")}
        orpha_entries = {disease_id for disease_id in disease_ids if disease_id.startswith("ORPHA:")}
        frequencies = frequencies_by_pair.get((gene_symbol, target_symptom), [])
        records.append(
            EvidenceRecord(
                perturbation_identifier=gene_symbol,
                symptom_identifier=target_symptom,
                relation="induces",
                evidence_class="monogenic",
                source="HPO genes_to_phenotype; diseases=" + ";".join(sorted(disease_ids)),
                independent_case_series_count=len(disease_ids),  # the version 0.3 proxy; see assign_evidence_grades
                has_omim_clinical_synopsis=bool(omim_entries),
                disease_identifiers=sorted(disease_ids),
                omim_entry_count=len(omim_entries),
                orpha_entry_count=len(orpha_entries),
                annotation_row_count=rows_by_pair[(gene_symbol, target_symptom)],
                max_annotation_frequency=max(frequencies) if frequencies else None,
                annotation_patient_count=denominators_by_pair.get((gene_symbol, target_symptom)) or None,
                evidence_available_date=dates_by_pair.get((gene_symbol, target_symptom)),
            )
        )
    return records
