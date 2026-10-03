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

References. An observation does not account for the number of publications behind it. The
reference column of phenotype.hpoa names the source of each (disease, term) annotation: a PubMed
identifier, the OMIM or Orphanet entry itself, occasionally a book or a URL. Most rows cite exactly
one source, and an OMIM or Orphanet citation stands for the curated synopsis rather than for any
count of studies, so a reference count measures curation style more than independent evidence.
The distinct references behind a (gene, symptom) pair are nevertheless recorded as descriptive
columns (distinct_reference_count, distinct_pubmed_reference_count) so multiplicity weighting can
be tested as an ablation; neither count enters the loss weight.

Reports. phenotype.hpoa is the per-annotation file: one row per (disease, term, reference, evidence code),
with the qualifier (NOT asserts absence), frequency, onset, sex and biocuration dates of that row.
genes_to_phenotype.txt collapses the rows of one (disease, term) key to a single row, keeps one of the
frequencies and does not mark the NOT qualifier, so it cannot tell presence from asserted absence.
monogenic_evidence_reports joins each genes_to_phenotype row back to its phenotype.hpoa rows on
(disease_id, hpo_id) and emits one EvidenceReport per phenotype.hpoa row per target symptom, with
report_value 0 for NOT-qualified or Excluded rows and for rows whose frequency is a patient fraction
with numerator 0 (0/35 is an observed absence in 35 patients, not a presence report; the fraction and
its denominator are kept so the sample-size factor applies to the absence claim) (evidence_reports.py).
The aggregated EvidenceRecord of monogenic_evidence_records remains for the Phase 1 counts and the older
tests, and there a 0/N row is still a frequency-0 row, which the grade weighting of the assembler keeps.

Output: one EvidenceRecord per (gene, target symptom) with evidence_class "monogenic", or one
EvidenceReport per annotation row behind such a pair.
"""
from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord
from mechanistic_pathway_learning.evidence.evidence_reports import UNJOINED_EVIDENCE_CODE, EvidenceReport, limitations_text

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
    inside_term_stanza = False
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line.startswith("["):  # any stanza header ([Term], [Typedef], [Instance]) ends the previous term
                inside_term_stanza = line == "[Term]"
                current_term = None
            elif inside_term_stanza and line.startswith("id: HP:"):
                current_term = line[len("id: "):]
                parents_by_term.setdefault(current_term, set())
            elif inside_term_stanza and line == "is_obsolete: true" and current_term is not None:
                parents_by_term.pop(current_term, None)  # obsolete terms are not part of the hierarchy
                current_term = None
            elif line.startswith("is_a:") and current_term is not None:
                parent_id = line[len("is_a:"):].strip().split(" ")[0]
                parents_by_term[current_term].add(parent_id)
    return dict(parents_by_term)


def load_hpo_term_names_from_obo(obo_path: Path) -> dict[str, str]:
    """Return term id -> name from an OBO file."""
    names_by_term: dict[str, str] = {}
    current_term: str | None = None
    inside_term_stanza = False
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line.startswith("["):
                inside_term_stanza = line == "[Term]"
                current_term = None
            elif inside_term_stanza and line.startswith("id: HP:"):
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


def is_zero_patient_fraction(raw_frequency: str | None) -> bool:
    """True for a patient fraction with numerator 0 and a positive denominator ("0/35"): an observed absence, not a presence report."""
    raw_frequency = (raw_frequency or "").strip()
    if "/" not in raw_frequency:
        return False
    try:
        numerator, denominator = raw_frequency.split("/", 1)
        return float(numerator) == 0.0 and float(denominator) > 0
    except ValueError:
        return False


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


def load_reference_publication_dates(publication_dates_path: Path) -> dict[str, date]:
    """PMID -> publication date from docs/hpo_reference_publication_dates.json (PubMed metadata; month or day may be missing)."""
    import json

    payload = json.loads(Path(publication_dates_path).read_text())
    entries = payload.get("publication_dates_by_pmid", payload)
    dates: dict[str, date] = {}
    for pmid, entry in entries.items():
        year = entry.get("year")
        if not year:
            continue
        dates[str(pmid)] = date(int(year), int(entry.get("month") or 1), int(entry.get("day") or 1))
    return dates


@dataclass
class HpoaAnnotationRow:
    """One row of phenotype.hpoa; ordinal is its 1-based position within its (disease_id, hpo_id) group in file order."""

    disease_id: str
    disease_name: str
    qualifier: str  # "" or "NOT"
    hpo_id: str
    references: list[str]  # prefixes kept: PMID:, OMIM:, ORPHA:, ISBN-13:, http
    evidence: str  # PCS, TAS, IEA or another code verbatim
    onset: str | None
    frequency: str  # raw frequency column ("" when empty)
    sex: str | None
    modifier: str
    aspect: str
    biocuration_dates: list[date]
    ordinal: int

    @property
    def is_negated(self) -> bool:
        """True when the row asserts absence: the NOT qualifier or the Excluded frequency qualifier."""
        return self.qualifier.strip().upper() == "NOT" or self.frequency.strip() == EXCLUDED_FREQUENCY_QUALIFIER


def load_hpo_annotation_rows(phenotype_hpoa_path: Path) -> dict[tuple[str, str], list[HpoaAnnotationRow]]:
    """(disease id, HPO term id) -> its phenotype.hpoa rows in file order, every provenance and qualifier kept."""
    rows_by_key: dict[tuple[str, str], list[HpoaAnnotationRow]] = defaultdict(list)
    with open(phenotype_hpoa_path, encoding="utf-8") as hpoa_file:
        header: list[str] | None = None
        for raw_line in hpoa_file:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            values = {name: fields[index].strip() if index < len(fields) else "" for index, name in enumerate(header)}
            key = (values.get("database_id", ""), values.get("hpo_id", ""))
            matches = BIOCURATION_DATE_PATTERN.findall(values.get("biocuration", ""))
            rows_by_key[key].append(HpoaAnnotationRow(
                disease_id=key[0], disease_name=values.get("disease_name", ""), qualifier=values.get("qualifier", ""), hpo_id=key[1],
                references=[reference.strip() for reference in values.get("reference", "").split(";") if reference.strip()],
                evidence=values.get("evidence", ""), onset=values.get("onset") or None, frequency=values.get("frequency", ""), sex=values.get("sex") or None,
                modifier=values.get("modifier", ""), aspect=values.get("aspect", ""),
                biocuration_dates=[date(int(year), int(month), int(day)) for year, month, day in matches],
                ordinal=len(rows_by_key[key]) + 1,
            ))
    return dict(rows_by_key)


def hpoa_row_availability_date(
    row: HpoaAnnotationRow,
    dated_provenance_prefixes: tuple[str, ...] = ("OMIM:",),
    publication_dates_by_pmid: Mapping[str, date] | None = None,
) -> date | None:
    """Earliest of the row's biocuration dates and the publication dates of its PubMed references, for dated provenances.

    Publication precedes curation, so with the lookup present the date is the publication date. Orphanet rows
    get None by default because their biocuration date is the release import date.
    """
    if not row.disease_id.startswith(dated_provenance_prefixes):
        return None
    candidates = list(row.biocuration_dates)
    if publication_dates_by_pmid:
        for reference in row.references:
            if reference.startswith("PMID:") and reference[5:] in publication_dates_by_pmid:
                candidates.append(publication_dates_by_pmid[reference[5:]])
    return min(candidates) if candidates else None


def load_hpo_annotation_dates(
    phenotype_hpoa_path: Path,
    dated_provenance_prefixes: tuple[str, ...] = ("OMIM:",),
    publication_dates_by_pmid: Mapping[str, date] | None = None,
) -> dict[tuple[str, str], date]:
    """(disease id, HPO term id) -> earliest availability date, for annotations from the dated provenances.

    The availability date of one annotation row is hpoa_row_availability_date; the key's date is the earliest
    over its rows. Rows with the NOT qualifier are skipped. Orphanet rows are excluded by default because
    their biocuration date is the release import date.
    """
    dates: dict[tuple[str, str], date] = {}
    for key, rows in load_hpo_annotation_rows(phenotype_hpoa_path).items():
        for row in rows:
            if row.qualifier.strip().upper() == "NOT":
                continue
            availability_date = hpoa_row_availability_date(row, dated_provenance_prefixes, publication_dates_by_pmid)
            if availability_date is not None and (key not in dates or availability_date < dates[key]):
                dates[key] = availability_date
    return dates


def load_hpo_annotation_references(phenotype_hpoa_path: Path) -> dict[tuple[str, str], set[str]]:
    """(disease id, HPO term id) -> the distinct reference strings cited for that annotation, all provenances.

    Rows with the NOT qualifier are skipped. Reference strings keep their prefix (PMID:, OMIM:, ORPHA:,
    ISBN-13:, http...) so callers can count the PubMed subset separately.
    """
    references: dict[tuple[str, str], set[str]] = {}
    for key, rows in load_hpo_annotation_rows(phenotype_hpoa_path).items():
        cited = {reference for row in rows if row.qualifier.strip().upper() != "NOT" for reference in row.references}
        if cited:
            references[key] = cited
    return references


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
    annotation_references: Mapping[tuple[str, str], set[str]] | None = None,
) -> list[EvidenceRecord]:
    """One record per (gene symbol, target symptom) for genes present in the physiology graph.

    Rows whose frequency qualifier is Excluded (HP:0040285) are dropped: they assert that the
    feature is absent in that disease. A pair whose rows are all Excluded yields no record.
    With annotation_dates (load_hpo_annotation_dates), each record carries the earliest date of
    the dated (disease, term) annotations behind it, or None when none of them is dated. With
    annotation_references (load_hpo_annotation_references), each record counts the distinct references
    and the distinct PubMed references cited across its (disease, term) annotations.
    """
    symptoms_by_term = expand_symptom_terms(target_symptom_to_hpo_ids, parents_by_term, excluded_hpo_ids_by_symptom)
    disease_ids_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    frequencies_by_pair: dict[tuple[str, str], list[float]] = defaultdict(list)
    denominators_by_pair: dict[tuple[str, str], int] = defaultdict(int)
    rows_by_pair: dict[tuple[str, str], int] = defaultdict(int)
    dates_by_pair: dict[tuple[str, str], date] = {}
    references_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
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
        cited_references = annotation_references.get((row.get("disease_id", ""), row.get("hpo_id", ""))) if annotation_references else None
        for target_symptom in target_symptoms:
            pair = (gene_symbol, target_symptom)
            disease_ids_by_pair[pair].add(row.get("disease_id", ""))
            rows_by_pair[pair] += 1
            if cited_references:
                references_by_pair[pair] |= cited_references
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
        cited = references_by_pair.get((gene_symbol, target_symptom), set())
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
                distinct_reference_count=len(cited),
                distinct_pubmed_reference_count=sum(1 for reference in cited if reference.startswith("PMID:")),
            )
        )
    return records


HPO_PROVENANCE_SOURCES: dict[str, str] = {"OMIM:": "HPO-OMIM", "ORPHA:": "HPO-Orphanet"}


def hpo_source_for_disease(disease_id: str) -> str | None:
    """Report source name for a disease identifier prefix, None for a prefix the report table does not know."""
    for prefix, source in HPO_PROVENANCE_SOURCES.items():
        if disease_id.startswith(prefix):
            return source
    return None


def monogenic_evidence_reports(
    annotation_rows: Iterable[dict[str, str]],
    target_symptom_to_hpo_ids: Mapping[str, Iterable[str]],
    parents_by_term: Mapping[str, set[str]],
    genes_in_graph: set[str],
    excluded_hpo_ids_by_symptom: Mapping[str, Iterable[str]] | None = None,
    hpoa_rows_by_key: Mapping[tuple[str, str], list[HpoaAnnotationRow]] | None = None,
    publication_dates_by_pmid: Mapping[str, date] | None = None,
    dropped_unknown_provenance: list[dict[str, str]] | None = None,
    association_types: Mapping[tuple[str, str], tuple[str, bool | None]] | None = None,
) -> list[EvidenceReport]:
    """One EvidenceReport per phenotype.hpoa row behind each (gene in the graph, target symptom), in file order.

    Duplicate (gene_symbol, disease_id, hpo_id) triples of genes_to_phenotype.txt are read once. Each is joined
    to its phenotype.hpoa rows on (disease_id, hpo_id); frequency, qualifier, evidence code, references, onset,
    sex, disease name and date come from the phenotype.hpoa row. When hpoa_rows_by_key is None or has no row
    for the key, the genes_to_phenotype row itself becomes one report with evidence_code "unjoined". A term in
    the expansion of two target symptoms yields one report per symptom. Rows whose disease prefix is neither
    OMIM: nor ORPHA: are dropped and, when a list is given, appended to dropped_unknown_provenance.
    perturbation_nodes is left empty for the assembler to fill.
    """
    symptoms_by_term = expand_symptom_terms(target_symptom_to_hpo_ids, parents_by_term, excluded_hpo_ids_by_symptom)
    seen_triples: set[tuple[str, str, str]] = set()
    reports: list[EvidenceReport] = []
    for row in annotation_rows:
        gene_symbol = row.get("gene_symbol", "")
        hpo_id = row.get("hpo_id", "")
        disease_id = row.get("disease_id", "")
        if gene_symbol not in genes_in_graph:
            continue
        target_symptoms = symptoms_by_term.get(hpo_id)
        if not target_symptoms:
            continue
        triple = (gene_symbol, disease_id, hpo_id)
        if triple in seen_triples:
            continue
        seen_triples.add(triple)
        source = hpo_source_for_disease(disease_id)
        if source is None:
            if dropped_unknown_provenance is not None:
                dropped_unknown_provenance.append(dict(row))
            continue
        hpo_name = row.get("hpo_name", "")
        joined_rows = (hpoa_rows_by_key or {}).get((disease_id, hpo_id)) or []
        if not joined_rows:
            raw_frequency = (row.get("frequency") or "").strip()
            joined_rows = [HpoaAnnotationRow(disease_id, disease_id, "", hpo_id, [], UNJOINED_EVIDENCE_CODE, None, "" if raw_frequency == "-" else raw_frequency, None, "", "", [], 0)]
        for hpoa_row in joined_rows:
            for target_symptom in target_symptoms:
                reports.append(evidence_report_from_hpoa_row(gene_symbol, target_symptom, source, hpo_name, hpoa_row, publication_dates_by_pmid, association_types=association_types))
    return reports


def evidence_report_from_hpoa_row(
    gene_symbol: str,
    target_symptom: str,
    source: str,
    hpo_name: str,
    hpoa_row: HpoaAnnotationRow,
    publication_dates_by_pmid: Mapping[str, date] | None = None,
    association_types: Mapping[tuple[str, str], tuple[str, bool | None]] | None = None,
) -> EvidenceReport:
    """Column derivation of one monogenic report from one phenotype.hpoa row (or the unjoined stand-in row).

    report_value is 0 when the row asserts absence: the NOT qualifier, the Excluded frequency qualifier, or a
    patient fraction with numerator 0 (module docstring, "Reports").
    """
    unjoined = hpoa_row.evidence == UNJOINED_EVIDENCE_CODE
    frequency = 0.0 if hpoa_row.frequency.strip() == EXCLUDED_FREQUENCY_QUALIFIER else parse_frequency_qualifier(hpoa_row.frequency)
    denominator = parse_frequency_denominator(hpoa_row.frequency)
    references = ";".join(hpoa_row.references)
    availability_date = None if unjoined else hpoa_row_availability_date(hpoa_row, publication_dates_by_pmid=publication_dates_by_pmid)
    entry_reference = hpoa_row.disease_id
    association_label, association_causal = (association_types or {}).get((gene_symbol, hpoa_row.disease_id), ("unknown", None))
    report = EvidenceReport(
        report_id=f"{source}|{gene_symbol}|{hpoa_row.disease_id}|{hpoa_row.hpo_id}|{target_symptom}|{hpoa_row.ordinal}",
        perturbation_id=gene_symbol, perturbation_type="gene", perturbation_label=gene_symbol,
        symptom=target_symptom, relation="induces", evidence_class="monogenic", source=source,
        report_value=0 if hpoa_row.is_negated or is_zero_patient_fraction(hpoa_row.frequency) else 1,
        source_record_id=hpoa_row.disease_id, source_record_label=hpoa_row.disease_name,
        source_term_id=hpoa_row.hpo_id, source_term_label=hpo_name,
        evidence_code=hpoa_row.evidence,
        model_description=f"human loss-of-function; {hpoa_row.disease_id} {hpoa_row.disease_name}",
        frequency=frequency, frequency_denominator=denominator, placebo_flag=False,
        onset=hpoa_row.onset, sex=hpoa_row.sex, references=references,
        pubmed_reference_count=sum(1 for reference in hpoa_row.references if reference.startswith("PMID:")),
        evidence_date=availability_date.isoformat() if availability_date else None,
        rubric_log_sample_size=float(math.log1p(denominator)) if denominator is not None else 0.0,
        rubric_frequency_known=1.0 if frequency is not None else 0.0,
        rubric_evidence_code_pcs=1.0 if hpoa_row.evidence == "PCS" else 0.0,
        rubric_evidence_code_tas=1.0 if hpoa_row.evidence == "TAS" else 0.0,
        rubric_evidence_code_iea=1.0 if hpoa_row.evidence == "IEA" else 0.0,
        rubric_placebo_controlled=0.0,
        rubric_curated_synopsis=1.0 if entry_reference in hpoa_row.references else 0.0,
        association_type=association_label, rubric_causal_association=1.0 if association_causal else 0.0,
    )
    report.limitations = limitations_text(report)
    return report
