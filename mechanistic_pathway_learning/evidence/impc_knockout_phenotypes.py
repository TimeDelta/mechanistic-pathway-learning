"""IMPC knockout phenotypes as cross-species evidence reports, positive and negative, per human gene and symptom.

What the source is. The International Mouse Phenotyping Consortium phenotypes single-gene knockouts on a standard
battery (open field, light-dark box, fear conditioning, sleep-wake, indirect calorimetry among them) and publishes
two things: the significant genotype effects as Mammalian Phenotype (MP) terms (genotype-phenotype-assertions), and
the full statistical results, significant or not (statistical-results-ALL). The second is what makes this source
different from every curated database here: a knockout tested on the open field with no locomotor effect is a
measured absence at that test's power, not silence.

What a report is. A mouse gene is mapped to its human gene through HGNC's MGI identifiers, one to one, and to the
graph's gene node; the human gene is the loss-of-function perturbation, as in the monogenic class. A significant MP
term that configs/impc_symptom_crosswalk.csv maps to a symptom is a positive (report_value 1) for that symptom and
relation; a procedure the crosswalk lists as a negative procedure for a symptom, on which the gene was tested with
no significant call among the symptom's terms, is one negative (report_value 0) per (gene, symptom, procedure code,
zygosity), carrying the parameters tested, the smallest genotype p-value and the mutant counts. The crosswalk's
validity column says how far the mouse readout stands for the human symptom; "weak" proxies are left out unless
asked for, and "endophenotype_only" and "none" are never read. Every row is evidence_class "animal", the design's
grade D: a soft prior, never a label and never an evaluation positive (docs/experiment_design.md section 4.2).

Licence: "IMPC data is freely available to download and use", Creative Commons Attribution 4.0
(https://www.mousephenotype.org/about-impc/accessing-the-data/, read 15 October 2026). Release 24.0 (the files of
16 March 2026), pinned in scripts/fetch_raw_sources.sh.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator

from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

SOURCE_NAME = "IMPC-knockout"
EVIDENCE_CLASS = "animal"
EVIDENCE_CODE_SIGNIFICANT = "impc_significant"
EVIDENCE_CODE_NOT_SIGNIFICANT = "impc_not_significant"
IMPC_RELEASE = "24.0"
IMPC_RELEASE_DATE = "2026-03-16"
DEFAULT_CROSSWALK_PATH = Path("configs/impc_symptom_crosswalk.csv")
VALIDITY_READ_BY_DEFAULT = ("strong", "moderate")
VALIDITY_READ_WITH_WEAK = ("strong", "moderate", "weak")
VALIDITY_RUBRIC = {"strong": 1.0, "moderate": 0.6, "weak": 0.3}


@dataclass(frozen=True)
class CrosswalkRow:
    symptom: str
    relation: str
    mp_term_ids: tuple[str, ...]
    negative_procedure_codes: tuple[str, ...]
    validity: str
    note: str


@dataclass
class GeneProcedureTests:
    """What one knockout's tests on one procedure code and zygosity showed."""
    parameters_tested: set[str] = field(default_factory=set)
    significant_terms: set[str] = field(default_factory=set)
    smallest_p_value: float | None = None
    mutant_count: int = 0
    control_count: int = 0
    procedure_names: set[str] = field(default_factory=set)
    centers: set[str] = field(default_factory=set)


def read_crosswalk(path: Path = DEFAULT_CROSSWALK_PATH) -> list[CrosswalkRow]:
    rows = []
    with open(path, encoding="utf-8") as crosswalk_file:
        for row in csv.DictReader(crosswalk_file):
            rows.append(CrosswalkRow(row["target_symptom"], row["relation"],
                                     tuple(term.strip() for term in row["mp_term_ids"].split(";") if term.strip()),
                                     tuple(code.strip() for code in row["negative_procedure_codes"].split(";") if code.strip()),
                                     row["validity"].strip(), row.get("note", "")))
    return rows


def procedure_code(procedure_stable_id: str) -> str | None:
    """The procedure family of an IMPC procedure id: IMPC_OFD_001, JAX_OFD_001 and HRWLLA_OFD_001 are all "OFD". A legacy
    numeric code (ESLIM_008_001, M-G-P_007_001) names no family and is left unread."""
    parts = str(procedure_stable_id).split("_")
    if len(parts) < 3:
        return None
    code = parts[-2]
    return code if code.isalpha() else None


def mouse_to_human_symbols(hgnc_path: Path) -> dict[str, str]:
    """MGI accession -> HGNC symbol, one to one; an MGI id listed on two HGNC rows, or a row listing two MGI ids, is dropped."""
    candidates: dict[str, set[str]] = defaultdict(set)
    with open(hgnc_path, encoding="utf-8") as hgnc_file:
        reader = csv.DictReader(hgnc_file, delimiter="\t")
        for row in reader:
            mgi_ids = [identifier.strip() for identifier in str(row.get("mgd_id") or "").split("|") if identifier.strip()]
            if len(mgi_ids) != 1:
                continue
            candidates[mgi_ids[0]].add(row["symbol"])
    return {mgi: next(iter(symbols)) for mgi, symbols in candidates.items() if len(symbols) == 1}


def read_statistical_results(path: Path, codes_of_interest: set[str]) -> Iterator[dict]:
    """Rows of statistical-results-ALL.csv.gz whose procedure family is of interest, streamed."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as results_file:
        for row in csv.DictReader(results_file):
            code = procedure_code(row.get("procedure_stable_id", ""))
            if code in codes_of_interest:
                yield row


def to_float(value) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(number) else number


def to_int(value) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def collect_tests(rows: Iterable[dict]) -> dict[tuple[str, str, str], GeneProcedureTests]:
    """(MGI accession, procedure code, zygosity) -> what the tests showed, over every successful row."""
    tests: dict[tuple[str, str, str], GeneProcedureTests] = defaultdict(GeneProcedureTests)
    for row in rows:
        if str(row.get("status", "")).strip() != "Successful":
            continue
        code = procedure_code(row.get("procedure_stable_id", ""))
        key = (row.get("marker_accession_id", ""), code or "", row.get("zygosity", ""))
        record = tests[key]
        record.parameters_tested.add(row.get("parameter_stable_id", ""))
        record.procedure_names.add(row.get("procedure_name", ""))
        record.centers.add(row.get("phenotyping_center", ""))
        p_value = to_float(row.get("p_value")) or to_float(row.get("genotype_effect_p_value"))
        if p_value is not None and (record.smallest_p_value is None or p_value < record.smallest_p_value):
            record.smallest_p_value = p_value
        mutants = to_int(row.get("male_mutant_count")) + to_int(row.get("female_mutant_count"))
        controls = to_int(row.get("male_control_count")) + to_int(row.get("female_control_count"))
        record.mutant_count = max(record.mutant_count, mutants)
        record.control_count = max(record.control_count, controls)
        if str(row.get("significant", "")).strip().lower() == "true" and row.get("mp_term_id"):
            record.significant_terms.add(row["mp_term_id"])
    return tests


def collect_assertions(path: Path) -> dict[str, dict[str, dict]]:
    """MGI accession -> MP term id -> the strongest assertion (smallest p) with its procedure, parameter, zygosity, effect."""
    assertions: dict[str, dict[str, dict]] = defaultdict(dict)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as assertions_file:
        for row in csv.DictReader(assertions_file):
            term = row.get("mp_term_id", "")
            if not term:
                continue
            p_value = to_float(row.get("p_value"))
            existing = assertions[row["marker_accession_id"]].get(term)
            if existing is None or (p_value is not None and (existing["p_value"] is None or p_value < existing["p_value"])):
                assertions[row["marker_accession_id"]][term] = {
                    "p_value": p_value, "procedure_stable_id": row.get("procedure_stable_id", ""), "procedure_name": row.get("procedure_name", ""),
                    "parameter_name": row.get("parameter_name", ""), "zygosity": row.get("zygosity", ""), "effect_size": to_float(row.get("effect_size")),
                    "percentage_change": row.get("percentage_change", ""), "marker_symbol": row.get("marker_symbol", ""),
                    "allele_symbol": row.get("allele_symbol", ""), "center": row.get("phenotyping_center", ""), "mp_term_name": row.get("mp_term_name", "")}
    return assertions


def base_row(human_symbol: str, node_id: str, symptom: str, relation: str, validity: str) -> dict:
    row = {column: "" for column in REPORT_COLUMNS}
    row.update({
        "perturbation_id": human_symbol, "perturbation_type": "gene", "perturbation_label": human_symbol, "symptom": symptom, "relation": relation,
        "evidence_class": EVIDENCE_CLASS, "source": SOURCE_NAME, "frequency": None, "frequency_denominator": None, "placebo_flag": False, "onset": None,
        "sex": None, "references": f"IMPC release {IMPC_RELEASE}", "pubmed_reference_count": 0, "evidence_date": IMPC_RELEASE_DATE, "evidence_date_source": "release",
        "rubric_frequency_known": 0.0, "rubric_evidence_code_pcs": 0.0, "rubric_evidence_code_tas": 0.0, "rubric_evidence_code_iea": 0.0,
        "rubric_placebo_controlled": 0.0, "rubric_curated_synopsis": 0.0, "association_type": "", "rubric_causal_association": 0.0,
        "perturbation_nodes": json.dumps([[node_id, -1.0, 1.0]]),
        # extra rubric features the reliability fit can read: a cross-species row, and how far the readout stands for the symptom
        "rubric_cross_species": 1.0, "rubric_proxy_validity": VALIDITY_RUBRIC.get(validity, 0.0),
    })
    return row


def build_reports(assertions: dict[str, dict[str, dict]], tests: dict[tuple[str, str, str], GeneProcedureTests], mouse_to_human: dict[str, str],
                  node_of_symbol: dict[str, str], crosswalk: list[CrosswalkRow], validities: tuple[str, ...] = VALIDITY_READ_BY_DEFAULT) -> tuple[list[dict], dict]:
    """Report rows (REPORT_COLUMNS plus two rubric columns) and counts.

    Positives: for every mapped MP term a knockout carries. Negatives: for every (gene, symptom, procedure code,
    zygosity) tested on a procedure the crosswalk lists for the symptom with no significant call among the symptom's
    terms (either relation), one row. A gene whose mouse symbol has no one-to-one human gene, or whose human gene
    has no graph node, produces nothing and is counted.
    """
    read_rows = [row for row in crosswalk if row.validity in validities]
    term_to_rows: dict[str, list[CrosswalkRow]] = defaultdict(list)
    for row in read_rows:
        for term in row.mp_term_ids:
            term_to_rows[term].append(row)
    negative_rows_by_code: dict[str, list[CrosswalkRow]] = defaultdict(list)
    for row in read_rows:
        for code in row.negative_procedure_codes:
            negative_rows_by_code[code].append(row)
    terms_of_symptom: dict[str, set[str]] = defaultdict(set)
    for row in crosswalk:  # every relation and validity: a weak-validity hit still blocks a negative for its symptom
        terms_of_symptom[row.symptom].update(row.mp_term_ids)

    rows: list[dict] = []
    counts = {"mouse_genes_with_assertions": len(assertions), "mouse_genes_tested_on_read_procedures": len({key[0] for key in tests}),
              "no_one_to_one_human_gene": 0, "human_gene_not_in_graph": 0, "positives": 0, "negatives": 0, "negatives_blocked_by_a_hit": 0}
    skipped_genes: set[str] = set()

    def human_node(mgi: str) -> tuple[str, str] | None:
        human = mouse_to_human.get(mgi)
        if human is None:
            if mgi not in skipped_genes:
                counts["no_one_to_one_human_gene"] += 1
                skipped_genes.add(mgi)
            return None
        node = node_of_symbol.get(human)
        if node is None:
            if mgi not in skipped_genes:
                counts["human_gene_not_in_graph"] += 1
                skipped_genes.add(mgi)
            return None
        return human, node

    for mgi, terms in sorted(assertions.items()):
        resolved = human_node(mgi)
        if resolved is None:
            continue
        human, node = resolved
        for term, assertion in sorted(terms.items()):
            for crosswalk_row in term_to_rows.get(term, []):
                row = base_row(human, node, crosswalk_row.symptom, crosswalk_row.relation, crosswalk_row.validity)
                row.update({
                    "report_id": f"{SOURCE_NAME}|{human}|{term}|{crosswalk_row.symptom}|{crosswalk_row.relation}|{assertion['zygosity']}",
                    "report_value": 1, "source_record_id": mgi, "source_record_label": f"{assertion['marker_symbol']} {assertion['allele_symbol']}".strip(),
                    "source_term_id": term, "source_term_label": assertion["mp_term_name"], "evidence_code": EVIDENCE_CODE_SIGNIFICANT,
                    "model_description": (f"mouse knockout ({assertion['zygosity']}, {assertion['center']}); {assertion['procedure_name']}, {assertion['parameter_name']}; "
                                          f"p {assertion['p_value']}, effect size {assertion['effect_size']}, change {assertion['percentage_change']}; human orthologue {human} by HGNC"),
                    "rubric_log_sample_size": 0.0,
                    "limitations": (f"cross-species: a mouse {assertion['mp_term_name']} read as {crosswalk_row.symptom} with {crosswalk_row.validity} validity ({crosswalk_row.note}); "
                                    f"a constitutive knockout on one background strain; the effect is the strongest colony's"),
                })
                rows.append(row)
                counts["positives"] += 1

    for (mgi, code, zygosity), record in sorted(tests.items()):
        if code not in negative_rows_by_code:
            continue
        resolved = human_node(mgi)
        if resolved is None:
            continue
        human, node = resolved
        for crosswalk_row in negative_rows_by_code[code]:
            if record.significant_terms & terms_of_symptom[crosswalk_row.symptom]:
                counts["negatives_blocked_by_a_hit"] += 1
                continue
            row = base_row(human, node, crosswalk_row.symptom, crosswalk_row.relation, crosswalk_row.validity)
            row.update({
                "report_id": f"{SOURCE_NAME}|{human}|{code}|{crosswalk_row.symptom}|absent|{zygosity}",
                "report_value": 0, "source_record_id": mgi, "source_record_label": ", ".join(sorted(record.procedure_names)),
                "source_term_id": "", "source_term_label": f"no significant call on {', '.join(sorted(record.procedure_names))}",
                "evidence_code": EVIDENCE_CODE_NOT_SIGNIFICANT,
                "model_description": (f"mouse knockout ({zygosity}; {', '.join(sorted(record.centers))}); {len(record.parameters_tested)} parameters of "
                                      f"{', '.join(sorted(record.procedure_names))} tested, no significant genotype effect on the symptom's terms; smallest p "
                                      f"{record.smallest_p_value}; mutants {record.mutant_count}, controls {record.control_count}; human orthologue {human} by HGNC"),
                "rubric_log_sample_size": math.log1p(record.mutant_count),
                "limitations": (f"cross-species: absence of a mouse readout taken as absence of {crosswalk_row.symptom} with {crosswalk_row.validity} validity "
                                f"({crosswalk_row.note}); an absence at the power of {record.mutant_count} mutants, not a demonstration of no effect"),
            })
            rows.append(row)
            counts["negatives"] += 1
    return rows, counts


def write_summary(path: Path, counts: dict, crosswalk: list[CrosswalkRow], validities: tuple[str, ...]) -> None:
    path.write_text(json.dumps({"source": SOURCE_NAME, "release": IMPC_RELEASE, "validities_read": list(validities), "counts": counts,
                                "crosswalk": [{"symptom": row.symptom, "relation": row.relation, "mp_term_ids": list(row.mp_term_ids),
                                               "negative_procedure_codes": list(row.negative_procedure_codes), "validity": row.validity} for row in crosswalk]},
                               indent=1) + "\n", encoding="utf-8")
