"""Evidence class E3, PubTator3 half: relations between a target-symptom MeSH descriptor and a graph gene or a SIDER drug,
one row per (PMID, relation), as reports in the evidence_reports.parquet schema (design sections 4.2 and 4.3;
docs/evidence_reports_spec.md section 1). Supersedes the stub load_literature_predications.py.

Literature reports are soft priors only (grades C to E): they never become evaluation positives, a pair that no
paper mentions is unlabelled rather than negative, and every row carries the exact extraction model behind it
(PubTator3 relation type, descriptor level, chemical match method) so the reliability model can weight it by its
limitations rather than by its count.

The bulk file relation2pubtator3.gz (ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/) has one tab-separated row per
(PMID, relation type, entity 1, entity 2) with entities written "Type|Identifier": Disease|MESH:D001007,
Gene|8813, Chemical|MESH:D005473, Variant and Species forms. It carries no publication date; dates come from
E-utilities esummary when the fetch script is asked for them.

Relation types are mapped to the repository's relations as follows and every drop is counted:
  cause, positive_correlate (positive_correlation)      -> induces, report_value 1
  treat, prevent, negative_correlate (negative_correlation) -> relieves, report_value 1
  associate (association)                               -> associated_with, report_value 1, undirected, its own relation
  stimulate, inhibit                                    chemical-gene relations, dropped with a count
  cotreat, compare, interact, drug_interact             not perturbation-symptom relations, dropped with a count
"""
from __future__ import annotations

import csv
import gzip
import json
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import SymptomMeshDescriptor
from mechanistic_pathway_learning.evidence.sider_chemical_matching import ChemicalMatch

PUBTATOR_SOURCE_NAME = "PubTator3"
LITERATURE_EVIDENCE_CLASS = "literature"
ASSOCIATED_WITH_RELATION = "associated_with"

RELATION_TYPE_TO_RELATION: dict[str, str] = {
    "cause": "induces",
    "positive_correlate": "induces",
    "positive_correlation": "induces",
    "treat": "relieves",
    "prevent": "relieves",
    "negative_correlate": "relieves",
    "negative_correlation": "relieves",
    "associate": ASSOCIATED_WITH_RELATION,
    "association": ASSOCIATED_WITH_RELATION,
}
CHEMICAL_GENE_RELATION_TYPES: tuple[str, ...] = ("stimulate", "inhibit")
DISEASE_ENTITY_TYPE = "Disease"
GENE_ENTITY_TYPE = "Gene"
CHEMICAL_ENTITY_TYPE = "Chemical"

# The columns of evidence_reports.parquet (docs/evidence_reports_spec.md section 1), in order; fixed there.
EVIDENCE_REPORT_COLUMNS: tuple[str, ...] = (
    "report_id", "perturbation_id", "perturbation_type", "perturbation_label", "symptom", "relation", "evidence_class", "source",
    "report_value", "source_record_id", "source_record_label", "source_term_id", "source_term_label", "evidence_code",
    "model_description", "frequency", "frequency_denominator", "placebo_flag", "onset", "sex", "references",
    "pubmed_reference_count", "evidence_date", "rubric_log_sample_size", "rubric_frequency_known", "rubric_evidence_code_pcs",
    "rubric_evidence_code_tas", "rubric_evidence_code_iea", "rubric_placebo_controlled", "rubric_curated_synopsis",
    "limitations", "perturbation_nodes",
)
# Literature-specific columns appended after the fixed set; the assembler may ignore them.
LITERATURE_EXTRA_COLUMNS: tuple[str, ...] = (
    "pmid", "relation_type", "entity_role", "perturbation_ncbi_gene_id", "chemical_mesh_id", "chemical_match_method",
    "mesh_descriptor", "mesh_descriptor_level", "publication_year", "direct_evidence",
)
LITERATURE_REPORT_COLUMNS: tuple[str, ...] = EVIDENCE_REPORT_COLUMNS + LITERATURE_EXTRA_COLUMNS


def map_relation_type(relation_type: str) -> tuple[str, int] | None:
    """(relation, report_value) for a PubTator3 relation type, or None when the type is not a perturbation-symptom relation."""
    relation = RELATION_TYPE_TO_RELATION.get(relation_type.strip().lower())
    return None if relation is None else (relation, 1)


def parse_entity(entity_field: str) -> tuple[str, list[str]]:
    """("Disease", ["D001007"]) from "Disease|MESH:D001007"; ("Gene", ["1234", "5678"]) from "Gene|1234;5678" (a mention
    normalized to several genes yields one identifier per gene). The MESH: prefix is stripped from disease and
    chemical identifiers; other identifier schemes are kept verbatim."""
    entity_type, separator, identifier = entity_field.partition("|")
    if not separator:
        return entity_field, []
    identifiers = [part[5:] if part.upper().startswith("MESH:") else part for part in identifier.split(";") if part]
    return entity_type, identifiers


def iterate_bulk_relation_rows(path_or_handle) -> Iterator[tuple[str, str, str, str]]:
    """(pmid, relation_type, entity1, entity2) per row of relation2pubtator3.gz, read as a stream."""
    if isinstance(path_or_handle, (str, Path)):
        with gzip.open(path_or_handle, "rt", encoding="utf-8", newline="") as handle:
            yield from _iterate_rows(handle)
    else:
        yield from _iterate_rows(path_or_handle)


def _iterate_rows(handle) -> Iterator[tuple[str, str, str, str]]:
    for row in csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
        if len(row) >= 4:
            yield row[0], row[1], row[2], row[3]


@dataclass(frozen=True)
class PubtatorRelationRow:
    pmid: str
    relation_type: str
    relation: str
    report_value: int
    entity_role: str  # "disease_first" when entity 1 is the disease, "disease_second" otherwise
    perturbation_type: str  # "gene" or "drug"
    perturbation_id: str  # gene symbol or STITCH flat id
    perturbation_label: str  # gene symbol or SIDER drug name
    perturbation_ncbi_gene_id: int | None
    chemical_mesh_id: str | None
    chemical_match_method: str | None
    symptom: str
    mesh_descriptor: str
    mesh_descriptor_level: str


def disease_rows_from_stream(rows: Iterable[tuple[str, str, str, str]], descriptors: dict[str, SymptomMeshDescriptor]) -> tuple[list[tuple[str, str, str, str]], Counter]:
    """First pass over the bulk stream: keep every row with a target descriptor on either side and count the rest.
    The kept rows (a few hundred thousand) fit in memory, so the partner filters can run after the chemical names
    have been resolved without a second pass over the file."""
    kept: list[tuple[str, str, str, str]] = []
    counts: Counter = Counter()
    for row in rows:
        counts["bulk_rows_total"] += 1
        counts[f"bulk_rows_by_type:{row[1]}"] += 1
        first_type, first_ids = parse_entity(row[2])
        second_type, second_ids = parse_entity(row[3])
        first_is_target = first_type == DISEASE_ENTITY_TYPE and any(identifier in descriptors for identifier in first_ids)
        second_is_target = second_type == DISEASE_ENTITY_TYPE and any(identifier in descriptors for identifier in second_ids)
        if first_is_target or second_is_target:
            kept.append(row)
            counts["bulk_rows_with_target_descriptor"] += 1
    return kept, counts


def filter_relation_rows(
    rows: Iterable[tuple[str, str, str, str]],
    descriptors: dict[str, SymptomMeshDescriptor],
    ncbi_gene_id_to_symbol: dict[int, str],
    chemical_matches_by_mesh_id: dict[str, list[ChemicalMatch]],
) -> tuple[list[PubtatorRelationRow], Counter]:
    """Rows with a target descriptor on one side and a graph gene or a matched SIDER drug on the other, with the
    relation type mapped; every excluded row is counted under the reason it was excluded."""
    kept: list[PubtatorRelationRow] = []
    counts: Counter = Counter()
    for pmid, relation_type, first_field, second_field in rows:
        first_type, first_ids = parse_entity(first_field)
        second_type, second_ids = parse_entity(second_field)
        first_targets = [identifier for identifier in first_ids if first_type == DISEASE_ENTITY_TYPE and identifier in descriptors]
        second_targets = [identifier for identifier in second_ids if second_type == DISEASE_ENTITY_TYPE and identifier in descriptors]
        if not first_targets and not second_targets:
            counts["dropped_no_target_descriptor"] += 1
            continue
        if first_targets and second_targets:
            counts["dropped_both_sides_target_descriptors"] += 1
            continue
        if first_targets:
            entity_role, disease_ids, partner_type, partner_ids = "disease_first", first_targets, second_type, second_ids
        else:
            entity_role, disease_ids, partner_type, partner_ids = "disease_second", second_targets, first_type, first_ids
        relation_type_lower = relation_type.strip().lower()
        if relation_type_lower in CHEMICAL_GENE_RELATION_TYPES:
            counts[f"dropped_chemical_gene_relation_type:{relation_type_lower}"] += 1
            continue
        mapped = map_relation_type(relation_type_lower)
        if mapped is None:
            counts[f"dropped_relation_type_not_mapped:{relation_type_lower}"] += 1
            continue
        relation, report_value = mapped
        if partner_type == GENE_ENTITY_TYPE:
            partners = _gene_partners(partner_ids, ncbi_gene_id_to_symbol, counts)
        elif partner_type == CHEMICAL_ENTITY_TYPE:
            partners = _chemical_partners(partner_ids, chemical_matches_by_mesh_id, counts)
        else:
            counts[f"dropped_partner_type:{partner_type}"] += 1
            continue
        if not partners:
            continue
        if len(partner_ids) > 1:
            counts["rows_with_several_partner_identifiers"] += 1
        for descriptor_id in disease_ids:
            descriptor = descriptors[descriptor_id]
            for perturbation_type, perturbation_id, perturbation_label, ncbi_gene_id, chemical_mesh_id, match_method in partners:
                kept.append(
                    PubtatorRelationRow(
                        pmid=pmid, relation_type=relation_type_lower, relation=relation, report_value=report_value, entity_role=entity_role,
                        perturbation_type=perturbation_type, perturbation_id=perturbation_id, perturbation_label=perturbation_label,
                        perturbation_ncbi_gene_id=ncbi_gene_id, chemical_mesh_id=chemical_mesh_id, chemical_match_method=match_method,
                        symptom=descriptor.target_symptom, mesh_descriptor=descriptor.mesh_descriptor, mesh_descriptor_level=descriptor.level,
                    )
                )
                counts["rows_kept"] += 1
    return kept, counts


def _gene_partners(partner_ids: list[str], ncbi_gene_id_to_symbol: dict[int, str], counts: Counter) -> list[tuple]:
    partners: list[tuple] = []
    for identifier in partner_ids:
        try:
            ncbi_gene_id = int(identifier)
        except ValueError:
            counts["dropped_gene_identifier_not_integer"] += 1
            continue
        symbol = ncbi_gene_id_to_symbol.get(ncbi_gene_id)
        if symbol is None:
            counts["dropped_gene_not_in_graph"] += 1
            continue
        partners.append(("gene", symbol, symbol, ncbi_gene_id, None, None))
    return partners


def _chemical_partners(partner_ids: list[str], chemical_matches_by_mesh_id: dict[str, list[ChemicalMatch]], counts: Counter) -> list[tuple]:
    partners: list[tuple] = []
    for identifier in partner_ids:
        matches = chemical_matches_by_mesh_id.get(identifier)
        if not matches:
            counts["dropped_chemical_not_matched_to_sider"] += 1
            continue
        for match in matches:
            partners.append(("drug", match.stitch_flat_id, match.drug_name, None, match.chemical_mesh_id, match.match_method))
    return partners


def publication_year_from_iso_date(iso_date: str | None) -> int | None:
    if not iso_date:
        return None
    try:
        return int(iso_date[:4])
    except ValueError:
        return None


def literature_report_record(
    *,
    source: str,
    source_record_id: str,
    source_term_id: str,
    source_term_label: str,
    ordinal: int,
    perturbation_type: str,
    perturbation_id: str,
    perturbation_label: str,
    symptom: str,
    relation: str,
    report_value: int,
    evidence_code: str,
    model_description: str,
    references: str,
    pubmed_reference_count: int,
    evidence_date: str | None,
    limitations: str,
    perturbation_nodes: str,
    extras: dict,
) -> dict:
    """One report row in the fixed schema plus the literature extras; the rubric one-hots are 0 because no literature
    row carries an HPO evidence code, a frequency or a placebo arm, and rubric_curated_synopsis is 0 because a paper
    is a named study, not a synopsis. The reliability model therefore treats every literature report with its
    "other" evidence-code factor until the document appraisal supplies rubric features."""
    record = {
        "report_id": f"{source}|{source_record_id}|{source_term_id}|{symptom}|{ordinal}",
        "perturbation_id": perturbation_id,
        "perturbation_type": perturbation_type,
        "perturbation_label": perturbation_label,
        "symptom": symptom,
        "relation": relation,
        "evidence_class": LITERATURE_EVIDENCE_CLASS,
        "source": source,
        "report_value": int(report_value),
        "source_record_id": source_record_id,
        "source_record_label": source_record_id,
        "source_term_id": source_term_id,
        "source_term_label": source_term_label,
        "evidence_code": evidence_code,
        "model_description": model_description,
        "frequency": None,
        "frequency_denominator": None,
        "placebo_flag": False,
        "onset": None,
        "sex": None,
        "references": references,
        "pubmed_reference_count": int(pubmed_reference_count),
        "evidence_date": evidence_date,
        "rubric_log_sample_size": 0.0,
        "rubric_frequency_known": 0.0,
        "rubric_evidence_code_pcs": 0.0,
        "rubric_evidence_code_tas": 0.0,
        "rubric_evidence_code_iea": 0.0,
        "rubric_placebo_controlled": 0.0,
        "rubric_curated_synopsis": 0.0,
        "limitations": limitations,
        "perturbation_nodes": perturbation_nodes,
    }
    for column in LITERATURE_EXTRA_COLUMNS:
        record[column] = extras.get(column)
    return record


def gene_perturbation_nodes_json(gene_symbol: str) -> str:
    """The graph node of a gene with the loss-of-function sign of the monogenic class. A literature relation does not
    say which direction the gene was perturbed, so the sign is a documented convention (model_description says so)
    that the assembler may override; drugs get an empty list and are joined to their ChEMBL targets by STITCH flat
    id by the assembler, as SIDER rows are."""
    return json.dumps([[f"GENE:{gene_symbol}", -1.0, 1.0]])


def pubtator_limitations_text(row: PubtatorRelationRow, dated: bool) -> str:
    clauses = ["machine-extracted relation from PubTator3, no study design read"]
    if row.relation == ASSOCIATED_WITH_RELATION:
        clauses.append("undirected association")
    if row.mesh_descriptor_level == "diagnosis":
        clauses.append("diagnosis-level descriptor, not the symptom itself")
    if row.perturbation_type == "gene":
        clauses.append("gene mention, perturbation direction not stated")
    else:
        clauses.append(f"drug matched to MeSH by {row.chemical_match_method}")
    if not dated:
        clauses.append("undated")
    clauses.append("soft prior only, never an evaluation positive")
    return "; ".join(clauses)


def relation_rows_to_reports(rows: list[PubtatorRelationRow], publication_dates: dict[str, str] | None = None) -> pd.DataFrame:
    """The report table for the PubTator3 rows in input order; ordinal numbers a row within its (pmid, descriptor,
    symptom) group so report_id is unique and stable for a pinned bulk file."""
    publication_dates = publication_dates or {}
    ordinals: Counter = Counter()
    records: list[dict] = []
    for row in rows:
        group_key = (row.pmid, row.mesh_descriptor, row.symptom)
        ordinals[group_key] += 1
        evidence_date = publication_dates.get(row.pmid)
        model = f"human gene mention; PubTator3 {row.relation_type}; sign set to loss-of-function by convention" if row.perturbation_type == "gene" else f"human; drug mention; PubTator3 {row.relation_type}; ChEMBL targets joined by the assembler"
        records.append(
            literature_report_record(
                source=PUBTATOR_SOURCE_NAME,
                source_record_id=f"PMID:{row.pmid}",
                source_term_id=f"MESH:{row.mesh_descriptor}",
                source_term_label=row.mesh_descriptor,
                ordinal=ordinals[group_key],
                perturbation_type=row.perturbation_type,
                perturbation_id=row.perturbation_id,
                perturbation_label=row.perturbation_label,
                symptom=row.symptom,
                relation=row.relation,
                report_value=row.report_value,
                evidence_code=f"pubtator3_{row.relation_type}",
                model_description=model,
                references=f"PMID:{row.pmid}",
                pubmed_reference_count=1,
                evidence_date=evidence_date,
                limitations=pubtator_limitations_text(row, evidence_date is not None),
                perturbation_nodes=gene_perturbation_nodes_json(row.perturbation_id) if row.perturbation_type == "gene" else "[]",
                extras={
                    "pmid": row.pmid,
                    "relation_type": row.relation_type,
                    "entity_role": row.entity_role,
                    "perturbation_ncbi_gene_id": row.perturbation_ncbi_gene_id,
                    "chemical_mesh_id": row.chemical_mesh_id,
                    "chemical_match_method": row.chemical_match_method,
                    "mesh_descriptor": row.mesh_descriptor,
                    "mesh_descriptor_level": row.mesh_descriptor_level,
                    "publication_year": publication_year_from_iso_date(evidence_date),
                    "direct_evidence": None,
                },
            )
        )
    return reports_dataframe(records)


def reports_dataframe(records: list[dict]) -> pd.DataFrame:
    """Stable dtypes for the literature report table, empty-safe."""
    table = pd.DataFrame(records, columns=list(LITERATURE_REPORT_COLUMNS))
    table["report_value"] = table.report_value.astype("int64")
    table["pubmed_reference_count"] = table.pubmed_reference_count.astype("int64")
    table["frequency"] = pd.to_numeric(table.frequency, errors="coerce").astype("float64")
    table["frequency_denominator"] = pd.array(table.frequency_denominator.tolist(), dtype="Int64")
    table["placebo_flag"] = table.placebo_flag.astype(bool)
    table["perturbation_ncbi_gene_id"] = pd.array(table.perturbation_ncbi_gene_id.tolist(), dtype="Int64")
    table["publication_year"] = pd.array(table.publication_year.tolist(), dtype="Int64")
    for column in ("rubric_log_sample_size", "rubric_frequency_known", "rubric_evidence_code_pcs", "rubric_evidence_code_tas", "rubric_evidence_code_iea", "rubric_placebo_controlled", "rubric_curated_synopsis"):
        table[column] = table[column].astype("float64")
    for column in ("onset", "sex", "evidence_date", "chemical_mesh_id", "chemical_match_method", "direct_evidence", "pmid"):
        table[column] = table[column].astype("object").where(table[column].notna(), None)
    return table


def summarize_reports(table: pd.DataFrame) -> dict:
    """Counts for the spec file: per symptom, relation, perturbation type, descriptor level and distinct PMIDs."""
    if table.empty:
        return {"rows": 0}
    return {
        "rows": int(len(table)),
        "distinct_pmids": int(table.pmid.dropna().nunique()),
        "distinct_perturbations": int(table.perturbation_id.nunique()),
        "distinct_genes": int(table.loc[table.perturbation_type == "gene", "perturbation_id"].nunique()),
        "distinct_drugs": int(table.loc[table.perturbation_type == "drug", "perturbation_id"].nunique()),
        "distinct_pairs": int(table[["perturbation_id", "symptom", "relation"]].drop_duplicates().shape[0]),
        "rows_by_symptom": {key: int(value) for key, value in table.symptom.value_counts().sort_index().items()},
        "rows_by_relation": {key: int(value) for key, value in table.relation.value_counts().sort_index().items()},
        "rows_by_relation_type": {key: int(value) for key, value in table.relation_type.value_counts().sort_index().items()} if table.relation_type.notna().any() else {},
        "rows_by_perturbation_type": {key: int(value) for key, value in table.perturbation_type.value_counts().sort_index().items()},
        "rows_by_descriptor_level": {key: int(value) for key, value in table.mesh_descriptor_level.value_counts().sort_index().items()},
        "rows_by_descriptor": {key: int(value) for key, value in table.mesh_descriptor.value_counts().sort_index().items()},
        "rows_by_chemical_match_method": {key: int(value) for key, value in table.chemical_match_method.dropna().value_counts().sort_index().items()},
        "rows_dated": int(table.evidence_date.notna().sum()),
        "rows_by_symptom_and_perturbation_type": {f"{symptom}:{perturbation_type}": int(count) for (symptom, perturbation_type), count in table.groupby(["symptom", "perturbation_type"]).size().items()},
    }


def load_pubtator_relations(path: Path) -> pd.DataFrame:
    """The filtered table written by experiments/fetch_pubtator_relations.py."""
    return pd.read_parquet(path)
