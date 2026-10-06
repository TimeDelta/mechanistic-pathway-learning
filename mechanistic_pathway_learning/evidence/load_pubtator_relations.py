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
  cause                                        -> induces, report_value 1 (on target rows the partner is a chemical or a variant)
  treat, prevent                               -> relieves, report_value 1
  stimulate, inhibit on a Disease-Gene row     -> induces, report_value 1, perturbation sign from the type: the bulk file
                                                  renders a gene's positive or negative correlation with a disease as
                                                  stimulate (gene activity rises with the disease, sign +1.0) or inhibit
                                                  (gene activity falls with the disease, sign -1.0, the loss-of-function
                                                  direction of the monogenic class); on any other row they are chemical-gene
                                                  relations and are dropped with a count
  positive_correlate, negative_correlate       -> induces, relieves; mapped for completeness, no row of the pinned bulk file
                                                  carries them against a Disease entity
  associate (association)                      -> associated_with, report_value 1, undirected, its own relation
  cotreat, compare, interact, drug_interact    not perturbation-symptom relations, dropped with a count
A directional relation (induces or relieves) on a descriptor of crosswalk level mixed_polarity_diagnosis (Bipolar
Disorder D001714, whose literature covers both mood poles) is demoted to associated_with with a count and a
limitation clause, because "treats bipolar disorder" says nothing about the direction of the mania symptom.

Source names are "PubTator3-<relation type>" so that the assembler's predication parser
(assemble_evidence_table.predication_type_for_literature_sources) and the reliability fit see one sensor per
extraction type.
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
MIXED_POLARITY_DESCRIPTOR_LEVEL = "mixed_polarity_diagnosis"  # crosswalk level whose directional relations are demoted
LOSS_OF_FUNCTION_SIGN = -1.0  # the perturbation sign of the monogenic class, the default for a gene mention without direction

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
# On a Disease-Gene row these types carry the sign of the gene's correlation with the disease (see the module docstring).
DISEASE_GENE_RELATION_TYPE_SIGNS: dict[str, float] = {"stimulate": 1.0, "inhibit": -1.0}
DISEASE_ENTITY_TYPE = "Disease"
GENE_ENTITY_TYPE = "Gene"
CHEMICAL_ENTITY_TYPE = "Chemical"

# The columns of evidence_reports.parquet (docs/evidence_reports_spec.md section 1), in order; fixed there.
try:  # the shared report dataclass is the authority on the column set; the literal below is the fallback
    from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS as EVIDENCE_REPORT_COLUMNS
except ImportError:  # pragma: no cover
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


def map_relation_type(relation_type: str, partner_entity_type: str | None = None) -> tuple[str, int] | None:
    """(relation, report_value) for a PubTator3 relation type, or None when the type is not a perturbation-symptom relation.

    stimulate and inhibit map to induces only on a Disease-Gene row (partner_entity_type "Gene"), where they encode the
    sign of the gene's correlation with the disease; with any other partner they are chemical-gene relations and None.
    """
    relation_type_lower = relation_type.strip().lower()
    if relation_type_lower in DISEASE_GENE_RELATION_TYPE_SIGNS:
        return ("induces", 1) if partner_entity_type == GENE_ENTITY_TYPE else None
    relation = RELATION_TYPE_TO_RELATION.get(relation_type_lower)
    return None if relation is None else (relation, 1)


def literature_source_name(source_prefix: str, relation_type: str) -> str:
    """The source string the assembler parses: "<prefix>-<relation type>", for example "PubTator3-cause" or "CTD-curated"."""
    return f"{source_prefix}-{relation_type.strip().lower()}"


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
    perturbation_sign: float = LOSS_OF_FUNCTION_SIGN  # sign of the gene perturbation; drugs carry none here
    direction_demoted: bool = False  # a directional relation demoted to associated_with on a mixed-polarity descriptor


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
        mapped = map_relation_type(relation_type_lower, partner_type)
        if mapped is None:
            if relation_type_lower in CHEMICAL_GENE_RELATION_TYPES:
                counts[f"dropped_chemical_gene_relation_type:{relation_type_lower}"] += 1
            else:
                counts[f"dropped_relation_type_not_mapped:{relation_type_lower}"] += 1
            continue
        relation, report_value = mapped
        perturbation_sign = DISEASE_GENE_RELATION_TYPE_SIGNS.get(relation_type_lower, LOSS_OF_FUNCTION_SIGN)
        if relation_type_lower in DISEASE_GENE_RELATION_TYPE_SIGNS:
            counts[f"disease_gene_relation_type_kept:{relation_type_lower}"] += 1
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
            descriptor_relation, demoted = demote_direction_on_mixed_polarity_descriptor(relation, descriptor.level)
            if demoted:
                counts[f"demoted_to_associated_with_on_mixed_polarity_descriptor:{relation_type_lower}"] += 1
            for perturbation_type, perturbation_id, perturbation_label, ncbi_gene_id, chemical_mesh_id, match_method in partners:
                kept.append(
                    PubtatorRelationRow(
                        pmid=pmid, relation_type=relation_type_lower, relation=descriptor_relation, report_value=report_value, entity_role=entity_role,
                        perturbation_type=perturbation_type, perturbation_id=perturbation_id, perturbation_label=perturbation_label,
                        perturbation_ncbi_gene_id=ncbi_gene_id, chemical_mesh_id=chemical_mesh_id, chemical_match_method=match_method,
                        symptom=descriptor.target_symptom, mesh_descriptor=descriptor.mesh_descriptor, mesh_descriptor_level=descriptor.level,
                        perturbation_sign=perturbation_sign, direction_demoted=demoted,
                    )
                )
                counts["rows_kept"] += 1
    return kept, counts


def demote_direction_on_mixed_polarity_descriptor(relation: str, descriptor_level: str) -> tuple[str, bool]:
    """(relation to record, demoted) for a relation on a descriptor of the given crosswalk level: induces and relieves on a
    mixed-polarity diagnosis (Bipolar Disorder) become associated_with, since the diagnosis has the opposite symptom as a
    cardinal feature and a treat or cause relation on it carries no sign for the mania symptom."""
    if descriptor_level == MIXED_POLARITY_DESCRIPTOR_LEVEL and relation in ("induces", "relieves"):
        return ASSOCIATED_WITH_RELATION, True
    return relation, False


def _gene_partners(partner_ids: list[str], ncbi_gene_id_to_symbol: dict[int, str | list[str]], counts: Counter) -> list[tuple]:
    """One partner per graph symbol of each NCBI id; a value that is a list (two graph nodes sharing one NCBI id, such as a
    gene present under its current and its previous symbol) yields one partner per node."""
    partners: list[tuple] = []
    for identifier in partner_ids:
        try:
            ncbi_gene_id = int(identifier)
        except ValueError:
            counts["dropped_gene_identifier_not_integer"] += 1
            continue
        symbols = ncbi_gene_id_to_symbol.get(ncbi_gene_id)
        if not symbols:
            counts["dropped_gene_not_in_graph"] += 1
            continue
        for symbol in (symbols if isinstance(symbols, list) else [symbols]):
            partners.append(("gene", symbol, symbol, ncbi_gene_id, None, None))
        if isinstance(symbols, list) and len(symbols) > 1:
            counts["rows_with_gene_id_shared_by_several_nodes"] += 1
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


def gene_perturbation_nodes_json(gene_symbol: str, perturbation_sign: float = LOSS_OF_FUNCTION_SIGN) -> str:
    """The graph node of a gene with its perturbation sign: the sign of the gene's correlation with the disease for a
    stimulate or inhibit row, else the loss-of-function sign of the monogenic class as a documented convention
    (model_description says which). Drugs get an empty list and are joined to their ChEMBL targets by STITCH flat id
    by the assembler, as SIDER rows are."""
    return json.dumps([[f"GENE:{gene_symbol}", float(perturbation_sign), 1.0]])


def pubtator_limitations_text(row: PubtatorRelationRow, dated: bool) -> str:
    clauses = ["machine-extracted relation from PubTator3, no study design read"]
    if row.relation == ASSOCIATED_WITH_RELATION:
        clauses.append("undirected association")
    if row.mesh_descriptor_level in ("diagnosis", MIXED_POLARITY_DESCRIPTOR_LEVEL):
        clauses.append("diagnosis-level descriptor, not the symptom itself")
    if row.direction_demoted:
        clauses.append("diagnosis covers both mood poles, direction not attributable to the symptom")
    if row.perturbation_type == "gene" and row.relation_type in DISEASE_GENE_RELATION_TYPE_SIGNS:
        clauses.append(f"gene-disease correlation, perturbation sign {row.perturbation_sign:+.0f} from the relation type")
    elif row.perturbation_type == "gene":
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
        if row.perturbation_type == "gene" and row.relation_type in DISEASE_GENE_RELATION_TYPE_SIGNS:
            model = f"human gene mention; PubTator3 {row.relation_type}; perturbation sign {row.perturbation_sign:+.0f} from the relation type"
        elif row.perturbation_type == "gene":
            model = f"human gene mention; PubTator3 {row.relation_type}; sign set to loss-of-function by convention"
        else:
            model = f"human; drug mention; PubTator3 {row.relation_type}; targets supplied by the assembler from the SIDER map"
        records.append(
            literature_report_record(
                source=literature_source_name(PUBTATOR_SOURCE_NAME, row.relation_type),
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
                perturbation_nodes=gene_perturbation_nodes_json(row.perturbation_id, row.perturbation_sign) if row.perturbation_type == "gene" else "[]",
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
    for column_name in EVIDENCE_REPORT_COLUMNS:  # columns added to the shared dataclass after this module was written get their defaults
        if records and column_name not in records[0]:
            default_value = 0.0 if column_name.startswith("rubric_") else ""
            for record in records:
                record.setdefault(column_name, default_value)
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
        "rows_by_source": {key: int(value) for key, value in table.source.value_counts().sort_index().items()},
        "rows_demoted_to_associated_with_on_mixed_polarity_descriptor": int(((table.relation == ASSOCIATED_WITH_RELATION) & table.relation_type.isin(["cause", "treat", "prevent", "stimulate", "inhibit", "marker/mechanism", "therapeutic"])).sum()) if table.relation_type.notna().any() else 0,
        "gene_rows_by_relation_type": {key: int(value) for key, value in table.loc[table.perturbation_type == "gene", "relation_type"].value_counts().sort_index().items()} if table.relation_type.notna().any() else {},
        "rows_by_descriptor_level": {key: int(value) for key, value in table.mesh_descriptor_level.value_counts().sort_index().items()},
        "rows_by_descriptor": {key: int(value) for key, value in table.mesh_descriptor.value_counts().sort_index().items()},
        "rows_by_chemical_match_method": {key: int(value) for key, value in table.chemical_match_method.dropna().value_counts().sort_index().items()},
        "rows_dated": int(table.evidence_date.notna().sum()),
        "rows_by_symptom_and_perturbation_type": {f"{symptom}:{perturbation_type}": int(count) for (symptom, perturbation_type), count in table.groupby(["symptom", "perturbation_type"]).size().items()},
    }


def load_pubtator_relations(path: Path) -> pd.DataFrame:
    """The filtered table written by experiments/fetch_pubtator_relations.py."""
    return pd.read_parquet(path)
