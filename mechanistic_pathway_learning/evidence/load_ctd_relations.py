"""Evidence class E3, CTD half: curated chemical-disease associations with direction, one report per cited paper
(design section 4.2; docs/evidence_reports_spec.md section 1). Supersedes the CTD stub in load_literature_predications.py.

CTD_chemicals_diseases.tsv.gz (ctdbase.org/reports/) has the columns ChemicalName, ChemicalID (MeSH), CasRN,
DiseaseName, DiseaseID (MESH: or OMIM:), DirectEvidence, InferenceGeneSymbol, InferenceScore, OmimIDs, PubMedIDs.
Only rows with DirectEvidence are curated statements: "marker/mechanism" (the chemical correlates with or plays a
role in the disease) maps to induces and "therapeutic" to relieves; a row may carry both, separated by "|", and then
yields one report per evidence value. Rows with empty DirectEvidence are inferred through a gene and are dropped
with a count. PubMedIDs is a "|" list; each cited paper becomes its own report so the reliability model can count
papers, and a direct-evidence row that cites no paper is kept as one undated report with no reference.

A direct-evidence statement on a descriptor of crosswalk level mixed_polarity_diagnosis (Bipolar Disorder D001714)
is demoted to associated_with with a count and a limitation clause, as in load_pubtator_relations. One SIDER drug is
kept per CTD chemical (sider_chemical_matching.collapse_matches_to_one_drug_per_mesh_id), so a paper is never counted
under a second STITCH id for the same chemical. The source is "CTD-curated", the name the assembler's predication
parser expects.

Literature reports are soft priors only (grades C to E): never evaluation positives, and an unobserved pair is
unlabelled, not negative.
"""
from __future__ import annotations

import csv
import gzip
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.load_pubtator_relations import demote_direction_on_mixed_polarity_descriptor, literature_report_record, literature_source_name, reports_dataframe
from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import SymptomMeshDescriptor
from mechanistic_pathway_learning.evidence.sider_chemical_matching import ChemicalMatch, collapse_matches_to_one_drug_per_mesh_id, match_ctd_chemical_name

CTD_SOURCE_NAME = "CTD"
CTD_CURATED_SOURCE = literature_source_name(CTD_SOURCE_NAME, "curated")  # "CTD-curated"
DIRECT_EVIDENCE_TO_RELATION: dict[str, str] = {"marker/mechanism": "induces", "therapeutic": "relieves"}
CTD_COLUMNS: tuple[str, ...] = ("ChemicalName", "ChemicalID", "CasRN", "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubMedIDs")


def iterate_ctd_rows(path_or_handle) -> Iterator[dict[str, str]]:
    """Data rows of the CTD file as dicts keyed by CTD_COLUMNS; the "#" header block is skipped."""
    if isinstance(path_or_handle, (str, Path)):
        with gzip.open(path_or_handle, "rt", encoding="utf-8", newline="") as handle:
            yield from _iterate_ctd_handle(handle)
    else:
        yield from _iterate_ctd_handle(path_or_handle)


def _iterate_ctd_handle(handle) -> Iterator[dict[str, str]]:
    for row in csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
        if not row or row[0].startswith("#"):
            continue
        padded = row + [""] * (len(CTD_COLUMNS) - len(row))
        yield dict(zip(CTD_COLUMNS, padded))


@dataclass(frozen=True)
class CtdRelationRow:
    chemical_mesh_id: str
    chemical_name: str
    stitch_flat_id: str
    drug_name: str
    chemical_match_method: str
    disease_mesh_id: str
    disease_name: str
    direct_evidence: str  # "marker/mechanism" or "therapeutic"
    relation: str
    pmid: str | None
    symptom: str
    mesh_descriptor_level: str
    direction_demoted: bool = False  # a directional statement demoted to associated_with on a mixed-polarity descriptor


def expand_ctd_row(row: dict[str, str], descriptor: SymptomMeshDescriptor, matches: list[ChemicalMatch], counts: Counter) -> list[CtdRelationRow]:
    """One CtdRelationRow per (SIDER match, direct-evidence value, PubMed id) of a CTD row already known to hit a target
    descriptor; empty DirectEvidence yields nothing and is counted as inferred."""
    evidence_values = [value.strip() for value in row["DirectEvidence"].split("|") if value.strip()]
    if not evidence_values:
        counts["dropped_inferred_no_direct_evidence"] += 1
        return []
    pmids = [pmid.strip() for pmid in row["PubMedIDs"].split("|") if pmid.strip()]
    if not pmids:
        counts["direct_evidence_rows_without_pubmed_id"] += 1
    expanded: list[CtdRelationRow] = []
    if len(evidence_values) > 1:
        counts["rows_with_both_direct_evidence_values"] += 1
    for evidence_value in evidence_values:
        relation = DIRECT_EVIDENCE_TO_RELATION.get(evidence_value)
        if relation is None:
            counts[f"dropped_direct_evidence_not_mapped:{evidence_value}"] += 1
            continue
        relation, demoted = demote_direction_on_mixed_polarity_descriptor(relation, descriptor.level)
        if demoted:
            counts[f"demoted_to_associated_with_on_mixed_polarity_descriptor:{evidence_value}"] += 1
        for match in matches:
            for pmid in pmids or [None]:
                expanded.append(
                    CtdRelationRow(
                        chemical_mesh_id=row["ChemicalID"], chemical_name=row["ChemicalName"], stitch_flat_id=match.stitch_flat_id,
                        drug_name=match.drug_name, chemical_match_method=match.match_method, disease_mesh_id=descriptor.mesh_descriptor,
                        disease_name=row["DiseaseName"], direct_evidence=evidence_value, relation=relation, pmid=pmid,
                        symptom=descriptor.target_symptom, mesh_descriptor_level=descriptor.level, direction_demoted=demoted,
                    )
                )
                counts["reports_kept"] += 1
    return expanded


def filter_ctd_rows(
    rows: Iterable[dict[str, str]],
    descriptors: dict[str, SymptomMeshDescriptor],
    sider_names_index: dict[str, list[tuple[str, str]]],
    chemical_matches_by_mesh_id: dict[str, list[ChemicalMatch]] | None = None,
) -> tuple[list[CtdRelationRow], Counter]:
    """Rows whose DiseaseID is a target descriptor and whose chemical is a SIDER drug, by ChemicalName first
    (ctd_chemical_name) and otherwise by a MeSH id already matched through PubTator3 (its own method is kept); one
    SIDER drug per chemical, the PubTator3 map's choice when it has one, else the collapse rule of sider_chemical_matching."""
    chemical_matches_by_mesh_id = chemical_matches_by_mesh_id or {}
    preferred_stitch_ids = {match.stitch_flat_id for matches in chemical_matches_by_mesh_id.values() for match in matches}
    kept: list[CtdRelationRow] = []
    counts: Counter = Counter()
    for row in rows:
        counts["ctd_rows_total"] += 1
        disease_id = row["DiseaseID"]
        descriptor = descriptors.get(disease_id[5:]) if disease_id.upper().startswith("MESH:") else None
        if descriptor is None:
            counts["dropped_disease_not_target_descriptor"] += 1
            continue
        counts["ctd_rows_with_target_descriptor"] += 1
        matches = match_ctd_chemical_name(row["ChemicalID"], row["ChemicalName"], sider_names_index)
        if not matches:
            matches = chemical_matches_by_mesh_id.get(row["ChemicalID"], [])
        if not matches:
            counts["dropped_chemical_not_matched_to_sider"] += 1
            continue
        if len(matches) > 1:
            matches, collapsed = collapse_matches_to_one_drug_per_mesh_id(matches, preferred_stitch_ids)
            counts["ctd_rows_with_chemical_collapsed_to_one_sider_drug"] += 1
        counts["ctd_rows_with_target_descriptor_and_sider_chemical"] += 1
        kept.extend(expand_ctd_row(row, descriptor, matches, counts))
    return kept, counts


def ctd_limitations_text(row: CtdRelationRow) -> str:
    clauses_prefix = ["diagnosis covers both mood poles, direction not attributable to the symptom"] if row.direction_demoted else []
    clauses = clauses_prefix + ["curated chemical-disease statement from CTD, study design not read"]
    if row.direct_evidence == "marker/mechanism":
        clauses.append("marker/mechanism covers correlation as well as causation")
    if row.mesh_descriptor_level == "diagnosis":
        clauses.append("diagnosis-level descriptor, not the symptom itself")
    clauses.append(f"drug matched to MeSH by {row.chemical_match_method}")
    if row.pmid is None:
        clauses.append("no PubMed reference")
    clauses.append("undated")
    clauses.append("soft prior only, never an evaluation positive")
    return "; ".join(clauses)


def ctd_rows_to_reports(rows: list[CtdRelationRow], publication_dates: dict[str, str] | None = None) -> pd.DataFrame:
    """The report table for the CTD rows; ordinal numbers a report within its (chemical, descriptor, symptom) group."""
    publication_dates = publication_dates or {}
    ordinals: Counter = Counter()
    records: list[dict] = []
    for row in rows:
        group_key = (row.chemical_mesh_id, row.disease_mesh_id, row.symptom)
        ordinals[group_key] += 1
        evidence_date = publication_dates.get(row.pmid) if row.pmid else None
        records.append(
            literature_report_record(
                source=CTD_CURATED_SOURCE,
                source_record_id=f"MESH:{row.chemical_mesh_id}",
                source_term_id=f"MESH:{row.disease_mesh_id}",
                source_term_label=row.disease_mesh_id,
                ordinal=ordinals[group_key],
                perturbation_type="drug",
                perturbation_id=row.stitch_flat_id,
                perturbation_label=row.drug_name,
                symptom=row.symptom,
                relation=row.relation,
                report_value=1,
                evidence_code=f"ctd_direct_evidence_{'marker_mechanism' if row.direct_evidence == 'marker/mechanism' else row.direct_evidence}",
                model_description=f"human or animal, unspecified; CTD curated {row.direct_evidence}; targets supplied by the assembler from the SIDER map",
                references=f"PMID:{row.pmid}" if row.pmid else "",
                pubmed_reference_count=1 if row.pmid else 0,
                evidence_date=evidence_date,
                limitations=ctd_limitations_text(row).replace("; undated", "") if evidence_date else ctd_limitations_text(row),
                perturbation_nodes="[]",
                extras={
                    "pmid": row.pmid,
                    "relation_type": row.direct_evidence,
                    "entity_role": "disease_second",
                    "perturbation_ncbi_gene_id": None,
                    "chemical_mesh_id": row.chemical_mesh_id,
                    "chemical_match_method": row.chemical_match_method,
                    "mesh_descriptor": row.disease_mesh_id,
                    "mesh_descriptor_level": row.mesh_descriptor_level,
                    "publication_year": int(evidence_date[:4]) if evidence_date else None,
                    "direct_evidence": row.direct_evidence,
                },
            )
        )
    return reports_dataframe(records)


def load_ctd_relations(path: Path) -> pd.DataFrame:
    """The filtered table written by experiments/fetch_ctd_chemical_disease.py."""
    return pd.read_parquet(path)
