"""Assemble the perturbation-symptom observation table (design sections 4.2 and 4.3).

Joins the monogenic records (HPO x Human-GEM) and the pharmacological records
(SIDER x ChEMBL) to the physiology graph and writes one row per observation:

  perturbation_id        gene symbol (E1) or SIDER STITCH flat id (E2)
  perturbation_type      "gene" or "drug"
  perturbation_label     readable name
  group_id               leakage group for splits: the gene itself, or the dominant ChEMBL target for a drug
  disease_cluster_id     second leakage group for genes: connected component of the gene-disease graph
                         over the target-symptom annotations (genes annotated to one disease share its
                         whole phenotype profile); drugs keep their group_id here
  symptom                target symptom from docs/symptom_crosswalk.csv
  relation               "induces" or "relieves"
  evidence_class         "monogenic" or "pharmacological"
  grade, weight          from assign_evidence_grades over the pair's positive reports plus its 0 of N patient-fraction
                         rows, which the version 0.4 loader counted as frequency-0 rows and which the grade weighting
                         keeps so its weights stay comparable (--weighting grade, the default), or
                         reliability_global_scale x reliability_posterior (--weighting reliability); 0 for a pair
                         with no such row, and 0 under reliability weighting for a pair with no positive report
  perturbation_nodes     JSON list of [node_id, sign, magnitude] in the graph
  label_frequency        SIDER frequency midpoint (E2) or largest HPO frequency midpoint (E1) when reported
  evidence_date          earliest availability date behind a monogenic pair: publication date of the cited PubMed
                         reference when docs/hpo_reference_publication_dates.json has it, else the OMIM biocuration
                         date (ISO string; null when Orphanet-only); the time-split axis of design section 6.1
  omim_entry_count, orpha_entry_count, annotation_row_count, annotation_patient_count, disease_identifiers
  distinct_reference_count, distinct_pubmed_reference_count  descriptive multiplicity of the HPO references behind a monogenic row
                         provenance of a monogenic record (null for drugs)
  source                 provenance string
  report_count, positive_report_count, negative_report_count
                         explicit reports behind the pair (evidence_reports.parquet), with value 1 and 0
  reliability_posterior  P(link real | reports) from the report-level Dawid-Skene fit, on every row whatever the weighting
  weighting              "grade" or "reliability", the same on every row

The table is an aggregation of evidence_reports.parquet, one row per report (evidence_reports.py), written
next to it: a pair is present when it has at least one explicit report, positive or negative, and absent
otherwise (unobserved pairs are unlabelled, not negative). Rows whose perturbation has no node in the
graph are written to a separate unmapped table so the gap is visible rather than silently dropped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import (
    DEFAULT_GRADE_A_POLICY,
    GRADE_A_POLICIES,
    EvidenceRecord,
    assign_evidence_grade,
    loss_weight_for_record,
)
from mechanistic_pathway_learning.evidence.load_open_targets_statements import (
    human_worded_statements,
    reliability_reports,
)
from mechanistic_pathway_learning.evidence.evidence_reliability_model import (
    ReportReliabilityFit,
    fit_weighted_dawid_skene,
    observation_weights_from_posterior,
    report_count_matrices,
)
from mechanistic_pathway_learning.evidence.evidence_reports import (
    SIDER_FREQUENCY_EVIDENCE_CODES,
    REPORT_COLUMNS,
    SOURCE_NAMES,
    UNJOINED_EVIDENCE_CODE,
    EvidenceReport,
    RubricWeightDefaults,
    load_rubric_weight_defaults,
    reports_to_dataframe,
    rubric_weights_for_table,
)
from mechanistic_pathway_learning.evidence.allele_mechanisms import (
    DEFAULT_GENE2PHENOTYPE_PATH,
    DEFAULT_MECHANISM_ADJUDICATIONS_PATH,
    DEFAULT_NEOMORPHIC_ALLELES_PATH,
    HYPERMORPH,
    NEOMORPH,
    UNTYPED,
    DiseaseMechanism,
    allele_of,
    apply_adjudications,
    disease_mechanisms,
    first_source_unless_absent,
    gene2phenotype_mechanisms,
    gene_of_perturbation,
    load_gene2phenotype,
    load_orphanet_omim_mappings,
    read_mechanism_adjudications,
    read_neomorphic_alleles,
)
from mechanistic_pathway_learning.evidence.gene_disease_association_types import genes_with_gain_of_function_associations_only
from mechanistic_pathway_learning.evidence.load_drug_label_events import NERVOUS_SYSTEM_ATC_PREFIXES, drug_label_reports, has_admitted_atc_code, load_sider_events, read_admitted_atc_prefixes
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    load_hpo_annotation_rows,
    load_hpo_is_a_parents_from_obo,
    load_reference_publication_dates,
    monogenic_evidence_reports,
    parse_genes_to_phenotype,
    read_crosswalk_hpo_terms,
)
from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import DEFAULT_BINDING_SITES_PATH, GTOPDB_DIRECTORY, ComplexTargetNarrowing
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH,
    NON_PROTEIN_TARGETS_PATH,
    GraphNodeLookup,
    check_graph_compounds_have_no_mechanism,
    drug_targets_for_pubchem_cid,
    graph_compound_targets,
    has_dominant_target,
    load_chembl_caches,
    load_drugs_acting_as_graph_compounds,
)


def gene_node_lookup(nodes: pd.DataFrame) -> dict[str, str]:
    genes = nodes[nodes.node_type == "gene"]
    return {symbol: node_id for symbol, node_id in zip(genes.gene_symbol, genes.node_id) if isinstance(symbol, str) and symbol}


DEFAULT_DISEASE_CLUSTER_MAX_GENES = 20


def disease_cluster_ids(records: list[EvidenceRecord], max_genes_per_linking_entry: int | None = DEFAULT_DISEASE_CLUSTER_MAX_GENES) -> dict[str, str]:
    """Gene -> cluster id, where genes sharing any disease identifier are in one cluster (union-find).

    Disease entries annotated to max_genes_per_linking_entry or more genes among the records do not link
    genes: on the full graph they are Orphanet group-level entries (non-specific early-onset epileptic
    encephalopathy, familial dilated cardiomyopathy, systemic lupus erythematosus and the like), which are
    heterogeneous groups rather than one disease, and chaining through them joined a third of all genes
    into one cluster. Their annotations stay in the evidence; only their linking is dropped, which is
    recorded in the summary (open question 11 of the design). The cluster is named after its
    alphabetically first gene. Genes with no shared disease form singletons.
    """
    parent: dict[str, str] = {}

    def find(item: str) -> str:
        while parent.setdefault(item, item) != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(first: str, second: str) -> None:
        root_first, root_second = find(first), find(second)
        if root_first != root_second:
            parent[max(root_first, root_second)] = min(root_first, root_second)

    genes_by_disease: dict[str, set[str]] = defaultdict(set)
    for record in records:
        find(record.perturbation_identifier)
        for disease_id in record.disease_identifiers:
            genes_by_disease[disease_id].add(record.perturbation_identifier)
    for disease_id, genes in genes_by_disease.items():
        if max_genes_per_linking_entry is not None and len(genes) >= max_genes_per_linking_entry:
            continue
        ordered = sorted(genes)
        for gene in ordered[1:]:
            union(ordered[0], gene)
    return {gene: "cluster:" + find(gene) for gene in list(parent)}


WEIGHTINGS = ("grade", "reliability")
DEFAULT_WEIGHTING = "grade"
POSTERIOR_QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)


@dataclass
class AssembledEvidence:
    observations: pd.DataFrame  # evidence_records.parquet
    unmapped: pd.DataFrame  # unmapped_records.parquet
    reports: pd.DataFrame  # evidence_reports.parquet, rubric_weight filled
    reliability_fit: ReportReliabilityFit
    reports_dropped_unknown_provenance: int = 0  # genes_to_phenotype rows whose disease prefix is neither OMIM: nor ORPHA:
    genes_seeded_as_gain_of_function: tuple[str, ...] = ()  # genes seeded with sign +1 (seed_gain_of_function_genes_positive), sorted
    allele_summary: dict | None = None  # allele_perturbations: the alleles by class and the basis of each mechanism (decision 9.54)
    complex_target_narrowing: dict | None = None  # decision 8.11: per drug, the ChEMBL complex targets narrowed to GtoPdb's subunits


def disease_cluster_ids_from_reports(reports: pd.DataFrame, max_genes_per_linking_entry: int | None = DEFAULT_DISEASE_CLUSTER_MAX_GENES) -> dict[str, str]:
    """disease_cluster_ids over the monogenic reports: every disease entry of a gene links it, positive and negative reports alike."""
    diseases_by_gene: dict[str, set[str]] = defaultdict(set)
    monogenic = reports[reports.evidence_class == "monogenic"] if len(reports) else reports
    for perturbation_id, disease_id in zip(monogenic.perturbation_id, monogenic.source_record_id):
        diseases_by_gene[gene_of_perturbation(perturbation_id)].add(disease_id)  # every allele of a gene links through the gene
    records = [EvidenceRecord(gene_symbol, "", "induces", "monogenic", disease_identifiers=sorted(disease_ids)) for gene_symbol, disease_ids in diseases_by_gene.items()]
    return disease_cluster_ids(records, max_genes_per_linking_entry)


# The classes a pair's reports may carry, best first. A pair's primary class is the best one present, whatever the
# order the reports arrived in, so a label pair stays a label pair however many soft sources join it later.
# "human_association" (grade C) and "animal" (grade D) are the soft classes the design names and the trial-arm,
# FAERS and IMPC sources use: they enter the reliability fit and the priors and never the grade of a label.
EVIDENCE_CLASS_RANK: dict[str, int] = {"monogenic": 0, "pharmacological": 1, "open_targets": 2, "human_association": 3, "animal": 4, "literature": 5}
SOFT_EVIDENCE_CLASSES: tuple[str, ...] = ("human_association", "animal")


def primary_report(pair_reports: pd.DataFrame) -> pd.Series:
    """The first report of the best class the pair carries (EVIDENCE_CLASS_RANK); the class that decides the grade."""
    ranks = pair_reports.evidence_class.map(lambda name: EVIDENCE_CLASS_RANK.get(str(name), len(EVIDENCE_CLASS_RANK)))
    return pair_reports.loc[ranks.idxmin()]


def grade_weighting_rows(pair_reports: pd.DataFrame) -> pd.DataFrame:
    """The reports the grade weighting reads for one pair: positive reports plus monogenic 0 of N patient-fraction rows.

    A phenotype.hpoa row stating that 0 of N patients had the symptom is an absence report (report_value 0)
    in the report table and the reliability fit. The version 0.4 loader (monogenic_evidence_records) read the
    same row as a presence row with frequency 0, which the frequency clip turned into the floor weight, and
    the default weighting must reproduce those weights exactly for the running sweeps, so the grade and the
    weight are still computed over these rows. Whether the pair is a positive label is decided by
    positive_report_count, which counts report_value 1 only. The rule is the monogenic loader's: a trial arm with
    no event, which also has frequency 0 and a denominator, is an absence and nothing else.
    """
    zero_fraction = (pair_reports.frequency == 0.0) & pair_reports.frequency_denominator.notna() & (pair_reports.evidence_class == "monogenic")
    return pair_reports[(pair_reports.report_value == 1) | zero_fraction]


def evidence_record_from_positive_reports(pair_reports: pd.DataFrame) -> EvidenceRecord:
    """The EvidenceRecord that assign_evidence_grade and loss_weight_for_record read for one pair, from its grade_weighting_rows.

    Monogenic: provenance counts, the clinical-synopsis flag and the version 0.3 case-series proxy over the
    positive reports' disease entries, the largest positive frequency, the summed patient denominators, the
    distinct references and the earliest date. Pharmacological: one label event per preferred term; the record
    carries the label frequency (mean of that term's frequency reports) of the term with the largest grade
    weight, because the downstream loader took the maximum weight over the one-row-per-term table of earlier
    revisions and the default weighting has to reproduce it. A pair with no positive report yields a record
    with no provenance (grade C) or, for a drug, no dominant target (grade E).
    """
    first = primary_report(pair_reports)
    positive = grade_weighting_rows(pair_reports)
    positive = positive[positive.evidence_class == first.evidence_class]  # the grade reads the pair's primary class only; external classes (literature) enter the reliability fit, not the grade
    if first.evidence_class in SOFT_EVIDENCE_CLASSES:
        # a pair carried by trial arms, FAERS or a knockout alone: grade C or D by class, never a label (experiment_data.DEFAULT_LABEL_GRADES)
        return EvidenceRecord(first.perturbation_id, first.symptom, first.relation, str(first.evidence_class),
                              source=";".join(sorted(set(positive.source))) if len(positive) else str(first.source))
    if first.evidence_class == "monogenic":
        disease_ids = sorted(set(positive.source_record_id))
        frequencies = positive.frequency.dropna()
        denominators = positive.frequency_denominator.dropna()
        cited = {reference for references in positive.references for reference in str(references).split(";") if reference}
        dates = positive.evidence_date.dropna()
        return EvidenceRecord(
            perturbation_identifier=first.perturbation_id, symptom_identifier=first.symptom, relation=first.relation, evidence_class="monogenic",
            source="HPO genes_to_phenotype; diseases=" + ";".join(disease_ids),
            independent_case_series_count=len(disease_ids), has_omim_clinical_synopsis=any(disease_id.startswith("OMIM:") for disease_id in disease_ids),
            disease_identifiers=disease_ids,
            omim_entry_count=sum(1 for disease_id in disease_ids if disease_id.startswith("OMIM:")),
            orpha_entry_count=sum(1 for disease_id in disease_ids if disease_id.startswith("ORPHA:")),
            annotation_row_count=int(len(positive)),
            max_annotation_frequency=float(frequencies.max()) if len(frequencies) else None,
            annotation_patient_count=int(denominators.sum()) if len(denominators) else None,
            association_type_known=bool((positive.association_type.fillna("").astype(str).isin(["", "unknown"]) == False).any()) if "association_type" in positive.columns else False,  # noqa: E712
            causal_association_count=int(positive.rubric_causal_association.fillna(0).sum()) if "rubric_causal_association" in positive.columns else 0,
            evidence_available_date=date.fromisoformat(min(dates)) if len(dates) else None,
            evidence_date_source=str(positive.loc[positive.evidence_date == min(dates), "evidence_date_source"].iloc[0]) if len(dates) and "evidence_date_source" in positive.columns else "",
            distinct_reference_count=len(cited), distinct_pubmed_reference_count=sum(1 for reference in cited if reference.startswith("PMID:")),
        )
    if first.evidence_class == "literature":
        return EvidenceRecord(first.perturbation_id, first.symptom, first.relation, "literature", source=";".join(sorted(set(positive.source))) if len(positive) else str(first.source),
                              predication_type=predication_type_for_literature_sources(set(positive.source), first.relation))
    source_label = pharmacological_source_label(positive.source) if len(positive) else "SIDER 4.1"
    if not len(positive):
        return EvidenceRecord(first.perturbation_id, first.symptom, first.relation, "pharmacological", source=source_label, cns_penetrant=False, has_dominant_target=False)
    best_record: EvidenceRecord | None = None
    best_weight = -1.0
    for _, term_reports in positive.groupby("source_term_id", sort=False):
        term_frequencies = [float(value) for value in term_reports[term_reports.evidence_code.isin(SIDER_FREQUENCY_EVIDENCE_CODES)].frequency if not pd.isna(value)]
        term_frequency = sum(term_frequencies) / len(term_frequencies) if term_frequencies else None
        record = EvidenceRecord(first.perturbation_id, first.symptom, first.relation, "pharmacological", source=source_label, cns_penetrant=True, has_dominant_target=True, label_event_frequency=term_frequency)
        weight = loss_weight_for_record(record)
        if weight > best_weight:
            best_record, best_weight = record, weight
    return best_record


def pharmacological_source_label(sources) -> str:
    """"SIDER 4.1" when every source is a SIDER table (the string earlier revisions wrote), otherwise the sorted source names joined."""
    names = sorted(set(str(source) for source in sources))
    return "SIDER 4.1" if names and all(name.startswith("SIDER") for name in names) else ";".join(names)


def predication_type_for_literature_sources(sources: set[str], relation: str) -> str:
    """The soft-prior weight class of assign_evidence_grades.DEFAULT_PREDICATION_WEIGHTS for a literature-only pair, read
    from the extraction type the source name carries ("PubTator3-cause", "CTD-curated"): a cause or a curated
    marker/mechanism statement is CAUSES, a treat, prevent or curated therapeutic statement AFFECTS, an extracted
    association or correlation (associate, stimulate, inhibit, the correlation types) ASSOCIATED_WITH; COMENTION is
    kept for a source that names no extraction type at all."""
    names = {source.split("-", 1)[1].lower() if "-" in source else source.lower() for source in sources}
    if "cause" in names or ("curated" in names and relation == "induces"):
        return "CAUSES"
    if names & {"treat", "prevent", "curated"}:
        return "AFFECTS"
    if names & {"associate", "association", "stimulate", "inhibit", "positive_correlation", "negative_correlation", "positive_correlate", "negative_correlate"}:
        return "ASSOCIATED_WITH"
    return "COMENTION"


def targets_from_report_descriptions(pair_reports: pd.DataFrame) -> str:
    """The "CHEMBLid:ACTION;..." target string an external drug report carries in its model_description, for drugs the SIDER map does not know."""
    for description in pair_reports.model_description.astype(str):
        marker = "ChEMBL targets "
        if marker in description:
            return description.split(marker, 1)[1].strip()
    return ""


def pharmacological_label_frequency(positive: pd.DataFrame) -> float | None:
    """Largest over preferred terms of the mean treatment-midpoint of that term's frequency reports (what the downstream maximum over the earlier one-row-per-term table gave)."""
    term_means = []
    for _, term_reports in positive.groupby("source_term_id", sort=False):
        values = [float(value) for value in term_reports[term_reports.evidence_code.isin(SIDER_FREQUENCY_EVIDENCE_CODES)].frequency if not pd.isna(value)]
        if values:
            term_means.append(sum(values) / len(values))
    return max(term_means) if term_means else None


def aggregate_reports_to_observations(
    reports: pd.DataFrame,
    grade_a_policy: str,
    cluster_by_gene: dict[str, str],
    metabolic_node_ids: set[str],
    drug_targets_by_perturbation: dict[str, str] | None = None,
) -> pd.DataFrame:
    """One observation row per (perturbation_id, symptom, relation) with at least one explicit report, in order of first report.

    Grade and weight come from evidence_record_from_positive_reports over grade_weighting_rows; a pair with no
    such row has weight 0 and grade C (monogenic) or E (pharmacological) and stays in the table so the absence
    claim is visible. Provenance counts and disease identifiers run over all explicit reports; label_frequency,
    patient count, references and dates over the grade weighting rows. positive_report_count counts
    report_value 1 only. drug_targets_by_perturbation gives, per drug, the sorted "CHEMBLid:ACTION;..." string
    that names its group and its source line.
    """
    rows: list[dict] = []
    if not len(reports):
        return pd.DataFrame(rows)
    for (perturbation_id, symptom, relation), pair_reports in reports.groupby(["perturbation_id", "symptom", "relation"], sort=False):
        first = primary_report(pair_reports)
        positive = pair_reports[pair_reports.report_value == 1]
        graded = grade_weighting_rows(pair_reports)
        record = evidence_record_from_positive_reports(pair_reports)
        grade = assign_evidence_grade(record, grade_a_policy)
        weight = loss_weight_for_record(record, grade_a_policy=grade_a_policy) if len(graded) else 0.0
        node_ids = [node_id for node_id, _, _ in json.loads(first.perturbation_nodes)] if first.perturbation_nodes else []
        row = {
            "perturbation_id": perturbation_id, "perturbation_type": first.perturbation_type, "perturbation_label": first.perturbation_label,
            "symptom": symptom, "relation": relation, "evidence_class": first.evidence_class, "grade": grade, "weight": float(weight),
            "in_metabolic_layer": any(node_id in metabolic_node_ids for node_id in node_ids),
            "report_count": int(len(pair_reports)), "positive_report_count": int(len(positive)), "negative_report_count": int(len(pair_reports) - len(positive)),
            "perturbation_nodes": first.perturbation_nodes,
        }
        if first.perturbation_type == "gene":
            # the disease identifiers and provenance counts read the monogenic reports only; a gene pair carried by a
            # knockout alone (class animal) has none and is grouped as a singleton cluster of its own gene
            monogenic_reports = pair_reports[pair_reports.evidence_class == "monogenic"]
            all_disease_ids = sorted(set(monogenic_reports.source_record_id))
            denominators = graded.frequency_denominator.dropna()
            other_sources = sorted(set(pair_reports.loc[pair_reports.evidence_class != "monogenic", "source"]))
            gene = gene_of_perturbation(perturbation_id)  # an allele's leakage group is its gene's (decision 9.54)
            row.update({
                "group_id": gene, "disease_cluster_id": cluster_by_gene.get(gene, "cluster:" + gene),
                "label_frequency": record.max_annotation_frequency,
                "source": ("HPO genes_to_phenotype; diseases=" + ";".join(all_disease_ids) if all_disease_ids else "")
                          + ("; " if all_disease_ids and other_sources else "") + ("; ".join(other_sources) if other_sources else ""),
                "omim_entry_count": sum(1 for disease_id in all_disease_ids if disease_id.startswith("OMIM:")),
                "orpha_entry_count": sum(1 for disease_id in all_disease_ids if disease_id.startswith("ORPHA:")),
                "annotation_row_count": int(pair_reports[["source_record_id", "source_term_id"]].drop_duplicates().shape[0]),
                "annotation_patient_count": int(denominators.sum()) if len(denominators) else None,
                "disease_identifiers": ";".join(all_disease_ids),
                "distinct_reference_count": record.distinct_reference_count, "distinct_pubmed_reference_count": record.distinct_pubmed_reference_count,
                "evidence_date": record.evidence_available_date.isoformat() if record.evidence_available_date else None,
                "evidence_date_source": record.evidence_date_source or None,
            })
        else:
            targets = (drug_targets_by_perturbation or {}).get(perturbation_id, "") or targets_from_report_descriptions(pair_reports)
            group_id = "|".join(part.split(":")[0] for part in targets.split(";") if part) if targets else perturbation_id
            row.update({
                "group_id": group_id, "disease_cluster_id": group_id,
                "label_frequency": pharmacological_label_frequency(graded),
                "source": pharmacological_source_label(pair_reports.source) + "; targets " + targets,
                "omim_entry_count": None, "orpha_entry_count": None, "annotation_row_count": None, "annotation_patient_count": None, "disease_identifiers": None,
                "distinct_reference_count": None, "distinct_pubmed_reference_count": None, "evidence_date": None, "evidence_date_source": None,
            })
        if "allele_class" in pair_reports.columns:
            classes = sorted({value for value in pair_reports.allele_class if isinstance(value, str) and value})
            row["allele_class"] = ";".join(classes)
        rows.append(row)
    column_order = ["perturbation_id", "perturbation_type", "perturbation_label", "group_id", "disease_cluster_id", "symptom", "relation", "evidence_class", "grade", "weight",
                    "label_frequency", "source", "in_metabolic_layer", "omim_entry_count", "orpha_entry_count", "annotation_row_count", "annotation_patient_count", "disease_identifiers",
                    "distinct_reference_count", "distinct_pubmed_reference_count", "evidence_date", "evidence_date_source", "perturbation_nodes", "report_count", "positive_report_count", "negative_report_count"]
    if "allele_class" in reports.columns:
        column_order.append("allele_class")
    return pd.DataFrame(rows, columns=column_order)


def append_external_reports(reports: pd.DataFrame, path: Path, defaults: RubricWeightDefaults) -> pd.DataFrame:
    """Append a report table produced outside the assembler (OnSIDES label statements, literature reports).

    Rows with a qualifies column are kept only when it is true. The shared columns are taken as they are, missing
    shared columns get their defaults, rubric_ feature columns are kept, and a rubric_weight the table already
    carries (the literature appraisal computes its own) is kept, otherwise the defaults apply. The appended rows
    enter the report table and the reliability fit as their own sources; whether they become labels is decided
    downstream by grade and label_grades.
    """
    external = pd.read_parquet(path)
    if "qualifies" in external.columns:
        external = external[external.qualifies == True]  # noqa: E712
    kept = [column for column in external.columns if column in REPORT_COLUMNS or column.startswith("rubric_")]
    external = external[kept].copy()
    for column in REPORT_COLUMNS:
        if column not in external.columns:
            external[column] = 0.0 if column.startswith("rubric_") else ("" if column not in ("frequency", "frequency_denominator", "onset", "sex", "evidence_date") else None)
    computed = rubric_weights_for_table(external, defaults)
    external["rubric_weight"] = external["rubric_weight"].where(external["rubric_weight"].notna(), computed) if "rubric_weight" in external.columns else computed
    external = external[[column for column in reports.columns if column in external.columns] + [column for column in external.columns if column not in reports.columns]]
    return pd.concat([reports, external], ignore_index=True, sort=False)


def fit_report_reliability(reports: pd.DataFrame, observations: pd.DataFrame, defaults: RubricWeightDefaults) -> ReportReliabilityFit:
    """Report-level Dawid-Skene fit with the observations as items (row order) and the sources that have reports, in SOURCE_NAMES order.

    Prevalence is fitted per evidence class (no source covers both gene and drug items, so a pooled value would
    hand the HPO agreement to every drug pair), and a source with no explicit negative report keeps its
    specificity at the prior mean, since silence is the only negative it has (evidence_reliability_model
    docstring). Both are recorded in per_source_summary and the summary JSON.
    """
    item_keys = [(perturbation_id, symptom, relation) for perturbation_id, symptom, relation in zip(observations.perturbation_id, observations.symptom, observations.relation)] if len(observations) else []
    item_groups = [str(evidence_class) for evidence_class in observations.evidence_class] if len(observations) else []
    present_sources = set(reports.source) if len(reports) else set()
    source_names = [source for source in SOURCE_NAMES if source in present_sources] + sorted(present_sources - set(SOURCE_NAMES))  # external sources (OnSIDES, literature) follow the built-in ones
    positive_weights, negative_weights, coverage = report_count_matrices(reports, item_keys, source_names, defaults.implicit_negative_weight)
    explicit_cells = explicit_report_cells(reports, item_keys, source_names)
    explicit_negative_weight = (negative_weights - defaults.implicit_negative_weight * coverage * (~explicit_cells)).sum(axis=0)
    fit = fit_weighted_dawid_skene(positive_weights, negative_weights, source_names, item_keys, item_groups=item_groups, estimate_specificity=explicit_negative_weight > 0)
    for column, source in enumerate(source_names):
        fit.per_source_summary[source] = {
            "items_covered": int(coverage[:, column].sum()),
            "items_with_explicit_reports": int(explicit_cells[:, column].sum()),
            "explicit_positive_weight": float(positive_weights[:, column].sum()),
            "explicit_negative_weight": float(explicit_negative_weight[column]),
            "implicit_negative_cells": int((coverage[:, column] * (~explicit_cells[:, column])).sum()),
            "sensitivity": float(fit.sensitivity[column]),
            "specificity": float(fit.specificity[column]),
            "specificity_held_at_prior": fit.specificity_hold_reasons.get(source),
            "degenerate": source in fit.degenerate_sources,
        }
    return fit


def reliability_weights(fit: ReportReliabilityFit, observations: pd.DataFrame, global_scale: float) -> np.ndarray:
    """Weights under --weighting reliability: global_scale x posterior for pairs with a positive report, 0 otherwise.

    Refused (ValueError) when the fit lists a degenerate source, because a positive report of such a source
    lowers the posterior and the weights would train against the evidence.
    """
    if fit.degenerate_sources:
        raise ValueError(f"--weighting reliability refused: fitted sensitivity + specificity <= 1 for {fit.degenerate_sources}; the fit is degenerate for these sources (see reliability.per_source in the summary)")
    return np.where(observations.positive_report_count > 0, observation_weights_from_posterior(fit.posterior, global_scale), 0.0)


def explicit_report_cells(reports: pd.DataFrame, item_keys: list[tuple[str, str, str]], source_names: list[str]) -> np.ndarray:
    """Boolean [num_items, num_sources]: the source has at least one explicit report on the item."""
    item_index = {key: index for index, key in enumerate(item_keys)}
    source_index = {name: index for index, name in enumerate(source_names)}
    explicit = np.zeros((len(item_keys), len(source_names)), dtype=bool)
    if len(reports):
        for perturbation_id, symptom, relation, source in zip(reports.perturbation_id, reports.symptom, reports.relation, reports.source):
            row, column = item_index.get((perturbation_id, symptom, relation)), source_index.get(source)
            if row is not None and column is not None:
                explicit[row, column] = True
    return explicit


def posterior_pair_categories(reports: pd.DataFrame, observations: pd.DataFrame, fit: ReportReliabilityFit) -> np.ndarray:
    """Category per observation row for the posterior summary: how the covering sources speak about the pair."""
    item_keys = [(perturbation_id, symptom, relation) for perturbation_id, symptom, relation in zip(observations.perturbation_id, observations.symptom, observations.relation)]
    positive_weights, negative_weights, coverage = report_count_matrices(reports, item_keys, fit.source_names, 0.0)
    categories = np.empty(len(item_keys), dtype=object)
    for row in range(len(item_keys)):
        has_positive, has_negative = positive_weights[row].sum() > 0, negative_weights[row].sum() > 0
        covering, positive_sources = int(coverage[row].sum()), int((positive_weights[row] > 0).sum())
        if has_positive and has_negative:
            categories[row] = "conflicting"
        elif not has_positive:
            categories[row] = "all_negative"
        elif covering >= 2 and positive_sources >= 2:
            categories[row] = "covered_by_two_positive_in_both"
        elif covering >= 2:
            categories[row] = "positive_in_one_silent_in_other"
        else:
            categories[row] = "covered_by_one_source_positive"
    return categories


def summarize_reports(reports: pd.DataFrame, fit: ReportReliabilityFit, defaults: RubricWeightDefaults | None = None,
                      observations: pd.DataFrame | None = None, reports_dropped_unknown_provenance: int = 0) -> dict:
    """The "reports" and "reliability" blocks of evidence_summary.json."""
    defaults = defaults or RubricWeightDefaults()
    pair_values = reports.groupby(["perturbation_id", "symptom", "relation"]).report_value.agg(["min", "max"]) if len(reports) else pd.DataFrame(columns=["min", "max"])
    posterior_quantiles: dict[str, dict] = {}
    if observations is not None and len(observations) and len(fit.posterior) == len(observations):
        categories = posterior_pair_categories(reports, observations, fit)
        for category in ("covered_by_two_positive_in_both", "positive_in_one_silent_in_other", "covered_by_one_source_positive", "conflicting", "all_negative"):
            values = fit.posterior[categories == category]
            posterior_quantiles[category] = {"pairs": int(len(values)), "mean": float(values.mean()) if len(values) else None,
                                             **{str(q): float(np.quantile(values, q)) for q in POSTERIOR_QUANTILES}} if len(values) else {"pairs": 0}
    return {
        "reports": {
            "rows": int(len(reports)),
            "by_source_and_value": {f"{source}|{value}": int(n) for (source, value), n in reports.groupby(["source", "report_value"]).size().items()} if len(reports) else {},
            "by_evidence_code": {code: int(n) for code, n in reports.groupby("evidence_code").size().items()} if len(reports) else {},
            "reports_without_hpoa_row": int((reports.evidence_code == UNJOINED_EVIDENCE_CODE).sum()) if len(reports) else 0,
            "reports_dropped_unknown_provenance": int(reports_dropped_unknown_provenance),
            "pairs_all_negative": int((pair_values["max"] == 0).sum()),
            "pairs_conflicting": int(((pair_values["min"] == 0) & (pair_values["max"] == 1)).sum()),
            "rubric_weight_quantiles": {str(q): float(reports.rubric_weight.quantile(q)) for q in POSTERIOR_QUANTILES} if len(reports) else {},
        },
        "reliability": {
            "source_names": list(fit.source_names),
            "prevalence": float(fit.prevalence),
            "prevalence_by_evidence_class": {group: float(value) for group, value in fit.prevalence_by_group.items()},
            "num_iterations": int(fit.num_iterations),
            "converged": bool(fit.converged),
            "weakly_identified": bool(fit.weakly_identified),
            "degenerate_sources": list(fit.degenerate_sources),
            "rubric_weight_defaults": defaults.as_dict(),
            "posterior_quantiles": posterior_quantiles,
            "per_source": fit.per_source_summary,
        },
    }


def bridge_chembl_parents_by_sider_cid(onsides_bridge_path: Path | None) -> dict[int, list[str]]:
    """SIDER PubChem CID -> the ChEMBL parents the OnSIDES bridge unified with it (empty without a bridge file).

    The SIDER route reaches ChEMBL through UniChem; for a CID UniChem leaves without a usable ChEMBL id, the parent
    the bridge found for the same drug by RxNorm, UNII or name is the fallback (docs/drug_targets_any_type.md).
    """
    if onsides_bridge_path is None or not Path(onsides_bridge_path).exists():
        return {}
    from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import load_ingredient_identifier_bridge

    parents_by_cid: dict[int, set[str]] = {}
    for bridge in load_ingredient_identifier_bridge(onsides_bridge_path).values():
        if bridge.unified_with_sider and bridge.sider_pubchem_cid is not None and bridge.chembl_parent:
            parents_by_cid.setdefault(int(bridge.sider_pubchem_cid), set()).add(bridge.chembl_parent)
    return {cid: sorted(parents) for cid, parents in parents_by_cid.items()}


def sider_drug_targets(pubchem_cid: int, caches, bridge_parents_by_cid: dict[int, list[str]], compounds_by_pubchem_cid: dict[int, list[str]],
                       compounds_by_chembl_id: dict[str, list[str]]) -> list:
    """Mechanism targets of a SIDER drug: the UniChem route, then the bridge parents when it finds no target, then the
    graph-compound fallback for a drug that is itself a graph metabolite."""
    drug_targets = drug_targets_for_pubchem_cid(pubchem_cid, caches)
    bridge_parents = bridge_parents_by_cid.get(int(pubchem_cid), [])
    if not drug_targets and bridge_parents:
        drug_targets = drug_targets_for_pubchem_cid(pubchem_cid, caches, extra_chembl_ids=bridge_parents)
    chembl_ids = list(caches.pubchem_to_chembl.get(str(pubchem_cid), [])) + bridge_parents
    return graph_compound_targets(drug_targets, int(pubchem_cid), chembl_ids, compounds_by_pubchem_cid, compounds_by_chembl_id)


def seed_monogenic_reports(monogenic_reports: list[EvidenceReport], node_by_symbol: dict[str, str], seed_gain_of_function_genes_positive: bool = False) -> set[str]:
    """Fill perturbation_nodes of every monogenic report with its gene's node, sign and magnitude 1; returns the genes seeded +1.

    Every gene is seeded with sign -1, a loss of function. With seed_gain_of_function_genes_positive, a gene whose every
    causal association behind these reports is Orphanet's gain-of-function type is seeded with sign +1 and its reports'
    model_description says so; a gene with causal associations of both kinds keeps sign -1, since a gene is one
    perturbation with one sign (docs/label_source_review.md, section 5).
    """
    gain_of_function_genes: set[str] = set()
    if seed_gain_of_function_genes_positive:
        causal_association_types_by_gene: dict[str, set[str]] = defaultdict(set)
        for report in monogenic_reports:
            if report.rubric_causal_association == 1.0:
                causal_association_types_by_gene[report.perturbation_id].add(report.association_type)
        gain_of_function_genes = genes_with_gain_of_function_associations_only(causal_association_types_by_gene)
    for report in monogenic_reports:
        seeded_as_gain = report.perturbation_id in gain_of_function_genes
        report.perturbation_nodes = json.dumps([[node_by_symbol[report.perturbation_id], 1.0 if seeded_as_gain else -1.0, 1.0]])
        if seeded_as_gain:
            report.model_description = report.model_description.replace("human loss-of-function", "human gain-of-function", 1)
    return gain_of_function_genes


def assign_alleles(monogenic_reports: list[EvidenceReport], mechanisms: dict[tuple[str, str], DiseaseMechanism], neomorphic_alleles: dict,
                   node_by_symbol: dict[str, str]) -> dict[str, tuple[str, str, str]]:
    """Move each monogenic report to the allele of its gene that its disease's mechanism names (decision 9.54,
    allele_mechanisms.allele_of) and seed it: a loss of function or an untyped allele -1 on the gene, a hypermorphic
    gain +1 on the gene, a neomorph on the targets its curated mechanism edges name (on a graph with entry nodes the
    loader seeds the allele's own node instead). Returns report_id -> (gene, allele class, basis of the mechanism)."""
    assigned: dict[str, tuple[str, str, str]] = {}
    for report in monogenic_reports:
        gene = report.perturbation_id
        mechanism = mechanisms.get((gene, report.source_record_id), DiseaseMechanism(UNTYPED, "untyped"))
        allele = allele_of(gene, mechanism, neomorphic_alleles)
        report.perturbation_id, report.perturbation_label = allele.perturbation_id, allele.label
        if allele.allele_class == NEOMORPH:
            curated = neomorphic_alleles[gene]
            seeds = [[node_by_symbol[edge.target_symbol], edge.sign, edge.magnitude] for edge in curated.edges if edge.target_symbol in node_by_symbol]
            if not seeds:
                raise ValueError(f"no mechanism target of {curated.allele_id} is a gene of this graph")
            report.perturbation_nodes = json.dumps(seeds)
            report.model_description = report.model_description.replace("human loss-of-function", f"human neomorphic allele, {curated.label}", 1)
        else:
            report.perturbation_nodes = json.dumps([[node_by_symbol[gene], allele.seed_sign, 1.0]])
            if allele.allele_class == HYPERMORPH:
                report.model_description = report.model_description.replace("human loss-of-function", "human gain-of-function", 1)
        assigned[report.report_id] = (gene, allele.allele_class, mechanism.basis)
    return assigned


def summarise_alleles(assigned: dict[str, tuple[str, str, str]], reports: list[EvidenceReport]) -> dict:
    """The alleles by class, the genes with more than one, and the basis behind each mechanism, over the reports."""
    alleles_by_gene: dict[str, set[str]] = defaultdict(set)
    class_of_allele: dict[str, str] = {}
    basis_counts: dict[str, int] = defaultdict(int)
    for report in reports:
        if report.report_id not in assigned:
            continue
        gene, allele_class, basis = assigned[report.report_id]
        alleles_by_gene[gene].add(report.perturbation_id)
        class_of_allele[report.perturbation_id] = allele_class
        basis_counts[basis.split(" ")[0].rstrip(";") if not basis.startswith("conflict") else "conflict"] += 1
    by_class: dict[str, int] = defaultdict(int)
    for allele_class in class_of_allele.values():
        by_class[allele_class] += 1
    several = sorted(gene for gene, alleles in alleles_by_gene.items() if len(alleles) > 1)
    return {"alleles_by_class": dict(sorted(by_class.items())), "genes": len(alleles_by_gene), "genes_with_several_alleles": len(several),
            "genes_with_several_alleles_list": several, "reports_by_mechanism_basis": dict(sorted(basis_counts.items())),
            "conflicting_pairs": sorted({f"{assigned[report.report_id][0]} {report.source_record_id}" for report in reports
                                         if report.report_id in assigned and assigned[report.report_id][2].startswith("conflict")}),
            "pairs_typed_by_the_second_source_only": len({(assigned[report.report_id][0], report.source_record_id) for report in reports
                                                          if report.report_id in assigned and assigned[report.report_id][2].startswith("gene2phenotype")}),
            "pairs_disputed_by_the_second_source": sorted({f"{assigned[report.report_id][0]} {report.source_record_id}: {assigned[report.report_id][2]}" for report in reports
                                                          if report.report_id in assigned and "; disputed by " in assigned[report.report_id][2]
                                                          and not assigned[report.report_id][2].startswith("adjudicated")}),
            "pairs_adjudicated": sorted({f"{assigned[report.report_id][0]} {report.source_record_id}: {assigned[report.report_id][1]}, {assigned[report.report_id][2]}" for report in reports
                                         if report.report_id in assigned and assigned[report.report_id][2].startswith("adjudicated")})}


def drug_target_description(drug_targets) -> str:
    return ";".join(f"{target.target_chembl_id}:{target.action_type}" for target in sorted(drug_targets, key=lambda target: target.target_chembl_id))


def assemble(
    crosswalk_path: Path,
    hpo_obo_path: Path,
    hpo_annotations_path: Path,
    graph_directory: Path,
    sider_directory: Path | None,
    chembl_directory: Path | None,
    max_drug_targets: int | None = 1,
    grade_a_policy: str = DEFAULT_GRADE_A_POLICY,
    phenotype_hpoa_path: Path | None = None,
    reference_publication_dates_path: Path | None = None,
    genes_to_disease_path: Path | None = None,
    orphadata_product6_path: Path | None = None,
    disease_cluster_max_genes: int | None = DEFAULT_DISEASE_CLUSTER_MAX_GENES,
    weighting: str = DEFAULT_WEIGHTING,
    rubric_weight_defaults: RubricWeightDefaults | None = None,
    reliability_global_scale: float = 1.0,
    extra_report_paths: list[Path] | None = None,
    onsides_bridge_path: Path | None = None,
    non_protein_targets_path: Path | None = NON_PROTEIN_TARGETS_PATH,
    drugs_acting_as_graph_compounds_path: Path | None = DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH,
    seed_gain_of_function_genes_positive: bool = False,
    admitted_atc_prefixes: tuple[str, ...] = NERVOUS_SYSTEM_ATC_PREFIXES,
    open_targets_scoped_rows_path: Path | None = None,
    allele_perturbations: bool = False,
    orphadata_product1_path: Path | None = None,
    neomorphic_alleles_path: Path | None = DEFAULT_NEOMORPHIC_ALLELES_PATH,
    gene2phenotype_path: Path | None = DEFAULT_GENE2PHENOTYPE_PATH,
    mechanism_adjudications_path: Path | None = DEFAULT_MECHANISM_ADJUDICATIONS_PATH,
    complex_target_narrowing=None,
) -> AssembledEvidence:
    """Build the report table, aggregate it to observations, fit the report-level reliability model and weight the observations.

    seed_gain_of_function_genes_positive seeds a gene with sign +1 when every causal association behind its reports is
    Orphanet's gain-of-function type (genes_with_gain_of_function_associations_only); every other gene keeps sign -1,
    as all genes do by default. admitted_atc_prefixes is the list of ATC groups and substances a SIDER drug must fall
    under, group N alone by default; the OnSIDES reports are built apart and take the same list there.

    allele_perturbations (decision 9.54) makes each gene's mechanism class its own perturbation instead: the reports of
    a disease typed gain of function (by Orphanet, or carried from the same disease's Orphanet record through Orphadata
    product 1, orphadata_product1_path) move to "<gene>|gain_of_function", seeded +1, or to the gene's curated neomorphic
    allele (neomorphic_alleles_path); every other report stays with the gene, seeded -1. Where Orphanet leaves a pair
    untyped, Gene2Phenotype's curated mechanism is the second source (gene2phenotype_path; None builds without it, which
    is recorded), and a pair the two dispute is settled by its cited adjudication where one exists
    (mechanism_adjudications_path). It replaces seed_gain_of_function_genes_positive, and the two are not combined.

    complex_target_narrowing (perturbation/gtopdb_drug_targets.ComplexTargetNarrowing, decision 8.11) narrows each
    SIDER drug's ChEMBL complex targets to the subunits the Guide to PHARMACOLOGY names before the drug is seeded;
    None seeds every subunit of a complex, as every build before it did."""
    if allele_perturbations and seed_gain_of_function_genes_positive:
        raise ValueError("allele_perturbations replaces seed_gain_of_function_genes_positive; pass one of them")
    if weighting not in WEIGHTINGS:
        raise ValueError(f"unknown weighting {weighting!r}; choose from {WEIGHTINGS}")
    rubric_weight_defaults = rubric_weight_defaults or RubricWeightDefaults()
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    node_by_symbol = gene_node_lookup(nodes)
    metabolic_genes = nodes[(nodes.node_type == "gene") & (nodes.get("in_metabolic_layer", False) == True)]  # noqa: E712
    metabolite_node_ids = set(nodes[nodes.node_type == "metabolite"].node_id)
    metabolic_node_ids = set(metabolic_genes.node_id) | metabolite_node_ids  # a drug seeded on metabolites acts in the metabolic layer
    metabolic_symbols = set(metabolic_genes.gene_symbol.dropna())
    node_lookup = GraphNodeLookup.from_nodes(nodes, non_protein_targets_path)
    crosswalk_terms = read_crosswalk_hpo_terms(crosswalk_path)
    symptom_to_hpo = {symptom: roots for symptom, (roots, _) in crosswalk_terms.items()}
    symptom_to_excluded = {symptom: excluded for symptom, (_, excluded) in crosswalk_terms.items()}
    unmapped: list[dict] = []

    parents = load_hpo_is_a_parents_from_obo(hpo_obo_path)
    publication_dates = load_reference_publication_dates(reference_publication_dates_path) if reference_publication_dates_path is not None and reference_publication_dates_path.exists() else None
    hpoa_rows_by_key = load_hpo_annotation_rows(phenotype_hpoa_path) if phenotype_hpoa_path is not None and phenotype_hpoa_path.exists() else None
    association_types = None
    orphanet_types = None
    if genes_to_disease_path is not None and genes_to_disease_path.exists():
        from mechanistic_pathway_learning.evidence.gene_disease_association_types import association_type_for, load_omim_association_types, load_orphanet_association_types

        omim_types = load_omim_association_types(genes_to_disease_path)
        orphanet_types = load_orphanet_association_types(orphadata_product6_path) if orphadata_product6_path is not None and orphadata_product6_path.exists() else None
        association_types = {key: association_type_for(key[0], key[1], omim_types, orphanet_types) for key in set(omim_types) | set(orphanet_types or {})}
    dropped_unknown_provenance: list[dict] = []
    report_list: list[EvidenceReport] = monogenic_evidence_reports(parse_genes_to_phenotype(hpo_annotations_path), symptom_to_hpo, parents, set(node_by_symbol), symptom_to_excluded,
                                                                   hpoa_rows_by_key, publication_dates, dropped_unknown_provenance, association_types=association_types)
    allele_assignment: dict[str, tuple[str, str, str]] = {}
    allele_summary = None
    if allele_perturbations:
        if orphanet_types is None or orphadata_product1_path is None or not Path(orphadata_product1_path).exists():
            raise ValueError("allele_perturbations needs Orphadata product 6 (the typed associations) and product 1 (the OMIM mappings)")
        mappings = load_orphanet_omim_mappings(Path(orphadata_product1_path))
        pairs = {(report.perturbation_id, report.source_record_id) for report in report_list}
        mechanisms = disease_mechanisms(orphanet_types, mappings, pairs)
        second_source_summary: dict = {"gene2phenotype": None}
        if gene2phenotype_path is not None:
            if not Path(gene2phenotype_path).exists():
                raise ValueError(f"{gene2phenotype_path} is missing: fetch it (scripts/fetch_raw_sources.sh) or pass gene2phenotype_path=None to build without the second source")
            gene2phenotype_records, record_counts = load_gene2phenotype(Path(gene2phenotype_path))
            mechanisms = first_source_unless_absent(mechanisms, gene2phenotype_mechanisms(gene2phenotype_records, mappings, pairs), pairs)
            second_source_summary = {"gene2phenotype": str(gene2phenotype_path), "gene2phenotype_sha256": hashlib.sha256(Path(gene2phenotype_path).read_bytes()).hexdigest(), **record_counts}
        if mechanism_adjudications_path is not None:
            adjudications = read_mechanism_adjudications(Path(mechanism_adjudications_path))
            mechanisms = apply_adjudications(mechanisms, adjudications)
            second_source_summary.update(mechanism_adjudications=str(mechanism_adjudications_path),
                                         mechanism_adjudications_sha256=hashlib.sha256(Path(mechanism_adjudications_path).read_bytes()).hexdigest(),
                                         adjudications_without_a_report=sorted(f"{gene} {disease}" for gene, disease in adjudications if (gene, disease) not in pairs))
        neomorphic_alleles = read_neomorphic_alleles(Path(neomorphic_alleles_path)) if neomorphic_alleles_path is not None else {}
        allele_assignment = assign_alleles(report_list, mechanisms, neomorphic_alleles, node_by_symbol)
        allele_summary = {**summarise_alleles(allele_assignment, report_list), **second_source_summary}
        gain_of_function_genes = {gene for gene, allele_class, _ in allele_assignment.values() if allele_class in (HYPERMORPH, NEOMORPH)}
    else:
        gain_of_function_genes = seed_monogenic_reports(report_list, node_by_symbol, seed_gain_of_function_genes_positive)

    drug_targets_by_perturbation: dict[str, str] = {}
    narrowing_by_perturbation: dict[str, dict] = {}
    if sider_directory is not None and chembl_directory is not None and (chembl_directory / "targets.json").exists():
        caches = load_chembl_caches(chembl_directory)
        bridge_parents_by_cid = bridge_chembl_parents_by_sider_cid(onsides_bridge_path)
        compounds_by_pubchem_cid, compounds_by_chembl_id = load_drugs_acting_as_graph_compounds(drugs_acting_as_graph_compounds_path) if drugs_acting_as_graph_compounds_path else ({}, {})
        check_graph_compounds_have_no_mechanism(compounds_by_chembl_id, caches.mechanisms_by_molecule)
        for event in load_sider_events(sider_directory, crosswalk_path):
            drug_targets = sider_drug_targets(event.pubchem_cid, caches, bridge_parents_by_cid, compounds_by_pubchem_cid, compounds_by_chembl_id)
            if not (has_dominant_target(drug_targets, max_drug_targets) and has_admitted_atc_code(event.atc_codes, admitted_atc_prefixes)):
                continue
            if complex_target_narrowing is not None:
                drug_targets, narrowing_record = complex_target_narrowing.narrow(drug_targets, event.pubchem_cid, event.drug_name)
                if narrowing_record["narrowed"] or narrowing_record["kept_whole"]:
                    narrowing_by_perturbation[event.stitch_flat_id] = {"drug": event.drug_name, **narrowing_record}
            perturbation_nodes = node_lookup.perturbation_nodes(drug_targets)
            targets = drug_target_description(drug_targets)
            if not perturbation_nodes:
                record = EvidenceRecord(event.stitch_flat_id, event.target_symptom, event.relation, "pharmacological", source=event.source,
                                        cns_penetrant=True, has_dominant_target=True, label_event_frequency=event.label_frequency)
                group_id = "|".join(sorted(target.target_chembl_id for target in drug_targets))
                unmapped.append({
                    "perturbation_id": event.stitch_flat_id, "perturbation_type": "drug", "perturbation_label": event.drug_name,
                    "group_id": group_id, "disease_cluster_id": group_id, "symptom": event.target_symptom, "relation": event.relation,
                    "evidence_class": "pharmacological", "grade": assign_evidence_grade(record), "weight": loss_weight_for_record(record),
                    "label_frequency": event.label_frequency, "source": f"{event.source}; targets " + targets,
                    "in_metabolic_layer": any(symbol in metabolic_symbols for target in drug_targets for symbol in target.gene_symbols),
                    "omim_entry_count": None, "orpha_entry_count": None, "annotation_row_count": None, "annotation_patient_count": None, "disease_identifiers": None, "distinct_reference_count": None, "distinct_pubmed_reference_count": None,
                    "evidence_date": None,
                })
                continue
            drug_targets_by_perturbation[event.stitch_flat_id] = targets
            model_description = ("human; drug indication; " if event.relation == "relieves" else "human; drug label; ") + "ChEMBL targets " + targets
            report_list.extend(drug_label_reports(event, model_description, json.dumps(perturbation_nodes)))

    reports = reports_to_dataframe(report_list)
    if allele_perturbations:
        reports["allele_class"] = reports.report_id.map(lambda report_id: allele_assignment[report_id][1] if report_id in allele_assignment else "")
        reports["allele_mechanism_basis"] = reports.report_id.map(lambda report_id: allele_assignment[report_id][2] if report_id in allele_assignment else "")
    reports["rubric_weight"] = rubric_weights_for_table(reports, rubric_weight_defaults)
    for extra_path in extra_report_paths or []:
        reports = append_external_reports(reports, extra_path, rubric_weight_defaults)
    # Amendment item 12d and decision 8.5 (the owner, 13 October 2026: "add it to the Dawid-Skene fit estimates";
    # option (a) chosen 14 October after the effect was measured). The class reaches the reliability fit and not
    # the grade, so its whole effect is on the fitted sensitivities and on each item's posterior.
    #
    # ONE GUARD IS INACTIVE HERE AND THAT IS THE CHOICE THE OWNER MADE. `reliability_reports` also skips a pair the
    # label selection has relabelled an evidenced negative, but that set is computed downstream in
    # evaluation/experiment_data.py and does not exist inside this function. Measured before choosing: the class
    # emits 14 reports, no existing source's sensitivity moves by as much as 1e-5, no item's posterior moves by
    # more than 0.0019, and exactly 4 of the touched items are pairs the selection flags `masks_a_negative`
    # (GABRA1 depressed_mood, HTR2A cognitive_impairment, HTR2A psychosis, SLC6A4 anxiety). Threading the set
    # across the module boundary would have guarded those 4 items and changed nothing measurable, so option (a)
    # was taken and the inactive guard is recorded here rather than left to be discovered.
    open_targets_report_count = 0
    if open_targets_scoped_rows_path is not None and Path(open_targets_scoped_rows_path).exists():
        items_with_a_positive_report = {
            (row.perturbation_id, row.symptom, row.relation)
            for row in reports[reports.report_value == 1].itertuples()}
        gene_perturbation_ids = set(reports[reports.perturbation_type == "gene"].perturbation_id)
        perturbation_of_target = {gene_symbol: gene_symbol for gene_symbol in gene_perturbation_ids}
        statements = human_worded_statements(pd.read_parquet(open_targets_scoped_rows_path))
        open_targets_rows = reliability_reports(statements, items_with_a_positive_report, perturbation_of_target)
        if open_targets_rows:
            open_targets_frame = pd.DataFrame(open_targets_rows)
            open_targets_frame["rubric_weight"] = rubric_weights_for_table(open_targets_frame, rubric_weight_defaults)
            reports = pd.concat([reports, open_targets_frame], ignore_index=True)
            open_targets_report_count = len(open_targets_frame)
        print(f"Open Targets target safety: {open_targets_report_count} reliability reports appended "
              f"from {len(statements)} curated statements (evidenced-negative guard inactive; decision 8.5 option a)")

    cluster_by_gene = disease_cluster_ids_from_reports(reports, disease_cluster_max_genes)
    observations = aggregate_reports_to_observations(reports, grade_a_policy, cluster_by_gene, metabolic_node_ids, drug_targets_by_perturbation)
    reliability_fit = fit_report_reliability(reports, observations, rubric_weight_defaults)
    if len(observations):
        observations["reliability_posterior"] = reliability_fit.posterior
        if weighting == "reliability":
            observations["weight"] = reliability_weights(reliability_fit, observations, reliability_global_scale)
        observations["weighting"] = weighting
    narrowing_summary = None
    if complex_target_narrowing is not None:
        narrowing_summary = {**complex_target_narrowing.provenance, "drugs_narrowed": sum(1 for record in narrowing_by_perturbation.values() if record["narrowed"]),
                             "by_drug": dict(sorted(narrowing_by_perturbation.items()))}
    return AssembledEvidence(observations, pd.DataFrame(unmapped), reports, reliability_fit, len(dropped_unknown_provenance), tuple(sorted(gain_of_function_genes)), allele_summary,
                             narrowing_summary)


def summarize_observations(observations: pd.DataFrame, unmapped: pd.DataFrame) -> dict:
    if not len(observations):
        return {"observations": 0, "unmapped": int(len(unmapped))}
    genes = observations[observations.perturbation_type == "gene"]
    return {
        "observations": int(len(observations)),
        "unmapped": int(len(unmapped)),
        "by_class_and_relation": {f"{cls}|{rel}": int(n) for (cls, rel), n in observations.groupby(["evidence_class", "relation"]).size().items()},
        "distinct_perturbations": int(observations.perturbation_id.nunique()),
        "distinct_groups": int(observations.group_id.nunique()),
        "distinct_disease_clusters": int(observations.disease_cluster_id.nunique()),
        "largest_disease_cluster_genes": int(genes.groupby("disease_cluster_id").perturbation_id.nunique().max()) if len(genes) else 0,
        "by_symptom": {symptom: int(n) for symptom, n in observations.groupby("symptom").size().items()},
        "by_grade": {grade: int(n) for grade, n in observations.groupby("grade").size().items()},
        "weight_quantiles": {str(q): float(observations.weight.quantile(q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)},
        "monogenic_rows_with_frequency": int(genes.label_frequency.notna().sum()) if len(genes) else 0,
        "monogenic_rows_with_date": int(genes.evidence_date.notna().sum()) if len(genes) and "evidence_date" in genes else 0,
        "monogenic_rows_with_zero_pubmed_references": int((genes.distinct_pubmed_reference_count == 0).sum()) if len(genes) and "distinct_pubmed_reference_count" in genes else 0,
        "monogenic_rows_with_two_or_more_pubmed_references": int((genes.distinct_pubmed_reference_count >= 2).sum()) if len(genes) and "distinct_pubmed_reference_count" in genes else 0,
        "monogenic_rows_dated_after_2015": int((genes.evidence_date.dropna() > "2015-12-31").sum()) if len(genes) and "evidence_date" in genes else 0,
        "disease_cluster_concentration_by_symptom": disease_cluster_concentration(genes) if len(genes) else {},
        "disease_entries_not_linking_clusters": non_linking_disease_entries(genes) if len(genes) else {},
        "pairs_all_negative": int((observations.positive_report_count == 0).sum()) if "positive_report_count" in observations else 0,
        "pairs_conflicting": int(((observations.positive_report_count > 0) & (observations.negative_report_count > 0)).sum()) if "negative_report_count" in observations else 0,
    }


def non_linking_disease_entries(gene_observations: pd.DataFrame, max_genes: int = DEFAULT_DISEASE_CLUSTER_MAX_GENES) -> dict[str, int]:
    """Disease entries annotated to max_genes or more of the genes in the table (the group-level entries that do not link clusters)."""
    genes_by_disease: dict[str, set[str]] = defaultdict(set)
    for gene, identifiers in zip(gene_observations.perturbation_id, gene_observations.disease_identifiers):
        for disease_id in str(identifiers or "").split(";"):
            if disease_id:
                genes_by_disease[disease_id].add(gene)
    return {disease_id: len(genes) for disease_id, genes in sorted(genes_by_disease.items(), key=lambda item: -len(item[1])) if len(genes) >= max_genes}


def disease_cluster_concentration(gene_observations: pd.DataFrame, top_clusters: int = 3) -> dict:
    """Per symptom: share of positive genes contributed by the largest disease clusters (assumption A7 audit).

    A symptom whose positives come mostly from one cluster (one disease or one group of diseases sharing
    genes) is learnable as disease identity rather than as mechanism, and the disease-cluster split is
    what keeps that from inflating the scores.
    """
    concentration: dict[str, dict] = {}
    for symptom, rows in gene_observations.groupby("symptom"):
        counts = rows.groupby("disease_cluster_id").perturbation_id.nunique().sort_values(ascending=False)
        total = int(counts.sum())
        top = counts.head(top_clusters)
        concentration[symptom] = {
            "positive_genes": total, "clusters": int(len(counts)),
            "largest_cluster_share": float(counts.iloc[0] / total) if total else 0.0,
            f"top_{top_clusters}_share": float(top.sum() / total) if total else 0.0,
            "top_clusters": {cluster: int(n) for cluster, n in top.items()},
        }
    return concentration


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--hpo-obo", type=Path, default=Path("data/raw/hpo/hp.obo"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--phenotype-hpoa", type=Path, default=Path("data/raw/hpo/phenotype.hpoa"), help="disease-level annotations with biocuration dates (time split)")
    parser.add_argument("--genes-to-disease", type=Path, default=Path("data/raw/hpo/genes_to_disease.txt"), help="OMIM association types (MENDELIAN, POLYGENIC); grade A needs a causal association when types are known")
    parser.add_argument("--orphadata-product6", type=Path, default=Path("data/raw/orphadata/en_product6.xml"), help="Orphadata product 6 gene-disease association types")
    parser.add_argument("--reference-publication-dates", type=Path, default=Path("docs/hpo_reference_publication_dates.json"), help="PMID -> publication date lookup; dates pairs by publication rather than curation")
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--max-drug-targets", type=int, default=1, help="most mechanism targets a qualifying SIDER drug may have; 0 lifts the single-target rule (pass the same value to build_onsides_reports.py)")
    parser.add_argument("--grade-a-policy", choices=GRADE_A_POLICIES, default=DEFAULT_GRADE_A_POLICY)
    parser.add_argument("--disease-cluster-max-genes", type=int, default=DEFAULT_DISEASE_CLUSTER_MAX_GENES,
                        help="disease entries annotated to this many genes or more (group-level Orphanet entries) do not link disease clusters; 0 disables the cap")
    parser.add_argument("--weighting", choices=WEIGHTINGS, default=DEFAULT_WEIGHTING,
                        help="grade: the grade table times the frequency clip (default); reliability: reliability_global_scale x the report-level Dawid-Skene posterior")
    parser.add_argument("--report-rubric-weights", type=Path, default=None, help="JSON overriding any field of RubricWeightDefaults (evidence-code factors, sample-size reference, minimum weight, implicit-negative weight)")
    parser.add_argument("--reliability-global-scale", type=float, default=1.0, help="multiplier on the posterior under --weighting reliability")
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--extra-reports", type=Path, nargs="*", default=[], help="report tables produced outside the assembler (OnSIDES label statements, literature reports) appended to the report table and the reliability fit")
    parser.add_argument("--onsides-bridge", type=Path, default=None, help="OnSIDES ingredient identifier bridge; its ChEMBL parents are the fallback for SIDER CIDs UniChem leaves without a mechanism")
    parser.add_argument("--non-protein-targets", type=Path, default=NON_PROTEIN_TARGETS_PATH, help="ChEMBL non-protein targets -> Human-GEM metabolites (docs/drug_targets_any_type.md)")
    parser.add_argument("--drugs-acting-as-graph-compounds", type=Path, default=DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH, help="drugs with no mechanism target that are themselves a graph metabolite")
    parser.add_argument("--seed-gain-of-function-genes-positive", action="store_true",
                        help="seed a gene with sign +1 when every causal association behind its reports is Orphanet's gain-of-function type; without it every gene has sign -1")
    parser.add_argument("--allele-perturbations", action="store_true",
                        help="decision 9.54: each gene's mechanism class is its own perturbation (loss of function or untyped: the gene, -1; gain of function: "
                             "<gene>|gain_of_function, +1, or the gene's curated neomorphic allele); replaces --seed-gain-of-function-genes-positive")
    parser.add_argument("--orphadata-product1", type=Path, default=Path("data/raw/orphadata/en_product1.xml"),
                        help="Orphadata product 1: how each Orphanet disorder maps to OMIM, which carries a typed mechanism to the same disease's OMIM record")
    parser.add_argument("--neomorphic-alleles", type=Path, default=DEFAULT_NEOMORPHIC_ALLELES_PATH,
                        help="curated neomorphic alleles, one row per cited mechanism edge (--allele-perturbations)")
    parser.add_argument("--mechanism-adjudications", default=str(DEFAULT_MECHANISM_ADJUDICATIONS_PATH),
                        help="cited adjudications of the pairs the two sources of a mechanism dispute (--allele-perturbations); none applies none")
    parser.add_argument("--gene2phenotype", default=str(DEFAULT_GENE2PHENOTYPE_PATH),
                        help="Gene2Phenotype's all-panel download, the second source of a disease's mechanism where Orphanet types none (--allele-perturbations); none builds without it")
    parser.add_argument("--admitted-atc-prefixes", type=Path, default=None,
                        help="written list of the ATC groups and substances admitted in place of group N alone (atc_prefix column); pass the same file to build_onsides_reports.py")
    parser.add_argument("--narrow-complex-targets", action=argparse.BooleanOptionalAction, default=True,
                        help="decision 8.11 (the owner, 13 October 2026): a drug's ChEMBL complex target is narrowed to the subunits the Guide to PHARMACOLOGY "
                             "names for it, with the binding site's obligate subunits restored; GtoPdb targets outside every ChEMBL row are not added. "
                             "Pass the same choice to build_onsides_reports.py")
    parser.add_argument("--gtopdb-dir", type=Path, default=GTOPDB_DIRECTORY, help="the Guide to PHARMACOLOGY release the narrowing reads")
    parser.add_argument("--binding-sites", type=Path, default=DEFAULT_BINDING_SITES_PATH, help="the drugs' binding sites, which name the obligate subunits")
    parser.add_argument("--open-targets-scoped-rows", type=Path, default=None,
                        help="data/processed/off_target_scoping/open_targets_symptom_rows.parquet (amendment item 12d, decision 8.5): the curated target-safety statements enter the Dawid-Skene reliability fit and reach no grade and no label. Omit the flag to assemble without the class")
    arguments = parser.parse_args()
    rubric_weight_defaults = load_rubric_weight_defaults(arguments.report_rubric_weights)
    assembled = assemble(arguments.crosswalk, arguments.hpo_obo, arguments.hpo_annotations, arguments.graph_dir, arguments.sider_dir, arguments.chembl_dir,
                         arguments.max_drug_targets or None, arguments.grade_a_policy, arguments.phenotype_hpoa, arguments.reference_publication_dates,
                         arguments.genes_to_disease, arguments.orphadata_product6, arguments.disease_cluster_max_genes or None, arguments.weighting, rubric_weight_defaults, arguments.reliability_global_scale,
                         extra_report_paths=list(arguments.extra_reports), onsides_bridge_path=arguments.onsides_bridge,
                         non_protein_targets_path=arguments.non_protein_targets, drugs_acting_as_graph_compounds_path=arguments.drugs_acting_as_graph_compounds,
                         seed_gain_of_function_genes_positive=arguments.seed_gain_of_function_genes_positive,
                         admitted_atc_prefixes=read_admitted_atc_prefixes(arguments.admitted_atc_prefixes) if arguments.admitted_atc_prefixes is not None else NERVOUS_SYSTEM_ATC_PREFIXES,
                         open_targets_scoped_rows_path=arguments.open_targets_scoped_rows,
                         allele_perturbations=arguments.allele_perturbations, orphadata_product1_path=arguments.orphadata_product1,
                         neomorphic_alleles_path=arguments.neomorphic_alleles,
                         gene2phenotype_path=None if str(arguments.gene2phenotype).lower() == "none" else Path(arguments.gene2phenotype),
                         mechanism_adjudications_path=None if str(arguments.mechanism_adjudications).lower() == "none" else Path(arguments.mechanism_adjudications),
                         complex_target_narrowing=ComplexTargetNarrowing.from_gtopdb(arguments.gtopdb_dir, arguments.binding_sites) if arguments.narrow_complex_targets else None)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in (("evidence_records.parquet", assembled.observations), ("unmapped_records.parquet", assembled.unmapped), ("evidence_reports.parquet", assembled.reports)):  # atomic replace for concurrent readers
        table.to_parquet(arguments.output_dir / (name + ".tmp"), index=False)
        os.replace(arguments.output_dir / (name + ".tmp"), arguments.output_dir / name)
    summary = summarize_observations(assembled.observations, assembled.unmapped)
    summary["grade_a_policy"] = arguments.grade_a_policy
    summary["weighting"] = arguments.weighting
    summary["reliability_global_scale"] = arguments.reliability_global_scale
    if arguments.seed_gain_of_function_genes_positive:  # recorded only when set, so a table built without it has the summary it always had
        summary["genes_seeded_as_gain_of_function"] = list(assembled.genes_seeded_as_gain_of_function)
    if assembled.allele_summary is not None:  # decision 9.54, recorded only when the alleles are perturbations
        summary["alleles"] = {**assembled.allele_summary, "orphadata_product1": str(arguments.orphadata_product1), "neomorphic_alleles": str(arguments.neomorphic_alleles)}
    summary["complex_target_narrowing"] = assembled.complex_target_narrowing  # decision 8.11; None when built with --no-narrow-complex-targets
    if arguments.admitted_atc_prefixes is not None:
        summary["admitted_atc_prefixes"] = {"file": str(arguments.admitted_atc_prefixes), "prefixes": list(read_admitted_atc_prefixes(arguments.admitted_atc_prefixes))}
    summary.update(summarize_reports(assembled.reports, assembled.reliability_fit, rubric_weight_defaults, assembled.observations, assembled.reports_dropped_unknown_provenance))
    (arguments.output_dir / "evidence_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
