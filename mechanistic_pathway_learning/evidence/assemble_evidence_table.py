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
  grade, weight          from assign_evidence_grades (fixed-grade fallback; the learned
                         reliability model replaces weight downstream)
  perturbation_nodes     JSON list of [node_id, sign, magnitude] in the graph
  label_frequency        SIDER frequency midpoint (E2) or largest HPO frequency midpoint (E1) when reported
  omim_entry_count, orpha_entry_count, annotation_row_count, annotation_patient_count, disease_identifiers
                         provenance of a monogenic record (null for drugs)
  source                 provenance string

Rows whose perturbation has no node in the graph are written to a separate
unmapped table so the gap is visible rather than silently dropped.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import (
    DEFAULT_GRADE_A_POLICY,
    GRADE_A_POLICIES,
    EvidenceRecord,
    assign_evidence_grade,
    loss_weight_for_record,
)
from mechanistic_pathway_learning.evidence.load_drug_label_events import is_nervous_system_atc, load_sider_events
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    load_hpo_is_a_parents_from_obo,
    monogenic_evidence_records,
    parse_genes_to_phenotype,
    read_crosswalk_hpo_terms,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    drug_targets_for_pubchem_cid,
    has_dominant_target,
    load_chembl_caches,
)


def gene_node_lookup(nodes: pd.DataFrame) -> dict[str, str]:
    genes = nodes[nodes.node_type == "gene"]
    return {symbol: node_id for symbol, node_id in zip(genes.gene_symbol, genes.node_id) if isinstance(symbol, str) and symbol}


def disease_cluster_ids(records: list[EvidenceRecord]) -> dict[str, str]:
    """Gene -> cluster id, where genes sharing any disease identifier are in one cluster (union-find).

    The cluster is named after its alphabetically first gene. Genes with no shared disease form singletons.
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
    for genes in genes_by_disease.values():
        ordered = sorted(genes)
        for gene in ordered[1:]:
            union(ordered[0], gene)
    return {gene: "cluster:" + find(gene) for gene in list(parent)}


def assemble(
    crosswalk_path: Path,
    hpo_obo_path: Path,
    hpo_annotations_path: Path,
    graph_directory: Path,
    sider_directory: Path | None,
    chembl_directory: Path | None,
    max_drug_targets: int = 1,
    grade_a_policy: str = DEFAULT_GRADE_A_POLICY,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    node_by_symbol = gene_node_lookup(nodes)
    metabolic_symbols = set(nodes[(nodes.node_type == "gene") & (nodes.get("in_metabolic_layer", False) == True)].gene_symbol.dropna())  # noqa: E712
    crosswalk_terms = read_crosswalk_hpo_terms(crosswalk_path)
    symptom_to_hpo = {symptom: roots for symptom, (roots, _) in crosswalk_terms.items()}
    symptom_to_excluded = {symptom: excluded for symptom, (_, excluded) in crosswalk_terms.items()}
    rows: list[dict] = []
    unmapped: list[dict] = []

    parents = load_hpo_is_a_parents_from_obo(hpo_obo_path)
    monogenic_records = monogenic_evidence_records(parse_genes_to_phenotype(hpo_annotations_path), symptom_to_hpo, parents, set(node_by_symbol), symptom_to_excluded)
    cluster_by_gene = disease_cluster_ids(monogenic_records)
    for record in monogenic_records:
        node_id = node_by_symbol.get(record.perturbation_identifier)
        base = {
            "perturbation_id": record.perturbation_identifier, "perturbation_type": "gene", "perturbation_label": record.perturbation_identifier,
            "group_id": record.perturbation_identifier, "disease_cluster_id": cluster_by_gene[record.perturbation_identifier],
            "symptom": record.symptom_identifier, "relation": record.relation,
            "evidence_class": "monogenic", "grade": assign_evidence_grade(record, grade_a_policy), "weight": loss_weight_for_record(record, grade_a_policy=grade_a_policy),
            "label_frequency": record.max_annotation_frequency, "source": record.source, "in_metabolic_layer": record.perturbation_identifier in metabolic_symbols,
            "omim_entry_count": record.omim_entry_count, "orpha_entry_count": record.orpha_entry_count, "annotation_row_count": record.annotation_row_count,
            "annotation_patient_count": record.annotation_patient_count, "disease_identifiers": ";".join(record.disease_identifiers),
        }
        if node_id is None:
            unmapped.append(base)
            continue
        base["perturbation_nodes"] = json.dumps([[node_id, -1.0, 1.0]])
        rows.append(base)

    if sider_directory is not None and chembl_directory is not None and (chembl_directory / "targets.json").exists():
        caches = load_chembl_caches(chembl_directory)
        for event in load_sider_events(sider_directory, crosswalk_path):
            drug_targets = drug_targets_for_pubchem_cid(event.pubchem_cid, caches)
            if not (has_dominant_target(drug_targets, max_drug_targets) and is_nervous_system_atc(event.atc_codes)):
                continue
            record = EvidenceRecord(event.stitch_flat_id, event.target_symptom, event.relation, "pharmacological", source=event.source,
                                    cns_penetrant=True, has_dominant_target=True, label_event_frequency=event.label_frequency)
            perturbation_nodes = []
            for target in drug_targets:
                mapped = [node_by_symbol[symbol] for symbol in target.gene_symbols if symbol in node_by_symbol]
                perturbation_nodes.extend([node, target.sign, 1.0 / len(mapped)] for node in mapped)
            group_id = "|".join(sorted(target.target_chembl_id for target in drug_targets))
            base = {
                "perturbation_id": event.stitch_flat_id, "perturbation_type": "drug", "perturbation_label": event.drug_name,
                "group_id": group_id, "disease_cluster_id": group_id, "symptom": event.target_symptom, "relation": event.relation,
                "evidence_class": "pharmacological", "grade": assign_evidence_grade(record), "weight": loss_weight_for_record(record),
                "label_frequency": event.label_frequency, "source": f"{event.source}; targets " + ";".join(f"{target.target_chembl_id}:{target.action_type}" for target in drug_targets),
                "in_metabolic_layer": any(symbol in metabolic_symbols for target in drug_targets for symbol in target.gene_symbols),
                "omim_entry_count": None, "orpha_entry_count": None, "annotation_row_count": None, "annotation_patient_count": None, "disease_identifiers": None,
            }
            if not perturbation_nodes:
                unmapped.append(base)
                continue
            base["perturbation_nodes"] = json.dumps(perturbation_nodes)
            rows.append(base)
    return pd.DataFrame(rows), pd.DataFrame(unmapped)


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
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--hpo-obo", type=Path, default=Path("data/raw/hpo/hp.obo"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--max-drug-targets", type=int, default=1)
    parser.add_argument("--grade-a-policy", choices=GRADE_A_POLICIES, default=DEFAULT_GRADE_A_POLICY)
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/evidence"))
    arguments = parser.parse_args()
    observations, unmapped = assemble(arguments.crosswalk, arguments.hpo_obo, arguments.hpo_annotations, arguments.graph_dir, arguments.sider_dir, arguments.chembl_dir,
                                      arguments.max_drug_targets, arguments.grade_a_policy)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    observations.to_parquet(arguments.output_dir / "evidence_records.parquet", index=False)
    unmapped.to_parquet(arguments.output_dir / "unmapped_records.parquet", index=False)
    summary = summarize_observations(observations, unmapped)
    summary["grade_a_policy"] = arguments.grade_a_policy
    (arguments.output_dir / "evidence_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
