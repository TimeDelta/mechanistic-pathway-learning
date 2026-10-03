"""Assemble the perturbation-symptom observation table (design sections 4.2 and 4.3).

Joins the monogenic records (HPO x Human-GEM) and the pharmacological records
(SIDER x ChEMBL) to the physiology graph and writes one row per observation:

  perturbation_id        gene symbol (E1) or SIDER STITCH flat id (E2)
  perturbation_type      "gene" or "drug"
  perturbation_label     readable name
  group_id               leakage group for splits: the gene itself, or the dominant ChEMBL target for a drug
  symptom                target symptom from docs/symptom_crosswalk.csv
  relation               "induces" or "relieves"
  evidence_class         "monogenic" or "pharmacological"
  grade, weight          from assign_evidence_grades (fixed-grade fallback; the learned
                         reliability model replaces weight downstream)
  perturbation_nodes     JSON list of [node_id, sign, magnitude] in the graph
  label_frequency        SIDER frequency midpoint when reported
  source                 provenance string

Rows whose perturbation has no node in the graph are written to a separate
unmapped table so the gap is visible rather than silently dropped.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord, assign_evidence_grade, loss_weight_for_record
from mechanistic_pathway_learning.evidence.load_drug_label_events import is_nervous_system_atc, load_sider_events
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    load_hpo_is_a_parents_from_obo,
    monogenic_evidence_records,
    parse_genes_to_phenotype,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (
    drug_targets_for_pubchem_cid,
    has_dominant_target,
    load_chembl_caches,
)


def read_crosswalk_hpo_ids(crosswalk_path: Path) -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = {}
    with open(crosswalk_path, encoding="utf-8") as crosswalk_file:
        for row in csv.DictReader(crosswalk_file):
            mapping[row["target_symptom"]] = [hpo_id.strip() for hpo_id in (row.get("hpo_ids") or "").split(";") if hpo_id.strip()]
    return mapping


def gene_node_lookup(nodes: pd.DataFrame) -> dict[str, str]:
    genes = nodes[nodes.node_type == "gene"]
    return {symbol: node_id for symbol, node_id in zip(genes.gene_symbol, genes.node_id) if isinstance(symbol, str) and symbol}


def assemble(
    crosswalk_path: Path,
    hpo_obo_path: Path,
    hpo_annotations_path: Path,
    graph_directory: Path,
    sider_directory: Path | None,
    chembl_directory: Path | None,
    max_drug_targets: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    node_by_symbol = gene_node_lookup(nodes)
    metabolic_symbols = set(nodes[(nodes.node_type == "gene") & (nodes.get("in_metabolic_layer", False) == True)].gene_symbol.dropna())  # noqa: E712
    symptom_to_hpo = read_crosswalk_hpo_ids(crosswalk_path)
    rows: list[dict] = []
    unmapped: list[dict] = []

    parents = load_hpo_is_a_parents_from_obo(hpo_obo_path)
    for record in monogenic_evidence_records(parse_genes_to_phenotype(hpo_annotations_path), symptom_to_hpo, parents, set(node_by_symbol)):
        node_id = node_by_symbol.get(record.perturbation_identifier)
        base = {
            "perturbation_id": record.perturbation_identifier, "perturbation_type": "gene", "perturbation_label": record.perturbation_identifier,
            "group_id": record.perturbation_identifier, "symptom": record.symptom_identifier, "relation": record.relation,
            "evidence_class": "monogenic", "grade": assign_evidence_grade(record), "weight": loss_weight_for_record(record),
            "label_frequency": None, "source": record.source, "in_metabolic_layer": record.perturbation_identifier in metabolic_symbols,
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
            base = {
                "perturbation_id": event.stitch_flat_id, "perturbation_type": "drug", "perturbation_label": event.drug_name,
                "group_id": "|".join(sorted(target.target_chembl_id for target in drug_targets)), "symptom": event.target_symptom, "relation": event.relation,
                "evidence_class": "pharmacological", "grade": assign_evidence_grade(record), "weight": loss_weight_for_record(record),
                "label_frequency": event.label_frequency, "source": f"{event.source}; targets " + ";".join(f"{target.target_chembl_id}:{target.action_type}" for target in drug_targets),
                "in_metabolic_layer": any(symbol in metabolic_symbols for target in drug_targets for symbol in target.gene_symbols),
            }
            if not perturbation_nodes:
                unmapped.append(base)
                continue
            base["perturbation_nodes"] = json.dumps(perturbation_nodes)
            rows.append(base)
    return pd.DataFrame(rows), pd.DataFrame(unmapped)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--hpo-obo", type=Path, default=Path("data/raw/hpo/hp.obo"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--max-drug-targets", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/evidence"))
    arguments = parser.parse_args()
    observations, unmapped = assemble(arguments.crosswalk, arguments.hpo_obo, arguments.hpo_annotations, arguments.graph_dir, arguments.sider_dir, arguments.chembl_dir, arguments.max_drug_targets)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    observations.to_parquet(arguments.output_dir / "evidence_records.parquet", index=False)
    unmapped.to_parquet(arguments.output_dir / "unmapped_records.parquet", index=False)
    summary = {
        "observations": int(len(observations)),
        "unmapped": int(len(unmapped)),
        "by_class_and_relation": {f"{cls}|{rel}": int(n) for (cls, rel), n in observations.groupby(["evidence_class", "relation"]).size().items()} if len(observations) else {},
        "distinct_perturbations": int(observations.perturbation_id.nunique()) if len(observations) else 0,
        "distinct_groups": int(observations.group_id.nunique()) if len(observations) else 0,
        "by_symptom": {symptom: int(n) for symptom, n in observations.groupby("symptom").size().items()} if len(observations) else {},
    }
    (arguments.output_dir / "evidence_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
