"""Phase 1 go/no-go counts (design section 9).

Evidence class E1 (monogenic) is computed from the HPO gene-to-phenotype file and the
Human-GEM gene table. Evidence class E2 (pharmacological) is computed from SIDER 4.1
joined to ChEMBL mechanisms (experiments/fetch_chembl_drug_targets.py) when those caches
exist. Writes data/processed/phase1_counts.json, docs/phase1_counts.md and prints the table.

Usage:
  python experiments/run_phase1_counts.py \
      --crosswalk docs/symptom_crosswalk.csv \
      --hpo-obo data/raw/hpo/hp.obo --hpo-annotations data/raw/hpo/genes_to_phenotype.txt \
      --human-gem-genes data/raw/Human-GEM/model/genes.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

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

GRADE_A_OR_B_PERTURBATIONS_PER_SYMPTOM = 15
MINIMUM_SYMPTOMS_PASSING = 8


def read_symptom_crosswalk(crosswalk_path: Path) -> dict[str, list[str]]:
    target_symptom_to_hpo_ids: dict[str, list[str]] = {}
    with open(crosswalk_path, encoding="utf-8") as crosswalk_file:
        for row in csv.DictReader(crosswalk_file):
            hpo_ids = [hpo_id.strip() for hpo_id in (row.get("hpo_ids") or "").split(";") if hpo_id.strip()]
            target_symptom_to_hpo_ids[row["target_symptom"]] = hpo_ids
    return target_symptom_to_hpo_ids


def read_human_gem_gene_symbols(genes_table_path: Path) -> set[str]:
    gene_symbols: set[str] = set()
    with open(genes_table_path, encoding="utf-8") as genes_file:
        for row in csv.DictReader(genes_file, delimiter="\t", quotechar='"'):
            for symbol in (row.get("geneSymbols") or "").split(";"):
                if symbol.strip():
                    gene_symbols.add(symbol.strip())
    return gene_symbols


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--hpo-obo", type=Path, default=Path("data/raw/hpo/hp.obo"))
    parser.add_argument("--hpo-annotations", type=Path, default=Path("data/raw/hpo/genes_to_phenotype.txt"))
    parser.add_argument("--human-gem-genes", type=Path, default=Path("data/raw/Human-GEM/model/genes.tsv"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/phase1_counts.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase1_counts.md"))
    arguments = parser.parse_args()

    target_symptom_to_hpo_ids = read_symptom_crosswalk(arguments.crosswalk)
    parents_by_term = load_hpo_is_a_parents_from_obo(arguments.hpo_obo)
    annotation_rows = parse_genes_to_phenotype(arguments.hpo_annotations)
    human_gem_gene_symbols = read_human_gem_gene_symbols(arguments.human_gem_genes)
    records = monogenic_evidence_records(annotation_rows, target_symptom_to_hpo_ids, parents_by_term, human_gem_gene_symbols)

    genes_by_symptom: dict[str, set[str]] = defaultdict(set)
    multi_disease_genes_by_symptom: dict[str, set[str]] = defaultdict(set)
    for record in records:
        genes_by_symptom[record.symptom_identifier].add(record.perturbation_identifier)
        if record.independent_case_series_count >= 2:
            multi_disease_genes_by_symptom[record.symptom_identifier].add(record.perturbation_identifier)

    e2_available = (arguments.sider_dir / "meddra_all_se.tsv.gz").exists() and (arguments.chembl_dir / "targets.json").exists()
    induces_by_symptom: dict[str, set[str]] = defaultdict(set)
    induces_in_gem_by_symptom: dict[str, set[str]] = defaultdict(set)
    relieves_by_symptom: dict[str, set[str]] = defaultdict(set)
    if e2_available:
        caches = load_chembl_caches(arguments.chembl_dir)
        for event in load_sider_events(arguments.sider_dir, arguments.crosswalk):
            drug_targets = drug_targets_for_pubchem_cid(event.pubchem_cid, caches)
            if not (has_dominant_target(drug_targets) and is_nervous_system_atc(event.atc_codes)):
                continue
            target_genes = {symbol for target in drug_targets for symbol in target.gene_symbols}
            if event.relation == "induces":
                induces_by_symptom[event.target_symptom].add(event.stitch_flat_id)
                if target_genes & human_gem_gene_symbols:
                    induces_in_gem_by_symptom[event.target_symptom].add(event.stitch_flat_id)
            else:
                relieves_by_symptom[event.target_symptom].add(event.stitch_flat_id)

    per_symptom = {}
    header = f"{'symptom':30s} {'E1 genes':>8s} {'>=2 dis':>8s} {'E2 induce':>9s} {'E2 in GEM':>9s} {'E2 relieve':>10s} {'A+B':>5s}  status"
    print(header)
    markdown_lines = ["# Phase 1 counts (generated by experiments/run_phase1_counts.py)", "",
                      "E1: Human-GEM genes with an HPO annotation in the symptom or its descendants. E2: SIDER 4.1 drugs with a single ChEMBL mechanism target and an ATC N code (centrally acting proxy); E2 in GEM: those whose target gene is in Human-GEM. A+B: E1 genes plus E2 inducing drugs.", "",
                      "| symptom | E1 genes | E1 genes with 2+ diseases | E2 inducing drugs | E2 inducing, target in Human-GEM | E2 relieving drugs | grade A+B perturbations | status |", "|---|---|---|---|---|---|---|---|"]
    symptoms_passing = 0
    for symptom in target_symptom_to_hpo_ids:
        gene_count = len(genes_by_symptom.get(symptom, set()))
        multi_count = len(multi_disease_genes_by_symptom.get(symptom, set()))
        induce_count = len(induces_by_symptom.get(symptom, set()))
        induce_gem_count = len(induces_in_gem_by_symptom.get(symptom, set()))
        relieve_count = len(relieves_by_symptom.get(symptom, set()))
        grade_a_or_b = gene_count + induce_count
        passes = grade_a_or_b >= GRADE_A_OR_B_PERTURBATIONS_PER_SYMPTOM
        symptoms_passing += int(passes)
        status = "pass" if passes else "below threshold"
        per_symptom[symptom] = {"e1_genes": gene_count, "e1_genes_with_two_or_more_diseases": multi_count, "e2_inducing_drugs": induce_count, "e2_inducing_drugs_target_in_human_gem": induce_gem_count, "e2_relieving_drugs": relieve_count, "grade_a_or_b_perturbations": grade_a_or_b, "passes_threshold": passes}
        print(f"{symptom:30s} {gene_count:8d} {multi_count:8d} {induce_count:9d} {induce_gem_count:9d} {relieve_count:10d} {grade_a_or_b:5d}  {status}")
        markdown_lines.append(f"| {symptom} | {gene_count} | {multi_count} | {induce_count} | {induce_gem_count} | {relieve_count} | {grade_a_or_b} | {status} |")
    all_genes = {record.perturbation_identifier for record in records}
    all_drugs = set().union(*induces_by_symptom.values()) if induces_by_symptom else set()
    verdict = symptoms_passing >= MINIMUM_SYMPTOMS_PASSING
    summary = (f"distinct Human-GEM genes with any target-symptom annotation: {len(all_genes)}; distinct qualifying E2 drugs: {len(all_drugs)}; "
               f"symptoms passing: {symptoms_passing} (need {MINIMUM_SYMPTOMS_PASSING}); go/no-go: {'go' if verdict else 'no-go'}")
    print("\n" + summary)
    markdown_lines += ["", summary, "", "Sources: HPO release and Human-GEM version as pinned in docs/data_sources.md; SIDER 4.1 (labels to 2015); ChEMBL mechanisms via experiments/fetch_chembl_drug_targets.py." if e2_available else "E2 columns are zero because the SIDER or ChEMBL caches were not found."]
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps({
        "human_gem_genes_in_table": len(human_gem_gene_symbols),
        "hpo_annotation_rows": len(annotation_rows),
        "distinct_human_gem_genes_with_target_symptom_annotation": len(all_genes),
        "distinct_qualifying_e2_drugs": len(all_drugs),
        "e2_available": e2_available,
        "per_symptom": per_symptom,
        "symptoms_passing": symptoms_passing,
        "go": verdict,
        "threshold_per_symptom": GRADE_A_OR_B_PERTURBATIONS_PER_SYMPTOM,
        "minimum_symptoms_passing": MINIMUM_SYMPTOMS_PASSING,
    }, indent=1))
    arguments.markdown_output.write_text("\n".join(markdown_lines) + "\n")


if __name__ == "__main__":
    main()
