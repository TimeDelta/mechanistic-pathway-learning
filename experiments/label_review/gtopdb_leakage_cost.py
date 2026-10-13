"""What the Guide to PHARMACOLOGY re-derivation costs in leakage, before any rebuild.

The owner, 13 October 2026: "go ahead and remeasure leakage." Changing a drug's targets changes which drugs share one,
and drugs sharing a target are merged into one leakage group (`merge_drugs_with_their_targets`), so a target change is
a leakage change. The target-limit change of 12 October was measured this way before it was accepted and this is the
same reading for the same reason.

Three arms, so the effect of the rule is separable from the effect of its threshold:

  - **as built**, the ChEMBL targets the graph carries now;
  - **re-derived**, the row-wise rule of `perturbation/gtopdb_drug_targets.py`: a single-gene ChEMBL row kept, a
    complex row narrowed only where GtoPdb names two or more of its subunits, a GtoPdb target outside every ChEMBL
    row added;
  - **re-derived, narrowing on one subunit**, the same with the threshold dropped to one, which is the rule the owner
    rejected over cocaine and is kept here to show what it would have cost.

Measurement only: no graph, evidence table or configuration is written.

Usage:
  OMP_NUM_THREADS=4 python experiments/label_review/gtopdb_leakage_cost.py \
      --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes \
      --evidence-dir data/processed/evidence_full_v4 \
      --label-selection data/processed/label_selection/better_v2_full_v4.parquet \
      --markdown-output docs/gtopdb_leakage_cost.md
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.label_review.rebuild_leakage_cost import kept_positive_counts, scheme  # noqa: E402
from experiments.measure_gtopdb_target_rederivation import chembl_rows_per_drug, drug_names  # noqa: E402
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data  # noqa: E402
from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import (  # noqa: E402
    ligand_id_of_drug,
    one_gene_per_row,
    pubchem_cid_from_perturbation_id,
    read_interactions,
    replace_complex_rows,
    targets_per_drug,
)


def rederived_seeds(data, graph_directory: Path, gtopdb_directory: Path, chembl_directory: Path,
                    subunits_needed: int, narrowing_only: bool = False) -> tuple[dict[str, list[int]], dict[str, int]]:
    """perturbation id -> the node indices its re-derived targets map to, and how many drugs changed.

    A gene the graph has no protein node for is dropped, as it is today; a drug with no GtoPdb row, or none the rule
    acts on, keeps the seeds the graph gives it."""
    names = drug_names(graph_directory)
    entry = pd.read_parquet(graph_directory / "drug_entry_nodes.parquet")
    drugs = sorted(set(entry.perturbation_id.astype(str)))
    cids = {drug: cid for drug in drugs if (cid := pubchem_cid_from_perturbation_id(drug)) is not None}
    ligand_of_drug, _ = ligand_id_of_drug(cids, {drug: names.get(drug, "") for drug in drugs}, gtopdb_directory)
    per_drug = targets_per_drug(ligand_of_drug, read_interactions(gtopdb_directory))
    chembl_rows = chembl_rows_per_drug(drugs, chembl_directory)
    node_of_gene = {str(node_id).split(":", 1)[1]: index for index, node_id in enumerate(data.node_ids)
                    if str(node_id).startswith("PROTEIN:")}
    seeds_of_drug, changed = {}, 0
    for drug in drugs:
        rows = per_drug.get(drug)
        if rows is None or drug not in chembl_rows:
            continue
        usable = one_gene_per_row(rows)
        if usable.empty:
            continue
        signs, reasons = replace_complex_rows(chembl_rows[drug], usable, subunits_needed=subunits_needed)
        if narrowing_only:  # keep GtoPdb's narrowing of a complex, drop its targets outside every ChEMBL row
            signs = {gene: sign for gene, sign in signs.items() if reasons[gene] != "gtopdb, outside every chembl row"}
        nodes = sorted({node_of_gene[gene] for gene in signs if gene in node_of_gene})
        if nodes:
            seeds_of_drug[drug] = nodes
            changed += 1
    return seeds_of_drug, changed


def reading_with_seeds(data, keep, kept_positive, seeds_of_drug: dict[str, list[int]]) -> dict:
    """The leakage reading with each named drug's seeds replaced; everything else as the table gives it."""
    original = list(data.perturbation_seeds)
    try:
        data.perturbation_seeds = [seeds_of_drug.get(perturbation_id, seeds)
                                   for perturbation_id, seeds in zip(data.perturbation_ids, original)]
        reading = scheme(data, keep, kept_positive)
    finally:
        data.perturbation_seeds = original
    return reading


def summarise(reading: dict) -> dict:
    return {key: reading[key] for key in ("drugs", "kept_positive_pairs", "drug_holding_groups",
                                          "two_largest_blocks_by_drugs", "largest_block_share_of_drugs",
                                          "largest_block_share_of_kept_pairs")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--label-selection", type=Path, required=True)
    parser.add_argument("--gtopdb-dir", type=Path, default=Path("data/raw/gtopdb_2026"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster_and_targets",
                                label_selection=arguments.label_selection)
    # kept_positive_counts reads this index, which rebuild_leakage_cost.py's own main() sets on the data it loads
    data.perturbation_index = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    kept_positive = kept_positive_counts(data, arguments.label_selection)
    keep = np.ones(len(data.perturbation_ids), dtype=bool)

    record = {"graph_dir": str(arguments.graph_dir), "evidence_dir": str(arguments.evidence_dir), "arms": {}}
    record["arms"]["as built (ChEMBL)"] = summarise(scheme(data, keep, kept_positive))
    for label, subunits_needed, narrowing_only in (
            ("re-derived, narrowing only (no GtoPdb target outside a ChEMBL row)", 2, True),
            ("re-derived, narrowing and additions", 2, False),
            ("re-derived, narrowing on one subunit and additions (the rejected rule)", 1, False)):
        seeds_of_drug, changed = rederived_seeds(data, arguments.graph_dir, arguments.gtopdb_dir,
                                                 arguments.chembl_dir, subunits_needed, narrowing_only)
        reading = summarise(reading_with_seeds(data, keep, kept_positive, seeds_of_drug))
        reading["drugs_whose_seeds_changed"] = changed
        record["arms"][label] = reading

    print(json.dumps(record, indent=2))
    if arguments.json_output:
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(json.dumps(record, indent=2) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.markdown_output.write_text(write_markdown(record))


def write_markdown(record: dict) -> str:
    lines = ["# What the GtoPdb re-derivation costs in leakage", "",
             "Measured by `experiments/label_review/gtopdb_leakage_cost.py` on "
             f"`{record['graph_dir']}` with `{record['evidence_dir']}`. Measurement only; nothing is written by it.",
             "",
             "| arm | drugs | kept positive pairs | drug-holding groups | two largest blocks | largest block, share of drugs | share of kept pairs |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for label, reading in record["arms"].items():
        lines.append(f"| {label} | {reading['drugs']} | {reading['kept_positive_pairs']} | "
                     f"{reading['drug_holding_groups']} | {reading['two_largest_blocks_by_drugs']} | "
                     f"{reading['largest_block_share_of_drugs']} | {reading['largest_block_share_of_kept_pairs']} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
