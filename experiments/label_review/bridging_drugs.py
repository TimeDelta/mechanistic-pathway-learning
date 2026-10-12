"""Which added drugs join the two large drug blocks of a base table, and what the groups are without them.

Part of the label source review (docs/label_source_review.md, section 12). Counting only: no label, table or
configuration of the study is changed.

Why. The registered leakage rule chains: perturbations that share a node are joined into one group. On the three-target
table the drugs fall mainly into two groups, called blocks here (114 drugs around the monoamine receptors and
transporters, 71 around the GABA-A receptor and the ion channels). A wider drug rule can join the two blocks into one,
and one drug whose targets lie in both is enough to do it. This script names those drugs and recounts the groups with
them left out.

Usage:
  PYTHONPATH=. python experiments/label_review/bridging_drugs.py
  PYTHONPATH=. python experiments/label_review/bridging_drugs.py --base approved --wider approved_no_cap approved_with_list_no_cap
"""
from __future__ import annotations

import argparse
import shutil
import sys
import tempfile
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from target_limit_comparison import CONFIRMATORY_GRAPH, REVIEW_DIRECTORY, load_variant

from mechanistic_pathway_learning.evaluation.experiment_data import ExperimentData, load_experiment_data


def seed_nodes_of(data: ExperimentData, position: int) -> set[str]:
    return {data.node_ids[seed] for seed in data.perturbation_seeds[position]}


def nodes_of_the_two_largest_drug_blocks(base_data: ExperimentData) -> list[tuple[int, set[str]]]:
    """(number of drugs, seed nodes) of the two leakage groups of the base table that hold the most drugs."""
    group_ids = np.array(base_data.group_ids)
    is_drug = np.array(base_data.perturbation_types) == "drug"
    blocks = []
    for group_id, drug_count in Counter(group_ids[is_drug]).most_common(2):
        nodes: set[str] = set()
        for position in np.flatnonzero(group_ids == group_id):
            nodes |= seed_nodes_of(base_data, position)
        blocks.append((drug_count, nodes))
    return blocks


def bridging_drugs(base_variant: str, wider_variant: str) -> pd.DataFrame:
    """The drugs the wider table adds whose targets lie in both blocks of the base table."""
    base_data, base_records, _ = load_variant(base_variant)
    wider_data, wider_records, wider_kept_drug_pairs = load_variant(wider_variant)
    blocks = nodes_of_the_two_largest_drug_blocks(base_data)
    base_drug_ids = set(base_records[base_records.perturbation_type == "drug"].perturbation_id)
    label_by_id = dict(zip(wider_records.perturbation_id, wider_records.perturbation_label))
    kept_pairs_by_drug = Counter(wider_kept_drug_pairs.perturbation_id)
    node_table = pd.read_parquet(CONFIRMATORY_GRAPH / "nodes.parquet")
    name_of_node = dict(zip(node_table.node_id, node_table.display_name))
    rows = []
    for position, perturbation_id in enumerate(wider_data.perturbation_ids):
        if wider_data.perturbation_types[position] != "drug" or perturbation_id in base_drug_ids:
            continue
        seed_nodes = seed_nodes_of(wider_data, position)
        shared_with_each_block = [seed_nodes & block_nodes for _, block_nodes in blocks]
        if all(shared_with_each_block):
            rows.append({"perturbation_id": perturbation_id, "drug": label_by_id.get(perturbation_id, perturbation_id), "kept_pairs": kept_pairs_by_drug.get(perturbation_id, 0),
                         **{f"targets_in_the_block_of_{drug_count}_drugs": ", ".join(sorted(str(name_of_node.get(node, node)) for node in shared)[:8])
                            for (drug_count, _), shared in zip(blocks, shared_with_each_block)}})
    return pd.DataFrame(rows)


def groups_without(wider_variant: str, left_out_ids: set[str]) -> dict:
    """Drugs, kept drug pairs and the two largest drug groups of the wider table with some perturbations left out."""
    source_directory = REVIEW_DIRECTORY / f"evidence_{wider_variant}"
    selection = pd.read_parquet(REVIEW_DIRECTORY / f"selection_{wider_variant}.parquet")
    with tempfile.TemporaryDirectory() as temporary_name:
        temporary_directory = Path(temporary_name) / "evidence"
        shutil.copytree(source_directory, temporary_directory)
        for table_path in temporary_directory.glob("*.parquet"):
            table = pd.read_parquet(table_path)
            if "perturbation_id" in table.columns:
                table[~table.perturbation_id.isin(left_out_ids)].to_parquet(table_path, index=False)
        selection_path = Path(temporary_name) / "selection.parquet"
        selection[~selection.perturbation_id.isin(left_out_ids)].to_parquet(selection_path, index=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            data = load_experiment_data(CONFIRMATORY_GRAPH, temporary_directory, label_selection=selection_path, group_by="disease_cluster_and_targets")
    kept_drug_pairs = selection[(selection.perturbation_type == "drug") & selection.keep.astype(bool) & ~selection.perturbation_id.isin(left_out_ids)]
    kept_pairs_by_drug = Counter(kept_drug_pairs.perturbation_id)
    perturbation_ids = np.array(data.perturbation_ids)
    group_ids = np.array(data.group_ids)
    is_drug = np.array(data.perturbation_types) == "drug"
    summary = {"drugs": int(is_drug.sum()), "kept_drug_pairs": len(kept_drug_pairs), "groups_holding_a_drug": len(set(group_ids[is_drug]))}
    for rank, (group_id, drug_count) in enumerate(Counter(group_ids[is_drug]).most_common(2), start=1):
        member_drug_ids = perturbation_ids[is_drug & (group_ids == group_id)]
        kept_pairs_in_group = sum(kept_pairs_by_drug[drug_id] for drug_id in member_drug_ids)
        summary |= {f"group_{rank}_drugs": drug_count, f"group_{rank}_share_of_kept_drug_pairs": round(kept_pairs_in_group / max(len(kept_drug_pairs), 1), 3)}
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", default="approved", help="the counting build whose two largest drug blocks are the reference")
    parser.add_argument("--wider", nargs="+", default=["approved_no_cap", "approved_with_list", "approved_with_list_no_cap"])
    arguments = parser.parse_args()
    pd.set_option("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 90)
    for wider_variant in arguments.wider:
        if not (REVIEW_DIRECTORY / f"selection_{wider_variant}.parquet").exists():
            print(f"{wider_variant}: not built, skipped")
            continue
        bridging = bridging_drugs(arguments.base, wider_variant)
        print(f"\n== {wider_variant}: {len(bridging)} added drug(s) with targets in both blocks of {arguments.base}")
        if len(bridging):
            print(bridging.drop(columns="perturbation_id").to_string(index=False))
        print("   with every drug kept:      ", groups_without(wider_variant, set()))
        if len(bridging):
            print("   with those drugs left out: ", groups_without(wider_variant, set(bridging.perturbation_id)))


if __name__ == "__main__":
    main()
