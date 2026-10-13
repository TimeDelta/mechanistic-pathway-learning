#!/usr/bin/env python
"""What the pending v4 rebuild costs in leakage, measured before it is run.

The owner, 12 October 2026: "Before the v4 rebuild, what does it cost in terms of leakage?"

The rebuild applies three registered changes to the drug set at once, and each can move the leakage blocks that the
held-out split depends on:

  1. the limit on a drug's mechanism targets is lifted (amendment item 1);
  2. orphenadrine is set aside by the stated bridging rule - a drug whose target set joins two otherwise separate
     leakage blocks, found by removing each drug in turn;
  3. metoprolol and propranolol enter under the central-action rule (amendment item 14).

The definitions are those of docs/label_source_review.md 26.1, so the numbers compare with the table there: leakage
groups under the registered chained rule (disease clusters with each drug joined to the genes it targets and to the
drugs sharing a target node), kept positive pairs from the matching label selection, and the leave-one-drug-out
bridging measurement. Change 3's drugs are read from the any-ATC table, where neither has more than three mechanism
targets, so the limit there does not bind and their rows are the rows the rebuild would give them.

Measurement only: no evidence table, selection, graph or model is changed, and no held-out set is drawn.

Usage:
  PYTHONPATH=. python experiments/label_review/rebuild_leakage_cost.py
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, merge_drugs_with_their_targets

_spec = importlib.util.spec_from_file_location("scoping", Path(__file__).resolve().parents[1] / "scope_off_target_binding.py")
scoping = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(scoping)

SET_ASIDE_BY_THE_BRIDGING_RULE = ("orphenadrine",)
ADMITTED_BY_THE_CENTRAL_ACTION_RULE = ("metoprolol", "propranolol")


def kept_positive_counts(data, selection_path: Path) -> np.ndarray:
    """Kept positive pairs per perturbation under a label selection, as 26.1 counted them."""
    outcomes = data.outcomes > 0
    if selection_path is None or not selection_path.exists():
        return outcomes.sum(axis=1)
    selection = pd.read_parquet(selection_path)
    kept = selection.keep.astype(bool)
    if "counts_as_negative" in selection.columns:
        kept = kept & ~selection.counts_as_negative.fillna(False).astype(bool)
    set_aside = set(zip(selection.perturbation_id[~kept], selection.symptom[~kept]))
    symptom_column = {symptom: index for index, symptom in enumerate(data.symptoms)}
    mask = np.ones_like(outcomes, dtype=bool)
    for perturbation_id, symptom in set_aside:
        if perturbation_id in data.perturbation_index and symptom in symptom_column:
            mask[data.perturbation_index[perturbation_id], symptom_column[symptom]] = False
    return (outcomes & mask).sum(axis=1)


def scheme(data, keep: np.ndarray, kept_positive: np.ndarray) -> dict:
    """The leakage reading of the rows where keep is True, under the registered chained grouping."""
    positions = np.flatnonzero(keep)
    ids = [data.perturbation_ids[p] for p in positions]
    kinds = [data.perturbation_types[p] for p in positions]
    clusters = [data.group_ids[p] for p in positions]
    seeds = [data.perturbation_seeds[p] for p in positions]
    groups = merge_drugs_with_their_targets(ids, kinds, clusters, seeds, data.node_ids)
    reading = scoping.largest_group(groups, kept_positive[positions], kinds)
    by_drugs = Counter()
    for group, kind in zip(groups, kinds):
        by_drugs[group] += kind == "drug"
    two_largest = [count for _, count in by_drugs.most_common(2)]
    kept_in_group = Counter()
    for group, count in zip(groups, kept_positive[positions]):
        kept_in_group[group] += int(count)
    largest_group_name = by_drugs.most_common(1)[0][0]
    return {"drugs": int(sum(kind == "drug" for kind in kinds)),
            "kept_positive_pairs": int(kept_positive[positions].sum()),
            "drug_holding_groups": reading["drug_holding_groups"],
            "two_largest_blocks_by_drugs": two_largest,
            "largest_block_share_of_drugs": reading["largest_by_drugs"]["share_of_drugs"],
            "largest_block_share_of_kept_pairs": round(kept_in_group[largest_group_name] / max(int(kept_positive[positions].sum()), 1), 3),
            "groups": groups, "ids": ids, "kinds": kinds}


def bridging(data, keep: np.ndarray, kept_positive: np.ndarray, label_of: dict[str, str], how_many: int = 8) -> list[dict]:
    """The drugs whose removal shrinks the largest block most: the stated rule for finding a bridge."""
    base = scheme(data, keep, kept_positive)
    whole = max(Counter(group for group, kind in zip(base["groups"], base["kinds"]) if kind == "drug").values())
    effects = []
    for position in np.flatnonzero(keep):
        if data.perturbation_types[position] != "drug":
            continue
        trial = keep.copy()
        trial[position] = False
        counts = Counter(group for group, kind in zip(*[scheme(data, trial, kept_positive)[key] for key in ("groups", "kinds")]) if kind == "drug")
        effects.append({"drug": label_of.get(data.perturbation_ids[position], data.perturbation_ids[position]),
                        "drugs_removed_from_the_largest_block": whole - max(counts.values())})
    return sorted(effects, key=lambda entry: -entry["drugs_removed_from_the_largest_block"])[:how_many]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--no-limit-dir", type=Path, default=Path("data/processed/label_review/evidence_approved_no_cap"))
    parser.add_argument("--no-limit-selection", type=Path, default=Path("data/processed/label_review/selection_approved_no_cap.parquet"))
    parser.add_argument("--three-target-dir", type=Path, default=Path("data/processed/label_review/evidence_approved"))
    parser.add_argument("--three-target-selection", type=Path, default=Path("data/processed/label_review/selection_approved.parquet"))
    parser.add_argument("--any-atc-dir", type=Path, default=Path("data/processed/label_review/evidence_any_atc_cap3"))
    parser.add_argument("--any-atc-selection", type=Path, default=Path("data/processed/label_review/selection_any_atc_cap3.parquet"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/rebuild_leakage_cost.json"))
    parser.add_argument("--document", type=Path, default=Path("docs/rebuild_leakage_cost.md"))
    arguments = parser.parse_args()

    # The rebuild's own table: the no-limit rows, plus the two drugs the central-action rule admits, read from the
    # any-ATC table where the limit does not bind for either of them.
    no_limit = pd.read_parquet(arguments.no_limit_dir / "evidence_records.parquet")
    any_atc = pd.read_parquet(arguments.any_atc_dir / "evidence_records.parquet")
    added = any_atc[any_atc.perturbation_label.str.lower().isin(ADMITTED_BY_THE_CENTRAL_ACTION_RULE)]
    too_many_targets = [row.perturbation_label for row in added.drop_duplicates("perturbation_id").itertuples(index=False)
                        if len(json.loads(row.perturbation_nodes)) > 3]
    if too_many_targets:
        raise ValueError(f"{too_many_targets} have more than three mechanism targets in {arguments.any_atc_dir}, so the limit there binds and their rows are not the rebuild's")
    combined_directory = arguments.output.parent / "evidence_rebuild_candidate"
    combined_directory.mkdir(parents=True, exist_ok=True)
    pd.concat([no_limit, added[no_limit.columns]], ignore_index=True).to_parquet(combined_directory / "evidence_records.parquet", index=False)
    selection = pd.concat([pd.read_parquet(arguments.no_limit_selection),
                           pd.read_parquet(arguments.any_atc_selection).pipe(lambda table: table[table.perturbation_id.isin(added.perturbation_id)])],
                          ignore_index=True)
    combined_selection = arguments.output.parent / "selection_rebuild_candidate.parquet"
    selection.to_parquet(combined_selection, index=False)

    readings, bridges = {}, {}
    for name, evidence_directory, selection_path, drop in (
            ("three targets, the table in use", arguments.three_target_dir, arguments.three_target_selection, ()),
            ("no limit, every drug", arguments.no_limit_dir, arguments.no_limit_selection, ()),
            ("no limit, orphenadrine set aside", arguments.no_limit_dir, arguments.no_limit_selection, SET_ASIDE_BY_THE_BRIDGING_RULE),
            ("the rebuild: no limit, orphenadrine out, the two beta blockers in", combined_directory, combined_selection, SET_ASIDE_BY_THE_BRIDGING_RULE)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            data = load_experiment_data(arguments.graph_dir, evidence_directory, group_by="disease_cluster")
        data.perturbation_index = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
        label_of = dict(zip(data.perturbation_ids, data.perturbation_labels))
        kept_positive = kept_positive_counts(data, selection_path)
        keep = np.array([label_of[p].lower() not in {name.lower() for name in drop} for p in data.perturbation_ids])
        reading = scheme(data, keep, kept_positive)
        readings[name] = {key: value for key, value in reading.items() if key not in ("groups", "ids", "kinds")}
        print(f"\n{name}: {readings[name]}")
        if name.startswith("the rebuild"):
            bridges[name] = bridging(data, keep, kept_positive, label_of)
            print("  removing one drug at a time, the drugs whose removal shrinks the largest block most:")
            for entry in bridges[name]:
                print(f"    {entry['drug']}: {entry['drugs_removed_from_the_largest_block']}")

    summary = {"readings": readings, "bridging_after_the_rebuild": bridges,
               "set_aside_by_the_bridging_rule": list(SET_ASIDE_BY_THE_BRIDGING_RULE),
               "admitted_by_the_central_action_rule": list(ADMITTED_BY_THE_CENTRAL_ACTION_RULE)}
    arguments.output.write_text(json.dumps(summary, indent=2) + "\n")
    write_document(arguments.document, summary)


def write_document(path: Path, summary: dict) -> None:
    provenance = ("Measured by `experiments/label_review/rebuild_leakage_cost.py` before the rebuild is run, with the "
                  "definitions of `docs/label_source_review.md` 26.1 so the numbers compare with the table there. "
                  "Nothing is rebuilt, no held-out set is drawn.")
    lines = ["# What the v4 rebuild costs in leakage", "", provenance, "",
             "| scheme | drugs | kept positive pairs | drug-holding groups | two largest blocks, by drugs | largest block's share of drugs / of kept pairs |",
             "| --- | --- | --- | --- | --- | --- |"]
    for name, reading in summary["readings"].items():
        blocks = ", ".join(str(count) for count in reading["two_largest_blocks_by_drugs"])
        lines.append(f"| {name} | {reading['drugs']} | {reading['kept_positive_pairs']} | {reading['drug_holding_groups']} | "
                     f"{blocks} | {reading['largest_block_share_of_drugs']:.0%} / {reading['largest_block_share_of_kept_pairs']:.0%} |")
    for name, entries in summary["bridging_after_the_rebuild"].items():
        how = ("Each row is the number of drugs the largest block loses when that drug alone is removed, which is the "
               "stated rule for finding a bridge.")
        lines += ["", f"## Bridging drugs after the rebuild ({name})", "", how, "",
                  "| drug | drugs removed from the largest block |", "| --- | --- |"]
        lines += [f"| {entry['drug']} | {entry['drugs_removed_from_the_largest_block']} |" for entry in entries]
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
