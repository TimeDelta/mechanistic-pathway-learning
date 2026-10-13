"""Splitting the held-out set on mechanism targets instead of on drugs.

The owner, 13 October 2026: "NOT splitting on drug itself but splitting on each target individually (so each drug can
have some labels in the test/validation and some in the held-out but all of the labels for a specific mechanism are
consistently split across test and held out) ... test this with full simulated redraws".

Splitting on drugs cannot work: measured over 200 random one-seventh draws, 84.5 percent of held-out drugs have every
target already present in some training drug, because the drugs genuinely share mechanisms (58 of the 192 targets are
held by six or more drugs). Grouping hard enough to remove that gives one block of 90 percent of the drugs. So the
question is whether the **target** is a better unit than the drug.

The sentence admits more than one reading and this measures all of them rather than picking one, because which is
meant changes the answer:

  - **all_targets_held**: a drug's labels are held out when every one of its targets is in the held-out target set. A
    promiscuous drug is then almost never held out.
  - **primary_target_held**: held out when the drug's best-bound target is, read from the occupancy table. This is the
    reading on which "each drug can have some labels in test and some in held-out" is true once labels are attributed
    per target.
  - **any_target_held**: held out when any of its targets is. A promiscuous drug is then almost always held out.
  - **mechanism_edge**: the unit is a (drug, target) mechanism edge rather than a drug. A drug keeps labels on both
    sides, and what is held out is the drug acting **through a particular target**.

For each, over full redraws, the reading that matters is the same one that condemned the drug split: how much of the
held-out side is already covered by the training side's mechanisms.

Measurement only; nothing is written to a graph, evidence table or configuration.

Usage:
  OMP_NUM_THREADS=4 python experiments/label_review/target_split_options.py \
      --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes \
      --evidence-dir data/processed/evidence_full_v4 --markdown-output docs/target_split_options.md
"""
from __future__ import annotations

import argparse
import json
import random
import statistics
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from experiments.measure_gtopdb_target_rederivation import chembl_rows_per_drug, drug_names  # noqa: E402
from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import (  # noqa: E402
    ligand_id_of_drug,
    one_gene_per_row,
    pubchem_cid_from_perturbation_id,
    read_interactions,
    replace_complex_rows,
    targets_per_drug,
)

SPLIT_RULES = ("all_targets_held", "primary_target_held", "any_target_held", "mechanism_edge")
DEFAULT_REDRAWS = 200
DEFAULT_HELD_OUT_FRACTION = 1 / 7


def rederived_targets(graph_directory: Path, gtopdb_directory: Path, chembl_directory: Path) -> dict[str, set[str]]:
    """Each drug's targets with the GtoPdb narrowing and every addition kept, which is the set the owner accepted."""
    entry = pd.read_parquet(graph_directory / "drug_entry_nodes.parquet")
    names = drug_names(graph_directory)
    drugs = sorted(set(entry.perturbation_id.astype(str)))
    cids = {drug: cid for drug in drugs if (cid := pubchem_cid_from_perturbation_id(drug)) is not None}
    ligand_of_drug, _ = ligand_id_of_drug(cids, {drug: names.get(drug, "") for drug in drugs}, gtopdb_directory)
    per_drug = targets_per_drug(ligand_of_drug, read_interactions(gtopdb_directory))
    chembl_rows = chembl_rows_per_drug(drugs, chembl_directory)
    targets: dict[str, set[str]] = {}
    for drug in drugs:
        base = {str(node).split(":", 1)[1] for node in entry[entry.perturbation_id == drug].target_node_id}
        if drug in per_drug and drug in chembl_rows:
            usable = one_gene_per_row(per_drug[drug])
            if not usable.empty:
                signs, _ = replace_complex_rows(chembl_rows[drug], usable)
                base = set(signs)
        targets[drug] = base
    return targets


def primary_target_of_drug(occupancy_path: Path) -> dict[str, str]:
    """The best-bound target per drug, by measured occupancy, for the primary_target_held rule."""
    if not occupancy_path.exists():
        return {}
    table = pd.read_csv(occupancy_path).dropna(subset=["occupancy_high"])
    best = table.sort_values("occupancy_high").groupby("perturbation_id").tail(1)
    return {str(row.perturbation_id): str(row.gene) for row in best.itertuples()}


def held_out_drugs(rule: str, targets: dict[str, set[str]], held_targets: set[str],
                   primary_of_drug: dict[str, str]) -> set[str]:
    if rule == "all_targets_held":
        return {drug for drug, own in targets.items() if own and own <= held_targets}
    if rule == "any_target_held":
        return {drug for drug, own in targets.items() if own & held_targets}
    if rule == "primary_target_held":
        return {drug for drug, own in targets.items() if primary_of_drug.get(drug) in held_targets}
    raise ValueError(rule)


def reading(rule: str, targets: dict[str, set[str]], primary_of_drug: dict[str, str],
            redraws: int, held_out_fraction: float, seed: int = 0) -> dict:
    """Over full redraws of the target partition: how large the held-out side is and how much of it the training side
    already covers mechanistically, which is the quantity that condemned the drug split."""
    generator = random.Random(seed)
    all_targets = sorted({target for own in targets.values() for target in own})
    drugs = sorted(targets)
    held_shares, covered_shares, unseen_shares, empty_draws = [], [], [], 0
    for _ in range(redraws):
        held_targets = set(generator.sample(all_targets, max(1, round(held_out_fraction * len(all_targets)))))
        if rule == "mechanism_edge":
            # the unit is a (drug, target) pair: every drug keeps labels, and what is held out is the drug acting
            # through a held-out target. "Covered" then asks whether that same target is still reachable in training
            # through another drug, which is the leak this rule is meant to close and closes by construction.
            held_pairs = [(drug, target) for drug, own in targets.items() for target in own & held_targets]
            training_pairs = [(drug, target) for drug, own in targets.items() for target in own - held_targets]
            held_shares.append(len(held_pairs) / max(1, len(held_pairs) + len(training_pairs)))
            training_targets = {target for _, target in training_pairs}
            covered_shares.append(sum(1 for _, target in held_pairs if target in training_targets) / max(1, len(held_pairs)))
            unseen_shares.append(1.0 - covered_shares[-1])
            continue
        held = held_out_drugs(rule, targets, held_targets, primary_of_drug)
        if not held:
            empty_draws += 1
            continue
        training_targets = {target for drug in drugs if drug not in held for target in targets[drug]}
        held_shares.append(len(held) / len(drugs))
        covered_shares.append(sum(1 for drug in held if targets[drug] and targets[drug] <= training_targets) / len(held))
        unseen_shares.append(sum(1 for drug in held if targets[drug] - training_targets) / len(held))
    return {"rule": rule, "draws_with_an_empty_held_out_side": empty_draws,
            "held_out_share": round(statistics.mean(held_shares), 4) if held_shares else None,
            "held_out_share_spread": round(statistics.pstdev(held_shares), 4) if len(held_shares) > 1 else None,
            "fully_covered_by_training_mechanisms": round(statistics.mean(covered_shares), 4) if covered_shares else None,
            "has_a_mechanism_training_never_saw": round(statistics.mean(unseen_shares), 4) if unseen_shares else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--gtopdb-dir", type=Path, default=Path("data/raw/gtopdb_2026"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--occupancy", type=Path, default=Path("data/processed/exposure/target_occupancy_v4.csv"))
    parser.add_argument("--redraws", type=int, default=DEFAULT_REDRAWS)
    parser.add_argument("--held-out-fraction", type=float, default=DEFAULT_HELD_OUT_FRACTION)
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()

    targets = rederived_targets(arguments.graph_dir, arguments.gtopdb_dir, arguments.chembl_dir)
    primary = primary_target_of_drug(arguments.occupancy)
    record = {"graph_dir": str(arguments.graph_dir), "drugs": len(targets),
              "distinct_targets": len({target for own in targets.values() for target in own}),
              "redraws": arguments.redraws, "held_out_fraction": round(arguments.held_out_fraction, 4),
              "drugs_with_a_measured_primary_target": sum(1 for drug in targets if drug in primary),
              "readings": [reading(rule, targets, primary, arguments.redraws, arguments.held_out_fraction)
                           for rule in SPLIT_RULES]}
    # the drug split, for comparison, drawn the same way
    generator = random.Random(0)
    drugs = sorted(targets)
    covered = []
    for _ in range(arguments.redraws):
        held = set(generator.sample(drugs, max(1, round(arguments.held_out_fraction * len(drugs)))))
        training_targets = {target for drug in drugs if drug not in held for target in targets[drug]}
        covered.append(sum(1 for drug in held if targets[drug] and targets[drug] <= training_targets) / len(held))
    record["readings"].insert(0, {"rule": "drug (as drawn today)", "draws_with_an_empty_held_out_side": 0,
                                  "held_out_share": round(arguments.held_out_fraction, 4), "held_out_share_spread": 0.0,
                                  "fully_covered_by_training_mechanisms": round(statistics.mean(covered), 4),
                                  "has_a_mechanism_training_never_saw": round(1 - statistics.mean(covered), 4)})

    print(json.dumps(record, indent=2))
    if arguments.json_output:
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(json.dumps(record, indent=2) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.markdown_output.write_text(write_markdown(record))


def write_markdown(record: dict) -> str:
    lines = ["# Splitting the held-out set on mechanism targets instead of drugs", "",
             f"Measured by `experiments/label_review/target_split_options.py` over {record['redraws']} full redraws at a "
             f"held-out fraction of {record['held_out_fraction']:.3f}, on {record['drugs']} drugs and "
             f"{record['distinct_targets']} distinct targets, with the GtoPdb narrowing and every addition kept.",
             "",
             "`fully covered` is the share of the held-out side whose mechanisms all appear on the training side: the "
             "quantity that condemned the drug split, where it is 84.5 percent. `unseen mechanism` is its complement, "
             "and is the share that supports a generalisation claim.",
             "",
             "| rule | held-out share | spread | fully covered by training mechanisms | has a mechanism training never saw |",
             "| --- | --- | --- | --- | --- |"]
    for entry in record["readings"]:
        def show(key, percent=True):
            value = entry[key]
            return "-" if value is None else (f"{value:.1%}" if percent else f"{value:.3f}")
        lines.append(f"| {entry['rule']} | {show('held_out_share')} | {show('held_out_share_spread', False)} | "
                     f"{show('fully_covered_by_training_mechanisms')} | {show('has_a_mechanism_training_never_saw')} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    main()
