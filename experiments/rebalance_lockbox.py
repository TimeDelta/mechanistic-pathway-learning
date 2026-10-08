"""Derive configs/lockbox_v2.json from configs/lockbox_v1.json: return lockbox groups holding a drug to development until
the lockbox holds the drug share of the whole data, then move development groups that no full-graph model has read into
the lockbox until it is as large as before (the user's decision of 8 October 2026; docs/preregistration.md).

Why. lockbox_v1 holds 64 of the 142 drugs (45 percent), because the draw's drug-stratum target counted cluster:ABCA3
(232 perturbations, kept in development for its size) while only the other groups could be drawn (draw_lockbox.py). All
34 development drugs with a kept positive are in cluster:ABCA3, so no development reading trained on drug positives and
then tested on drug positives elsewhere, and a lockbox run's validation set held no drug kept positive. The lockbox's
drug groups share no target node with any development perturbation (merge_drugs_with_their_targets joins a drug to every
perturbation targeting one of its nodes), so whole drug groups can move without a group straddling the two sides.

Rule (fixed before it was run; membership, perturbation types and group sizes only, no label and no model output):
1. Lockbox groups holding a drug are put in a seeded random order (numpy default_rng(--seed) over the sorted group ids).
   In that order a group returns to development when the lockbox still holds at least round(--drug-share x all drugs)
   drugs after it leaves; a group that would take the lockbox below that is skipped.
2. Development groups are eligible to enter the lockbox when they hold no drug and none of their perturbations has been
   read by a full-graph model or by the slice runs: not in the test fold or the early-stopping validation set of the
   development pilots (lockbox_v1 development set, grouped fold 0 of 5, seed 0, validation fraction 0.15 with the large
   groups kept in training, as experiments/run_main_model.py drew them) and not a perturbation of the slice evidence
   (data/processed/evidence). Eligible groups are put in a seeded random order and enter the lockbox, in that order,
   until it holds at least as many perturbations as lockbox_v1.
3. As in draw_lockbox.py, the record lists the label counts inside the new lockbox, which fix the symptoms of its macro
   average (five kept positives or more) before anything is scored.
Drug groups that return to development were never read by a model either (they were in the lockbox), so the
development set gains drug positives outside cluster:ABCA3 without anything entering the lockbox that a decision has
seen. What remains: the swapped-in genes were training perturbations of the fold-0 development pilots (fitted, not
scored), and 85 lockbox_v1 perturbations that stay are slice perturbations, as before.

Usage:
  python experiments/rebalance_lockbox.py --output configs/lockbox_v2.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from draw_lockbox import MINIMUM_POSITIVES_TO_SCORE  # noqa: E402
from run_main_model import early_stopping_validation  # noqa: E402

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, restrict_to_perturbations  # noqa: E402
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds  # noqa: E402


def file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def groups_returned_to_development(drug_groups: list[str], drugs_in_group: dict[str, int], lockbox_drugs: int, minimum_drugs: int, seed: int) -> list[str]:
    """Rule step 1: lockbox drug groups, in a seeded order, that leave while the lockbox keeps at least minimum_drugs drugs."""
    order = np.random.default_rng(seed).permutation(len(drug_groups))
    returned = []
    for position in order:
        group = sorted(drug_groups)[position]
        if lockbox_drugs - drugs_in_group[group] >= minimum_drugs:
            returned.append(group)
            lockbox_drugs -= drugs_in_group[group]
    return returned


def groups_entering_lockbox(eligible_groups: list[str], size_of_group: dict[str, int], current_size: int, target_size: int, seed: int) -> list[str]:
    """Rule step 2: eligible development groups, in a seeded order, until the lockbox holds target_size perturbations."""
    order = np.random.default_rng(seed + 1).permutation(len(eligible_groups))
    entering = []
    for position in order:
        if current_size >= target_size:
            break
        group = sorted(eligible_groups)[position]
        entering.append(group)
        current_size += size_of_group[group]
    return entering


def read_by_development_pilots(data, in_lockbox: np.ndarray, num_folds: int = 5, fold: int = 0, seed: int = 0, validation_fraction: float = 0.15) -> set[str]:
    """Perturbation ids in the test fold and the early-stopping validation set of the lockbox_v1 development pilots."""
    development = restrict_to_perturbations(data, ~in_lockbox)
    fold_of = assign_grouped_folds(development.perturbation_ids, development.group_ids, num_folds, seed)
    test_ids = {p for p in development.perturbation_ids if fold_of[p] == fold}
    pool = np.array([i for i, p in enumerate(development.perturbation_ids) if fold_of[p] != fold], dtype=int)
    _, validation = early_stopping_validation(development, pool, SimpleNamespace(validation_fraction=validation_fraction, keep_large_groups_in_training=True, seed=seed))
    return test_ids | {development.perturbation_ids[i] for i in validation}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v2.parquet"))
    parser.add_argument("--slice-evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--slice-graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v1.json"))
    parser.add_argument("--drug-share", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--output", type=Path, default=Path("configs/lockbox_v2.json"))
    parser.add_argument("--overwrite", action="store_true", help="replace an existing lockbox file (only before anything was scored on it)")
    arguments = parser.parse_args()
    if arguments.output.exists() and not arguments.overwrite:
        raise SystemExit(f"{arguments.output} exists; a lockbox is drawn once (pass --overwrite only before anything was scored on it)")
    previous = json.loads(arguments.lockbox.read_text())
    group_by = previous["group_by"]
    if file_sha256(arguments.evidence_dir / "evidence_records.parquet") != previous["evidence_records_sha256"]:
        raise SystemExit(f"{arguments.lockbox} was drawn on other evidence than {arguments.evidence_dir}")
    if file_sha256(Path(previous["label_selection"])) != previous["label_selection_sha256"]:
        raise SystemExit(f"the label selection {previous['label_selection']} changed since {arguments.lockbox} was drawn")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=group_by, label_selection=arguments.label_selection)
    previous_data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=group_by, label_selection=Path(previous["label_selection"]))
    if previous_data.perturbation_ids != data.perturbation_ids or previous_data.group_ids != data.group_ids or not np.array_equal(previous_data.outcomes, data.outcomes):
        raise SystemExit("the new selection loads other perturbations, groups or positives than the one the lockbox was drawn on")
    in_previous = np.array([p in set(previous["perturbation_ids"]) for p in data.perturbation_ids])
    is_drug = np.array([kind == "drug" for kind in data.perturbation_types])
    size_of_group = Counter(data.group_ids)
    drugs_in_group = Counter(group for group, drug in zip(data.group_ids, is_drug) if drug)
    lockbox_groups = set(previous["groups"])
    if {group for group, inside in zip(data.group_ids, in_previous) if inside} != lockbox_groups:
        raise SystemExit(f"{arguments.lockbox}: its groups do not match its perturbations under group_by {group_by}")
    lockbox_drug_target_nodes = {data.node_ids[s] for i in np.flatnonzero(in_previous & is_drug) for s in data.perturbation_seeds[i]}
    development_nodes = {data.node_ids[s] for i in np.flatnonzero(~in_previous) for s in data.perturbation_seeds[i]}
    if lockbox_drug_target_nodes & development_nodes:
        raise SystemExit("a lockbox drug shares a target node with a development perturbation; groups cannot move whole")

    minimum_drugs = int(round(arguments.drug_share * is_drug.sum()))
    lockbox_drug_groups = sorted(group for group in lockbox_groups if drugs_in_group[group] > 0)
    returned = groups_returned_to_development(lockbox_drug_groups, drugs_in_group, int((in_previous & is_drug).sum()), minimum_drugs, arguments.seed)

    already_read = read_by_development_pilots(data, in_previous)
    slice_ids = set(load_experiment_data(arguments.slice_graph_dir, arguments.slice_evidence_dir, group_by="disease_cluster").perturbation_ids)
    members_of = defaultdict(list)
    for perturbation_id, group in zip(data.perturbation_ids, data.group_ids):
        members_of[group].append(perturbation_id)
    eligible = sorted(group for group in set(data.group_ids) - lockbox_groups
                      if drugs_in_group[group] == 0 and not (set(members_of[group]) & (already_read | slice_ids)))
    size_after_return = int(in_previous.sum()) - sum(size_of_group[group] for group in returned)
    entering = groups_entering_lockbox(eligible, dict(size_of_group), size_after_return, int(in_previous.sum()), arguments.seed)

    chosen = (lockbox_groups - set(returned)) | set(entering)
    in_lockbox = np.array([group in chosen for group in data.group_ids])
    kept = data.outcomes * data.label_mask
    positives_by_symptom = {symptom: int(kept[in_lockbox, column].sum()) for column, symptom in enumerate(data.symptoms)}
    record = {
        "drawn_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rule": {"derived_from": str(arguments.lockbox), "derived_from_sha256": file_sha256(arguments.lockbox), "drug_share": arguments.drug_share,
                 "minimum_lockbox_drugs": minimum_drugs, "seed": arguments.seed,
                 "step_1": "lockbox groups holding a drug, in numpy default_rng(seed).permutation order over their sorted ids, return to development while the lockbox keeps at least minimum_lockbox_drugs drugs",
                 "step_2": "development groups with no drug and no perturbation in the development pilots' fold-0 test or validation set or in the slice evidence, in numpy default_rng(seed + 1).permutation order over their sorted ids, enter until the lockbox is as large as before",
                 "read_by_development_pilots": len(already_read), "slice_perturbations": len(slice_ids), "eligible_groups": len(eligible),
                 "eligible_perturbations": sum(size_of_group[group] for group in eligible)},
        "group_by": group_by,
        "evidence_records_sha256": previous["evidence_records_sha256"],
        "label_selection": str(arguments.label_selection),
        "label_selection_sha256": file_sha256(arguments.label_selection),
        "groups_kept_in_development_for_size": previous["groups_kept_in_development_for_size"],
        "groups_returned_to_development": sorted(returned),
        "groups_entering_from_development": sorted(entering),
        "num_perturbations": int(in_lockbox.sum()), "num_development_perturbations": int((~in_lockbox).sum()),
        "num_genes": int((in_lockbox & ~is_drug).sum()), "num_drugs": int((in_lockbox & is_drug).sum()),
        "num_drugs_with_a_kept_positive": int(((kept[in_lockbox & is_drug]).sum(axis=1) > 0).sum()),
        "kept_positive_pairs": int(kept[in_lockbox].sum()), "kept_positive_pairs_development": int(kept[~in_lockbox].sum()),
        "kept_positives_by_symptom": positives_by_symptom,
        "symptoms_in_the_macro_average": sorted(symptom for symptom, count in positives_by_symptom.items() if count >= MINIMUM_POSITIVES_TO_SCORE),
        "groups": sorted(chosen),
        "perturbation_ids": sorted(p for p, inside in zip(data.perturbation_ids, in_lockbox) if inside),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({key: value for key, value in record.items() if key not in ("groups", "perturbation_ids")}, indent=1))


if __name__ == "__main__":
    main()
