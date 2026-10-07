"""Draw the lockbox: leakage groups set aside before any full-graph model is scored and scored once, by the named
confirmatory models, after every pilot is finished (docs/preregistration.md).

Rule (fixed before the draw):
- Units are the leakage groups of --group-by (disease_cluster_and_targets: a gene with the genes sharing a disease entry,
  a drug with the genes it targets and the drugs sharing a target), so no group straddles the lockbox.
- Groups holding more than --max-group-share of the perturbations stay in development: one of them would make the lockbox
  mostly one disease cluster (the largest group on evidence_full_v2 holds 232 of 1,539 perturbations).
- Two strata: groups holding a drug with a kept positive pair, and the rest. Within each stratum the groups are put in a
  random order (numpy default_rng(--seed)) and taken in that order until the lockbox holds at least --share of the
  stratum's perturbations, so drugs reach the lockbox in proportion.
- Nothing about any model or baseline enters the draw. The output records the label counts inside the lockbox, which
  decide which symptoms the macro average will take (five positives or more), so that list is fixed before scoring.

Writes configs/lockbox_v1.json (committed: perturbation and group ids, the rule, the seed and the SHA-256 of the
evidence table and label selection it was drawn on).

Usage:
  python experiments/draw_lockbox.py --output configs/lockbox_v1.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

MINIMUM_POSITIVES_TO_SCORE = 5


def draw_groups(group_ids: list[str], stratum_of_group: dict[str, str], size_of_group: dict[str, int], share: float, seed: int,
                excluded_groups: set[str]) -> set[str]:
    """Groups taken, per stratum, in a seeded random order until the stratum's lockbox reaches share of its perturbations."""
    generator = np.random.default_rng(seed)
    chosen: set[str] = set()
    for stratum in sorted(set(stratum_of_group.values())):
        groups = sorted(group for group in set(group_ids) if stratum_of_group[group] == stratum)
        target = share * sum(size_of_group[group] for group in groups)
        taken = 0
        for position in generator.permutation(len(groups)):
            group = groups[position]
            if taken >= target:
                break
            if group in excluded_groups:
                continue
            chosen.add(group)
            taken += size_of_group[group]
    return chosen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v1_full_v2.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--share", type=float, default=0.20)
    parser.add_argument("--max-group-share", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20261007)
    parser.add_argument("--output", type=Path, default=Path("configs/lockbox_v1.json"))
    parser.add_argument("--overwrite", action="store_true", help="replace an existing lockbox file (a lockbox is drawn once)")
    arguments = parser.parse_args()
    if arguments.output.exists() and not arguments.overwrite:
        raise SystemExit(f"{arguments.output} exists; a lockbox is drawn once (pass --overwrite only before anything was scored on it)")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    kept = data.outcomes * data.label_mask
    size_of_group: dict[str, int] = defaultdict(int)
    stratum_of_group: dict[str, str] = defaultdict(lambda: "no_drug_with_a_kept_positive")
    for row, (group, kind) in enumerate(zip(data.group_ids, data.perturbation_types)):
        size_of_group[group] += 1
        if kind == "drug" and kept[row].sum() > 0:
            stratum_of_group[group] = "drug_with_a_kept_positive"
    for group in size_of_group:
        stratum_of_group[group] = stratum_of_group[group]
    large_groups = {group for group, size in size_of_group.items() if size > arguments.max_group_share * len(data.perturbation_ids)}
    chosen = draw_groups(data.group_ids, dict(stratum_of_group), dict(size_of_group), arguments.share, arguments.seed, large_groups)
    in_lockbox = np.array([group in chosen for group in data.group_ids])
    lockbox_ids = sorted(p for p, inside in zip(data.perturbation_ids, in_lockbox) if inside)
    is_drug = np.array([kind == "drug" for kind in data.perturbation_types])
    positives_by_symptom = {symptom: int(kept[in_lockbox, column].sum()) for column, symptom in enumerate(data.symptoms)}
    record = {
        "drawn_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rule": {"group_by": arguments.group_by, "share": arguments.share, "max_group_share": arguments.max_group_share, "seed": arguments.seed,
                 "strata": "groups holding a drug with a kept positive pair; the rest", "order": "numpy default_rng(seed).permutation over the sorted group ids of each stratum"},
        "group_by": arguments.group_by,
        "evidence_records_sha256": hashlib.sha256((arguments.evidence_dir / "evidence_records.parquet").read_bytes()).hexdigest(),
        "label_selection": str(arguments.label_selection),
        "label_selection_sha256": hashlib.sha256(arguments.label_selection.read_bytes()).hexdigest(),
        "groups_kept_in_development_for_size": sorted(large_groups),
        "num_perturbations": int(in_lockbox.sum()), "num_development_perturbations": int((~in_lockbox).sum()),
        "num_genes": int((in_lockbox & ~is_drug).sum()), "num_drugs": int((in_lockbox & is_drug).sum()),
        "num_drugs_with_a_kept_positive": int(((kept[in_lockbox & is_drug]).sum(axis=1) > 0).sum()),
        "kept_positive_pairs": int(kept[in_lockbox].sum()), "kept_positive_pairs_development": int(kept[~in_lockbox].sum()),
        "kept_positives_by_symptom": positives_by_symptom,
        "symptoms_in_the_macro_average": sorted(symptom for symptom, count in positives_by_symptom.items() if count >= MINIMUM_POSITIVES_TO_SCORE),
        "groups": sorted(chosen),
        "perturbation_ids": lockbox_ids,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({key: value for key, value in record.items() if key not in ("groups", "perturbation_ids")}, indent=1))


if __name__ == "__main__":
    main()
