"""Evaluation splits (experiment design, section 6.1).

Three regimes, each stricter than random edge splits:

- Grouped perturbation-wise folds. Every perturbation (a gene or a drug) is in
  exactly one fold, and all perturbations sharing a group (for drugs, the
  dominant target; for genes, the gene itself) land in the same fold. This is the
  disjoint-group protocol that removes the optimism of random splits.
- Pathway-wise folds. All perturbations anchored in one curated module are held
  out together, testing generalization to an unseen mechanism.
- Time split. Observations are divided by the date the evidence became available.
"""
from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Mapping, Sequence
from datetime import date


def assign_grouped_folds(
    perturbation_identifiers: Sequence[str],
    group_identifiers: Sequence[str],
    num_folds: int,
    random_seed: int = 0,
) -> dict[str, int]:
    """Return perturbation identifier -> fold index, keeping each group inside one fold.

    Groups are assigned greedily, largest first with ties shuffled, to the fold that
    currently holds the fewest perturbations, which balances fold sizes.
    """
    if len(perturbation_identifiers) != len(group_identifiers):
        raise ValueError("perturbation_identifiers and group_identifiers must have the same length")
    if num_folds < 2:
        raise ValueError("num_folds must be at least 2")
    members_by_group: dict[str, list[str]] = defaultdict(list)
    for perturbation_identifier, group_identifier in zip(perturbation_identifiers, group_identifiers):
        members_by_group[group_identifier].append(perturbation_identifier)
    group_order = list(members_by_group)
    random.Random(random_seed).shuffle(group_order)
    group_order.sort(key=lambda group_identifier: -len(members_by_group[group_identifier]))
    fold_sizes = [0] * num_folds
    fold_by_perturbation: dict[str, int] = {}
    for group_identifier in group_order:
        smallest_fold = min(range(num_folds), key=lambda fold_index: fold_sizes[fold_index])
        for perturbation_identifier in members_by_group[group_identifier]:
            fold_by_perturbation[perturbation_identifier] = smallest_fold
        fold_sizes[smallest_fold] += len(members_by_group[group_identifier])
    return fold_by_perturbation


def pathway_wise_holdout_sets(perturbation_to_curated_module: Mapping[str, str]) -> dict[str, list[str]]:
    """Return curated module identifier -> perturbations held out together for that module."""
    held_out_by_module: dict[str, list[str]] = defaultdict(list)
    for perturbation_identifier, module_identifier in perturbation_to_curated_module.items():
        held_out_by_module[module_identifier].append(perturbation_identifier)
    return dict(held_out_by_module)


def time_split(
    observation_identifiers: Sequence[str],
    evidence_available_dates: Sequence[date],
    cutoff_date: date,
) -> tuple[list[str], list[str]]:
    """Return (training identifiers, test identifiers) split by evidence availability date.

    Observations dated on or before the cutoff train; later ones test. Observations
    with no date are excluded from both sets, since their position is unknown.
    """
    if len(observation_identifiers) != len(evidence_available_dates):
        raise ValueError("observation_identifiers and evidence_available_dates must have the same length")
    training_identifiers: list[str] = []
    test_identifiers: list[str] = []
    for observation_identifier, evidence_date in zip(observation_identifiers, evidence_available_dates):
        if evidence_date is None:
            continue
        if evidence_date <= cutoff_date:
            training_identifiers.append(observation_identifier)
        else:
            test_identifiers.append(observation_identifier)
    return training_identifiers, test_identifiers
