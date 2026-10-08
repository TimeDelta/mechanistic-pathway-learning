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


def read_curated_modules(curated_modules_path) -> dict[str, list[str]]:
    """docs/curated_pathway_modules.csv -> module id -> gene symbols (the scaffolding of design section 3.2)."""
    import csv

    genes_by_module: dict[str, list[str]] = {}
    with open(curated_modules_path, encoding="utf-8") as modules_file:
        for row in csv.DictReader(modules_file):
            genes_by_module[row["module_id"]] = [symbol.strip() for symbol in (row.get("genes") or "").split(";") if symbol.strip()]
    return genes_by_module


def perturbations_anchored_in_module(perturbation_seed_node_indices: Sequence, module_node_indices: set[int]) -> list[bool]:
    """True for every perturbation that writes onto at least one node of the module (a module gene, or a drug targeting one).

    This is the pathway-wise hold-out set: all of them leave the training data together.
    """
    return [any(int(node_index) in module_node_indices for node_index in seeds) for seeds in perturbation_seed_node_indices]


def training_mask_without_group_partners(test_mask: Sequence[bool], group_ids: Sequence[str]) -> list[bool]:
    """Training rows for a hold-out: every perturbation outside it whose leakage group has no member inside it.

    A grouped fold or the lockbox holds whole groups, so there this is the complement of the test mask. A pathway-wise
    hold-out selects perturbations by seed gene instead, and a perturbation outside it that shares a disease cluster or
    a drug target with one inside would carry the held-out labels into training; such partners are neither trained on
    nor scored.
    """
    held_out_groups = {group for group, held_out in zip(group_ids, test_mask) if held_out}
    return [not held_out and group not in held_out_groups for group, held_out in zip(group_ids, test_mask)]


CATCH_ALL_SUBSYSTEMS: frozenset[str] = frozenset({
    "Transport reactions", "Exchange/demand reactions", "Isolated", "Miscellaneous", "Artificial reactions", "Pool reactions",
})


def primary_subsystem_by_gene_node(
    node_subsystem: Sequence[str],
    edge_source: Sequence[int],
    edge_target: Sequence[int],
    edge_relation: Sequence[int],
    catalyzed_by_relation_index: int,
    ignored_subsystems: frozenset[str] = CATCH_ALL_SUBSYSTEMS,
) -> dict[int, str]:
    """Gene node index -> the subsystem most of its catalyzed reactions belong to (ties broken alphabetically).

    Reconstruction subsystems (Human-GEM) partition the reactions into pathways such as "Tryptophan
    metabolism"; holding out every gene whose primary subsystem is S is a pathway-wise split that is
    populated far more evenly than the sixteen hand-curated modules of design section 3.2. Catch-all
    subsystems (transport, exchange, isolated reactions) are not pathways: a gene is assigned to its most
    frequent pathway subsystem when it has one, and to the catch-all only when it has nothing else, so
    that a transporter of the urea cycle is held out with the urea cycle and not with every transporter.
    """
    from collections import Counter, defaultdict

    counts_by_gene: dict[int, Counter] = defaultdict(Counter)
    for source, target, relation in zip(edge_source, edge_target, edge_relation):
        if relation != catalyzed_by_relation_index:
            continue
        subsystem = node_subsystem[target]
        if subsystem:
            counts_by_gene[int(source)][subsystem] += 1
    primary: dict[int, str] = {}
    for gene_node, counts in counts_by_gene.items():
        pathway_counts = {name: count for name, count in counts.items() if name not in ignored_subsystems} or dict(counts)
        best_count = max(pathway_counts.values())
        primary[gene_node] = sorted(name for name, count in pathway_counts.items() if count == best_count)[0]
    return primary


def subsystem_holdout_masks(
    perturbation_seed_node_indices: Sequence,
    primary_subsystem: Mapping[int, str],
    outcomes_per_perturbation: Sequence[float],
    min_holdout_positives: int,
) -> dict[str, list[bool]]:
    """Subsystem -> mask of perturbations whose seed nodes include a gene of that primary subsystem, keeping subsystems with enough positives.

    With a label selection, outcomes_per_perturbation should count the kept positives only (outcomes x label mask): the
    pairs set aside are not scored, so they cannot make a hold-out scorable."""
    nodes_by_subsystem: dict[str, set[int]] = defaultdict(set)
    for gene_node, subsystem in primary_subsystem.items():
        nodes_by_subsystem[subsystem].add(gene_node)
    masks: dict[str, list[bool]] = {}
    for subsystem, node_set in sorted(nodes_by_subsystem.items()):
        if subsystem in CATCH_ALL_SUBSYSTEMS:
            continue
        mask = perturbations_anchored_in_module(perturbation_seed_node_indices, node_set)
        if sum(outcome for outcome, held_out in zip(outcomes_per_perturbation, mask) if held_out) >= min_holdout_positives:
            masks[subsystem] = mask
    return masks
