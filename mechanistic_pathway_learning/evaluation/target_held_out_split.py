"""Drawing the held-out set by holding out mechanism targets (preregistration amendment item 25, 13 October 2026).

The rule, and what is principled about it. A drug is held out **whole** when any of its mechanism targets is in the
held-out target partition, because the claim under test is about mechanism and not about drug identity, and this is
what puts a mechanism the model never trained on into every held-out drug. The owner, 13 October 2026: "NOT splitting
on drug itself but splitting on each target individually"; and of the drug split it replaces, "i have to just accept
that there will be some leakage" - which this is what avoids.

Why the drug split could not stand, measured over 200 draws: 84.4 percent of held-out drugs had **every** target
already present in a training drug, because the drugs genuinely share mechanisms, and grouping hard enough to remove
that gives one block holding 90 percent of them.

What is a convention, said plainly rather than dressed up (the owner: "that's not principled; it's just pushing the
arbitrariness further down the chain"). The **size** of the held-out set is a convention. Inside a one-shot lockbox
design nothing can make it otherwise; only cross-validation removes the choice, and leave-one-target-out, which has no
convention in it at all, costs 192 folds against the registered 60 runs for the same split quality. One seventh is
legitimate because it was fixed long before this rule existed, not because it is derived.

A gene perturbation of a held-out target goes with it. Otherwise the model learns what that node does from the gene
and is credited for it on the drug, which is the leak `merge_drugs_with_their_targets` was written to stop and which
this rule inherits rather than drops.
"""
from __future__ import annotations

import random

REGISTERED_HELD_OUT_FRACTION = 1 / 7
# a draw stops growing once it is this far past the registered fraction, since no later prefix can be nearer
STOP_SEARCHING_PAST = 0.05


def targets_of_drugs(perturbation_ids, perturbation_types, perturbation_seeds, node_ids) -> dict[str, set[str]]:
    """perturbation id -> the node ids it seeds, for the drug perturbations only."""
    names = list(node_ids)
    return {str(perturbation_id): {names[int(seed)] for seed in seeds if 0 <= int(seed) < len(names)}
            for perturbation_id, kind, seeds in zip(perturbation_ids, perturbation_types, perturbation_seeds)
            if kind == "drug"}


def held_out_targets(targets_of_drug: dict[str, set[str]], seed: int,
                     held_out_fraction: float = REGISTERED_HELD_OUT_FRACTION) -> set[str]:
    """The held-out target partition: targets in the seed's order, kept to the prefix whose held-out share of drugs
    lands nearest the registered fraction.

    The count of targets is a consequence of that fraction and of how promiscuous the drugs happen to be, so it moves
    on its own when a rebuild changes the target sets. Which targets, and in which order, is the seed's alone."""
    drugs = sorted(targets_of_drug)
    if not drugs:
        return set()
    every_target = sorted({target for own in targets_of_drug.values() for target in own})
    order = list(every_target)
    random.Random(seed).shuffle(order)
    held, best_gap, best = set(), abs(0.0 - held_out_fraction), set()
    for target in order:
        held = held | {target}
        share = sum(1 for drug in drugs if targets_of_drug[drug] & held) / len(drugs)
        gap = abs(share - held_out_fraction)
        if gap < best_gap:
            best_gap, best = gap, set(held)
        elif share > held_out_fraction + STOP_SEARCHING_PAST:
            break  # no longer prefix can be nearer, since the share only grows
    return best


def held_out_perturbations(perturbation_ids, perturbation_types, perturbation_seeds, node_ids, seed: int,
                           held_out_fraction: float = REGISTERED_HELD_OUT_FRACTION,
                           gene_node_of_perturbation: dict[str, str] | None = None) -> tuple[set[str], set[str]]:
    """(the held-out perturbation ids, the held-out target node ids).

    A drug is held out when any of its targets is. A gene perturbation whose own node is a held-out target is held out
    with it, so the model cannot learn that node's effect from the gene and be credited for it on the drug."""
    targets_of_drug = targets_of_drugs(perturbation_ids, perturbation_types, perturbation_seeds, node_ids)
    held_targets = held_out_targets(targets_of_drug, seed, held_out_fraction)
    held = {drug for drug, own in targets_of_drug.items() if own & held_targets}
    names = list(node_ids)
    for perturbation_id, kind, seeds in zip(perturbation_ids, perturbation_types, perturbation_seeds):
        if kind == "drug":
            continue
        own = {names[int(index)] for index in seeds if 0 <= int(index) < len(names)}
        if gene_node_of_perturbation and str(perturbation_id) in gene_node_of_perturbation:
            own = own | {gene_node_of_perturbation[str(perturbation_id)]}
        if own & held_targets:
            held.add(str(perturbation_id))
    return held, held_targets


def mechanisms_never_trained(targets_of_drug: dict[str, set[str]], held: set[str]) -> dict[str, set[str]]:
    """Per held-out drug, the targets no training drug has: what the unseen-mechanism stratum of item 25c is scored on.

    Under this rule every held-out drug has at least one, by construction. A drug for which it is **all** of them is in
    the stratum; the rest have part of their mechanism in training, which the registration states as the price of
    keeping the label exact."""
    trained_targets = {target for drug, own in targets_of_drug.items() if drug not in held for target in own}
    return {drug: targets_of_drug[drug] - trained_targets for drug in held if drug in targets_of_drug}
