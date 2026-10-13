"""Drawing the held-out set by mechanism target (preregistration amendment item 25; the owner, 13 October 2026: "NOT
splitting on drug itself but splitting on each target individually")."""
import pytest

from mechanistic_pathway_learning.evaluation.target_held_out_split import (
    REGISTERED_HELD_OUT_FRACTION,
    held_out_perturbations,
    held_out_targets,
    mechanisms_never_trained,
    targets_of_drugs,
)

# six drugs over five targets; d1 and d2 share PROTEIN:A, which is what the drug split would have merged
NODES = ["PROTEIN:A", "PROTEIN:B", "PROTEIN:C", "PROTEIN:D", "PROTEIN:E", "GENE:A"]
IDS = ["d1", "d2", "d3", "d4", "d5", "d6", "gA"]
KINDS = ["drug"] * 6 + ["gene"]
SEEDS = [[0, 1], [0], [2], [3], [4], [1, 4], [5]]


def test_a_drug_is_held_out_whole_when_any_of_its_targets_is():
    targets = targets_of_drugs(IDS, KINDS, SEEDS, NODES)
    assert targets["d1"] == {"PROTEIN:A", "PROTEIN:B"} and "gA" not in targets  # genes are not drugs here
    held, held_targets = held_out_perturbations(IDS, KINDS, SEEDS, NODES, seed=0, held_out_fraction=0.34)
    for drug, own in targets.items():
        assert (drug in held) == bool(own & held_targets), drug  # the rule, exactly
    assert held_targets and held


def test_every_held_out_drug_carries_a_mechanism_training_never_saw():
    """The property the rule exists for, and the one the drug split could not give: 84.4 percent of its held-out drugs
    had every target already in training."""
    targets = targets_of_drugs(IDS, KINDS, SEEDS, NODES)
    for seed in range(25):
        held, _ = held_out_perturbations(IDS, KINDS, SEEDS, NODES, seed=seed, held_out_fraction=0.34)
        unseen = mechanisms_never_trained(targets, held)
        assert all(missing for missing in unseen.values()), seed


def test_a_gene_perturbation_of_a_held_out_target_goes_with_it():
    """Otherwise the model learns the node's effect from the gene and is credited for it on the drug, which is the leak
    merge_drugs_with_their_targets was written to stop."""
    nodes = ["PROTEIN:A", "PROTEIN:B"]
    ids, kinds, seeds = ["d1", "d2", "gB"], ["drug", "drug", "gene"], [[0], [1], [1]]
    held, held_targets = held_out_perturbations(ids, kinds, seeds, nodes, seed=0, held_out_fraction=0.5)
    assert ("gB" in held) == ("PROTEIN:B" in held_targets)
    # and through an explicit gene-node map, for a split graph where the gene seeds its own node
    held_mapped, held_mapped_targets = held_out_perturbations(
        ids, kinds, [[0], [1], []], nodes, seed=0, held_out_fraction=0.5,
        gene_node_of_perturbation={"gB": "PROTEIN:B"})
    assert ("gB" in held_mapped) == ("PROTEIN:B" in held_mapped_targets)


def test_the_held_out_share_lands_near_the_registered_fraction_and_the_target_count_follows():
    """The size is the registered convention; the number of targets is whatever delivers it, and is not chosen."""
    targets = {f"d{index}": {f"t{index % 17}", f"t{(index * 7) % 17}"} for index in range(140)}
    shares, counts = [], []
    for seed in range(60):
        held = held_out_targets(targets, seed)
        counts.append(len(held))
        shares.append(sum(1 for own in targets.values() if own & held) / len(targets))
    mean_share = sum(shares) / len(shares)
    assert abs(mean_share - REGISTERED_HELD_OUT_FRACTION) < 0.05, mean_share
    assert len(set(counts)) > 1  # the count varies with the draw rather than being fixed in advance


def test_the_draw_is_the_seeds_alone_and_is_reproducible():
    first = held_out_targets(targets_of_drugs(IDS, KINDS, SEEDS, NODES), seed=4, held_out_fraction=0.34)
    again = held_out_targets(targets_of_drugs(IDS, KINDS, SEEDS, NODES), seed=4, held_out_fraction=0.34)
    other = held_out_targets(targets_of_drugs(IDS, KINDS, SEEDS, NODES), seed=5, held_out_fraction=0.34)
    assert first == again and first != other


def test_the_unseen_mechanism_stratum_is_the_drugs_with_nothing_left_in_training():
    targets = {"d1": {"a"}, "d2": {"a", "b"}, "d3": {"c"}}
    unseen = mechanisms_never_trained(targets, {"d1", "d2"})
    assert unseen["d1"] == {"a"}  # d1's only target is held out with it, so nothing of it is trained
    assert unseen["d2"] == {"a", "b"}
    partly = mechanisms_never_trained({"d1": {"a"}, "d2": {"a", "b"}}, {"d2"})
    assert partly["d2"] == {"b"}  # "a" is still trained, through d1


def test_an_empty_table_and_a_drug_with_no_target_are_handled():
    assert held_out_targets({}, seed=0) == set()
    held, held_targets = held_out_perturbations(["d1"], ["drug"], [[]], ["PROTEIN:A"], seed=0)
    assert held == set() and held_targets == set()
