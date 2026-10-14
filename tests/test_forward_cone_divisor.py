"""Decision 9.16, reversed to per seed by the owner on 14 October 2026 ("oh ... do the per-seed cone divisor then").

The chain graph in the first test is the smallest case that shows the point: a destination fed by one live edge and
one dead edge divides the live message by two under item 28 and by one under the cone divisor.
"""
import torch

from mechanistic_pathway_learning.models.forward_cone_divisor import (
    cone_restricted_in_degree,
    divisor_inflation_against_the_union,
    forward_reachable_nodes,
    union_of_seed_cones,
)


def test_the_cone_holds_the_seed_and_what_it_reaches_and_nothing_else():
    # 0 -> 1 -> 2, and an orphan 3 -> 1 that no seed can reach
    edge_source = torch.tensor([0, 1, 3])
    edge_target = torch.tensor([1, 2, 1])
    reachable = forward_reachable_nodes(edge_source, edge_target, torch.tensor([0]), num_graph_nodes=4)
    assert reachable.tolist() == [True, True, True, False]


def test_a_dead_edge_leaves_the_divisor_so_the_live_message_is_not_halved():
    edge_source = torch.tensor([0, 1, 3])
    edge_target = torch.tensor([1, 2, 1])
    edge_weight = torch.ones(3)
    reachable = forward_reachable_nodes(edge_source, edge_target, torch.tensor([0]), num_graph_nodes=4)
    cone_divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight, reachable, num_graph_nodes=4)
    every_edge_divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight,
                                                   torch.ones(4, dtype=torch.bool), num_graph_nodes=4)
    assert cone_divisor[1].item() == 1.0, "only the edge from inside the cone counts"
    assert every_edge_divisor[1].item() == 2.0, "item 28 counts the orphan edge too"


def test_the_row_sum_stays_at_most_one_so_the_contraction_is_untouched():
    edge_source = torch.tensor([0, 0, 1, 3])
    edge_target = torch.tensor([1, 2, 2, 2])
    edge_weight = torch.tensor([1.0, 2.0, 3.0, 4.0])
    reachable = forward_reachable_nodes(edge_source, edge_target, torch.tensor([0]), num_graph_nodes=4)
    divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight, reachable, num_graph_nodes=4)
    source_is_in_the_cone = reachable[edge_source]
    normalised = torch.where(source_is_in_the_cone, edge_weight.abs(), torch.zeros_like(edge_weight)) / divisor[edge_target]
    row_sums = torch.zeros(4).index_add_(0, edge_target, normalised)
    assert bool((row_sums <= 1.0 + 1e-6).all()), f"a row sum exceeded one: {row_sums.tolist()}"


def test_the_union_and_the_per_seed_divisors_are_not_identical():
    """Two seeds with disjoint cones meeting at one destination: the union divides by both, per seed by one each."""
    edge_source = torch.tensor([0, 1])
    edge_target = torch.tensor([2, 2])
    edge_weight = torch.ones(2)
    seeds = [torch.tensor([0]), torch.tensor([1])]
    union_reachable = union_of_seed_cones(edge_source, edge_target, seeds, num_graph_nodes=3)
    union_divisor = cone_restricted_in_degree(edge_source, edge_target, edge_weight, union_reachable, num_graph_nodes=3)
    first_seed_divisor = cone_restricted_in_degree(
        edge_source, edge_target, edge_weight,
        forward_reachable_nodes(edge_source, edge_target, seeds[0], num_graph_nodes=3), num_graph_nodes=3)
    assert union_divisor[2].item() == 2.0
    assert first_seed_divisor[2].item() == 1.0

    report = divisor_inflation_against_the_union(edge_source, edge_target, edge_weight, seeds, num_graph_nodes=3)
    assert report["median_union_over_per_seed"] == 2.0
    assert report["share_of_destinations_where_they_agree"] < 1.0


def test_seeds_that_share_a_cone_make_the_two_divisors_agree():
    edge_source = torch.tensor([0, 1, 2])
    edge_target = torch.tensor([2, 2, 3])
    edge_weight = torch.ones(3)
    seeds = [torch.tensor([0, 1]), torch.tensor([0, 1])]
    report = divisor_inflation_against_the_union(edge_source, edge_target, edge_weight, seeds, num_graph_nodes=4)
    assert report["share_of_destinations_where_they_agree"] == 1.0
