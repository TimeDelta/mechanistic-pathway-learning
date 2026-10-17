"""experiments/measure_perturbation_reach.py: which perturbations' limiting-step routes stop at their seeds, and why."""
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import torch

from experiments.measure_perturbation_reach import (REASON_ENCODER_DOES_NOT_WALK, REASON_INTO_WEIGHTLESS_NODES, REASON_LOADER_DROPS,
                                                    REASON_NO_EDGE_IN_FILE, REASON_ONLY_INTO_DRUG_NODES, REASON_SEED_WITHOUT_WEIGHT,
                                                    measure, reachable_nodes, reasons_routes_stop, walkable_adjacency)
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder


def chain_adjacency(num_nodes, edges, weighted=None, entry=None):
    source = np.array([edge[0] for edge in edges], dtype=int)
    target = np.array([edge[1] for edge in edges], dtype=int)
    weighted = np.ones(num_nodes, dtype=bool) if weighted is None else np.asarray(weighted, dtype=bool)
    entry = np.zeros(num_nodes, dtype=bool) if entry is None else np.asarray(entry, dtype=bool)
    return walkable_adjacency(num_nodes, source, target, weighted, entry), weighted, entry


def test_a_route_runs_along_walked_edges_and_counts_its_seed():
    adjacency, weighted, _ = chain_adjacency(4, [(0, 1), (1, 2)])
    assert reachable_nodes(adjacency, [0], weighted) == {0, 1, 2}
    assert reachable_nodes(adjacency, [3], weighted) == {3}, "a seed with no edge out reaches itself"


def test_a_node_without_class_weight_ends_a_route_and_a_seed_without_it_starts_none():
    adjacency, weighted, _ = chain_adjacency(4, [(0, 1), (1, 2), (2, 3)], weighted=[True, True, False, True])
    assert reachable_nodes(adjacency, [0], weighted) == {0, 1}, "the recurrence's min with w = 0 stops the route at node 2"
    assert reachable_nodes(adjacency, [2], weighted) == set(), "b starts at the seed's class weight, 0 here"


def test_an_entry_node_is_reached_only_as_a_seed():
    """A drug node exists only in the perturbation that seeds it, so another perturbation's route cannot pass through it."""
    adjacency, weighted, _ = chain_adjacency(4, [(0, 1), (1, 2), (2, 3)], entry=[False, False, True, False])
    assert reachable_nodes(adjacency, [0], weighted) == {0, 1}
    assert reachable_nodes(adjacency, [2], weighted) == {2, 3}, "a drug seeded on its own node walks its mechanism edges"


def test_the_reasons_name_each_way_a_route_can_stop():
    weighted = np.array([True, True, False, True, True])
    entry = np.array([False, False, False, True, False])
    no_edges = reasons_routes_stop({0}, [0], weighted, entry, {}, {}, {})
    assert no_edges == ([REASON_NO_EDGE_IN_FILE], {})
    dropped = reasons_routes_stop({0}, [0], weighted, entry, {0: [(4, "binds")]}, {}, {})
    assert dropped[0] == [REASON_LOADER_DROPS] and dropped[1][REASON_LOADER_DROPS] == Counter({"binds": 1})
    not_walked = reasons_routes_stop({0}, [0], weighted, entry, {0: [(4, "activates")]}, {0: [(4, "activates")]}, {})
    assert not_walked[0] == [REASON_ENCODER_DOES_NOT_WALK] and not_walked[1][REASON_ENCODER_DOES_NOT_WALK] == Counter({"activates": 1})
    into_drug = reasons_routes_stop({0}, [0], weighted, entry, {0: [(3, "sequesters")]}, {0: [(3, "sequesters")]}, {0: {3}})
    assert into_drug == ([REASON_ONLY_INTO_DRUG_NODES], {})
    into_weightless = reasons_routes_stop({0}, [0], weighted, entry, {0: [(2, "activates")]}, {0: [(2, "activates")]}, {0: {2}})
    assert into_weightless == ([REASON_INTO_WEIGHTLESS_NODES], {})
    seed_without_weight = reasons_routes_stop(set(), [2], weighted, entry, {}, {}, {})
    assert seed_without_weight == ([REASON_SEED_WITHOUT_WEIGHT], {})


def test_an_edge_to_a_reached_node_is_no_reason():
    """The protein's edge back to its gene is walked and ends inside the reached set: no route leaves through it, and
    it is not why the routes stop."""
    weighted = np.ones(2, dtype=bool)
    entry = np.zeros(2, dtype=bool)
    reasons, _ = reasons_routes_stop({0, 1}, [0], weighted, entry, {0: [(1, "encodes")], 1: [(0, "binds")]},
                                     {0: [(1, "encodes")], 1: [(0, "binds")]}, {0: {1}, 1: {0}})
    assert reasons == [], "nothing blocks a route here; the reach is the whole component"


def toy_graph(tmp_path):
    """Nodes: 0 gene A, 1 protein A (no edge out), 2 gene B, 3 protein B, 4 a currency metabolite B makes, 5 a reaction
    the currency metabolite feeds, 6 drug D (an entry node), 7 a protein D targets, 8 gene C, 9 protein C, whose only
    edge in the file the loader did not load, and 10 a node without class weight downstream of protein C's other edge."""
    node_ids = ["gene:A", "protein:A", "gene:B", "protein:B", "metabolite:atp", "reaction:R", "drug:D", "protein:T", "gene:C", "protein:C",
                "protein:W"]
    node_types = ["gene", "protein", "gene", "protein", "metabolite", "reaction", "drug", "protein", "gene", "protein", "protein"]
    relation_types = ["encodes", "activates", "product_of", "substrate_of", "targets"]
    loaded_edges = [(0, 1, "encodes"), (2, 3, "encodes"), (3, 4, "activates"), (4, 5, "substrate_of"), (6, 7, "targets"), (8, 9, "encodes"),
                    (9, 10, "activates")]
    file_edges = loaded_edges + [(9, 3, "binds")]
    pd.DataFrame({"source_id": [node_ids[s] for s, _, _ in file_edges], "target_id": [node_ids[t] for _, t, _ in file_edges],
                  "relation_type": [relation for _, _, relation in file_edges]}).to_parquet(tmp_path / "edges.parquet")
    edge_source = torch.tensor([s for s, _, _ in loaded_edges])
    edge_target = torch.tensor([t for _, t, _ in loaded_edges])
    edge_relation = torch.tensor([relation_types.index(relation) for _, _, relation in loaded_edges])
    num_nodes = len(node_ids)
    class_weights = torch.ones((num_nodes, 2))
    class_weights[10] = 0.0
    entry = torch.zeros(num_nodes, dtype=torch.bool)
    entry[6] = True
    currency = torch.zeros(num_nodes, dtype=torch.bool)
    currency[4] = True
    torch.manual_seed(0)
    encoder = LinearResponseEncoder(num_nodes, relation_types, edge_source=edge_source, edge_target=edge_target, edge_relation=edge_relation,
                                    edge_sign=torch.ones(len(loaded_edges)), node_features=torch.eye(num_nodes), node_state_dim=4,
                                    non_propagating_nodes=currency, cell_class_weights=class_weights, class_weight_rule="limiting_step",
                                    entry_node_mask=entry, entry_relation_names=("targets",), edge_weight=torch.ones(len(loaded_edges)))
    perturbation_ids = ["A", "B", "D", "C"]
    data = SimpleNamespace(
        node_ids=node_ids, node_index={node: index for index, node in enumerate(node_ids)}, node_types=np.array(node_types),
        edge_source=edge_source.numpy(), edge_target=edge_target.numpy(), edge_relation=edge_relation.numpy(), relation_types=relation_types,
        perturbation_ids=perturbation_ids, perturbation_types=["gene", "gene", "drug", "gene"],
        perturbation_seeds=[np.array([0]), np.array([2]), np.array([6]), np.array([8])],
        outcomes=np.array([[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]]), label_mask=np.array([[True, True], [True, True], [True, True], [False, True]]))
    return data, encoder


def test_measure_reads_the_encoders_own_bottleneck_edges(tmp_path):
    data, encoder = toy_graph(tmp_path)
    table, summary = measure(data, encoder, tmp_path, dead_end_reach=2)
    by_id = table.set_index("perturbation_id")
    assert by_id.loc["A", "reachable_nodes"] == 2 and by_id.loc["A", "dead_end"]
    assert by_id.loc["A", "reasons"] == REASON_NO_EDGE_IN_FILE, "protein A has no edge out in the file"
    assert by_id.loc["B", "reachable_nodes"] == 3 and not by_id.loc["B", "dead_end"], "gene B, protein B and the currency metabolite"
    assert by_id.loc["D", "reachable_nodes"] == 2 and by_id.loc["D", "dead_end"], "the drug's mechanism edge is a route (decision 9.66)"
    assert REASON_NO_EDGE_IN_FILE in by_id.loc["D", "reasons"]
    assert by_id.loc["C", "dead_end"]
    assert set(by_id.loc["C", "reasons"].split("; ")) == {REASON_LOADER_DROPS, REASON_INTO_WEIGHTLESS_NODES}
    assert "binds:1" in by_id.loc["C", "relations_by_reason"]
    assert by_id.loc["C", "positive_pairs"] == 0, "a masked pair is not a positive"
    assert summary["dead_ends"] == 3 and summary["by_type"]["gene"]["dead_ends"] == 2
    assert summary["by_type"]["gene"]["positive_pairs_on_dead_ends"] == 1


def test_a_currency_nodes_edges_are_named_as_edges_the_encoder_does_not_walk(tmp_path):
    data, encoder = toy_graph(tmp_path)
    table, _ = measure(data, encoder, tmp_path, dead_end_reach=3)
    reasons = table.set_index("perturbation_id").loc["B", "reasons"]
    assert REASON_ENCODER_DOES_NOT_WALK in reasons, "a currency metabolite does not propagate, so its substrate_of edge is not a route"


def test_measure_needs_the_limiting_step_rules_class_weights(tmp_path):
    data, _ = toy_graph(tmp_path)
    encoder = SimpleNamespace(cell_class_node_weight=None)
    with pytest.raises(ValueError, match="limiting-step"):
        measure(data, encoder, tmp_path, dead_end_reach=2)
