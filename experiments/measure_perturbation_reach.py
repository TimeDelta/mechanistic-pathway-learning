"""Which perturbations' limiting-step routes reach past their own seeds, and why the ones that stop there do.

Under the registered cell-class rule (class_weight_rule "limiting_step"), the encoder scales each node's settled response
by b, the weakest class weight on the widest route from the perturbed nodes (LinearResponseEncoder.class_bottleneck), so
a node with no route from the seeds carries no response, and a perturbation whose routes stop at its seeds is read
through the seeds' own states and nothing else. This script counts, for every perturbation, the nodes a route can reach,
and for each perturbation that reaches no more than --dead-end-reach nodes it states why the routes stop, for every
node they reach:

  - "no edge out in the graph file": no edge leaves the node in edges.parquet;
  - "edges the loader drops": edges leave it in the file, and load_experiment_data does not load them;
  - "edges the encoder does not walk": loaded edges leave it, and the bottleneck's edge set does not hold them (a
    currency node's edges, which do not propagate, or a relation the adjacency leaves out);
  - "edges only into drug nodes": the edges it walks end at drug entry nodes, which exist only in the perturbation that
    seeds them;
  - "edges into nodes without class weight": the edges it walks end at nodes whose class weight is 0 in every class.

The count is an upper bound on the nodes with a nonzero bottleneck, since a route here may change class at any node and
the recurrence lets it change class only at an extracellular pool; a perturbation at or below the bound has no larger
exact reach, so every dead end found here is one.

It takes the trainer's own arguments, so the encoder it reads is the one a run with those arguments trains, then its
own:

  python experiments/measure_perturbation_reach.py <the trainer's flags> \\
      --reach-output runs/perturbation_reach.csv --reach-summary-output runs/perturbation_reach.json
"""
from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.sparse.csgraph import breadth_first_order

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

DEFAULT_DEAD_END_REACH = 2  # a gene perturbation's own gene and the protein it encodes

REASON_NO_EDGE_IN_FILE = "no edge out in the graph file"
REASON_LOADER_DROPS = "edges the loader drops"
REASON_ENCODER_DOES_NOT_WALK = "edges the encoder does not walk"
REASON_ONLY_INTO_DRUG_NODES = "edges only into drug nodes"
REASON_INTO_WEIGHTLESS_NODES = "edges into nodes without class weight"
REASON_SEED_WITHOUT_WEIGHT = "the seed has no class weight"


def walkable_adjacency(num_nodes: int, walked_edge_source: np.ndarray, walked_edge_target: np.ndarray,
                       node_carries_class_weight: np.ndarray, entry_node_mask: np.ndarray) -> sparse.csr_matrix:
    """The edges a route can take: walked by the bottleneck, both ends with a class weight above 0 in some class, and
    not into an entry node, since an entry node is present only in the perturbation that seeds it."""
    route_can_take_edge = (node_carries_class_weight[walked_edge_source] & node_carries_class_weight[walked_edge_target]
                           & ~entry_node_mask[walked_edge_target])
    return sparse.csr_matrix((np.ones(int(route_can_take_edge.sum()), dtype=np.int8),
                              (walked_edge_source[route_can_take_edge], walked_edge_target[route_can_take_edge])),
                             shape=(num_nodes, num_nodes))


def reachable_nodes(adjacency: sparse.csr_matrix, seed_nodes: list[int], node_carries_class_weight: np.ndarray) -> set[int]:
    """Every node a route from the seeds can reach, the seeds included; a seed without class weight starts no route,
    because the recurrence starts it at its class weight."""
    reached = set()
    for seed_node in seed_nodes:
        if not node_carries_class_weight[seed_node]:
            continue
        reached.update(int(node) for node in breadth_first_order(adjacency, seed_node, directed=True, return_predecessors=False))
    return reached


def reasons_routes_stop(reached: set[int], seed_nodes: list[int], node_carries_class_weight: np.ndarray, entry_node_mask: np.ndarray,
                        file_out_edges: dict[int, list[tuple[int, str]]], loaded_out_edges: dict[int, list[tuple[int, str]]],
                        walked_out_targets: dict[int, set[int]]) -> tuple[list[str], dict[str, Counter]]:
    """Why no route leaves the reached set: the reasons that hold at some reached node, in the order first met, and
    the relations of the edges behind each reason that concerns edges.

    file_out_edges: per node, (target, relation) of each edge leaving it in the graph file, target -1 when the loader
    has no such node; loaded_out_edges: the same for the loaded edges; walked_out_targets: per node, the targets of
    the edges the bottleneck walks from it."""
    if not any(node_carries_class_weight[seed_node] for seed_node in seed_nodes):
        return [REASON_SEED_WITHOUT_WEIGHT], {}
    reasons: list[str] = []
    relations_by_reason: dict[str, Counter] = {}

    def add(reason: str, relations: Counter | None = None) -> None:
        if reason not in reasons:
            reasons.append(reason)
        if relations:
            relations_by_reason.setdefault(reason, Counter()).update(relations)

    for node in sorted(reached):
        file_edges = file_out_edges.get(node, [])
        loaded_edges = loaded_out_edges.get(node, [])
        walked_targets = walked_out_targets.get(node, set())
        if not file_edges and not loaded_edges and not walked_targets:
            add(REASON_NO_EDGE_IN_FILE)
            continue
        loaded_edge_set = set(loaded_edges)
        dropped_at_load = Counter(relation for target, relation in file_edges if (target, relation) not in loaded_edge_set)
        if dropped_at_load:
            add(REASON_LOADER_DROPS, dropped_at_load)
        not_walked = Counter(relation for target, relation in loaded_edges if target not in walked_targets)
        if not_walked:
            add(REASON_ENCODER_DOES_NOT_WALK, not_walked)
        blocked_targets = [target for target in walked_targets if target not in reached]
        if blocked_targets:
            if all(entry_node_mask[target] for target in blocked_targets):
                add(REASON_ONLY_INTO_DRUG_NODES)
            else:
                add(REASON_INTO_WEIGHTLESS_NODES)
    return reasons, relations_by_reason


def edges_out_of_each_node(edge_source: np.ndarray, edge_target: np.ndarray, edge_relation_names: np.ndarray) -> dict[int, list[tuple[int, str]]]:
    edges_by_node: dict[int, list[tuple[int, str]]] = {}
    for source_node, target_node, relation_name in zip(edge_source.tolist(), edge_target.tolist(), edge_relation_names.tolist()):
        edges_by_node.setdefault(int(source_node), []).append((int(target_node), str(relation_name)))
    return edges_by_node


def targets_out_of_each_node(edge_source: np.ndarray, edge_target: np.ndarray) -> dict[int, set[int]]:
    targets_by_node: dict[int, set[int]] = {}
    for source_node, target_node in zip(edge_source.tolist(), edge_target.tolist()):
        targets_by_node.setdefault(int(source_node), set()).add(int(target_node))
    return targets_by_node


def describe_relations(relations: Counter) -> str:
    return ";".join(f"{name}:{count}" for name, count in sorted(relations.items()))


def git_commit() -> str:
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], capture_output=True, text=True, check=True).stdout.strip()
        return commit + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def measure(data, encoder, graph_dir: Path, dead_end_reach: int) -> tuple[pd.DataFrame, dict]:
    """One row per perturbation with its reach and, for a dead end, why its routes stop; and a summary."""
    num_nodes = len(data.node_ids)
    if encoder.cell_class_node_weight is None:
        raise ValueError("the reach is defined under the limiting-step rule, which needs cell-class weights")
    node_carries_class_weight = (encoder.cell_class_node_weight.max(dim=1).values > 0).cpu().numpy()
    entry_node_mask = np.zeros(num_nodes, dtype=bool) if encoder.entry_node_mask is None else encoder.entry_node_mask.cpu().numpy().astype(bool)
    walked_edge_source = encoder.bottleneck_edge_source.cpu().numpy()
    walked_edge_target = encoder.bottleneck_edge_target.cpu().numpy()
    adjacency = walkable_adjacency(num_nodes, walked_edge_source, walked_edge_target, node_carries_class_weight, entry_node_mask)

    file_edges = pd.read_parquet(Path(graph_dir) / "edges.parquet", columns=["source_id", "target_id", "relation_type"])
    file_source = file_edges["source_id"].map(data.node_index)
    file_target = file_edges["target_id"].map(data.node_index).fillna(-1)
    source_is_loaded = file_source.notna().to_numpy()
    file_out_edges = edges_out_of_each_node(file_source[source_is_loaded].astype(int).to_numpy(), file_target[source_is_loaded].astype(int).to_numpy(),
                                            file_edges.loc[source_is_loaded, "relation_type"].to_numpy())
    loaded_relation_names = np.asarray(data.relation_types, dtype=object)[np.asarray(data.edge_relation, dtype=int)]
    loaded_out_edges = edges_out_of_each_node(np.asarray(data.edge_source, dtype=int), np.asarray(data.edge_target, dtype=int), loaded_relation_names)
    walked_out_targets = targets_out_of_each_node(walked_edge_source, walked_edge_target)

    positives = np.asarray(data.outcomes) > 0
    if data.label_mask is not None:
        positives &= np.asarray(data.label_mask, dtype=bool)
    rows = []
    reason_counts: Counter = Counter()
    relation_counts_by_reason: dict[str, Counter] = {}
    for row, perturbation_id in enumerate(data.perturbation_ids):
        seed_nodes = [int(node) for node in data.perturbation_seeds[row]]
        reached = reachable_nodes(adjacency, seed_nodes, node_carries_class_weight)
        is_dead_end = len(reached) <= dead_end_reach
        reasons, relations_by_reason = ([], {})
        if is_dead_end:
            reasons, relations_by_reason = reasons_routes_stop(reached, seed_nodes, node_carries_class_weight, entry_node_mask,
                                                               file_out_edges, loaded_out_edges, walked_out_targets)
            reason_counts.update(reasons)
            for reason, relations in relations_by_reason.items():
                relation_counts_by_reason.setdefault(reason, Counter()).update(relations)
        rows.append({
            "perturbation_id": perturbation_id,
            "perturbation_type": data.perturbation_types[row],
            "seed_ids": ";".join(str(data.node_ids[node]) for node in seed_nodes),
            "seed_node_types": ";".join(str(data.node_types[node]) for node in seed_nodes),
            "reachable_nodes": len(reached),
            "positive_pairs": int(positives[row].sum()),
            "dead_end": is_dead_end,
            "reasons": "; ".join(reasons),
            "reached_node_ids": ";".join(str(data.node_ids[node]) for node in sorted(reached)) if is_dead_end else "",
            "relations_out_of_reached_nodes_in_file": describe_relations(Counter(relation for node in reached for _, relation in file_out_edges.get(node, [])))
                                                      if is_dead_end else "",
            "relations_by_reason": "; ".join(f"{reason}: {describe_relations(relations)}" for reason, relations in relations_by_reason.items()),
        })
    table = pd.DataFrame(rows)
    by_type = {}
    for kind, group in table.groupby("perturbation_type"):
        by_type[kind] = {"perturbations": int(len(group)), "dead_ends": int(group["dead_end"].sum()),
                         "median_reachable_nodes": float(group["reachable_nodes"].median()),
                         "positive_pairs": int(group["positive_pairs"].sum()),
                         "positive_pairs_on_dead_ends": int(group.loc[group["dead_end"], "positive_pairs"].sum())}
    summary = {
        "script": "experiments/measure_perturbation_reach.py",
        "commit": git_commit(),
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "graph_dir": str(graph_dir),
        "nodes": num_nodes,
        "nodes_with_class_weight": int(node_carries_class_weight.sum()),
        "edges_walked": int(len(walked_edge_source)),
        "edges_a_route_can_take": int(adjacency.nnz),
        "dead_end_reach": dead_end_reach,
        "perturbations": int(len(table)),
        "dead_ends": int(table["dead_end"].sum()),
        "by_type": by_type,
        "reasons": dict(reason_counts),
        "relations_by_reason": {reason: dict(relations) for reason, relations in relation_counts_by_reason.items()},
    }
    return table, summary


def main() -> None:
    import torch

    from experiments.run_main_model import build_argument_parser, build_models, resolve_registered_defaults
    from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

    parser = build_argument_parser()
    parser.add_argument("--reach-output", type=Path, required=True, help="CSV, one row per perturbation")
    parser.add_argument("--reach-summary-output", type=Path, required=True, help="JSON summary")
    parser.add_argument("--dead-end-reach", type=int, default=DEFAULT_DEAD_END_REACH,
                        help="a perturbation whose routes reach at most this many nodes is a dead end (default 2: a gene and its protein)")
    arguments = parser.parse_args()
    resolve_registered_defaults(arguments)
    if arguments.encoder != "linear_response" or arguments.class_weight_rule != "limiting_step":
        raise ValueError("the reach is the limiting-step rule's: use --encoder linear_response --class-weight-rule limiting_step")
    torch.manual_seed(arguments.seed)
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only,
                                group_by=arguments.group_by, label_grades=tuple(arguments.label_grades) if arguments.label_grades else None,
                                label_selection=arguments.label_selection, drug_entry=arguments.drug_entry or "nodes",
                                combined_perturbations=arguments.combined_perturbations)
    encoder, _ = build_models(data, arguments, torch.device("cpu"))
    table, summary = measure(data, encoder, arguments.graph_dir, arguments.dead_end_reach)
    arguments.reach_output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.reach_output, index=False)
    arguments.reach_summary_output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in ("perturbations", "dead_ends", "by_type", "reasons")}, indent=2))


if __name__ == "__main__":
    main()
