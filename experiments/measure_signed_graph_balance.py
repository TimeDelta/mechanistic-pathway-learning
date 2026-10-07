"""Frustration of the signed graph, which decides whether the edge signs carry structure.

Why this matters. The linear-response encoder propagates over a signed adjacency whose entry
[destination, source] is the edge sign divided by a normaliser. Switch the node basis by a
diagonal of plus and minus ones, D, and every relation's matrix becomes D S_r D, so the fixed
point (I - sum_r gain_r S_r)^-1 u becomes D (I - sum_r gain_r S'_r)^-1 D u. If some D makes
every entry non-negative, signed propagation equals unsigned propagation up to a per-node sign
on the input and the readout, which a linear head absorbs. The signs then carry no structure
of their own, and a signed encoder cannot beat an unsigned walk for that reason alone.

Such a D exists exactly when the signed graph is balanced: the switched entry is
sign * node_sign[source] * node_sign[target], so every edge imposes
node_sign[source] * node_sign[target] = sign. That is the balance condition on the underlying
undirected signed graph (Harary's theorem), reached here by breadth-first two-colouring per
connected component. The violated-edge count under one breadth-first spanning forest is an
upper bound on the frustration index, not the index itself, which is NP-hard; a count of zero
is exact and proves balance.

The measurement runs on the relation stack the encoder builds, not on the stored edge table: the
metabolic graph stores no negative edge, so a balance measurement taken from the table reports zero
frustration for the trivial reason that every sign is the same.

Usage: PYTHONPATH=. python3 experiments/measure_signed_graph_balance.py --graph-dir data/processed/graph
"""
from __future__ import annotations

import argparse
import json
from collections import deque
from pathlib import Path

import pandas as pd


def measure_frustration(edges: pd.DataFrame) -> dict:
    """Two-colour the nodes so each edge wants node_sign[source] * node_sign[target] == sign."""
    adjacency: dict[str, list[tuple[str, int]]] = {}
    for source, target, sign in zip(edges["source_id"], edges["target_id"], edges["sign"]):
        required_product = 1 if float(sign) >= 0 else -1
        adjacency.setdefault(source, []).append((target, required_product))
        adjacency.setdefault(target, []).append((source, required_product))

    node_sign: dict[str, int] = {}
    components = 0
    for start_node in adjacency:
        if start_node in node_sign:
            continue
        components += 1
        node_sign[start_node] = 1
        frontier = deque([start_node])
        while frontier:
            current = frontier.popleft()
            for neighbour, required_product in adjacency[current]:
                if neighbour not in node_sign:
                    node_sign[neighbour] = node_sign[current] * required_product
                    frontier.append(neighbour)

    frustrated_by_relation: dict[str, int] = {}
    total_by_relation: dict[str, int] = {}
    self_loops_negative = 0
    for source, target, sign, relation in zip(edges["source_id"], edges["target_id"],
                                              edges["sign"], edges["relation_type"]):
        total_by_relation[relation] = total_by_relation.get(relation, 0) + 1
        required_product = 1 if float(sign) >= 0 else -1
        if source == target:
            if required_product == -1:
                self_loops_negative += 1
                frustrated_by_relation[relation] = frustrated_by_relation.get(relation, 0) + 1
            continue
        if node_sign[source] * node_sign[target] != required_product:
            frustrated_by_relation[relation] = frustrated_by_relation.get(relation, 0) + 1

    frustrated_total = sum(frustrated_by_relation.values())
    return {
        "num_edges": int(len(edges)),
        "num_nodes_in_signed_subgraph": len(node_sign),
        "num_connected_components": components,
        "num_negative_edges": int((edges["sign"] < 0).sum()),
        "num_negative_self_loops": self_loops_negative,
        "frustrated_edges_upper_bound": frustrated_total,
        "frustrated_fraction_upper_bound": frustrated_total / max(len(edges), 1),
        "is_balanced": frustrated_total == 0,
        "per_relation": {relation: {"edges": total_by_relation[relation],
                                    "frustrated": frustrated_by_relation.get(relation, 0)}
                         for relation in sorted(total_by_relation)},
    }


def add_derived_depletion_relation(edges: pd.DataFrame) -> pd.DataFrame:
    """The reverse of every substrate_of edge with sign -1, as signed_stacked_adjacency appends it.

    Without this the measurement is vacuous on the metabolic graph, whose stored edge table holds no
    negative edge at all: every negative sign in the model is derived here.
    """
    substrate_edges = edges[edges["relation_type"] == "substrate_of"]
    derived = pd.DataFrame({"source_id": substrate_edges["target_id"].values,
                            "target_id": substrate_edges["source_id"].values,
                            "relation_type": "depletes_substrate",
                            "sign": -1.0})
    return pd.concat([edges[["source_id", "target_id", "relation_type", "sign"]], derived], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", default="data/processed/graph")
    parser.add_argument("--output-json", default=None)
    parser.add_argument("--stored-edges-only", action="store_true",
                        help="measure the stored edge table alone, which carries no negative edge on the "
                             "metabolic graph and so reports a vacuous zero")
    arguments = parser.parse_args()

    edges = pd.read_parquet(Path(arguments.graph_dir) / "edges.parquet")
    if not arguments.stored_edges_only:
        edges = add_derived_depletion_relation(edges)
    result = measure_frustration(edges)
    result["graph_dir"] = arguments.graph_dir

    print(f"graph {arguments.graph_dir}")
    print(f"  edges {result['num_edges']}, negative {result['num_negative_edges']}, "
          f"nodes touched {result['num_nodes_in_signed_subgraph']}, components {result['num_connected_components']}")
    print(f"  frustrated edges (upper bound) {result['frustrated_edges_upper_bound']} "
          f"= {100 * result['frustrated_fraction_upper_bound']:.3f}% of edges")
    print(f"  balanced: {result['is_balanced']}")
    for relation, counts in result["per_relation"].items():
        if counts["frustrated"]:
            print(f"    {relation}: {counts['frustrated']} of {counts['edges']} frustrated")
    if arguments.output_json:
        Path(arguments.output_json).write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
