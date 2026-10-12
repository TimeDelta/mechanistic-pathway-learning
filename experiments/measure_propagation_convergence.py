"""How many propagation steps the linear response needs to settle, which is the number that would replace 8.

The owner, 12 October 2026: "Don't make it a third encoder. Use the convergence instead of 8 steps. 8 steps is
arbitrary, not physiological."

Why this can be done for the linear-response encoder without a new encoder. Two properties it already has: it holds
the perturbation on at every step (the recurrence adds the sustained input each round, h(t+1) = (1 - damping) h(t) +
damping (sum_r gain_r S_r h(t) + u)), and its gains are bounded by the contraction, so the transition is a
contraction and the iteration has a fixed point. Iterating to a tolerance instead of to a fixed count is therefore a
change of stopping rule, not of architecture. Message passing has neither property: it writes the perturbation once
and gives every layer its own weights, so running it to convergence is a different model.

This script measures the stopping rule's cost: the relative change per step, and the step at which it falls below a
tolerance, at several gains. The gain is the thing that decides the answer, so it is swept rather than assumed.

Structural only: every relation takes the same gain, because the gains are learned and no model has been trained on
this table. Nothing here changes a label, a graph, a model or a model input.

Usage:
  PYTHONPATH=. python experiments/measure_propagation_convergence.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REGISTERED_STEPS = 8
REGISTERED_DAMPING = 0.5
REGISTERED_CONTRACTION = 0.9
GAINS = (0.9, 0.7, 0.5, 0.3)
TOLERANCES = (1e-2, 1e-3, 1e-4)


def normalised_adjacency(nodes: pd.DataFrame, edges: pd.DataFrame) -> sparse.csr_matrix:
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    sources = np.array([position_of_node[node] for node in edges.source_id], dtype=np.int64)
    targets = np.array([position_of_node[node] for node in edges.target_id], dtype=np.int64)
    adjacency = sparse.coo_matrix((np.ones(len(sources), dtype=np.float32), (targets, sources)),
                                  shape=(len(nodes), len(nodes))).tocsr()
    in_degree = np.asarray(adjacency.sum(axis=1)).ravel()
    in_degree[in_degree == 0] = 1.0
    return sparse.diags(1.0 / in_degree).tocsr() @ adjacency


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--maximum-steps", type=int, default=400)
    parser.add_argument("--damping", type=float, default=REGISTERED_DAMPING)
    parser.add_argument("--output", type=Path, default=Path("data/processed/tissue_saturation/convergence.json"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id"])
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet", columns=["source_id", "target_id"])
    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet",
                              columns=["perturbation_type", "perturbation_label", "perturbation_nodes"])
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    columns = []
    for row in records[records.perturbation_type == "drug"].drop_duplicates(subset=["perturbation_label"]).itertuples(index=False):
        column = np.zeros(len(nodes), dtype=np.float32)
        for node_id, sign, magnitude in json.loads(row.perturbation_nodes):
            index = position_of_node.get(node_id)
            if index is not None:
                column[index] += float(sign) * float(magnitude)
        if np.abs(column).sum():
            columns.append(column)
    sustained = np.stack(columns, axis=1)
    base_adjacency = normalised_adjacency(nodes, edges)
    print(f"graph {len(nodes)} nodes, {len(edges)} edges; {sustained.shape[1]} drug perturbations; "
          f"damping {arguments.damping}")

    summary = {}
    for gain in GAINS:
        adjacency = gain * base_adjacency
        state = sustained.copy()
        steps_to, previous = {}, None
        reached_at_registered = None
        for step in range(1, arguments.maximum_steps + 1):
            previous = state
            state = (1.0 - arguments.damping) * state + arguments.damping * (adjacency @ state + sustained)
            change = np.linalg.norm(state - previous) / max(np.linalg.norm(state), 1e-12)
            for tolerance in TOLERANCES:
                if tolerance not in steps_to and change < tolerance:
                    steps_to[tolerance] = step
            if step == REGISTERED_STEPS:
                reached_at_registered = state.copy()
            if len(steps_to) == len(TOLERANCES):
                break
        settled = state
        share = (np.linalg.norm(reached_at_registered) / max(np.linalg.norm(settled), 1e-12)
                 if reached_at_registered is not None else float("nan"))
        summary[str(gain)] = {"steps_to_tolerance": {str(tolerance): steps_to.get(tolerance) for tolerance in TOLERANCES},
                              "norm_at_registered_steps_over_settled": float(share)}
        print(f"\ngain {gain}{' (the contraction bound)' if gain == REGISTERED_CONTRACTION else ''}:")
        for tolerance in TOLERANCES:
            print(f"  relative change below {tolerance:g} at step {steps_to.get(tolerance, f'>{arguments.maximum_steps}')}")
        print(f"  the field at the registered {REGISTERED_STEPS} steps is {share:.3f} of the settled field by norm")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(summary, indent=1, sort_keys=True))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
