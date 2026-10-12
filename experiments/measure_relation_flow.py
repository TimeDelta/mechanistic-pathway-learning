"""How much training evidence passes through each relation, and the rank that evidence can support.

The owner, 12 October 2026: "I was actually thinking that the rank of each relation type should be estimated based
on the number of training examples that pass through the relation type and how much that can support."

That is a better criterion than the edge count of 20.3. An edge is structure; what identifies a parameter is
gradient, and gradient reaches a relation only through labelled examples whose perturbation actually sends signal
along its edges. A relation with 41,638 edges that no training perturbation reaches is as unidentified as one with
three.

What is measured, per relation r:

  flow(p, r)      the message magnitude a perturbation p sends along r's edges at the settled state,
                  sum over edges of |h_source| / in_degree(target), which is the quantity the encoder aggregates
  examples(r)     the labelled training pairs behind that flow: each development perturbation counted once per kept
                  pair it holds, and counted for r only when its flow through r clears a share of its own total
  share(r)        flow(p, r) weighted by each perturbation's kept pairs, normalised over relations

Then the rank. Parameters have to be bought out of one budget, not allocated per relation independently: the same
2,975 pairs flow through every relation at once, so a scheme that gives each relation observations-per-parameter in
isolation spends the same evidence many times. The budget here is the kept training pairs divided by a stated number
of pairs per parameter, split between relations in proportion to share(r), and converted to a rank by
rank = parameters / (2 * node_state_dim) for the factorisation W_r = U_r V_r transposed.

Structural: every relation takes the same gain, since the gains are learned and no model has been trained on this
table. Nothing here changes a label, a graph, a model or a model input.

Usage:
  PYTHONPATH=. python experiments/measure_relation_flow.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

REGISTERED_CONTRACTION = 0.9
REGISTERED_DAMPING = 0.5
SETTLE_STEPS = 60  # past the 45 at which the relative change falls below 1e-4 at the contraction bound
PAIRS_PER_PARAMETER = (5, 10, 20)
TOUCHED = 0.001  # a perturbation counts as reaching a relation when that relation carries this share of its flow


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path,
                        default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v2.json"))
    parser.add_argument("--node-state-dim", type=int, default=32)
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/relation_flow.json"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id"])
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet", columns=["source_id", "target_id", "relation_type"])
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    source = np.array([position_of_node[node] for node in edges.source_id], dtype=np.int64)
    target = np.array([position_of_node[node] for node in edges.target_id], dtype=np.int64)
    relations = sorted(edges.relation_type.unique())
    relation_of_edge = edges.relation_type.to_numpy()

    # The encoder normalises by in-degree within a relation, so the divisor is counted that way here too.
    in_degree = np.zeros((len(relations), len(nodes)), dtype=np.float64)
    for index, relation in enumerate(relations):
        of_relation = relation_of_edge == relation
        np.add.at(in_degree[index], target[of_relation], 1.0)
    in_degree[in_degree == 0] = 1.0

    whole = sparse.coo_matrix((np.ones(len(source), dtype=np.float32), (target, source)),
                              shape=(len(nodes), len(nodes))).tocsr()
    whole_in_degree = np.asarray(whole.sum(axis=1)).ravel()
    whole_in_degree[whole_in_degree == 0] = 1.0
    adjacency = REGISTERED_CONTRACTION * (sparse.diags(1.0 / whole_in_degree).tocsr() @ whole)

    selection = pd.read_parquet(arguments.label_selection)
    kept_pairs_of_perturbation = selection[selection.keep.astype(bool)].perturbation_id.value_counts().to_dict()
    held_out = set(json.loads(arguments.lockbox.read_text())["perturbation_ids"])
    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_nodes"]).drop_duplicates(
        subset=["perturbation_id"])

    seeded, weights, identifiers = [], [], []
    for row in records.itertuples(index=False):
        if row.perturbation_id in held_out:
            continue
        pairs = kept_pairs_of_perturbation.get(row.perturbation_id, 0)
        if not pairs:
            continue
        column = np.zeros(len(nodes), dtype=np.float32)
        for node_id, sign, magnitude in json.loads(row.perturbation_nodes):
            index = position_of_node.get(node_id)
            if index is not None:
                column[index] += float(sign) * float(magnitude)
        if np.abs(column).sum():
            seeded.append(column)
            weights.append(pairs)
            identifiers.append(row.perturbation_id)
    sustained = np.stack(seeded, axis=1)
    pair_weight = np.asarray(weights, dtype=np.float64)
    print(f"development perturbations with a kept pair: {sustained.shape[1]}; kept pairs behind them: "
          f"{int(pair_weight.sum()):,}")

    state = sustained.copy()
    for _ in range(SETTLE_STEPS):
        state = (1.0 - REGISTERED_DAMPING) * state + REGISTERED_DAMPING * (adjacency @ state + sustained)
    magnitude = np.abs(state)  # [nodes, perturbations]

    # flow[relation, perturbation]: what each perturbation sends along that relation's edges
    flow = np.zeros((len(relations), sustained.shape[1]), dtype=np.float64)
    for index, relation in enumerate(relations):
        of_relation = relation_of_edge == relation
        divisor = in_degree[index][target[of_relation]][:, None]
        flow[index] = (magnitude[source[of_relation]] / divisor).sum(axis=0)

    total_per_perturbation = flow.sum(axis=0)
    total_per_perturbation[total_per_perturbation == 0] = 1.0
    share_per_perturbation = flow / total_per_perturbation
    examples = (share_per_perturbation >= TOUCHED) @ pair_weight
    weighted_flow = flow @ pair_weight
    share = weighted_flow / weighted_flow.sum()

    edge_counts = edges.relation_type.value_counts().to_dict()
    order = np.argsort(-share)
    print(f"\n{'relation':<30} {'edges':>8} {'examples':>10} {'flow share':>11}")
    for index in order:
        print(f"{relations[index]:<30} {edge_counts[relations[index]]:>8,} {examples[index]:>10,.0f} "
              f"{share[index]:>10.3%}")

    print("\nRanks the evidence supports, with one budget split by flow share "
          f"(a rank costs 2 x {arguments.node_state_dim} parameters):")
    results = {}
    for pairs_per_parameter in PAIRS_PER_PARAMETER:
        budget = pair_weight.sum() / pairs_per_parameter
        ranks = {}
        for index, relation in enumerate(relations):
            allowed = share[index] * budget / (2 * arguments.node_state_dim)
            ranks[relation] = int(np.clip(np.floor(allowed), 0, arguments.node_state_dim))
        spent = sum(2 * arguments.node_state_dim * rank for rank in ranks.values())
        results[str(pairs_per_parameter)] = {"budget_parameters": int(budget), "ranks": ranks, "spent": spent}
        above_zero = {relation: rank for relation, rank in ranks.items() if rank}
        print(f"\n  {pairs_per_parameter} pairs per parameter: budget {budget:,.0f} parameters, spent {spent:,}")
        print(f"    relations with rank 1 or more: {len(above_zero)} of {len(relations)}")
        for relation, rank in sorted(above_zero.items(), key=lambda item: -item[1]):
            print(f"      {relation:<30} rank {rank}")
        if len(above_zero) < len(relations):
            print(f"    rank 0 (the evidence supports no free transform; a diagonal gain is the floor): "
                  f"{sorted(set(relations) - set(above_zero))}")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(
        {"relations": relations, "edges": edge_counts,
         "examples": {relation: float(examples[index]) for index, relation in enumerate(relations)},
         "flow_share": {relation: float(share[index]) for index, relation in enumerate(relations)},
         "ranks": results}, indent=1, sort_keys=True))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
