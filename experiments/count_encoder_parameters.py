"""What each encoder's parameters are spent on, and what the standard reductions would cost.

The owner, 12 October 2026: "there's gotta be a way to parameterize that up-projection so that it doesn't require
that many params. Otherwise, I think I'm gonna drop the MP encoder variant because its just getting too parameter
hungry."

The message-passing encoder's relation transforms are one [node_state_dim, node_state_dim] matrix per relation per
layer, which is where its parameters are. This script counts them against the confirmatory graph and against the
labelled pairs they are fitted on, and prices four reductions, three of which are the standard ones for a relational
graph network (Schlichtkrull et al. 2018 introduced basis and block-diagonal decomposition for exactly this reason).

Counting only: no label, table, graph or model input is changed, and nothing is built.

Usage:
  PYTHONPATH=. python experiments/count_encoder_parameters.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

PERTURBATION_FEATURE_DIM = 3  # relational_message_passing_encoder.PERTURBATION_FEATURE_DIM


def relation_parameter_schemes(num_relations: int, node_state_dim: int, num_layers: int) -> list[dict]:
    """Each scheme: what the relation transforms would cost, and what one relation owns under it."""
    full_per_layer = num_relations * node_state_dim * node_state_dim
    schemes = [
        {"scheme": "registered: a full matrix per relation per layer",
         "total": num_layers * full_per_layer,
         "one_relation": num_layers * node_state_dim * node_state_dim,
         "note": "what the encoder does now"},
        {"scheme": "the same matrices shared across layers",
         "total": full_per_layer,
         "one_relation": node_state_dim * node_state_dim,
         "note": "what convergence needs anyway, since an iterated update cannot have per-round weights"},
    ]
    for num_bases in (2, 4, 8):
        schemes.append({
            "scheme": f"basis decomposition, {num_bases} bases, shared across layers",
            "total": num_bases * node_state_dim * node_state_dim + num_relations * num_bases,
            "one_relation": num_bases,
            "note": "W_r = sum_b a_rb V_b; the bases are shared and a relation owns only its mixing coefficients"})
    for num_blocks in (2, 4, 8):
        schemes.append({
            "scheme": f"block diagonal, {num_blocks} blocks, shared across layers",
            "total": num_relations * num_blocks * (node_state_dim // num_blocks) ** 2,
            "one_relation": num_blocks * (node_state_dim // num_blocks) ** 2,
            "note": "each relation mixes only within blocks of the node state"})
    for rank in (2, 4, 8):
        schemes.append({
            "scheme": f"per-relation low rank {rank}, shared across layers",
            "total": num_relations * 2 * node_state_dim * rank,
            "one_relation": 2 * node_state_dim * rank,
            "note": "W_r = U_r V_r^T, each relation with its own rank-r factors; no relation shares a subspace with "
                    "another, so it assumes only that each transform is low rank"})
    schemes.append({
        "scheme": "diagonal, shared across layers",
        "total": num_relations * node_state_dim,
        "one_relation": node_state_dim,
        "note": "one gain per relation per channel, which is what the linear-response encoder already is"})
    return schemes


def rank_matched_to_evidence(edge_counts: dict[str, int], node_state_dim: int,
                             ranks_by_edges: tuple[tuple[int, int], ...]) -> dict:
    """A per-relation rank chosen by how many edges the relation has, so capacity follows evidence.

    ranks_by_edges is (at least this many edges, rank), largest threshold first."""
    total, per_relation = 0, {}
    for relation, edges in edge_counts.items():
        rank = next(rank for threshold, rank in ranks_by_edges if edges >= threshold)
        per_relation[relation] = {"edges": edges, "rank": rank, "parameters": 2 * node_state_dim * rank}
        total += 2 * node_state_dim * rank
    return {"total": total, "per_relation": per_relation}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path,
                        default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--node-state-dim", type=int, default=32)
    parser.add_argument("--num-layers", type=int, default=3)
    parser.add_argument("--propagation-channels", type=int, default=4)
    parser.add_argument("--descriptor-columns", type=int, default=138)
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/encoder_parameters.json"))
    arguments = parser.parse_args()

    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet", columns=["relation_type"])
    num_relations = edges.relation_type.nunique()
    selection = pd.read_parquet(arguments.label_selection)
    kept_pairs = int(selection.keep.astype(bool).sum())
    state, channels = arguments.node_state_dim, arguments.propagation_channels
    print(f"graph {arguments.graph_dir.name}: {num_relations} relations; node state {state}; "
          f"{arguments.num_layers} message-passing layers; {channels} propagation channels")
    print(f"labelled pairs the encoders are fitted on (kept in the selection): {kept_pairs}")

    message_passing = {
        "relation transforms": arguments.num_layers * num_relations * state * state,
        "self transforms": arguments.num_layers * state * state,
        "layer biases": arguments.num_layers * state,
        "feature projection": arguments.descriptor_columns * state + state,
        "perturbation injection": PERTURBATION_FEATURE_DIM * state + state,
    }
    linear_response = {
        "relation gains": num_relations * channels,
        "expansion to the node state": channels * state + state,
        "feature projection": arguments.descriptor_columns * state + state,
        "perturbation injection": PERTURBATION_FEATURE_DIM * state + state,
    }
    for name, parts in (("message passing", message_passing), ("linear response", linear_response)):
        total = sum(parts.values())
        print(f"\n{name}: {total:,} encoder parameters, {total / kept_pairs:.1f} per labelled pair")
        for part, count in sorted(parts.items(), key=lambda item: -item[1]):
            print(f"    {part}: {count:,} ({count / total:.0%})")

    print(f"\nWhat the relation transforms would cost under each scheme "
          f"(the registered {message_passing['relation transforms']:,} is the first row):")
    schemes = relation_parameter_schemes(num_relations, state, arguments.num_layers)
    print(f"  {'scheme':<58} {'total':>9} {'one relation':>13} {'per pair':>9}")
    for scheme in schemes:
        print(f"  {scheme['scheme']:<58} {scheme['total']:>9,} {scheme['one_relation']:>13,} "
              f"{scheme['total'] / kept_pairs:>9.2f}")
    print()
    for scheme in schemes:
        print(f"  {scheme['scheme']}: {scheme['note']}")

    # A rank matched to the evidence, since the relations differ in size by four orders of magnitude.
    edge_counts = edges.relation_type.value_counts().to_dict()
    tiers = ((10_000, 8), (1_000, 4), (100, 2), (0, 1))
    matched = rank_matched_to_evidence(edge_counts, state, tiers)
    print(f"\nA rank matched to each relation's evidence, with tiers {tiers}:")
    print(f"  total {matched['total']:,} parameters, {matched['total'] / kept_pairs:.2f} per labelled pair")
    for relation, entry in sorted(matched["per_relation"].items(), key=lambda item: -item[1]["edges"]):
        print(f"    {relation:<30} {entry['edges']:>7,} edges -> rank {entry['rank']}, {entry['parameters']:>5,} parameters")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(
        {"relations": num_relations, "kept_pairs": kept_pairs, "message_passing": message_passing,
         "linear_response": linear_response, "schemes": schemes, "rank_matched_to_evidence": matched},
        indent=1, sort_keys=True))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
