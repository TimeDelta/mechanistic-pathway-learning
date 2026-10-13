"""What the limiting-step cell-class rule costs and what it changes, measured on the confirmatory graph.

The owner, 13 October 2026: "min along the path makes WAY more sense to me. do that so it aligns with the physiology".
The rule is in place (linear_response_encoder.CLASS_WEIGHT_RULES, "limiting_step"), and this script is the measurement
that was promised before recommending it, because a running minimum changes what the gains mean and a max-min pass over
the graph is not free. It reports:

  1. the size of the pass: edges in the adjacency it walks, and the rounds the bottleneck needs before it stops growing;
  2. its cost: wall time for the bottleneck alone and for a whole response under each rule, at a realistic batch;
  3. what it restores: the per-node attenuation b against the attenuation the product rule applies, and the ratio of the
     two settled fields, which is the signal the rule gives back;
  4. what it gives up: under the limiting-step rule the class weights leave the propagation, so inside the loop the
     classes differ only in their learned per-channel gains. The script measures that directly, by how closely the
     per-class response patterns agree under each rule and under an all-ones encoder that has only the gains to
     separate its classes, so the loss of class-specific routing is a number and not an assertion;
  5. whether it can be computed once: the bottleneck depends on no parameter, so it could be cached per perturbation.
     The script counts the distinct seed sets and seed nodes and gives the memory each cache would take.

Usage, once the v4 graph is rebuilt with a cell-class weight table that covers every node:
  OMP_NUM_THREADS=4 python experiments/measure_limiting_step_rule.py \
      --graph-dir <the rebuilt v4 graph> --evidence-dir data/processed/evidence_full_v4 \
      --cell-class-weights <its weight table> --markdown-output docs/limiting_step_rule_cost.md
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.models.linear_response_encoder import LinearResponseEncoder

BATCH_SIZE = 32
TIMED_REPEATS = 3
QUANTILES = (0.05, 0.25, 0.5, 0.75, 0.95)


def build_encoder(data, cell_class_weights: torch.Tensor, class_weight_rule: str, fixed_point_tolerance: float | None,
                  max_propagation_steps: int, channels_per_cell_class: int) -> LinearResponseEncoder:
    """The same seed under both rules, so the gains and input weights are identical and only the rule differs."""
    torch.manual_seed(0)
    return LinearResponseEncoder(
        len(data.node_ids), data.relation_types,
        edge_source=torch.as_tensor(data.edge_source), edge_target=torch.as_tensor(data.edge_target),
        edge_relation=torch.as_tensor(data.edge_relation), edge_sign=torch.as_tensor(data.edge_sign),
        node_features=torch.as_tensor(data.structural_node_features()), node_state_dim=32,
        non_propagating_nodes=torch.as_tensor(data.is_currency),
        cell_class_weights=cell_class_weights, channels_per_cell_class=channels_per_cell_class,
        class_weight_rule=class_weight_rule, fixed_point_tolerance=fixed_point_tolerance,
        max_propagation_steps=max_propagation_steps,
    )


def batch_of_perturbations(data, rows: list[int]) -> tuple[torch.Tensor, torch.Tensor]:
    widest = max(len(data.perturbation_seeds[row]) for row in rows)
    node_index = torch.full((len(rows), widest), -1, dtype=torch.long)
    sign_and_magnitude = torch.zeros((len(rows), widest, 2))
    for slot, row in enumerate(rows):
        seeds = data.perturbation_seeds[row]
        node_index[slot, : len(seeds)] = torch.as_tensor(seeds)
        sign_and_magnitude[slot, : len(seeds), 0] = torch.as_tensor(data.perturbation_signs[row], dtype=torch.float32)
        sign_and_magnitude[slot, : len(seeds), 1] = torch.as_tensor(data.perturbation_magnitudes[row], dtype=torch.float32)
    return node_index, sign_and_magnitude


def rounds_until_the_bottleneck_settles(encoder: LinearResponseEncoder, node_index: torch.Tensor, cap: int) -> dict:
    """The round at which the pass stops changing, and the reach it has gained by each round: a bottleneck that is
    still growing at the propagation's step count means the rule's ceiling is set by the reach and not by the graph."""
    reach_by_round, previous = [], None
    for rounds in range(1, cap + 1):
        with torch.no_grad():
            bottleneck = encoder.class_bottleneck(node_index, rounds=rounds)
        reach_by_round.append({"rounds": rounds, "nonzero_fraction": float((bottleneck > 0).to(torch.float32).mean()),
                               "mean_where_nonzero": float(bottleneck[bottleneck > 0].mean()) if bool((bottleneck > 0).any()) else 0.0})
        if previous is not None and torch.equal(bottleneck, previous):
            return {"settled_at": rounds - 1, "reach_by_round": reach_by_round}
        previous = bottleneck
    return {"settled_at": None, "reach_by_round": reach_by_round}


def timed(call, repeats: int = TIMED_REPEATS) -> float:
    call()  # once to pay any one-off allocation
    started = time.perf_counter()
    for _ in range(repeats):
        call()
    return (time.perf_counter() - started) / repeats


def product_rule_attenuation(encoder: LinearResponseEncoder, node_index: torch.Tensor, rounds: int) -> torch.Tensor:
    """[N, B, K]: the product of the class weights along the route the bottleneck chose, computed by the same
    recurrence with the minimum replaced by multiplication, so the two numbers are comparable route by route. This is
    what the product rule charges a response that travels the widest route; the propagation charges it over every
    route at once, which is why the fields are compared separately below."""
    with torch.no_grad():
        weights = encoder.cell_class_node_weight
        batch_size, num_classes = node_index.shape[0], weights.shape[1]
        seeded = torch.zeros((batch_size, encoder.num_graph_nodes + 1), dtype=torch.bool)
        seeded.scatter_(1, torch.where(node_index >= 0, node_index, torch.full_like(node_index, encoder.num_graph_nodes)), True)
        seeded = seeded[:, : encoder.num_graph_nodes].permute(1, 0)[:, :, None]
        along = torch.where(seeded, weights[:, None, :], torch.zeros_like(weights[:, None, :]))
        flat = along.reshape(encoder.num_graph_nodes, batch_size * num_classes)
        node_limit = weights[:, None, :].expand(-1, batch_size, -1).reshape(encoder.num_graph_nodes, batch_size * num_classes).contiguous()
        source, target = encoder.bottleneck_edge_source, encoder.bottleneck_edge_target
        for _ in range(rounds):
            arriving = torch.zeros_like(flat)
            for start in range(0, flat.shape[1], 32):
                arriving[:, start : start + 32] = arriving[:, start : start + 32].index_reduce(
                    0, target, flat[:, start : start + 32][source], "amax", include_self=True)
            updated = torch.maximum(flat, node_limit * arriving)
            if torch.equal(updated, flat):
                break
            flat = updated
        return flat.reshape(encoder.num_graph_nodes, batch_size, num_classes)


def between_class_pattern_agreement(field: torch.Tensor, num_classes: int, channels_per_class: int) -> dict:
    """How much of the per-class structure is the route and how much is the scale. field is [B, N, K * C], class-major.

    For each perturbation the K class blocks are flattened to one vector per class and the correlation is taken between
    every pair: a rule that only rescales one shared pattern gives 1, a rule whose classes carry the response over
    different routes gives less. The magnitude of the correlation is what is reported, because each channel carries the
    input with its own learned weight and the sign of that weight is arbitrary at initialisation - a signed correlation
    returns -1 for two classes that hold the same pattern with opposite input signs, which is a statement about the
    initialisation and not about the routing. (A first version of this script reported the signed value and read -0.99
    under both rules for that reason.)"""
    batch_size, num_nodes, _ = field.shape
    by_class = field.reshape(batch_size, num_nodes, num_classes, channels_per_class).mean(dim=3)  # [B, N, K]
    off_diagonal_mask = ~torch.eye(num_classes, dtype=torch.bool)
    agreements = []
    for row in range(batch_size):
        columns = by_class[row].transpose(0, 1)  # [K, N]
        if float(columns.abs().max()) == 0.0:
            continue
        centred = columns - columns.mean(dim=1, keepdim=True)
        unit = centred / centred.norm(dim=1, keepdim=True).clamp_min(1e-30)
        agreements.append(float((unit @ unit.transpose(0, 1)).abs()[off_diagonal_mask].median()))
    return {"perturbations_measured": len(agreements),
            "median_between_class_agreement": float(np.median(agreements)) if agreements else float("nan"),
            "smallest_between_class_agreement": float(np.min(agreements)) if agreements else float("nan")}


def cache_sizes(data, num_nodes: int, num_classes: int) -> dict:
    """The bottleneck depends on no parameter, so it is the same tensor on every round of training for the same
    perturbation. These are the two caches that could hold it: one entry per perturbation, or one per seed node with
    the entries combined by maximum (b over a set of seeds is the maximum of b over each seed)."""
    seed_sets = {tuple(sorted(seeds)) for seeds in data.perturbation_seeds}
    seed_nodes = {node for seeds in data.perturbation_seeds for node in seeds}
    bytes_per_entry = num_nodes * num_classes * 4
    return {"perturbations": len(data.perturbation_seeds), "distinct_seed_sets": len(seed_sets),
            "distinct_seed_nodes": len(seed_nodes), "megabytes_per_entry_float32": bytes_per_entry / 1e6,
            "gigabytes_for_every_seed_set": len(seed_sets) * bytes_per_entry / 1e9,
            "gigabytes_for_every_seed_node": len(seed_nodes) * bytes_per_entry / 1e9}


def quantiles_of(values: torch.Tensor) -> dict:
    flat = values.reshape(-1).to(torch.float32)
    return {f"q{int(quantile * 100):02d}": float(torch.quantile(flat, quantile)) for quantile in QUANTILES}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # no defaults: the v1 and v2 graphs are retired and the v4 graph is not yet rebuilt, so the caller names the
    # graph, the evidence and the weight table that belong together
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    parser.add_argument("--cell-class-weights", type=Path, required=True)
    parser.add_argument("--channels-per-cell-class", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    parser.add_argument("--fixed-point-tolerance", type=float, default=1e-3)
    parser.add_argument("--max-propagation-steps", type=int, default=60)
    parser.add_argument("--settling-cap", type=int, default=24, help="rounds to try when looking for the round the bottleneck settles at")
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster_and_targets")
    weight_table = pd.read_parquet(arguments.cell_class_weights)
    missing = set(data.node_ids) - set(weight_table.index)
    if missing:
        raise ValueError(f"{len(missing)} graph nodes have no row in {arguments.cell_class_weights}")
    cell_class_weights = torch.as_tensor(weight_table.loc[data.node_ids].to_numpy(dtype=np.float32))

    encoders = {rule: build_encoder(data, cell_class_weights, rule, arguments.fixed_point_tolerance,
                                    arguments.max_propagation_steps, arguments.channels_per_cell_class)
                for rule in ("product", "limiting_step")}
    limiting = encoders["limiting_step"]
    rows = list(range(min(arguments.batch_size, len(data.perturbation_seeds))))
    node_index, sign_and_magnitude = batch_of_perturbations(data, rows)

    record: dict = {
        "graph_dir": str(arguments.graph_dir), "evidence_dir": str(arguments.evidence_dir),
        "cell_class_weights": str(arguments.cell_class_weights),
        "nodes": len(data.node_ids), "edge_rows": int(len(data.edge_source)),
        "adjacency_edges_walked": int(limiting.bottleneck_edge_source.numel()),
        "cell_classes": int(cell_class_weights.shape[1]), "channels": int(limiting.propagation_channels),
        "batch_size": len(rows), "varying_weight_nodes": int((cell_class_weights.min(dim=1).values < cell_class_weights.max(dim=1).values).sum()),
        "weight_quantiles": quantiles_of(cell_class_weights),
    }

    with torch.no_grad():
        product_field = encoders["product"].response(node_index, sign_and_magnitude)
        product_steps = encoders["product"].propagation_steps_taken
        limiting_field = limiting.response(node_index, sign_and_magnitude)
        limiting_steps = limiting.propagation_steps_taken
    record["propagation_steps"] = {"product": product_steps, "limiting_step": limiting_steps}

    record["settling"] = rounds_until_the_bottleneck_settles(limiting, node_index, arguments.settling_cap)

    with torch.no_grad():
        bottleneck = limiting.class_bottleneck(node_index, rounds=limiting_steps)
        along_the_route = product_rule_attenuation(limiting, node_index, rounds=limiting_steps)
    reached = bottleneck > 0
    record["attenuation"] = {
        "nodes_reached_fraction": float(reached.to(torch.float32).mean()),
        "bottleneck_quantiles_where_reached": quantiles_of(bottleneck[reached]),
        "product_along_the_route_quantiles_where_reached": quantiles_of(along_the_route[reached]),
        "median_ratio_bottleneck_over_product_along_the_route": float(torch.median(bottleneck[reached] / along_the_route[reached].clamp_min(1e-30))),
        "product_along_the_route_below_1e-6_fraction": float((along_the_route[reached] < 1e-6).to(torch.float32).mean()),
        "bottleneck_below_1e-6_fraction": float((bottleneck[reached] < 1e-6).to(torch.float32).mean()),
    }

    product_magnitude, limiting_magnitude = product_field.abs(), limiting_field.abs()
    both_nonzero = (product_magnitude > 0) & (limiting_magnitude > 0)
    record["field"] = {
        "product_nonzero_fraction": float((product_magnitude > 0).to(torch.float32).mean()),
        "limiting_step_nonzero_fraction": float((limiting_magnitude > 0).to(torch.float32).mean()),
        "product_median_magnitude": float(product_magnitude[product_magnitude > 0].median()) if bool((product_magnitude > 0).any()) else 0.0,
        "limiting_step_median_magnitude": float(limiting_magnitude[limiting_magnitude > 0].median()) if bool((limiting_magnitude > 0).any()) else 0.0,
        "median_ratio_where_both_nonzero": float((limiting_magnitude[both_nonzero] / product_magnitude[both_nonzero]).median()) if bool(both_nonzero.any()) else float("nan"),
        "sign_agreement_where_both_nonzero": float((torch.sign(product_field[both_nonzero]) == torch.sign(limiting_field[both_nonzero])).to(torch.float32).mean()) if bool(both_nonzero.any()) else float("nan"),
    }

    # the floor for the limiting-step rule: with every weight at 1 nothing attenuates, so what is left between the
    # classes is the learned per-channel gain alone. If the limiting-step rule sits at this floor, its classes differ
    # only in scale and the rule has given up class-specific routing entirely.
    gain_only = build_encoder(data, torch.ones_like(cell_class_weights), "product", arguments.fixed_point_tolerance,
                              arguments.max_propagation_steps, arguments.channels_per_cell_class)
    with torch.no_grad():
        gain_only_field = gain_only.response(node_index, sign_and_magnitude)
    record["routing"] = {
        rule: between_class_pattern_agreement(field, int(cell_class_weights.shape[1]), arguments.channels_per_cell_class)
        for rule, field in (("product", product_field), ("limiting_step", limiting_field), ("gain_only_floor", gain_only_field))
    }

    with torch.no_grad():
        record["seconds"] = {
            "bottleneck_alone": timed(lambda: limiting.class_bottleneck(node_index, rounds=limiting_steps)),
            "response_product": timed(lambda: encoders["product"].response(node_index, sign_and_magnitude)),
            "response_limiting_step": timed(lambda: limiting.response(node_index, sign_and_magnitude)),
        }
    record["seconds"]["limiting_step_over_product"] = record["seconds"]["response_limiting_step"] / record["seconds"]["response_product"]
    record["cache"] = cache_sizes(data, len(data.node_ids), int(cell_class_weights.shape[1]))

    print(json.dumps(record, indent=2))
    if arguments.json_output:
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(json.dumps(record, indent=2) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.markdown_output.write_text(write_markdown(record, arguments))


def write_markdown(record: dict, arguments) -> str:
    settling = record["settling"]
    attenuation, field, seconds, cache = record["attenuation"], record["field"], record["seconds"], record["cache"]
    lines = [
        "# What the limiting-step cell-class rule costs",
        "",
        "Measured by `experiments/measure_limiting_step_rule.py` on "
        f"`{record['graph_dir']}` ({record['nodes']:,} nodes, {record['edge_rows']:,} edge rows, "
        f"{record['adjacency_edges_walked']:,} adjacency entries), with "
        f"`{record['cell_class_weights']}` ({record['cell_classes']} classes, "
        f"{record['varying_weight_nodes']:,} nodes whose weight varies by class) at a batch of {record['batch_size']}.",
        "",
        "## The pass",
        "",
        f"- The propagation settles in {record['propagation_steps']['limiting_step']} rounds at a tolerance of "
        f"{arguments.fixed_point_tolerance}, and the bottleneck is taken over that many rounds.",
        f"- The bottleneck itself stops changing at round "
        f"{settling['settled_at'] if settling['settled_at'] is not None else f'more than {arguments.settling_cap}'}.",
        f"- It reaches {attenuation['nodes_reached_fraction']:.1%} of the node-by-class entries.",
        "",
        "## What it restores",
        "",
        "| quantity | q05 | q25 | median | q75 | q95 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for name, key in (("bottleneck b", "bottleneck_quantiles_where_reached"),
                      ("product of the weights along the same route", "product_along_the_route_quantiles_where_reached")):
        quantiles = attenuation[key]
        lines.append(f"| {name} | " + " | ".join(f"{quantiles[column]:.3g}" for column in ("q05", "q25", "q50", "q75", "q95")) + " |")
    lines += [
        "",
        f"The median ratio of the two is {attenuation['median_ratio_bottleneck_over_product_along_the_route']:.3g}. "
        f"{attenuation['product_along_the_route_below_1e-6_fraction']:.1%} of the reached entries have a product below "
        f"1e-6 against {attenuation['bottleneck_below_1e-6_fraction']:.1%} of the bottlenecks.",
        "",
        f"On the settled fields, the limiting-step rule holds a non-zero value at "
        f"{field['limiting_step_nonzero_fraction']:.1%} of the entries against {field['product_nonzero_fraction']:.1%} "
        f"for the product rule, its median magnitude is {field['limiting_step_median_magnitude']:.3g} against "
        f"{field['product_median_magnitude']:.3g}, and where both are non-zero the median ratio is "
        f"{field['median_ratio_where_both_nonzero']:.3g} with "
        f"{field['sign_agreement_where_both_nonzero']:.1%} sign agreement.",
        "",
        "## What it gives up",
        "",
        "The class weights leave the propagation, so inside the loop the classes differ only in their learned "
        "per-channel gains and the class's expression enters afterwards as a per-node ceiling. The agreement between "
        "the per-class response patterns is the measure of that: a rule that only rescales one shared pattern gives 1. "
        "`gain_only_floor` is the same encoder with every weight set to 1, so it is what the classes agree to when "
        "nothing but the gains separates them.",
        "",
        "| rule | median between-class agreement | smallest |",
        "| --- | --- | --- |",
    ]
    for rule in ("product", "limiting_step", "gain_only_floor"):
        routing = record["routing"][rule]
        lines.append(f"| {rule} | {routing['median_between_class_agreement']:.4f} | {routing['smallest_between_class_agreement']:.4f} |")
    lines += [
        "",
        "## The cost",
        "",
        f"- the bottleneck alone: {seconds['bottleneck_alone']:.3f} s per call",
        f"- a whole response under the product rule: {seconds['response_product']:.3f} s",
        f"- a whole response under the limiting-step rule: {seconds['response_limiting_step']:.3f} s, "
        f"{seconds['limiting_step_over_product']:.2f} times the product rule",
        "",
        "The bottleneck depends on no parameter, only on the graph, the class weights and which nodes are perturbed, so "
        "it is the same tensor on every round of training for the same perturbation and can be computed once. "
        f"The evidence holds {cache['perturbations']:,} perturbations, {cache['distinct_seed_sets']:,} distinct seed "
        f"sets and {cache['distinct_seed_nodes']:,} distinct seed nodes; an entry is "
        f"{cache['megabytes_per_entry_float32']:.2f} MB in float32, so a cache of every seed set is "
        f"{cache['gigabytes_for_every_seed_set']:.1f} GB and one of every seed node "
        f"{cache['gigabytes_for_every_seed_node']:.1f} GB (b over a set of seeds is the maximum of b over each seed, so "
        "the per-seed cache is enough).",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
