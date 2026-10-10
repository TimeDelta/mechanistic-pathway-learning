"""Probe, with nothing trained, what the message-passing encoder would do if its layers were reused for more rounds.

Why. The registered message-passing encoder runs 3 layers, each with its own weights, so a perturbation's difference
field is exactly zero more than 3 edges from its seeds. docs/drug_entry_nodes_measured.md, section 4, gives what that
reaches: about 1 percent of the nodes a gene knockout can reach and about 35 percent of what a drug can, against 96 to
98 percent at the 8 steps of the linear-response encoder. The user asked why message passing does not carry the signal
until every reachable node has been reached, with parameters reused across layers. This script measures the three
things that question turns on before anything is built:

1. the size of the difference field and the nodes it covers when one layer's weights are reused for 3, 8, 18 and 26
   rounds, at the registered initialisation and with the weights halved (whether a reused layer is stable);
2. the time and peak memory of one forward and backward pass at 3 and at 6 layers, batch 16, each in its own process
   (what a round costs);
3. how far the linear response is from its fixed point at its registered number of steps (whether that encoder
   already is the converged case).

The perturbations are a fixed sample of gene knockouts and drugs, drawn by --seed. Labels are not read.

Usage:
  python experiments/probe_message_passing_depth.py
  python experiments/probe_message_passing_depth.py --document-only
"""
from __future__ import annotations

import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.measure_drug_entry_nodes import batches, linear_response_encoder, markdown_table, message_passing_encoder, padded  # noqa: E402
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data  # noqa: E402

MEASUREMENTS_FILE = "message_passing_depth_probe.json"
ROUNDS = (3, 8, 18, 26)
WEIGHT_SCALES = (1.0, 0.5)
COST_LAYERS = (3, 6)
COST_BATCH_SIZE = 16
LINEAR_RESPONSE_STEPS = (8, 16, 32, 64, 128)


def load(arguments):
    return load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster_and_targets", label_selection=arguments.label_selection)


def cost_of_one_pass(arguments) -> None:
    """One forward and backward pass at --cost-layers; prints one JSON line. Run in its own process so that the peak
    memory is this pass's."""
    torch.manual_seed(arguments.seed)
    data = load(arguments)
    memory_after_loading = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    encoder, adjacencies = message_passing_encoder(data, arguments.cost_layers, entry=False)
    positions = np.random.default_rng(arguments.seed).choice(len(data.perturbation_ids), COST_BATCH_SIZE, replace=False)
    node_index, sign_and_magnitude = padded([data.perturbation_seeds[position] for position in positions], [data.perturbation_signs[position] for position in positions],
                                            [data.perturbation_magnitudes[position] for position in positions])
    started = time.time()
    field = encoder.perturbation_difference_field(node_index, sign_and_magnitude, adjacencies)
    forward_seconds = time.time() - started
    field.abs().mean().backward()
    print(json.dumps({"layers": arguments.cost_layers, "batch_size": COST_BATCH_SIZE, "forward_seconds": forward_seconds, "forward_and_backward_seconds": time.time() - started,
                      "peak_memory_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, "memory_after_loading_gb": memory_after_loading}))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3_parkinsonism"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"))
    parser.add_argument("--perturbations-per-kind", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/message_passing_depth_probe.md"))
    parser.add_argument("--document-only", action="store_true", help="rewrite the document from the measurements saved beside the graph")
    parser.add_argument("--cost-layers", type=int, default=0, help="internal: measure one forward and backward pass at this many layers and print it")
    arguments = parser.parse_args()
    if arguments.cost_layers:
        cost_of_one_pass(arguments)
        return
    if arguments.document_only:
        write_document(json.loads((arguments.graph_dir / MEASUREMENTS_FILE).read_text()), arguments)
        return

    torch.manual_seed(arguments.seed)
    generator = np.random.default_rng(arguments.seed)
    data = load(arguments)
    is_drug = np.array(data.perturbation_types) == "drug"
    sample = np.concatenate([generator.choice(np.flatnonzero(~is_drug), arguments.perturbations_per_kind, replace=False),
                             generator.choice(np.flatnonzero(is_drug), arguments.perturbations_per_kind, replace=False)])
    is_drug_row = is_drug[sample]
    seeds = [data.perturbation_seeds[position] for position in sample]
    signs = [data.perturbation_signs[position] for position in sample]
    magnitudes = [data.perturbation_magnitudes[position] for position in sample]
    results: dict = {"graph": str(arguments.graph_dir), "evidence": str(arguments.evidence_dir), "seed": arguments.seed, "graph_nodes": len(data.node_ids),
                     "gene_knockouts": int((~is_drug_row).sum()), "drugs": int(is_drug_row.sum()), "reused_layer": [], "linear_response": [], "cost": []}

    # 1. one layer's weights reused for every round: the layer is drawn once, in an encoder of one layer, and copied
    # into every layer of each deeper encoder, so every row of the table reads the same shared layer
    torch.manual_seed(arguments.seed)
    one_layer, _ = message_passing_encoder(data, 1, entry=False)
    for weight_scale in WEIGHT_SCALES:
        for rounds in ROUNDS:
            encoder, adjacencies = message_passing_encoder(data, rounds, entry=False)
            with torch.no_grad():
                encoder.feature_projection.load_state_dict(one_layer.feature_projection.state_dict())
                encoder.perturbation_injection.load_state_dict(one_layer.perturbation_injection.state_dict())
                for name in ("relation_weight", "self_weight", "layer_bias"):
                    shared = getattr(one_layer, name) * (weight_scale if name != "layer_bias" else 1.0)
                    getattr(encoder, name).copy_(shared.expand_as(getattr(encoder, name)).clone())
            norms, covered = [], []
            with torch.no_grad():
                for node_index, sign_and_magnitude in batches(seeds, signs, magnitudes):
                    field = encoder.perturbation_difference_field(node_index, sign_and_magnitude, adjacencies)
                    norms.extend(field.flatten(1).norm(dim=1).tolist())
                    covered.extend((field.abs().amax(dim=2) > 0).sum(dim=1).tolist())
            norms, covered = np.array(norms), np.array(covered)
            results["reused_layer"].append({
                "weight_scale": weight_scale, "rounds": rounds, "all_finite": bool(np.isfinite(norms).all()),
                "median_field_norm_gene_knockouts": float(np.median(norms[~is_drug_row])), "median_field_norm_drugs": float(np.median(norms[is_drug_row])),
                "median_nodes_with_a_nonzero_field_gene_knockouts": float(np.median(covered[~is_drug_row])),
                "median_nodes_with_a_nonzero_field_drugs": float(np.median(covered[is_drug_row]))})
            print(results["reused_layer"][-1], flush=True)

    # 3. the linear response against its own fixed point
    previous = None
    for steps in LINEAR_RESPONSE_STEPS:
        torch.manual_seed(arguments.seed)  # the same gains and input weights at every number of steps
        encoder = linear_response_encoder(data, steps, entry=False)
        with torch.no_grad():
            response = torch.cat([encoder.response(node_index, sign_and_magnitude) for node_index, sign_and_magnitude in batches(seeds, signs, magnitudes)])
        if previous is not None:
            relative = ((response - previous).flatten(1).norm(dim=1) / response.flatten(1).norm(dim=1).clamp_min(1e-30)).numpy()
            results["linear_response"].append({"steps": steps, "against_steps": steps // 2, "median_relative_change": float(np.median(relative)), "largest_relative_change": float(relative.max())})
            print(results["linear_response"][-1], flush=True)
        previous = response

    # 2. the cost of a pass, each in its own process
    for layers in COST_LAYERS:
        command = [sys.executable, __file__, "--graph-dir", str(arguments.graph_dir), "--evidence-dir", str(arguments.evidence_dir), "--label-selection", str(arguments.label_selection),
                   "--seed", str(arguments.seed), "--cost-layers", str(layers)]
        output = subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1]
        results["cost"].append(json.loads(output))
        print(results["cost"][-1], flush=True)

    (arguments.graph_dir / MEASUREMENTS_FILE).write_text(json.dumps(results, indent=1))
    write_document(results, arguments)


def write_document(results: dict, arguments) -> None:
    first_cost, second_cost = results["cost"]
    layers_between = second_cost["layers"] - first_cost["layers"]
    seconds_per_round = (second_cost["forward_and_backward_seconds"] - first_cost["forward_and_backward_seconds"]) / layers_between
    memory_per_round = (second_cost["peak_memory_gb"] - first_cost["peak_memory_gb"]) / layers_between

    def extrapolated(rounds: int) -> str:
        more = rounds - second_cost["layers"]
        return f"{second_cost['forward_and_backward_seconds'] + more * seconds_per_round:.0f} s and {second_cost['peak_memory_gb'] + more * memory_per_round:.1f} GB at {rounds} rounds"

    lines = [
        "# Message passing with its layers reused for more rounds: probes with nothing trained",
        "",
        "Generated by `experiments/probe_message_passing_depth.py`. Do not edit by hand.",
        "",
        f"Graph `{results['graph']}` ({results['graph_nodes']:,} nodes), perturbations of `{results['evidence']}`: {results['gene_knockouts']} gene knockouts and "
        f"{results['drugs']} drugs seeded on their targets, drawn with seed {results['seed']}. No label is read and no parameter is fitted. The reach of a perturbation by number of "
        "hops, over every perturbation, is in `docs/drug_entry_nodes_measured.md`, section 4.",
        "",
        "## 1. One layer's weights reused for every round",
        "",
        "The first layer's maps and bias copied into every layer, so the encoder applies one shared layer; the weights at the registered initialisation and halved. "
        "The field norm is the Euclidean norm of a perturbation's difference field over all nodes and channels.",
        "",
        *markdown_table(["weights", "rounds", "median field norm, gene knockouts", "drugs", "median nodes with a nonzero field, gene knockouts", "drugs", "all values finite"],
                        [[f"x {row['weight_scale']:g}", row["rounds"], f"{row['median_field_norm_gene_knockouts']:.3g}", f"{row['median_field_norm_drugs']:.3g}",
                          f"{row['median_nodes_with_a_nonzero_field_gene_knockouts']:,.0f}", f"{row['median_nodes_with_a_nonzero_field_drugs']:,.0f}", "yes" if row["all_finite"] else "no"]
                         for row in results["reused_layer"]]),
        "",
        "Read the norm columns down each block: if the size of the field changes by orders of magnitude with the number of rounds, a reused layer with the registered update "
        "is not stable on its own, and if perturbations stopped at different rounds the size of the field would carry the stop round.",
        "",
        "## 2. What a round costs",
        "",
        f"One forward and backward pass of the registered encoder (own weights per layer), batch {first_cost['batch_size']}, each in its own process. Shared weights cost the same "
        "per round, because a round stores the same activations whichever weights it reads.",
        "",
        *markdown_table(["layers", "forward, s", "forward and backward, s", "peak memory, GB", "memory after loading the data, GB"],
                        [[row["layers"], f"{row['forward_seconds']:.1f}", f"{row['forward_and_backward_seconds']:.1f}", f"{row['peak_memory_gb']:.2f}", f"{row['memory_after_loading_gb']:.2f}"]
                         for row in results["cost"]]),
        "",
        f"Per added round: {seconds_per_round:.1f} s and {memory_per_round:.2f} GB. Extrapolated linearly, not measured: {extrapolated(9)}, {extrapolated(18)}, {extrapolated(26)}. "
        "The times were taken on two shared cores with other jobs running, so they compare rounds with each other and are not a benchmark. Seed masking propagates a reference "
        "beside every perturbation and doubles both.",
        "",
        "## 3. The linear response against its own fixed point",
        "",
        "The relative change of the response when the number of steps is doubled, over the same perturbations:",
        "",
        *markdown_table(["steps", "against", "median relative change", "largest"],
                        [[row["steps"], row["against_steps"], f"{row['median_relative_change']:.2e}", f"{row['largest_relative_change']:.2e}"] for row in results["linear_response"]]),
        "",
        "The registered linear-response encoder runs 8 steps.",
    ]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.markdown_output}")


if __name__ == "__main__":
    main()
