"""Measure what drug entry nodes change before anything is trained on them (docs/drug_entry_nodes.md).

Nothing here fits a parameter or reads an outcome to choose anything. Every encoder is at its random start with a fixed
seed, and two encoders that are compared share every parameter they have in common, so a difference between them comes
from the graph and not from the draw. Outcomes are read once, to count kept positives per drug.

What is measured:

  1. the loader: the variant read with drug_entry "targets" against the source graph, field by field, and what
     drug_entry "nodes" leaves unchanged (leakage groups, outcomes, label mask, degrees);
  2. the drug targets as a graph structure: how many study drugs share a target node;
  3. presence: the field of a gene knockout with every drug node always present and its mechanism edges averaged like
     any relation, against the field with a drug node present only where it is seeded; both against the source graph;
  4. the hop: how many nodes a drug's signal reaches from its targets and from its own node, by breadth-first search
     along the edges each encoder propagates on;
  5. the linear-response encoder's two arms: the response of a drug seeded on its node against its response seeded on
     its targets and scaled by the entry gain;
  6. the carriers' sequestration edges: what they change in a drug's own field under each encoder, and how the drugs
     with a known carrier differ from the others in kept positives;
  7. the drug classes of the study by ATC group, for the pharmacokinetic question of docs/drug_entry_nodes.md.

Writes docs/drug_entry_nodes_measured.md and a JSON beside the variant graph. Reads no lockbox: the label selection
and the evidence table are read whole, as the trainer reads them before any split.

Usage:
  OMP_NUM_THREADS=2 PYTHONPATH=. python experiments/measure_drug_entry_nodes.py
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.stats import mannwhitneyu

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.graph.drug_entry_nodes import DRUG_MECHANISM_RELATION, ENTRY_RELATIONS, SEQUESTERS_RELATION
from mechanistic_pathway_learning.models.linear_response_encoder import DEPLETES_SUBSTRATE_RELATION, SUBSTRATE_RELATION, LinearResponseEncoder
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

OUTPUT_DOCUMENT = Path("docs/drug_entry_nodes_measured.md")
MEASUREMENTS_FILE = "drug_entry_measurements.json"
NODE_STATE_DIM = 16
TOLERANCE = 1e-6
ATC_GROUP_NAMES = {"N01": "anaesthetics", "N02": "analgesics", "N03": "antiepileptics", "N04": "anti-Parkinson drugs", "N05": "psycholeptics",
                   "N06": "psychoanaleptics", "N07": "other nervous system drugs"}


STARTED = time.time()


def progress(message: str) -> None:
    print(f"[{time.time() - STARTED:6.0f} s] {message}", file=sys.stderr, flush=True)


def same_value(first, second) -> bool:
    if first is None or second is None:
        return first is None and second is None
    if isinstance(first, np.ndarray):
        return first.shape == second.shape and bool(np.array_equal(first, second, equal_nan=first.dtype.kind == "f"))
    if isinstance(first, list) and first and isinstance(first[0], np.ndarray):
        return len(first) == len(second) and all(np.array_equal(one, other) for one, other in zip(first, second))
    return first == second


def padded(seeds: list[np.ndarray], signs: list[np.ndarray], magnitudes: list[np.ndarray]) -> tuple[torch.Tensor, torch.Tensor]:
    width = max(1, max(len(one) for one in seeds))
    node_index = torch.full((len(seeds), width), -1, dtype=torch.long)
    sign_and_magnitude = torch.zeros((len(seeds), width, 2))
    for row, (one_seeds, one_signs, one_magnitudes) in enumerate(zip(seeds, signs, magnitudes)):
        node_index[row, : len(one_seeds)] = torch.as_tensor(one_seeds, dtype=torch.long)
        sign_and_magnitude[row, : len(one_seeds), 0] = torch.as_tensor(one_signs, dtype=torch.float32)
        sign_and_magnitude[row, : len(one_seeds), 1] = torch.as_tensor(one_magnitudes, dtype=torch.float32)
    return node_index, sign_and_magnitude


def edge_tensors(data) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    return (torch.as_tensor(data.edge_source, dtype=torch.long), torch.as_tensor(data.edge_target, dtype=torch.long),
            torch.as_tensor(data.edge_relation, dtype=torch.long), torch.as_tensor(data.edge_sign, dtype=torch.float32))


def without_relation(data, relation: str):
    """A copy of data without the edges of one relation (the relation list is unchanged)."""
    keep = np.asarray(data.edge_relation) != data.relation_types.index(relation)
    return dataclasses.replace(data, edge_source=np.asarray(data.edge_source)[keep], edge_target=np.asarray(data.edge_target)[keep],
                               edge_relation=np.asarray(data.edge_relation)[keep], edge_sign=np.asarray(data.edge_sign)[keep],
                               edge_weight=None if data.edge_weight is None else np.asarray(data.edge_weight)[keep])


def message_passing_encoder(data, num_layers: int, entry: bool, before_first_layer: bool = False,
                            sequestration_carries_change: bool = False) -> tuple[RelationalMessagePassingEncoder, list]:
    source, target, relation, sign = edge_tensors(data)
    keywords = {}
    if sequestration_carries_change:
        sequestration_index = data.relation_types.index(SEQUESTERS_RELATION)
        is_sequestration = relation == sequestration_index
        keywords.update(change_only_edge_index=torch.stack([source[is_sequestration], target[is_sequestration]]), change_only_relation_index=sequestration_index)
    if entry:
        is_entry = np.isin(np.asarray(data.edge_relation), [data.relation_types.index(name) for name in ENTRY_RELATIONS])
        keywords |= {"entry_before_first_layer": before_first_layer, "entry_node_mask": torch.as_tensor(data.entry_node_mask), "entry_edge_index": torch.stack([source[is_entry], target[is_entry]]),
                    "entry_edge_sign": sign[is_entry], "entry_edge_weight": torch.as_tensor(np.asarray(data.edge_weight)[is_entry], dtype=torch.float32),
                    "entry_relation_indices": tuple(data.relation_types.index(name) for name in ENTRY_RELATIONS)}
    encoder = RelationalMessagePassingEncoder(len(data.node_ids), len(data.relation_types), NODE_STATE_DIM, num_layers,
                                              node_features=torch.as_tensor(data.structural_node_features()), **keywords)
    adjacencies = RelationalMessagePassingEncoder.build_relation_adjacencies(torch.stack([source, target]), relation, len(data.node_ids), len(data.relation_types))
    return encoder, adjacencies


def share_message_passing_parameters(reference: RelationalMessagePassingEncoder, reference_data, other: RelationalMessagePassingEncoder, other_data) -> None:
    """Copy into other every parameter it shares with reference: the feature map by feature column, the relation
    maps by relation name, the rest whole. Layer biases are drawn once and shared, nonzero, because a zero bias would
    hide the state an always-present node gets from the bias alone."""
    reference_types, other_types = sorted(set(reference_data.node_types.tolist())), sorted(set(other_data.node_types.tolist()))
    num_other_columns = other.feature_projection.weight.shape[1]
    column_of_reference_type = {name: index for index, name in enumerate(reference_types)}
    # structural columns: one per node type, then the columns every graph has in the same order
    reference_columns = [column_of_reference_type[name] for name in other_types] + list(range(len(reference_types), reference.feature_projection.weight.shape[1]))
    if len(reference_columns) != num_other_columns:
        raise ValueError("the two graphs differ in more than their node types")
    with torch.no_grad():
        other.feature_projection.weight.copy_(reference.feature_projection.weight[:, reference_columns])
        other.feature_projection.bias.copy_(reference.feature_projection.bias)
        other.perturbation_injection.load_state_dict(reference.perturbation_injection.state_dict())
        shared_layers = min(other.relation_weight.shape[0], reference.relation_weight.shape[0])  # an encoder with one more layer keeps its own last one
        for index, name in enumerate(other_data.relation_types):
            other.relation_weight[:shared_layers, index] = reference.relation_weight[:shared_layers, reference_data.relation_types.index(name)]
        other.self_weight[:shared_layers] = reference.self_weight[:shared_layers]
        other.layer_bias[:shared_layers] = reference.layer_bias[:shared_layers]
        if other.entry_unsigned_weight is not None and reference.entry_unsigned_weight is not None:
            layers_in_common = min(other.entry_unsigned_weight.shape[0], reference.entry_unsigned_weight.shape[0])
            other.entry_unsigned_weight[:layers_in_common] = reference.entry_unsigned_weight[:layers_in_common]
            other.entry_signed_weight[:layers_in_common] = reference.entry_signed_weight[:layers_in_common]


def linear_response_encoder(data, num_steps: int, entry: bool) -> LinearResponseEncoder:
    keywords = {}
    if entry:
        keywords = {"entry_node_mask": torch.as_tensor(data.entry_node_mask), "entry_relation_names": ENTRY_RELATIONS, "edge_weight": torch.as_tensor(data.edge_weight)}
    return LinearResponseEncoder(len(data.node_ids), data.relation_types, *edge_tensors(data), torch.as_tensor(data.structural_node_features()), NODE_STATE_DIM,
                                 non_propagating_nodes=torch.as_tensor(data.is_currency), num_propagation_steps=num_steps, **keywords)


def share_linear_response_parameters(reference: LinearResponseEncoder, other: LinearResponseEncoder) -> None:
    with torch.no_grad():
        for index, name in enumerate(other.relation_names):
            other.gain_logit[index] = reference.gain_logit[reference.relation_names.index(name)]
        other.input_weight.copy_(reference.input_weight)
        if other.entry_gain_logit is not None and reference.entry_gain_logit is not None:
            other.entry_gain_logit.copy_(reference.entry_gain_logit)


def batches(seeds, signs, magnitudes, batch_size: int = 4):
    """Padded (node index, sign and magnitude) batches. The fields are compared batch by batch and dropped: one field
    of a few hundred perturbations over 50,000 nodes is more than a gigabyte, and three of them ended an earlier run of
    this script at the container's memory limit."""
    for start in range(0, len(seeds), batch_size):
        yield padded(seeds[start : start + batch_size], signs[start : start + batch_size], magnitudes[start : start + batch_size])


class RowDifferences:
    """Per perturbation, how far two fields are apart, accumulated over batches."""

    def __init__(self) -> None:
        self.largest_difference, self.scale, self.nodes_changed = [], [], []

    def add(self, reference: torch.Tensor, other: torch.Tensor) -> None:
        difference = (reference - other).abs().amax(dim=2)  # [rows, nodes]
        self.largest_difference.append(difference.amax(dim=1))
        self.scale.append(reference.abs().amax(dim=(1, 2)))
        self.nodes_changed.append((difference > TOLERANCE).sum(dim=1))

    def summary(self, rows: np.ndarray | None = None) -> dict:
        """Over every perturbation, or over those where rows is True."""
        largest, scale, nodes_changed = torch.cat(self.largest_difference), torch.cat(self.scale), torch.cat(self.nodes_changed)
        if rows is not None:
            keep = torch.as_tensor(rows, dtype=torch.bool)
            largest, scale, nodes_changed = largest[keep], scale[keep], nodes_changed[keep]
        changed = largest > TOLERANCE * scale.clamp_min(1.0)
        relative = largest / scale.clamp_min(1e-30)
        return {"rows": int(len(largest)), "rows_changed": int(changed.sum()), "largest_absolute_difference": float(largest.max()),
                "median_relative_difference_of_changed_rows": float(relative[changed].median()) if bool(changed.any()) else 0.0,
                "median_nodes_changed_per_changed_row": float(nodes_changed[changed].float().median()) if bool(changed.any()) else 0.0}


def reach_counts(num_nodes: int, source: np.ndarray, target: np.ndarray, seed_lists: list[np.ndarray], max_hops: int) -> np.ndarray:
    """[len(seed_lists), max_hops + 1]: nodes within k directed hops of the seeds, the seeds included, for k = 0..max_hops."""
    adjacency = sparse.csr_matrix((np.ones(len(source), dtype=np.float32), (target, source)), shape=(num_nodes, num_nodes))
    counts = np.zeros((len(seed_lists), max_hops + 1), dtype=int)
    for row, seeds in enumerate(seed_lists):
        reached = np.zeros(num_nodes, dtype=bool)
        reached[np.asarray(seeds, dtype=int)] = True
        counts[row, 0] = int(reached.sum())
        for hop in range(1, max_hops + 1):
            reached |= adjacency.dot(reached.astype(np.float32)) > 0
            counts[row, hop] = int(reached.sum())
    return counts


def sequestration_as_change(as_nodes, arguments) -> dict:
    """What a sequestration edge does when it carries the carrier's change from rest and not the carrier's state.

    Message passing: the drugs whose own field differs from the field on the graph without the edges (the edge is meant
    to be silent while the carrier is at rest). Linear response, which carries the change by construction: each drug
    with a known carrier is simulated alone, with its carriers lowered beside it and the carriers lowered alone; the
    part of the response that passes through the drug node is both - drug - carriers alone, and its size against the
    drug's own response is the factor by which one unit less carrier scales the drug's effect. Nothing is trained, so
    the factor's sign is the graph's and its size is the untrained gains'."""
    torch.manual_seed(arguments.seed)
    layers = arguments.message_passing_layers
    is_drug = np.array(as_nodes.perturbation_types) == "drug"
    drug_positions = np.flatnonzero(is_drug)
    seeds = [as_nodes.perturbation_seeds[position] for position in drug_positions]
    signs = [as_nodes.perturbation_signs[position] for position in drug_positions]
    magnitudes = [as_nodes.perturbation_magnitudes[position] for position in drug_positions]
    sequestration = np.asarray(as_nodes.edge_relation) == as_nodes.relation_types.index(SEQUESTERS_RELATION)
    carriers_of_node: dict[int, list[int]] = {}
    for carrier, drug_node in zip(np.asarray(as_nodes.edge_source)[sequestration].tolist(), np.asarray(as_nodes.edge_target)[sequestration].tolist()):
        carriers_of_node.setdefault(drug_node, []).append(carrier)
    has_carrier = np.array([int(seed_list[0]) in carriers_of_node for seed_list in seeds])

    without_carriage = without_relation(as_nodes, SEQUESTERS_RELATION)
    with_change, with_change_adjacencies = message_passing_encoder(as_nodes, layers, entry=True, sequestration_carries_change=True)
    without_edges, without_edges_adjacencies = message_passing_encoder(without_carriage, layers, entry=True)
    share_message_passing_parameters(with_change, as_nodes, without_edges, without_carriage)
    differences = []
    with torch.no_grad():
        for node_index, sign_and_magnitude in batches(seeds, signs, magnitudes):
            differences.append((with_change.perturbation_difference_field(node_index, sign_and_magnitude, with_change_adjacencies)
                                - without_edges.perturbation_difference_field(node_index, sign_and_magnitude, without_edges_adjacencies)).abs().amax(dim=(1, 2)))
    difference = torch.cat(differences).numpy()
    block: dict = {"message_passing": {"drugs_with_a_carrier_whose_field_changes": int((difference[has_carrier] > TOLERANCE).sum()),
                                        "drugs_without_a_carrier_whose_field_changes": int((difference[~has_carrier] > TOLERANCE).sum()),
                                        "largest_absolute_difference": float(difference.max())}}
    progress("sequestration as change: message passing measured")

    rows = []
    for steps in (arguments.propagation_steps, 4 * arguments.propagation_steps):
        encoder = linear_response_encoder(as_nodes, steps, entry=True)
        factors, cosines = [], []
        with torch.no_grad():
            for seed_list in (seed_list for seed_list, known in zip(seeds, has_carrier) if known):
                drug_node = int(seed_list[0])
                carriers = carriers_of_node[drug_node]
                width = 1 + len(carriers)
                node_index = torch.full((3, width), -1, dtype=torch.long)
                sign_and_magnitude = torch.zeros((3, width, 2))
                node_index[0, 0], node_index[1, 0] = drug_node, drug_node
                sign_and_magnitude[0, 0], sign_and_magnitude[1, 0] = torch.tensor([1.0, 1.0]), torch.tensor([1.0, 1.0])
                for column, carrier in enumerate(carriers):  # the carriers lowered: beside the drug in row 1, alone in row 2
                    node_index[1, 1 + column], node_index[2, column] = carrier, carrier
                    sign_and_magnitude[1, 1 + column], sign_and_magnitude[2, column] = torch.tensor([-1.0, 1.0]), torch.tensor([-1.0, 1.0])
                response = encoder.response(node_index, sign_and_magnitude)  # [3, nodes, channels]
                drug, through_the_drug = response[0], response[1] - response[0] - response[2]
                per_channel_factor = (through_the_drug * drug).sum(dim=0) / (drug * drug).sum(dim=0).clamp_min(1e-30)
                per_channel_cosine = (through_the_drug * drug).sum(dim=0) / (through_the_drug.norm(dim=0) * drug.norm(dim=0)).clamp_min(1e-30)
                factors.append(per_channel_factor.numpy())
                cosines.append(per_channel_cosine.numpy())
        factors, cosines = np.array(factors), np.array(cosines)
        has_response = np.abs(factors).max(axis=1) > 0  # a drug whose mechanism signs are all 0 has no response to scale
        rows.append({"steps": steps, "drugs": int(len(factors)), "drugs_with_a_response": int(has_response.sum()),
                     "drugs_whose_effect_rises_in_every_channel": int((factors[has_response] > 0).all(axis=1).sum()),
                     "median_factor": float(np.median(factors[has_response])), "smallest_factor": float(factors[has_response].min()), "largest_factor": float(factors[has_response].max()),
                     "median_cosine_with_the_drug_s_own_response": float(np.median(cosines[has_response]))})
        progress(f"sequestration as change: linear response at {steps} steps measured")
    block["linear_response_with_the_carriers_lowered"] = rows
    return block


def sequestration_as_change_lines(block: dict | None) -> list[str]:
    if not block:
        return []
    mp = block["message_passing"]
    lines = [
        "### The edge carrying the carrier's change, not its state",
        "",
        "`--sequestration-carries change` (the trainer's default on a graph with the relation) sends along a sequestration edge how far the carrier has moved from its "
        "unperturbed state. In message passing, against the graph without the edges:",
        "",
        f"- drugs with a carrier whose own field changes: **{mp['drugs_with_a_carrier_whose_field_changes']}**; drugs without: {mp['drugs_without_a_carrier_whose_field_changes']} "
        f"(largest absolute difference {mp['largest_absolute_difference']:.1e}).",
        "",
        "The linear response is a deviation from rest, so it carries the change by construction. Each drug with a known carrier simulated alone, with its carriers lowered "
        "by one unit beside it, and the carriers lowered alone; the part that passes through the drug node is the second minus the first minus the third, and the factor is "
        "its size along the drug's own response, per channel. Nothing is trained: the sign of the factor is the graph's and its size is the untrained gains'.",
        "",
        *markdown_table(["steps", "drugs with a carrier", "with a nonzero response", "effect rises in every channel", "median factor", "smallest", "largest", "median cosine with the drug's own response"],
                        [[row["steps"], row["drugs"], row["drugs_with_a_response"], row["drugs_whose_effect_rises_in_every_channel"], f"{row['median_factor']:.3f}",
                          f"{row['smallest_factor']:.3f}", f"{row['largest_factor']:.3f}", f"{row['median_cosine_with_the_drug_s_own_response']:.4f}"]
                         for row in block["linear_response_with_the_carriers_lowered"]]),
        "",
    ]
    return lines


SATURATION_HOPS = 60  # the breadth-first search stops growing well before this on the confirmatory graph; the block records whether it did
LINEAR_RESPONSE_STEPS_FOR_REACH = 8  # the linear-response default, for the share a search of that many hops reaches


def gene_knockout_reach(on_source, layers: int) -> dict:
    """Nodes within the registered number of hops of a gene knockout's seeds, beside the same count for the drugs seeded
    on their targets: the depth a drug node is to be compared with is the one a gene knockout already has, since both
    then sit one edge before the protein."""
    is_drug = np.array(on_source.perturbation_types) == "drug"
    source, target = np.asarray(on_source.edge_source), np.asarray(on_source.edge_target)
    block = {"hops": layers, "largest_hop_searched": SATURATION_HOPS, "linear_response_steps": LINEAR_RESPONSE_STEPS_FOR_REACH}
    for name, positions in (("gene_knockouts", np.flatnonzero(~is_drug)), ("drugs_seeded_on_their_targets", np.flatnonzero(is_drug))):
        counts = reach_counts(len(on_source.node_ids), source, target, [on_source.perturbation_seeds[position] for position in positions], SATURATION_HOPS)
        seed_prefixes = sorted({on_source.node_ids[int(seed)].split(":")[0] for position in positions for seed in on_source.perturbation_seeds[position]})
        final = counts[:, -1]
        hops_to_everything = np.array([int(np.argmax(row == row[-1])) for row in counts])
        hops_to_most = np.array([int(np.argmax(row >= 0.95 * row[-1])) for row in counts])
        block[name] = {"perturbations": int(len(positions)), "seed_node_kinds": seed_prefixes,
                       "median_nodes_by_hops": [float(value) for value in np.median(counts[:, : layers + 2], axis=0)],
                       "quartiles_at_the_registered_depth": [float(value) for value in np.percentile(counts[:, layers], [25, 75])],
                       "still_growing_at_the_largest_hop": int((counts[:, -1] > counts[:, -2]).sum()),
                       "median_nodes_reachable": float(np.median(final)),
                       "median_hops_to_everything_reachable": float(np.median(hops_to_everything)), "largest_hops_to_everything_reachable": int(hops_to_everything.max()),
                       "median_hops_to_95_percent": float(np.median(hops_to_most)),
                       "median_share_reached_at_the_registered_depth": float(np.median(counts[:, layers] / final)),
                       "median_share_reached_at_the_linear_response_steps": float(np.median(counts[:, LINEAR_RESPONSE_STEPS_FOR_REACH] / final))}
    return block


def linear_response_edges(data) -> tuple[np.ndarray, np.ndarray]:
    """The edges the linear-response encoder propagates on: the stored ones and a reverse edge per substrate edge
    (depletes_substrate), without the edges that leave a currency metabolite."""
    source, target, relation = np.asarray(data.edge_source), np.asarray(data.edge_target), np.asarray(data.edge_relation)
    is_substrate = relation == data.relation_types.index(SUBSTRATE_RELATION)
    all_source, all_target = np.concatenate([source, target[is_substrate]]), np.concatenate([target, source[is_substrate]])
    keep = ~np.asarray(data.is_currency)[all_source]
    return all_source[keep], all_target[keep]


def markdown_table(header: list[str], rows: list[list]) -> list[str]:
    return ["| " + " | ".join(header) + " |", "| " + " | ".join("---" for _ in header) + " |"] + ["| " + " | ".join(str(cell) for cell in row) + " |" for row in rows]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source-graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3_parkinsonism"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"))
    parser.add_argument("--carrier-table", type=Path, default=Path("configs/drug_plasma_carriers.csv"))
    parser.add_argument("--onsides-bridge", type=Path, default=Path("data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json"))
    parser.add_argument("--num-gene-knockouts", type=int, default=200, help="gene perturbations sampled for the presence measurement")
    parser.add_argument("--message-passing-layers", type=int, default=3, help="the confirmatory configuration's depth")
    parser.add_argument("--propagation-steps", type=int, default=8, help="the linear-response default")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--document-only", action="store_true", help="rewrite the document from the measurements saved beside the graph, measuring nothing")
    parser.add_argument("--sequestration-as-change-only", action="store_true",
                        help="measure only what the sequestration edges do when they carry the carrier's change, add it to the saved measurements and rewrite the document")
    parser.add_argument("--gene-knockout-reach-only", action="store_true",
                        help="measure only the reach of the gene knockouts (a breadth-first search, no encoder), add it to the saved measurements and rewrite the document")
    arguments = parser.parse_args()
    if arguments.document_only:
        write_document(json.loads((arguments.graph_dir / MEASUREMENTS_FILE).read_text()), arguments)
        return
    if arguments.sequestration_as_change_only:
        saved = json.loads((arguments.graph_dir / MEASUREMENTS_FILE).read_text())
        as_nodes = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, drug_entry="nodes", group_by="disease_cluster_and_targets", label_selection=arguments.label_selection)
        saved["carriers"]["sequestration_as_change"] = sequestration_as_change(as_nodes, arguments)
        (arguments.graph_dir / MEASUREMENTS_FILE).write_text(json.dumps(saved, indent=1))
        write_document(saved, arguments)
        return
    if arguments.gene_knockout_reach_only:
        saved = json.loads((arguments.graph_dir / MEASUREMENTS_FILE).read_text())
        on_source = load_experiment_data(arguments.source_graph_dir, arguments.evidence_dir, group_by="disease_cluster_and_targets", label_selection=arguments.label_selection)
        saved["reach"]["gene_knockouts_beside_drugs"] = gene_knockout_reach(on_source, arguments.message_passing_layers)
        (arguments.graph_dir / MEASUREMENTS_FILE).write_text(json.dumps(saved, indent=1))
        write_document(saved, arguments)
        return
    torch.manual_seed(arguments.seed)
    generator = np.random.default_rng(arguments.seed)

    loader_keywords = {"group_by": "disease_cluster_and_targets", "label_selection": arguments.label_selection}
    on_source = load_experiment_data(arguments.source_graph_dir, arguments.evidence_dir, **loader_keywords)
    as_targets = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, drug_entry="targets", **loader_keywords)
    as_nodes = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, drug_entry="nodes", **loader_keywords)
    build_summary = json.loads((arguments.graph_dir / "drug_entry_summary.json").read_text())
    results: dict = {"graph": str(arguments.graph_dir), "source_graph": str(arguments.source_graph_dir), "evidence": str(arguments.evidence_dir),
                     "label_selection": str(arguments.label_selection), "seed": arguments.seed, "build": {key: build_summary[key] for key in (
                         "drug_nodes_added", "mechanism_edges_added", "mechanism_edges_by_sign", "drugs_given_a_sequestration_edge", "sequestration_edges_added",
                         "sequestration_edges_by_carrier_gene", "drugs_attached_to_a_graph_compound", "nodes", "edges", "nodes_sha256", "edges_sha256")}}

    progress("loaded")
    # 1. the loader
    differing_fields = [field.name for field in dataclasses.fields(on_source) if not same_value(getattr(on_source, field.name), getattr(as_targets, field.name))]
    is_drug = np.array(as_nodes.perturbation_types) == "drug"
    kept_positive = (on_source.outcomes > 0) & np.asarray(on_source.label_mask, dtype=bool)
    results["loader"] = {
        "fields_compared": len(dataclasses.fields(on_source)), "fields_differing_under_targets": differing_fields,
        "structural_features_identical_under_targets": bool(np.array_equal(on_source.structural_node_features(), as_targets.structural_node_features())),
        "perturbations": len(on_source.perturbation_ids), "drug_perturbations": int(is_drug.sum()), "gene_perturbations": int((~is_drug).sum()),
        "identical_under_nodes": {name: same_value(getattr(on_source, name), getattr(as_nodes, name)) for name in ("perturbation_ids", "group_ids", "outcomes", "weights", "label_mask")},
        "degrees_for_strata_identical_under_nodes": bool(np.array_equal(on_source.perturbation_degrees_for_strata, as_nodes.perturbation_degrees_for_strata)),
        "gene_seeds_identical_under_nodes": all(np.array_equal(one, other) for one, other, drug in zip(on_source.perturbation_seeds, as_nodes.perturbation_seeds, is_drug) if not drug),
        "kept_positive_pairs": int(kept_positive.sum()), "drugs_with_a_kept_positive": int((kept_positive.sum(axis=1)[is_drug] > 0).sum())}

    progress("loader compared")
    # 2. the targets as a structure
    mechanism = np.asarray(as_nodes.edge_relation) == as_nodes.relation_types.index(DRUG_MECHANISM_RELATION)
    drugs_per_target = pd.Series(np.asarray(as_nodes.edge_target)[mechanism]).value_counts()
    most_shared = [(as_nodes.node_ids[int(node)], int(count)) for node, count in drugs_per_target.head(5).items()]
    relations_feeding = {}
    for node, relation in zip(np.asarray(on_source.edge_target), np.asarray(on_source.edge_relation)):
        relations_feeding.setdefault(int(node), set()).add(int(relation))
    feeding_counts = pd.Series([len(relations_feeding.get(int(node), ())) for node in drugs_per_target.index])
    results["targets"] = {"target_nodes": int(len(drugs_per_target)), "target_nodes_shared_by_two_or_more_drugs": int((drugs_per_target >= 2).sum()),
                          "median_drugs_per_target": float(drugs_per_target.median()), "most_drugs_on_one_target": int(drugs_per_target.max()), "most_shared_targets": most_shared,
                          "mechanism_edges_at_a_shared_target": int(drugs_per_target[drugs_per_target >= 2].sum()),
                          "relations_feeding_a_target_in_the_source_graph": {str(count): int(number) for count, number in feeding_counts.value_counts().sort_index().items()}}

    progress("targets counted")
    # 3. presence
    gene_positions = np.flatnonzero(~is_drug)
    sampled = np.sort(generator.choice(gene_positions, size=min(arguments.num_gene_knockouts, len(gene_positions)), replace=False))
    carrier_gene_nodes = {int(node) for node in np.asarray(as_nodes.edge_source)[np.asarray(as_nodes.edge_relation) == as_nodes.relation_types.index(SEQUESTERS_RELATION)]}
    target_nodes = set(int(node) for node in drugs_per_target.index)
    # the knockouts a drug node could reach back to: genes whose protein is a drug target or a carrier, read from the seeds' encoded proteins
    encodes = np.asarray(on_source.edge_relation) == on_source.relation_types.index("encodes")
    protein_of_gene_node = {}
    for gene_node, protein_node in zip(np.asarray(on_source.edge_source)[encodes], np.asarray(on_source.edge_target)[encodes]):
        protein_of_gene_node.setdefault(int(gene_node), set()).add(int(protein_node))
    touches_a_drug = np.array([any(protein_of_gene_node.get(int(seed), set()) & (target_nodes | carrier_gene_nodes) for seed in on_source.perturbation_seeds[position]) for position in gene_positions])
    near_drugs = gene_positions[touches_a_drug]
    measured_genes = np.unique(np.concatenate([sampled, near_drugs]))
    is_near_a_drug = np.isin(measured_genes, near_drugs)
    gene_seeds = [on_source.perturbation_seeds[position] for position in measured_genes]
    gene_signs = [on_source.perturbation_signs[position] for position in measured_genes]
    gene_magnitudes = [on_source.perturbation_magnitudes[position] for position in measured_genes]
    num_source_nodes = len(on_source.node_ids)

    entry_mp, entry_adjacencies = message_passing_encoder(as_nodes, arguments.message_passing_layers, entry=True)
    with torch.no_grad():
        entry_mp.layer_bias.copy_(torch.randn_like(entry_mp.layer_bias) * 0.1)
    naive_mp, naive_adjacencies = message_passing_encoder(as_nodes, arguments.message_passing_layers, entry=False)
    share_message_passing_parameters(entry_mp, as_nodes, naive_mp, as_nodes)
    source_mp, source_adjacencies = message_passing_encoder(on_source, arguments.message_passing_layers, entry=False)
    share_message_passing_parameters(entry_mp, as_nodes, source_mp, on_source)
    seeded_against_source, always_against_source, largest_at_absent_drug_node, knockouts_reaching_a_drug_node = RowDifferences(), RowDifferences(), 0.0, 0
    with torch.no_grad():
        for node_index, sign_and_magnitude in batches(gene_seeds, gene_signs, gene_magnitudes):
            source_field = source_mp.perturbation_difference_field(node_index, sign_and_magnitude, source_adjacencies)
            entry_field = entry_mp.perturbation_difference_field(node_index, sign_and_magnitude, entry_adjacencies)
            naive_field = naive_mp.perturbation_difference_field(node_index, sign_and_magnitude, naive_adjacencies)
            seeded_against_source.add(source_field, entry_field[:, :num_source_nodes])
            always_against_source.add(source_field, naive_field[:, :num_source_nodes])
            largest_at_absent_drug_node = max(largest_at_absent_drug_node, float(entry_field[:, num_source_nodes:].abs().max()))
            knockouts_reaching_a_drug_node += int((naive_field[:, num_source_nodes:].abs().amax(dim=(1, 2)) > TOLERANCE).sum())
    results["presence_message_passing"] = {
        "gene_knockouts_measured": int(len(measured_genes)), "of_which_a_drug_target_or_a_carrier": int(len(near_drugs)), "layers": arguments.message_passing_layers,
        "present_only_where_seeded_against_source": seeded_against_source.summary(),
        "largest_field_at_a_drug_node_present_only_where_seeded": largest_at_absent_drug_node,
        "always_present_against_source": always_against_source.summary(),
        "always_present_against_source_drug_target_or_carrier_knockouts": always_against_source.summary(is_near_a_drug),
        "always_present_against_source_other_knockouts": always_against_source.summary(~is_near_a_drug),
        "knockouts_with_a_nonzero_field_at_some_drug_node_always_present": knockouts_reaching_a_drug_node}
    progress("message-passing presence measured")
    entry_lr = linear_response_encoder(as_nodes, arguments.propagation_steps, entry=True)
    naive_lr = linear_response_encoder(as_nodes, arguments.propagation_steps, entry=False)
    share_linear_response_parameters(entry_lr, naive_lr)
    source_lr = linear_response_encoder(on_source, arguments.propagation_steps, entry=False)
    share_linear_response_parameters(entry_lr, source_lr)
    seeded_against_source, always_against_source, largest_at_absent_drug_node, knockouts_reaching_a_drug_node = RowDifferences(), RowDifferences(), 0.0, 0
    with torch.no_grad():
        for node_index, sign_and_magnitude in batches(gene_seeds, gene_signs, gene_magnitudes, batch_size=8):
            source_response, entry_response, naive_response = (encoder.response(node_index, sign_and_magnitude) for encoder in (source_lr, entry_lr, naive_lr))
            seeded_against_source.add(source_response, entry_response[:, :num_source_nodes])
            always_against_source.add(source_response, naive_response[:, :num_source_nodes])
            largest_at_absent_drug_node = max(largest_at_absent_drug_node, float(entry_response[:, num_source_nodes:].abs().max()))
            knockouts_reaching_a_drug_node += int((naive_response[:, num_source_nodes:].abs().amax(dim=(1, 2)) > 0).sum())
    results["presence_linear_response"] = {
        "gene_knockouts_measured": int(len(measured_genes)), "steps": arguments.propagation_steps,
        "present_only_where_seeded_against_source": seeded_against_source.summary(),
        "largest_response_at_a_drug_node_present_only_where_seeded": largest_at_absent_drug_node,
        "always_present_against_source": always_against_source.summary(),
        "always_present_against_source_drug_target_or_carrier_knockouts": always_against_source.summary(is_near_a_drug),
        "always_present_against_source_other_knockouts": always_against_source.summary(~is_near_a_drug),
        "knockouts_with_a_nonzero_response_at_some_drug_node_always_present": knockouts_reaching_a_drug_node}
    progress("linear-response presence measured")
    # 4. the hop
    drug_positions = np.flatnonzero(is_drug)
    target_seed_lists = [on_source.perturbation_seeds[position] for position in drug_positions]
    node_seed_lists = [as_nodes.perturbation_seeds[position] for position in drug_positions]
    layers, steps = arguments.message_passing_layers, arguments.propagation_steps
    mp_from_targets = reach_counts(num_source_nodes, np.asarray(on_source.edge_source), np.asarray(on_source.edge_target), target_seed_lists, layers + 1)
    mp_keep = np.asarray(as_nodes.edge_relation) != as_nodes.relation_types.index(SEQUESTERS_RELATION)
    mp_from_nodes = reach_counts(len(as_nodes.node_ids), np.asarray(as_nodes.edge_source)[mp_keep], np.asarray(as_nodes.edge_target)[mp_keep], node_seed_lists, layers + 1)
    lr_source, lr_target = linear_response_edges(on_source)
    lr_from_targets = reach_counts(num_source_nodes, lr_source, lr_target, target_seed_lists, steps)
    lr_variant_source, lr_variant_target = linear_response_edges(without_relation(as_nodes, SEQUESTERS_RELATION))
    lr_from_nodes = reach_counts(len(as_nodes.node_ids), lr_variant_source, lr_variant_target, node_seed_lists, steps)

    def reach_row(counts_targets: np.ndarray, counts_nodes: np.ndarray, hops: int) -> dict:
        from_targets = counts_targets[:, hops]
        from_node = counts_nodes[:, hops] - 1  # the drug node itself is not a node of the source graph
        return {"hops": hops, "median_nodes_from_the_targets": float(np.median(from_targets)), "median_nodes_from_the_drug_node": float(np.median(from_node)),
                "median_share_kept": float(np.median(from_node / np.maximum(from_targets, 1))), "drugs_that_reach_fewer_nodes": int((from_node < from_targets).sum())}

    results["reach"] = {
        "message_passing": [reach_row(mp_from_targets, mp_from_nodes, layers)],
        "message_passing_one_more_layer_from_the_drug_node": {
            "hops_from_the_drug_node": layers + 1, "drugs_reaching_exactly_the_nodes_count_of_the_targets_at_the_registered_depth": int((mp_from_nodes[:, layers + 1] - 1 == mp_from_targets[:, layers]).sum()),
            "median_nodes": float(np.median(mp_from_nodes[:, layers + 1] - 1))},
        "linear_response": [reach_row(lr_from_targets, lr_from_nodes, steps)],
        "graph_nodes": num_source_nodes,
        "gene_knockouts_beside_drugs": gene_knockout_reach(on_source, layers)}

    # the same question asked of the encoders: nodes of the source graph where a drug's difference field is not zero
    without_carriage_for_reach = without_relation(as_nodes, SEQUESTERS_RELATION)
    reach_encoders = {"seeded on the targets": (source_mp, source_adjacencies, target_seed_lists, on_source)}
    for name, num_layers, before_first_layer in (("drug node", layers, False), ("drug node, mechanism also before the first layer", layers, True), ("drug node, one more layer", layers + 1, False)):
        encoder, encoder_adjacencies = message_passing_encoder(without_carriage_for_reach, num_layers, entry=True, before_first_layer=before_first_layer)
        share_message_passing_parameters(entry_mp, as_nodes, encoder, without_carriage_for_reach)
        reach_encoders[name] = (encoder, encoder_adjacencies, node_seed_lists, as_nodes)
    nonzero_nodes: dict[str, list] = {name: [] for name in reach_encoders}
    with torch.no_grad():
        for name, (encoder, encoder_adjacencies, seed_lists, seed_data) in reach_encoders.items():
            seed_signs = [seed_data.perturbation_signs[position] for position in drug_positions]
            seed_magnitudes = [seed_data.perturbation_magnitudes[position] for position in drug_positions]
            for node_index, sign_and_magnitude in batches(seed_lists, seed_signs, seed_magnitudes):
                field = encoder.perturbation_difference_field(node_index, sign_and_magnitude, encoder_adjacencies)[:, :num_source_nodes]
                nonzero_nodes[name].extend((field.abs().amax(dim=2) > 0).sum(dim=1).tolist())
    reference_counts = np.array(nonzero_nodes["seeded on the targets"])
    results["reach"]["message_passing_field"] = [
        {"arm": name, "median_nodes_with_a_nonzero_field": float(np.median(counts)), "median_share_of_the_target_seeded_count": float(np.median(np.array(counts) / np.maximum(reference_counts, 1))),
         "drugs_with_as_many_nodes_as_seeded_on_the_targets": int((np.array(counts) >= reference_counts).sum())}
        for name, counts in nonzero_nodes.items()]
    progress("reach counted")
    # 5. the linear response of the two arms
    target_signs = [on_source.perturbation_signs[position] for position in drug_positions]
    target_magnitudes = [on_source.perturbation_magnitudes[position] for position in drug_positions]
    drug_signs = [as_nodes.perturbation_signs[position] for position in drug_positions]
    drug_magnitudes = [as_nodes.perturbation_magnitudes[position] for position in drug_positions]
    arm_rows, zero_signed = [], 0
    for num_steps in (steps, 4 * steps):
        entry_lr.num_propagation_steps, source_lr.num_propagation_steps = num_steps, num_steps
        cosines, relatives, nonzero_drugs = [], [], 0
        with torch.no_grad():
            for (node_index, node_values), (target_index, target_values) in zip(batches(node_seed_lists, drug_signs, drug_magnitudes, batch_size=8),
                                                                                 batches(target_seed_lists, target_signs, target_magnitudes, batch_size=8)):
                on_node = entry_lr.response(node_index, node_values)[:, :num_source_nodes]
                on_targets = source_lr.response(target_index, target_values) * torch.sigmoid(entry_lr.entry_gain_logit)[0]
                nonzero = on_targets.abs().amax(dim=(1, 2)) > 0
                nonzero_drugs += int(nonzero.sum())
                if bool(nonzero.any()):
                    flat_node, flat_targets = on_node[nonzero].flatten(1), on_targets[nonzero].flatten(1)
                    cosines.append(torch.nn.functional.cosine_similarity(flat_node, flat_targets, dim=1))
                    relatives.append((flat_node - flat_targets).abs().amax(dim=1) / flat_targets.abs().amax(dim=1))
        cosine, relative = torch.cat(cosines), torch.cat(relatives)
        arm_rows.append({"steps": num_steps, "drugs_with_a_nonzero_response": nonzero_drugs, "median_cosine_similarity": float(cosine.median()),
                         "smallest_cosine_similarity": float(cosine.min()), "median_largest_relative_difference": float(relative.median())})
        zero_signed = len(drug_positions) - nonzero_drugs
    entry_lr.num_propagation_steps, source_lr.num_propagation_steps = steps, steps
    results["linear_response_arms"] = {"rows": arm_rows, "drugs_with_a_zero_response_under_both_arms": zero_signed}
    progress("arms compared")
    # 6. the carriers
    carriers = pd.read_csv(arguments.carrier_table, dtype=str).fillna("")
    carrier_of = dict(zip(carriers.perturbation_id, carriers.carrier))
    has_carrier = np.array([as_nodes.perturbation_ids[position] in carrier_of for position in drug_positions])
    without_carriage = without_relation(as_nodes, SEQUESTERS_RELATION)
    lr_without = linear_response_encoder(without_carriage, steps, entry=True)
    share_linear_response_parameters(entry_lr, lr_without)
    mp_without, mp_without_adjacencies = message_passing_encoder(without_carriage, layers, entry=True)
    share_message_passing_parameters(entry_mp, as_nodes, mp_without, without_carriage)
    lr_differences, mp_differences, mp_scales = [], [], []
    with torch.no_grad():
        for node_index, sign_and_magnitude in batches(node_seed_lists, drug_signs, drug_magnitudes):
            lr_differences.append((entry_lr.response(node_index, sign_and_magnitude) - lr_without.response(node_index, sign_and_magnitude)).abs().amax(dim=(1, 2)))
            with_field = entry_mp.perturbation_difference_field(node_index, sign_and_magnitude, entry_adjacencies)
            without_field = mp_without.perturbation_difference_field(node_index, sign_and_magnitude, mp_without_adjacencies)
            mp_differences.append((with_field - without_field).abs().amax(dim=(1, 2)))
            mp_scales.append(without_field.abs().amax(dim=(1, 2)).clamp_min(1e-30))
    lr_difference, mp_difference, mp_scale = torch.cat(lr_differences), torch.cat(mp_differences), torch.cat(mp_scales)
    kept_per_drug = kept_positive.sum(axis=1)[drug_positions]
    with_carrier, without_carrier = kept_per_drug[has_carrier], kept_per_drug[~has_carrier]
    test = mannwhitneyu(with_carrier, without_carrier, alternative="two-sided")
    by_carrier = {}
    for name in ("albumin", "orosomucoid", "both"):
        of_kind = np.array([carrier_of.get(as_nodes.perturbation_ids[position]) == name for position in drug_positions])
        by_carrier[name] = {"drugs": int(of_kind.sum()), "kept_positives": int(kept_per_drug[of_kind].sum()), "drugs_with_a_kept_positive": int((kept_per_drug[of_kind] > 0).sum())}
    results["carriers"] = {
        "drugs_with_a_known_carrier": int(has_carrier.sum()), "drugs_without": int((~has_carrier).sum()),
        "linear_response": {"drugs_whose_response_the_sequestration_edges_change": int((lr_difference > 0).sum()), "largest_absolute_difference": float(lr_difference.max())},
        "message_passing": {"drugs_with_a_carrier_whose_field_changes": int((mp_difference[has_carrier] > TOLERANCE).sum()),
                            "drugs_without_a_carrier_whose_field_changes": int((mp_difference[~has_carrier] > TOLERANCE).sum()),
                            "median_relative_change_with_a_carrier": float((mp_difference / mp_scale)[has_carrier].median())},
        "kept_positives": {"with_a_carrier": {"drugs": int(len(with_carrier)), "total": int(with_carrier.sum()), "mean": float(with_carrier.mean()), "median": float(np.median(with_carrier)),
                                              "drugs_with_at_least_one": int((with_carrier > 0).sum())},
                           "without": {"drugs": int(len(without_carrier)), "total": int(without_carrier.sum()), "mean": float(without_carrier.mean()), "median": float(np.median(without_carrier)),
                                       "drugs_with_at_least_one": int((without_carrier > 0).sum())},
                           "mann_whitney_u": float(test.statistic), "two_sided_p": float(test.pvalue), "by_carrier": by_carrier}}

    results["carriers"]["sequestration_as_change"] = sequestration_as_change(as_nodes, arguments)
    progress("carriers measured")
    # 7. the drug classes
    atc_of_perturbation: dict[str, set[str]] = {}
    if arguments.onsides_bridge.exists():
        for ingredient in json.loads(arguments.onsides_bridge.read_text())["ingredients"].values():
            atc_of_perturbation.setdefault(str(ingredient.get("perturbation_id")), set()).update(str(code) for code in ingredient.get("atc_codes") or [])
    drug_ids = [as_nodes.perturbation_ids[position] for position in drug_positions]
    groups_per_drug = [{code[:3] for code in atc_of_perturbation.get(drug_id, set()) if code.startswith("N")} for drug_id in drug_ids]
    group_counts = {group: {"drugs": int(sum(group in groups for groups in groups_per_drug)),
                            "kept_positives": int(sum(int(count) for count, groups in zip(kept_per_drug, groups_per_drug) if group in groups)),
                            "with_a_known_carrier": int(sum(carried for carried, groups in zip(has_carrier, groups_per_drug) if group in groups))}
                    for group in sorted(ATC_GROUP_NAMES)}
    opioid = [any(code.startswith("N02A") for code in atc_of_perturbation.get(drug_id, set())) for drug_id in drug_ids]
    results["drug_classes"] = {"drugs_with_a_nervous_system_atc_code": int(sum(bool(groups) for groups in groups_per_drug)), "by_group": group_counts,
                               "opioids_n02a": {"drugs": int(sum(opioid)), "kept_positives": int(kept_per_drug[np.array(opioid)].sum())},
                               "anaesthetics_or_opioids": {"drugs": int(sum(flag or "N01" in groups for flag, groups in zip(opioid, groups_per_drug))),
                                                           "kept_positives": int(sum(int(count) for count, flag, groups in zip(kept_per_drug, opioid, groups_per_drug) if flag or "N01" in groups))},
                               "kept_positives_of_all_drugs": int(kept_per_drug.sum())}

    results["num_gene_knockouts_drawn"] = arguments.num_gene_knockouts
    (arguments.graph_dir / MEASUREMENTS_FILE).write_text(json.dumps(results, indent=1) + "\n")
    write_document(results, arguments)
    print(json.dumps({key: value for key, value in results.items() if key not in ("build",)}, indent=1)[:6000])


def gene_knockout_reach_lines(block: dict | None) -> list[str]:
    """The paragraph and table that set the drugs' reach beside the gene knockouts'; empty for measurements saved before
    the block existed."""
    if not block:
        return []
    genes, drugs, hops = block["gene_knockouts"], block["drugs_seeded_on_their_targets"], block["hops"]

    def row(name: str, of_kind: dict) -> list:
        return [name, of_kind["perturbations"], ", ".join(of_kind["seed_node_kinds"]), *[f"{value:,.0f}" for value in of_kind["median_nodes_by_hops"][1:]],
                " to ".join(f"{value:,.0f}" for value in of_kind["quartiles_at_the_registered_depth"])]

    return [
        "What the registered depth gives a gene knockout, for comparison. A gene knockout is seeded on a gene node, one `encodes` edge before its protein, and a drug "
        f"seeded on its targets is seeded on the protein nodes themselves, so at {hops} layers the drug is one edge ahead. A drug seeded on its own node is one edge "
        "before the protein, where the gene knockout is. Median nodes within k hops of the seeds, by breadth-first search on the source graph:",
        "",
        *markdown_table(["perturbations", "number", "seeded on", *[f"k = {hop}" for hop in range(1, hops + 2)], f"lower to upper quartile at k = {hops}"],
                        [row("gene knockouts", genes), row("drugs seeded on their targets", drugs)]),
        "",
        *saturation_lines(block),
    ]


def saturation_lines(block: dict) -> list[str]:
    """How far the registered depth is from carrying a signal everywhere it can go; empty for measurements saved before
    the search ran to saturation."""
    genes, drugs = block["gene_knockouts"], block["drugs_seeded_on_their_targets"]
    if "median_hops_to_everything_reachable" not in genes:
        return []
    hops, steps = block["hops"], block["linear_response_steps"]

    def row(name: str, of_kind: dict) -> list:
        return [name, f"{of_kind['median_nodes_reachable']:,.0f}", f"{of_kind['median_hops_to_everything_reachable']:.0f}", of_kind["largest_hops_to_everything_reachable"],
                f"{of_kind['median_hops_to_95_percent']:.0f}", f"{of_kind['median_share_reached_at_the_registered_depth']:.3f}",
                f"{of_kind['median_share_reached_at_the_linear_response_steps']:.3f}"]

    growing = genes["still_growing_at_the_largest_hop"] + drugs["still_growing_at_the_largest_hop"]
    return [
        f"How many hops a signal needs to reach every node it can reach, by the same search run to {block['largest_hop_searched']} hops "
        f"({growing} perturbations were still gaining nodes at the last hop). The share columns are of each perturbation's own reachable set, at the "
        f"message-passing depth of {hops} and at the {steps} steps of the linear-response default:",
        "",
        *markdown_table(["perturbations", "median nodes reachable", "median hops to all of them", "largest", "median hops to 95 percent",
                         f"median share at {hops} hops", f"median share at {steps} hops"],
                        [row("gene knockouts", genes), row("drugs seeded on their targets", drugs)]),
        "",
    ]


def write_document(results: dict, arguments) -> None:
    build, loader, targets, reach, carriers, classes = results["build"], results["loader"], results["targets"], results["reach"], results["carriers"], results["drug_classes"]
    presence_mp, presence_lr, arms = results["presence_message_passing"], results["presence_linear_response"], results["linear_response_arms"]
    kept = carriers["kept_positives"]
    mp_reach, lr_reach = reach["message_passing"][0], reach["linear_response"][0]
    more_layer = reach["message_passing_one_more_layer_from_the_drug_node"]

    def presence_rows(block: dict, unit: str) -> list[list]:
        seeded, always = block["present_only_where_seeded_against_source"], block["always_present_against_source"]
        return [["present only where seeded", f"{seeded['rows_changed']} of {seeded['rows']}", f"{seeded['largest_absolute_difference']:.1e}", "-", "-"],
                ["always present, averaged like any relation", f"{always['rows_changed']} of {always['rows']}", f"{always['largest_absolute_difference']:.1e}",
                 f"{always['median_relative_difference_of_changed_rows']:.1e}", f"{always['median_nodes_changed_per_changed_row']:.0f}"]]

    lines = [
        "# Drug entry nodes: what they change before anything is trained",
        "",
        "Generated by `experiments/measure_drug_entry_nodes.py`. Do not edit by hand.",
        "",
        f"Graph `{results['graph']}`, built from `{results['source_graph']}` by `experiments/build_drug_entry_node_variant.py`; evidence `{results['evidence']}`; "
        f"label selection `{results['label_selection']}`. Every encoder is at its random start (seed {results['seed']}, {NODE_STATE_DIM} state channels, structural node "
        "features only), and two encoders that are compared share every parameter they have in common. No parameter is fitted and no outcome chooses anything.",
        "",
        "## 1. What the build adds, and what the loader keeps",
        "",
        f"- drug nodes: **{build['drug_nodes_added']}**, one per drug perturbation; drugs attached to a graph compound instead: {len(build['drugs_attached_to_a_graph_compound'])}",
        f"- mechanism edges (`targets`): **{build['mechanism_edges_added']}**, by sign " + ", ".join(f"{sign}: {count}" for sign, count in build["mechanism_edges_by_sign"].items()),
        f"- sequestration edges (`sequesters`, sign -1): **{build['sequestration_edges_added']}** for {build['drugs_given_a_sequestration_edge']} drugs, by carrier gene "
        + ", ".join(f"{gene}: {count}" for gene, count in build["sequestration_edges_by_carrier_gene"].items()),
        f"- the graph: {build['nodes']:,} nodes and {build['edges']:,} edges (nodes sha256 `{build['nodes_sha256'][:16]}`, edges `{build['edges_sha256'][:16]}`)",
        "",
        f"Read with `drug_entry=\"targets\"`, the ablation arm, the variant gives the source graph's arrays: of {loader['fields_compared']} fields of the loaded data, "
        f"**{len(loader['fields_differing_under_targets'])}** differ" + (f" ({', '.join(loader['fields_differing_under_targets'])})" if loader["fields_differing_under_targets"] else "")
        + f", and the structural feature matrix is {'identical' if loader['structural_features_identical_under_targets'] else 'different'}.",
        "",
        f"Read with `drug_entry=\"nodes\"`, each of the {loader['drug_perturbations']} drug perturbations is seeded on its own node and the {loader['gene_perturbations']:,} gene "
        "perturbations keep their seeds. Identical to the source graph: "
        + ", ".join(f"{name} ({'yes' if same else 'NO'})" for name, same in loader["identical_under_nodes"].items())
        + f", degrees for strata ({'yes' if loader['degrees_for_strata_identical_under_nodes'] else 'NO'}). So the two arms share their folds, their label mask and their degree strata.",
        "",
        "## 2. How many study drugs share a target",
        "",
        f"The {build['mechanism_edges_added']} mechanism edges land on **{targets['target_nodes']}** target nodes, and **{targets['target_nodes_shared_by_two_or_more_drugs']}** of those "
        f"are the target of two or more study drugs ({targets['mechanism_edges_at_a_shared_target']} of the edges). The median target has {targets['median_drugs_per_target']:.0f} drugs "
        f"and the most shared has {targets['most_drugs_on_one_target']}: " + ", ".join(f"{node} ({count})" for node, count in targets["most_shared_targets"]) + ".",
        "",
        "This count is a property of which drugs the study contains. An encoder that averages a relation's messages over the edges entering a node divides each drug's signal at "
        "a shared target by it.",
        "",
        "## 3. Presence: a gene knockout's field with and without the drug nodes",
        "",
        f"{presence_mp['gene_knockouts_measured']} gene knockouts: {results.get('num_gene_knockouts_drawn', arguments.num_gene_knockouts)} drawn at random and every knockout of a gene whose protein is a drug target or a "
        f"carrier ({presence_mp['of_which_a_drug_target_or_a_carrier']}). Each row compares the field on the variant, restricted to the source graph's nodes, with the field on the source graph.",
        "",
        f"Message passing, {presence_mp['layers']} layers, difference field:",
        "",
        *markdown_table(["drug nodes", "knockouts whose field changes", "largest absolute difference", "median relative difference of those", "median nodes changed in those"],
                        presence_rows(presence_mp, "field")),
        "",
        f"Linear response, {presence_lr['steps']} steps:",
        "",
        *markdown_table(["drug nodes", "knockouts whose response changes", "largest absolute difference", "median relative difference of those", "median nodes changed in those"],
                        presence_rows(presence_lr, "response")),
        "",
        "A relative difference is the largest difference in a knockout's field divided by the largest entry of that field on the source graph.",
        "",
        "With the drug nodes always present, by whether the knocked-out gene's protein is a drug target or a carrier:",
        "",
        *markdown_table(["encoder", "knockouts", "whose field changes", "median relative difference of those", "largest absolute difference"],
                        [[encoder_name, group_name, f"{block['rows_changed']} of {block['rows']}", f"{block['median_relative_difference_of_changed_rows']:.1e}", f"{block['largest_absolute_difference']:.1e}"]
                         for encoder_name, presence in (("message passing", presence_mp), ("linear response", presence_lr))
                         for group_name, key in (("drug target or carrier", "always_present_against_source_drug_target_or_carrier_knockouts"), ("other", "always_present_against_source_other_knockouts"))
                         if key in presence for block in [presence[key]]]),
        "",
        f"With every drug node always present, knockouts that leave a nonzero value at some drug node: {presence_mp['knockouts_with_a_nonzero_field_at_some_drug_node_always_present']} "
        f"of {presence_mp['gene_knockouts_measured']} under message passing and {presence_lr['knockouts_with_a_nonzero_response_at_some_drug_node_always_present']} of "
        f"{presence_lr['gene_knockouts_measured']} under the linear response, which reaches a carrier within its {presence_lr['steps']} steps and follows the sequestration edge into "
        "the drugs that carrier binds. The model then computes a change in a drug nobody took and passes it on to that drug's targets. Present only where seeded, the largest value "
        f"at a drug node is {presence_mp['largest_field_at_a_drug_node_present_only_where_seeded']:.1e} and {presence_lr['largest_response_at_a_drug_node_present_only_where_seeded']:.1e}.",
        "",
        "## 4. The hop a drug node costs",
        "",
        "Nodes of the source graph that a drug's signal can reach, by breadth-first search along the edges each encoder propagates on, over the "
        f"{loader['drug_perturbations']} drugs. The graph has {reach['graph_nodes']:,} nodes.",
        "",
        *markdown_table(["encoder", "hops", "median nodes reached from the targets", "median from the drug node", "median share kept", "drugs that reach fewer nodes"],
                        [["message passing", mp_reach["hops"], f"{mp_reach['median_nodes_from_the_targets']:,.0f}", f"{mp_reach['median_nodes_from_the_drug_node']:,.0f}",
                          f"{mp_reach['median_share_kept']:.2f}", f"{mp_reach['drugs_that_reach_fewer_nodes']} of {loader['drug_perturbations']}"],
                         ["linear response", lr_reach["hops"], f"{lr_reach['median_nodes_from_the_targets']:,.0f}", f"{lr_reach['median_nodes_from_the_drug_node']:,.0f}",
                          f"{lr_reach['median_share_kept']:.2f}", f"{lr_reach['drugs_that_reach_fewer_nodes']} of {loader['drug_perturbations']}"]]),
        "",
        f"With one more message-passing layer ({more_layer['hops_from_the_drug_node']} from the drug node) the median reach is {more_layer['median_nodes']:,.0f} nodes, and "
        f"{more_layer['drugs_reaching_exactly_the_nodes_count_of_the_targets_at_the_registered_depth']} of {loader['drug_perturbations']} drugs reach as many nodes as their targets "
        f"reach at the registered depth of {mp_reach['hops']}.",
        "",
        "The same count read off the message-passing encoder itself, as the nodes of the source graph where a drug's difference field is not zero (sequestration edges left out, "
        "so the arms differ only in how the drug enters):",
        "",
        *markdown_table(["arm", "median nodes with a nonzero field", "median share of the target-seeded count", "drugs with at least as many nodes"],
                        [[row["arm"], f"{row['median_nodes_with_a_nonzero_field']:,.0f}", f"{row['median_share_of_the_target_seeded_count']:.2f}",
                          f"{row['drugs_with_as_many_nodes_as_seeded_on_the_targets']} of {loader['drug_perturbations']}"] for row in reach.get("message_passing_field", [])]),
        "",
        *gene_knockout_reach_lines(reach.get("gene_knockouts_beside_drugs")),
        "## 5. The linear-response arms are one model up to a gain",
        "",
        "The response of a drug seeded on its node, against its response seeded on its targets multiplied by the entry gain of each channel:",
        "",
        *markdown_table(["steps", "drugs with a nonzero response", "median cosine similarity", "smallest", "median largest relative difference"],
                        [[row["steps"], row["drugs_with_a_nonzero_response"], f"{row['median_cosine_similarity']:.6f}", f"{row['smallest_cosine_similarity']:.6f}",
                          f"{row['median_largest_relative_difference']:.2e}"] for row in arms["rows"]]),
        "",
        f"{arms['drugs_with_a_zero_response_under_both_arms']} drugs have a zero response under both arms, because every sign of their mechanism is 0 (a binding agent in ChEMBL's "
        "action types) and the linear response is sign times magnitude.",
        "",
        "## 6. The carriers' sequestration edges",
        "",
        f"{carriers['drugs_with_a_known_carrier']} drugs have a known carrier and {carriers['drugs_without']} do not. What the {build['sequestration_edges_added']} edges change in "
        "a drug's own field, the variant with them against the variant without them, when the edge is averaged in like any other relation and so carries the carrier's "
        "state (`--sequestration-carries presence` in the trainer):",
        "",
        f"- linear response: **{carriers['linear_response']['drugs_whose_response_the_sequestration_edges_change']}** drugs change (largest absolute difference "
        f"{carriers['linear_response']['largest_absolute_difference']:.1e});",
        f"- message passing: **{carriers['message_passing']['drugs_with_a_carrier_whose_field_changes']}** of the {carriers['drugs_with_a_known_carrier']} drugs with a carrier change "
        f"(median relative change {carriers['message_passing']['median_relative_change_with_a_carrier']:.3f}) and "
        f"{carriers['message_passing']['drugs_without_a_carrier_whose_field_changes']} of the {carriers['drugs_without']} without.",
        "",
        "Kept positive pairs per drug, by whether a carrier is known:",
        "",
        *markdown_table(["drugs", "number", "kept positives", "mean per drug", "median", "drugs with at least one"],
                        [["with a known carrier", kept["with_a_carrier"]["drugs"], kept["with_a_carrier"]["total"], f"{kept['with_a_carrier']['mean']:.2f}", f"{kept['with_a_carrier']['median']:.0f}",
                          kept["with_a_carrier"]["drugs_with_at_least_one"]],
                         ["without", kept["without"]["drugs"], kept["without"]["total"], f"{kept['without']['mean']:.2f}", f"{kept['without']['median']:.0f}", kept["without"]["drugs_with_at_least_one"]]]),
        "",
        f"Mann-Whitney U {kept['mann_whitney_u']:.0f}, two-sided p = {kept['two_sided_p']:.3g}. By carrier: "
        + "; ".join(f"{name} {block['drugs']} drugs, {block['kept_positives']} kept positives, {block['drugs_with_a_kept_positive']} with at least one" for name, block in kept["by_carrier"].items()) + ".",
        "",
        *sequestration_as_change_lines(carriers.get("sequestration_as_change")),
        "## 7. The drug classes of the study",
        "",
        f"{classes['drugs_with_a_nervous_system_atc_code']} of the {loader['drug_perturbations']} drugs carry an ATC nervous-system code in the identifier bridge; a drug with codes in two "
        f"groups is counted in both. The drugs hold {classes['kept_positives_of_all_drugs']} kept positives in all.",
        "",
        *markdown_table(["ATC group", "drugs", "kept positives", "with a known carrier"],
                        [[f"{group} {ATC_GROUP_NAMES[group]}", block["drugs"], block["kept_positives"], block["with_a_known_carrier"]] for group, block in classes["by_group"].items()]),
        "",
        f"Opioids (N02A): {classes['opioids_n02a']['drugs']} drugs, {classes['opioids_n02a']['kept_positives']} kept positives. Anaesthetics or opioids together: "
        f"{classes['anaesthetics_or_opioids']['drugs']} drugs, {classes['anaesthetics_or_opioids']['kept_positives']} kept positives.",
    ]
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")


if __name__ == "__main__":
    main()
