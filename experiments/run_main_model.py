"""Train and evaluate the proposed model B6, or the sigmoid-head baseline B3, on one split (design sections 5 and 6).

Resumable for a runtime-limited cluster: a checkpoint is written every --checkpoint-every-minutes
and on SIGUSR1; --resume continues from it; a DONE marker is written when the split finishes, which
slurm/train_resumable.sbatch waits for before it stops requeueing. One split per invocation lets
folds, seeds and module hold-outs run as array jobs.

Splits: --fold k of the grouped perturbation-wise split (leakage groups by --group-by), or
--holdout-module <module_id> for the pathway-wise split (every perturbation writing onto a gene of
that curated module is test data). A grouped validation subset of the training perturbations
(--validation-fraction) drives early stopping on macro AUPRC; the best validation state is restored
before the test evaluation.

Encoding: --field difference (default) reads forward(perturbed) - forward(unperturbed), so the head
sees only what the perturbation changed; --field absolute reads the raw field (base node states
included), the version 0.3 behaviour kept as an ablation. --pooling sum (default) or mean selects
how a module reads its gated support.

Objective per batch: evidence-weighted BCE (positives smoothed to --positive-target; unobserved pairs
are weak negatives at --negative-weight) + --description-length-coefficient * model cost (expected
active support nodes and nonzero links; zero for the sigmoid head).
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import signal
import subprocess
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.negative_controls import permute_symptom_labels_within_degree_strata
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    assign_grouped_folds,
    perturbations_anchored_in_module,
    primary_subsystem_by_gene_node,
    read_curated_modules,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    expected_calibration_error,
    hits_at_k,
    mean_reciprocal_rank,
    per_symptom_auprc,
    per_symptom_auroc,
)
from sklearn.metrics import average_precision_score, roc_auc_score
from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
from mechanistic_pathway_learning.models.baselines.local_descriptor_encoder import LocalDescriptorEncoder
from mechanistic_pathway_learning.models.baselines.zero_field_encoder import ZeroFieldEncoder
from mechanistic_pathway_learning.models.laboratory_readout import LaboratoryLabelIndex, LaboratoryReadout, laboratory_sign_loss
from mechanistic_pathway_learning.graph.cofactor_edges import CARRIER_RULES, cofactor_edge_mask
from mechanistic_pathway_learning.models.linear_response_encoder import (
    CROSS_RELATION_AGGREGATORS,
    EDGE_SIGNS,
    LinearResponseEncoder,
    MIXTURE_WEIGHTINGS,
    NORMALISATIONS,
)
from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from mechanistic_pathway_learning.models.soft_constraint_losses import evidence_weighted_binary_cross_entropy

CHECKPOINT_REQUESTED = False
MINIMUM_POSITIVES_TO_SCORE = 5


def request_checkpoint(signal_number, frame) -> None:  # noqa: ARG001
    global CHECKPOINT_REQUESTED
    CHECKPOINT_REQUESTED = True


def pad_perturbations(data, indices: np.ndarray) -> tuple[torch.Tensor, torch.Tensor]:
    max_nodes = max(1, max(len(data.perturbation_seeds[i]) for i in indices))
    node_index = torch.full((len(indices), max_nodes), -1, dtype=torch.long)
    sign_and_magnitude = torch.zeros((len(indices), max_nodes, 2))
    for row, i in enumerate(indices):
        seeds = data.perturbation_seeds[i]
        node_index[row, : len(seeds)] = torch.as_tensor(seeds, dtype=torch.long)
        sign_and_magnitude[row, : len(seeds), 0] = torch.as_tensor(data.perturbation_signs[i], dtype=torch.float32)
        sign_and_magnitude[row, : len(seeds), 1] = torch.as_tensor(data.perturbation_magnitudes[i], dtype=torch.float32)
    return node_index, sign_and_magnitude


def evaluate_laboratory_labels(encoder, readout, label_index, data, indices: np.ndarray, adjacencies, arguments, device) -> dict:
    """Sign agreement and AUROC of the readout's predicted metabolite changes against the measured directions of the
    given (held-out) perturbations."""
    from sklearn.metrics import roc_auc_score

    encoder.eval()
    scores, directions = [], []
    with torch.no_grad():
        for start in range(0, len(indices), arguments.batch_size):
            batch = indices[start : start + arguments.batch_size]
            if not label_index.count(batch):
                continue
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
            positions, nodes, slots, batch_directions = label_index.batch(batch, device)
            scores.append(readout.label_scores(field, positions, nodes, slots, len(batch_directions)).cpu().numpy())
            directions.append(batch_directions.cpu().numpy())
    encoder.train()
    if not scores:
        return {"labels": 0}
    scores, directions = np.concatenate(scores), np.concatenate(directions)
    reached = scores != 0
    return {"labels": int(len(scores)), "reached": int(reached.sum()),
            "sign_agreement": float(np.mean(np.sign(scores[reached]) == directions[reached])) if reached.any() else None,
            "majority_direction_rate": float(np.mean(directions == 1)),
            "auroc": float(roc_auc_score(directions == 1, scores)) if len(set(directions.tolist())) == 2 else None}


def perturbation_covariate(data) -> np.ndarray:
    """Standardised log1p of the summed degree of each perturbation's nodes: the hub signal that popularity and the
    random walk exploit, given to the heads explicitly under --degree-offset. It is a property of the input graph, not
    of the labels, so standardising over every perturbation leaks nothing."""
    log_degree = np.log1p(data.perturbation_degrees)
    return ((log_degree - log_degree.mean()) / max(log_degree.std(), 1e-8)).astype(np.float32)


def covariate_of(data, batch: np.ndarray, arguments, device):
    if not arguments.degree_offset:
        return None
    return torch.as_tensor(perturbation_covariate(data)[batch], device=device)


def node_feature_matrix(data, arguments) -> np.ndarray:
    """The fixed structural features, followed by the node descriptors of --node-descriptors when given (one block per
    node type, zero outside it: mechanistic_pathway_learning/graph/node_descriptors.py)."""
    structural = data.structural_node_features()
    if not getattr(arguments, "node_descriptors", None):
        return structural
    descriptors = pd.read_parquet(arguments.node_descriptors)
    missing = set(data.node_ids) - set(descriptors.index)
    if missing:
        raise ValueError(f"{len(missing)} graph nodes have no row in {arguments.node_descriptors}")
    return np.concatenate([structural, descriptors.loc[data.node_ids].to_numpy(dtype=np.float32)], axis=1)


def file_sha256(path) -> str | None:
    if not path:
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_models(data, arguments, device):
    if arguments.encoder == "none":
        encoder = ZeroFieldEncoder(len(data.node_ids), arguments.node_state_dim).to(device)
    elif arguments.encoder == "local_descriptors":
        encoder = LocalDescriptorEncoder(len(data.node_ids), torch.as_tensor(node_feature_matrix(data, arguments)), arguments.node_state_dim).to(device)
    elif arguments.encoder == "linear_response":
        if arguments.field != "difference":
            raise ValueError("the linear-response encoder is linear in its input, so its field is a difference field; use --field difference")
        cofactor_edges = None
        if arguments.cofactor_relations:
            cofactor_edges = torch.as_tensor(cofactor_edge_mask(data.edge_source, data.edge_target, data.edge_relation, data.relation_types,
                                                                data.node_base_metabolite_id, data.node_display_name, data.is_currency,
                                                                carrier_rule=arguments.carrier_rule))
            print(f"carrier edges given their own relations: {int(cofactor_edges.sum())} (rule {arguments.carrier_rule})")
        distinct_node_types = sorted(set(data.node_types.tolist()))
        node_type_index = torch.as_tensor([distinct_node_types.index(node_type) for node_type in data.node_types], dtype=torch.long)
        encoder = LinearResponseEncoder(len(data.node_ids), data.relation_types, torch.as_tensor(data.edge_source), torch.as_tensor(data.edge_target),
                                        torch.as_tensor(data.edge_relation), torch.as_tensor(data.edge_sign), torch.as_tensor(node_feature_matrix(data, arguments)),
                                        arguments.node_state_dim, non_propagating_nodes=torch.as_tensor(data.is_currency), cofactor_edges=cofactor_edges,
                                        num_propagation_steps=arguments.propagation_steps,
                                        propagation_channels=arguments.propagation_channels, response_scale=arguments.response_scale, damping=arguments.propagation_damping,
                                        normalisation=arguments.normalisation, cross_relation_aggregator=arguments.cross_relation_aggregator,
                                        mixture_weighting=arguments.mixture_weighting, node_type_index=node_type_index,
                                        edge_signs=arguments.edge_signs).to(device)
        print(f"linear-response encoder: {arguments.normalisation} normalisation, {arguments.cross_relation_aggregator} across relations"
              + (f" weighted by {arguments.mixture_weighting}" if arguments.cross_relation_aggregator == "softmax_mixture" else "")
              + f", {len(distinct_node_types)} node types ({', '.join(distinct_node_types)})"
              + (", every edge sign +1" if arguments.edge_signs == "all_positive" else ""))
    else:
        if arguments.node_descriptors and arguments.node_features != "typed":
            raise ValueError("--node-descriptors extends the typed node features; use --node-features typed")
        node_features = torch.as_tensor(node_feature_matrix(data, arguments)) if arguments.node_features == "typed" else None
        encoder = RelationalMessagePassingEncoder(len(data.node_ids), len(data.relation_types), arguments.node_state_dim, arguments.num_layers, node_features=node_features).to(device)
    if arguments.head == "sigmoid":
        head = RelationalGnnSigmoidHead(arguments.node_state_dim, len(data.symptoms), hidden_dim=arguments.sigmoid_hidden_dim, pooling=arguments.pooling,
                                        degree_offset=arguments.degree_offset).to(device)
    else:
        head = NoisyOrPathwayModuleHead(len(data.node_ids), arguments.node_state_dim, arguments.num_modules, len(data.symptoms),
                                        gate_initial_log_alpha=arguments.gate_initial_log_alpha, pooling=arguments.pooling,
                                        gate_initial_log_alpha_noise=arguments.gate_init_noise, initial_readout_bias=arguments.module_bias_init,
                                        degree_offset=arguments.degree_offset).to(device)
    return encoder, head


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def git_provenance() -> dict:
    """The commit the code runs from and any tracked files changed relative to it, so a result names the exact code
    that produced it (a run started from a tree with tracked changes is not reproducible from the commit alone)."""
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPOSITORY_ROOT, capture_output=True, text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPOSITORY_ROOT, capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "tracked_changes": None}
    return {"commit": commit, "tracked_changes": [line[3:] for line in status.splitlines() if line.strip()]}


def optimizer_parameter_groups(encoder, head, arguments, extra_parameters=()) -> list[dict]:
    """One group at the main learning rate, plus one group for each noisy-OR time scale (links, leaks, module biases,
    gates) whose learning rate is set; Adam moves a parameter by about one learning rate per step, so the rate is the
    time scale on which that parameter can change. extra_parameters (the laboratory readout) join the main group."""
    rates = {"links": arguments.link_learning_rate, "leaks": arguments.leak_learning_rate, "module_biases": arguments.module_bias_learning_rate,
             "gates": getattr(arguments, "gate_learning_rate", 0.0)}
    if not hasattr(head, "time_scale_parameter_groups") or not any(rates.values()):
        return [{"params": list(encoder.parameters()) + list(head.parameters()) + list(extra_parameters)}]
    separate_groups = [{"params": parameters, "lr": rates[name], "weight_decay": 0.0}
                       for name, parameters in head.time_scale_parameter_groups().items() if rates[name]]
    separated_ids = {id(parameter) for group in separate_groups for parameter in group["params"]}
    other_parameters = list(encoder.parameters()) + [parameter for parameter in head.parameters() if id(parameter) not in separated_ids] + list(extra_parameters)
    return [{"params": other_parameters}, *separate_groups]


def encode(encoder, node_index, sign_and_magnitude, adjacencies, field_kind: str):
    if field_kind == "difference":
        return encoder.perturbation_difference_field(node_index, sign_and_magnitude, adjacencies)
    return encoder(node_index, sign_and_magnitude, relation_adjacencies=adjacencies)


def predict(encoder, head, data, indices: np.ndarray, adjacencies, arguments, device) -> np.ndarray:
    encoder.eval()
    head.eval()
    predictions = []
    with torch.no_grad():
        for start in range(0, len(indices), arguments.batch_size):
            batch = indices[start : start + arguments.batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
            predictions.append(head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device)).symptom_probability.cpu().numpy())
    encoder.train()
    head.train()
    return np.concatenate(predictions, axis=0) if predictions else np.zeros((0, len(data.symptoms)))


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray) -> float:
    values = [per_symptom_auprc(predictions, outcomes, s) for s in range(outcomes.shape[1]) if MINIMUM_POSITIVES_TO_SCORE <= outcomes[:, s].sum() < outcomes.shape[0]]
    return float(np.mean(values)) if values else float("nan")


def split_indices(data, arguments) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Return (train, validation, test) index arrays and a name for the split directory."""
    all_indices = np.arange(len(data.perturbation_ids))
    if arguments.time_split_cutoff is not None:
        test_mask = np.ones(len(all_indices), dtype=bool)  # every perturbation is scored at the pair level; training uses pre-cutoff positives only
        split_name = f"time_{arguments.time_split_cutoff.isoformat()}_seed{arguments.seed}"
        train_pool = all_indices
        validation = np.array([], dtype=int)
        if arguments.validation_fraction > 0:
            num_validation_folds = max(2, int(round(1.0 / arguments.validation_fraction)))
            validation_fold = assign_grouped_folds(data.perturbation_ids, data.group_ids, num_validation_folds, arguments.seed + 1000)
            validation = np.array([i for i in all_indices if validation_fold[data.perturbation_ids[i]] == 0])
            train_pool = np.array([i for i in all_indices if validation_fold[data.perturbation_ids[i]] != 0])
        return train_pool, validation, all_indices[test_mask], split_name
    if arguments.holdout_module:
        module_genes = read_curated_modules(arguments.curated_modules)[arguments.holdout_module]
        module_nodes = {data.node_index[f"GENE:{symbol}"] for symbol in module_genes if f"GENE:{symbol}" in data.node_index}
        test_mask = np.array(perturbations_anchored_in_module(data.perturbation_seeds, module_nodes))
        split_name = f"module_{arguments.holdout_module}_seed{arguments.seed}"
    elif arguments.holdout_subsystem:
        if data.node_subsystem is None:
            raise ValueError("the graph has no subsystem column; rebuild it with the Human-GEM yml so reaction nodes carry subsystems")
        primary = primary_subsystem_by_gene_node(data.node_subsystem, data.edge_source, data.edge_target, data.edge_relation, data.relation_types.index("catalyzed_by"))
        subsystem_nodes = {gene_node for gene_node, subsystem in primary.items() if subsystem == arguments.holdout_subsystem}
        if not subsystem_nodes:
            raise ValueError(f"no gene has primary subsystem {arguments.holdout_subsystem!r}")
        test_mask = np.array(perturbations_anchored_in_module(data.perturbation_seeds, subsystem_nodes))
        safe_name = "".join(character if character.isalnum() else "_" for character in arguments.holdout_subsystem)
        split_name = f"subsystem_{safe_name}_seed{arguments.seed}"
    else:
        fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
        test_mask = np.array([fold_by_perturbation[p] == arguments.fold for p in data.perturbation_ids])
        split_name = f"fold{arguments.fold}_seed{arguments.seed}"
    train_pool = all_indices[~test_mask]
    validation = np.array([], dtype=int)
    if arguments.validation_fraction > 0 and len(train_pool) >= 20:
        num_validation_folds = max(2, int(round(1.0 / arguments.validation_fraction)))
        validation_fold = assign_grouped_folds([data.perturbation_ids[i] for i in train_pool], [data.group_ids[i] for i in train_pool], num_validation_folds, arguments.seed + 1000)
        validation = np.array([i for i in train_pool if validation_fold[data.perturbation_ids[i]] == 0])
        train_pool = np.array([i for i in train_pool if validation_fold[data.perturbation_ids[i]] != 0])
    return train_pool, validation, all_indices[test_mask], split_name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--run-dir", type=Path, default=Path("runs/main_model"))
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--holdout-module", type=str, default="", help="pathway-wise split: curated module id to hold out instead of a grouped fold")
    parser.add_argument("--holdout-subsystem", type=str, default="", help="pathway-wise split: Human-GEM subsystem whose genes (by primary subsystem) are held out")
    parser.add_argument("--permute-labels", action="store_true", help="negative control: permute outcome rows within degree strata before training and testing")
    parser.add_argument("--time-split-cutoff", type=date.fromisoformat, default=None,
                        help="monogenic time split (design 6.1): train on pairs dated on or before this day across all perturbations; score the pairs that could still become positive")
    parser.add_argument("--group-by", choices=["gene", "disease_cluster"], default="gene")
    parser.add_argument("--label-grades", nargs="*", default=["A", "B"], help="evidence grades that count as positive labels; pass A B C to keep grade C rows as the version 0.3 ablation did")
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--selection-metric", choices=["loss", "auprc"], default="loss",
                        help="early stopping on the validation evidence-weighted BCE (smooth on small validation sets) or on validation macro AUPRC")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--head", choices=["noisy_or", "sigmoid"], default="noisy_or")
    parser.add_argument("--field", choices=["difference", "absolute"], default="difference")
    parser.add_argument("--node-features", choices=["identity", "typed"], default="identity",
                        help="identity: a learned embedding per node; typed: fixed structural features only (type, compartment, degree, flags), the inductive variant")
    parser.add_argument("--pooling", choices=["sum", "mean"], default="sum")
    parser.add_argument("--encoder", choices=["message_passing", "linear_response", "none", "local_descriptors"], default="message_passing",
                        help="route 1 encoder: L layers of message passing, the time-invariant signed linear-response state space (linear_response_encoder.py), none (a zero field: with --degree-offset, the degree-only control) or local_descriptors (the perturbed node's own features and nothing from the graph: the descriptors-only control)")
    parser.add_argument("--laboratory-label-weight", type=float, default=0.0,
                        help="weight of the auxiliary loss on measured metabolite directions of the training genes (models/laboratory_readout.py); 0 leaves it out")
    parser.add_argument("--laboratory-labels", type=Path, default=Path("data/processed/laboratory_labels.parquet"),
                        help="label table of experiments/check_laboratory_label_coverage.py --labels-output")
    parser.add_argument("--node-descriptors", type=Path, default=None,
                        help="parquet of fixed node descriptors indexed by node_id (experiments/build_node_descriptors.py), appended to the structural node features")
    parser.add_argument("--degree-offset", action="store_true",
                        help="give the head the standardised log degree of each perturbation: a degree-dependent leak (noisy-OR) or logit offset (sigmoid), so the field only has to explain what degree does not")
    parser.add_argument("--propagation-steps", type=int, default=8, help="linear-response encoder: steps of the shared transition (the reach in edges)")
    parser.add_argument("--propagation-channels", type=int, default=4, help="linear-response encoder: channels propagated with their own gains (time scales), expanded linearly to --node-state-dim")
    parser.add_argument("--cofactor-relations", action="store_true",
                        help="linear-response encoder: give carrier edges (cofactor_edges.py) relations of their own, so their coupling gets learned gains")
    parser.add_argument("--carrier-rule", choices=CARRIER_RULES, default="all",
                        help="which clause of the carrier rule applies (cofactor_edges.py): all is the curated list and the recurring-pair "
                             "heuristic, names drops the heuristic, neuronal keeps only the transmitter-synthesis and oxidative carriers")
    parser.add_argument("--response-scale", choices=["linear", "signed_log"], default="linear",
                        help="linear-response encoder: read the response as it is, or through sign(h) log(1 + |h| / s) with a learned scale, so changes many edges away stay readable")
    parser.add_argument("--normalisation", choices=NORMALISATIONS, default="in_degree",
                        help="linear-response encoder: divide each message by the number of edges feeding the node (in_degree), or apply one global scale set by the spectral radius so the stoichiometric counts survive (spectral)")
    parser.add_argument("--edge-signs", choices=EDGE_SIGNS, default="graph",
                        help="linear-response encoder: the graph's signs, or every edge +1 with every gain positive (the sign ablation)")
    parser.add_argument("--cross-relation-aggregator", choices=CROSS_RELATION_AGGREGATORS, default="mean",
                        help="linear-response encoder: how the per-relation messages into a node combine; mean cancels exactly when a positive and a negative relation carry equal weight, softmax_mixture learns a weighting over order statistics instead")
    parser.add_argument("--mixture-weighting", choices=MIXTURE_WEIGHTINGS, default="sparsemax",
                        help="linear-response encoder: how --cross-relation-aggregator softmax_mixture turns its logits into weights; sparsemax can place exactly zero weight on a statistic, softmax cannot")
    parser.add_argument("--propagation-damping", type=float, default=0.5, help="linear-response encoder: weight of the new state per step (sets the transient, not the fixed point)")
    parser.add_argument("--node-state-dim", type=int, default=32)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--num-modules", type=int, default=8)
    parser.add_argument("--sigmoid-hidden-dim", type=int, default=64)
    parser.add_argument("--gate-initial-log-alpha", type=float, default=-1.0)
    parser.add_argument("--init-leak-from-base-rate", action="store_true", help="noisy-OR head: start each symptom's leak at its training base rate")
    parser.add_argument("--module-bias-init", type=float, default=0.0, help="noisy-OR head: initial readout bias of every module; negative values make modules off by default")
    parser.add_argument("--gate-init-noise", type=float, default=0.01, help="noisy-OR head: standard deviation of the per-gate noise added to the initial log-alpha (symmetry breaking between modules)")
    parser.add_argument("--link-learning-rate", type=float, default=0.0,
                        help="noisy-OR head: learning rate of the module-to-symptom links (0 = the main rate); links start at logit -3 and Adam moves them about one rate per step, so at 0.002 they need about 1,500 steps to reach 0.5")
    parser.add_argument("--leak-learning-rate", type=float, default=0.0, help="noisy-OR head: learning rate of the symptom leaks (0 = the main rate); slow, so a leak started at the base rate stays there")
    parser.add_argument("--module-bias-learning-rate", type=float, default=0.0, help="noisy-OR head: learning rate of the module readout biases (0 = the main rate)")
    parser.add_argument("--gate-learning-rate", type=float, default=0.0, help="noisy-OR head: learning rate of the support gate log-alphas (0 = the main rate); fast, so gates no field reaches can close within the epochs early stopping allows")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--learning-rate", type=float, default=0.002)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--max-epochs", type=int, default=60)
    parser.add_argument("--positive-target", type=float, default=0.99)
    parser.add_argument("--positive-target-from-frequency", action="store_true",
                        help="open question 8: use the reported HPO or label frequency (floored at --minimum-frequency-target) as the target of a positive pair instead of --positive-target")
    parser.add_argument("--minimum-frequency-target", type=float, default=0.05)
    parser.add_argument("--negative-weight", type=float, default=0.2)
    parser.add_argument("--description-length-coefficient", type=float, default=1e-6)
    parser.add_argument("--checkpoint-every-minutes", type=float, default=20.0)
    parser.add_argument("--time-budget-seconds", type=float, default=0.0, help="stop training after this many seconds (0 = no limit)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    arguments = parser.parse_args()
    signal.signal(signal.SIGUSR1, request_checkpoint)
    torch.manual_seed(arguments.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by, label_grades=tuple(arguments.label_grades) if arguments.label_grades else None)
    if arguments.permute_labels:
        data.outcomes = permute_symptom_labels_within_degree_strata(data.outcomes, data.perturbation_degrees, random_seed=arguments.seed)
    time_split = None
    if arguments.time_split_cutoff is not None:
        if data.evidence_dates is None or (data.evidence_dates > 0).sum() == 0:
            raise ValueError("the evidence table carries no evidence_date column; rebuild it with phenotype.hpoa present")
        cutoff_ordinal = arguments.time_split_cutoff.toordinal()
        dated = data.evidence_dates > 0
        full_outcomes = data.outcomes.copy()
        time_split = {
            "new_positive": ((full_outcomes > 0) & dated & (data.evidence_dates > cutoff_ordinal)).astype(float),
            "undated_positive": (full_outcomes > 0) & ~dated,
        }
        data.outcomes = ((full_outcomes > 0) & dated & (data.evidence_dates <= cutoff_ordinal)).astype(float)  # the model only ever sees pre-cutoff positives
        time_split["scored_pairs"] = ~(data.outcomes > 0) & ~time_split["undated_positive"]
    train_indices, validation_indices, test_indices, split_name = split_indices(data, arguments)
    if arguments.permute_labels:
        split_name += "_permuted"
    encoder, head = build_models(data, arguments, device)
    if arguments.init_leak_from_base_rate and arguments.head == "noisy_or":
        head.initialize_leak_from_base_rates(torch.as_tensor(data.outcomes[train_indices].mean(axis=0), dtype=torch.float32, device=device))
    adjacencies = None  # the linear-response encoder builds its signed adjacency from the edges at construction
    if arguments.encoder == "message_passing":
        adjacencies = [adjacency.to(device) if adjacency is not None else None for adjacency in RelationalMessagePassingEncoder.build_relation_adjacencies(
            torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long), torch.as_tensor(data.edge_relation, dtype=torch.long), len(data.node_ids), len(data.relation_types))]
    laboratory_index, laboratory_readout = None, None
    if arguments.laboratory_label_weight > 0:
        laboratory_index = LaboratoryLabelIndex(pd.read_parquet(arguments.laboratory_labels), data.perturbation_ids, data.node_base_metabolite_id)
        laboratory_readout = LaboratoryReadout(arguments.node_state_dim).to(device)
        print(f"laboratory labels: {laboratory_index.count(train_indices)} on training perturbations, {laboratory_index.count(test_indices)} held out")
    optimizer = torch.optim.AdamW(optimizer_parameter_groups(encoder, head, arguments, list(laboratory_readout.parameters()) if laboratory_readout is not None else ()),
                                  lr=arguments.learning_rate, weight_decay=arguments.weight_decay)
    split_directory = arguments.run_dir / split_name
    if arguments.resume and (split_directory / "DONE").exists():
        print(f"{split_name}: DONE marker present; skipping (delete the marker to retrain)")
        return
    split_directory.mkdir(parents=True, exist_ok=True)
    checkpoint_path = split_directory / "checkpoint.pt"
    state = {"epoch": 0, "history": [], "best_validation_auprc": float("-inf"), "best_validation_loss": float("inf"), "best_epoch": -1, "epochs_without_improvement": 0, "best_encoder": None, "best_head": None, "code_provenance": []}
    if arguments.resume and checkpoint_path.exists():
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        encoder.load_state_dict(checkpoint["encoder"])
        head.load_state_dict(checkpoint["head"])
        if laboratory_readout is not None and "laboratory_readout" in checkpoint:
            laboratory_readout.load_state_dict(checkpoint["laboratory_readout"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        state = checkpoint["state"]
        print(f"resumed at epoch {state['epoch']}")
    if "code_provenance" not in state:  # a checkpoint written before commits were recorded: its epochs ran under unrecorded code
        state["code_provenance"] = [{"started_at_epoch": 0, "commit": None, "tracked_changes": None}] if state["epoch"] > 0 else []
    provenance = git_provenance()
    if provenance["tracked_changes"]:
        print(f"warning: tracked files differ from commit {provenance['commit']}: {provenance['tracked_changes']}")
    state["code_provenance"].append({"started_at_epoch": state["epoch"], **provenance})  # one entry per process, so a resumed run lists every commit it ran under

    def save_checkpoint() -> None:
        torch.save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "state": state,
                    **({"laboratory_readout": laboratory_readout.state_dict()} if laboratory_readout is not None else {})}, checkpoint_path)

    outcomes = torch.as_tensor(data.outcomes, dtype=torch.float32)
    weights = torch.as_tensor(np.where(data.outcomes > 0, np.maximum(data.weights, 1e-3), arguments.negative_weight), dtype=torch.float32)
    positive_targets = None
    if arguments.positive_target_from_frequency:
        frequency = np.where(np.isnan(data.frequencies), arguments.positive_target, np.maximum(data.frequencies, arguments.minimum_frequency_target))
        positive_targets = torch.as_tensor(np.where(data.outcomes > 0, frequency, 0.0), dtype=torch.float32)
        weights = torch.as_tensor(np.where(data.outcomes > 0, 1.0, arguments.negative_weight), dtype=torch.float32)  # the frequency is the target, not the weight
    node_cost = float(np.log(len(data.node_ids)))
    started = time.time()
    last_checkpoint = time.time()
    generator = np.random.default_rng(arguments.seed + state["epoch"])
    stopped_early = False
    for epoch in range(state["epoch"], arguments.max_epochs):
        order = generator.permutation(train_indices)
        epoch_loss, epoch_bce, epoch_penalty, epoch_laboratory = 0.0, 0.0, 0.0, 0.0
        for start in range(0, len(order), arguments.batch_size):
            batch = order[start : start + arguments.batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
            output = head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
            if positive_targets is not None:
                bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, positive_targets[batch].to(device), weights[batch].to(device), positive_target=1.0)
            else:
                bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device), positive_target=arguments.positive_target)
            penalty = arguments.description_length_coefficient * head.description_length_penalty(node_cost=node_cost)
            loss = bce + penalty
            if laboratory_readout is not None and laboratory_index.count(batch):
                positions, nodes, slots, directions = laboratory_index.batch(batch, device)
                laboratory_loss = laboratory_sign_loss(laboratory_readout.label_scores(field, positions, nodes, slots, len(directions)), directions)
                loss = loss + arguments.laboratory_label_weight * laboratory_loss
                epoch_laboratory += laboratory_loss.item() * len(batch)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(batch)
            epoch_bce += bce.item() * len(batch)
            epoch_penalty += float(penalty) * len(batch)
            if CHECKPOINT_REQUESTED or (time.time() - last_checkpoint) / 60.0 >= arguments.checkpoint_every_minutes:
                save_checkpoint()
                last_checkpoint = time.time()
                if CHECKPOINT_REQUESTED:
                    print("checkpoint written on signal; exiting for requeue")
                    return
        entry = {"epoch": epoch, "train_loss": epoch_loss / max(1, len(order)), "train_bce": epoch_bce / max(1, len(order)), "train_penalty": epoch_penalty / max(1, len(order)), "elapsed_seconds": time.time() - started}
        if laboratory_readout is not None:
            entry["train_laboratory_loss"] = epoch_laboratory / max(1, len(order))
        if len(validation_indices):
            validation_predictions = predict(encoder, head, data, validation_indices, adjacencies, arguments, device)
            entry["validation_macro_auprc"] = macro_auprc(validation_predictions, data.outcomes[validation_indices])
            entry["validation_loss"] = float(evidence_weighted_binary_cross_entropy(torch.as_tensor(validation_predictions, dtype=torch.float32), outcomes[validation_indices], weights[validation_indices], positive_target=arguments.positive_target))
            improved = (entry["validation_loss"] < state["best_validation_loss"] - 1e-5) if arguments.selection_metric == "loss" else (entry["validation_macro_auprc"] > state["best_validation_auprc"] + 1e-4)
            if improved:
                state.update(best_validation_auprc=entry["validation_macro_auprc"], best_validation_loss=entry["validation_loss"], best_epoch=epoch, epochs_without_improvement=0,
                             best_encoder=copy.deepcopy(encoder.state_dict()), best_head=copy.deepcopy(head.state_dict()))
                if laboratory_readout is not None:
                    state["best_laboratory_readout"] = copy.deepcopy(laboratory_readout.state_dict())
            else:
                state["epochs_without_improvement"] += 1
        state["history"].append(entry)
        state["epoch"] = epoch + 1
        print(f"epoch {epoch} loss {entry['train_loss']:.4f} (bce {entry['train_bce']:.4f}, penalty {entry['train_penalty']:.4f})" + (f" validation loss {entry['validation_loss']:.4f} macro AUPRC {entry['validation_macro_auprc']:.3f}" if "validation_loss" in entry else "") + f" elapsed {entry['elapsed_seconds']:.0f}s", flush=True)
        save_checkpoint()
        if len(validation_indices) and state["epochs_without_improvement"] >= arguments.patience:
            stopped_early = True
            print(f"early stop at epoch {epoch}; best validation epoch {state['best_epoch']}")
            break
        if arguments.time_budget_seconds and time.time() - started > arguments.time_budget_seconds:
            print("time budget reached")
            break
    if state["best_encoder"] is not None:
        encoder.load_state_dict(state["best_encoder"])
        head.load_state_dict(state["best_head"])
        if laboratory_readout is not None and state.get("best_laboratory_readout") is not None:
            laboratory_readout.load_state_dict(state["best_laboratory_readout"])

    predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
    test_outcomes = data.outcomes[test_indices]
    laboratory_results = None
    if laboratory_readout is not None:
        laboratory_results = evaluate_laboratory_labels(encoder, laboratory_readout, laboratory_index, data, test_indices, adjacencies, arguments, device)
        print(f"held-out laboratory labels: {laboratory_results}")
    per_symptom = {}
    time_split_results = None
    if time_split is not None:
        generator = np.random.default_rng(arguments.seed)
        permuted = time_split["new_positive"].copy()
        for symptom_index in range(permuted.shape[1]):
            rows = np.where(time_split["scored_pairs"][:, symptom_index])[0]
            permuted[rows, symptom_index] = time_split["new_positive"][generator.permutation(rows), symptom_index]
        tables = {}
        for label_name, labels in (("observed", time_split["new_positive"]), ("permuted", permuted)):
            table = {}
            for symptom_index, symptom in enumerate(data.symptoms):
                mask = time_split["scored_pairs"][:, symptom_index]
                positives = labels[mask, symptom_index].sum()
                if positives < MINIMUM_POSITIVES_TO_SCORE or positives == mask.sum():
                    continue
                table[symptom] = {"scored_pairs": int(mask.sum()), "new_positives": int(positives), "base_rate": float(positives / mask.sum()),
                                  "auprc": float(average_precision_score(labels[mask, symptom_index], predictions[mask, symptom_index])),
                                  "auroc": float(roc_auc_score(labels[mask, symptom_index], predictions[mask, symptom_index]))}
            tables[label_name] = table
        time_split_results = {
            "cutoff": arguments.time_split_cutoff.isoformat(), "training_positive_pairs": int(data.outcomes.sum()), "new_positive_pairs": int(time_split["new_positive"].sum()),
            "undated_positive_pairs": int(time_split["undated_positive"].sum()), "scored_pairs": int(time_split["scored_pairs"].sum()),
            "per_symptom": tables["observed"], "per_symptom_permuted": tables["permuted"],
            "macro_auprc": float(np.mean([e["auprc"] for e in tables["observed"].values()])) if tables["observed"] else float("nan"),
            "macro_auroc": float(np.mean([e["auroc"] for e in tables["observed"].values()])) if tables["observed"] else float("nan"),
            "macro_auprc_permuted": float(np.mean([e["auprc"] for e in tables["permuted"].values()])) if tables["permuted"] else float("nan"),
            "macro_auroc_permuted": float(np.mean([e["auroc"] for e in tables["permuted"].values()])) if tables["permuted"] else float("nan"),
        }
        test_outcomes = time_split["new_positive"]  # ranking metrics below rank the new positives
    for symptom_index, symptom in enumerate(data.symptoms):
        if time_split is not None:
            break
        positives = test_outcomes[:, symptom_index].sum()
        if positives < MINIMUM_POSITIVES_TO_SCORE or positives == len(test_indices):
            continue
        per_symptom[symptom] = {
            "positives": int(positives),
            "auprc": bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap).__dict__,
            "auroc": bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap).__dict__,
            "base_rate": float(test_outcomes[:, symptom_index].mean()),
        }
    results = {
        "split": split_name, "fold": None if (arguments.holdout_module or arguments.holdout_subsystem) else arguments.fold, "holdout_module": arguments.holdout_module or None,
        "holdout_subsystem": arguments.holdout_subsystem or None, "labels_permuted": bool(arguments.permute_labels), "seed": arguments.seed,
        "arguments": {key: (value if isinstance(value, (int, float, str, bool, list, type(None))) else str(value)) for key, value in vars(arguments).items()},
        "num_train": int(len(train_indices)), "num_validation": int(len(validation_indices)), "num_test": int(len(test_indices)),
        "epochs_completed": state["epoch"], "best_epoch": state["best_epoch"], "stopped_early": stopped_early, "history": state["history"],
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, test_outcomes), "hits_at_3": hits_at_k(predictions, test_outcomes, 3),
        "expected_calibration_error": expected_calibration_error(predictions, test_outcomes), "per_symptom": per_symptom, "symptoms": data.symptoms,
        "test_perturbation_ids": [data.perturbation_ids[i] for i in test_indices],
        "time_split": time_split_results,
        "code_provenance": state["code_provenance"],
        "node_descriptors_sha256": file_sha256(getattr(arguments, "node_descriptors", None)),
        "laboratory_labels": laboratory_results,
    }
    if time_split_results is not None:
        results["macro_auprc"], results["macro_auroc"] = time_split_results["macro_auprc"], time_split_results["macro_auroc"]
    np.save(split_directory / "test_predictions.npy", predictions)
    if arguments.head == "noisy_or":
        with torch.no_grad():
            head.eval()
            support = head.module_support().cpu().numpy()
            expected_support = head.support_gate.expected_active_node_count().cpu().numpy()
            links = head.link_probability(0).cpu().numpy()
            leaks = head.leak_probability(0).cpu().numpy()  # at covariate zero (the mean log degree) under --degree-offset
            test_field_activations, test_leaks = [], []
            for start in range(0, len(test_indices), arguments.batch_size):
                batch = test_indices[start : start + arguments.batch_size]
                node_index, sign_and_magnitude = pad_perturbations(data, batch)
                output = head(encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field), perturbation_covariate=covariate_of(data, batch, arguments, device))
                test_field_activations.append(output.module_activation.cpu().numpy())
                test_leaks.append(np.broadcast_to(output.leak_probability.cpu().numpy(), (len(batch), output.leak_probability.shape[-1])))
            if arguments.degree_offset and test_leaks:
                np.save(split_directory / "test_leaks.npy", np.concatenate(test_leaks, axis=0))  # per-perturbation leaks for the sufficiency test
        results.update({
            "module_support_sizes": [int((support[k] > 0.5).sum()) for k in range(support.shape[0])],
            "module_expected_support_sizes": [float(x) for x in expected_support],
            "module_support_mass": [float(support[k].sum()) for k in range(support.shape[0])],
            "module_symptom_links": links.round(3).tolist(), "symptom_leaks": leaks.round(3).tolist(),
            "test_module_activation_mean": np.concatenate(test_field_activations, axis=0).mean(axis=0).round(3).tolist() if test_field_activations else [],
        })
        np.save(split_directory / "module_support.npy", support)
        np.save(split_directory / "test_module_activations.npy", np.concatenate(test_field_activations, axis=0) if test_field_activations else np.zeros((0, support.shape[0])))
        print("module support sizes:", results["module_support_sizes"])
    (split_directory / "results.json").write_text(json.dumps(results, indent=1))
    print(f"{split_name}: macro AUPRC {results['macro_auprc']:.3f} macro AUROC {results['macro_auroc']:.3f} MRR {results['mean_reciprocal_rank']:.3f} hits@3 {results['hits_at_3']:.3f} ECE {results['expected_calibration_error']:.3f}")
    (split_directory / "DONE").write_text("done\n")


if __name__ == "__main__":
    main()
