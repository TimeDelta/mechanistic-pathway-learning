"""Train and evaluate the proposed model B6, or the sigmoid-head baseline B3, on one split (design sections 5 and 6).

Resumable for a runtime-limited cluster: a checkpoint is written every --checkpoint-every-minutes
and on SIGUSR1; --resume continues from it; a DONE marker is written when the split finishes, which
slurm/train_resumable.sbatch waits for before it stops requeueing. One split per invocation lets
folds, seeds and module hold-outs run as array jobs.

Splits: --fold k of the grouped perturbation-wise split (leakage groups by --group-by), or
--holdout-module <module_id> or --holdout-subsystem <name> for the pathway-wise split (every perturbation writing onto
a gene of that curated module or primary subsystem is test data, and the perturbations sharing a leakage group with
one of them are left out of training). A grouped validation subset of the training perturbations
(--validation-fraction) drives early stopping (--selection-metric, the validation loss by default); the best
validation state is restored before the test evaluation. With --keep-large-groups-in-training a leakage group larger
than half the expected validation set never forms that subset (early_stopping_validation). With --refit-on-validation the
early-stopping run only chooses the number of epochs: a fresh model is then trained on the training and validation
perturbations together for best epoch + 1 epochs and scores the test (refit_checkpoint.pt, resumable), so it fits on
the same perturbations as the baselines, which need no validation set.

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
from collections import Counter
import hashlib
import json
import os
import signal
import subprocess
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.negative_controls import degree_stratified_row_permutation, duplicate_edge_count, fast_degree_preserving_rewiring, reciprocated_relations
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    assign_grouped_folds,
    perturbations_anchored_in_module,
    primary_subsystem_by_gene_node,
    read_curated_modules,
    training_mask_without_group_partners,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    expected_calibration_error,
    hits_at_k,
    mean_reciprocal_rank,
    micro_auprc,
    per_symptom_auprc,
    per_symptom_auroc,
    scorable_symptom,
    scored_rows,
)
from sklearn.metrics import average_precision_score, roc_auc_score
from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
from mechanistic_pathway_learning.models.baselines.local_descriptor_encoder import LocalDescriptorEncoder
from mechanistic_pathway_learning.models.baselines.zero_field_encoder import ZeroFieldEncoder
from mechanistic_pathway_learning.models.laboratory_readout import LaboratoryLabelIndex, LaboratoryReadout, laboratory_sign_loss
from mechanistic_pathway_learning.graph.brain_expression_weights import ALL_CELLS_CLASS, EXTRACELLULAR_COMPARTMENT
from mechanistic_pathway_learning.graph.cofactor_edges import CARRIER_RULES, cofactor_edge_mask
from mechanistic_pathway_learning.graph.node_descriptors import DESCRIPTOR_BLOCK_PREFIXES, descriptor_blocks
from mechanistic_pathway_learning.models.linear_response_encoder import (
    CROSS_RELATION_AGGREGATORS,
    EDGE_SIGNS,
    RELATION_GAINS,
    LinearResponseEncoder,
    MIXTURE_WEIGHTINGS,
    NORMALISATIONS,
)
from mechanistic_pathway_learning.models.descriptor_treatments import DESCRIPTOR_TREATMENTS
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
    dropped_blocks = getattr(arguments, "drop_descriptor_blocks", None) or []
    if dropped_blocks:
        blocks = descriptor_blocks(descriptors.columns)
        empty = [block for block in dropped_blocks if not blocks[block]]
        if empty:
            raise ValueError(f"--drop-descriptor-blocks names blocks with no column in {arguments.node_descriptors}: {empty}")
        descriptors = descriptors.drop(columns=[column for block in dropped_blocks for column in blocks[block]])
        if descriptors.shape[1] == 0:
            raise ValueError("--drop-descriptor-blocks leaves no descriptor column; leave out --node-descriptors instead")
    return np.concatenate([structural, descriptors.loc[data.node_ids].to_numpy(dtype=np.float32)], axis=1)


def descriptor_treatment_inputs(data, arguments) -> dict:
    """The encoder keywords of --descriptor-treatment: the treatment and the number of descriptor columns at the end of
    node_feature_matrix."""
    treatment = getattr(arguments, "descriptor_treatment", "plain")
    if treatment == "plain":
        return {}
    if not getattr(arguments, "node_descriptors", None):
        raise ValueError(f"--descriptor-treatment {treatment} needs --node-descriptors")
    if arguments.encoder not in ("message_passing", "linear_response"):
        raise ValueError(f"--descriptor-treatment {treatment} is defined for the message passing and linear-response encoders")
    if treatment == "zero_init_slow" and not getattr(arguments, "descriptor_learning_rate", 0.0):
        raise ValueError("--descriptor-treatment zero_init_slow needs --descriptor-learning-rate")
    if treatment != "zero_init_slow" and getattr(arguments, "descriptor_learning_rate", 0.0):
        raise ValueError("--descriptor-learning-rate is defined for --descriptor-treatment zero_init_slow only")
    num_structural_columns = data.structural_node_features().shape[1]
    num_descriptor_columns = node_feature_matrix(data, arguments).shape[1] - num_structural_columns
    return {"descriptor_treatment": treatment, "num_descriptor_columns": num_descriptor_columns}


def cell_class_inputs(data, arguments) -> tuple[torch.Tensor | None, torch.Tensor | None, torch.Tensor | None]:
    """Per-node class weights [N, K] of --cell-class-weights (experiments/build_cell_class_weights.py), and, with
    --extracellular-coupling, the extracellular metabolites as one pool shared by every class but the all-cells one
    (that class stays the propagation without cell classes)."""
    if not getattr(arguments, "cell_class_weights", None):
        if getattr(arguments, "extracellular_coupling", False):
            raise ValueError("--extracellular-coupling needs --cell-class-weights")
        return None, None, None
    table = pd.read_parquet(arguments.cell_class_weights)
    missing = set(data.node_ids) - set(table.index)
    if missing:
        raise ValueError(f"{len(missing)} graph nodes have no row in {arguments.cell_class_weights}")
    weights = torch.as_tensor(table.loc[data.node_ids].to_numpy(dtype=np.float32))
    if not getattr(arguments, "extracellular_coupling", False):
        return weights, None, None
    compartments = data.node_compartment if data.node_compartment is not None else np.full(len(data.node_ids), "")
    pool_nodes = torch.as_tensor((data.node_types == "metabolite") & (np.asarray(compartments, dtype=object) == EXTRACELLULAR_COMPARTMENT))
    shared_classes = torch.as_tensor([column != ALL_CELLS_CLASS for column in table.columns])
    return weights, pool_nodes, shared_classes


# arguments a resumed process may change without changing what the checkpoint was trained for
RESUME_CONTROL_ARGUMENTS = {"run_dir", "resume", "refit_on_validation", "max_epochs", "checkpoint_every_minutes", "time_budget_seconds", "timing_batches", "num_bootstrap"}


def array_sha256(values) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    return hashlib.sha256(str(array.dtype).encode() + str(array.shape).encode() + array.tobytes()).hexdigest()


def configuration_fingerprint(arguments, data, label_mask) -> dict:
    """The training arguments and digests of the inputs as loaded, recorded in the checkpoint so that a resume under
    other arguments or rebuilt inputs is reported (the split signature covers only the split and the controls)."""
    fingerprint = {f"argument:{name}": str(value) for name, value in sorted(vars(arguments).items()) if name not in RESUME_CONTROL_ARGUMENTS}
    fingerprint.update({
        "input:outcomes": array_sha256(data.outcomes), "input:weights": array_sha256(data.weights),
        "input:label_mask": None if label_mask is None else array_sha256(label_mask),
        "input:edges": array_sha256(np.stack([np.asarray(data.edge_source), np.asarray(data.edge_target), np.asarray(data.edge_relation)])),
        "input:seeds": hashlib.sha256(repr([(list(map(int, seeds)), list(map(float, signs)), list(map(float, magnitudes)))
                                            for seeds, signs, magnitudes in zip(data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes)]).encode()).hexdigest()})
    for name in ("node_descriptors", "cell_class_weights", "laboratory_labels"):
        value = getattr(arguments, name, None)
        fingerprint[f"file:{name}"] = file_sha256(value) if value and Path(value).is_file() else None
    return fingerprint


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
        cell_class_weights, pool_nodes, shared_pool_classes = cell_class_inputs(data, arguments)
        if cell_class_weights is not None:
            print(f"cell-class channels: {cell_class_weights.shape[1]} classes x {arguments.channels_per_cell_class} channels"
                  + (f", {int(pool_nodes.sum())} extracellular metabolites pooled across {int(shared_pool_classes.sum())} classes" if pool_nodes is not None else ""))
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
                                        edge_signs=arguments.edge_signs, relation_gains=arguments.relation_gains,
                                        cell_class_weights=cell_class_weights, channels_per_cell_class=arguments.channels_per_cell_class,
                                        extracellular_pool_nodes=pool_nodes, shared_pool_classes=shared_pool_classes,
                                        **descriptor_treatment_inputs(data, arguments)).to(device)
        print(f"linear-response encoder: {arguments.normalisation} normalisation, {arguments.cross_relation_aggregator} across relations"
              + (f" weighted by {arguments.mixture_weighting}" if arguments.cross_relation_aggregator == "softmax_mixture" else "")
              + f", {len(distinct_node_types)} node types ({', '.join(distinct_node_types)})"
              + {"all_positive": ", every edge sign +1", "permuted": ", edge signs permuted among edges"}.get(arguments.edge_signs, "")
              + (", one gain shared by every relation" if arguments.relation_gains == "shared" else ""))
    elif getattr(arguments, "cell_class_weights", None):
        raise ValueError("--cell-class-weights is defined for the linear-response encoder only (other encoders read cell-class expression through --node-descriptors)")
    else:
        if arguments.node_descriptors and arguments.node_features != "typed":
            raise ValueError("--node-descriptors extends the typed node features; use --node-features typed")
        node_features = torch.as_tensor(node_feature_matrix(data, arguments)) if arguments.node_features == "typed" else None
        encoder = RelationalMessagePassingEncoder(len(data.node_ids), len(data.relation_types), arguments.node_state_dim, arguments.num_layers, node_features=node_features,
                                                  **descriptor_treatment_inputs(data, arguments)).to(device)
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


def atomic_torch_save(payload: dict, path: Path) -> None:
    """torch.save to a temporary file, then rename it over the checkpoint, so a kill during the write leaves the
    previous checkpoint intact instead of a truncated one that every later resume fails to load."""
    temporary_path = path.with_name(path.name + ".partial")
    torch.save(payload, temporary_path)
    os.replace(temporary_path, path)


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
    time scale on which that parameter can change. extra_parameters (the laboratory readout) join the main group. Under
    --descriptor-treatment zero_init_slow the encoder's descriptor map is one more group, at --descriptor-learning-rate
    with the main weight decay."""
    rates = {"links": arguments.link_learning_rate, "leaks": arguments.leak_learning_rate, "module_biases": arguments.module_bias_learning_rate,
             "gates": getattr(arguments, "gate_learning_rate", 0.0)}
    descriptor_parameters = encoder.descriptor_parameters() if hasattr(encoder, "descriptor_parameters") else []
    descriptor_groups = [{"params": descriptor_parameters, "lr": arguments.descriptor_learning_rate}] if descriptor_parameters else []
    descriptor_ids = {id(parameter) for parameter in descriptor_parameters}
    encoder_parameters = [parameter for parameter in encoder.parameters() if id(parameter) not in descriptor_ids]
    if not hasattr(head, "time_scale_parameter_groups") or not any(rates.values()):
        return [{"params": encoder_parameters + list(head.parameters()) + list(extra_parameters)}, *descriptor_groups]
    separate_groups = [{"params": parameters, "lr": rates[name], "weight_decay": 0.0}
                       for name, parameters in head.time_scale_parameter_groups().items() if rates[name]]
    separated_ids = {id(parameter) for group in separate_groups for parameter in group["params"]}
    other_parameters = encoder_parameters + [parameter for parameter in head.parameters() if id(parameter) not in separated_ids] + list(extra_parameters)
    return [{"params": other_parameters}, *separate_groups, *descriptor_groups]


def encode(encoder, node_index, sign_and_magnitude, adjacencies, field_kind: str):
    if field_kind == "difference":
        return encoder.perturbation_difference_field(node_index, sign_and_magnitude, adjacencies)
    return encoder(node_index, sign_and_magnitude, relation_adjacencies=adjacencies)


def predict(encoder, head, data, indices: np.ndarray, adjacencies, arguments, device, with_log_parts: bool = False):
    """Symptom probabilities [len(indices), S]; with_log_parts, also log P and log(1 - P) from the head (for the
    validation loss under --bce-in-log-space)."""
    encoder.eval()
    head.eval()
    predictions, log_probabilities, log_complements = [], [], []
    with torch.no_grad():
        for start in range(0, len(indices), arguments.batch_size):
            batch = indices[start : start + arguments.batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
            output = head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
            predictions.append(output.symptom_probability.cpu().numpy())
            if with_log_parts:
                log_probabilities.append(output.symptom_log_probability.cpu().numpy())
                log_complements.append(output.symptom_log_complement.cpu().numpy())
    encoder.train()
    head.train()
    empty = np.zeros((0, len(data.symptoms)))
    probabilities = np.concatenate(predictions, axis=0) if predictions else empty
    if not with_log_parts:
        return probabilities
    return probabilities, (np.concatenate(log_probabilities, axis=0) if log_probabilities else empty), (np.concatenate(log_complements, axis=0) if log_complements else empty)


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray | None = None) -> float:
    values = [per_symptom_auprc(predictions, outcomes, s, mask) for s in range(outcomes.shape[1]) if scorable_symptom(outcomes, s, MINIMUM_POSITIVES_TO_SCORE, mask)]
    return float(np.mean(values)) if values else float("nan")


def initialise_leaks_from_base_rates(head, data, fit_indices: np.ndarray, label_mask, arguments, device) -> None:
    """--init-leak-from-base-rate: start each noisy-OR leak at its symptom's base rate over the perturbations the model
    is fitted on (over the labelled pairs only under a label selection)."""
    if not (arguments.init_leak_from_base_rate and arguments.head == "noisy_or"):
        return
    if label_mask is None:
        base_rates = data.outcomes[fit_indices].mean(axis=0)
    else:  # over the labelled pairs only
        base_rates = (data.outcomes[fit_indices] * label_mask[fit_indices]).sum(axis=0) / np.maximum(label_mask[fit_indices].sum(axis=0), 1)
    head.initialize_leak_from_base_rates(torch.as_tensor(base_rates, dtype=torch.float32, device=device))


def early_stopping_validation(data, pool: np.ndarray, arguments) -> tuple[np.ndarray, np.ndarray]:
    """Split a training pool into (training, validation): fold 0 of a grouped split of the pool into
    round(1 / --validation-fraction) folds, so no leakage group is on both sides.

    assign_grouped_folds places the largest group first, in fold 0, so without --keep-large-groups-in-training a group
    larger than the expected validation set becomes the whole of it. On the full development data with the lockbox
    removed that is one group of 232 of 1,222 perturbations holding 27 percent of the positives and 63 drugs, which a
    lockbox run (whose pool is every development perturbation) would never train on, for every seed. With the flag,
    groups larger than half the expected validation set (--validation-fraction x pool size / 2) stay in training and the
    folds are drawn over the remaining groups; when no group is that large the split is the one without the flag."""
    num_validation_folds = max(2, int(round(1.0 / arguments.validation_fraction)))
    if getattr(arguments, "validation_draw", "fold0") == "rotated_stratified":
        return rotated_stratified_validation(data, pool, arguments, num_validation_folds)
    eligible = pool
    if getattr(arguments, "keep_large_groups_in_training", False):
        group_size_limit = arguments.validation_fraction * len(pool) / 2
        group_sizes = Counter(data.group_ids[i] for i in pool)
        eligible = np.array([i for i in pool if group_sizes[data.group_ids[i]] <= group_size_limit], dtype=int)
        kept = {group: size for group, size in group_sizes.items() if size > group_size_limit}
        if kept:
            print(f"early-stopping validation: {len(kept)} leakage group(s) larger than {group_size_limit:.0f} perturbations kept in training "
                  f"({sum(kept.values())} perturbations, largest {max(kept.values())})")
        if len(eligible) == 0:
            raise ValueError("every leakage group is larger than half the expected validation set; lower --validation-fraction or drop --keep-large-groups-in-training")
    validation_fold = assign_grouped_folds([data.perturbation_ids[i] for i in eligible], [data.group_ids[i] for i in eligible], num_validation_folds, arguments.seed + 1000)
    validation = np.array([i for i in eligible if validation_fold[data.perturbation_ids[i]] == 0], dtype=int)
    in_validation = set(validation.tolist())
    return np.array([i for i in pool if i not in in_validation], dtype=int), validation


VALIDATION_PARTITION_SEED = 1000
VALIDATION_GROUP_SIZE_DIVISOR = 4


def rotated_stratified_validation(data, pool: np.ndarray, arguments, num_validation_folds: int) -> tuple[np.ndarray, np.ndarray]:
    """--validation-draw rotated_stratified (the user's decisions of 8 October 2026): the pool's groups are split once into
    num_validation_folds grouped folds, separately for groups holding a drug and the rest, so each fold has about the pool's
    share of each; the validation set is fold (seed mod num_validation_folds), so seeds 0 to 4 stop on five disjoint
    validation sets instead of one. Groups larger than a quarter of the expected validation set (--validation-fraction x
    pool size / 4) stay in training, so no one leakage group is more than about a quarter of a validation set; the fold-0
    draw allowed half, and its validation set was 41 percent one disease cluster (cluster:ARNT2, 58 of 142). The
    partition uses a fixed seed (VALIDATION_PARTITION_SEED), not the run's seed, so that the folds and hence the five
    validation sets are disjoint."""
    group_size_limit = arguments.validation_fraction * len(pool) / VALIDATION_GROUP_SIZE_DIVISOR
    group_sizes = Counter(data.group_ids[i] for i in pool)
    eligible = [i for i in pool if group_sizes[data.group_ids[i]] <= group_size_limit]
    if not eligible:
        raise ValueError("every leakage group is larger than a quarter of the expected validation set; lower --validation-fraction")
    holds_a_drug = {data.group_ids[i] for i in eligible if data.perturbation_types[i] == "drug"}
    fold_of: dict[str, int] = {}
    for offset, stratum in enumerate((True, False)):  # the offset keeps the largest group of each stratum out of one shared fold
        members = [i for i in eligible if (data.group_ids[i] in holds_a_drug) == stratum]
        if members:
            stratum_folds = assign_grouped_folds([data.perturbation_ids[i] for i in members], [data.group_ids[i] for i in members],
                                                 num_validation_folds, VALIDATION_PARTITION_SEED)
            fold_of.update({perturbation_id: (fold + offset) % num_validation_folds for perturbation_id, fold in stratum_folds.items()})
    validation_fold = arguments.seed % num_validation_folds
    validation = np.array([i for i in eligible if fold_of[data.perturbation_ids[i]] == validation_fold], dtype=int)
    in_validation = set(validation.tolist())
    kept = sorted(size for size in group_sizes.values() if size > group_size_limit)
    print(f"early-stopping validation: fold {validation_fold} of {num_validation_folds} (seed {arguments.seed}), drawn per stratum (groups with a drug, the rest); "
          f"{len(kept)} group(s) larger than {group_size_limit:.0f} perturbations kept in training" + (f" (largest {kept[-1]})" if kept else ""))
    return np.array([i for i in pool if i not in in_validation], dtype=int), validation


def split_indices(data, arguments, in_lockbox: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Return (train, validation, test) index arrays and a name for the split directory. With in_lockbox (--score-lockbox)
    the test set is the lockbox and training and validation come from the development perturbations."""
    all_indices = np.arange(len(data.perturbation_ids))
    if arguments.time_split_cutoff is not None:
        test_mask = np.ones(len(all_indices), dtype=bool)  # every perturbation is scored at the pair level; training uses pre-cutoff positives only
        split_name = f"time_{arguments.time_split_cutoff.isoformat()}_seed{arguments.seed}"
        train_pool = all_indices
        validation = np.array([], dtype=int)
        if arguments.validation_fraction > 0:
            train_pool, validation = early_stopping_validation(data, all_indices, arguments)
        return train_pool, validation, all_indices[test_mask], split_name
    if in_lockbox is not None:
        test_mask = np.asarray(in_lockbox, dtype=bool)
        split_name = f"lockbox_seed{arguments.seed}"
    elif arguments.holdout_module:
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
    if arguments.holdout_module or arguments.holdout_subsystem:  # the hold-out is chosen by seed gene, not by leakage group
        training_mask = np.array(training_mask_without_group_partners(test_mask, data.group_ids))
        num_group_partners = int((~test_mask & ~training_mask).sum())
        if num_group_partners:
            print(f"pathway-wise hold-out: {num_group_partners} perturbations share a leakage group with the held-out ones and are left out of training")
        train_pool = all_indices[training_mask]
    else:  # whole groups: the lockbox (read_lockbox refuses a straddling group) and grouped folds
        train_pool = all_indices[~test_mask]
    validation = np.array([], dtype=int)
    if arguments.validation_fraction > 0 and len(train_pool) >= 20:
        train_pool, validation = early_stopping_validation(data, train_pool, arguments)
    return train_pool, validation, all_indices[test_mask], split_name


def time_training_steps(data, encoder, head, optimizer, train_indices, validation_indices, adjacencies, label_mask, arguments, device) -> None:
    """Seconds per training step (forward, backward, optimizer step) and per prediction batch; nothing is scored."""
    outcomes = torch.as_tensor(data.outcomes, dtype=torch.float32)
    weights = torch.as_tensor(np.where(data.outcomes > 0, np.maximum(data.weights, 1e-3), arguments.negative_weight), dtype=torch.float32)
    if label_mask is not None:
        weights = weights * torch.as_tensor(label_mask, dtype=torch.float32)
    order = np.random.default_rng(arguments.seed).permutation(train_indices)
    step_seconds = []
    for step in range(arguments.timing_batches):
        batch = order[(step * arguments.batch_size) % len(order):][: arguments.batch_size]
        started = time.time()
        node_index, sign_and_magnitude = pad_perturbations(data, batch)
        field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
        output = head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
        loss = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device), positive_target=arguments.positive_target)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        step_seconds.append(time.time() - started)
    started = time.time()
    with torch.no_grad():
        encoder.eval()
        head.eval()
        node_index, sign_and_magnitude = pad_perturbations(data, order[: arguments.batch_size])
        head(encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field), perturbation_covariate=covariate_of(data, order[: arguments.batch_size], arguments, device))
    prediction_seconds = time.time() - started
    steady = float(np.median(step_seconds[1:] if len(step_seconds) > 1 else step_seconds))
    steps_per_epoch = int(np.ceil(len(train_indices) / arguments.batch_size))
    prediction_batches = int(np.ceil(len(validation_indices) / arguments.batch_size))
    print(json.dumps({"timing_batches": arguments.timing_batches, "first_step_seconds": round(step_seconds[0], 2), "median_step_seconds": round(steady, 2),
                      "prediction_batch_seconds": round(prediction_seconds, 2), "train_perturbations": int(len(train_indices)), "steps_per_epoch": steps_per_epoch,
                      "epoch_minutes_estimate": round((steady * steps_per_epoch + prediction_seconds * prediction_batches) / 60.0, 1),
                      "torch_threads": torch.get_num_threads()}))


def build_argument_parser() -> argparse.ArgumentParser:
    """The trainer's arguments; experiments/evaluate_on_rewired_graph.py rebuilds a finished run from them."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--run-dir", type=Path, default=Path("runs/main_model"))
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--holdout-module", type=str, default="", help="pathway-wise split: curated module id to hold out instead of a grouped fold")
    parser.add_argument("--holdout-subsystem", type=str, default="", help="pathway-wise split: Human-GEM subsystem whose genes (by primary subsystem) are held out")
    parser.add_argument("--permute-labels", action="store_true",
                        help="negative control: permute perturbation rows (outcomes with their pair weights, frequencies and label mask) within degree strata, and within the lockbox and the development set apart, before training and testing")
    parser.add_argument("--lockbox", type=Path, default=None,
                        help="lockbox file from experiments/draw_lockbox.py; its perturbations are removed before the split, so pilot folds never see them")
    parser.add_argument("--score-lockbox", action="store_true",
                        help="with --lockbox: train on every development perturbation (validation drawn from them) and score the lockbox once; the confirmatory runs (docs/preregistration.md)")
    parser.add_argument("--rewire-swaps-per-edge", type=int, default=0,
                        help="negative control: degree-preserving rewiring of the graph within each relation before the model is built (attempted swaps per edge; 0 = the real graph); seeded by --seed")
    parser.add_argument("--keep-reciprocated-relations-symmetric", action="store_true",
                        help="with --rewire-swaps-per-edge: a relation stored in both directions (binds) is rewired as undirected edges, so it stays symmetric "
                             "(negative_controls.reciprocated_relations, fast_degree_preserving_rewiring); off by default, which reproduces earlier rewired runs")
    parser.add_argument("--time-split-cutoff", type=date.fromisoformat, default=None,
                        help="monogenic time split (design 6.1): train on pairs dated on or before this day across all perturbations; score the pairs that could still become positive")
    parser.add_argument("--group-by", choices=["gene", "disease_cluster", "disease_cluster_and_targets"], default="gene")
    parser.add_argument("--label-grades", nargs="*", default=["A", "B"], help="evidence grades that count as positive labels; pass A B C to keep grade C rows as the version 0.3 ablation did")
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="parquet of (perturbation_id, symptom, keep) from experiments/build_label_selection.py; positive pairs with keep False are masked out of the loss and every metric")
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--refit-on-validation", action="store_true",
                        help="after early stopping, train a fresh model on the training and validation perturbations together for best epoch + 1 epochs "
                             "and score the test with it (the early-stopped model's test predictions are kept in test_predictions_early_stopped.npy)")
    parser.add_argument("--keep-large-groups-in-training", action="store_true",
                        help="leakage groups larger than half the expected validation set stay in training and never form the early-stopping validation set "
                             "(without it the largest group goes there first; see early_stopping_validation)")
    parser.add_argument("--validation-draw", choices=["fold0", "rotated_stratified"], default="fold0",
                        help="fold0: the early-stopping validation set is fold 0 of the grouped split of the pool (the same set for every seed); "
                             "rotated_stratified: fold seed mod 7 of a split drawn per stratum (groups with a drug, the rest), groups above a quarter of the "
                             "expected validation set kept in training (rotated_stratified_validation)")
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
    parser.add_argument("--drop-descriptor-blocks", nargs="*", choices=list(DESCRIPTOR_BLOCK_PREFIXES), default=[],
                        help="with --node-descriptors, leave these blocks out (graph/node_descriptors.py DESCRIPTOR_BLOCK_PREFIXES); the descriptors of gene "
                             "nodes are protein and gene_brain, and reaction_brain carries gene expression onto reactions")
    parser.add_argument("--descriptor-treatment", choices=list(DESCRIPTOR_TREATMENTS), default="plain",
                        help="with --node-descriptors, how they enter the message passing base state or the linear-response output gate (models/descriptor_treatments.py): "
                             "plain (one map of all the columns), seed_masked (each perturbation's perturbed nodes read their structural columns only) or "
                             "zero_init_slow (a separate descriptor map, initialised at zero and trained at --descriptor-learning-rate)")
    parser.add_argument("--descriptor-learning-rate", type=float, default=0.0,
                        help="--descriptor-treatment zero_init_slow: learning rate of the descriptor map (required there)")
    parser.add_argument("--cell-class-weights", type=Path, default=None,
                        help="linear-response encoder: parquet of per-node cell-class weights (experiments/build_cell_class_weights.py); the response then propagates once per class")
    parser.add_argument("--channels-per-cell-class", type=int, default=1, help="linear-response encoder with --cell-class-weights: channels (time scales) per class")
    parser.add_argument("--extracellular-coupling", action="store_true",
                        help="linear-response encoder with --cell-class-weights: extracellular metabolites are one pool shared by every class but all_cells")
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
                        help="linear-response encoder: the graph's signs, every edge +1 with every gain positive (the sign ablation), "
                             "or the signs shuffled among edges with a fixed seed (the sign permutation)")
    parser.add_argument("--relation-gains", choices=RELATION_GAINS, default="per_relation",
                        help="linear-response encoder: a learned gain per relation, or one gain shared by every relation (the relation-typing ablation)")
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
    parser.add_argument("--bce-in-log-space", action="store_true",
                        help="training and validation loss from the head's log P and log(1 - P) instead of P clamped to [1e-6, 1 - 1e-6], under which a "
                             "pair beyond the clamp has a constant loss and no gradient (a sigmoid logit below -13.8 for a positive); off by default, "
                             "which reproduces every run before 8 October 2026")
    parser.add_argument("--positive-target-from-frequency", action="store_true",
                        help="open question 8: use the reported HPO or label frequency (floored at --minimum-frequency-target) as the target of a positive pair instead of --positive-target")
    parser.add_argument("--minimum-frequency-target", type=float, default=0.05)
    parser.add_argument("--negative-weight", type=float, default=0.2)
    parser.add_argument("--description-length-coefficient", type=float, default=1e-6)
    parser.add_argument("--checkpoint-every-minutes", type=float, default=20.0)
    parser.add_argument("--time-budget-seconds", type=float, default=0.0, help="stop training after this many seconds (0 = no limit)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--timing-batches", type=int, default=0,
                        help="run this many training steps on the training perturbations, print seconds per step and per epoch, and exit: no split directory, validation or test score")
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    return parser


def main() -> None:
    arguments = build_argument_parser().parse_args()
    signal.signal(signal.SIGUSR1, request_checkpoint)
    torch.manual_seed(arguments.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by, label_grades=tuple(arguments.label_grades) if arguments.label_grades else None,
                                label_selection=arguments.label_selection)
    if data.label_selection_summary is not None:
        print(f"label selection: {data.label_selection_summary}")
        if arguments.time_split_cutoff is not None:
            raise ValueError("--label-selection is not implemented with --time-split-cutoff")
    if arguments.permute_labels and arguments.time_split_cutoff is not None:
        raise ValueError("--permute-labels is not defined with --time-split-cutoff (the time split permutes its new positives itself)")
    in_lockbox, lockbox_summary = None, None
    if arguments.score_lockbox and arguments.lockbox is None:
        raise ValueError("--score-lockbox needs --lockbox")
    if arguments.lockbox is not None:
        if arguments.time_split_cutoff is not None or (arguments.score_lockbox and (arguments.holdout_module or arguments.holdout_subsystem)):
            raise ValueError("--score-lockbox is defined for the grouped split only, and --lockbox not with the time split; "
                             "a module or subsystem hold-out runs on the development perturbations (--lockbox without --score-lockbox)")
        in_lockbox = read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir)
        lockbox_summary = {"path": str(arguments.lockbox), "sha256": file_sha256(arguments.lockbox), "num_lockbox_perturbations": int(in_lockbox.sum()),
                           "role": "scored" if arguments.score_lockbox else "removed before the split"}
        print(f"lockbox: {lockbox_summary}")
        if not arguments.score_lockbox:
            data = restrict_to_perturbations(data, ~in_lockbox)
            in_lockbox = None
    permutation_source_rows_sha256 = None
    if arguments.permute_labels:  # rows move whole, so a pair keeps its weight, frequency and label-mask entry
        source_row = degree_stratified_row_permutation(data.perturbation_degrees, random_seed=arguments.seed, partition=in_lockbox)
        permutation_source_rows_sha256 = array_sha256(source_row)  # experiments/score_confirmatory.py checks it against its own draw
        data.outcomes, data.weights = data.outcomes[source_row], data.weights[source_row]
        if data.frequencies is not None:
            data.frequencies = data.frequencies[source_row]
        if data.label_mask is not None:
            data.label_mask = data.label_mask[source_row]
    label_mask = data.label_mask  # None without --label-selection: every pair is labelled
    if arguments.encoder == "linear_response":  # its input is sign x magnitude, so a perturbation whose seeds all have sign 0 gives a zero field
        zero_input = [perturbation_id for perturbation_id, signs, magnitudes in zip(data.perturbation_ids, data.perturbation_signs, data.perturbation_magnitudes)
                      if not np.any(np.asarray(signs, dtype=float) * np.asarray(magnitudes, dtype=float))]
        if zero_input:
            print(f"warning: {len(zero_input)} perturbations have no seed with a nonzero sign and magnitude; the linear-response encoder gives them "
                  f"the prediction of no perturbation: {zero_input}")
    rewiring_summary = None
    if arguments.rewire_swaps_per_edge > 0:
        original_edges = np.stack([data.edge_source, data.edge_target])
        undirected_relations = reciprocated_relations(original_edges, data.edge_relation) if arguments.keep_reciprocated_relations_symmetric else []
        rewired = fast_degree_preserving_rewiring(original_edges, data.edge_relation, num_swaps_per_edge=arguments.rewire_swaps_per_edge, random_seed=arguments.seed,
                                                  undirected_relations=undirected_relations)
        rewiring_summary = {"swaps_per_edge": arguments.rewire_swaps_per_edge, "seed": arguments.seed, "num_edges": int(rewired.shape[1]),
                            "share_of_edges_unchanged": float((rewired == original_edges).all(axis=0).mean()),
                            "duplicate_edges_before": duplicate_edge_count(original_edges, data.edge_relation), "duplicate_edges_after": duplicate_edge_count(rewired, data.edge_relation)}
        if arguments.keep_reciprocated_relations_symmetric:  # only then, so the split signature of earlier rewired runs is unchanged
            rewiring_summary["undirected_relations"] = [data.relation_types[relation] for relation in undirected_relations]
        data.edge_source, data.edge_target = rewired[0], rewired[1]  # before the encoder and the adjacencies are built
        print(f"rewired graph: {rewiring_summary}")
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
    train_indices, validation_indices, test_indices, split_name = split_indices(data, arguments, in_lockbox)
    if arguments.refit_on_validation and (len(validation_indices) == 0 or arguments.laboratory_label_weight > 0 or arguments.time_split_cutoff is not None):
        raise ValueError("--refit-on-validation needs a validation set and is not implemented with the laboratory readout or the time split")
    if arguments.permute_labels:
        split_name += "_permuted"
    if arguments.rewire_swaps_per_edge > 0:
        split_name += "_rewired"
    encoder, head = build_models(data, arguments, device)
    initialise_leaks_from_base_rates(head, data, train_indices, label_mask, arguments, device)
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
    if arguments.timing_batches > 0:
        time_training_steps(data, encoder, head, optimizer, train_indices, validation_indices, adjacencies, label_mask, arguments, device)
        return
    split_directory = arguments.run_dir / split_name
    split_signature = {"test_perturbation_ids_sha256": hashlib.sha256("\n".join(data.perturbation_ids[i] for i in test_indices).encode()).hexdigest(),
                       "train_perturbation_ids_sha256": hashlib.sha256("\n".join(data.perturbation_ids[i] for i in train_indices).encode()).hexdigest(),
                       "lockbox_sha256": (lockbox_summary or {}).get("sha256"), "rewiring": rewiring_summary, "labels_permuted": bool(arguments.permute_labels)}
    refit_after_done = False
    if arguments.resume and (split_directory / "DONE").exists():
        finished = json.loads((split_directory / "results.json").read_text()) if (split_directory / "results.json").exists() else {}
        if "test_perturbation_ids" in finished and finished["test_perturbation_ids"] != [data.perturbation_ids[i] for i in test_indices]:
            raise SystemExit(f"{split_directory} is DONE for other test perturbations (another lockbox or grouping); use another --run-dir")
        if not (arguments.refit_on_validation and finished.get("refit") is None and (split_directory / "checkpoint.pt").exists()):
            print(f"{split_name}: DONE marker present; skipping (delete the marker to retrain)")
            return
        # an early-stopping run finished before --refit-on-validation was added: keep its results and train only the refit
        refit_after_done = True
        if not (split_directory / "results_early_stopped.json").exists():
            (split_directory / "results_early_stopped.json").write_text((split_directory / "results.json").read_text())
        print(f"{split_name}: DONE without a refit; refitting from the stored early-stopping checkpoint (its results kept in results_early_stopped.json)")
    split_directory.mkdir(parents=True, exist_ok=True)
    checkpoint_path = split_directory / "checkpoint.pt"
    state = {"epoch": 0, "history": [], "best_validation_auprc": float("-inf"), "best_validation_loss": float("inf"), "best_epoch": -1, "epochs_without_improvement": 0, "best_encoder": None, "best_head": None, "code_provenance": []}
    checkpoint_torch_rng_state = None
    fingerprint = configuration_fingerprint(arguments, data, label_mask)
    configuration_changes = []
    if arguments.resume and checkpoint_path.exists():
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        # a persistent buffer (the message-passing node features) comes back from the checkpoint, not from the current inputs
        buffers_from_current_inputs = {name: buffer.detach().clone() for name, buffer in encoder.named_buffers() if name in checkpoint["encoder"]}
        encoder.load_state_dict(checkpoint["encoder"])
        for name, buffer in encoder.named_buffers():
            if name in buffers_from_current_inputs and (buffer.shape != buffers_from_current_inputs[name].shape or not torch.equal(buffer, buffers_from_current_inputs[name])):
                print(f"warning: the checkpoint's encoder buffer {name} differs from the one built from the current inputs; "
                      "the run continues with the checkpoint's values, while results.json records the current files' hashes")
        head.load_state_dict(checkpoint["head"])
        if laboratory_readout is not None and "laboratory_readout" in checkpoint:
            laboratory_readout.load_state_dict(checkpoint["laboratory_readout"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        state = checkpoint["state"]
        checkpoint_torch_rng_state = checkpoint.get("torch_rng_state")  # absent in checkpoints written before 8 October 2026
        if state.get("split_signature", split_signature) != split_signature:
            raise SystemExit(f"{checkpoint_path} was written for another split (signature {state['split_signature']}); use another --run-dir")
        previous_fingerprint = state.get("configuration_fingerprint")  # absent in checkpoints written before 8 October 2026
        if previous_fingerprint is not None:
            configuration_changes = sorted(key for key in set(previous_fingerprint) | set(fingerprint) if previous_fingerprint.get(key) != fingerprint.get(key))
            if configuration_changes:
                print(f"warning: {checkpoint_path} was written under other arguments or inputs; resuming anyway. Changed: {configuration_changes}")
        if refit_after_done:
            state.update(training_finished=True, stopped_early=bool(finished.get("stopped_early")))
        print(f"resumed at epoch {state['epoch']}")
    if "code_provenance" not in state:  # a checkpoint written before commits were recorded: its epochs ran under unrecorded code
        state["code_provenance"] = [{"started_at_epoch": 0, "commit": None, "tracked_changes": None}] if state["epoch"] > 0 else []
    provenance = git_provenance()
    if provenance["tracked_changes"]:
        print(f"warning: tracked files differ from commit {provenance['commit']}: {provenance['tracked_changes']}")
    state["split_signature"] = split_signature
    state["configuration_fingerprint"] = fingerprint
    state["code_provenance"].append({"started_at_epoch": state["epoch"], **provenance, **({"configuration_changes": configuration_changes} if configuration_changes else {})})  # one entry per process, so a resumed run lists every commit it ran under

    def save_checkpoint() -> None:
        atomic_torch_save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "state": state,
                           "torch_rng_state": torch.get_rng_state(),  # the noisy-OR gate noise draws from it
                           **({"laboratory_readout": laboratory_readout.state_dict()} if laboratory_readout is not None else {})}, checkpoint_path)

    outcomes = torch.as_tensor(data.outcomes, dtype=torch.float32)
    weights = torch.as_tensor(np.where(data.outcomes > 0, np.maximum(data.weights, 1e-3), arguments.negative_weight), dtype=torch.float32)
    if label_mask is not None:
        weights = weights * torch.as_tensor(label_mask, dtype=torch.float32)  # a pair set aside by the selection is neither positive nor negative
    positive_targets = None
    if arguments.positive_target_from_frequency:
        frequency = np.where(np.isnan(data.frequencies), arguments.positive_target, np.maximum(data.frequencies, arguments.minimum_frequency_target))
        positive_targets = torch.as_tensor(np.where(data.outcomes > 0, frequency, 0.0), dtype=torch.float32)
        weights = torch.as_tensor(np.where(data.outcomes > 0, 1.0, arguments.negative_weight), dtype=torch.float32)  # the frequency is the target, not the weight
        if label_mask is not None:
            weights = weights * torch.as_tensor(label_mask, dtype=torch.float32)
    node_cost = float(np.log(len(data.node_ids)))

    def training_step(model_encoder, model_head, model_laboratory_readout, model_optimizer, batch) -> tuple[float, float, float, float | None]:
        """One optimizer step on a batch; returns the loss, its BCE and penalty parts and the laboratory loss (None without it)."""
        node_index, sign_and_magnitude = pad_perturbations(data, batch)
        field = encode(model_encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
        output = model_head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
        log_parts = ({"log_probability": output.symptom_log_probability, "log_complement": output.symptom_log_complement}
                     if arguments.bce_in_log_space else {})
        if positive_targets is not None:
            bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, positive_targets[batch].to(device), weights[batch].to(device), positive_target=1.0, **log_parts)
        else:
            bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device), positive_target=arguments.positive_target, **log_parts)
        penalty = arguments.description_length_coefficient * model_head.description_length_penalty(node_cost=node_cost)
        loss = bce + penalty
        laboratory_value = None
        if model_laboratory_readout is not None and laboratory_index.count(batch):
            positions, nodes, slots, directions = laboratory_index.batch(batch, device)
            laboratory_loss = laboratory_sign_loss(model_laboratory_readout.label_scores(field, positions, nodes, slots, len(directions)), directions)
            loss = loss + arguments.laboratory_label_weight * laboratory_loss
            laboratory_value = laboratory_loss.item()
        model_optimizer.zero_grad()
        loss.backward()
        model_optimizer.step()
        return loss.item(), bce.item(), float(penalty), laboratory_value

    started = time.time()
    last_checkpoint = time.time()
    # A resumed run continues the epoch-order generator and the torch generator where the checkpoint left them, and a
    # mid-epoch checkpoint finishes its epoch with the same order from the next batch, so a resumed run trains as an
    # uninterrupted one does. Checkpoints written before 8 October 2026 hold neither: they reseed the order with
    # seed + epoch and restart a mid-epoch checkpoint's epoch from its first batch, as that code did.
    generator = np.random.default_rng(arguments.seed + state["epoch"])
    if state.get("epoch_order_generator_state") is not None:
        generator.bit_generator.state = state["epoch_order_generator_state"]
    if checkpoint_torch_rng_state is not None:
        torch.set_rng_state(checkpoint_torch_rng_state)
    stopped_early = bool(state.get("stopped_early", False))
    time_budget_reached = False
    # a checkpoint written after early stopping ended (training_finished) is not trained further on resume
    for epoch in range(arguments.max_epochs if state.get("training_finished") else state["epoch"], arguments.max_epochs):
        in_progress = state.pop("epoch_in_progress", None)
        if in_progress is not None and in_progress["epoch"] == epoch:
            order = np.asarray(in_progress["order"], dtype=int)
            first_batch_start = in_progress["next_batch_start"]
            epoch_loss, epoch_bce, epoch_penalty, epoch_laboratory = in_progress["sums"]
            print(f"continuing epoch {epoch} at batch start {first_batch_start}")
        else:
            order = generator.permutation(train_indices)
            first_batch_start = 0
            epoch_loss, epoch_bce, epoch_penalty, epoch_laboratory = 0.0, 0.0, 0.0, 0.0
        for start in range(first_batch_start, len(order), arguments.batch_size):
            batch = order[start : start + arguments.batch_size]
            loss_value, bce_value, penalty_value, laboratory_value = training_step(encoder, head, laboratory_readout, optimizer, batch)
            if laboratory_value is not None:
                epoch_laboratory += laboratory_value * len(batch)
            epoch_loss += loss_value * len(batch)
            epoch_bce += bce_value * len(batch)
            epoch_penalty += penalty_value * len(batch)
            if CHECKPOINT_REQUESTED or (time.time() - last_checkpoint) / 60.0 >= arguments.checkpoint_every_minutes:
                state["epoch_in_progress"] = {"epoch": epoch, "order": order.tolist(), "next_batch_start": start + arguments.batch_size,
                                              "sums": [epoch_loss, epoch_bce, epoch_penalty, epoch_laboratory]}
                state["epoch_order_generator_state"] = generator.bit_generator.state
                save_checkpoint()
                state.pop("epoch_in_progress")
                last_checkpoint = time.time()
                if CHECKPOINT_REQUESTED:
                    print("checkpoint written on signal; exiting for requeue")
                    return
        entry = {"epoch": epoch, "train_loss": epoch_loss / max(1, len(order)), "train_bce": epoch_bce / max(1, len(order)), "train_penalty": epoch_penalty / max(1, len(order)), "elapsed_seconds": time.time() - started}
        if laboratory_readout is not None:
            entry["train_laboratory_loss"] = epoch_laboratory / max(1, len(order))
        if len(validation_indices):
            validation_log_parts = {}
            if arguments.bce_in_log_space:
                validation_predictions, validation_log_probability, validation_log_complement = predict(encoder, head, data, validation_indices, adjacencies, arguments, device, with_log_parts=True)
                validation_log_parts = {"log_probability": torch.as_tensor(validation_log_probability, dtype=torch.float32),
                                        "log_complement": torch.as_tensor(validation_log_complement, dtype=torch.float32)}
            else:
                validation_predictions = predict(encoder, head, data, validation_indices, adjacencies, arguments, device)
            entry["validation_macro_auprc"] = macro_auprc(validation_predictions, data.outcomes[validation_indices], None if label_mask is None else label_mask[validation_indices])
            # the validation loss takes the training loss's targets: with --positive-target-from-frequency each positive's
            # target is its frequency (it took 0.99 for every positive before 8 October, with the frequency mode's weights)
            validation_targets, validation_positive_target = ((positive_targets[validation_indices], 1.0) if positive_targets is not None
                                                              else (outcomes[validation_indices], arguments.positive_target))
            entry["validation_loss"] = float(evidence_weighted_binary_cross_entropy(torch.as_tensor(validation_predictions, dtype=torch.float32), validation_targets, weights[validation_indices],
                                                                                    positive_target=validation_positive_target, **validation_log_parts))
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
        state["epoch_order_generator_state"] = generator.bit_generator.state
        print(f"epoch {epoch} loss {entry['train_loss']:.4f} (bce {entry['train_bce']:.4f}, penalty {entry['train_penalty']:.4f})" + (f" validation loss {entry['validation_loss']:.4f} macro AUPRC {entry['validation_macro_auprc']:.3f}" if "validation_loss" in entry else "") + f" elapsed {entry['elapsed_seconds']:.0f}s", flush=True)
        save_checkpoint()
        if len(validation_indices) and state["epochs_without_improvement"] >= arguments.patience:
            stopped_early = True
            print(f"early stop at epoch {epoch}; best validation epoch {state['best_epoch']}")
            break
        if arguments.time_budget_seconds and time.time() - started > arguments.time_budget_seconds:
            print("time budget reached")
            time_budget_reached = True
            break
    # a resume after early stopping (a crash while scoring, or during the refit) must not train the early-stopping model
    # further; a run that used up --max-epochs stays resumable with a larger --max-epochs, except under the refit
    if stopped_early or time_budget_reached or arguments.refit_on_validation:
        state["training_finished"] = True
    state["stopped_early"] = stopped_early
    save_checkpoint()
    if arguments.refit_on_validation:
        (split_directory / "DONE").unlink(missing_ok=True)  # written again when the refit has scored the test
    if state["best_encoder"] is not None:
        encoder.load_state_dict(state["best_encoder"])
        head.load_state_dict(state["best_head"])
        if laboratory_readout is not None and state.get("best_laboratory_readout") is not None:
            laboratory_readout.load_state_dict(state["best_laboratory_readout"])

    refit_summary = None
    if arguments.refit_on_validation:
        # Goodfellow, Bengio and Courville, Deep Learning (2016), section 7.8, algorithm 7.2: initialise again and train on
        # all the training data for the number of epochs early stopping chose; same epochs, not same optimizer steps
        early_stopped_predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
        np.save(split_directory / "test_predictions_early_stopped.npy", early_stopped_predictions)
        early_stopped_test_mask = None if label_mask is None else label_mask[test_indices]
        early_stopped_scores = {"macro_auprc": macro_auprc(early_stopped_predictions, data.outcomes[test_indices], early_stopped_test_mask),
                                "micro_auprc": micro_auprc(early_stopped_predictions, data.outcomes[test_indices], early_stopped_test_mask)}
        if state["best_epoch"] < 0:
            raise ValueError("no validation epoch improved (the validation loss was never finite), so the refit has no number of epochs to train for")
        refit_indices = np.sort(np.concatenate([train_indices, validation_indices]))
        refit_checkpoint_path = split_directory / "refit_checkpoint.pt"
        torch.manual_seed(arguments.seed)
        encoder, head = build_models(data, arguments, device)
        initialise_leaks_from_base_rates(head, data, refit_indices, label_mask, arguments, device)
        optimizer = torch.optim.AdamW(optimizer_parameter_groups(encoder, head, arguments), lr=arguments.learning_rate, weight_decay=arguments.weight_decay)
        refit_state = {"epoch": 0, "next_batch_start": 0, "epoch_sums": [0.0, 0.0, 0.0], "history": [], "epochs": state["best_epoch"] + 1,
                       "num_perturbations": int(len(refit_indices)), "code_provenance": []}
        if arguments.resume and refit_checkpoint_path.exists():
            refit_checkpoint = torch.load(refit_checkpoint_path, map_location=device, weights_only=False)
            encoder.load_state_dict(refit_checkpoint["encoder"])
            head.load_state_dict(refit_checkpoint["head"])
            optimizer.load_state_dict(refit_checkpoint["optimizer"])
            refit_state = refit_checkpoint["refit_state"]
            torch.set_rng_state(refit_checkpoint["torch_rng_state"])
            print(f"refit resumed at epoch {refit_state['epoch']}, batch start {refit_state['next_batch_start']}")
        refit_state["code_provenance"].append({"started_at_epoch": refit_state["epoch"], **provenance})

        def save_refit_checkpoint() -> None:
            atomic_torch_save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "refit_state": refit_state,
                               "torch_rng_state": torch.get_rng_state()}, refit_checkpoint_path)

        print(f"refit on {len(refit_indices)} training and validation perturbations for {refit_state['epochs']} epochs (best validation epoch {state['best_epoch']})")
        last_checkpoint = time.time()
        for epoch in range(refit_state["epoch"], refit_state["epochs"]):
            order = np.random.default_rng(arguments.seed + 100_000 + epoch).permutation(refit_indices)  # per epoch, so a resume sees the same order
            for start in range(refit_state["next_batch_start"], len(order), arguments.batch_size):
                batch = order[start : start + arguments.batch_size]
                loss_value, bce_value, penalty_value, _ = training_step(encoder, head, None, optimizer, batch)
                refit_state["epoch_sums"] = [refit_state["epoch_sums"][0] + loss_value * len(batch), refit_state["epoch_sums"][1] + bce_value * len(batch),
                                             refit_state["epoch_sums"][2] + penalty_value * len(batch)]
                refit_state["next_batch_start"] = start + arguments.batch_size
                if CHECKPOINT_REQUESTED or (time.time() - last_checkpoint) / 60.0 >= arguments.checkpoint_every_minutes:
                    save_refit_checkpoint()
                    last_checkpoint = time.time()
                    if CHECKPOINT_REQUESTED:
                        print("refit checkpoint written on signal; exiting for requeue")
                        return
            loss_sum, bce_sum, penalty_sum = refit_state["epoch_sums"]
            refit_state["history"].append({"epoch": epoch, "train_loss": loss_sum / len(order), "train_bce": bce_sum / len(order), "train_penalty": penalty_sum / len(order)})
            refit_state.update(epoch=epoch + 1, next_batch_start=0, epoch_sums=[0.0, 0.0, 0.0])
            save_refit_checkpoint()
            print(f"refit epoch {epoch} loss {refit_state['history'][-1]['train_loss']:.4f}", flush=True)
        refit_summary = {"epochs": refit_state["epochs"], "num_perturbations": refit_state["num_perturbations"], "history": refit_state["history"],
                         "code_provenance": refit_state["code_provenance"], "early_stopped_test_predictions": "test_predictions_early_stopped.npy",
                         "early_stopped_scores": early_stopped_scores}
        if not arguments.score_lockbox:
            print(f"early-stopped model on the test: macro AUPRC {early_stopped_scores['macro_auprc']:.3f}, micro AUPRC {early_stopped_scores['micro_auprc']:.3f}")

    predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
    test_outcomes = data.outcomes[test_indices]
    test_mask = None if label_mask is None else label_mask[test_indices]
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
        if not scorable_symptom(test_outcomes, symptom_index, MINIMUM_POSITIVES_TO_SCORE, test_mask):
            continue
        labelled_rows = scored_rows(len(test_outcomes), symptom_index, test_mask)
        per_symptom[symptom] = {
            "positives": int(test_outcomes[labelled_rows, symptom_index].sum()),
            "auprc": bootstrap_interval(lambda p, y, m=None: per_symptom_auprc(p, y, symptom_index, m), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap, mask=test_mask).__dict__,
            "auroc": bootstrap_interval(lambda p, y, m=None: per_symptom_auroc(p, y, symptom_index, m), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap, mask=test_mask).__dict__,
            "base_rate": float(test_outcomes[labelled_rows, symptom_index].mean()),
            **({"masked_pairs": int((~labelled_rows).sum())} if test_mask is not None else {}),
        }
    results = {
        "split": split_name, "fold": None if (arguments.holdout_module or arguments.holdout_subsystem or arguments.score_lockbox) else arguments.fold, "holdout_module": arguments.holdout_module or None,
        "holdout_subsystem": arguments.holdout_subsystem or None, "labels_permuted": bool(arguments.permute_labels), "seed": arguments.seed,
        "arguments": {key: (value if isinstance(value, (int, float, str, bool, list, type(None))) else str(value)) for key, value in vars(arguments).items()},
        "num_train": int(len(train_indices)), "num_validation": int(len(validation_indices)), "num_test": int(len(test_indices)),
        "epochs_completed": state["epoch"], "best_epoch": state["best_epoch"], "stopped_early": stopped_early, "history": state["history"],
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "micro_auprc": micro_auprc(predictions, test_outcomes, test_mask) if time_split is None else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, test_outcomes, test_mask), "hits_at_3": hits_at_k(predictions, test_outcomes, 3, test_mask),
        "expected_calibration_error": expected_calibration_error(predictions, test_outcomes, mask=test_mask), "per_symptom": per_symptom, "symptoms": data.symptoms,
        "test_perturbation_ids": [data.perturbation_ids[i] for i in test_indices],
        "time_split": time_split_results,
        "code_provenance": state["code_provenance"],
        "graph_files_sha256": {name: file_sha256(Path(arguments.graph_dir) / f"{name}.parquet") for name in ("nodes", "edges")},
        "evidence_records_sha256": file_sha256(Path(arguments.evidence_dir) / "evidence_records.parquet"),
        "permutation_source_rows_sha256": permutation_source_rows_sha256,
        "configuration_fingerprint": state.get("configuration_fingerprint"),
        "node_descriptors_sha256": file_sha256(getattr(arguments, "node_descriptors", None)),
        "cell_class_weights_sha256": file_sha256(getattr(arguments, "cell_class_weights", None)),
        "label_selection": data.label_selection_summary, "label_selection_sha256": file_sha256(getattr(arguments, "label_selection", None)),
        "laboratory_labels": laboratory_results,
        "lockbox": lockbox_summary, "rewiring": rewiring_summary, "refit": refit_summary,
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
    if arguments.score_lockbox:  # the lockbox is read once, by experiments/score_confirmatory.py; logs read by check-ins must not show it
        print(f"{split_name}: lockbox predictions written to {split_directory} (scores not printed)")
    else:
        print(f"{split_name}: macro AUPRC {results['macro_auprc']:.3f} macro AUROC {results['macro_auroc']:.3f} MRR {results['mean_reciprocal_rank']:.3f} hits@3 {results['hits_at_3']:.3f} ECE {results['expected_calibration_error']:.3f}")
    (split_directory / "DONE").write_text("done\n")


if __name__ == "__main__":
    main()
