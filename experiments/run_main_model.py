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
import hashlib
import json
import os
import signal
import subprocess
import time
from collections import Counter
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from mechanistic_pathway_learning.evaluation.experiment_data import (
    DRUG_ENTRY_MODES,
    load_experiment_data,
    read_lockbox,
    restrict_to_perturbations,
)
from mechanistic_pathway_learning.evaluation.negative_controls import (
    degree_stratified_row_permutation,
    duplicate_edge_count,
    fast_degree_preserving_rewiring,
    reciprocated_relations,
)
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
    macro_auprc_by_degree_bin,
    mean_reciprocal_rank,
    micro_auprc,
    per_symptom_auprc,
    per_symptom_auroc,
    scorable_symptom,
    scored_rows,
)
from mechanistic_pathway_learning.evidence.allele_mechanisms import gene_of_perturbation
from mechanistic_pathway_learning.evidence.calibrated_reliability import (
    MINIMUM_CALIBRATION_CASES,
    calibrated_probability_matrix,
    cell_table,
    feature_frame,
    fit_calibration,
    permute_cell_rows,
    summarise_calibration,
)
from mechanistic_pathway_learning.evidence.soft_priors import (
    DEFAULT_COUNT_CAP,
    DEFAULT_RECORD_GRADES,
    DEFAULT_TEMPERING_EXPONENT,
    DEFAULT_TIERS_PATH,
    apply_soft_priors,
    eligible_cells,
    literature_popularity,
    permute_soft_priors,
    read_tiers,
    resolve_pairs,
    soft_evidence_from_records,
    soft_evidence_from_reports,
    soft_prior_matrices,
    summarise as summarise_soft_priors,
)
from mechanistic_pathway_learning.evidence.combined_perturbations import (
    COMBINATION_PERTURBATION_TYPE,
    DRUG_ROW_KINDS,
    PERTURBATION_TYPES_NEVER_SCORED,
)
from mechanistic_pathway_learning.graph.brain_expression_weights import ALL_CELLS_CLASS, EXTRACELLULAR_COMPARTMENT, MEMBER_RELATION
from mechanistic_pathway_learning.graph.cofactor_edges import CARRIER_RULES, cofactor_edge_mask
from mechanistic_pathway_learning.graph.compartment_ph import compartment_ph_inputs
from mechanistic_pathway_learning.graph.drug_entry_nodes import CHANGE_CARRYING_RELATIONS, ENTRY_RELATIONS, SEQUESTERS_RELATION
from mechanistic_pathway_learning.graph.edge_evidence import (
    DEFAULT_COUNT_CAP as DEFAULT_EDGE_COUNT_CAP,
    DEFAULT_LABEL_MINIMUM_HUMAN_RECORDS,
    DEFAULT_LABEL_MINIMUM_RECORDS,
    DEFAULT_TEMPERING as DEFAULT_EDGE_TEMPERING,
    DEFAULT_TIERS_PATH as DEFAULT_EDGE_TIERS_PATH,
    curated_edge_labels,
    edge_priors,
    permute_edge_priors,
    read_edge_tiers,
    summarise as summarise_edge_priors,
)
from mechanistic_pathway_learning.graph.gene_protein_split import ENCODES_RELATION
from mechanistic_pathway_learning.graph.node_descriptors import DESCRIPTOR_BLOCK_PREFIXES, columns_a_seed_still_reads, descriptor_blocks
from mechanistic_pathway_learning.models.baselines.local_descriptor_encoder import LocalDescriptorEncoder
from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
from mechanistic_pathway_learning.models.baselines.zero_field_encoder import ZeroFieldEncoder
from mechanistic_pathway_learning.models.descriptor_treatments import DESCRIPTOR_TREATMENTS
from mechanistic_pathway_learning.models.edge_gate_passes import (
    DEFAULT_SURROGATE_WIDTH,
    binary_pass_step,
    closure_evidence,
    estimated_closure_evidence,
    exact_flip_checks,
)
from mechanistic_pathway_learning.models.edge_gates import EDGE_GATE_KINDS, HARD_CONCRETE_BETA, EdgeGates, HardConcreteEdgeGates, penalty_denominator
from mechanistic_pathway_learning.models.laboratory_readout import LaboratoryLabelIndex, LaboratoryReadout, laboratory_sign_loss
from mechanistic_pathway_learning.models.linear_response_encoder import (
    CLASS_WEIGHT_RULES,
    CROSS_RELATION_AGGREGATORS,
    DEFAULT_FIXED_POINT_TOLERANCE,
    DEFAULT_MAX_PROPAGATION_STEPS,
    EDGE_SIGNS,
    INPUT_WEIGHT_FORMS,
    MIXTURE_WEIGHTINGS,
    NORMALISATIONS,
    RELATION_GAINS,
    LinearResponseEncoder,
)
from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
from mechanistic_pathway_learning.models.relational_message_passing_encoder import (
    DEFAULT_RELATION_FACTOR_RANK,
    RELATION_PARAMETERISATIONS,
    CONJUNCTION_AGGREGATIONS,
    DEFAULT_CONJUNCTION_RELATIONS,
    DEFAULT_SOFT_MINIMUM_TEMPERATURE,
    RelationalMessagePassingEncoder,
)
from mechanistic_pathway_learning.models.soft_constraint_losses import (
    PROBABILITY_EPSILON,
    equal_drug_loss_shares,
    evidence_weighted_binary_cross_entropy,
)

MESSAGE_PASSING_EDGE_SIGNS: tuple[str, ...] = ("graph", "unsigned")
TRAINED_RELATION = "induces"  # the relation load_experiment_data reads by default and the soft priors are drawn on
CHECKPOINT_REQUESTED = False
MINIMUM_POSITIVES_TO_SCORE = 5


def request_checkpoint(signal_number, frame) -> None:
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
    return {"labels": len(scores), "reached": int(reached.sum()),
            "sign_agreement": float(np.mean(np.sign(scores[reached]) == directions[reached])) if reached.any() else None,
            "majority_direction_rate": float(np.mean(directions == 1)),
            "auroc": float(roc_auc_score(directions == 1, scores)) if len(set(directions.tolist())) == 2 else None}


def perturbation_covariate(data) -> np.ndarray:
    """Standardised log1p of the summed degree of each perturbation's nodes: the hub signal that popularity and the
    random walk exploit, given to the heads explicitly under --degree-offset. It is a property of the input graph, not
    of the labels, so standardising over every perturbation leaks nothing."""
    log_degree = np.log1p(data.perturbation_degrees_for_strata)
    return ((log_degree - log_degree.mean()) / max(log_degree.std(), 1e-8)).astype(np.float32)


def covariate_of(data, batch: np.ndarray, arguments, device):
    if not arguments.degree_offset:
        return None
    return torch.as_tensor(perturbation_covariate(data)[batch], device=device)


def descriptor_column_names(data, arguments) -> list[str]:
    """The descriptor table's columns in the order node_feature_matrix appends them, after --drop-descriptor-blocks."""
    if not getattr(arguments, "node_descriptors", None):
        return []
    descriptors = data.node_descriptor_table if getattr(data, "node_descriptor_table", None) is not None else pd.read_parquet(arguments.node_descriptors)
    dropped_blocks = getattr(arguments, "drop_descriptor_blocks", None) or []
    if dropped_blocks:
        blocks = descriptor_blocks(descriptors.columns)
        descriptors = descriptors.drop(columns=[column for block in dropped_blocks for column in blocks[block]])
    return [str(column) for column in descriptors.columns]


def node_feature_matrix(data, arguments) -> np.ndarray:
    """The fixed structural features, followed by the node descriptors of --node-descriptors when given (one block per
    node type, zero outside it: mechanistic_pathway_learning/graph/node_descriptors.py)."""
    structural = data.structural_node_features()
    if not getattr(arguments, "node_descriptors", None):
        return structural
    descriptors = data.node_descriptor_table if getattr(data, "node_descriptor_table", None) is not None else pd.read_parquet(arguments.node_descriptors)
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


# The nine relations amendment item 8 names, measured on development training pairs: 2,100 to 2,468 pairs reach each,
# against 451 or fewer for the other eleven.
REGISTERED_FACTORISED_RELATIONS = ("encodes", "catalyzed_by", "substrate_of", "activates", "product_of", "binds",
                                   "regulates_transcription_of", "inhibits", "member_of")


def relation_parameterisation_inputs(data, arguments) -> dict:
    """The encoder keywords of --relation-parameterisation (amendment item 8). A named relation the graph does not have
    is left out and reported, rather than shifting the indices of the ones it does have."""
    if getattr(arguments, "relation_parameterisation", "per_layer_dense") != "shared_low_rank":
        return {}
    named = list(getattr(arguments, "factorised_relations", REGISTERED_FACTORISED_RELATIONS) or ())
    present = [relation for relation in named if relation in set(data.relation_types)]
    missing = [relation for relation in named if relation not in set(data.relation_types)]
    indices = tuple(list(data.relation_types).index(relation) for relation in present)
    parameters = len(data.relation_types) * arguments.node_state_dim + len(indices) * 2 * arguments.node_state_dim * arguments.relation_factor_rank
    print(f"relation transforms: one shared set at every round, {len(data.relation_types)} diagonals and "
          f"{len(present)} rank-{arguments.relation_factor_rank} factorisations, {parameters:,} parameters"
          + (f"; named but not in this graph: {missing}" if missing else ""))
    inputs = {"relation_parameterisation": "shared_low_rank", "factorised_relation_indices": indices,
              "relation_factor_rank": arguments.relation_factor_rank}
    if arguments.fixed_point_tolerance is not None:  # amendment item 7, which item 8's shared weights made possible
        inputs["fixed_point_tolerance"] = arguments.fixed_point_tolerance
        inputs["max_message_passing_rounds"] = arguments.max_propagation_steps
        print(f"message passing runs to a fixed point at tolerance {arguments.fixed_point_tolerance}, "
              f"at most {arguments.max_propagation_steps} rounds")
    return inputs


def descriptor_treatment_inputs(data, arguments) -> dict:
    """The encoder keywords of --descriptor-treatment: the treatment and the number of descriptor columns at the end of
    node_feature_matrix."""
    resolve_registered_defaults(arguments)
    treatment = arguments.descriptor_treatment
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
    inputs = {"descriptor_treatment": treatment, "num_descriptor_columns": num_descriptor_columns}
    if treatment == "seed_masked":
        kept = columns_a_seed_still_reads(descriptor_column_names(data, arguments))
        if len(kept) != num_descriptor_columns:
            raise ValueError(f"{sum(1 for _ in kept)} descriptor names against {num_descriptor_columns} descriptor columns")
        inputs["descriptor_columns_kept"] = torch.as_tensor(kept)
        print(f"seed masking leaves {sum(kept)} of {num_descriptor_columns} descriptor columns readable at a perturbed "
              f"node: the perturbed gene's regional and whole-tissue expression (the owner, 13 October 2026)")
        partner_pairs = seed_mask_partner_pairs(data)
        if partner_pairs.shape[1]:
            inputs["seed_mask_partner_index"] = torch.as_tensor(partner_pairs, dtype=torch.long)
    return inputs


def seed_mask_partner_pairs(data) -> np.ndarray:
    """[2, num_pairs] of (node, partner): a perturbation seeded on the node has its partner masked with it under
    seed_masked. On a split graph the partners of a seed are its encodes partners, both directions. A drug entry node's
    partners are the targets of its mechanism edges and their encodes partners, one direction only, so a drug seeded on
    its node hides the descriptors of the nodes it hides when it is seeded on its targets, and a seed on a target does
    not reach back to the drugs that act on it."""
    nodes_of_pairs, partners_of_pairs = [], []
    encodes_partners: dict[int, list[int]] = {}
    if ENCODES_RELATION in data.relation_types:  # a split graph: the mask covers each seed's encodes partners
        encodes = np.asarray(data.edge_relation) == data.relation_types.index(ENCODES_RELATION)
        genes, proteins = np.asarray(data.edge_source)[encodes], np.asarray(data.edge_target)[encodes]
        nodes_of_pairs += [genes, proteins]
        partners_of_pairs += [proteins, genes]
        for gene, protein in zip(genes.tolist(), proteins.tolist()):
            encodes_partners.setdefault(gene, []).append(protein)
            encodes_partners.setdefault(protein, []).append(gene)
        print(f"seed masking covers each seed's encodes partners ({int(encodes.sum())} encodes edges)")
    if getattr(data, "entry_node_mask", None) is not None:
        is_entry_edge, _, _, _ = entry_edge_arrays(data)
        entry_pairs = set()
        for entry_node, target in zip(np.asarray(data.edge_source)[is_entry_edge].tolist(), np.asarray(data.edge_target)[is_entry_edge].tolist()):
            entry_pairs.add((entry_node, target))
            entry_pairs.update((entry_node, partner) for partner in encodes_partners.get(target, []))
        ordered_entry_pairs = np.array(sorted(entry_pairs), dtype=int).reshape(-1, 2)
        nodes_of_pairs.append(ordered_entry_pairs[:, 0])
        partners_of_pairs.append(ordered_entry_pairs[:, 1])
        print(f"seed masking covers each drug node's targets and their encodes partners ({len(ordered_entry_pairs)} pairs)")
    if not nodes_of_pairs:
        return np.zeros((2, 0), dtype=int)
    return np.stack([np.concatenate(nodes_of_pairs), np.concatenate(partners_of_pairs)])


def cell_class_inputs(data, arguments) -> tuple[torch.Tensor | None, torch.Tensor | None, torch.Tensor | None]:
    """Per-node class weights [N, K] of --cell-class-weights (experiments/build_cell_class_weights.py), and, with
    --extracellular-coupling, the extracellular metabolites as one pool shared by every class but the all-cells one
    (that class stays the propagation without cell classes)."""
    if not getattr(arguments, "cell_class_weights", None):
        if getattr(arguments, "extracellular_coupling", False):
            raise ValueError("--extracellular-coupling needs --cell-class-weights")
        return None, None, None
    table = data.cell_class_weight_table if getattr(data, "cell_class_weight_table", None) is not None else pd.read_parquet(arguments.cell_class_weights)
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


def recompute_reaction_expression(data, edges_before: np.ndarray, edges_after: np.ndarray, arguments) -> dict:
    """--recompute-reaction-expression-on-rewired-graph: rewrite each reaction's gene rule onto its rewired catalysts,
    recompute the reaction_brain descriptor block and the reactions' cell-class weights from the rewritten rules, and
    recompute each protein complex's cell-class gate from its rewired members (graph/rewired_expression_features.py),
    so a rewired run's node inputs follow its wiring (the user's decision of 8 October 2026, extended to the complex
    gates by decision 4.2b on 13 October: "Rewire the complex gates for that control"). Gene and metabolite rows and
    the intrinsic blocks stay as they are.

    The flag keeps its name although it now reaches the complexes too, because experiments/
    evaluate_on_rewired_graph.py rebuilds a finished run from its recorded arguments and a rename would break that.
    What a rebuild of a run made before 13 October does change is its complex gates, which is the point of 4.2b and
    is recorded in the returned counts."""
    from mechanistic_pathway_learning.graph.rewired_expression_features import (
        rewired_complex_cell_class_weights,
        rewired_gene_reaction_rules,
        rewired_reaction_brain_block,
        rewired_reaction_cell_class_weights,
    )
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    if list(nodes.node_id) != list(data.node_ids):
        raise ValueError(f"{arguments.graph_dir}/nodes.parquet is not in the node order of the loaded graph")
    catalysis = np.asarray(data.edge_relation) == data.relation_types.index("catalyzed_by")
    rules, counts = rewired_gene_reaction_rules(nodes, edges_before, edges_after, catalysis)
    record = {**counts, "gene_expression_table_sha256": file_sha256(arguments.gene_expression_table), "class_expression_table_sha256": None}
    if getattr(arguments, "node_descriptors", None):
        descriptors = pd.read_parquet(arguments.node_descriptors).copy()
        block = rewired_reaction_brain_block(nodes, rules, pd.read_parquet(arguments.gene_expression_table))
        missing = [column for column in block.columns if column not in descriptors.columns]
        if missing:
            raise ValueError(f"{arguments.node_descriptors} has no column {missing[0]}; it was not built with the brain expression blocks")
        descriptors.loc[block.index, block.columns] = block.to_numpy()
        data.node_descriptor_table = descriptors
    if getattr(arguments, "cell_class_weights", None):
        weights = pd.read_parquet(arguments.cell_class_weights).copy()
        class_expression = pd.read_parquet(arguments.class_expression_table)
        carries_all_cells = ALL_CELLS_CLASS in weights.columns
        recomputed = rewired_reaction_cell_class_weights(nodes, rules, class_expression, include_all_cells_class=carries_all_cells)
        weights.loc[recomputed.index, weights.columns] = recomputed[weights.columns].to_numpy()
        membership = np.asarray(data.edge_relation) == data.relation_types.index(MEMBER_RELATION) \
            if MEMBER_RELATION in data.relation_types else np.zeros(len(data.edge_relation), dtype=bool)
        complex_weights, complex_counts = rewired_complex_cell_class_weights(
            nodes, edges_after, membership, class_expression, include_all_cells_class=carries_all_cells)
        weights.loc[complex_weights.index, weights.columns] = complex_weights[weights.columns].to_numpy()
        record.update(complex_counts)
        data.cell_class_weight_table = weights
        record["class_expression_table_sha256"] = file_sha256(arguments.class_expression_table)
    return record


# The registered v4 model is the trainer's default (the owner, 13 October 2026: "all of the v4 stuff should be ON by
# default. the other models were shit quite frankly because they didn't even follow the physiology very well").
# experiments/run_main_model_batch.py states the same model as explicit flags (V4_MESSAGE_PASSING_MODEL,
# V4_LINEAR_RESPONSE_MODEL) and tests/test_registered_defaults.py holds the two equal, so a default cannot drift from the
# registration without a test failing. A run made before a default changed is rebuilt (experiments/
# evaluate_on_rewired_graph.py) with the value it was trained under, which is the value below for any argument its
# results.json does not record, because the argument did not exist yet.
VALUES_BEFORE_THE_V4_DEFAULTS = {
    "group_by": "gene", "refit_on_validation": False, "keep_large_groups_in_training": False, "node_features": "identity",
    "descriptor_treatment": "plain", "extracellular_coupling": False, "class_weight_rule": "product", "input_weights": "signed",
    "fixed_point_tolerance": None, "cofactor_relations": False, "num_layers": 2, "start_at_weighted_optimum": False,
    "module_bias_init": 0.0, "gate_init_noise": 0.01, "link_learning_rate": 0.0, "leak_learning_rate": 0.0, "gate_learning_rate": 0.0,
    "placebo_negative_weight": None, "placebo_negative_weights_file": None, "equal_drug_shares": False, "normalisation": "in_degree", "message_passing_edge_signs": "unsigned", "combined_perturbations": None,
    "compartment_ph": False, "relation_parameterisation": "per_layer_dense", "reversible_reactions": False,
    "soft_prior_strength": 0.0, "soft_prior_reports": [], "permute_soft_priors": False, "soft_prior_probability": "tier",
    "edge_prior_strength": 0.0, "permute_edge_priors": False,
}


def path_or_none(value: str) -> Path | None:
    """A file path, or 'none' (or the empty string) for no file."""
    return None if value.strip().lower() in ("", "none") else Path(value)


def tolerance_or_none(value: str) -> float | None:
    """--fixed-point-tolerance: a positive number, or 'none' for a fixed number of steps."""
    if value.strip().lower() == "none":
        return None
    tolerance = float(value)
    if tolerance <= 0:
        raise argparse.ArgumentTypeError("the tolerance must be positive; use 'none' for a fixed number of steps")
    return tolerance


def resolve_registered_defaults(arguments) -> None:
    """Fill the arguments whose registered value depends on the rest of the run (left at None by the parser). Runs in
    place and is idempotent, so a caller that builds a model without main() reaches the same model."""
    linear_response = getattr(arguments, "encoder", None) == "linear_response"
    if getattr(arguments, "descriptor_treatment", None) is None:
        masks = getattr(arguments, "node_descriptors", None) and getattr(arguments, "encoder", None) in ("message_passing", "linear_response")
        arguments.descriptor_treatment = "seed_masked" if masks else "plain"
    if getattr(arguments, "cofactor_relations", None) is None:
        arguments.cofactor_relations = linear_response
    if getattr(arguments, "compartment_ph", None) is None:
        arguments.compartment_ph = linear_response
    if getattr(arguments, "reversible_reactions", None) is None:
        arguments.reversible_reactions = linear_response
    if getattr(arguments, "extracellular_coupling", None) is None:
        arguments.extracellular_coupling = linear_response and bool(getattr(arguments, "cell_class_weights", None))
    if getattr(arguments, "refit_on_validation", None) is None:
        not_implemented = getattr(arguments, "laboratory_label_weight", 0) > 0 or getattr(arguments, "time_split_cutoff", None) is not None
        arguments.refit_on_validation = not not_implemented
        if not_implemented:
            print("--refit-on-validation is off: it is not implemented with the laboratory readout or the time split")


# arguments a resumed process may change without changing what the checkpoint was trained for
RESUME_CONTROL_ARGUMENTS = {"run_dir", "resume", "refit_on_validation", "max_epochs", "checkpoint_every_minutes", "time_budget_seconds", "timing_batches", "num_bootstrap"}
# flags added after runs had started: off, they leave the fingerprint as it was, so a resumed run reports no change
ARGUMENTS_RECORDED_ONLY_WHEN_SET = {"start_at_weighted_optimum", "rewire_encodes", "equal_drug_shares", "normalise_drug_input", "drug_entry",
                                    "drug_mechanism_before_first_layer", "sequestration_carries", "permute_drug_targets", "placebo_negative_weight",
                                    "soft_prior_strength", "edge_prior_strength"}
# the settings of the soft prior mean nothing while its strength is 0, so they enter the fingerprint only when it is on
SOFT_PRIOR_SETTINGS = {"soft_prior_tiers", "soft_prior_tempering", "soft_prior_count_cap", "soft_prior_record_grades", "soft_prior_reports", "permute_soft_priors",
                       "soft_prior_probability", "soft_prior_calibration_minimum_cases"}
# and the edge priors' settings mean nothing while their strength is 0 (decision 9.57), nor do those of the gates they
# make (9.62) and of the curated edges taken as known (9.63); until 17 October 2026 the gates' settings of 9.62 were
# recorded whatever the strength, so a run without edge priors resumed under them reported changes that were not
EDGE_PRIOR_SETTINGS = {"edge_prior_tiers", "edge_prior_tempering", "edge_prior_count_cap", "permute_edge_priors", "edge_gate_learning_rate",
                       "edge_gate_kind", "edge_gate_beta", "edge_gate_binary_epochs", "edge_gate_binary_learning_rate", "edge_gate_surrogate_width",
                       "edge_gate_binary_full_batch", "edge_gate_closure_evidence", "edge_gate_flip_checks", "edge_gate_binary_learning_rate_decay",
                       "curated_edge_labels", "curated_label_minimum_records", "curated_label_minimum_human_records", "curated_label_minimum_sign_records",
                       "curated_label_exclude_disputed_signs", "curated_label_bulk_papers"}
EDGE_EVIDENCE_FILE = "edge_evidence.parquet"  # experiments/build_edge_evidence.py writes it beside the graph's edges.parquet
EDGE_GATED_ENCODERS: tuple[str, ...] = ("linear_response", "message_passing")


def json_ready(value):
    """An argument value as results.json records it: basic types as they are, a list element by element, a path or
    anything else as its string."""
    if isinstance(value, (int, float, str, bool, type(None))):
        return value
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return str(value)


def array_sha256(values) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    return hashlib.sha256(str(array.dtype).encode() + str(array.shape).encode() + array.tobytes()).hexdigest()


def rewiring_fixed_relations(relation_types: list[str], rewire_encodes: bool, has_entry_nodes: bool = False) -> list[int]:
    """Relation indices the rewiring leaves as they are: encodes on a split graph unless --rewire-encodes (an open
    decision of the user, docs/gene_protein_split.md); none on a merged graph. With drug entry nodes, also a drug's
    mechanism edges, its carriers' sequestration edges and its enzymes' clearance edges: they are how the perturbation enters, as its seeds are
    without the nodes, and seeds are not rewired. A rewiring that moved them would give a drug other targets, which
    tests something else than whether the graph beyond the targets carries the signal."""
    fixed = [relation_types.index(ENCODES_RELATION)] if ENCODES_RELATION in relation_types and not rewire_encodes else []
    if has_entry_nodes:
        fixed += [relation_types.index(relation) for relation in (*ENTRY_RELATIONS, *CHANGE_CARRYING_RELATIONS) if relation in relation_types]
    return fixed


def entry_node_inputs(data) -> dict:
    """The encoder keywords that say which nodes are drug entry nodes; empty on a graph without them or under
    --drug-entry targets, so every earlier run builds the encoder it always built."""
    if getattr(data, "entry_node_mask", None) is None:
        return {}
    return {"entry_node_mask": torch.as_tensor(data.entry_node_mask)}


SEQUESTRATION_MODES = ("change", "presence")


def sequestration_inputs(data, arguments) -> dict:
    """The message-passing keywords for the edges that end on a drug and say "more of the source, less of the drug": a
    carrier's sequestration edges and, on a graph that holds them, an enzyme's clearance edges (drug_entry_nodes.py,
    CHANGE_CARRYING_RELATIONS). Under "change" (the default on a graph that holds such a relation) an edge carries how
    far its source has moved from its unperturbed state, so it is silent in a perturbation that leaves the carrier or
    the enzyme at rest; under "presence" it is averaged like any other relation and carries the source's state itself,
    which gives every drug with a known carrier the same extra input (docs/drug_entry_nodes_measured.md, section 6).
    Empty on a graph without the relations. With sequestration alone the keywords are the ones they always were."""
    present = [relation for relation in CHANGE_CARRYING_RELATIONS if relation in data.relation_types]
    if not present:
        return {}
    relation_indices = [data.relation_types.index(relation) for relation in present]
    edge_relation = np.asarray(data.edge_relation)
    carries_change = np.isin(edge_relation, relation_indices)
    mode = getattr(arguments, "sequestration_carries", None) or "change"
    counts = ", ".join(f"{relation} {int((edge_relation == index).sum())}" for relation, index in zip(present, relation_indices))
    print(f"edges from a carrier or a clearing enzyme to a drug ({counts}), carrying the source's {'change from its unperturbed state' if mode == 'change' else 'state (presence)'}")
    if mode == "presence" or not carries_change.any():
        return {}
    edge_index = torch.as_tensor(np.stack([np.asarray(data.edge_source)[carries_change], np.asarray(data.edge_target)[carries_change]]), dtype=torch.long)
    if present == [SEQUESTERS_RELATION]:
        return {"change_only_edge_index": edge_index, "change_only_relation_index": relation_indices[0]}
    return {"change_only_edge_index": edge_index, "change_only_edge_relation": torch.as_tensor(edge_relation[carries_change], dtype=torch.long)}


def entry_edge_arrays(data) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(is entry edge, relation indices of the entry relations, edge signs, edge weights) of a graph with entry nodes."""
    entry_relation_indices = np.array([data.relation_types.index(relation) for relation in ENTRY_RELATIONS if relation in data.relation_types], dtype=int)
    is_entry_edge = np.isin(np.asarray(data.edge_relation), entry_relation_indices)
    weights = np.ones(len(is_entry_edge)) if data.edge_weight is None else np.asarray(data.edge_weight, dtype=float)
    return is_entry_edge, entry_relation_indices, np.asarray(data.edge_sign, dtype=float), weights


def perturbations_with_a_zero_signed_input(data) -> list[str]:
    """Perturbations whose input carries no sign: every seed has sign 0 or magnitude 0, or, for a drug seeded on its
    entry node, every mechanism edge of that node does. The linear-response encoder gives them a zero field."""
    signed_input_of_entry_node = None
    if getattr(data, "entry_node_mask", None) is not None:
        is_entry_edge, _, signs, weights = entry_edge_arrays(data)
        signed_input_of_entry_node = np.zeros(len(data.node_ids))
        np.add.at(signed_input_of_entry_node, np.asarray(data.edge_source)[is_entry_edge], np.abs(signs[is_entry_edge] * weights[is_entry_edge]))
    zero_input = []
    for perturbation_id, seeds, signs, magnitudes in zip(data.perturbation_ids, data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes):
        seed_input = np.abs(np.asarray(signs, dtype=float) * np.asarray(magnitudes, dtype=float))
        if signed_input_of_entry_node is not None:  # a seed on an entry node passes on only what its mechanism edges carry
            on_entry_node = np.asarray(data.entry_node_mask)[np.asarray(seeds, dtype=int)]
            seed_input = np.where(on_entry_node, seed_input * signed_input_of_entry_node[np.asarray(seeds, dtype=int)], seed_input)
        if not np.any(seed_input):
            zero_input.append(perturbation_id)
    return zero_input


def configuration_fingerprint(arguments, data, label_mask) -> dict:
    """The training arguments and digests of the inputs as loaded, recorded in the checkpoint so that a resume under
    other arguments or rebuilt inputs is reported (the split signature covers only the split and the controls)."""
    soft_priors_on = bool(getattr(arguments, "soft_prior_strength", 0.0))
    edge_priors_on = bool(getattr(arguments, "edge_prior_strength", 0.0))
    fingerprint = {f"argument:{name}": str(value) for name, value in sorted(vars(arguments).items())
                   if name not in RESUME_CONTROL_ARGUMENTS and not (name in ARGUMENTS_RECORDED_ONLY_WHEN_SET and not value)
                   and not (name in SOFT_PRIOR_SETTINGS and not soft_priors_on) and not (name in EDGE_PRIOR_SETTINGS and not edge_priors_on)}
    fingerprint.update({
        "input:outcomes": array_sha256(data.outcomes), "input:weights": array_sha256(data.weights),
        "input:label_mask": None if label_mask is None else array_sha256(label_mask),
        "input:edges": array_sha256(np.stack([np.asarray(data.edge_source), np.asarray(data.edge_target), np.asarray(data.edge_relation)])),
        "input:seeds": hashlib.sha256(repr([(list(map(int, seeds)), list(map(float, signs)), list(map(float, magnitudes)))
                                            for seeds, signs, magnitudes in zip(data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes)]).encode()).hexdigest()})
    for name in ("node_descriptors", "cell_class_weights", "laboratory_labels"):
        value = getattr(arguments, name, None)
        fingerprint[f"file:{name}"] = file_sha256(value) if value and Path(value).is_file() else None
    if soft_priors_on:
        fingerprint["file:soft_prior_tiers"] = file_sha256(arguments.soft_prior_tiers)
        for path in getattr(arguments, "soft_prior_reports", None) or []:
            fingerprint[f"file:soft_prior_reports:{Path(path).name}"] = file_sha256(path) if Path(path).is_file() else None
    if edge_priors_on:
        fingerprint["file:edge_prior_tiers"] = file_sha256(arguments.edge_prior_tiers)
        evidence_path = Path(arguments.graph_dir) / EDGE_EVIDENCE_FILE
        fingerprint["file:edge_evidence"] = file_sha256(evidence_path) if evidence_path.is_file() else None
    return fingerprint


def file_sha256(path) -> str | None:
    if not path:
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compartment_ph_encoder_inputs(data) -> dict:
    """compartment_ph_inputs as the linear-response encoder takes them, the index arrays as tensors. Until 16 October
    2026 the numpy arrays went to the encoder as they were, and its constructor failed (".long()" on an array) on every
    graph that names a proton node in a compartment, so the registered model, which turns compartment pH on, could not
    be built on the real graph; the test graphs have no compartments and passed."""
    return {name: torch.as_tensor(value, dtype=torch.long) if isinstance(value, np.ndarray) else value for name, value in compartment_ph_inputs(data).items()}


def build_models(data, arguments, device):
    """The encoder and the head. With edge priors on the data (build_edge_priors), the linear-response encoder takes its
    edge gates here; the message-passing encoder takes them once its adjacencies exist (attach_edge_gates)."""
    resolve_registered_defaults(arguments)
    if getattr(data, "edge_prior_gated", None) is not None and arguments.encoder not in EDGE_GATED_ENCODERS:
        raise ValueError(f"the edge priors gate the edges an encoder propagates along; the {arguments.encoder} encoder reads none (use one of {EDGE_GATED_ENCODERS})")
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
                                        reversible_reactions=(torch.as_tensor(data.is_reversible_reaction)
                                                              if arguments.reversible_reactions and data.is_reversible_reaction is not None else None),
                                        num_propagation_steps=arguments.propagation_steps,
                                        fixed_point_tolerance=arguments.fixed_point_tolerance,
                                        max_propagation_steps=arguments.max_propagation_steps,
                                class_weight_rule=arguments.class_weight_rule, input_weights=arguments.input_weights,
                                        **(compartment_ph_encoder_inputs(data) if arguments.compartment_ph else {}),
                                        propagation_channels=arguments.propagation_channels, response_scale=arguments.response_scale, damping=arguments.propagation_damping,
                                        normalisation=arguments.normalisation, cross_relation_aggregator=arguments.cross_relation_aggregator,
                                        mixture_weighting=arguments.mixture_weighting, node_type_index=node_type_index,
                                        edge_signs=arguments.edge_signs, relation_gains=arguments.relation_gains,
                                        cell_class_weights=cell_class_weights, channels_per_cell_class=arguments.channels_per_cell_class,
                                        extracellular_pool_nodes=pool_nodes, shared_pool_classes=shared_pool_classes,
                                        **descriptor_treatment_inputs(data, arguments),
                                        **({"entry_relation_names": ENTRY_RELATIONS, "edge_weight": torch.as_tensor(data.edge_weight), **entry_node_inputs(data)}
                                           if entry_node_inputs(data) else {}),
                                        edge_gates=edge_gates_of(data, arguments)).to(device)
        # Item 20, built 14 October 2026: the stoichiometric demand of the limiting-step pass. Before this
        # bottleneck_edge_demand was None and every demand was 1, so a reaction consuming two molecules counted as one.
        # A graph built before the column existed carries no coefficient and keeps the old behaviour.
        stoichiometric_demand_counts = None
        if data.edge_stoichiometric_coefficient is not None and arguments.stoichiometric_demand:
            stoichiometric_demand_counts = encoder.set_bottleneck_edge_demand_from_edges(
                torch.as_tensor(data.edge_source), torch.as_tensor(data.edge_target),
                torch.as_tensor(data.edge_relation), torch.as_tensor(data.edge_stoichiometric_coefficient))
            print(f"stoichiometric demand: {stoichiometric_demand_counts['entries_with_a_stated_coefficient']} of "
                  f"{stoichiometric_demand_counts['bottleneck_entries']} bottleneck entries carry a stated coefficient, "
                  f"{stoichiometric_demand_counts['entries_above_one']} above 1 (largest "
                  f"{stoichiometric_demand_counts['largest_demand']:.0f}), "
                  f"{stoichiometric_demand_counts['entries_clamped_up_from_a_fraction']} clamped up from a fraction")
        elif data.edge_stoichiometric_coefficient is None:
            print("stoichiometric demand: the graph states no coefficients, so every demand stays 1 (item 20 unbuilt "
                  "for this graph; it needs the v4 rebuild)")
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
        conjunction_relations = [relation for relation in arguments.conjunction_relations if relation in set(data.relation_types)]
        missing_relations = [relation for relation in arguments.conjunction_relations if relation not in set(data.relation_types)]
        conjunction_indices = tuple(list(data.relation_types).index(relation) for relation in conjunction_relations)
        entry_inputs = entry_node_inputs(data)
        if entry_inputs:
            is_entry_edge, entry_relation_indices, edge_signs, edge_weights = entry_edge_arrays(data)
            entry_inputs.update(entry_edge_index=torch.as_tensor(np.stack([np.asarray(data.edge_source)[is_entry_edge], np.asarray(data.edge_target)[is_entry_edge]]), dtype=torch.long),
                                entry_edge_sign=torch.as_tensor(edge_signs[is_entry_edge], dtype=torch.float32),
                                entry_edge_weight=torch.as_tensor(edge_weights[is_entry_edge], dtype=torch.float32),
                                entry_relation_indices=tuple(int(index) for index in entry_relation_indices),
                                entry_before_first_layer=bool(getattr(arguments, "drug_mechanism_before_first_layer", False)))
            print(f"drug entry nodes: {int(np.asarray(data.entry_node_mask).sum())} nodes, {int(is_entry_edge.sum())} mechanism edges, present only where seeded"
                  + (", mechanism edges also delivered before the first layer" if entry_inputs["entry_before_first_layer"] else ""))
        entry_inputs.update(sequestration_inputs(data, arguments))
        encoder = RelationalMessagePassingEncoder(len(data.node_ids), len(data.relation_types), arguments.node_state_dim, arguments.num_layers, node_features=node_features,
                                                  conjunction_relation_indices=conjunction_indices,
                                                  conjunction_aggregation=arguments.conjunction_aggregation,
                                                  soft_minimum_temperature=arguments.soft_minimum_temperature,
                                                  **relation_parameterisation_inputs(data, arguments),
                                                  **descriptor_treatment_inputs(data, arguments), **entry_inputs).to(device)
        if arguments.conjunction_aggregation == "mean":
            print("conjunction relations: none; every relation is aggregated by its in-degree mean")
        elif conjunction_indices:
            edges_of_relation = int(np.isin(data.edge_relation, np.array(conjunction_indices)).sum())
            print(f"conjunction relations ({arguments.conjunction_aggregation}"
                  + (f", temperature {arguments.soft_minimum_temperature:g}" if arguments.conjunction_aggregation == "soft_minimum" else "")
                  + f"): {', '.join(conjunction_relations)}, {edges_of_relation} edges"
                  + (f"; absent from this graph: {', '.join(missing_relations)}" if missing_relations else ""))
        else:
            print(f"conjunction relations: none of {', '.join(arguments.conjunction_relations)} is in this graph, "
                  "so every relation is aggregated by its in-degree mean")
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


def run_edge_gate_binary_pass(encoder, head, gates: HardConcreteEdgeGates, fit_rows: np.ndarray, fitted: dict, data, adjacencies, arguments, device,
                              edge_penalty_denominator: float) -> dict:
    """The binary pass over hard edge gates (decision 9.62; models/edge_gate_passes.py), after the fit and before the
    test is scored: the gates alone, every other parameter frozen, the encoder and head in evaluation mode so every
    forward reads the binary read-out, and the data loss the fit minimised (the weighted BCE of the final fit's targets,
    without the head's description-length penalty, which no gate moves). Then, with --edge-gate-flip-checks, the exact
    checks of the most promising flips."""
    batches = [fit_rows[start : start + arguments.batch_size] for start in range(0, len(fit_rows), arguments.batch_size)]

    def data_loss(batch) -> torch.Tensor:
        node_index, sign_and_magnitude = pad_perturbations(data, batch)
        field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
        output = head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
        log_parts = ({"log_probability": output.symptom_log_probability, "log_complement": output.symptom_log_complement} if arguments.bce_in_log_space else {})
        if fitted["positive_targets"] is not None:
            return evidence_weighted_binary_cross_entropy(output.symptom_probability, fitted["positive_targets"][batch].to(device), fitted["weights"][batch].to(device),
                                                          positive_target=1.0, **log_parts)
        return evidence_weighted_binary_cross_entropy(output.symptom_probability, fitted["outcomes"][batch].to(device), fitted["weights"][batch].to(device),
                                                      positive_target=arguments.positive_target, **log_parts)

    frozen = [parameter for parameter in list(encoder.parameters()) + list(head.parameters()) if parameter is not gates.log_alpha and parameter.requires_grad]
    for parameter in frozen:
        parameter.requires_grad_(False)
    encoder.eval()
    head.eval()
    gates.set_phase("binary")
    started = time.time()
    closed_before = int((gates.read_out() < 0.5).sum())
    optimizer = torch.optim.Adam([gates.log_alpha], lr=arguments.edge_gate_binary_learning_rate)
    learning_rate_decay = float(getattr(arguments, "edge_gate_binary_learning_rate_decay", 1.0))  # absent from a run recorded before 17 October 2026
    history = []
    try:
        for epoch in range(arguments.edge_gate_binary_epochs):
            flips = 0
            if arguments.edge_gate_binary_full_batch:
                flips += binary_pass_step(encoder, gates, data_loss, None, edge_penalty_denominator, optimizer, arguments.edge_gate_surrogate_width, batches=batches)["flipped"]
            else:
                for batch in batches:
                    flips += binary_pass_step(encoder, gates, data_loss, batch, edge_penalty_denominator, optimizer, arguments.edge_gate_surrogate_width)["flipped"]
            history.append({"epoch": epoch, "flips": flips, "closed_at_read_out": int((gates.read_out() < 0.5).sum()), "elapsed_seconds": time.time() - started,
                            "learning_rate": optimizer.param_groups[0]["lr"]})
            print(f"edge-gate binary pass epoch {epoch}: {flips} flips, {history[-1]['closed_at_read_out']} closed", flush=True)
            for group in optimizer.param_groups:  # 1 keeps the rate; below 1 anneals it, so gates near the threshold stop crossing it
                group["lr"] *= learning_rate_decay
        checks = None
        if arguments.edge_gate_flip_checks > 0:
            checks = exact_flip_checks(encoder, gates, data_loss, batches, edge_penalty_denominator, arguments.edge_gate_flip_checks, arguments.edge_gate_surrogate_width)
        evidence = None
        budget = int(getattr(arguments, "edge_gate_closure_evidence", 0))
        if budget != 0:  # -1: every closed edge (decision 9.63 (c)); N: the N against the strongest priors
            evidence = closure_evidence(encoder, gates, data_loss, batches, edge_penalty_denominator, None if budget < 0 else budget,
                                        arguments.edge_gate_surrogate_width)
        # every closure, with its first-order evidence from one more pass, for the review list (decision 9.63 (c))
        estimated = estimated_closure_evidence(encoder, gates, data_loss, batches, edge_penalty_denominator)
        closed_edges = [{"edge": edge, **{key: record[key] for key in ("prior", "estimated_evidence_nats", "estimated_posterior_false")}}
                        for edge, record in estimated.items()]
    finally:
        for parameter in frozen:
            parameter.requires_grad_(True)
    return {"fit_rows": int(len(fit_rows)), "epochs": arguments.edge_gate_binary_epochs, "learning_rate": arguments.edge_gate_binary_learning_rate,
            "learning_rate_decay": learning_rate_decay, "surrogate_width": arguments.edge_gate_surrogate_width, "full_batch": bool(arguments.edge_gate_binary_full_batch),
            "closed_before": closed_before, "closed_after": int((gates.read_out() < 0.5).sum()), "history": history, "flip_checks": checks,
            "closure_evidence": evidence, "closed_edges": closed_edges, "seconds": time.time() - started}


CLOSURE_EVIDENCE_TABLE_COLUMNS = ["family", "kind", "records", "contested", "species_known_records", "human_records", "animal_only_records", "pubmed_ids",
                                  "pubmed_species", "resources"]
CLOSURE_REVIEW_COLUMNS = ["edge", "source_id", "target_id", "relation", "prior", "evidence_nats", "posterior_false", "estimated_evidence_nats",
                          "estimated_posterior_false", *CLOSURE_EVIDENCE_TABLE_COLUMNS]


def write_closure_review(binary_pass: dict, data, directory: Path) -> dict:
    """The review list of decision 9.63 (c): every edge the fit closed ("Fit decides": it stays closed), each with what
    the build knows of it, written to edge_gate_closures.csv, against the strongest priors first. Per edge: its graph
    edge (source, target, relation), its prior, its exact closure evidence in nats and P(false) where measured (empty
    where --edge-gate-closure-evidence left it out), the first-order estimate of both for every closure (one pass for
    all; models/edge_gate_passes.estimated_closure_evidence), its family and kind, its independent records with their PubMed ids,
    how many iCite scored, how many are human evidence and how many animal-only, each paper's species scores, and the
    resources that state the edge. Returns the summary results.json keeps in place of the rows."""
    measured = {record["edge"]: record for record in (binary_pass.get("closure_evidence") or {}).get("records", [])}
    evidence_table = getattr(data, "edge_evidence", None)
    rows = []
    for closed in sorted(binary_pass["closed_edges"], key=lambda entry: (-entry["prior"], entry["edge"])):
        edge = closed["edge"]
        record = measured.get(edge, {})
        row = {"edge": edge, "source_id": str(data.node_ids[int(data.edge_source[edge])]), "target_id": str(data.node_ids[int(data.edge_target[edge])]),
               "relation": str(data.relation_types[int(data.edge_relation[edge])]), "prior": closed["prior"],
               "evidence_nats": record.get("evidence_nats"), "posterior_false": record.get("posterior_false"),
               "estimated_evidence_nats": closed.get("estimated_evidence_nats"), "estimated_posterior_false": closed.get("estimated_posterior_false")}
        if evidence_table is not None:
            for column in CLOSURE_EVIDENCE_TABLE_COLUMNS:
                if column in evidence_table.columns:
                    value = evidence_table[column].iat[edge]
                    row[column] = value.item() if hasattr(value, "item") else value
        rows.append(row)
    table = pd.DataFrame(rows, columns=CLOSURE_REVIEW_COLUMNS)
    table.to_csv(directory / "edge_gate_closures.csv", index=False)
    evidences = pd.to_numeric(table.evidence_nats, errors="coerce")
    posteriors = pd.to_numeric(table.posterior_false, errors="coerce")
    return {"closed": len(rows), "evaluated": int(evidences.notna().sum()), "file": "edge_gate_closures.csv",
            "evidence_at_least_1_nat": int((evidences >= 1.0).sum()), "posterior_false_at_least_0_9": int((posteriors >= 0.9).sum()),
            "by_kind": ({f"{family}:{kind}": int(count) for (family, kind), count in table.groupby(["family", "kind"]).size().items()}
                        if evidence_table is not None and len(table) else {})}


def optimizer_parameter_groups(encoder, head, arguments, extra_parameters=()) -> list[dict]:
    """One group at the main learning rate, plus one group for each noisy-OR time scale (links, leaks, module biases,
    gates) whose learning rate is set; Adam moves a parameter by about one learning rate per step, so the rate is the
    time scale on which that parameter can change. extra_parameters (the laboratory readout) join the main group. Under
    --descriptor-treatment zero_init_slow the encoder's descriptor map is one more group, at --descriptor-learning-rate
    with the main weight decay. With edge priors the edge gates' logits are a group of their own, at
    --edge-gate-learning-rate (the main rate when unset) and without weight decay: decay pulls a logit toward 0, every
    gate toward 0.5, which is a prior of its own that the evidence does not state."""
    rates = {"links": arguments.link_learning_rate, "leaks": arguments.leak_learning_rate, "module_biases": arguments.module_bias_learning_rate,
             "gates": getattr(arguments, "gate_learning_rate", 0.0)}
    descriptor_parameters = encoder.descriptor_parameters() if hasattr(encoder, "descriptor_parameters") else []
    descriptor_groups = [{"params": descriptor_parameters, "lr": arguments.descriptor_learning_rate}] if descriptor_parameters else []
    edge_gate_parameters = list(encoder.edge_gates.parameters()) if getattr(encoder, "edge_gates", None) is not None else []
    edge_gate_groups = ([{"params": edge_gate_parameters, "lr": getattr(arguments, "edge_gate_learning_rate", None) or arguments.learning_rate, "weight_decay": 0.0}]
                        if edge_gate_parameters else [])
    grouped_apart = {id(parameter) for parameter in descriptor_parameters + edge_gate_parameters}
    encoder_parameters = [parameter for parameter in encoder.parameters() if id(parameter) not in grouped_apart]
    if not hasattr(head, "time_scale_parameter_groups") or not any(rates.values()):
        return [{"params": encoder_parameters + list(head.parameters()) + list(extra_parameters)}, *descriptor_groups, *edge_gate_groups]
    separate_groups = [{"params": parameters, "lr": rates[name], "weight_decay": 0.0}
                       for name, parameters in head.time_scale_parameter_groups().items() if rates[name]]
    separated_ids = {id(parameter) for group in separate_groups for parameter in group["params"]}
    other_parameters = encoder_parameters + [parameter for parameter in head.parameters() if id(parameter) not in separated_ids] + list(extra_parameters)
    return [{"params": other_parameters}, *separate_groups, *descriptor_groups, *edge_gate_groups]


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


REGISTERED_PLACEBO_NEGATIVE_WEIGHTS = Path("configs/placebo_negative_weights.csv")
PLACEBO_WEIGHT_FILE_COLUMNS = ("perturbation_id", "symptom", "weight")


def placebo_negative_weight_matrix(path, perturbation_ids: list[str], symptoms: list[str]) -> tuple[np.ndarray, dict]:
    """The per-pair placebo negative weights of amendment item 3, as a [perturbation, symptom] array of NaN except
    where the file names the pair, with a record of what matched.

    The file (configs/placebo_negative_weights.csv) holds one row per never-above-placebo pair: weight
    n_effective / (n_effective + 214) where the drug's label states both arm sizes, and the registered default of 0.5
    where it does not. A row naming a perturbation or a symptom this data set does not hold is counted rather than
    applied, because the file is written against one evidence table and a rebuild may drop a drug; the counts go into
    the run record so a file that has gone stale is visible instead of silently weighing fewer pairs than registered.
    A file that matches no pair at all is refused, since that is the wrong file rather than a stale one."""
    table = pd.read_csv(path)
    missing_columns = [column for column in PLACEBO_WEIGHT_FILE_COLUMNS if column not in table.columns]
    if missing_columns:
        raise ValueError(f"{path} has no column {missing_columns[0]}; it is not a placebo negative weight file")
    below_zero = table[table.weight < 0]
    if len(below_zero):
        raise ValueError(f"{path} gives {below_zero.iloc[0].perturbation_id} / {below_zero.iloc[0].symptom} a weight of "
                         f"{below_zero.iloc[0].weight}; a loss weight below zero inverts the gradient for that pair")
    row_of_perturbation = {perturbation_id: row for row, perturbation_id in enumerate(perturbation_ids)}
    column_of_symptom = {symptom: column for column, symptom in enumerate(symptoms)}
    matrix = np.full((len(perturbation_ids), len(symptoms)), np.nan)
    matched, unmatched = 0, []
    for entry in table.itertuples():
        row, column = row_of_perturbation.get(entry.perturbation_id), column_of_symptom.get(entry.symptom)
        if row is None or column is None:
            unmatched.append(f"{entry.perturbation_id}/{entry.symptom}")
            continue
        matrix[row, column] = float(entry.weight)
        matched += 1
    if matched == 0:
        raise ValueError(f"{path} names {len(table)} pairs and none is in this evidence table; it is the wrong file")
    return matrix, {"file": str(path), "sha256": file_sha256(path), "pairs_in_file": int(len(table)),
                    "pairs_matched": matched, "pairs_not_in_this_evidence_table": unmatched}


def placebo_weight_pairs_in_force(matrix: np.ndarray, outcomes: np.ndarray, evidenced_negative) -> int:
    """Pairs a per-pair placebo weight actually acts on: named by the file, negative, and marked evidenced negative
    by the selection. Matching an id and a symptom is not enough, because the weight reaches only the pairs the
    selection relabelled; a selection that stops relabelling a pair leaves the file matching 61 of 61 while no
    weight applies to it, which is the drift the run record has to make visible."""
    if evidenced_negative is None:
        return 0
    return int((np.isfinite(matrix) & (np.asarray(outcomes) <= 0) & np.asarray(evidenced_negative, dtype=bool)).sum())


def genes_with_several_alleles(perturbation_ids) -> set[str]:
    """The genes that appear as more than one perturbation: an allele is "<gene>|<allele>" (decision 9.54)."""
    alleles_of_gene: dict[str, int] = {}
    for perturbation_id in perturbation_ids:
        gene = gene_of_perturbation(perturbation_id)
        alleles_of_gene[gene] = alleles_of_gene.get(gene, 0) + 1
    return {gene for gene, count in alleles_of_gene.items() if count > 1 and any(identifier.startswith(gene + "|") for identifier in perturbation_ids)}


def build_soft_priors(data, arguments) -> tuple[np.ndarray | None, np.ndarray | None, dict | None, pd.DataFrame | None, pd.DataFrame | None]:
    """The soft-prior matrices for data's perturbations and symptoms (evidence/soft_priors.py): None while the strength
    is 0. Also returns the literature report table found among --soft-prior-reports, for the popularity bins, and under
    --soft-prior-probability calibrated the prior cells with their features (evidence/calibrated_reliability.py), which
    each fit calibrates on its own rows."""
    if not arguments.soft_prior_strength:
        return None, None, None, None, None
    tiers = read_tiers(arguments.soft_prior_tiers)
    records = pd.read_parquet(Path(arguments.evidence_dir) / "evidence_records.parquet")
    frames = [soft_evidence_from_records(records, tuple(arguments.soft_prior_record_grades), TRAINED_RELATION)]
    literature = None
    admission_counts: dict = {}
    allele_genes = genes_with_several_alleles(data.perturbation_ids)
    set_aside_for_alleles = 0
    for path in arguments.soft_prior_reports:
        reports = pd.read_parquet(path)
        if allele_genes and "evidence_class" in reports.columns:
            # a literature statement names a gene, not one of its alleles: where a gene has several (decision 9.54) it
            # cannot be given to one of them, so it is set aside rather than credited to the gene's default allele
            ambiguous = (reports.evidence_class == "literature") & reports.perturbation_id.astype(str).isin(allele_genes)
            set_aside_for_alleles += int(ambiguous.sum())
            reports = reports[~ambiguous]
        frames.append(soft_evidence_from_reports(reports, TRAINED_RELATION, admission_counts))
        if "evidence_class" in reports.columns and (reports.evidence_class == "literature").any():
            literature = reports if literature is None else pd.concat([literature, reports], ignore_index=True)
    resolved = resolve_pairs(pd.concat(frames, ignore_index=True), tiers, arguments.soft_prior_tempering, arguments.soft_prior_count_cap)
    probability, strength = soft_prior_matrices(resolved, data.perturbation_ids, data.symptoms)
    eligible = eligible_cells(data.outcomes, data.label_mask, data.evidenced_negative, data.masked_for_soft_evidence)
    summary = summarise_soft_priors(resolved, probability, strength, eligible)
    summary.update({"strength": arguments.soft_prior_strength, "tempering": arguments.soft_prior_tempering, "count_cap": arguments.soft_prior_count_cap, "admission": admission_counts,
                    "record_grades": list(arguments.soft_prior_record_grades), "tiers": str(arguments.soft_prior_tiers), "reports": [str(path) for path in arguments.soft_prior_reports],
                    "probability_source": arguments.soft_prior_probability})
    if allele_genes:  # recorded only on a table with alleles, so an earlier run's summary is as it was
        summary["literature_reports_set_aside_for_genes_with_several_alleles"] = set_aside_for_alleles
    calibration_cells = None
    if arguments.soft_prior_probability == "calibrated":
        popularity = {} if literature is None else dict(zip(data.perturbation_ids, np.log1p(literature_popularity(data.perturbation_ids, literature))))
        features = feature_frame(resolved.pairs, dict(zip(data.perturbation_ids, data.perturbation_types)), popularity)
        calibration_cells = cell_table(resolved.pairs, features, data.perturbation_ids, data.symptoms)
    return probability, strength, summary, literature, calibration_cells


def negative_share_scale(weights: np.ndarray, eligible: np.ndarray, negative_weight: float) -> np.ndarray:
    """Per row, how much --equal-drug-shares rescaled its unlabelled negatives: the weight of one such cell divided by
    --negative-weight (1 for a gene row, or a row with no such cell). A soft prior on a drug row is a pseudo-observation
    of that drug and takes the same share, else a drug with few labels would weigh its priors less than its labels."""
    scale = np.ones(len(weights))
    for row in range(len(weights)):
        cells = np.flatnonzero(eligible[row] & (weights[row] > 0))
        if len(cells) and negative_weight > 0:
            scale[row] = float(weights[row, cells[0]]) / negative_weight
    return scale


def fit_targets_with_soft_priors(outcomes: torch.Tensor, weights: torch.Tensor, positive_targets: torch.Tensor | None, fit_rows: np.ndarray,
                                 soft_prior_probability: np.ndarray | None, soft_prior_strength: np.ndarray | None, eligible: np.ndarray | None,
                                 global_strength: float, row_scale: np.ndarray | None = None) -> tuple[dict, dict | None]:
    """The target and weight tensors a fit trains on: the labels, and on the rows it fits the soft priors written into
    the eligible cells (evidence/soft_priors.py). Without priors the tensors are returned as they are, so the
    validation loss and every score, which read the plain matrices, see what the fit saw. With priors the targets
    tensor in use (positive_targets under --positive-target-from-frequency, else outcomes) and the weights are new
    tensors; the inputs are not modified. row_scale multiplies each row's prior strengths (negative_share_scale)."""
    targets = {"outcomes": outcomes, "weights": weights, "positive_targets": positive_targets}
    if soft_prior_probability is None:
        return targets, None
    base_targets = positive_targets.numpy() if positive_targets is not None else outcomes.numpy()
    strength = soft_prior_strength if row_scale is None else soft_prior_strength * np.asarray(row_scale, dtype=float)[:, None]
    fitted_targets, fitted_weights, applied = apply_soft_priors(base_targets, weights.numpy(), soft_prior_probability, strength, eligible, fit_rows, global_strength)
    targets["weights"] = torch.as_tensor(fitted_weights, dtype=torch.float32)
    if positive_targets is not None:
        targets["positive_targets"] = torch.as_tensor(fitted_targets, dtype=torch.float32)
    else:
        targets["outcomes"] = torch.as_tensor(fitted_targets, dtype=torch.float32)
    return targets, applied


def loss_pair_weights(data, label_mask, arguments, fit_indices: np.ndarray, frequency_targets: bool = False) -> np.ndarray:
    """The loss's weight per (perturbation, symptom): a positive at its evidence weight (at 1 when the frequency is the
    target), a negative at --negative-weight, a pair set aside by the selection at zero; under --equal-drug-shares each
    drug's positives and negatives are rescaled to the mean drug counts over fit_indices (equal_drug_loss_shares).

    A negative whose absence rests on a controlled comparison rather than on silence (the never-above-placebo pairs of a
    selection built with --never-above-placebo-as-negative, marked in data.evidenced_negative) takes
    --placebo-negative-weight, which defaults to --negative-weight so that nothing changes unless it is set. The user,
    11 October 2026, asked for a better-evidenced absence to count more in the loss share.

    Amendment item 3 makes that weight per pair where the drug's label states the two arm sizes: the pair takes
    n_effective / (n_effective + 214) rather than the default, read from --placebo-negative-weights-file into
    data.placebo_negative_weight_matrix (NaN where the file has no row for the pair, which then keeps the default).
    """
    positive = data.outcomes > 0
    weights = np.where(positive, 1.0 if frequency_targets else np.maximum(data.weights, 1e-3), arguments.negative_weight)
    evidenced_negative = getattr(data, "evidenced_negative", None)
    if evidenced_negative is not None and getattr(arguments, "placebo_negative_weight", None) is not None:
        is_evidenced_negative = ~positive & np.asarray(evidenced_negative, dtype=bool)
        weights = np.where(is_evidenced_negative, arguments.placebo_negative_weight, weights)
        per_pair = getattr(data, "placebo_negative_weight_matrix", None)
        if per_pair is not None:
            per_pair = np.asarray(per_pair, dtype=float)
            weights = np.where(is_evidenced_negative & np.isfinite(per_pair), per_pair, weights)
    if label_mask is not None:
        weights = weights * label_mask
    if getattr(arguments, "equal_drug_shares", False):
        weights = equal_drug_loss_shares(weights, positive, np.array([kind in DRUG_ROW_KINDS for kind in data.perturbation_types]), fit_indices,
                                         drug_loss_slots(data))
    return weights


def clamped_probability_report(predictions: np.ndarray, data, indices: np.ndarray, label_mask) -> dict:
    """How many test probabilities sit beyond the loss's clamp, and how many labelled pairs are among them.

    Registered as item 19 of the amendment of 12 October 2026 (the owner: "Leave the clamp off for now but make sure
    you check for the issue in the dev runs"). The training loss takes the logarithm of P clamped to
    [PROBABILITY_EPSILON, 1 - PROBABILITY_EPSILON], so a pair beyond the clamp has a constant loss and no gradient: a
    positive below it cannot be pulled back. The evidence that raised this was measured on the sigmoid head, which the
    confirmatory family no longer holds, so for the noisy-OR head the clamp is untested rather than harmless. A
    labelled positive below the clamp in any development run is the finding that reopens --bce-in-log-space.
    """
    probabilities = np.asarray(predictions, dtype=float)
    outcomes = np.asarray(data.outcomes)[indices] > 0
    labelled = np.ones_like(outcomes, dtype=bool) if label_mask is None else np.asarray(label_mask)[indices].astype(bool)
    below, above = probabilities < PROBABILITY_EPSILON, probabilities > 1.0 - PROBABILITY_EPSILON
    return {"epsilon": PROBABILITY_EPSILON, "pairs": int(probabilities.size),
            "below_the_clamp": int(below.sum()), "above_the_clamp": int(above.sum()),
            "share_below_the_clamp": round(float(below.mean()), 6),
            "labelled_positives_below_the_clamp": int((below & outcomes & labelled).sum()),
            "labelled_negatives_above_the_clamp": int((above & ~outcomes & labelled).sum()),
            "reopens_bce_in_log_space": bool((below & outcomes & labelled).any())}


def message_passing_edge_arrays(data, edge_signs: str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None, torch.Tensor | None]:
    """(edge index, relation, sign or None, weight or None): what the message-passing adjacencies are built from, and
    what its edge gates are laid out on (attach_edge_gates), so the two cannot read different signs or weights."""
    if edge_signs not in MESSAGE_PASSING_EDGE_SIGNS:
        raise ValueError(f"--message-passing-edge-signs must be one of {MESSAGE_PASSING_EDGE_SIGNS}, not {edge_signs!r}")
    edge_sign = None if edge_signs == "unsigned" or data.edge_sign is None else torch.as_tensor(np.asarray(data.edge_sign), dtype=torch.float32)
    edge_weight = None if edge_signs == "unsigned" or data.edge_weight is None else torch.as_tensor(np.asarray(data.edge_weight), dtype=torch.float32)
    return (torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long), torch.as_tensor(data.edge_relation, dtype=torch.long),
            edge_sign, edge_weight)


def message_passing_adjacencies(data, edge_signs: str) -> list:
    """The message-passing encoder's per-relation adjacencies, signed from the graph under edge_signs "graph".

    Registered as item 9 of the amendment of 12 October 2026 (the owner: "The edge signs are also settled: the MP has
    to adhere to the edge signs the same way the linear encoder does"), so "graph" is the default and "unsigned"
    reproduces every run made before that date. Under "graph" an edge also carries the weight the graph gives it, which
    today is the driving-force edges of the computed membrane potential; without such a column every weight is 1 and
    only the signs change."""
    edge_index, edge_relation, edge_sign, edge_weight = message_passing_edge_arrays(data, edge_signs)
    return RelationalMessagePassingEncoder.build_relation_adjacencies(edge_index, edge_relation, len(data.node_ids), len(data.relation_types),
                                                                      edge_sign=edge_sign, edge_weight=edge_weight)


def build_edge_priors(data, arguments) -> dict | None:
    """The edge priors of decision 9.57, set on data aligned to its loaded edges (edge_prior_gated, edge_prior_probability,
    edge_prior_strength) with their summary; at --edge-prior-strength 0, the registered model, the three are None and so
    is the summary (set to None in either case, so a data set reused across runs never keeps another run's priors).

    The table is <graph-dir>/edge_evidence.parquet (experiments/build_edge_evidence.py), one row per row of the graph's
    edges.parquet with each edge's source_id, target_id and relation_type. It is refused unless those are the graph's,
    and read through the loader's edge rows (ExperimentData.edge_row), which skip what drug_entry "targets" dropped. Read
    before any rewiring, which then keeps each prior at its edge position as it keeps the signs (main).

    With --curated-edge-labels (decision 9.63 (b), the default with edge priors) a curated edge with enough independent
    records, enough of them human evidence, and no IntAct absence is taken as known: ungated, so always open, as the
    simulation's --ungated-tiers did. The permutation control then shuffles the priors of the edges still gated; the
    labels stay on their edges. data.edge_evidence keeps the table's rows for the loaded edges, which the review list of
    the closures reads (9.63 (c))."""
    data.edge_prior_gated = data.edge_prior_probability = data.edge_prior_strength = None
    data.edge_evidence = None
    strength = float(getattr(arguments, "edge_prior_strength", 0.0) or 0.0)
    if strength == 0.0:
        return None
    if strength < 0.0:
        raise ValueError("--edge-prior-strength multiplies pseudo-counts and cannot be negative")
    if arguments.encoder not in EDGE_GATED_ENCODERS:
        raise ValueError(f"--edge-prior-strength gates the edges an encoder propagates along; the {arguments.encoder} encoder reads none (use one of {EDGE_GATED_ENCODERS})")
    path = Path(arguments.graph_dir) / EDGE_EVIDENCE_FILE
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing; build it with: python experiments/build_edge_evidence.py --graph-dir {arguments.graph_dir}")
    if getattr(data, "edge_row", None) is None:
        raise ValueError("these data carry no edge rows (ExperimentData.edge_row); load them with load_experiment_data")
    endpoint_columns = ["source_id", "target_id", "relation_type"]
    evidence = pd.read_parquet(path)
    graph_edges = pd.read_parquet(Path(arguments.graph_dir) / "edges.parquet", columns=endpoint_columns)
    missing = [column for column in endpoint_columns if column not in evidence.columns]
    if missing or len(evidence) != len(graph_edges) or not all(np.array_equal(evidence[column].to_numpy(dtype=object), graph_edges[column].to_numpy(dtype=object))
                                                             for column in endpoint_columns):
        raise ValueError(f"{path} does not hold this graph's edges row for row" + (" (it was written before the endpoint columns)" if missing else "")
                         + f"; rebuild it with: python experiments/build_edge_evidence.py --graph-dir {arguments.graph_dir}")
    rows = np.asarray(data.edge_row, dtype=np.int64)
    # the loaded edges are those rows: a rewiring moves targets only and runs after this, so sources and relations must agree
    node_ids, relation_names = np.asarray(data.node_ids, dtype=object), np.asarray(data.relation_types, dtype=object)
    if not (np.array_equal(node_ids[np.asarray(data.edge_source)], evidence.source_id.to_numpy(dtype=object)[rows])
            and np.array_equal(relation_names[np.asarray(data.edge_relation)], evidence.relation_type.to_numpy(dtype=object)[rows])):
        raise ValueError(f"the loaded edges are not the rows of {path} the loader names; the edge table and the graph disagree")
    evidence = evidence.iloc[rows].reset_index(drop=True)
    gated, probability, pseudo_count = edge_priors(evidence, read_edge_tiers(arguments.edge_prior_tiers), global_strength=strength,
                                                   tempering=arguments.edge_prior_tempering, count_cap=arguments.edge_prior_count_cap)
    labels_summary = None
    if getattr(arguments, "curated_edge_labels", False):
        minimum_records = getattr(arguments, "curated_label_minimum_records", DEFAULT_LABEL_MINIMUM_RECORDS)
        minimum_human_records = getattr(arguments, "curated_label_minimum_human_records", DEFAULT_LABEL_MINIMUM_HUMAN_RECORDS)
        minimum_sign_records = getattr(arguments, "curated_label_minimum_sign_records", None)  # absent from a run recorded before 17 October 2026
        exclude_disputed_signs = bool(getattr(arguments, "curated_label_exclude_disputed_signs", False))  # likewise
        bulk_papers = getattr(arguments, "curated_label_bulk_papers", "count")  # likewise
        try:
            labels = curated_edge_labels(evidence, minimum_records=minimum_records, minimum_human_records=minimum_human_records,
                                         minimum_sign_records=minimum_sign_records, exclude_disputed_signs=exclude_disputed_signs,
                                         bulk_papers=bulk_papers)
        except ValueError as error:
            raise ValueError(f"{path}: {error} (or run with --no-curated-edge-labels)") from None
        labels_summary = {"edges": int(labels.sum()), "were_gated": int((labels & gated).sum()), "minimum_records": minimum_records,
                          "minimum_human_records": minimum_human_records, "human_evidence": "NIH iCite human score above 0 (experiments/fetch_icite_species.py)",
                          "minimum_sign_records": minimum_sign_records, "exclude_disputed_signs": exclude_disputed_signs, "bulk_papers": bulk_papers}
        gated, probability, pseudo_count = gated & ~labels, np.where(labels, np.nan, probability), np.where(labels, 0.0, pseudo_count)
    summary = {"path": str(path), "sha256": file_sha256(path), "tiers": str(arguments.edge_prior_tiers), "tiers_sha256": file_sha256(arguments.edge_prior_tiers),
               "strength": strength, "tempering": arguments.edge_prior_tempering, "count_cap": arguments.edge_prior_count_cap,
               **summarise_edge_priors(evidence, gated, probability, pseudo_count), "curated_edge_labels": labels_summary, "permutation": None}
    if getattr(arguments, "permute_edge_priors", False):
        probability, pseudo_count = permute_edge_priors(gated, probability, pseudo_count, np.asarray(data.edge_relation), random_seed=arguments.seed)
        summary["permutation"] = {"seed": arguments.seed, "within": "relation", "priors_sha256": array_sha256(np.nan_to_num(probability, nan=-1.0))}
    data.edge_prior_gated, data.edge_prior_probability, data.edge_prior_strength = gated, probability, pseudo_count
    data.edge_evidence = evidence
    return summary


def edge_gates_of(data, arguments=None) -> EdgeGates | HardConcreteEdgeGates | None:
    """Fresh edge gates at the data's edge priors, or None without them (the registered model). --edge-gate-kind picks
    the hard gates decided on 16 October 2026 (decision 9.62 (e), the default) or the soft gates of 9.57."""
    if getattr(data, "edge_prior_gated", None) is None:
        return None
    kind = getattr(arguments, "edge_gate_kind", "soft") if arguments is not None else "soft"
    if kind == "hard_concrete":
        return HardConcreteEdgeGates(data.edge_prior_gated, data.edge_prior_probability, data.edge_prior_strength,
                                     beta=getattr(arguments, "edge_gate_beta", HARD_CONCRETE_BETA))
    return EdgeGates(data.edge_prior_gated, data.edge_prior_probability, data.edge_prior_strength)


def attach_edge_gates(encoder, data, adjacencies, arguments) -> None:
    """The message-passing encoder's edge gates, laid out on the adjacencies it reads, with the signs and weights those
    were built from. A no-op for any other encoder (the linear-response one takes its gates in build_models) and without
    edge priors. Call it on every freshly built encoder before its optimizer and before a state is loaded into it."""
    if not isinstance(encoder, RelationalMessagePassingEncoder):
        return
    gates = edge_gates_of(data, arguments)
    if gates is None:
        return
    edge_index, edge_relation, edge_sign, edge_weight = message_passing_edge_arrays(data, getattr(arguments, "message_passing_edge_signs", "graph"))
    encoder.attach_edge_gates(gates, edge_index, edge_relation, adjacencies, edge_sign=edge_sign, edge_weight=edge_weight)


def drug_loss_slots(data) -> np.ndarray | None:
    """One slot id per perturbation for the equal drug shares, or None for a table with no combination in it.

    A combination of a drug and a lowered gene carries that drug's labels, so on its own slot it would give the drug
    two shares of the loss where every drug without a statement has one (models/soft_constraint_losses.py). The first
    drug part is the drug the statement is about, which for every kind built so far is the only drug part. None leaves
    every row its own drug, which is the rule the shares were registered with."""
    if getattr(data, "perturbation_parts", None) is None:
        return None
    slots = list(data.perturbation_ids)
    for position, kind in enumerate(data.perturbation_types):
        if kind != COMBINATION_PERTURBATION_TYPE:
            continue
        drug_parts = [part.part_id for part in data.perturbation_parts[position] if part.part_kind == "drug"]
        if drug_parts:
            slots[position] = drug_parts[0]
    return np.array(slots, dtype=object)


def permuted_drug_targets(perturbation_types: list[str], perturbation_seeds: list[np.ndarray], perturbation_signs: list[np.ndarray],
                          perturbation_magnitudes: list[np.ndarray], random_seed: int) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray], int]:
    """Each drug's seeds exchanged with those of another drug that has the same number of them; genes untouched.

    The third of the three controls for the input-size objection (docs/label_source_review.md, 13.1). Every drug keeps
    its number of targets and its total input magnitude and loses which targets they are, so a model that scores as
    well here as on the real targets is reading the size of a drug's input and not what its targets do. The exchange is
    inside a group of equal count, which is what makes it a control of identity alone: a drug with three targets gets
    three other targets, never one or five.

    A count whose group holds one drug cannot be exchanged and is left as it is; the returned number of drugs moved
    says how many were. The permutation is a derangement where the group allows one (a cycle of the group, so no drug
    keeps its own targets), seeded by random_seed.

    A combination is left as it is as well, and for a reason rather than an omission: its seeds are its drug part's
    targets together with its lowered gene, and exchanging the whole list would move the lowered gene too, which is not
    an input-size confound and is the part of the example the statement is about. The control is a control of a drug's
    target identity, so it is applied where a row is a drug.

    Returns (seeds, signs, magnitudes, drugs moved); the three lists are new and the caller's are not modified."""
    generator = np.random.default_rng(random_seed)
    seeds, signs, magnitudes = list(perturbation_seeds), list(perturbation_signs), list(perturbation_magnitudes)
    positions_by_count: dict[int, list[int]] = {}
    for position, kind in enumerate(perturbation_types):
        if kind == "drug":
            positions_by_count.setdefault(len(perturbation_seeds[position]), []).append(position)
    moved = 0
    for count in sorted(positions_by_count):
        group = positions_by_count[count]
        if len(group) < 2:
            continue
        order = list(generator.permutation(len(group)))  # a cycle of a shuffled order: every member takes another's seeds
        taken_from = [group[order[(index + 1) % len(order)]] for index in range(len(order))]
        for position, source in zip([group[index] for index in order], taken_from):
            seeds[position], signs[position], magnitudes[position] = perturbation_seeds[source], perturbation_signs[source], perturbation_magnitudes[source]
            moved += 1
    return seeds, signs, magnitudes, moved


def drug_input_normalised(perturbation_types: list[str], perturbation_magnitudes: list[np.ndarray],
                          perturbation_parts: list[list] | None = None, perturbation_seed_parts: list[np.ndarray] | None = None) -> list[np.ndarray]:
    """Seed magnitudes with every drug's rescaled to sum to 1 in absolute value; genes are unchanged. A drug's magnitude is
    otherwise 1 per target (map_drug_targets_to_graph_nodes.py), so a drug with k targets pushes k times the field of a
    drug with one, and a model could learn the count of targets instead of what the targets do.

    In a combination the rescaling is per part: each drug part's seeds are rescaled to sum to 1 among themselves and
    each gene part is left alone, which is what makes the combination's drug input the same size as that drug's input
    on its own. Without the parts a combination is left unchanged rather than rescaled as a whole, because rescaling
    the whole list would make the lowered gene's magnitude depend on how many targets the drug has."""
    normalised = []
    for position, (kind, magnitudes) in enumerate(zip(perturbation_types, perturbation_magnitudes)):
        magnitudes = np.asarray(magnitudes, dtype=float)
        if kind == COMBINATION_PERTURBATION_TYPE:
            if perturbation_parts is None or perturbation_seed_parts is None:
                normalised.append(magnitudes)
                continue
            rescaled = magnitudes.copy()
            seed_part = np.asarray(perturbation_seed_parts[position], dtype=int)
            for part_position, part in enumerate(perturbation_parts[position]):
                of_this_part = seed_part == part_position
                total = np.abs(magnitudes[of_this_part]).sum()
                if part.part_kind == "drug" and total > 0:
                    rescaled[of_this_part] = magnitudes[of_this_part] / total
            normalised.append(rescaled)
            continue
        total = np.abs(magnitudes).sum()
        normalised.append(magnitudes / total if kind == "drug" and total > 0 else magnitudes)
    return normalised


def weighted_constant_optimum(data, fit_indices: np.ndarray, label_mask, arguments) -> np.ndarray:
    """Per symptom, the constant probability that minimises the training loss over fit_indices: the weighted mean of the
    loss's targets, sum(w t) / sum(w), with the loss's weights (a positive at its evidence weight, a negative at
    --negative-weight, a pair set aside by the selection at zero) and targets (--positive-target, or the frequency under
    --positive-target-from-frequency). With negatives at 0.2 it is about three times the raw base rate (0.204 against
    0.064 on the full-graph development pool, docs/best_epoch_zero.md)."""
    outcomes = data.outcomes[fit_indices]
    positive = outcomes > 0
    if arguments.positive_target_from_frequency:
        frequencies = data.frequencies[fit_indices]
        targets = np.where(positive, np.where(np.isnan(frequencies), arguments.positive_target, np.maximum(frequencies, arguments.minimum_frequency_target)), 0.0)
    else:
        targets = np.where(positive, arguments.positive_target, 0.0)
    weights = loss_pair_weights(data, label_mask, arguments, fit_indices, frequency_targets=arguments.positive_target_from_frequency)[fit_indices]
    return (weights * targets).sum(axis=0) / np.maximum(weights.sum(axis=0), 1e-12)


def initialise_rates_at_weighted_optimum(head, data, fit_indices: np.ndarray, label_mask, arguments, device) -> None:
    """--start-at-weighted-optimum: start each symptom's noisy-OR leak (in place of the raw base rate of
    --init-leak-from-base-rate) or the sigmoid head's output bias at the loss-optimal constant of weighted_constant_optimum,
    so that the first epochs are not spent on calibration. Measured on the full-graph pilots (docs/best_epoch_zero.md): a
    leak started at the raw base rate stays there at --leak-learning-rate 0.0002, and the modules make up the gap to the
    optimum; a sigmoid head starts at P = 0.5 and reaches the optimum within the first epoch, which is then the best one."""
    if not getattr(arguments, "start_at_weighted_optimum", False):
        return
    optimum = torch.as_tensor(weighted_constant_optimum(data, fit_indices, label_mask, arguments), dtype=torch.float32, device=device)
    if arguments.head == "noisy_or":
        head.initialize_leak_from_base_rates(optimum)
        return
    clamped = optimum.clamp(PROBABILITY_EPSILON, 1.0 - PROBABILITY_EPSILON)
    with torch.no_grad():
        head.readout[-1].bias.copy_(torch.log(clamped) - torch.log1p(-clamped))


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


def never_scored_rows(data) -> np.ndarray:
    """Rows of a perturbation type that is never scored (PERTURBATION_TYPES_NEVER_SCORED, item 12c). Data without a
    perturbation type at all has no such row: the splits that hold out a module or a subsystem are built from seeds
    and group ids and are given data that carries neither types nor outcomes."""
    types = getattr(data, "perturbation_types", None)
    if types is None:
        return np.zeros(len(data.perturbation_ids), dtype=bool)
    return np.array([kind in PERTURBATION_TYPES_NEVER_SCORED for kind in types], dtype=bool)


def withhold_never_scored_rows(data, test_mask: np.ndarray, train_pool: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Take the never-scored rows out of the test set, and out of training where their leakage group is held out.

    "Training only" cannot mean "always trained", which would be the leak the leakage groups exist to stop: training
    the receptor-activation row of a gene whose knockout sits in the held-out set shows the model that node perturbed.
    So such a row is **dropped** - neither scored nor trained - whenever its group is held out, and trained
    otherwise. The held-out draw over every other perturbation is untouched either way, because these rows leave the
    test set whatever group they are in."""
    never_scored = never_scored_rows(data)
    if not never_scored.any():
        return test_mask, train_pool, 0, 0
    group_ids = np.asarray(data.group_ids)
    held_out_groups = set(group_ids[test_mask])
    taken_out_of_test = int((test_mask & never_scored).sum())
    test_mask = test_mask & ~never_scored
    in_a_held_out_group = never_scored & np.array([group in held_out_groups for group in group_ids])
    kept_in_training = np.ones(len(group_ids), dtype=bool)
    kept_in_training[in_a_held_out_group] = False
    dropped = int(in_a_held_out_group.sum())
    return test_mask, train_pool[kept_in_training[train_pool]], taken_out_of_test, dropped


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
    test_mask, train_pool, taken_out_of_test, dropped = withhold_never_scored_rows(data, test_mask, train_pool)
    if taken_out_of_test or dropped:
        print(f"never-scored perturbation types: {taken_out_of_test} rows taken out of the test set, "
              f"{dropped} dropped from training because their leakage group is held out")
    validation = np.array([], dtype=int)
    if arguments.validation_fraction > 0 and len(train_pool) >= 20:
        train_pool, validation = early_stopping_validation(data, train_pool, arguments)
        # The validation set picks the stopping epoch by macro AUPRC, which is a score, so a never-scored row must
        # not reach it either: item 12c keeps these rows out of every reading precisely because macro AUPRC is taken
        # per symptom and 43 percent of their labels are one symptom. They go back to training instead of being lost.
        never_scored = never_scored_rows(data)
        if never_scored[validation].any():
            moved = validation[never_scored[validation]]
            validation = validation[~never_scored[validation]]
            train_pool = np.sort(np.concatenate([train_pool, moved]))
            print(f"never-scored perturbation types: {len(moved)} rows moved out of the validation set into training")
    return train_pool, validation, all_indices[test_mask], split_name


def time_training_steps(data, encoder, head, optimizer, train_indices, validation_indices, adjacencies, label_mask, arguments, device) -> None:
    """Seconds per training step (forward, backward, optimizer step) and per prediction batch; nothing is scored."""
    outcomes = torch.as_tensor(data.outcomes, dtype=torch.float32)
    weights = torch.as_tensor(loss_pair_weights(data, label_mask, arguments, train_indices), dtype=torch.float32)
    order = np.random.default_rng(arguments.seed).permutation(train_indices)
    step_seconds = []
    gates = getattr(encoder, "edge_gates", None)
    for step in range(arguments.timing_batches):
        batch = order[(step * arguments.batch_size) % len(order):][: arguments.batch_size]
        started = time.time()
        if isinstance(gates, HardConcreteEdgeGates) and gates.phase == "stochastic":
            gates.resample()  # as training_step does: a sample is drawn per batch, and a kept one would carry the last step's freed graph
        node_index, sign_and_magnitude = pad_perturbations(data, batch)
        field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
        output = head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
        loss = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device), positive_target=arguments.positive_target)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        step_seconds.append(time.time() - started)
    prediction_times = []
    with torch.no_grad():
        encoder.eval()
        head.eval()
        for repeat in range(2):  # the first reads after the train-mode steps (caches rebuilt for the read-out); the second is the steady cost
            batch = order[repeat * arguments.batch_size:][: arguments.batch_size]
            started = time.time()
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            head(encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field), perturbation_covariate=covariate_of(data, batch, arguments, device))
            prediction_times.append(time.time() - started)
    prediction_seconds = prediction_times[-1]
    steady = float(np.median(step_seconds[1:] if len(step_seconds) > 1 else step_seconds))
    steps_per_epoch = int(np.ceil(len(train_indices) / arguments.batch_size))
    prediction_batches = int(np.ceil(len(validation_indices) / arguments.batch_size))
    print(json.dumps({"timing_batches": arguments.timing_batches, "first_step_seconds": round(step_seconds[0], 2), "median_step_seconds": round(steady, 2),
                      "prediction_batch_seconds": round(prediction_seconds, 2), "first_prediction_batch_seconds": round(prediction_times[0], 2),
                      "train_perturbations": len(train_indices), "steps_per_epoch": steps_per_epoch,
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
    parser.add_argument("--recompute-reaction-expression-on-rewired-graph", action="store_true",
                        help="with --rewire-swaps-per-edge: rewrite each reaction's gene rule onto its rewired catalysts and recompute the reaction_brain "
                             "descriptor block and the reactions' cell-class weights from it (recompute_reaction_expression); off by default, which keeps "
                             "the real graph's reaction expression as every rewired run so far did")
    parser.add_argument("--gene-expression-table", type=Path, default=Path("data/processed/brain_expression/gene_expression_for_descriptors.parquet"),
                        help="the Ensembl-indexed table the brain descriptors were built from (experiments/write_expression_tables.py)")
    parser.add_argument("--class-expression-table", type=Path, default=Path("data/processed/brain_expression/class_expression_for_weights.parquet"),
                        help="the Ensembl-indexed class nCPM table the cell-class weights were built from (experiments/write_expression_tables.py)")
    parser.add_argument("--rewire-encodes", action="store_true",
                        help="with --rewire-swaps-per-edge on a split graph: rewire the encodes relation too (gene -> protein); by default it is held "
                             "fixed, so a knockout keeps its own protein (docs/gene_protein_split.md)")
    parser.add_argument("--keep-reciprocated-relations-symmetric", action="store_true",
                        help="with --rewire-swaps-per-edge: a relation stored in both directions (binds) is rewired as undirected edges, so it stays symmetric "
                             "(negative_controls.reciprocated_relations, fast_degree_preserving_rewiring); off by default, which reproduces earlier rewired runs")
    parser.add_argument("--time-split-cutoff", type=date.fromisoformat, default=None,
                        help="monogenic time split (design 6.1): train on pairs dated on or before this day across all perturbations; score the pairs that could still become positive")
    parser.add_argument("--group-by", choices=["gene", "disease_cluster", "disease_cluster_and_targets"], default="disease_cluster_and_targets")
    parser.add_argument("--label-grades", nargs="*", default=["A", "B"], help="evidence grades that count as positive labels; pass A B C to keep grade C rows as the version 0.3 ablation did")
    parser.add_argument("--message-passing-edge-signs", choices=MESSAGE_PASSING_EDGE_SIGNS, default="graph",
                        help="whether the message-passing encoder's adjacency carries the graph's edge signs and weights. \"graph\" is the registered "
                             "behaviour from 12 October 2026 (amendment item 9); \"unsigned\" reproduces every run made before it, when an inhibition "
                             "and an activation of the same relation arrived alike")
    parser.add_argument("--combined-perturbations", type=Path, default=None,
                        help="parquet of combination rows from experiments/build_combined_perturbations.py: a drug carried with a lowered gene, or with a second drug, "
                             "in one example (mechanistic_pathway_learning/evidence/combined_perturbations.py). Without it the table holds no combination and the run is unchanged")
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="parquet of (perturbation_id, symptom, keep) from experiments/build_label_selection.py; positive pairs with keep False are masked out of the loss and every metric")
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--refit-on-validation", action=argparse.BooleanOptionalAction, default=None,
                        help="after early stopping, train a fresh model on the training and validation perturbations together for best epoch + 1 epochs "
                             "and score the test with it (the early-stopped model's test predictions are kept in test_predictions_early_stopped.npy). "
                             "On by default (registered); off by default only where it is not implemented, the laboratory readout and the time split")
    parser.add_argument("--keep-large-groups-in-training", action=argparse.BooleanOptionalAction, default=True,
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
    parser.add_argument("--node-features", choices=["identity", "typed"], default="typed",
                        help="identity: a learned embedding per node; typed: fixed structural features only (type, compartment, degree, flags), the inductive variant")
    parser.add_argument("--pooling", choices=["sum", "mean"], default="sum")
    parser.add_argument("--conjunction-relations", nargs="*", default=list(DEFAULT_CONJUNCTION_RELATIONS),
                        help="relations whose messages combine by a per-channel minimum rather than a mean, because their "
                             "target is present only to the extent that every source is (message-passing encoder only)")
    parser.add_argument("--conjunction-aggregation", choices=list(CONJUNCTION_AGGREGATIONS), default="minimum",
                        help="minimum: the hard per-channel minimum; soft_minimum: its temperature surrogate, the control "
                             "that separates a conjunction from a hard minimum; mean: the in-degree mean, the arm before "
                             "this option existed")
    parser.add_argument("--soft-minimum-temperature", type=float, default=DEFAULT_SOFT_MINIMUM_TEMPERATURE,
                        help="temperature of --conjunction-aggregation soft_minimum; the hard minimum is its limit at 0")
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
    parser.add_argument("--descriptor-treatment", choices=list(DESCRIPTOR_TREATMENTS), default=None,
                        help="with --node-descriptors, how they enter the message passing base state or the linear-response output gate (models/descriptor_treatments.py): "
                             "plain (one map of all the columns), seed_masked (each perturbation's perturbed nodes read their structural columns only) or "
                             "zero_init_slow (a separate descriptor map, initialised at zero and trained at --descriptor-learning-rate). Left out: seed_masked "
                             "with --node-descriptors (the owner, 13 October 2026), plain without them, since there is nothing to mask")
    parser.add_argument("--descriptor-learning-rate", type=float, default=0.0,
                        help="--descriptor-treatment zero_init_slow: learning rate of the descriptor map (required there)")
    parser.add_argument("--cell-class-weights", type=Path, default=None,
                        help="linear-response encoder: parquet of per-node cell-class weights (experiments/build_cell_class_weights.py); the response then propagates once per class")
    parser.add_argument("--channels-per-cell-class", type=int, default=1, help="linear-response encoder with --cell-class-weights: channels (time scales) per class")
    parser.add_argument("--extracellular-coupling", action=argparse.BooleanOptionalAction, default=None,
                        help="linear-response encoder with --cell-class-weights: extracellular metabolites are one pool shared by every class but all_cells. "
                             "Left out, on whenever the linear-response encoder has cell-class weights (registered)")
    parser.add_argument("--compartment-ph", action=argparse.BooleanOptionalAction, default=None,
                        help="linear-response encoder: the settled proton deviation of each compartment becomes a throughput ceiling on its "
                             "nodes (graph/compartment_ph.py). Left out, on for the linear-response encoder wherever the graph names a proton "
                             "node in a single compartment (registered)")
    parser.add_argument("--class-weight-rule", choices=list(CLASS_WEIGHT_RULES), default="limiting_step",
                        help="linear-response encoder with --cell-class-weights: how a class's node weights attenuate the response. limiting_step (the owner, "
                             "13 October 2026), the weakest weight on the widest route; product, the weights multiplied along the path, every run before")
    parser.add_argument("--input-weights", choices=list(INPUT_WEIGHT_FORMS), default="positive",
                        help="linear-response encoder: positive (the owner, 13 October 2026: each class consistently signed), every channel's input weight "
                             "positive and starting equal; signed, one free weight per channel drawn with a random sign, every run before")
    parser.add_argument("--degree-offset", action="store_true",
                        help="give the head the standardised log degree of each perturbation: a degree-dependent leak (noisy-OR) or logit offset (sigmoid), so the field only has to explain what degree does not")
    parser.add_argument("--propagation-steps", type=int, default=8, help="linear-response encoder: steps of the shared transition (the reach in edges); ignored under --fixed-point-tolerance")
    parser.add_argument("--fixed-point-tolerance", type=tolerance_or_none, default=DEFAULT_FIXED_POINT_TOLERANCE,
                        help="run the propagation until the relative change of the state falls below this, in place of a fixed number of steps "
                             "(amendment item 7 of 12 October 2026 registers 1e-3, the default). 'none' takes --propagation-steps rounds, which "
                             "reproduces every run before that date")
    parser.add_argument("--max-propagation-steps", type=int, default=DEFAULT_MAX_PROPAGATION_STEPS,
                        help="the cap on the rounds under --fixed-point-tolerance, so a graph that does not settle still terminates; a run that "
                             "reaches it has not settled and the run record says so")
    parser.add_argument("--propagation-channels", type=int, default=4, help="linear-response encoder: channels propagated with their own gains (time scales), expanded linearly to --node-state-dim")
    parser.add_argument("--cofactor-relations", action=argparse.BooleanOptionalAction, default=None,
                        help="linear-response encoder: give carrier edges (cofactor_edges.py) relations of their own, so their coupling gets learned gains. "
                             "Left out, on for the linear-response encoder (registered)")
    parser.add_argument("--reversible-reactions", action=argparse.BooleanOptionalAction, default=None,
                        help="linear-response encoder: derive product -> reaction at sign -1 for every reversible reaction (product_opposes, "
                             "preregistration item 27), so a perturbation that accumulates a product drives the reaction backwards. Needs a graph "
                             "built without the old mirror edges and refuses one that has them. Left out, on for the linear-response encoder (registered)")
    parser.add_argument("--carrier-rule", choices=CARRIER_RULES, default="all",
                        help="which clause of the carrier rule applies (cofactor_edges.py): all is the curated list and the recurring-pair "
                             "heuristic, names drops the heuristic, neuronal keeps only the transmitter-synthesis and oxidative carriers")
    parser.add_argument("--response-scale", choices=["linear", "signed_log"], default="linear",
                        help="linear-response encoder: read the response as it is, or through sign(h) log(1 + |h| / s) with a learned scale, so changes many edges away stay readable")
    parser.add_argument("--no-stoichiometric-demand", dest="stoichiometric_demand", action="store_false",
                        help="leave every bottleneck demand at 1, the behaviour before item 20 was built")
    parser.add_argument("--normalisation", choices=NORMALISATIONS, default="total_in_degree",
                        help="linear-response encoder: how the signed adjacency is scaled so the iteration contracts. total_in_degree (registered, amendment item 28) "
                             "divides by the destination's in-degree summed over every relation, so the unsigned row sums are one and the stoichiometric counts "
                             "survive; in_degree divides per relation as well, which makes every metabolite a mean against a mean and zeroes 92 percent of them; "
                             "spectral applies one global scale, which also keeps the counts but shrinks the field nearly everywhere")
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
    parser.add_argument("--num-layers", type=int, default=3, help="message passing: layers; 3 on a graph that splits genes from proteins (registered)")
    parser.add_argument("--relation-parameterisation", choices=list(RELATION_PARAMETERISATIONS), default="shared_low_rank",
                        help="message passing: how a relation's transform is parameterised (amendment item 8). shared_low_rank, the registered "
                             "form, applies one set at every round, each relation owning a diagonal and the relations enough training examples "
                             "reach a rank-r factorisation, with a norm bound; per_layer_dense is a full matrix per relation per layer, which "
                             "every run before 13 October 2026 built")
    parser.add_argument("--relation-factor-rank", type=int, default=DEFAULT_RELATION_FACTOR_RANK,
                        help="message passing with --relation-parameterisation shared_low_rank: the factorisation's rank. Amendment item 8 "
                             "registers 2 and names 4 as the only step up that may be taken")
    parser.add_argument("--factorised-relations", nargs="*", default=list(REGISTERED_FACTORISED_RELATIONS),
                        help="message passing with --relation-parameterisation shared_low_rank: the relations that own a factorisation beside "
                             "their diagonal; the registered nine are the default, and a relation the graph does not have is left out")
    parser.add_argument("--num-modules", type=int, default=8)
    parser.add_argument("--sigmoid-hidden-dim", type=int, default=64)
    parser.add_argument("--gate-initial-log-alpha", type=float, default=-1.0)
    parser.add_argument("--init-leak-from-base-rate", action="store_true", help="noisy-OR head: start each symptom's leak at its training base rate")
    parser.add_argument("--start-at-weighted-optimum", action=argparse.BooleanOptionalAction, default=True,
                        help="start the noisy-OR leaks (in place of --init-leak-from-base-rate) or the sigmoid head's output biases at each symptom's "
                             "loss-optimal constant over the fitted perturbations, the weighted mean of the loss's targets (docs/best_epoch_zero.md)")
    parser.add_argument("--module-bias-init", type=float, default=-3.0, help="noisy-OR head: initial readout bias of every module; negative values make modules off by default")
    parser.add_argument("--gate-init-noise", type=float, default=0.5, help="noisy-OR head: standard deviation of the per-gate noise added to the initial log-alpha (symmetry breaking between modules)")
    parser.add_argument("--link-learning-rate", type=float, default=0.02,
                        help="noisy-OR head: learning rate of the module-to-symptom links (0 = the main rate); links start at logit -3 and Adam moves them about one rate per step, so at 0.002 they need about 1,500 steps to reach 0.5")
    parser.add_argument("--leak-learning-rate", type=float, default=0.0002, help="noisy-OR head: learning rate of the symptom leaks (0 = the main rate); slow, so a leak started at the base rate stays there")
    parser.add_argument("--module-bias-learning-rate", type=float, default=0.0, help="noisy-OR head: learning rate of the module readout biases (0 = the main rate)")
    parser.add_argument("--gate-learning-rate", type=float, default=0.05, help="noisy-OR head: learning rate of the support gate log-alphas (0 = the main rate); fast, so gates no field reaches can close within the epochs early stopping allows")
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
    parser.add_argument("--placebo-negative-weight", type=float, default=0.5,
                        help="loss weight of a negative that rests on a controlled comparison rather than on silence: the never-above-placebo pairs of a selection "
                             "built with --never-above-placebo-as-negative. 0.5 is registered (amendment item 3 of 12 October 2026), and it is the weight of a "
                             "pair the file below does not name")
    parser.add_argument("--placebo-negative-weights-file", type=path_or_none, default=REGISTERED_PLACEBO_NEGATIVE_WEIGHTS,
                        help="per-pair weights for those negatives (amendment item 3): n_effective / (n_effective + 214) where the drug's label states both arm "
                             "sizes, the default above where it does not. 'none' falls back to --placebo-negative-weight for every pair")
    parser.add_argument("--soft-prior-strength", type=float, default=0.0,
                        help="global multiplier on the soft priors of evidence/soft_priors.py (the power-prior exponent): a grade-C record or a report of the "
                             "trial-arm, FAERS, IMPC or literature tables becomes a target and a weight on an unlabelled pair of a fitted row. 0 (the default, and "
                             "the registered model) turns them off and leaves every run as it was")
    parser.add_argument("--soft-prior-reports", type=Path, nargs="*", default=[],
                        help="report tables (evidence_reports schema) whose rows become soft priors with --soft-prior-strength > 0, beside the grade-C rows of the "
                             "evidence table; a literature table among them also gives the popularity bins of the test scores")
    parser.add_argument("--soft-prior-tiers", type=Path, default=DEFAULT_TIERS_PATH, help="per evidence class and direction: prior probability and unit strength")
    parser.add_argument("--soft-prior-tempering", type=float, default=DEFAULT_TEMPERING_EXPONENT, help="exponent on the capped count of independent records")
    parser.add_argument("--soft-prior-count-cap", type=int, default=DEFAULT_COUNT_CAP, help="records beyond this count add nothing")
    parser.add_argument("--soft-prior-record-grades", nargs="*", default=list(DEFAULT_RECORD_GRADES), help="grades of the evidence table that become soft priors")
    parser.add_argument("--soft-prior-probability", choices=("tier", "calibrated"), default="tier",
                        help="what a prior's target probability is: the tier table's fixed value, or (decision 9.59) a logistic regression of the label outcome on the "
                             "prior's features, fitted on the rows each fit trains on and read at the median popularity; the tiers' strengths apply either way")
    parser.add_argument("--soft-prior-calibration-minimum-cases", type=int, default=MINIMUM_CALIBRATION_CASES,
                        help="below this many positive or negative calibration cells on the fit rows the tiers are kept")
    parser.add_argument("--permute-soft-priors", action="store_true",
                        help="negative control: the soft priors moved among eligible cells so that every perturbation and symptom keeps its number of them (the "
                             "layout is lost, the popularity is kept); the real arm must beat it for the prior to carry pair-specific information")
    parser.add_argument("--edge-prior-strength", type=float, default=0.0,
                        help="global multiplier on the edge priors of decision 9.57 (graph/edge_evidence.py, models/edge_gates.py): every edge whose evidence "
                             "has a tier gets a learnable presence gate, started at its prior probability and pulled back to it with this many pseudo-counts "
                             "per unit of its tier's strength, one pseudo-count weighing as one unit of label weight. Needs <graph-dir>/edge_evidence.parquet. "
                             "0 (the default, and the registered model) builds no gate and leaves every run as it was")
    parser.add_argument("--edge-prior-tiers", type=Path, default=DEFAULT_EDGE_TIERS_PATH,
                        help="per edge family and evidence kind: prior probability that the edge is present and unit strength")
    parser.add_argument("--edge-prior-tempering", type=float, default=DEFAULT_EDGE_TEMPERING, help="exponent on an edge's capped count of distinct records")
    parser.add_argument("--edge-prior-count-cap", type=int, default=DEFAULT_EDGE_COUNT_CAP, help="an edge's records beyond this count add nothing")
    parser.add_argument("--permute-edge-priors", action="store_true",
                        help="negative control: the priors of the gated edges shuffled among the gated edges of each relation (the wiring and each relation's "
                             "distribution of priors kept, which edge holds which lost); the real arm must beat it for the evidence to carry edge-specific information")
    parser.add_argument("--edge-gate-learning-rate", type=float, default=None,
                        help="learning rate of the edge gates' logits, which take no weight decay; the main learning rate when unset")
    parser.add_argument("--edge-gate-kind", choices=EDGE_GATE_KINDS, default="hard_concrete",
                        help="with edge priors: hard_concrete (decision 9.62 (e), the owner's choice of 16 October 2026), gates that are exactly 0 or 1 at the read-out, "
                             "sampled as hard concrete while training, with the prior sigmoid(n logit p) as each edge's cost; or soft, the deterministic gates of 9.57")
    parser.add_argument("--edge-gate-beta", type=float, default=HARD_CONCRETE_BETA, help="the hard concrete temperature (Louizos et al. 2018 use 2/3)")
    parser.add_argument("--edge-gate-binary-epochs", type=int, default=8,
                        help="hard gates: epochs of the binary pass after the fit, over the gates alone, with the read-out in the forward pass and the flips' "
                             "estimated effect as the proxy gradient (models/edge_gate_passes.py); 8 is the setting the simulations of "
                             "docs/edge_gate_simulation.md tested (the pass had not settled by then: 2 to 5 flips in its last epoch, and 32 epochs at a "
                             "constant rate did not settle it either); 0 skips it")
    parser.add_argument("--edge-gate-binary-learning-rate", type=float, default=0.1, help="hard gates: the binary pass's learning rate on the gates' logits")
    parser.add_argument("--edge-gate-binary-learning-rate-decay", type=float, default=1.0,
                        help="hard gates: the binary pass's learning rate is multiplied by this after each epoch; 1 keeps it, the setting the trainer has "
                             "always had. In the simulations (docs/edge_gate_simulation.md, batch 4 of 17 October 2026) 0.8 over 32 epochs settled the "
                             "pass to no flips in all five worlds, while a constant rate left 2 to 9 flips an epoch; the closures' mean precision (0.637, "
                             "against 0.636 at a constant rate and 0.633 at 8 epochs) and the held-out loss (within 0.0006 nats in every world) barely "
                             "moved. A settled pass makes the closed set independent of the epoch it stops at")
    parser.add_argument("--edge-gate-surrogate-width", type=float, default=DEFAULT_SURROGATE_WIDTH,
                        help="hard gates: the width, in open log-odds, of the steep sigmoid the binary pass's proxy gradient passes through")
    parser.add_argument("--edge-gate-binary-full-batch", action=argparse.BooleanOptionalAction, default=True,
                        help="hard gates: one binary-pass step per epoch on the flips' effect averaged over every fit batch (no minibatch noise), or one per batch")
    parser.add_argument("--edge-gate-closure-evidence", type=int, default=20,
                        help="hard gates: after the binary pass and the checks, measure for this many closed edges, against the strongest priors first, the "
                             "exact evidence for the closure, the rise in the data loss times W when the edge alone is reopened, and P(false) = sigmoid(evidence "
                             "- logit prior); -1 every closed edge; 0 none. Each edge is a forward pass over the fit rows with the limiting-step cache rebuilt, "
                             "about 46 minutes on two CPU threads for the linear-response model on v5 (17 October 2026), hence the bound. Every closed edge "
                             "is listed in edge_gate_closures.csv whatever this is, with its first-order estimate from one more pass (decision 9.63 (c))")
    parser.add_argument("--curated-edge-labels", action=argparse.BooleanOptionalAction, default=True,
                        help="with edge priors: a curated edge with at least --curated-label-minimum-records independent records, at least "
                             "--curated-label-minimum-human-records of them human evidence (NIH iCite human score above 0) and no IntAct absence is taken as "
                             "known, with no gate (decision 9.63 (b), the owner, 16 October 2026: \"Curated edges for other species should not be considered "
                             "labels unless there's corroborating evidence from humans\"); --no-curated-edge-labels gates it at its prior, as 9.57 did")
    parser.add_argument("--curated-label-minimum-records", type=int, default=DEFAULT_LABEL_MINIMUM_RECORDS,
                        help="decision 9.63 (b): the independent records (distinct PubMed ids) a curated edge needs to be taken as known")
    parser.add_argument("--curated-label-minimum-human-records", type=int, default=DEFAULT_LABEL_MINIMUM_HUMAN_RECORDS,
                        help="decision 9.63 (b): how many of those records must be human evidence")
    parser.add_argument("--curated-label-minimum-sign-records", type=int, default=None,
                        help="decision 9.70, off unless given: a signed signalling or transcription edge is taken as known only if its sign itself rests on this "
                             "many PubMed ids and, for a regulation, CollecTRI decided the sign from a paper rather than by default (graph/sign_evidence.py; "
                             "the edge evidence table needs its sign columns, from the OmniPath evidences downloads)")
    parser.add_argument("--curated-label-exclude-disputed-signs", action="store_true",
                        help="decision 9.72, off unless given: a regulation whose two CollecTRI releases (the 2023 release and CollecTRI2) state opposite "
                             "signs, or one of which states both, is not taken as known; it stays gated at its prior (graph/sign_evidence.py, sign_disputed)")
    parser.add_argument("--curated-label-bulk-papers", choices=("count", "one", "none"), default="count",
                        help="decision 9.70 (3): how a bulk paper (a PubMed id behind more than 100 signalling or transcription edges; 114 on v5) counts "
                             "toward --curated-label-minimum-records: one record each as before (count), all of an edge's together as one (one), or "
                             "not at all (none)")
    parser.add_argument("--edge-gate-flip-checks", type=int, default=0,
                        help="hard gates: after the binary pass, flip this many of the most promising edges one at a time and keep a flip only when the exact "
                             "objective falls; each check is one pass over the fit rows")
    parser.add_argument("--equal-drug-shares", action=argparse.BooleanOptionalAction, default=True,
                        help="rescale each drug's positives and, apart, its negatives so every drug holds the mean count of each over the training drugs; evidence weights still differ between drugs; gene rows unchanged")
    parser.add_argument("--normalise-drug-input", action="store_true",
                        help="rescale each drug's seed magnitudes to sum to 1, so a drug with more targets does not push more field; genes unchanged")
    parser.add_argument("--permute-drug-targets", action="store_true",
                        help="exchange each drug's target seeds with those of another drug that has the same number of them, so every drug keeps its input size and "
                             "loses which targets it acts on (the third control for the input-size objection, docs/label_source_review.md 13.1). A negative control: "
                             "a run of it is never a result of the study")
    parser.add_argument("--drug-entry", choices=DRUG_ENTRY_MODES, default=None,
                        help="on a graph with drug entry nodes (docs/drug_entry_nodes.md): nodes seeds each drug on its own node, which is the default there; "
                             "targets drops the nodes and seeds the drug's targets, the ablation arm and the inputs of the source graph. No effect on a graph without drug nodes")
    parser.add_argument("--sequestration-carries", choices=SEQUESTRATION_MODES, default=None,
                        help="message passing on a graph with sequestration edges (a plasma carrier to the drug it binds) or clearance edges (an enzyme to the drug it clears): "
                             "change, the default there, sends how far the carrier or enzyme has moved from its unperturbed state, so the edge is silent while its source is at rest; "
                             "presence averages the source's state in like any other relation, which marks every drug that has a known carrier. The linear-response encoder carries the change by construction")
    parser.add_argument("--drug-mechanism-before-first-layer", action="store_true",
                        help="message passing with drug entry nodes: also deliver each drug's mechanism edges once before the first layer, so the drug's signal reaches as far "
                             "beyond its targets as a seed on the targets does (docs/drug_entry_nodes_measured.md, the hop)")
    parser.add_argument("--description-length-coefficient", type=float, default=1e-6)
    parser.add_argument("--checkpoint-every-minutes", type=float, default=20.0)
    parser.add_argument("--time-budget-seconds", type=float, default=0.0, help="stop training after this many seconds (0 = no limit)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--timing-batches", type=int, default=0,
                        help="run this many training steps on the training perturbations, print seconds per step and per epoch, and exit: no split directory, validation or test score")
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    return parser


def refuse_a_negative_loss_weight(arguments) -> None:
    """A loss weight below zero inverts the gradient for the pairs it applies to, so it is refused rather than run.

    The weights multiply a per-pair cross-entropy, so a weight of -0.5 makes the optimiser maximise the loss on those
    pairs: it would train the model to predict the symptom for exactly the pairs a placebo-controlled comparison
    found no excess for. That is never what is meant, and a leading minus is easy to type by accident, so it stops
    the run instead of producing a model that is wrong in a way no metric would show."""
    for option, value in (("--negative-weight", arguments.negative_weight),
                          ("--placebo-negative-weight", getattr(arguments, "placebo_negative_weight", None))):
        # the per-pair file is checked when it is read (placebo_negative_weight_matrix), because its weights are
        # not on the command line
        if value is not None and value < 0:
            raise ValueError(f"{option} is {value}; a loss weight below zero inverts the gradient for those pairs, "
                             f"so the run would train the model towards the outcome the evidence argues against. "
                             f"Use a weight of zero or more (zero sets the pairs aside)")


def main() -> None:
    arguments = build_argument_parser().parse_args()
    resolve_registered_defaults(arguments)
    refuse_a_negative_loss_weight(arguments)
    signal.signal(signal.SIGUSR1, request_checkpoint)
    torch.manual_seed(arguments.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by, label_grades=tuple(arguments.label_grades) if arguments.label_grades else None,
                                label_selection=arguments.label_selection, drug_entry=arguments.drug_entry or "nodes",
                                combined_perturbations=arguments.combined_perturbations)
    if data.drug_entry == "nodes":
        if arguments.normalise_drug_input:
            raise ValueError("--normalise-drug-input rescales seed magnitudes, and a drug seeded on its entry node has one seed of magnitude 1; "
                             "its targets' shares are the mechanism edges' weights. Use --drug-entry targets with it")
        print(f"drug entry: {int(np.asarray(data.entry_node_mask).sum())} drug nodes in the graph, each drug perturbation seeded on its own")
    placebo_negative_weights_summary = None
    if getattr(arguments, "placebo_negative_weights_file", None) is not None:
        if data.evidenced_negative is None:
            print(f"--placebo-negative-weights-file {arguments.placebo_negative_weights_file} is not read: this label selection marks no "
                  f"never-above-placebo negative, so there is no pair for a per-pair weight to act on")
        else:
            data.placebo_negative_weight_matrix, placebo_negative_weights_summary = placebo_negative_weight_matrix(
                arguments.placebo_negative_weights_file, list(data.perturbation_ids), list(data.symptoms))
            placebo_negative_weights_summary["pairs_the_weight_acts_on"] = placebo_weight_pairs_in_force(
                data.placebo_negative_weight_matrix, data.outcomes, data.evidenced_negative)
            print(f"placebo negative weights: {placebo_negative_weights_summary['pairs_matched']} of "
                  f"{placebo_negative_weights_summary['pairs_in_file']} pairs of {arguments.placebo_negative_weights_file} are in this evidence table, "
                  f"and the weight acts on {placebo_negative_weights_summary['pairs_the_weight_acts_on']} of them"
                  + (f"; not in the table: {', '.join(placebo_negative_weights_summary['pairs_not_in_this_evidence_table'])}"
                     if placebo_negative_weights_summary["pairs_not_in_this_evidence_table"] else ""))
    combined_perturbations_summary = None
    if arguments.combined_perturbations is not None:
        is_combination = np.array([kind == COMBINATION_PERTURBATION_TYPE for kind in data.perturbation_types])
        combined_perturbations_summary = {"path": str(arguments.combined_perturbations), "sha256": file_sha256(arguments.combined_perturbations),
                                          "combinations": int(is_combination.sum()),
                                          "positive_pairs": int((data.outcomes[is_combination] > 0).sum())}
        print(f"combined perturbations: {combined_perturbations_summary}")
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
    # the soft priors are built on the perturbations that remain (the lockbox never receives one) and travel with the
    # label rows under --permute-labels, so a permuted-labels run sees no pair-specific information from them either
    soft_prior_probability, soft_prior_strength, soft_priors_summary, literature_reports, calibration_cells = build_soft_priors(data, arguments)
    # the edge priors (decision 9.57) are read for the graph's own edges, before any rewiring below, which keeps each at its
    # edge position as it keeps the signs: a relation keeps its distribution of priors and its sources their priors
    edge_priors_summary = build_edge_priors(data, arguments)
    if edge_priors_summary is not None:
        print(f"edge priors: {edge_priors_summary['gated']:,} of {edge_priors_summary['edges']:,} edges gated at strength {edge_priors_summary['strength']:g}, "
              f"{edge_priors_summary['contested']} with a curated absence" + (", priors permuted within each relation" if edge_priors_summary["permutation"] else ""))
    permutation_source_rows_sha256 = None
    if arguments.permute_labels:  # rows move whole, so a pair keeps its weight, frequency and label-mask entry
        source_row = degree_stratified_row_permutation(data.perturbation_degrees_for_strata, random_seed=arguments.seed, partition=in_lockbox)
        permutation_source_rows_sha256 = array_sha256(source_row)  # experiments/score_confirmatory.py checks it against its own draw
        data.outcomes, data.weights = data.outcomes[source_row], data.weights[source_row]
        if data.frequencies is not None:
            data.frequencies = data.frequencies[source_row]
        if data.label_mask is not None:
            data.label_mask = data.label_mask[source_row]
        if data.masked_for_soft_evidence is not None:
            data.masked_for_soft_evidence = data.masked_for_soft_evidence[source_row]
        if data.evidenced_negative is not None:
            data.evidenced_negative = data.evidenced_negative[source_row]
        if soft_prior_probability is not None:
            soft_prior_probability, soft_prior_strength = soft_prior_probability[source_row], soft_prior_strength[source_row]
        if calibration_cells is not None:
            calibration_cells = permute_cell_rows(calibration_cells, source_row)
    label_mask = data.label_mask  # None without --label-selection: every pair is labelled
    permuted_drug_targets_summary = None
    if arguments.permute_drug_targets:
        if data.drug_entry == "nodes":
            raise ValueError("--permute-drug-targets exchanges a drug's seeds, and a drug seeded on its entry node has one seed, its own node; "
                             "its targets are its mechanism edges. Use --drug-entry targets with it")
        seeds, signs, magnitudes, moved = permuted_drug_targets(data.perturbation_types, data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes, arguments.seed)
        data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes = seeds, signs, magnitudes
        permuted_drug_targets_summary = {"drugs_moved": moved, "drugs": int(sum(kind == "drug" for kind in data.perturbation_types)), "seed": arguments.seed,
                                         "seeds_sha256": array_sha256(np.concatenate([np.asarray(group, dtype=np.int64) for group in seeds]))}
        print(f"permuted drug targets: {permuted_drug_targets_summary}")
    if arguments.normalise_drug_input:
        data.perturbation_magnitudes = drug_input_normalised(data.perturbation_types, data.perturbation_magnitudes, data.perturbation_parts, data.perturbation_seed_parts)
    if arguments.encoder == "linear_response":  # its input is sign x magnitude, so a perturbation whose seeds all have sign 0 gives a zero field
        zero_input = perturbations_with_a_zero_signed_input(data)
        if zero_input:
            print(f"warning: {len(zero_input)} perturbations have no seed with a nonzero sign and magnitude; the linear-response encoder gives them "
                  f"the prediction of no perturbation: {zero_input}")
    rewiring_summary = None
    if arguments.rewire_swaps_per_edge > 0:
        original_edges = np.stack([data.edge_source, data.edge_target])
        undirected_relations = reciprocated_relations(original_edges, data.edge_relation) if arguments.keep_reciprocated_relations_symmetric else []
        fixed_relations = rewiring_fixed_relations(data.relation_types, getattr(arguments, "rewire_encodes", False), has_entry_nodes=data.entry_node_mask is not None)
        rewired = fast_degree_preserving_rewiring(original_edges, data.edge_relation, num_swaps_per_edge=arguments.rewire_swaps_per_edge, random_seed=arguments.seed,
                                                  undirected_relations=undirected_relations, fixed_relations=fixed_relations)
        rewiring_summary = {"swaps_per_edge": arguments.rewire_swaps_per_edge, "seed": arguments.seed, "num_edges": int(rewired.shape[1]),
                            "share_of_edges_unchanged": float((rewired == original_edges).all(axis=0).mean()),
                            "duplicate_edges_before": duplicate_edge_count(original_edges, data.edge_relation), "duplicate_edges_after": duplicate_edge_count(rewired, data.edge_relation)}
        if arguments.keep_reciprocated_relations_symmetric:  # only then, so the split signature of earlier rewired runs is unchanged
            rewiring_summary["undirected_relations"] = [data.relation_types[relation] for relation in undirected_relations]
        if fixed_relations:  # split graphs only (docs/gene_protein_split.md), so the summary of earlier rewired runs is unchanged
            rewiring_summary["fixed_relations"] = [data.relation_types[relation] for relation in fixed_relations]
        if arguments.recompute_reaction_expression_on_rewired_graph:  # only then, so the split signature of earlier rewired runs is unchanged
            rewiring_summary["reaction_expression_recomputed"] = recompute_reaction_expression(data, original_edges, rewired, arguments)
        # The signs stay at their edge positions. The rewiring exchanges targets between two edges of one relation, so
        # position i keeps its source and its sign and gets a new target: the control randomises which node an edge
        # points at, not whether that edge activates or represses, which keeps the sign distribution of every source and
        # every relation intact and randomises only the wiring. The question did not exist before 12 October 2026,
        # when the message-passing adjacency began to read the signs (amendment item 9), and the linear-response
        # encoder has always read them this way.
        data.edge_source, data.edge_target = rewired[0], rewired[1]  # before the encoder and the adjacencies are built
        print(f"rewired graph: {rewiring_summary}")
    elif arguments.recompute_reaction_expression_on_rewired_graph:
        raise ValueError("--recompute-reaction-expression-on-rewired-graph needs --rewire-swaps-per-edge")
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
    initialise_rates_at_weighted_optimum(head, data, train_indices, label_mask, arguments, device)
    adjacencies = None  # the linear-response encoder builds its signed adjacency from the edges at construction
    if arguments.encoder == "message_passing":
        adjacencies = [adjacency.to(device) if adjacency is not None else None
                       for adjacency in message_passing_adjacencies(data, arguments.message_passing_edge_signs)]
    attach_edge_gates(encoder, data, adjacencies, arguments)
    if getattr(encoder, "edge_gates", None) is not None:
        print(f"edge gates: {encoder.edge_gates.summary()['gated_edges_the_encoder_reads']:,} of the gated edges are read by the {arguments.encoder} encoder")
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
                       "lockbox_sha256": (lockbox_summary or {}).get("sha256"), "rewiring": rewiring_summary, "labels_permuted": bool(arguments.permute_labels),
                       "permuted_drug_targets": permuted_drug_targets_summary}
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
    # a pair set aside by the selection is neither positive nor negative; under --positive-target-from-frequency the frequency is the target, not the weight
    weights = torch.as_tensor(loss_pair_weights(data, label_mask, arguments, train_indices, frequency_targets=arguments.positive_target_from_frequency), dtype=torch.float32)
    positive_targets = None
    if arguments.positive_target_from_frequency:
        frequency = np.where(np.isnan(data.frequencies), arguments.positive_target, np.maximum(data.frequencies, arguments.minimum_frequency_target))
        positive_targets = torch.as_tensor(np.where(data.outcomes > 0, frequency, 0.0), dtype=torch.float32)
    node_cost = float(np.log(len(data.node_ids)))
    soft_prior_eligible, soft_prior_permutation = None, None
    tier_probability, tier_strength = soft_prior_probability, soft_prior_strength  # before any permutation: what a phase's calibration writes on
    if soft_prior_probability is not None:
        soft_prior_eligible = eligible_cells(data.outcomes, label_mask, data.evidenced_negative, data.masked_for_soft_evidence)
        if arguments.permute_soft_priors:
            soft_prior_probability, soft_prior_strength, soft_prior_permutation = permute_soft_priors(soft_prior_probability, soft_prior_strength, soft_prior_eligible,
                                                                                                     np.arange(len(data.perturbation_ids)), random_seed=arguments.seed)
            print(f"soft priors permuted: {soft_prior_permutation}")

    calibration_summaries: dict[str, dict] = {}

    def fit_targets(fit_rows: np.ndarray, phase: str = "fit") -> tuple[dict, dict | None]:
        row_scale = None
        probability = soft_prior_probability
        if calibration_cells is not None and soft_prior_probability is not None:
            # the regression is fitted on this phase's rows alone (the lockbox is already out, a held-out fold never in)
            calibration = fit_calibration(calibration_cells, data.outcomes, label_mask, fit_rows, arguments.soft_prior_calibration_minimum_cases)
            probability = calibrated_probability_matrix(calibration, calibration_cells, tier_probability)
            if arguments.permute_soft_priors:  # the same swaps as the tiers' matrix took: the prior cells and the seed are the same
                probability, _, _ = permute_soft_priors(probability, tier_strength, soft_prior_eligible, np.arange(len(data.perturbation_ids)), random_seed=arguments.seed)
            calibration_summaries[phase] = summarise_calibration(calibration)
            print(f"soft-prior calibration ({phase}): {calibration_summaries[phase]}")
        if soft_prior_probability is not None and arguments.equal_drug_shares:
            fit_weights = loss_pair_weights(data, label_mask, arguments, fit_rows, frequency_targets=arguments.positive_target_from_frequency)
            row_scale = negative_share_scale(fit_weights, soft_prior_eligible, arguments.negative_weight)
        return fit_targets_with_soft_priors(outcomes, weights, positive_targets, fit_rows, probability, soft_prior_strength, soft_prior_eligible,
                                            arguments.soft_prior_strength, row_scale)

    fitted, soft_priors_applied_in_fit = fit_targets(train_indices)
    if soft_priors_applied_in_fit is not None:
        print(f"soft priors in the fit: {soft_priors_applied_in_fit}")
    # the edge priors' penalty is divided by the total loss weight of the pairs the fit trains on, because the BCE is a
    # weighted mean: then one pseudo-count of an edge prior weighs as one unit of label weight (models/edge_gates.py)
    edge_penalty_denominators: dict[str, float] = {}
    edge_penalty_denominator = 1.0
    if edge_priors_summary is not None:
        edge_penalty_denominator = edge_penalty_denominators["fit"] = penalty_denominator(fitted["weights"].numpy(), train_indices)

    def training_step(model_encoder, model_head, model_laboratory_readout, model_optimizer, batch) -> tuple[float, float, float, float | None]:
        """One optimizer step on a batch; returns the loss, its BCE and penalty parts and the laboratory loss (None without it)."""
        model_gates = getattr(model_encoder, "edge_gates", None)
        if isinstance(model_gates, HardConcreteEdgeGates) and model_gates.phase == "stochastic":
            model_gates.resample()  # one hard concrete sample per batch, as Louizos et al. (2018) drew it
        node_index, sign_and_magnitude = pad_perturbations(data, batch)
        field = encode(model_encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
        output = model_head(field, relation_index=0, perturbation_covariate=covariate_of(data, batch, arguments, device))
        log_parts = ({"log_probability": output.symptom_log_probability, "log_complement": output.symptom_log_complement}
                     if arguments.bce_in_log_space else {})
        if fitted["positive_targets"] is not None:
            bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, fitted["positive_targets"][batch].to(device), fitted["weights"][batch].to(device), positive_target=1.0, **log_parts)
        else:
            bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, fitted["outcomes"][batch].to(device), fitted["weights"][batch].to(device), positive_target=arguments.positive_target, **log_parts)
        penalty = arguments.description_length_coefficient * model_head.description_length_penalty(node_cost=node_cost)
        if getattr(model_encoder, "edge_gates", None) is not None:  # decision 9.57: the gates' pull back to their evidence priors
            penalty = penalty + model_encoder.edge_gates.kl_penalty() / edge_penalty_denominator
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
        if getattr(encoder, "edge_gates", None) is not None:  # the edge priors' part of train_penalty, at the end of the epoch
            entry["edge_prior_penalty"] = float(encoder.edge_gates.kl_penalty()) / edge_penalty_denominator
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
    soft_priors_applied_in_refit = None
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
        fitted, soft_priors_applied_in_refit = fit_targets(refit_indices, "refit")  # the validation rows join the fit, and so do their soft priors
        if edge_priors_summary is not None:
            edge_penalty_denominator = edge_penalty_denominators["refit"] = penalty_denominator(fitted["weights"].numpy(), refit_indices)
        refit_checkpoint_path = split_directory / "refit_checkpoint.pt"
        torch.manual_seed(arguments.seed)
        encoder, head = build_models(data, arguments, device)
        attach_edge_gates(encoder, data, adjacencies, arguments)  # fresh gates at the priors, as every other parameter starts again
        initialise_leaks_from_base_rates(head, data, refit_indices, label_mask, arguments, device)
        initialise_rates_at_weighted_optimum(head, data, refit_indices, label_mask, arguments, device)
        optimizer = torch.optim.AdamW(optimizer_parameter_groups(encoder, head, arguments), lr=arguments.learning_rate, weight_decay=arguments.weight_decay)
        refit_state = {"epoch": 0, "next_batch_start": 0, "epoch_sums": [0.0, 0.0, 0.0], "history": [], "epochs": state["best_epoch"] + 1,
                       "num_perturbations": len(refit_indices), "code_provenance": []}
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
            if getattr(encoder, "edge_gates", None) is not None:
                refit_state["history"][-1]["edge_prior_penalty"] = float(encoder.edge_gates.kl_penalty()) / edge_penalty_denominator
            refit_state.update(epoch=epoch + 1, next_batch_start=0, epoch_sums=[0.0, 0.0, 0.0])
            save_refit_checkpoint()
            print(f"refit epoch {epoch} loss {refit_state['history'][-1]['train_loss']:.4f}", flush=True)
        refit_summary = {"epochs": refit_state["epochs"], "num_perturbations": refit_state["num_perturbations"], "history": refit_state["history"],
                         "code_provenance": refit_state["code_provenance"], "early_stopped_test_predictions": "test_predictions_early_stopped.npy",
                         "early_stopped_scores": early_stopped_scores}
        if not arguments.score_lockbox:
            print(f"early-stopped model on the test: macro AUPRC {early_stopped_scores['macro_auprc']:.3f}, micro AUPRC {early_stopped_scores['micro_auprc']:.3f}")

    edge_gate_binary_pass = None
    final_gates = getattr(encoder, "edge_gates", None)
    if isinstance(final_gates, HardConcreteEdgeGates) and (arguments.edge_gate_binary_epochs > 0 or arguments.edge_gate_flip_checks > 0
                                                           or arguments.edge_gate_closure_evidence != 0):
        final_fit_rows = np.sort(refit_indices) if arguments.refit_on_validation else np.sort(train_indices)
        edge_gate_binary_pass = run_edge_gate_binary_pass(encoder, head, final_gates, final_fit_rows, fitted, data, adjacencies, arguments, device, edge_penalty_denominator)
        print(f"edge gates after the binary pass: {final_gates.summary()}")
        # every closure, with its evidence, goes beside the read-out for review (decision 9.63 (c)); results.json keeps the summary
        edge_gate_binary_pass["closure_evidence"] = write_closure_review(edge_gate_binary_pass, data, split_directory)
        edge_gate_binary_pass.pop("closed_edges")
        # the scored model is the one after the pass, so it is what the saved state holds (experiments/evaluate_on_rewired_graph.py rebuilds from it)
        if arguments.refit_on_validation:
            save_refit_checkpoint()
        else:
            state["best_encoder"], state["best_head"] = copy.deepcopy(encoder.state_dict()), copy.deepcopy(head.state_dict())
            save_checkpoint()
    predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
    clamped_probabilities = clamped_probability_report(predictions, data, test_indices, label_mask)
    print(f"clamped probabilities: {clamped_probabilities}")
    if clamped_probabilities["reopens_bce_in_log_space"]:
        print(f"warning: {clamped_probabilities['labelled_positives_below_the_clamp']} labelled positives sit below the loss's clamp of "
              f"{PROBABILITY_EPSILON}, so they carry no gradient; amendment item 19 of 12 October 2026 asks that "
              f"--bce-in-log-space be reconsidered before the confirmatory runs if this happens")
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
        "arguments": {key: json_ready(value) for key, value in vars(arguments).items()},
        "num_train": len(train_indices), "num_validation": len(validation_indices), "num_test": len(test_indices),
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
        "placebo_negative_weights": placebo_negative_weights_summary,
        "combined_perturbations": combined_perturbations_summary,
        "clamped_probabilities": clamped_probabilities,
        "propagation_steps_taken": getattr(encoder, "propagation_steps_taken", None),
        "compartment_ph_sensitivity": getattr(encoder, "compartment_ph_sensitivity_report", None),
        "message_passing_rounds_taken": getattr(encoder, "rounds_taken", None),
        "propagation_reached_the_cap": bool(getattr(encoder, "fixed_point_tolerance", None) is not None
                                            and getattr(encoder, "propagation_steps_taken", 0) >= getattr(encoder, "max_propagation_steps", 0)),
        "laboratory_labels": laboratory_results,
        "lockbox": lockbox_summary, "rewiring": rewiring_summary, "refit": refit_summary,
        "soft_priors": None if soft_priors_summary is None else {**soft_priors_summary, "applied_in_fit": soft_priors_applied_in_fit, "applied_in_refit": soft_priors_applied_in_refit,
                                                                  "calibration": calibration_summaries or None,
                                                                 "permutation": soft_prior_permutation},
        "edge_priors": None if edge_priors_summary is None else {**edge_priors_summary, "penalty_denominator": edge_penalty_denominators,
                                                                  "gates": encoder.edge_gates.summary(), "gate_values": "edge_gate_values.npy",
                                                                  **({"binary_pass": edge_gate_binary_pass, "gate_read_out": "edge_gate_read_out.npy"}
                                                                     if isinstance(encoder.edge_gates, HardConcreteEdgeGates) else {})},
    }
    if literature_reports is not None and time_split is None:
        # the popularity axis the degree bins do not measure: distinct PMIDs per perturbation in the literature table
        popularity = literature_popularity(data.perturbation_ids, literature_reports)[test_indices]
        bins = macro_auprc_by_degree_bin(predictions, test_outcomes, popularity, num_bins=5, mask=test_mask)
        results["macro_auprc_by_literature_popularity_bin"] = {label.replace("degree_bin", "popularity_bin"): value for label, value in bins.items()}
        print(f"macro AUPRC by literature popularity bin: {results['macro_auprc_by_literature_popularity_bin']}")
    if time_split_results is not None:
        results["macro_auprc"], results["macro_auroc"] = time_split_results["macro_auprc"], time_split_results["macro_auroc"]
    np.save(split_directory / "test_predictions.npy", predictions)
    if edge_priors_summary is not None:  # [num_edges] the scored model's gate per edge, NaN where an edge has none
        np.save(split_directory / "edge_gate_values.npy", encoder.edge_gates.values_per_edge())
        if isinstance(encoder.edge_gates, HardConcreteEdgeGates):  # [num_edges] 1 open, 0 closed at the read-out the test was scored with
            np.save(split_directory / "edge_gate_read_out.npy", encoder.edge_gates.open_edges().cpu().numpy().astype(np.int8))
        print(f"edge gates: {results['edge_priors']['gates']}")
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
