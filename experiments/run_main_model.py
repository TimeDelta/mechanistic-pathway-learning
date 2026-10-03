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
import json
import signal
import time
from pathlib import Path

import numpy as np
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
from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
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


def build_models(data, arguments, device):
    encoder = RelationalMessagePassingEncoder(len(data.node_ids), len(data.relation_types), arguments.node_state_dim, arguments.num_layers).to(device)
    if arguments.head == "sigmoid":
        head = RelationalGnnSigmoidHead(arguments.node_state_dim, len(data.symptoms), hidden_dim=arguments.sigmoid_hidden_dim, pooling=arguments.pooling).to(device)
    else:
        head = NoisyOrPathwayModuleHead(len(data.node_ids), arguments.node_state_dim, arguments.num_modules, len(data.symptoms),
                                        gate_initial_log_alpha=arguments.gate_initial_log_alpha, pooling=arguments.pooling).to(device)
    return encoder, head


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
            predictions.append(head(field, relation_index=0).symptom_probability.cpu().numpy())
    encoder.train()
    head.train()
    return np.concatenate(predictions, axis=0) if predictions else np.zeros((0, len(data.symptoms)))


def macro_auprc(predictions: np.ndarray, outcomes: np.ndarray) -> float:
    values = [per_symptom_auprc(predictions, outcomes, s) for s in range(outcomes.shape[1]) if MINIMUM_POSITIVES_TO_SCORE <= outcomes[:, s].sum() < outcomes.shape[0]]
    return float(np.mean(values)) if values else float("nan")


def split_indices(data, arguments) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    """Return (train, validation, test) index arrays and a name for the split directory."""
    all_indices = np.arange(len(data.perturbation_ids))
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
    parser.add_argument("--group-by", choices=["gene", "disease_cluster"], default="gene")
    parser.add_argument("--validation-fraction", type=float, default=0.15)
    parser.add_argument("--patience", type=int, default=8)
    parser.add_argument("--selection-metric", choices=["loss", "auprc"], default="loss",
                        help="early stopping on the validation evidence-weighted BCE (smooth on small validation sets) or on validation macro AUPRC")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--head", choices=["noisy_or", "sigmoid"], default="noisy_or")
    parser.add_argument("--field", choices=["difference", "absolute"], default="difference")
    parser.add_argument("--pooling", choices=["sum", "mean"], default="sum")
    parser.add_argument("--node-state-dim", type=int, default=32)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--num-modules", type=int, default=8)
    parser.add_argument("--sigmoid-hidden-dim", type=int, default=64)
    parser.add_argument("--gate-initial-log-alpha", type=float, default=-1.0)
    parser.add_argument("--init-leak-from-base-rate", action="store_true", help="noisy-OR head: start each symptom's leak at its training base rate")
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
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by)
    if arguments.permute_labels:
        data.outcomes = permute_symptom_labels_within_degree_strata(data.outcomes, data.perturbation_degrees, random_seed=arguments.seed)
    train_indices, validation_indices, test_indices, split_name = split_indices(data, arguments)
    if arguments.permute_labels:
        split_name += "_permuted"
    encoder, head = build_models(data, arguments, device)
    if arguments.init_leak_from_base_rate and arguments.head == "noisy_or":
        head.initialize_leak_from_base_rates(torch.as_tensor(data.outcomes[train_indices].mean(axis=0), dtype=torch.float32, device=device))
    adjacencies = [adjacency.to(device) if adjacency is not None else None for adjacency in RelationalMessagePassingEncoder.build_relation_adjacencies(
        torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long), torch.as_tensor(data.edge_relation, dtype=torch.long), len(data.node_ids), len(data.relation_types))]
    optimizer = torch.optim.AdamW(list(encoder.parameters()) + list(head.parameters()), lr=arguments.learning_rate, weight_decay=arguments.weight_decay)
    split_directory = arguments.run_dir / split_name
    split_directory.mkdir(parents=True, exist_ok=True)
    checkpoint_path = split_directory / "checkpoint.pt"
    state = {"epoch": 0, "history": [], "best_validation_auprc": float("-inf"), "best_validation_loss": float("inf"), "best_epoch": -1, "epochs_without_improvement": 0, "best_encoder": None, "best_head": None}
    if arguments.resume and checkpoint_path.exists():
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        encoder.load_state_dict(checkpoint["encoder"])
        head.load_state_dict(checkpoint["head"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        state = checkpoint["state"]
        print(f"resumed at epoch {state['epoch']}")

    def save_checkpoint() -> None:
        torch.save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "state": state}, checkpoint_path)

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
        epoch_loss, epoch_bce, epoch_penalty = 0.0, 0.0, 0.0
        for start in range(0, len(order), arguments.batch_size):
            batch = order[start : start + arguments.batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)
            output = head(field, relation_index=0)
            if positive_targets is not None:
                bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, positive_targets[batch].to(device), weights[batch].to(device), positive_target=1.0)
            else:
                bce = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device), positive_target=arguments.positive_target)
            penalty = arguments.description_length_coefficient * head.description_length_penalty(node_cost=node_cost)
            loss = bce + penalty
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
        if len(validation_indices):
            validation_predictions = predict(encoder, head, data, validation_indices, adjacencies, arguments, device)
            entry["validation_macro_auprc"] = macro_auprc(validation_predictions, data.outcomes[validation_indices])
            entry["validation_loss"] = float(evidence_weighted_binary_cross_entropy(torch.as_tensor(validation_predictions, dtype=torch.float32), outcomes[validation_indices], weights[validation_indices], positive_target=arguments.positive_target))
            improved = (entry["validation_loss"] < state["best_validation_loss"] - 1e-5) if arguments.selection_metric == "loss" else (entry["validation_macro_auprc"] > state["best_validation_auprc"] + 1e-4)
            if improved:
                state.update(best_validation_auprc=entry["validation_macro_auprc"], best_validation_loss=entry["validation_loss"], best_epoch=epoch, epochs_without_improvement=0,
                             best_encoder=copy.deepcopy(encoder.state_dict()), best_head=copy.deepcopy(head.state_dict()))
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

    predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
    test_outcomes = data.outcomes[test_indices]
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
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
        "arguments": {key: (str(value) if isinstance(value, Path) else value) for key, value in vars(arguments).items()},
        "num_train": int(len(train_indices)), "num_validation": int(len(validation_indices)), "num_test": int(len(test_indices)),
        "epochs_completed": state["epoch"], "best_epoch": state["best_epoch"], "stopped_early": stopped_early, "history": state["history"],
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, test_outcomes), "hits_at_3": hits_at_k(predictions, test_outcomes, 3),
        "expected_calibration_error": expected_calibration_error(predictions, test_outcomes), "per_symptom": per_symptom, "symptoms": data.symptoms,
        "test_perturbation_ids": [data.perturbation_ids[i] for i in test_indices],
    }
    np.save(split_directory / "test_predictions.npy", predictions)
    if arguments.head == "noisy_or":
        with torch.no_grad():
            head.eval()
            support = head.module_support().cpu().numpy()
            expected_support = head.support_gate.expected_active_node_count().cpu().numpy()
            links = head.link_probability(0).cpu().numpy()
            leaks = head.leak_probability(0).cpu().numpy()
            test_field_activations = []
            for start in range(0, len(test_indices), arguments.batch_size):
                batch = test_indices[start : start + arguments.batch_size]
                node_index, sign_and_magnitude = pad_perturbations(data, batch)
                test_field_activations.append(head(encode(encoder, node_index.to(device), sign_and_magnitude.to(device), adjacencies, arguments.field)).module_activation.cpu().numpy())
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
