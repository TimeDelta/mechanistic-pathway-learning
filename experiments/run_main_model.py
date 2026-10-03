"""Train and evaluate the proposed model B6 under grouped cross-validation (design sections 5 and 6).

Resumable for a runtime-limited cluster: a checkpoint is written every
--checkpoint-every-minutes and on SIGUSR1; --resume continues from it; a DONE
marker is written when every fold finishes, which slurm/train_resumable.sbatch
waits for before it stops requeueing. One fold per invocation (--fold) lets
folds and seeds run as array jobs.

Training objective per batch of perturbations:
  evidence-weighted BCE (positives smoothed to 0.99; unobserved pairs are weak negatives)
  + description-length coefficient * model cost (expected active support nodes and links)
"""
from __future__ import annotations

import argparse
import json
import signal
import time
from pathlib import Path

import numpy as np
import torch

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    expected_calibration_error,
    hits_at_k,
    mean_reciprocal_rank,
    per_symptom_auprc,
    per_symptom_auroc,
)
from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder
from mechanistic_pathway_learning.models.soft_constraint_losses import evidence_weighted_binary_cross_entropy

CHECKPOINT_REQUESTED = False


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
    head = NoisyOrPathwayModuleHead(len(data.node_ids), arguments.node_state_dim, arguments.num_modules, len(data.symptoms), gate_initial_log_alpha=arguments.gate_initial_log_alpha).to(device)
    return encoder, head


def predict(encoder, head, data, indices: np.ndarray, adjacencies, batch_size: int, device) -> np.ndarray:
    encoder.eval()
    head.eval()
    predictions = []
    with torch.no_grad():
        for start in range(0, len(indices), batch_size):
            batch = indices[start : start + batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encoder(node_index.to(device), sign_and_magnitude.to(device), relation_adjacencies=adjacencies)
            predictions.append(head(field, relation_index=0).symptom_probability.cpu().numpy())
    encoder.train()
    head.train()
    return np.concatenate(predictions, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--run-dir", type=Path, default=Path("runs/main_model"))
    parser.add_argument("--fold", type=int, default=0)
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--node-state-dim", type=int, default=32)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--num-modules", type=int, default=8)
    parser.add_argument("--gate-initial-log-alpha", type=float, default=-1.0)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--learning-rate", type=float, default=0.002)
    parser.add_argument("--max-epochs", type=int, default=30)
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
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_perturbation = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
    train_indices = np.where(fold_of_perturbation != arguments.fold)[0]
    test_indices = np.where(fold_of_perturbation == arguments.fold)[0]
    encoder, head = build_models(data, arguments, device)
    adjacencies = [adjacency.to(device) if adjacency is not None else None for adjacency in RelationalMessagePassingEncoder.build_relation_adjacencies(
        torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long), torch.as_tensor(data.edge_relation, dtype=torch.long), len(data.node_ids), len(data.relation_types))]
    optimizer = torch.optim.AdamW(list(encoder.parameters()) + list(head.parameters()), lr=arguments.learning_rate, weight_decay=1e-4)
    fold_directory = arguments.run_dir / f"fold{arguments.fold}_seed{arguments.seed}"
    fold_directory.mkdir(parents=True, exist_ok=True)
    checkpoint_path = fold_directory / "checkpoint.pt"
    start_epoch = 0
    if arguments.resume and checkpoint_path.exists():
        checkpoint = torch.load(checkpoint_path, map_location=device)
        encoder.load_state_dict(checkpoint["encoder"])
        head.load_state_dict(checkpoint["head"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        start_epoch = checkpoint["epoch"] + 1
        print(f"resumed at epoch {start_epoch}")
    outcomes = torch.as_tensor(data.outcomes, dtype=torch.float32)
    weights = torch.as_tensor(np.where(data.outcomes > 0, np.maximum(data.weights, 1e-3), arguments.negative_weight), dtype=torch.float32)
    started = time.time()
    last_checkpoint = time.time()
    generator = np.random.default_rng(arguments.seed)
    history = []
    for epoch in range(start_epoch, arguments.max_epochs):
        order = generator.permutation(train_indices)
        epoch_loss = 0.0
        for start in range(0, len(order), arguments.batch_size):
            batch = order[start : start + arguments.batch_size]
            node_index, sign_and_magnitude = pad_perturbations(data, batch)
            field = encoder(node_index.to(device), sign_and_magnitude.to(device), relation_adjacencies=adjacencies)
            output = head(field, relation_index=0)
            loss = evidence_weighted_binary_cross_entropy(output.symptom_probability, outcomes[batch].to(device), weights[batch].to(device))
            loss = loss + arguments.description_length_coefficient * head.description_length_penalty(node_cost=float(np.log(len(data.node_ids))))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item() * len(batch)
            elapsed_minutes = (time.time() - last_checkpoint) / 60.0
            if CHECKPOINT_REQUESTED or elapsed_minutes >= arguments.checkpoint_every_minutes:
                torch.save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "epoch": epoch - 1}, checkpoint_path)
                last_checkpoint = time.time()
                if CHECKPOINT_REQUESTED:
                    print("checkpoint written on signal; exiting for requeue")
                    return
            if arguments.time_budget_seconds and time.time() - started > arguments.time_budget_seconds:
                break
        history.append({"epoch": epoch, "train_loss": epoch_loss / max(1, len(order))})
        print(f"epoch {epoch} loss {history[-1]['train_loss']:.4f} elapsed {time.time() - started:.0f}s")
        torch.save({"encoder": encoder.state_dict(), "head": head.state_dict(), "optimizer": optimizer.state_dict(), "epoch": epoch}, checkpoint_path)
        if arguments.time_budget_seconds and time.time() - started > arguments.time_budget_seconds:
            print("time budget reached")
            break
    predictions = predict(encoder, head, data, test_indices, adjacencies, arguments.batch_size, device)
    test_outcomes = data.outcomes[test_indices]
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        if test_outcomes[:, symptom_index].sum() < 5:
            continue
        per_symptom[symptom] = {
            "auprc": bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap).__dict__,
            "auroc": bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, test_outcomes, num_bootstrap=arguments.num_bootstrap).__dict__,
            "base_rate": float(test_outcomes[:, symptom_index].mean()),
        }
    with torch.no_grad():
        head.eval()
        support = head.module_support().cpu().numpy()
        links = head.link_probability(0).cpu().numpy()
    results = {
        "fold": arguments.fold, "seed": arguments.seed, "num_train": int(len(train_indices)), "num_test": int(len(test_indices)),
        "epochs_completed": len(history) + start_epoch, "history": history,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, test_outcomes), "hits_at_3": hits_at_k(predictions, test_outcomes, 3),
        "expected_calibration_error": expected_calibration_error(predictions, test_outcomes), "per_symptom": per_symptom,
        "module_support_sizes": [int((support[k] > 0.5).sum()) for k in range(support.shape[0])],
        "module_symptom_links": links.round(3).tolist(), "symptoms": data.symptoms,
    }
    (fold_directory / "results.json").write_text(json.dumps(results, indent=1))
    np.save(fold_directory / "module_support.npy", support)
    print(f"fold {arguments.fold}: macro AUPRC {results['macro_auprc']:.3f} macro AUROC {results['macro_auroc']:.3f} MRR {results['mean_reciprocal_rank']:.3f} hits@3 {results['hits_at_3']:.3f}")
    print("module support sizes:", results["module_support_sizes"])
    (fold_directory / "DONE").write_text("done\n")


if __name__ == "__main__":
    main()
