"""Diagnose why noisy-OR modules go unused: optimizer travel of the links, gate movement and the gap between the
training gate and the evaluation gate, read from the saved checkpoints and test activations of finished runs.

For each split it reports how far the link logits, gate log-alphas and readout biases moved from their initial values
against the largest distance Adam could move them (about one learning rate per optimizer step up to the best epoch),
the evaluation and expected training gate at the median log-alpha, and how much the test-set module activations vary
across perturbations and agree between modules.

Usage:
  python experiments/diagnose_noisy_or_modules.py --run-dirs runs/b6_default_disease_cluster runs/b6_mechanistic_disease_cluster \
      --markdown-output docs/b6_module_diagnosis.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import HardConcreteNodeGate

INITIAL_LINK_LOGIT = -3.0  # NoisyOrPathwayModuleHead default


def gate_values_at(log_alpha: float) -> tuple[float, float]:
    """(Louizos test-time gate, expected training gate) at one log-alpha."""
    gate = HardConcreteNodeGate(1, 1, initial_log_alpha=log_alpha, initial_log_alpha_noise=0.0)
    with torch.no_grad():
        return float(gate.louizos_test_time_gate()[0, 0]), float(gate.expected_gate()[0, 0])


def evaluation_gate_of_run(results: dict) -> tuple[str, str]:
    """(evaluation gate, commits) of a finished run. Runs that record their commit evaluate with the expected training
    gate; runs from before that, without a commit, used the Louizos estimator unless their arguments say otherwise
    (the diagnosis runs of commit f6985fc passed --gate-evaluation expected)."""
    commits = sorted({entry["commit"][:7] for entry in results.get("code_provenance", []) if entry.get("commit")})
    if commits or results.get("arguments", {}).get("gate_evaluation") == "expected":
        return "expected", ", ".join(commits) or "f6985fc"
    return "louizos", "not recorded"


def diagnose_split(split_directory: Path) -> dict:
    checkpoint = torch.load(split_directory / "checkpoint.pt", map_location="cpu", weights_only=False)
    results = json.loads((split_directory / "results.json").read_text())
    configuration = results.get("arguments", {})
    state = checkpoint["state"]
    best_head = state["best_head"]
    learning_rate = float(configuration.get("learning_rate", 0.002))
    scalar_learning_rate = float(configuration.get("link_learning_rate") or configuration.get("head_scalar_learning_rate") or learning_rate)  # head_scalar_learning_rate: runs of commits f6985fc to 1cf0dcb
    batch_size = int(configuration.get("batch_size", 16))
    train_size = int(results.get("num_train") or 0)
    steps_per_epoch = int(np.ceil(train_size / batch_size)) if train_size else None
    optimizer_steps_to_best = steps_per_epoch * (state["best_epoch"] + 1) if steps_per_epoch else None
    link_logit = best_head["module_symptom_link_logit"][0].numpy()
    log_alpha = best_head["support_gate.log_alpha"].numpy()
    readout_bias = best_head["module_readout_bias"].numpy()
    activations = np.load(split_directory / "test_module_activations.npy")
    off_diagonal = ~np.eye(activations.shape[1], dtype=bool)
    median_log_alpha = float(np.median(log_alpha))
    return {
        "split": split_directory.name,
        "evaluation_gate": evaluation_gate_of_run(results)[0],
        "code_commits": evaluation_gate_of_run(results)[1],
        "best_epoch": state["best_epoch"],
        "optimizer_steps_to_best": optimizer_steps_to_best,
        "largest_link_travel_possible": None if optimizer_steps_to_best is None else scalar_learning_rate * optimizer_steps_to_best,
        "largest_link_travel": float(np.abs(link_logit - INITIAL_LINK_LOGIT).max()),
        "largest_link_probability": float(1.0 / (1.0 + np.exp(-link_logit.max()))),
        "readout_bias_mean": float(readout_bias.mean()),
        "median_log_alpha": median_log_alpha,
        "log_alpha_99th_percentile": float(np.percentile(log_alpha, 99)),
        "louizos_gate_at_median": gate_values_at(median_log_alpha)[0],
        "expected_training_gate_at_median": gate_values_at(median_log_alpha)[1],
        "activation_mean": float(activations.mean()),
        "activation_standard_deviation_across_perturbations": float(activations.std(axis=0).mean()),
        "activation_correlation_between_modules": float(np.mean(np.corrcoef(activations.T)[off_diagonal])) if activations.std(axis=0).min() > 0 else float("nan"),
        "test_macro_auprc": results.get("macro_auprc"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--markdown-output", type=Path, default=None)
    arguments = parser.parse_args()
    rows = []
    for run_directory in arguments.run_dirs:
        for split_directory in sorted(path for path in run_directory.iterdir() if (path / "DONE").exists()):
            rows.append({"run": run_directory.name, **diagnose_split(split_directory)})
    columns = list(rows[0]) if rows else []
    lines = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    for row in rows:
        lines.append("| " + " | ".join(f"{value:.3f}" if isinstance(value, float) else str(value) for value in row.values()) + " |")
    table = "\n".join(lines)
    print(table)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(
            "# Noisy-OR module diagnosis (generated by experiments/diagnose_noisy_or_modules.py)\n\n"
            "Links start at logit -3 and Adam moves a parameter by about one learning rate per step, so largest_link_travel_possible "
            "bounds how far any link could have moved by the best epoch. louizos_gate_at_median is the test-time estimator of Louizos et al. "
            "at the median log-alpha and expected_training_gate_at_median the mean gate the training passes saw there; evaluation_gate "
            "names the one the run evaluated with, and code_commits the commits it ran under (runs before commit recording used the "
            "Louizos estimator).\n\n" + table + "\n")


if __name__ == "__main__":
    main()
