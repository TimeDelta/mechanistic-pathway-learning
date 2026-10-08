"""Graph reliance of trained runs: each finished fold's test perturbations scored with its trained weights on the real
graph and on a degree-preserving rewiring of it, without retraining.

The rewiring control of the confirmatory experiment (docs/preregistration.md) retrains on a rewired graph and asks
whether the real wiring is worth more than a degree-matched one. This asks a narrower question of one trained model:
how much of its score depends on the edges it was trained on. A model that reads a perturbation's own descriptors and
little else keeps its score when the edges move; a model whose predictions come through the propagation loses it. The
drop, macro and micro AUPRC on the real graph minus on the rewired graph, is the reading the descriptor treatments are
judged by (docs/preregistration.md, amendment of 8 October 2026 on the descriptor treatments). It measures dependence,
not usefulness: a model can depend on edges and still predict badly, so it is read beside the AUPRC itself.

The run is rebuilt from the arguments in its results.json (defaults of experiments/run_main_model.py for arguments
added later), its best weights are loaded from checkpoint.pt (or the refitted weights from refit_checkpoint.pt for a run
with --refit-on-validation), and the test predictions on the real graph are compared
with the stored test_predictions.npy (the reproduction check in the output). The rewiring is
fast_degree_preserving_rewiring within each relation, as --rewire-swaps-per-edge in the trainer. Each fold writes
rewired_graph_evaluation.json beside its results and is skipped when that file exists for the same swaps and seed.
Runs on permuted labels, rewired graphs, the time split or the lockbox are refused (the lockbox is blinded until the
confirmatory scoring); the development pilots of the confirmatory configurations (*_development, the lockbox removed
before the split) are allowed.

    OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/evaluate_on_rewired_graph.py \\
      --run-dirs runs/encoder/b3_typed_nodes_descriptors_brain_disease_cluster runs/b3_typed_nodes_disease_cluster \\
      --markdown-output docs/graph_reliance.md --json-output runs/graph_reliance.json
"""
from __future__ import annotations

import argparse
import json
from argparse import Namespace
from pathlib import Path

import numpy as np
import torch

from experiments.run_main_model import (build_argument_parser, build_models, macro_auprc, predict)
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.negative_controls import fast_degree_preserving_rewiring
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import micro_auprc
from mechanistic_pathway_learning.models.relational_message_passing_encoder import RelationalMessagePassingEncoder

OUTPUT_NAME = "rewired_graph_evaluation.json"
REPRODUCTION_TOLERANCE = 1e-4


def stored_arguments(results: dict) -> Namespace:
    """The trainer's arguments of a finished run: parser defaults, overwritten by the values in its results.json, with
    path arguments turned back into paths."""
    parser = build_argument_parser()
    defaults = vars(parser.parse_args([]))
    path_arguments = {action.dest for action in parser._actions if action.type is Path}
    stored = {key: (Path(value) if key in path_arguments and value is not None else value) for key, value in results["arguments"].items()}
    return Namespace(**{**defaults, **stored})


def refuse_unsupported(results: dict, arguments: Namespace, fold_directory: Path) -> None:
    if arguments.score_lockbox or fold_directory.parent.name.endswith("_confirmatory") or "lockbox" in fold_directory.name:
        raise SystemExit(f"{fold_directory}: lockbox runs are blinded until the confirmatory scoring")
    if results.get("labels_permuted") or results.get("rewiring") or arguments.time_split_cutoff is not None:
        raise SystemExit(f"{fold_directory}: graph reliance is defined for runs on the real labels, the real graph and the grouped split")
    if arguments.laboratory_label_weight > 0:
        print(f"{fold_directory}: the laboratory readout is not scored here; only the symptom predictions are")


def experiment_data(arguments: Namespace):
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by,
                                label_grades=tuple(arguments.label_grades) if arguments.label_grades else None, label_selection=arguments.label_selection)
    if arguments.lockbox is not None:
        data = restrict_to_perturbations(data, ~read_lockbox(arguments.lockbox, data))
    return data


def trained_models(data, arguments: Namespace, checkpoint: dict, device):
    """A fresh encoder and head on data's current edges with the run's best weights (a fresh encoder, so no adjacency
    cached for other edges is reused)."""
    encoder, head = build_models(data, arguments, device)
    state = checkpoint["state"]
    encoder.load_state_dict(state["best_encoder"] if state.get("best_encoder") is not None else checkpoint["encoder"])
    head.load_state_dict(state["best_head"] if state.get("best_head") is not None else checkpoint["head"])
    adjacencies = None
    if arguments.encoder == "message_passing":
        adjacencies = [adjacency.to(device) if adjacency is not None else None for adjacency in RelationalMessagePassingEncoder.build_relation_adjacencies(
            torch.as_tensor(np.stack([data.edge_source, data.edge_target]), dtype=torch.long), torch.as_tensor(data.edge_relation, dtype=torch.long),
            len(data.node_ids), len(data.relation_types))]
    return encoder, head, adjacencies


def evaluate_fold(fold_directory: Path, swaps_per_edge: int, rewiring_seed: int, data_cache: dict) -> dict:
    output_path = fold_directory / OUTPUT_NAME
    if output_path.exists():
        existing = json.loads(output_path.read_text())
        if existing.get("swaps_per_edge") == swaps_per_edge and existing.get("rewiring_seed") == rewiring_seed:
            return existing
    results = json.loads((fold_directory / "results.json").read_text())
    arguments = stored_arguments(results)
    refuse_unsupported(results, arguments, fold_directory)
    device = torch.device("cpu")
    data_key = (str(arguments.graph_dir), str(arguments.evidence_dir), arguments.group_by, tuple(arguments.label_grades or ()), str(arguments.label_selection),
                str(arguments.lockbox), arguments.metabolic_layer_only)
    if data_key not in data_cache:
        data_cache.clear()  # one data set in memory at a time
        data_cache[data_key] = experiment_data(arguments)
    data = data_cache[data_key]
    position = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    test_indices = np.array([position[perturbation_id] for perturbation_id in results["test_perturbation_ids"]])
    test_outcomes = data.outcomes[test_indices]
    test_mask = None if data.label_mask is None else data.label_mask[test_indices]
    if results.get("refit"):  # the scored model is the refit (--refit-on-validation), whose last state is the scored one
        refit_checkpoint = torch.load(fold_directory / "refit_checkpoint.pt", map_location=device, weights_only=False)
        checkpoint = {"state": {}, "encoder": refit_checkpoint["encoder"], "head": refit_checkpoint["head"]}
    else:
        checkpoint = torch.load(fold_directory / "checkpoint.pt", map_location=device, weights_only=False)
    torch.manual_seed(arguments.seed)

    encoder, head, adjacencies = trained_models(data, arguments, checkpoint, device)
    real_predictions = predict(encoder, head, data, test_indices, adjacencies, arguments, device)
    stored_predictions = np.load(fold_directory / "test_predictions.npy")
    reproduction_difference = float(np.abs(real_predictions - stored_predictions).max())

    original_source, original_target = data.edge_source.copy(), data.edge_target.copy()
    original_edges = np.stack([original_source, original_target])
    rewired = fast_degree_preserving_rewiring(original_edges, data.edge_relation, num_swaps_per_edge=swaps_per_edge, random_seed=rewiring_seed)
    try:
        data.edge_source, data.edge_target = rewired[0], rewired[1]
        rewired_encoder, rewired_head, rewired_adjacencies = trained_models(data, arguments, checkpoint, device)
        rewired_predictions = predict(rewired_encoder, rewired_head, data, test_indices, rewired_adjacencies, arguments, device)
    finally:
        data.edge_source, data.edge_target = original_source, original_target  # the cached data stays the real graph

    evaluation = {
        "swaps_per_edge": swaps_per_edge, "rewiring_seed": rewiring_seed, "share_of_edges_unchanged": float((rewired == original_edges).all(axis=0).mean()),
        "num_test": int(len(test_indices)), "reproduction_max_abs_difference": reproduction_difference,
        "reproduced": reproduction_difference <= REPRODUCTION_TOLERANCE,
        "macro_auprc_real": macro_auprc(real_predictions, test_outcomes, test_mask), "micro_auprc_real": micro_auprc(real_predictions, test_outcomes, test_mask),
        "macro_auprc_rewired": macro_auprc(rewired_predictions, test_outcomes, test_mask), "micro_auprc_rewired": micro_auprc(rewired_predictions, test_outcomes, test_mask),
        "stored_macro_auprc": results["macro_auprc"], "stored_micro_auprc": results.get("micro_auprc"),
    }
    evaluation["macro_reliance"] = evaluation["macro_auprc_real"] - evaluation["macro_auprc_rewired"]
    evaluation["micro_reliance"] = evaluation["micro_auprc_real"] - evaluation["micro_auprc_rewired"]
    output_path.write_text(json.dumps(evaluation, indent=2))
    return evaluation


def summarise(configuration_directory: Path, evaluations: list[dict]) -> dict:
    def mean_and_sd(key: str) -> tuple[float, float]:
        values = np.array([evaluation[key] for evaluation in evaluations], dtype=float)
        return float(values.mean()), float(values.std(ddof=1)) if len(values) > 1 else float("nan")

    summary = {"configuration": configuration_directory.name, "folds": len(evaluations), "all_reproduced": all(evaluation["reproduced"] for evaluation in evaluations),
               "largest_reproduction_difference": max(evaluation["reproduction_max_abs_difference"] for evaluation in evaluations)}
    for key in ("macro_auprc_real", "macro_auprc_rewired", "macro_reliance", "micro_auprc_real", "micro_auprc_rewired", "micro_reliance"):
        summary[key] = mean_and_sd(key)
    summary["per_fold_macro_reliance"] = [evaluation["macro_reliance"] for evaluation in evaluations]
    summary["per_fold_micro_reliance"] = [evaluation["micro_reliance"] for evaluation in evaluations]
    return summary


def markdown_table(summaries: list[dict], swaps_per_edge: int, rewiring_seed: int) -> str:
    lines = [
        "# Graph reliance",
        "",
        "Generated by experiments/evaluate_on_rewired_graph.py. Each fold's test perturbations scored with the trained weights on the real graph and on "
        f"a degree-preserving rewiring of it within each relation ({swaps_per_edge} attempted swaps per edge, seed {rewiring_seed}), without retraining. "
        "Reliance is real minus rewired; mean ± standard deviation over folds. It measures how much a trained model's score depends on its edges, not "
        "whether the real edges are better than degree-matched ones (that is the retrained rewiring control).",
        "",
        "The reproduction column is the largest absolute difference between a stored test prediction and the same prediction rebuilt on the real graph; "
        f"above {REPRODUCTION_TOLERANCE:g} the run was trained under earlier code whose arithmetic differs (for example the summation order of the propagation, changed on 8 October), "
        "and both of its scores here come from the current code.",
        "",
        "| configuration | folds | reproduction | macro real | macro rewired | macro reliance | micro real | micro rewired | micro reliance |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for summary in summaries:
        cells = [f"{summary[key][0]:.3f} ± {summary[key][1]:.3f}" for key in
                 ("macro_auprc_real", "macro_auprc_rewired", "macro_reliance", "micro_auprc_real", "micro_auprc_rewired", "micro_reliance")]
        lines.append(f"| {summary['configuration']} | {summary['folds']} | {summary['largest_reproduction_difference']:.0e} | " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True, help="configuration directories holding fold*_seed0 subdirectories")
    parser.add_argument("--swaps-per-edge", type=int, default=50, help="attempted swaps per edge (50, as the confirmatory rewiring control)")
    parser.add_argument("--rewiring-seed", type=int, default=0)
    parser.add_argument("--markdown-output", type=Path, default=None)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()
    summaries, data_cache = [], {}
    for configuration_directory in arguments.run_dirs:
        fold_directories = sorted(path for path in configuration_directory.glob("fold*_seed0") if (path / "DONE").exists())
        if not fold_directories:
            print(f"{configuration_directory}: no finished folds; skipped")
            continue
        evaluations = []
        for fold_directory in fold_directories:
            evaluation = evaluate_fold(fold_directory, arguments.swaps_per_edge, arguments.rewiring_seed, data_cache)
            print(f"{fold_directory}: macro {evaluation['macro_auprc_real']:.3f} -> {evaluation['macro_auprc_rewired']:.3f}, "
                  f"micro {evaluation['micro_auprc_real']:.3f} -> {evaluation['micro_auprc_rewired']:.3f}, reproduction difference {evaluation['reproduction_max_abs_difference']:.1e}",
                  flush=True)
            evaluations.append(evaluation)
        summaries.append(summarise(configuration_directory, evaluations))
    if arguments.json_output:
        arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_output.write_text(json.dumps(summaries, indent=2))
    table = markdown_table(summaries, arguments.swaps_per_edge, arguments.rewiring_seed)
    if arguments.markdown_output:
        arguments.markdown_output.write_text(table)
    print(table)


if __name__ == "__main__":
    main()
