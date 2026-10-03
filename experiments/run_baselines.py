"""Phase 2: baselines B0 (popularity, degree-scaled popularity) and B1 (random walk with restart)
under grouped perturbation-wise cross-validation (design sections 5.6 and 6.1).

Out-of-fold predictions are pooled and scored per symptom (AUPRC, AUROC with bootstrap
intervals) plus mean reciprocal rank and hits-at-3 per perturbation. Writes
runs/baselines/results.json and docs/phase2_baselines.md.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    hits_at_k,
    mean_reciprocal_rank,
    per_symptom_auprc,
    per_symptom_auroc,
)
from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline
from mechanistic_pathway_learning.models.baselines.random_walk_with_restart_baseline import (
    RandomWalkWithRestartBaseline,
    build_normalized_adjacency,
)


def out_of_fold_predictions(data, fold_of_perturbation: np.ndarray, num_folds: int, model_name: str, restart_probability: float) -> np.ndarray:
    predictions = np.zeros_like(data.outcomes)
    degrees = data.perturbation_degrees
    normalized_adjacency = None
    if model_name == "random_walk_with_restart":
        normalized_adjacency = build_normalized_adjacency(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0])
    for fold in range(num_folds):
        train = fold_of_perturbation != fold
        test = ~train
        if model_name in ("popularity", "degree_popularity"):
            model = PopularityBaseline(scale_by_degree=model_name == "degree_popularity").fit(data.outcomes[train])
            predictions[test] = model.predict(int(test.sum()), degrees[test])
        elif model_name == "random_walk_with_restart":
            model = RandomWalkWithRestartBaseline(restart_probability=restart_probability).fit(
                normalized_adjacency, [data.perturbation_seeds[i] for i in np.where(train)[0]], data.outcomes[train]
            )
            predictions[test] = model.predict([data.perturbation_seeds[i] for i in np.where(test)[0]])
        else:
            raise ValueError(model_name)
    return predictions


def score(data, predictions: np.ndarray, num_bootstrap: int) -> dict:
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        if data.outcomes[:, symptom_index].sum() < 5:
            continue
        auprc = bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, data.outcomes, num_bootstrap=num_bootstrap)
        auroc = bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, data.outcomes, num_bootstrap=num_bootstrap)
        per_symptom[symptom] = {"positives": int(data.outcomes[:, symptom_index].sum()), "auprc": auprc.__dict__, "auroc": auroc.__dict__, "base_rate": float(data.outcomes[:, symptom_index].mean())}
    return {
        "per_symptom": per_symptom,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, data.outcomes),
        "hits_at_3": hits_at_k(predictions, data.outcomes, 3),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--restart-probability", type=float, default=0.3)
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/baselines"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase2_baselines.md"))
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_perturbation = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
    results = {"num_perturbations": len(data.perturbation_ids), "num_symptoms": len(data.symptoms), "symptoms": data.symptoms, "models": {}}
    for model_name in ("popularity", "degree_popularity", "random_walk_with_restart"):
        predictions = out_of_fold_predictions(data, fold_of_perturbation, arguments.num_folds, model_name, arguments.restart_probability)
        results["models"][model_name] = score(data, predictions, arguments.num_bootstrap)
        print(f"{model_name:26s} macro AUPRC {results['models'][model_name]['macro_auprc']:.3f}  macro AUROC {results['models'][model_name]['macro_auroc']:.3f}  MRR {results['models'][model_name]['mean_reciprocal_rank']:.3f}  hits@3 {results['models'][model_name]['hits_at_3']:.3f}")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "results.json").write_text(json.dumps(results, indent=1))
    lines = ["# Phase 2 baselines (generated by experiments/run_baselines.py)", "",
             f"{len(data.perturbation_ids)} perturbations ({sum(t == 'gene' for t in data.perturbation_types)} genes, {sum(t == 'drug' for t in data.perturbation_types)} drugs), {len(data.symptoms)} symptoms, relation induces, {arguments.num_folds}-fold grouped cross-validation (drugs grouped by dominant target), out-of-fold predictions pooled. Unobserved pairs count as negatives (positive-unlabelled convention).", "",
             "| model | macro AUPRC | macro AUROC | MRR | hits@3 |", "|---|---|---|---|---|"]
    for model_name, entry in results["models"].items():
        lines.append(f"| {model_name} | {entry['macro_auprc']:.3f} | {entry['macro_auroc']:.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} |")
    lines += ["", "Per-symptom AUPRC (point estimate with 95 percent bootstrap interval over perturbations); base rate in parentheses.", "", "| symptom | " + " | ".join(results["models"]) + " |", "|---|" + "---|" * len(results["models"])]
    for symptom in data.symptoms:
        cells = []
        for entry in results["models"].values():
            s = entry["per_symptom"].get(symptom)
            cells.append("n/a" if s is None else f"{s['auprc']['point']:.3f} [{s['auprc']['lower']:.3f}, {s['auprc']['upper']:.3f}] ({s['base_rate']:.3f})")
        lines.append(f"| {symptom} | " + " | ".join(cells) + " |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
