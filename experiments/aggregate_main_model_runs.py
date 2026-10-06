"""Aggregate run directories of experiments/run_main_model.py into one comparison table (design sections 6 and 7).

A run directory (runs/<configuration>_<group_by>/) holds one split per subdirectory with results.json,
test_predictions.npy and the test perturbation ids. For grouped folds the out-of-fold predictions are
pooled and scored exactly as the baseline harness scores B0 and B1 (pooled per-symptom AUPRC and AUROC
with bootstrap intervals, per-fold mean and standard deviation of the macro metrics, MRR, hits-at-3,
expected calibration error). For hold-out splits (curated modules or reconstruction subsystems) the
held-out sets are pooled the same way. Module supports and link matrices of noisy-OR runs are
summarized (support sizes, number of symptoms per module above the link threshold).

Usage:
  python experiments/aggregate_main_model_runs.py --run-dirs runs/b6_default_disease_cluster runs/b3_sigmoid_disease_cluster \
      --baseline-results runs/baselines_disease_cluster/results.json --markdown-output docs/phase3_main_model.md
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    expected_calibration_error,
    hits_at_k,
    macro_auprc_by_degree_bin,
    mean_reciprocal_rank,
    paired_bootstrap_macro_difference,
    per_symptom_auprc,
    per_symptom_auroc,
)

MINIMUM_POSITIVES_TO_SCORE = 5


def macro_scores(predictions: np.ndarray, outcomes: np.ndarray) -> tuple[float, float]:
    auprcs, aurocs = [], []
    for symptom_index in range(outcomes.shape[1]):
        positives = outcomes[:, symptom_index].sum()
        if positives < MINIMUM_POSITIVES_TO_SCORE or positives == outcomes.shape[0]:
            continue
        auprcs.append(per_symptom_auprc(predictions, outcomes, symptom_index))
        aurocs.append(per_symptom_auroc(predictions, outcomes, symptom_index))
    return (float(np.mean(auprcs)) if auprcs else float("nan"), float(np.mean(aurocs)) if aurocs else float("nan"))


def aggregate_run_directory(run_directory: Path, data, num_bootstrap: int) -> dict | None:
    split_directories = sorted(path for path in run_directory.iterdir() if path.is_dir() and (path / "DONE").exists() and (path / "results.json").exists())
    split_directories = [path for path in split_directories if json.loads((path / "results.json").read_text()).get("time_split") is None]
    if not split_directories:
        return None
    position_of = {perturbation_id: index for index, perturbation_id in enumerate(data.perturbation_ids)}
    pooled_predictions = np.zeros_like(data.outcomes)
    scored = np.zeros(len(data.perturbation_ids), dtype=bool)
    per_split = []
    support_sizes, expected_support_sizes, symptoms_per_module, epochs = [], [], [], []
    code_commits: set[str] = set()
    for split_directory in split_directories:
        results = json.loads((split_directory / "results.json").read_text())
        predictions = np.load(split_directory / "test_predictions.npy")
        rows = np.array([position_of[p] for p in results["test_perturbation_ids"] if p in position_of])
        if len(rows) != len(results["test_perturbation_ids"]):
            raise ValueError(f"{split_directory}: test perturbations not all present in the current evidence table")
        pooled_predictions[rows] = predictions
        scored[rows] = True
        macro_auprc, macro_auroc = macro_scores(predictions, data.outcomes[rows])
        per_split.append({"split": results["split"], "num_test": len(rows), "macro_auprc": macro_auprc, "macro_auroc": macro_auroc,
                          "mean_reciprocal_rank": mean_reciprocal_rank(predictions, data.outcomes[rows]), "hits_at_3": hits_at_k(predictions, data.outcomes[rows], 3),
                          "epochs_completed": results.get("epochs_completed"), "best_epoch": results.get("best_epoch")})
        epochs.append(results.get("epochs_completed", 0))
        code_commits.update(entry["commit"][:7] for entry in results.get("code_provenance", []) if entry.get("commit"))
        if "module_symptom_links" in results:
            links = np.array(results["module_symptom_links"])
            support_sizes.append(results["module_support_sizes"])
            expected_support_sizes.append(results.get("module_expected_support_sizes", []))
            symptoms_per_module.append((links > 0.5).sum(axis=1).tolist())
    outcomes, predictions = data.outcomes[scored], pooled_predictions[scored]
    np.save(run_directory / "pooled_predictions.npy", pooled_predictions)
    np.save(run_directory / "pooled_scored_rows.npy", scored)
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        positives = outcomes[:, symptom_index].sum()
        if positives < MINIMUM_POSITIVES_TO_SCORE or positives == outcomes.shape[0]:
            continue
        per_symptom[symptom] = {
            "positives": int(positives), "base_rate": float(outcomes[:, symptom_index].mean()),
            "auprc": bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap).__dict__,
            "auroc": bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap).__dict__,
        }
    fold_auprcs = [entry["macro_auprc"] for entry in per_split if not np.isnan(entry["macro_auprc"])]
    fold_aurocs = [entry["macro_auroc"] for entry in per_split if not np.isnan(entry["macro_auroc"])]
    first_results = json.loads((split_directories[0] / "results.json").read_text())
    return {
        "run_directory": str(run_directory), "num_splits": len(split_directories), "num_scored_perturbations": int(scored.sum()),
        "code_commits": sorted(code_commits),
        "arguments": {key: first_results["arguments"].get(key) for key in ("head", "field", "pooling", "num_modules", "group_by", "description_length_coefficient", "learning_rate", "node_state_dim", "num_layers", "permute_labels")},
        "per_symptom": per_symptom,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "per_fold_macro_auprc_mean": float(np.mean(fold_auprcs)) if fold_auprcs else float("nan"), "per_fold_macro_auprc_sd": float(np.std(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auroc_mean": float(np.mean(fold_aurocs)) if fold_aurocs else float("nan"), "per_fold_macro_auroc_sd": float(np.std(fold_aurocs)) if fold_aurocs else float("nan"),
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, outcomes), "hits_at_3": hits_at_k(predictions, outcomes, 3),
        "expected_calibration_error": expected_calibration_error(predictions, outcomes),
        "macro_auprc_by_degree_bin": macro_auprc_by_degree_bin(predictions, outcomes, data.perturbation_degrees[scored]),
        "per_split": per_split, "mean_epochs": float(np.mean(epochs)) if epochs else float("nan"),
        "module_support_sizes": support_sizes, "module_expected_support_sizes": expected_support_sizes, "symptoms_per_module_above_half": symptoms_per_module,
    }


def collect_time_split_runs(run_directories: list[Path]) -> list[dict]:
    rows = []
    for run_directory in run_directories:
        for split_directory in sorted(path for path in run_directory.iterdir() if path.is_dir() and (path / "DONE").exists()):
            results = json.loads((split_directory / "results.json").read_text())
            if results.get("time_split"):
                rows.append({"run": run_directory.name, "split": results["split"], **results["time_split"]})
    return rows


def paired_comparisons(aggregated: dict, run_directories: list[Path], baseline_directory: Path | None, baseline_split: str, data, num_bootstrap: int, include_baseline_pairs: bool = True) -> list[dict]:
    """Paired bootstrap of the pooled macro AUPRC and AUROC difference between every run and every baseline (and between runs) on the rows both scored."""
    predictions_by_name: dict[str, np.ndarray] = {}
    rows_by_name: dict[str, np.ndarray] = {}
    for run_directory in run_directories:
        if run_directory.name in aggregated and (run_directory / "pooled_predictions.npy").exists():
            predictions_by_name[run_directory.name] = np.load(run_directory / "pooled_predictions.npy")
            rows_by_name[run_directory.name] = np.load(run_directory / "pooled_scored_rows.npy")
    if baseline_directory is not None and (baseline_directory / "perturbation_ids.json").exists():
        baseline_ids = json.loads((baseline_directory / "perturbation_ids.json").read_text())
        if baseline_ids == data.perturbation_ids:
            rows_path = baseline_directory / f"scored_rows_{baseline_split}.npy"
            for prediction_path in sorted(baseline_directory.glob(f"predictions_{baseline_split}_*.npy")):
                name = prediction_path.stem[len(f"predictions_{baseline_split}_"):]
                predictions_by_name[name] = np.load(prediction_path)
                rows_by_name[name] = np.load(rows_path) if rows_path.exists() else np.ones(len(data.perturbation_ids), dtype=bool)
    comparisons = []
    names = list(predictions_by_name)
    for position, first in enumerate(names):
        for second in names[position + 1:]:
            if first not in aggregated and second not in aggregated and not include_baseline_pairs:
                continue
            if first not in aggregated and second in aggregated:
                first, second = second, first  # a run-versus-baseline pair is reported once, with the run as A
            rows = rows_by_name[first] & rows_by_name[second]
            auprc = paired_bootstrap_macro_difference(predictions_by_name[first][rows], predictions_by_name[second][rows], data.outcomes[rows], per_symptom_auprc, num_bootstrap)
            auroc = paired_bootstrap_macro_difference(predictions_by_name[first][rows], predictions_by_name[second][rows], data.outcomes[rows], per_symptom_auroc, num_bootstrap)
            comparisons.append({"a": first, "b": second, "rows": int(rows.sum()), "macro_auprc": auprc, "macro_auroc": auroc})
    return comparisons


def summary_row(name: str, entry: dict) -> str:
    return (f"| {name} | {entry['macro_auprc']:.3f} | {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f} | {entry['macro_auroc']:.3f} | "
            f"{entry['per_fold_macro_auroc_mean']:.3f} ± {entry['per_fold_macro_auroc_sd']:.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} | "
            f"{entry.get('expected_calibration_error', float('nan')):.3f} | {entry['num_scored_perturbations']} |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--group-by", choices=["gene", "disease_cluster"], default="disease_cluster")
    parser.add_argument("--baseline-results", type=Path, default=None, help="results.json of run_baselines.py for the same split, to put B0 and B1 in the same table")
    parser.add_argument("--baseline-split", default="grouped", help="key under 'splits' in the baseline results to show")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--title", default="Phase 3: proposed model and sigmoid-head baseline")
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase3_main_model.md"))
    parser.add_argument("--json-output", type=Path, default=Path("runs/phase3_aggregate.json"))
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by)
    aggregated = {}
    for run_directory in arguments.run_dirs:
        entry = aggregate_run_directory(run_directory, data, arguments.num_bootstrap)
        if entry is None:
            print(f"{run_directory}: no finished splits")
            continue
        aggregated[run_directory.name] = entry
        print(f"{run_directory.name:40s} splits {entry['num_splits']}  macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  "
              f"macro AUROC {entry['macro_auroc']:.3f}  MRR {entry['mean_reciprocal_rank']:.3f}  hits@3 {entry['hits_at_3']:.3f}  ECE {entry['expected_calibration_error']:.3f}")
    baseline_entries = {}
    baseline_directory = None
    if arguments.baseline_results and arguments.baseline_results.exists():
        baseline_entries = json.loads(arguments.baseline_results.read_text())["splits"].get(arguments.baseline_split, {})
        baseline_directory = arguments.baseline_results.parent
    comparisons = paired_comparisons(aggregated, arguments.run_dirs, baseline_directory, arguments.baseline_split, data, arguments.num_bootstrap)
    for comparison in comparisons:
        print(f"{comparison['a']:34s} vs {comparison['b']:28s} macro AUPRC diff {comparison['macro_auprc']['difference']:+.3f} [{comparison['macro_auprc']['lower']:+.3f}, {comparison['macro_auprc']['upper']:+.3f}]  "
              f"macro AUROC diff {comparison['macro_auroc']['difference']:+.3f} [{comparison['macro_auroc']['lower']:+.3f}, {comparison['macro_auroc']['upper']:+.3f}]")
    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps({"runs": aggregated, "paired_comparisons": comparisons}, indent=1))

    lines = [f"# {arguments.title} (generated by experiments/aggregate_main_model_runs.py)", "",
             f"{len(data.perturbation_ids)} perturbations, {len(data.symptoms)} symptoms, leakage groups by {arguments.group_by}. Out-of-split predictions pooled; per-fold values are mean ± standard deviation over splits. "
             "Baseline rows come from experiments/run_baselines.py on the same split. ECE: expected calibration error over all perturbation-symptom pairs.", "",
             "| model | macro AUPRC (pooled) | macro AUPRC (per fold) | macro AUROC (pooled) | macro AUROC (per fold) | MRR | hits@3 | ECE | scored perturbations |", "|---|---|---|---|---|---|---|---|---|"]
    for name, entry in baseline_entries.items():
        lines.append(summary_row(name, {**entry, "expected_calibration_error": float("nan")}))
    for name, entry in aggregated.items():
        lines.append(summary_row(name, entry))
    if aggregated:
        degree_bins = list(next(iter(aggregated.values()))["macro_auprc_by_degree_bin"])
        lines += ["", "Macro AUPRC by perturbation degree tercile (pooled predictions):", "", "| model | " + " | ".join(degree_bins) + " |", "|---|" + "---|" * len(degree_bins)]
        for name, entry in list(baseline_entries.items()) + list(aggregated.items()):
            bins = entry.get("macro_auprc_by_degree_bin", {})
            lines.append(f"| {name} | " + " | ".join(f"{bins[b]:.3f}" if b in bins else "n/a" for b in degree_bins) + " |")
    if comparisons:
        lines += ["", "## Paired bootstrap comparisons (design section 7)", "", "Difference in pooled macro AUPRC and macro AUROC between two models on the rows both scored; 95 percent percentile interval of the paired bootstrap over perturbations. The pre-registered primary endpoint asks for a difference of at least 0.05 with an interval excluding zero.", "",
                  "| A | B | rows | macro AUPRC A - B [95% CI] | resamples favoring A | macro AUROC A - B [95% CI] |", "|---|---|---|---|---|---|"]
        for comparison in comparisons:
            a, r = comparison["macro_auprc"], comparison["macro_auroc"]
            lines.append(f"| {comparison['a']} | {comparison['b']} | {comparison['rows']} | {a['difference']:+.3f} [{a['lower']:+.3f}, {a['upper']:+.3f}] | {a['fraction_resamples_favoring_a']:.2f} | {r['difference']:+.3f} [{r['lower']:+.3f}, {r['upper']:+.3f}] |")
    lines += ["", "## Configurations", "", "Commits are those the splits ran under (recorded from 6 October 2026; 'not recorded' marks earlier runs, which evaluated noisy-OR gates with the Louizos test-time estimator rather than the expected training gate).", "",
              "| run | head | field | pooling | modules | description-length coefficient | learning rate | state dim | layers | labels permuted | splits | mean epochs | commits |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, entry in aggregated.items():
        a = entry["arguments"]
        commits = ", ".join(entry.get("code_commits") or []) or "not recorded"
        lines.append(f"| {name} | {a.get('head')} | {a.get('field')} | {a.get('pooling')} | {a.get('num_modules')} | {a.get('description_length_coefficient')} | {a.get('learning_rate')} | {a.get('node_state_dim')} | {a.get('num_layers')} | {a.get('permute_labels')} | {entry['num_splits']} | {entry['mean_epochs']:.1f} | {commits} |")
    lines += ["", "## Per-symptom AUPRC (pooled; 95 percent bootstrap interval over perturbations; base rate in parentheses)", "",
              "| symptom | " + " | ".join(list(baseline_entries) + list(aggregated)) + " |", "|---|" + "---|" * (len(baseline_entries) + len(aggregated))]
    for symptom in data.symptoms:
        cells = []
        for entry in list(baseline_entries.values()) + list(aggregated.values()):
            s = entry["per_symptom"].get(symptom)
            cells.append("n/a" if s is None else f"{s['auprc']['point']:.3f} [{s['auprc']['lower']:.3f}, {s['auprc']['upper']:.3f}] ({s['base_rate']:.3f})")
        lines.append(f"| {symptom} | " + " | ".join(cells) + " |")
    for name, entry in aggregated.items():
        if entry["module_support_sizes"]:
            lines += ["", f"## Module structure, {name}", "", "Support size counts gates above 0.5 in the evaluation gate (see the commits note under Configurations); expected support is the hard-concrete expected number of nonzero gates. Symptoms per module counts links above 0.5.", "",
                      "| split | support sizes | expected support sizes | symptoms per module above 0.5 |", "|---|---|---|---|"]
            for split_entry, sizes, expected, per_module in zip(entry["per_split"], entry["module_support_sizes"], entry["module_expected_support_sizes"], entry["symptoms_per_module_above_half"]):
                lines.append(f"| {split_entry['split']} | {sizes} | {[round(x) for x in expected] if expected else 'n/a'} | {per_module} |")
    time_split_rows = collect_time_split_runs(arguments.run_dirs)
    if time_split_rows:
        lines += ["", "## Time split (monogenic pairs by OMIM biocuration date)", "", "Training positives are the pairs dated on or before the cutoff; scored pairs are those that could still become positive; permuted: new-positive labels permuted within each symptom's scored pairs.", "",
                  "| run | split | training positives | new positives | scored pairs | macro AUPRC | macro AUPRC (permuted) | macro AUROC | macro AUROC (permuted) |", "|---|---|---|---|---|---|---|---|---|"]
        for row in time_split_rows:
            lines.append(f"| {row['run']} | {row['split']} | {row['training_positive_pairs']} | {row['new_positive_pairs']} | {row['scored_pairs']} | {row['macro_auprc']:.3f} | {row['macro_auprc_permuted']:.3f} | {row['macro_auroc']:.3f} | {row['macro_auroc_permuted']:.3f} |")
    lines += ["", "## Per-split macro metrics", "", "| run | split | test size | macro AUPRC | macro AUROC | MRR | hits@3 | epochs | best epoch |", "|---|---|---|---|---|---|---|---|---|"]
    for name, entry in aggregated.items():
        for split_entry in entry["per_split"]:
            lines.append(f"| {name} | {split_entry['split']} | {split_entry['num_test']} | {split_entry['macro_auprc']:.3f} | {split_entry['macro_auroc']:.3f} | {split_entry['mean_reciprocal_rank']:.3f} | {split_entry['hits_at_3']:.3f} | {split_entry['epochs_completed']} | {split_entry['best_epoch']} |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
