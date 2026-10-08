"""Module health of noisy-OR runs: whether the pathway modules add to some perturbations only or add the same amount to
every perturbation, a second leak (docs/best_epoch_zero.md; rule in docs/module_health.md, written before the module-fix
runs it reads). Reads each fold's results.json, test_module_activations.npy, test_predictions.npy and module_support.npy;
trains nothing.

Per fold:
- responding_modules: modules whose activation over the test perturbations has a standard deviation above 0.05.
- floor_share: per symptom, the modules' part of P for each test perturbation, m = 1 - (1 - P) / (1 - leak); the 10th
  percentile of m over the perturbations divided by its mean; the median over the symptoms whose mean m is above 0.01.
  Near 1: the modules add about the same to every perturbation. Near 0: they add to some perturbations only.
- median_link_correlation: the median Pearson correlation between the link vectors (over symptoms) of the pairs of
  modules with a nonempty hard support.
- nonempty_modules: modules with at least one node of hard support (gate above 0.5); median_support_size over those.
- macro_auprc over the symptoms with five or more kept positives in the fold's test set, micro_auprc over the symptoms
  with five or more kept positives in the data set.

Usage:
  python experiments/module_health.py --run-root runs/module_fix \\
      --pairs b6_mechanistic_gate_time_scales_weighted_start:b6_mechanistic_gate_time_scales \\
              b6_linear_response_gate_time_scales_cofactors_weighted_start:b6_linear_response_gate_time_scales_cofactors \\
      --markdown-output docs/module_health_results.md
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_confirmatory import MINIMUM_POSITIVES_TO_SCORE, macro_auprc, micro_auprc

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

RESPONDING_STANDARD_DEVIATION = 0.05
FLOOR_PERCENTILE = 10
MINIMUM_MEAN_MODULE_PART = 0.01
READINGS = ("responding_modules", "floor_share", "median_link_correlation", "nonempty_modules", "median_support_size", "macro_auprc", "micro_auprc")


def floor_share(predictions: np.ndarray, leaks: np.ndarray) -> float:
    module_part = 1.0 - (1.0 - predictions) / (1.0 - leaks[None, :])
    mean_part = module_part.mean(axis=0)
    shares = [np.percentile(module_part[:, symptom], FLOOR_PERCENTILE) / mean_part[symptom]
              for symptom in range(module_part.shape[1]) if mean_part[symptom] > MINIMUM_MEAN_MODULE_PART]
    return float(np.median(shares)) if shares else float("nan")


def median_link_correlation(links: np.ndarray, nonempty: np.ndarray) -> float:
    rows = links[nonempty]
    if len(rows) < 2:
        return float("nan")
    correlation = np.corrcoef(rows)
    upper = correlation[np.triu_indices(len(rows), k=1)]
    return float(np.nanmedian(upper))


def fold_health(split_directory: Path, data, mask: np.ndarray, micro_columns: list[int]) -> dict[str, float]:
    results = json.loads((split_directory / "results.json").read_text())
    activations = np.load(split_directory / "test_module_activations.npy")
    predictions = np.load(split_directory / "test_predictions.npy")
    support = np.load(split_directory / "module_support.npy")
    links = np.asarray(results["module_symptom_links"], dtype=float)
    links = links[0] if links.ndim == 3 else links
    leaks = np.asarray(results["symptom_leaks"], dtype=float)
    leaks = leaks[0] if leaks.ndim == 2 else leaks
    support_sizes = (support > 0.5).sum(axis=1)
    nonempty = support_sizes > 0
    row_of = {perturbation_id: row for row, perturbation_id in enumerate(data.perturbation_ids)}
    rows = np.array([row_of[p] for p in results["test_perturbation_ids"]])
    outcomes, fold_mask = data.outcomes[rows], mask[rows]
    macro_columns = [c for c in range(outcomes.shape[1]) if (outcomes[:, c] * fold_mask[:, c]).sum() >= MINIMUM_POSITIVES_TO_SCORE]
    return {"responding_modules": float((activations.std(axis=0) > RESPONDING_STANDARD_DEVIATION).sum()),
            "floor_share": floor_share(predictions, leaks),
            "median_link_correlation": median_link_correlation(links, nonempty),
            "nonempty_modules": float(nonempty.sum()),
            "median_support_size": float(np.median(support_sizes[nonempty])) if nonempty.any() else 0.0,
            "macro_auprc": macro_auprc(predictions, outcomes, fold_mask, macro_columns),
            "micro_auprc": micro_auprc(predictions, outcomes, fold_mask, micro_columns)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--pairs", nargs="+", required=True, help="fix_configuration:twin_configuration")
    parser.add_argument("--group-by", default="disease_cluster")
    parser.add_argument("--folds", type=int, nargs="*", default=[0, 1, 2, 3, 4])
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    arguments = parser.parse_args()

    loaded = {}
    record, lines = {}, ["# Module health of the module-fix runs (generated by experiments/module_health.py)", "",
                         ("Rule and thresholds: docs/module_health.md. Each cell: five-fold mean (fix, twin) and the paired mean difference fix minus "
                          "twin with the number of folds where the fix is higher."), "",
                         "| pair | reading | fix | twin | fix minus twin | folds fix higher |", "|---|---|---|---|---|---|"]
    for pair in arguments.pairs:
        fix, twin = pair.split(":")
        per_configuration = {}
        for configuration in (fix, twin):
            readings = {reading: [] for reading in READINGS}
            for fold in arguments.folds:
                split_directory = arguments.run_root / f"{configuration}_{arguments.group_by}" / f"fold{fold}_seed{arguments.seed}"
                if not (split_directory / "DONE").exists():
                    raise SystemExit(f"{split_directory} is not done")
                run_arguments = json.loads((split_directory / "results.json").read_text())["arguments"]
                key = (run_arguments["graph_dir"], run_arguments["evidence_dir"], run_arguments.get("label_selection"))
                if key not in loaded:
                    data = load_experiment_data(Path(key[0]), Path(key[1]), group_by=arguments.group_by, label_selection=Path(key[2]) if key[2] else None)
                    mask = data.label_mask if data.label_mask is not None else np.ones_like(data.outcomes, dtype=bool)
                    kept = data.outcomes * mask
                    loaded[key] = (data, mask, [c for c in range(kept.shape[1]) if kept[:, c].sum() >= MINIMUM_POSITIVES_TO_SCORE])
                for reading, value in fold_health(split_directory, *loaded[key]).items():
                    readings[reading].append(value)
            per_configuration[configuration] = readings
        record[pair] = per_configuration
        for reading in READINGS:
            fix_values, twin_values = np.array(per_configuration[fix][reading]), np.array(per_configuration[twin][reading])
            differences = fix_values - twin_values
            lines.append(f"| {pair} | {reading} | {np.nanmean(fix_values):.3f} | {np.nanmean(twin_values):.3f} | {np.nanmean(differences):+.3f} | "
                         f"{int((differences > 0).sum())} of {len(differences)} |")
    if arguments.json_output:
        arguments.json_output.write_text(json.dumps(record, indent=1) + "\n")
    if arguments.markdown_output:
        arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
