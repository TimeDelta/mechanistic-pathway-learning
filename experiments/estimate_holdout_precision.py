"""Precision of a comparison scored on one grouped hold-out the size of a lockbox, from baseline-against-baseline differences.

The confirmatory comparison is scored once, on a lockbox of leakage groups that no pilot sees. How wide its interval will be
depends on the lockbox's size, the label counts and how much two models' scores co-vary, so it is estimated here before the
lockbox is drawn, from baselines only (the baselines are not the hypothesis, so this spends nothing). Each test fold of a
k-fold grouped split is one hold-out of share 1/k scored by models fitted without it, which is what a lockbox is; for every
fold and every pair of baselines, the paired bootstrap over the fold's perturbations gives the 95 percent half-width of the
macro AUPRC difference (symptoms with five or more positives in the fold) and of the micro AUPRC difference. The median over
folds is the expected half-width at that share; the macro reading also reports how many symptoms reach five positives.

Usage:
  python experiments/estimate_holdout_precision.py --baseline-dir runs/full/baselines_v2_disease_cluster_and_targets --num-folds 5
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    paired_bootstrap_macro_difference, paired_bootstrap_micro_difference, per_symptom_auprc, scorable_symptom)

MINIMUM_POSITIVES_TO_SCORE = 5


def half_width(reading: dict) -> float:
    return (reading["upper"] - reading["lower"]) / 2.0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v1_full_v2.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--baseline-split", default="grouped")
    parser.add_argument("--num-folds", type=int, default=5, help="the fold count the baselines were run with; the hold-out share is 1 / num_folds")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=1000)
    parser.add_argument("--json-output", type=Path, default=None)
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    if json.loads((arguments.baseline_dir / "perturbation_ids.json").read_text()) != data.perturbation_ids:
        raise ValueError(f"{arguments.baseline_dir}: baselines scored another perturbation set")
    folds = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_row = np.array([folds[p] for p in data.perturbation_ids])
    predictions = {path.stem[len(f"predictions_{arguments.baseline_split}_"):]: np.load(path)
                   for path in sorted(arguments.baseline_dir.glob(f"predictions_{arguments.baseline_split}_*.npy"))}
    readings = []
    for fold in range(arguments.num_folds):
        rows = fold_of_row == fold
        outcomes, mask = data.outcomes[rows], data.label_mask[rows]
        symptoms_scored = sum(scorable_symptom(outcomes, s, MINIMUM_POSITIVES_TO_SCORE, mask) for s in range(outcomes.shape[1]))
        for first, second in itertools.combinations(sorted(predictions), 2):
            macro = paired_bootstrap_macro_difference(predictions[first][rows], predictions[second][rows], outcomes, per_symptom_auprc, arguments.num_bootstrap, mask=mask)
            micro = paired_bootstrap_micro_difference(predictions[first][rows], predictions[second][rows], outcomes, arguments.num_bootstrap, mask=mask)
            readings.append({"fold": fold, "perturbations": int(rows.sum()), "kept_positive_pairs": int((outcomes * mask).sum()), "symptoms_in_macro": int(symptoms_scored),
                             "a": first, "b": second, "macro_difference": macro["difference"], "macro_half_width": half_width(macro),
                             "micro_difference": micro["difference"], "micro_half_width": half_width(micro)})
            print(json.dumps(readings[-1]), flush=True)
    summary = {"hold_out_share": 1.0 / arguments.num_folds, "num_folds": arguments.num_folds, "group_by": arguments.group_by,
               "median_macro_half_width": float(np.median([r["macro_half_width"] for r in readings])),
               "max_macro_half_width": float(np.max([r["macro_half_width"] for r in readings])),
               "median_micro_half_width": float(np.median([r["micro_half_width"] for r in readings])),
               "max_micro_half_width": float(np.max([r["micro_half_width"] for r in readings])),
               "symptoms_in_macro_per_fold": sorted({r["fold"]: r["symptoms_in_macro"] for r in readings}.items()), "readings": readings}
    print(json.dumps({key: value for key, value in summary.items() if key != "readings"}, indent=1))
    output = arguments.json_output or arguments.baseline_dir / "holdout_precision.json"
    output.write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
