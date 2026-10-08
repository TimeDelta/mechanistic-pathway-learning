"""Choose one head per encoder on the development set before any lockbox run, so that two models are tested instead of
four (the user's decision of 8 October 2026: fewer tested models, so a smaller chance that one passes by luck).

Rule (written before any of the runs it reads exists; docs/preregistration.md):
- Inputs: for each encoder, the noisy-OR and the sigmoid configuration run on the five grouped development folds with
  seed 0 (run_main_model_batch.py --lockbox <lockbox> --folds 0 1 2 3 4, into <configuration>_<group-by>_development),
  and the development baselines (run_baselines.py --lockbox <lockbox>, the same five grouped folds and seed).
- Readings per fold, as experiments/score_confirmatory.py reads the lockbox: macro AUPRC over the symptoms with five or
  more kept positives in the fold's test set, micro AUPRC over the symptoms with five or more kept positives in the
  development set, and both again after ranking inside degree strata. Each is the model minus the best baseline of that
  reading, the baseline with the highest five-fold mean score, chosen without reference to any model.
- A head's margin is the smallest of: mean macro difference minus the macro floor, mean micro difference minus the
  micro floor, and the two mean within-strata differences (means over the five folds). These are the conditions of H1,
  so the head chosen is the one whose weakest condition stands furthest above its threshold. Margins equal to 0.001
  choose the noisy-OR head, the design's mechanistic head.
- A missing fold for either head stops the script; no choice is made on part of the folds.
The choice is a development reading: it spends nothing of the lockbox. Writes configs/head_choice.json and
docs/head_choice.md.

Usage:
  python experiments/choose_heads.py --lockbox configs/lockbox_v2.json --pairs \\
      confirmatory_v2_message_passing_noisy_or:confirmatory_v2_message_passing_sigmoid \\
      confirmatory_v2_linear_response_noisy_or:confirmatory_v2_linear_response_sigmoid
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_confirmatory import BASELINES, MINIMUM_POSITIVES_TO_SCORE, macro_auprc, micro_auprc  # noqa: E402

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations  # noqa: E402
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds  # noqa: E402
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import degree_strata, rank_normalise_within_groups  # noqa: E402

READINGS = ("macro", "micro", "macro_within_degree_strata", "micro_within_degree_strata")
TIE_TOLERANCE = 0.001


def fold_readings(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray, strata: np.ndarray, micro_columns: list[int]) -> dict[str, float]:
    kept = outcomes * mask
    macro_columns = [column for column in range(outcomes.shape[1]) if kept[:, column].sum() >= MINIMUM_POSITIVES_TO_SCORE]
    ranked = rank_normalise_within_groups(predictions, strata)
    return {"macro": macro_auprc(predictions, outcomes, mask, macro_columns), "micro": micro_auprc(predictions, outcomes, mask, micro_columns),
            "macro_within_degree_strata": macro_auprc(ranked, outcomes, mask, macro_columns),
            "micro_within_degree_strata": micro_auprc(ranked, outcomes, mask, micro_columns)}


def head_margin(mean_differences: dict[str, float], macro_floor: float, micro_floor: float) -> float:
    """The smallest of H1's conditions, each as its distance above its threshold."""
    return min(mean_differences["macro"] - macro_floor, mean_differences["micro"] - micro_floor,
               mean_differences["macro_within_degree_strata"], mean_differences["micro_within_degree_strata"])


def choose(noisy_or_margin: float, sigmoid_margin: float) -> str:
    return "sigmoid" if sigmoid_margin > noisy_or_margin + TIE_TOLERANCE else "noisy_or"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v2.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--lockbox", type=Path, required=True)
    parser.add_argument("--pairs", nargs="+", required=True, help="noisy_or_configuration:sigmoid_configuration, one per encoder")
    parser.add_argument("--run-root", type=Path, default=Path("runs/full"))
    parser.add_argument("--baseline-dir", type=Path, required=True, help="run_baselines.py --lockbox <lockbox> output on the development set")
    parser.add_argument("--baselines", nargs="*", default=list(BASELINES))
    parser.add_argument("--folds", type=int, nargs="*", default=[0, 1, 2, 3, 4])
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--minimum-macro-difference", type=float, default=0.032, help="the scorer's macro floor (docs/preregistration.md)")
    parser.add_argument("--minimum-micro-difference", type=float, default=0.043, help="the scorer's micro floor")
    parser.add_argument("--json-output", type=Path, default=Path("configs/head_choice.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/head_choice.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    in_lockbox = read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir)
    strata_all = degree_strata(data.perturbation_degrees)
    development = restrict_to_perturbations(data, ~in_lockbox)
    strata = strata_all[~in_lockbox]
    mask = development.label_mask if development.label_mask is not None else np.ones_like(development.outcomes, dtype=bool)
    kept = development.outcomes * mask
    micro_columns = [column for column in range(kept.shape[1]) if kept[:, column].sum() >= MINIMUM_POSITIVES_TO_SCORE]
    row_of = {perturbation_id: row for row, perturbation_id in enumerate(development.perturbation_ids)}
    fold_of = assign_grouped_folds(development.perturbation_ids, development.group_ids, arguments.num_folds, arguments.seed)
    if json.loads((arguments.baseline_dir / "perturbation_ids.json").read_text()) != development.perturbation_ids:
        raise SystemExit(f"{arguments.baseline_dir} scored other perturbations than this development set")
    baseline_predictions = {name: np.load(arguments.baseline_dir / f"predictions_grouped_{name}.npy") for name in arguments.baselines}

    fold_rows = {fold: np.array([row_of[p] for p in development.perturbation_ids if fold_of[p] == fold]) for fold in arguments.folds}
    baseline_scores = {name: {fold: fold_readings(predictions[rows], development.outcomes[rows], mask[rows], strata[rows], micro_columns)
                              for fold, rows in fold_rows.items()} for name, predictions in baseline_predictions.items()}
    best_baseline = {reading: max(arguments.baselines, key=lambda name: np.mean([baseline_scores[name][fold][reading] for fold in arguments.folds])) for reading in READINGS}

    choices, details, missing = {}, {}, []
    for pair in arguments.pairs:
        noisy_or_configuration, sigmoid_configuration = pair.split(":")
        per_head = {}
        for head, configuration in (("noisy_or", noisy_or_configuration), ("sigmoid", sigmoid_configuration)):
            differences = {reading: [] for reading in READINGS}
            for fold, rows in fold_rows.items():
                split_directory = arguments.run_root / f"{configuration}_{arguments.group_by}_development" / f"fold{fold}_seed{arguments.seed}"
                if not (split_directory / "DONE").exists():
                    missing.append(str(split_directory))
                    continue
                results = json.loads((split_directory / "results.json").read_text())
                test_rows = np.array([row_of[p] for p in results["test_perturbation_ids"]])
                if not np.array_equal(np.sort(test_rows), np.sort(rows)):
                    raise SystemExit(f"{split_directory} tested other perturbations than development fold {fold}")
                predictions = np.empty_like(development.outcomes)
                predictions[test_rows] = np.load(split_directory / "test_predictions.npy")
                model = fold_readings(predictions[rows], development.outcomes[rows], mask[rows], strata[rows], micro_columns)
                for reading in READINGS:
                    differences[reading].append(model[reading] - baseline_scores[best_baseline[reading]][fold][reading])
            mean_differences = {reading: float(np.mean(values)) if values else float("nan") for reading, values in differences.items()}
            per_head[head] = {"configuration": configuration, "per_fold": differences, "mean": mean_differences,
                              "margin": head_margin(mean_differences, arguments.minimum_macro_difference, arguments.minimum_micro_difference)}
        details[pair] = per_head
        choices[pair] = per_head[choose(per_head["noisy_or"]["margin"], per_head["sigmoid"]["margin"])]["configuration"]
    if missing:
        raise SystemExit("no choice: development runs missing:\n" + "\n".join(missing))

    record = {"lockbox": str(arguments.lockbox), "folds": arguments.folds, "seed": arguments.seed, "best_baseline": best_baseline,
              "floors": {"macro": arguments.minimum_macro_difference, "micro": arguments.minimum_micro_difference}, "chosen": choices, "details": details}
    arguments.json_output.write_text(json.dumps(record, indent=1) + "\n")
    lines = ["# Head choice per encoder on the development set (generated by experiments/choose_heads.py)", "",
             f"Lockbox {arguments.lockbox}; development folds {arguments.folds}, seed {arguments.seed}. Best baseline per reading (five-fold mean, "
             "no model enters the choice): " + "; ".join(f"{reading} {name}" for reading, name in best_baseline.items()) + ".", "",
             "Each cell: five-fold mean of model minus best baseline. Margin: the smallest of macro minus its floor, micro minus its floor and the two "
             "within-strata differences; the head with the larger margin is tested on the lockbox (equal to 0.001: noisy-OR).", "",
             "| encoder pair | head | macro | micro | macro within strata | micro within strata | margin | chosen |", "|---|---|---|---|---|---|---|---|"]
    for pair, per_head in details.items():
        for head, entry in per_head.items():
            cells = " | ".join(f"{entry['mean'][reading]:+.3f}" for reading in READINGS)
            lines.append(f"| {pair} | {head} | {cells} | {entry['margin']:+.3f} | {'yes' if choices[pair] == entry['configuration'] else ''} |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
