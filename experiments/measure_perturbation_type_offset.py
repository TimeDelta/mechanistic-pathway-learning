"""How much a perturbation-type offset alone earns on the confirmatory readings, measured on development data only.

Drugs and genes differ in base rate (over labelled development pairs, drugs 0.159 and genes 0.063), and the lockbox
holds a larger share of drugs than the development set (64 of 317 against 78 of 1,222). A predictor that knows only
whether a perturbation is a drug can then rank drugs above genes on every symptom, and none of the four confirmatory
baselines is fitted per perturbation type. This script draws pseudo-lockboxes from the development perturbations with
the lockbox's drug share, fits popularity and degree_popularity (the trainer's PopularityBaseline) and two graph-free
type-aware predictors on the rest, and scores the four H1 readings the way experiments/score_confirmatory.py does:
macro AUPRC over the symptoms with five or more kept positives in the drawn rows, micro AUPRC over the symptoms with five
or more kept positives in the development set, and both after ranking every score inside degree strata (degree_strata
over all perturbations; degrees are not outcomes). Lockbox membership is used only to leave lockbox rows out; no lockbox
outcome is read.

  type_popularity         per-symptom base rate among training perturbations of the same type
  type_degree_popularity  type_popularity times log(1 + perturbation degree), as degree_popularity scales popularity

With --pilot-runs, it also splits each pilot's test-fold macro AUROC by perturbation type, which shows how much of it
comes from ranking drugs against genes.

  OMP_NUM_THREADS=1 python experiments/measure_perturbation_type_offset.py \\
      --pilot-runs runs/full/confirmatory_*_disease_cluster_and_targets_development
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from score_confirmatory import MINIMUM_POSITIVES_TO_SCORE, macro_auprc, micro_auprc

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    degree_strata,
    per_symptom_auroc,
    rank_normalise_within_groups,
    scorable_symptom,
)
from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline

READINGS = ("macro", "micro", "macro_within_degree_strata", "micro_within_degree_strata")
TYPE_AWARE = ("type_popularity", "type_degree_popularity")


def base_rates(outcomes: np.ndarray, mask: np.ndarray, rows: np.ndarray) -> np.ndarray:
    return (outcomes[rows] * mask[rows]).sum(axis=0) / np.maximum(mask[rows].sum(axis=0), 1)


def readings(predictions: np.ndarray, outcomes: np.ndarray, mask: np.ndarray, strata: np.ndarray, macro_columns: list[int], micro_columns: list[int]) -> dict[str, float]:
    ranked = rank_normalise_within_groups(predictions, strata)
    return {"macro": macro_auprc(predictions, outcomes, mask, macro_columns), "micro": micro_auprc(predictions, outcomes, mask, micro_columns),
            "macro_within_degree_strata": macro_auprc(ranked, outcomes, mask, macro_columns),
            "micro_within_degree_strata": micro_auprc(ranked, outcomes, mask, micro_columns)}


def pseudo_lockbox_draws(data, strata: np.ndarray, num_draws: int, drugs_per_draw: int, genes_per_draw: int, seed: int) -> dict[str, dict[str, np.ndarray]]:
    kinds = np.asarray(data.perturbation_types)
    outcomes = (data.outcomes > 0.5).astype(float)
    mask = np.asarray(data.label_mask, dtype=bool) if data.label_mask is not None else np.ones_like(outcomes, dtype=bool)
    degrees = data.perturbation_degrees
    micro_columns = [column for column in range(outcomes.shape[1]) if (outcomes[:, column] * mask[:, column]).sum() >= MINIMUM_POSITIVES_TO_SCORE]
    drugs, genes = np.flatnonzero(kinds == "drug"), np.flatnonzero(kinds == "gene")
    generator = np.random.default_rng(seed)
    scores: dict[str, dict[str, list[float]]] = {}
    for _ in range(num_draws):
        test = np.concatenate([generator.choice(drugs, drugs_per_draw, replace=False), generator.choice(genes, genes_per_draw, replace=False)])
        train = np.setdiff1d(np.arange(len(kinds)), test)
        macro_columns = [column for column in range(outcomes.shape[1]) if (outcomes[test, column] * mask[test, column]).sum() >= MINIMUM_POSITIVES_TO_SCORE]
        typed = np.stack([base_rates(outcomes, mask, train[kinds[train] == kinds[row]]) for row in test])
        predictions = {
            "popularity": PopularityBaseline().fit(outcomes[train], training_label_mask=mask[train]).predict(len(test), degrees[test]),
            "degree_popularity": PopularityBaseline(scale_by_degree=True).fit(outcomes[train], training_label_mask=mask[train]).predict(len(test), degrees[test]),
            "type_popularity": typed,
            "type_degree_popularity": typed * np.log1p(degrees[test])[:, None],
        }
        for name, prediction in predictions.items():
            for reading, value in readings(prediction, outcomes[test], mask[test], strata[test], macro_columns, micro_columns).items():
                scores.setdefault(name, {}).setdefault(reading, []).append(value)
    return {name: {reading: np.array(values) for reading, values in by_reading.items()} for name, by_reading in scores.items()}


def pilot_rows(data, pilot_runs: list[Path]) -> list[dict]:
    """Test-fold macro AUROC of each pilot over all test perturbations and inside each perturbation type."""
    kinds = np.asarray(data.perturbation_types)
    outcomes = (data.outcomes > 0.5).astype(float)
    mask = np.asarray(data.label_mask, dtype=bool) if data.label_mask is not None else np.ones_like(outcomes, dtype=bool)
    index_of = {perturbation: index for index, perturbation in enumerate(data.perturbation_ids)}
    rows = []
    for run in pilot_runs:
        for split in sorted(path for path in run.iterdir() if (path / "results.json").exists() and (path / "test_predictions.npy").exists()):
            results = json.loads((split / "results.json").read_text())
            test = np.array([index_of[perturbation] for perturbation in results["test_perturbation_ids"]])
            predictions = np.load(split / "test_predictions.npy")
            row = {"run": run.name, "split": split.name, "test": len(test), "test_drugs": int((kinds[test] == "drug").sum()),
                   "training_drugs": int((kinds == "drug").sum() - (kinds[test] == "drug").sum())}
            for label, selected in (("all", np.ones(len(test), dtype=bool)), ("genes", kinds[test] == "gene"), ("drugs", kinds[test] == "drug")):
                y, m, p = outcomes[test][selected], mask[test][selected], predictions[selected]
                values = [per_symptom_auroc(p, y, column, m) for column in range(y.shape[1]) if scorable_symptom(y, column, 1, m)]
                values = [value for value in values if not np.isnan(value)]
                row[f"macro_auroc_{label}"] = float(np.mean(values)) if values else float("nan")
                row[f"mean_prediction_{label}"] = float(p[m].mean()) if m.any() else float("nan")
                row[f"base_rate_{label}"] = float(y[m].mean()) if m.any() else float("nan")
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v1_full_v2.parquet"))
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v1.json"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--num-draws", type=int, default=200)
    parser.add_argument("--drugs-per-draw", type=int, default=39, help="39 drugs and 154 genes give the lockbox's drug share (64 of 317) from half the development drugs")
    parser.add_argument("--genes-per-draw", type=int, default=154)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--pilot-runs", type=Path, nargs="*", default=[])
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/perturbation_type_offset.md"))
    parser.add_argument("--json-output", type=Path, default=Path("runs/perturbation_type_offset.json"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by, label_selection=arguments.label_selection)
    strata_of_all = degree_strata(data.perturbation_degrees)
    in_lockbox = read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir)
    data = restrict_to_perturbations(data, ~in_lockbox)
    strata = strata_of_all[~in_lockbox]
    kinds = np.asarray(data.perturbation_types)
    outcomes = (data.outcomes > 0.5).astype(float)
    mask = np.asarray(data.label_mask, dtype=bool) if data.label_mask is not None else np.ones_like(outcomes, dtype=bool)
    rate_of = {kind: float((outcomes[kinds == kind] * mask[kinds == kind]).sum() / mask[kinds == kind].sum()) for kind in ("drug", "gene")}
    symptoms_drug_higher = int((base_rates(outcomes, mask, np.flatnonzero(kinds == "drug")) > base_rates(outcomes, mask, np.flatnonzero(kinds == "gene"))).sum())
    stratum_counts = {int(stratum): {"drugs": int(((strata == stratum) & (kinds == "drug")).sum()), "genes": int(((strata == stratum) & (kinds == "gene")).sum())}
                      for stratum in np.unique(strata)}

    draws = pseudo_lockbox_draws(data, strata, arguments.num_draws, arguments.drugs_per_draw, arguments.genes_per_draw, arguments.seed)
    best_type_blind = {reading: np.maximum(draws["popularity"][reading], draws["degree_popularity"][reading]) for reading in READINGS}
    differences = {name: {reading: draws[name][reading] - best_type_blind[reading] for reading in READINGS} for name in TYPE_AWARE}
    pilots = pilot_rows(data, arguments.pilot_runs)

    summary = {"development_perturbations": len(kinds), "development_drugs": int((kinds == "drug").sum()), "base_rate": rate_of,
               "symptoms_with_higher_drug_rate": symptoms_drug_higher, "num_symptoms": len(data.symptoms), "drugs_and_genes_per_degree_stratum": stratum_counts,
               "draws": arguments.num_draws, "drugs_per_draw": arguments.drugs_per_draw, "genes_per_draw": arguments.genes_per_draw,
               "mean_reading": {name: {reading: float(values.mean()) for reading, values in by_reading.items()} for name, by_reading in draws.items()},
               "difference_from_best_type_blind_popularity": {name: {reading: {"mean": float(values.mean()), "lower": float(np.percentile(values, 2.5)), "upper": float(np.percentile(values, 97.5))}
                                                                   for reading, values in by_reading.items()} for name, by_reading in differences.items()},
               "pilots": pilots}
    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps(summary, indent=1))

    lines = ["# Perturbation-type offset on the confirmatory readings (generated by experiments/measure_perturbation_type_offset.py)", "",
             f"Development set: {len(kinds)} perturbations, {summary['development_drugs']} of them drugs. Base rate over labelled pairs: drugs {rate_of['drug']:.3f}, "
             f"genes {rate_of['gene']:.3f}; the drug rate is higher on {symptoms_drug_higher} of {len(data.symptoms)} symptoms. No confirmatory baseline is fitted per "
             "perturbation type, so a predictor that knows only whether a perturbation is a drug can rank drugs above genes on those symptoms.", "",
             "Drugs and genes per degree stratum (degree_strata over all perturbations, development rows shown): drugs sit in the upper strata but share them "
             "with genes, so ranking inside strata does not remove the offset.", "",
             "| stratum | drugs | genes |", "|---|---|---|"]
    lines += [f"| {stratum} | {counts['drugs']} | {counts['genes']} |" for stratum, counts in stratum_counts.items()]
    lines += ["", f"## Pseudo-lockboxes ({arguments.num_draws} draws of {arguments.drugs_per_draw} drugs and {arguments.genes_per_draw} genes from the development set, "
              "the lockbox's drug share; baselines fitted on the remaining development perturbations)", "",
              "Mean reading per predictor:", "", "| predictor | " + " | ".join(READINGS) + " |", "|---|" + "---|" * len(READINGS)]
    lines += [f"| {name} | " + " | ".join(f"{draws[name][reading].mean():.3f}" for reading in READINGS) + " |" for name in draws]
    lines += ["", "Type-aware predictor minus the better of popularity and degree_popularity in the same draw (mean, and 2.5 and 97.5 percentiles over draws). "
              "The H1 floors are 0.041 (macro) and 0.028 (micro).", "",
              "| predictor | " + " | ".join(READINGS) + " |", "|---|" + "---|" * len(READINGS)]
    lines += [f"| {name} | " + " | ".join(f"{values[reading].mean():+.3f} [{np.percentile(values[reading], 2.5):+.3f}, {np.percentile(values[reading], 97.5):+.3f}]" for reading in READINGS) + " |"
              for name, values in differences.items()]
    lines += ["", "The random walk and TransE baselines are not drawn here (each needs the full graph per draw); whether either already carries the offset is open."]
    if pilots:
        lines += ["", "## Pilot test folds split by perturbation type", "",
                  "Macro AUROC over symptoms with a positive and a negative among the rows; mean prediction and base rate over labelled pairs.", "",
                  "| run | split | test | test drugs | training drugs | AUROC all | AUROC genes | AUROC drugs | mean prediction genes | mean prediction drugs | base rate genes | base rate drugs |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        lines += [f"| {row['run']} | {row['split']} | {row['test']} | {row['test_drugs']} | {row['training_drugs']} | {row['macro_auroc_all']:.3f} | {row['macro_auroc_genes']:.3f} | "
                  f"{row['macro_auroc_drugs']:.3f} | {row['mean_prediction_genes']:.3f} | {row['mean_prediction_drugs']:.3f} | {row['base_rate_genes']:.3f} | {row['base_rate_drugs']:.3f} |" for row in pilots]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
