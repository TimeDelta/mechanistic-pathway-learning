"""Phase 2: baselines B0 (popularity, degree-scaled popularity) and B1 (random walk with restart)
under the grouped perturbation-wise split and the pathway-wise split (design sections 5.6 and 6.1),
with the label-permutation negative control (section 6.3).

Scoring. Out-of-fold predictions are pooled and scored per symptom (AUPRC, AUROC with bootstrap
intervals over perturbations) plus mean reciprocal rank and hits-at-3 per perturbation. Pooling
has an artifact for constant-per-fold predictors: a fold with a low training base rate has a
high test base rate, so pooled AUROC of the popularity baseline falls below 0.5. Per-fold
metrics (mean and standard deviation over folds of the macro AUPRC and AUROC) are reported
alongside the pooled numbers; the pooled numbers stay because they are what the bootstrap
intervals and the symptom-ranking metrics are computed on.

Pathway-wise splits. Two definitions of a pathway are held out in turn: the curated modules of
design section 3.2 (docs/curated_pathway_modules.csv) and the reconstruction subsystems of
Human-GEM (every gene whose primary subsystem is S). A held-out set is every perturbation that
writes onto a gene of the pathway (the gene itself, or a drug targeting it); pathways with fewer
than --min-holdout-positives positive pairs are skipped. Predictions are pooled over the held-out
sets. Held-out sets are small and symptom profiles inside a pathway are homogeneous, so per-symptom
AUPRC inside one held-out set is not meaningful: per-split macro metrics are only computed for sets
of at least --min-fold-size-for-macro perturbations, per-split MRR and hits-at-3 are always reported,
and the pooled metrics are compared with the same split run on permuted labels.

Negative controls. (1) Outcome rows are permuted among perturbations within degree strata (hub
structure kept, biology destroyed) and the whole cross-validation is rerun on the permuted
labels; a model whose score survives the permutation is reading degree, not biology. (2) The
graph is rewired with degree-preserving double-edge swaps within each relation type and the
random walk is rerun on it under the grouped split (--rewiring-swaps-per-edge; 0 skips it); the
popularity baselines do not use the graph and are unchanged.

Time split (monogenic, design section 6.1). Pairs are dated by the earliest OMIM biocuration date
behind them (evidence_date; Orphanet-only pairs are undated). Training positives are the pairs dated
on or before --time-split-cutoff; the test question is which of the remaining pairs became positive
afterwards, so the scored pairs are every (perturbation, symptom) pair that is neither a training
positive nor an undated positive, labelled 1 when dated after the cutoff and 0 otherwise. Per-symptom
AUPRC and AUROC are computed over those pairs, with the base rate of new positives shown, and the same
evaluation is run with the new-positive labels permuted among the scored pairs of each symptom.

Writes runs/baselines/results.json and docs/phase2_baselines.md.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.negative_controls import degree_preserving_rewiring, permute_symptom_labels_within_degree_strata
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    assign_grouped_folds,
    perturbations_anchored_in_module,
    primary_subsystem_by_gene_node,
    read_curated_modules,
    subsystem_holdout_masks,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval,
    hits_at_k,
    macro_auprc_by_degree_bin,
    mean_reciprocal_rank,
    per_symptom_auprc,
    per_symptom_auroc,
)
from sklearn.metrics import average_precision_score, roc_auc_score
from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline
from mechanistic_pathway_learning.models.baselines.random_walk_with_restart_baseline import (
    RandomWalkWithRestartBaseline,
    build_normalized_adjacency,
)

BASELINE_NAMES = ("popularity", "degree_popularity", "random_walk_with_restart")
MINIMUM_POSITIVES_TO_SCORE = 5


def fit_and_predict(data, outcomes: np.ndarray, train: np.ndarray, test: np.ndarray, model_name: str, restart_probability: float, normalized_adjacency) -> np.ndarray:
    degrees = data.perturbation_degrees
    if model_name in ("popularity", "degree_popularity"):
        model = PopularityBaseline(scale_by_degree=model_name == "degree_popularity").fit(outcomes[train])
        return model.predict(int(test.sum()), degrees[test])
    if model_name == "random_walk_with_restart":
        model = RandomWalkWithRestartBaseline(restart_probability=restart_probability).fit(
            normalized_adjacency, [data.perturbation_seeds[i] for i in np.where(train)[0]], outcomes[train]
        )
        return model.predict([data.perturbation_seeds[i] for i in np.where(test)[0]])
    raise ValueError(model_name)


def macro_scores(predictions: np.ndarray, outcomes: np.ndarray) -> tuple[float, float]:
    """Macro AUPRC and AUROC over symptoms with enough positives and at least one negative."""
    auprcs, aurocs = [], []
    for symptom_index in range(outcomes.shape[1]):
        positives = outcomes[:, symptom_index].sum()
        if positives < MINIMUM_POSITIVES_TO_SCORE or positives == outcomes.shape[0]:
            continue
        auprcs.append(per_symptom_auprc(predictions, outcomes, symptom_index))
        aurocs.append(per_symptom_auroc(predictions, outcomes, symptom_index))
    return (float(np.mean(auprcs)) if auprcs else float("nan"), float(np.mean(aurocs)) if aurocs else float("nan"))


def run_split(data, outcomes: np.ndarray, test_masks: list[np.ndarray], model_name: str, restart_probability: float, normalized_adjacency,
              min_fold_size_for_macro: int = 20, split_labels: list[str] | None = None) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """Fit on the complement of each test mask, predict the mask; return pooled predictions, the scored-row mask and per-fold scores."""
    predictions = np.zeros_like(outcomes)
    scored = np.zeros(outcomes.shape[0], dtype=bool)
    per_fold = []
    for fold_index, test in enumerate(test_masks):
        train = ~test
        if test.sum() == 0 or train.sum() == 0:
            continue
        predictions[test] = fit_and_predict(data, outcomes, train, test, model_name, restart_probability, normalized_adjacency)
        scored |= test
        macro_auprc, macro_auroc = macro_scores(predictions[test], outcomes[test]) if test.sum() >= min_fold_size_for_macro else (float("nan"), float("nan"))
        per_fold.append({"fold": fold_index, "label": split_labels[fold_index] if split_labels else str(fold_index), "num_test": int(test.sum()), "num_positive_pairs": int(outcomes[test].sum()),
                         "macro_auprc": macro_auprc, "macro_auroc": macro_auroc,
                         "mean_reciprocal_rank": mean_reciprocal_rank(predictions[test], outcomes[test]), "hits_at_3": hits_at_k(predictions[test], outcomes[test], 3)})
    return predictions, scored, per_fold


def score(data, predictions: np.ndarray, outcomes: np.ndarray, rows: np.ndarray, per_fold: list[dict], num_bootstrap: int) -> dict:
    degrees = data.perturbation_degrees[rows]
    predictions, outcomes = predictions[rows], outcomes[rows]
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        positives = outcomes[:, symptom_index].sum()
        if positives < MINIMUM_POSITIVES_TO_SCORE or positives == outcomes.shape[0]:
            continue
        auprc = bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap)
        auroc = bootstrap_interval(lambda p, y: per_symptom_auroc(p, y, symptom_index), predictions, outcomes, num_bootstrap=num_bootstrap)
        per_symptom[symptom] = {"positives": int(positives), "auprc": auprc.__dict__, "auroc": auroc.__dict__, "base_rate": float(outcomes[:, symptom_index].mean())}
    fold_auprcs = [entry["macro_auprc"] for entry in per_fold if not np.isnan(entry["macro_auprc"])]
    fold_aurocs = [entry["macro_auroc"] for entry in per_fold if not np.isnan(entry["macro_auroc"])]
    return {
        "num_scored_perturbations": int(rows.sum()),
        "per_symptom": per_symptom,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "per_fold_macro_auprc_mean": float(np.mean(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auprc_sd": float(np.std(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auroc_mean": float(np.mean(fold_aurocs)) if fold_aurocs else float("nan"),
        "per_fold_macro_auroc_sd": float(np.std(fold_aurocs)) if fold_aurocs else float("nan"),
        "per_fold_mrr_mean": float(np.mean([entry["mean_reciprocal_rank"] for entry in per_fold])) if per_fold else float("nan"),
        "per_fold_hits_at_3_mean": float(np.mean([entry["hits_at_3"] for entry in per_fold])) if per_fold else float("nan"),
        "per_fold": per_fold,
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, outcomes),
        "hits_at_3": hits_at_k(predictions, outcomes, 3),
        "macro_auprc_by_degree_bin": macro_auprc_by_degree_bin(predictions, outcomes, degrees),
    }


def pathway_wise_test_masks(data, curated_modules_path: Path, min_holdout_positives: int) -> tuple[list[np.ndarray], list[str]]:
    masks, module_ids = [], []
    for module_id, gene_symbols in read_curated_modules(curated_modules_path).items():
        module_nodes = {data.node_index[f"GENE:{symbol}"] for symbol in gene_symbols if f"GENE:{symbol}" in data.node_index}
        mask = np.array(perturbations_anchored_in_module(data.perturbation_seeds, module_nodes))
        if mask.sum() == 0 or data.outcomes[mask].sum() < min_holdout_positives:
            continue
        masks.append(mask)
        module_ids.append(module_id)
    return masks, module_ids


def subsystem_test_masks(data, min_holdout_positives: int) -> tuple[list[np.ndarray], list[str]]:
    if data.node_subsystem is None:
        return [], []
    primary = primary_subsystem_by_gene_node(data.node_subsystem, data.edge_source, data.edge_target, data.edge_relation, data.relation_types.index("catalyzed_by"))
    masks_by_subsystem = subsystem_holdout_masks(data.perturbation_seeds, primary, data.outcomes.sum(axis=1), min_holdout_positives)
    return [np.array(mask) for mask in masks_by_subsystem.values()], list(masks_by_subsystem)


def holdout_section(title: str, description: str, entries: dict, control_entries: dict, labels: list[str], masks: list[np.ndarray], outcomes: np.ndarray) -> list[str]:
    header = ["| model | macro AUPRC (pooled) | macro AUROC (pooled) | MRR (pooled) | hits@3 (pooled) | MRR (mean over hold-outs) | hits@3 (mean over hold-outs) | scored perturbations |", "|---|---|---|---|---|---|---|---|"]

    def row(name: str, entry: dict) -> str:
        return (f"| {name} | {entry['macro_auprc']:.3f} | {entry['macro_auroc']:.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} | "
                f"{entry['per_fold_mrr_mean']:.3f} | {entry['per_fold_hits_at_3_mean']:.3f} | {entry['num_scored_perturbations']} |")

    lines = [f"## {title}", "", description, "",
             "Held out (perturbations, positive pairs): " + "; ".join(f"{label} ({int(mask.sum())}, {int(outcomes[mask].sum())})" for label, mask in zip(labels, masks)) + ".", ""]
    lines += header + [row(name, entry) for name, entry in entries.items()] + [""]
    if control_entries:
        lines += ["Same split with labels permuted within degree strata:", ""] + header + [row(name, entry) for name, entry in control_entries.items()] + [""]
    lines += ["| hold-out | perturbations | positive pairs | " + " | ".join(f"{name} MRR / hits@3" for name in entries) + " |", "|---|---|---|" + "---|" * len(entries)]
    for position, label in enumerate(labels):
        cells = []
        for entry in entries.values():
            fold_entry = next((f for f in entry["per_fold"] if f.get("label") == label), None)
            cells.append("n/a" if fold_entry is None else f"{fold_entry['mean_reciprocal_rank']:.2f} / {fold_entry['hits_at_3']:.2f}")
        lines.append(f"| {label} | {int(masks[position].sum())} | {int(outcomes[masks[position]].sum())} | " + " | ".join(cells) + " |")
    return lines + [""]


def time_split_evaluation(data, cutoff: date, model_name: str, restart_probability: float, normalized_adjacency, seed: int) -> dict | None:
    """Train on pairs dated on or before the cutoff; score the pairs that could still become positive."""
    if data.evidence_dates is None or (data.evidence_dates > 0).sum() == 0:
        return None
    cutoff_ordinal = cutoff.toordinal()
    dated = data.evidence_dates > 0
    training_outcomes = ((data.outcomes > 0) & dated & (data.evidence_dates <= cutoff_ordinal)).astype(float)
    new_positive = ((data.outcomes > 0) & dated & (data.evidence_dates > cutoff_ordinal)).astype(float)
    undated_positive = (data.outcomes > 0) & ~dated
    scored_pairs = ~(training_outcomes > 0) & ~undated_positive
    all_rows = np.ones(len(data.perturbation_ids), dtype=bool)
    predictions = fit_and_predict(data, training_outcomes, all_rows, all_rows, model_name, restart_probability, normalized_adjacency)
    generator = np.random.default_rng(seed)

    def per_symptom_table(labels: np.ndarray) -> dict:
        table = {}
        for symptom_index, symptom in enumerate(data.symptoms):
            mask = scored_pairs[:, symptom_index]
            positives = labels[mask, symptom_index].sum()
            if positives < MINIMUM_POSITIVES_TO_SCORE or positives == mask.sum():
                continue
            table[symptom] = {"scored_pairs": int(mask.sum()), "new_positives": int(positives), "base_rate": float(positives / mask.sum()),
                              "auprc": float(average_precision_score(labels[mask, symptom_index], predictions[mask, symptom_index])),
                              "auroc": float(roc_auc_score(labels[mask, symptom_index], predictions[mask, symptom_index]))}
        return table

    permuted = new_positive.copy()
    for symptom_index in range(permuted.shape[1]):
        rows = np.where(scored_pairs[:, symptom_index])[0]
        permuted[rows, symptom_index] = new_positive[generator.permutation(rows), symptom_index]
    observed_table, permuted_table = per_symptom_table(new_positive), per_symptom_table(permuted)
    return {
        "cutoff": cutoff.isoformat(), "training_positive_pairs": int(training_outcomes.sum()), "new_positive_pairs": int(new_positive.sum()),
        "undated_positive_pairs": int(undated_positive.sum()), "scored_pairs": int(scored_pairs.sum()),
        "per_symptom": observed_table, "per_symptom_permuted": permuted_table,
        "macro_auprc": float(np.mean([e["auprc"] for e in observed_table.values()])) if observed_table else float("nan"),
        "macro_auroc": float(np.mean([e["auroc"] for e in observed_table.values()])) if observed_table else float("nan"),
        "macro_auprc_permuted": float(np.mean([e["auprc"] for e in permuted_table.values()])) if permuted_table else float("nan"),
        "macro_auroc_permuted": float(np.mean([e["auroc"] for e in permuted_table.values()])) if permuted_table else float("nan"),
    }


def summary_row(name: str, entry: dict) -> str:
    return (f"| {name} | {entry['macro_auprc']:.3f} | {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f} | {entry['macro_auroc']:.3f} | "
            f"{entry['per_fold_macro_auroc_mean']:.3f} ± {entry['per_fold_macro_auroc_sd']:.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} | {entry['num_scored_perturbations']} |")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--restart-probability", type=float, default=0.3)
    parser.add_argument("--group-by", choices=["gene", "disease_cluster"], default="gene")
    parser.add_argument("--min-holdout-positives", type=int, default=10)
    parser.add_argument("--min-fold-size-for-macro", type=int, default=20)
    parser.add_argument("--skip-permutation-control", action="store_true")
    parser.add_argument("--rewiring-swaps-per-edge", type=int, default=2, help="degree-preserving rewiring control for the random walk; 0 skips it")
    parser.add_argument("--time-split-cutoff", type=date.fromisoformat, default=date(2015, 12, 31), help="monogenic time split: pairs dated on or before this day train")
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/baselines"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase2_baselines.md"))
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
    fold_of_perturbation = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
    grouped_masks = [fold_of_perturbation == fold for fold in range(arguments.num_folds)]
    pathway_masks, pathway_module_ids = pathway_wise_test_masks(data, arguments.curated_modules, arguments.min_holdout_positives)
    subsystem_masks, subsystem_labels = subsystem_test_masks(data, arguments.min_holdout_positives)
    normalized_adjacency = build_normalized_adjacency(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0])
    permuted_outcomes = permute_symptom_labels_within_degree_strata(data.outcomes, data.perturbation_degrees, random_seed=arguments.seed)

    num_genes = sum(t == "gene" for t in data.perturbation_types)
    num_drugs = sum(t == "drug" for t in data.perturbation_types)
    results = {
        "num_perturbations": len(data.perturbation_ids), "num_genes": num_genes, "num_drugs": num_drugs, "num_symptoms": len(data.symptoms), "symptoms": data.symptoms,
        "group_by": arguments.group_by, "num_groups": len(set(data.group_ids)), "num_folds": arguments.num_folds,
        "pathway_wise_modules": pathway_module_ids, "pathway_wise_holdout_sizes": [int(mask.sum()) for mask in pathway_masks],
        "pathway_wise_positives": [int(data.outcomes[mask].sum()) for mask in pathway_masks],
        "subsystem_wise_subsystems": subsystem_labels, "subsystem_wise_holdout_sizes": [int(mask.sum()) for mask in subsystem_masks],
        "subsystem_wise_positives": [int(data.outcomes[mask].sum()) for mask in subsystem_masks],
        "splits": {"grouped": {}, "pathway_wise": {}, "subsystem_wise": {}, "grouped_label_permutation": {}, "pathway_wise_label_permutation": {}, "subsystem_wise_label_permutation": {}, "grouped_rewired_graph": {}},
    }
    if arguments.rewiring_swaps_per_edge > 0:
        rewired = degree_preserving_rewiring(np.stack([data.edge_source, data.edge_target]), data.edge_relation, num_swaps_per_edge=arguments.rewiring_swaps_per_edge, random_seed=arguments.seed)
        rewired_adjacency = build_normalized_adjacency(len(data.node_ids), rewired[0], rewired[1], np.where(data.is_currency)[0])
        predictions, rows, per_fold = run_split(data, data.outcomes, grouped_masks, "random_walk_with_restart", arguments.restart_probability, rewired_adjacency, arguments.min_fold_size_for_macro)
        results["splits"]["grouped_rewired_graph"]["random_walk_with_restart"] = score(data, predictions, data.outcomes, rows, per_fold, arguments.num_bootstrap)
        entry = results["splits"]["grouped_rewired_graph"]["random_walk_with_restart"]
        print(f"{'grouped_rewired_graph':26s} {'random_walk_with_restart':26s} macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  macro AUROC {entry['macro_auroc']:.3f}")
    split_plan = [("grouped", grouped_masks, [str(fold) for fold in range(arguments.num_folds)]), ("pathway_wise", pathway_masks, pathway_module_ids), ("subsystem_wise", subsystem_masks, subsystem_labels)]
    for model_name in BASELINE_NAMES:
        for split_name, masks, labels in split_plan:
            if not masks:
                continue
            predictions, rows, per_fold = run_split(data, data.outcomes, masks, model_name, arguments.restart_probability, normalized_adjacency, arguments.min_fold_size_for_macro, labels)
            results["splits"][split_name][model_name] = score(data, predictions, data.outcomes, rows, per_fold, arguments.num_bootstrap)
            if not arguments.skip_permutation_control:
                predictions, rows, per_fold = run_split(data, permuted_outcomes, masks, model_name, arguments.restart_probability, normalized_adjacency, arguments.min_fold_size_for_macro, labels)
                results["splits"][f"{split_name}_label_permutation"][model_name] = score(data, predictions, permuted_outcomes, rows, per_fold, arguments.num_bootstrap)
        for split_name, entries in results["splits"].items():
            if model_name in entries:
                entry = entries[model_name]
                print(f"{split_name:26s} {model_name:26s} macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  "
                      f"macro AUROC {entry['macro_auroc']:.3f} (per fold {entry['per_fold_macro_auroc_mean']:.3f})  MRR {entry['mean_reciprocal_rank']:.3f}  hits@3 {entry['hits_at_3']:.3f}")
    results["time_split"] = {}
    for model_name in BASELINE_NAMES:
        entry = time_split_evaluation(data, arguments.time_split_cutoff, model_name, arguments.restart_probability, normalized_adjacency, arguments.seed)
        if entry is not None:
            results["time_split"][model_name] = entry
            print(f"{'time_split':26s} {model_name:26s} macro AUPRC {entry['macro_auprc']:.3f} (permuted {entry['macro_auprc_permuted']:.3f})  macro AUROC {entry['macro_auroc']:.3f} (permuted {entry['macro_auroc_permuted']:.3f})  new positives {entry['new_positive_pairs']} of {entry['scored_pairs']} scored pairs")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "results.json").write_text(json.dumps(results, indent=1))

    lines = ["# Phase 2 baselines (generated by experiments/run_baselines.py)", "",
             f"{len(data.perturbation_ids)} perturbations ({num_genes} genes, {num_drugs} drugs), {len(data.symptoms)} symptoms, relation induces. "
             f"Grouped split: {arguments.num_folds}-fold cross-validation with leakage groups by {arguments.group_by} ({len(set(data.group_ids))} groups); out-of-fold predictions pooled. "
             "Unobserved pairs count as negatives (positive-unlabelled convention). Pooled metrics are computed on the pooled out-of-fold predictions; per-fold values are the mean and standard deviation over folds of the macro metric, which is immune to the base-rate artifact that pulls pooled AUROC of a constant-per-fold predictor below 0.5.", ""]
    header = ["| model | macro AUPRC (pooled) | macro AUPRC (per fold) | macro AUROC (pooled) | macro AUROC (per fold) | MRR | hits@3 | scored perturbations |", "|---|---|---|---|---|---|---|---|"]
    lines += ["## Grouped perturbation-wise split", ""] + header + [summary_row(name, entry) for name, entry in results["splits"]["grouped"].items()] + [""]
    degree_bins = list(next(iter(results["splits"]["grouped"].values()))["macro_auprc_by_degree_bin"])
    lines += ["Macro AUPRC by perturbation degree tercile (pooled out-of-fold predictions; bins are [low degree, high degree] with the number of perturbations):", "",
              "| model | " + " | ".join(degree_bins) + " |", "|---|" + "---|" * len(degree_bins)]
    for name, entry in results["splits"]["grouped"].items():
        lines.append(f"| {name} | " + " | ".join(f"{entry['macro_auprc_by_degree_bin'][b]:.3f}" for b in degree_bins) + " |")
    lines.append("")
    if results["splits"]["pathway_wise"]:
        lines += holdout_section("Pathway-wise split, curated modules (each module of design section 3.2 held out in turn)",
                                 "Per-symptom AUPRC inside one held-out module is not meaningful (sets of 4 to 8 genes with homogeneous symptom profiles), so only pooled per-symptom metrics and per-hold-out ranking metrics are shown.",
                                 results["splits"]["pathway_wise"], results["splits"]["pathway_wise_label_permutation"], pathway_module_ids, pathway_masks, data.outcomes)
    if results["splits"]["subsystem_wise"]:
        lines += holdout_section("Pathway-wise split, Human-GEM subsystems (every gene whose primary subsystem is the held-out one)",
                                 "Subsystems are the reconstruction's own pathway partition; they cover many more annotated genes than the curated modules and are the candidate definition of the pre-registered pathway-wise split.",
                                 results["splits"]["subsystem_wise"], results["splits"]["subsystem_wise_label_permutation"], subsystem_labels, subsystem_masks, data.outcomes)
    if results["splits"]["grouped_label_permutation"]:
        lines += ["## Negative control: labels permuted within degree strata (grouped split)", ""] + header + [summary_row(name, entry) for name, entry in results["splits"]["grouped_label_permutation"].items()] + [""]
    if results["splits"]["grouped_rewired_graph"]:
        lines += [f"## Negative control: degree-preserving rewiring of the graph ({arguments.rewiring_swaps_per_edge} swaps per edge; grouped split, random walk only)", ""] + header + [summary_row(name, entry) for name, entry in results["splits"]["grouped_rewired_graph"].items()] + [""]
    lines += ["## Per-symptom AUPRC, grouped split", "", "Point estimate with 95 percent bootstrap interval over perturbations; base rate in parentheses.", "",
              "| symptom | " + " | ".join(results["splits"]["grouped"]) + " |", "|---|" + "---|" * len(results["splits"]["grouped"])]
    for symptom in data.symptoms:
        cells = []
        for entry in results["splits"]["grouped"].values():
            s = entry["per_symptom"].get(symptom)
            cells.append("n/a" if s is None else f"{s['auprc']['point']:.3f} [{s['auprc']['lower']:.3f}, {s['auprc']['upper']:.3f}] ({s['base_rate']:.3f})")
        lines.append(f"| {symptom} | " + " | ".join(cells) + " |")
    if results["time_split"]:
        first = next(iter(results["time_split"].values()))
        lines += ["", f"## Time split, monogenic pairs by OMIM biocuration date (cutoff {first['cutoff']})", "",
                  f"Training positives: {first['training_positive_pairs']} pairs dated on or before the cutoff; new positives: {first['new_positive_pairs']} pairs dated after it; undated (Orphanet-only) positives excluded: {first['undated_positive_pairs']}; scored pairs: {first['scored_pairs']} (every pair that is neither a training positive nor an undated positive). "
                  "The biocuration date is when the HPO team recorded the annotation, an upper bound on publication, so this split measures prediction of later-curated links. Permuted: new-positive labels permuted among the scored pairs of each symptom.", "",
                  "| model | macro AUPRC | macro AUPRC (permuted) | macro AUROC | macro AUROC (permuted) |", "|---|---|---|---|---|"]
        for name, entry in results["time_split"].items():
            lines.append(f"| {name} | {entry['macro_auprc']:.3f} | {entry['macro_auprc_permuted']:.3f} | {entry['macro_auroc']:.3f} | {entry['macro_auroc_permuted']:.3f} |")
        lines += ["", "| symptom | scored pairs | new positives (base rate) | " + " | ".join(f"{name} AUPRC / AUROC" for name in results["time_split"]) + " |", "|---|---|---|" + "---|" * len(results["time_split"])]
        for symptom in data.symptoms:
            cells = []
            for entry in results["time_split"].values():
                s = entry["per_symptom"].get(symptom)
                cells.append("n/a" if s is None else f"{s['auprc']:.3f} / {s['auroc']:.3f}")
            s = first["per_symptom"].get(symptom)
            if s is not None:
                lines.append(f"| {symptom} | {s['scored_pairs']} | {s['new_positives']} ({s['base_rate']:.3f}) | " + " | ".join(cells) + " |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
