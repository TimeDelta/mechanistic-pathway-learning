"""Phase 2: baselines B0 (popularity, degree-scaled popularity), B1 (random walk with restart) and,
with --with-kg-embedding, B2 (TransE knowledge-graph embedding over the graph plus training evidence
triples) and, with --with-type-popularity, popularity and degree_popularity fitted per perturbation type
(type_popularity, type_degree_popularity; docs/perturbation_type_offset.md) under the grouped perturbation-wise split and the pathway-wise splits (design sections 5.6
and 6.1), with the label-permutation negative control (section 6.3). B2 is trained once per split
(seconds to a minute on CPU) and, unless --kg-embedding-all-splits is given, only for the grouped
split, its permutation control and the time split.

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
than --min-holdout-positives positive pairs (kept pairs only, under a label selection) are skipped.
Perturbations outside a held-out set that share a leakage group with one inside it are left out of
that fit, since a shared disease cluster or drug target carries the held-out labels. Predictions are
pooled over the held-out sets. Held-out sets are small and symptom profiles inside a pathway are homogeneous, so per-symptom
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

Lockbox (docs/preregistration.md). --lockbox removes the lockbox perturbations of configs/lockbox_v1.json before
anything else, so every split above runs on the development set. --lockbox with --score-lockbox fits once on the
development set and scores the lockbox (split "lockbox", with its label-permutation and rewired-graph controls); the
pathway-wise splits and the time split are not run then. The labels are permuted within the lockbox and within the
development set apart, so no lockbox label reaches the fit.

Writes runs/baselines/results.json and docs/phase2_baselines.md.
"""
from __future__ import annotations

import argparse
import json
from datetime import date
from functools import partial
from pathlib import Path

import numpy as np

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.measurement_coverage_strata import measurement_coverage_strata, read_measurement_coverage
from mechanistic_pathway_learning.evaluation.negative_controls import degree_preserving_rewiring, degree_stratified_row_permutation, duplicate_edge_count, fast_degree_preserving_rewiring, rewire_walk_graph
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    assign_grouped_folds,
    perturbations_anchored_in_module,
    primary_subsystem_by_gene_node,
    read_curated_modules,
    subsystem_holdout_masks,
    training_mask_without_group_partners,
)
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    rank_normalise_within_groups,
    stratified_auroc,
    bootstrap_interval,
    hits_at_k,
    macro_auprc_by_degree_bin,
    macro_auprc_by_stratum,
    mean_reciprocal_rank,
    micro_auprc,
    per_symptom_auprc,
    per_symptom_auroc,
    scorable_symptom,
    scored_rows,
)
from sklearn.metrics import average_precision_score, roc_auc_score
from mechanistic_pathway_learning.evaluation.result_cache import array_digest, cached_result, source_digest, sparse_digest
from mechanistic_pathway_learning.models.baselines.knowledge_graph_embedding_baseline import KnowledgeGraphEmbeddingBaseline
from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline
from mechanistic_pathway_learning.models.baselines.random_walk_with_restart_baseline import (
    RandomWalkWithRestartBaseline,
    build_normalized_adjacency,
)

BASELINE_NAMES = ("popularity", "degree_popularity", "random_walk_with_restart")
CACHED_BASELINES = ("random_walk_with_restart", "knowledge_graph_embedding_transe")  # the slow fits; popularity takes milliseconds
BASELINE_SOURCES = [Path(__file__), *sorted((Path(__file__).resolve().parents[1] / "mechanistic_pathway_learning" / "models" / "baselines").glob("*.py"))]
_DIGEST_MEMO: dict[int, tuple[object, str]] = {}
KG_EMBEDDING_NAME = "knowledge_graph_embedding_transe"
TYPE_POPULARITY_NAMES = ("type_popularity", "type_degree_popularity")
KG_EMBEDDING_SETTINGS = {"model_name": "TransE", "embedding_dim": 64, "num_epochs": 30, "batch_size": 4096, "learning_rate": 0.01, "margin": 1.0, "negatives_per_positive": 4}
MINIMUM_POSITIVES_TO_SCORE = 5


def fit_and_predict(data, outcomes: np.ndarray, train: np.ndarray, test: np.ndarray, model_name: str, restart_probability: float, normalized_adjacency,
                    random_seed: int = 0, label_mask: np.ndarray | None = None) -> np.ndarray:
    """Fit one baseline on the training rows and predict the test rows; random_seed sets the TransE initialisation and
    negative sampling, and label_mask (a label selection) keeps set-aside pairs out of the popularity base rates."""
    degrees = data.perturbation_degrees_for_strata
    if model_name in ("popularity", "degree_popularity"):
        model = PopularityBaseline(scale_by_degree=model_name == "degree_popularity").fit(outcomes[train], training_label_mask=None if label_mask is None else label_mask[train])
        return model.predict(int(test.sum()), degrees[test])
    if model_name in TYPE_POPULARITY_NAMES:  # popularity or degree_popularity fitted on the training perturbations of the test perturbation's type
        kinds = np.asarray(data.perturbation_types)
        predictions = np.zeros((int(test.sum()), outcomes.shape[1]))
        for kind in np.unique(kinds[test]):
            training_rows = train & (kinds == kind)
            if not training_rows.any():  # no training perturbation of this type: fall back to every training perturbation
                training_rows = train
            model = PopularityBaseline(scale_by_degree=model_name == "type_degree_popularity").fit(outcomes[training_rows], training_label_mask=None if label_mask is None else label_mask[training_rows])
            test_rows_of_kind = kinds[test] == kind
            predictions[test_rows_of_kind] = model.predict(int(test_rows_of_kind.sum()), degrees[test][test_rows_of_kind])
        return predictions
    if model_name == "random_walk_with_restart":
        model = RandomWalkWithRestartBaseline(restart_probability=restart_probability).fit(
            normalized_adjacency, [data.perturbation_seeds[i] for i in np.where(train)[0]], outcomes[train]
        )
        return model.predict([data.perturbation_seeds[i] for i in np.where(test)[0]])
    if model_name == KG_EMBEDDING_NAME:
        model = KnowledgeGraphEmbeddingBaseline(**KG_EMBEDDING_SETTINGS, random_seed=random_seed).fit(
            len(data.node_ids), data.edge_source, data.edge_target, data.edge_relation, len(data.relation_types),
            [data.perturbation_seeds[i] for i in np.where(train)[0]], outcomes[train]
        )
        return model.predict([data.perturbation_seeds[i] for i in np.where(test)[0]])
    raise ValueError(model_name)


def macro_scores(predictions: np.ndarray, outcomes: np.ndarray, label_mask: np.ndarray | None = None) -> tuple[float, float]:
    """Macro AUPRC and AUROC over symptoms with enough positives and at least one negative among the labelled pairs."""
    auprcs, aurocs = [], []
    for symptom_index in range(outcomes.shape[1]):
        if not scorable_symptom(outcomes, symptom_index, MINIMUM_POSITIVES_TO_SCORE, label_mask):
            continue
        auprcs.append(per_symptom_auprc(predictions, outcomes, symptom_index, label_mask))
        aurocs.append(per_symptom_auroc(predictions, outcomes, symptom_index, label_mask))
    return (float(np.mean(auprcs)) if auprcs else float("nan"), float(np.mean(aurocs)) if aurocs else float("nan"))


def memoised_digest(value, compute) -> str:
    """Digest of an object that lives for the whole run (the graph, an adjacency), computed once; the reference kept in
    the memo stops its id from being reused."""
    if id(value) not in _DIGEST_MEMO:
        _DIGEST_MEMO[id(value)] = (value, compute(value))
    return _DIGEST_MEMO[id(value)][1]


def fold_cache_key(data, training_outcomes: np.ndarray, train: np.ndarray, test: np.ndarray, model_name: str, restart_probability: float,
                   normalized_adjacency, random_seed: int, label_mask: np.ndarray | None) -> dict:
    """Every input of one baseline fit: the rows, the training labels, the graph the fit reads and the code."""
    graph = memoised_digest(data, lambda d: array_digest(np.array([len(d.node_ids)]), d.edge_source, d.edge_target, d.edge_relation,
                                                         np.array([len(seeds) for seeds in d.perturbation_seeds]),
                                                         np.concatenate([np.asarray(seeds, dtype=np.int64) for seeds in d.perturbation_seeds]) if d.perturbation_seeds else None))
    return {"what": "baseline fold", "model": model_name, "restart_probability": restart_probability, "random_seed": random_seed,
            "kg_embedding_settings": KG_EMBEDDING_SETTINGS if model_name == KG_EMBEDDING_NAME else None,
            "rows": array_digest(train, test), "training_outcomes": array_digest(training_outcomes), "label_mask": array_digest(label_mask),
            "graph": graph, "adjacency": memoised_digest(normalized_adjacency, sparse_digest) if model_name == "random_walk_with_restart" else None,
            "source": memoised_digest(BASELINE_SOURCES, lambda paths: source_digest(*paths))}


def run_split(data, outcomes: np.ndarray, test_masks: list[np.ndarray], model_name: str, restart_probability: float, normalized_adjacency,
              min_fold_size_for_macro: int = 20, split_labels: list[str] | None = None, label_mask: np.ndarray | None = None,
              random_seed: int = 0, group_ids: list[str] | None = None, cache_directory: Path | None = None) -> tuple[np.ndarray, np.ndarray, list[dict], np.ndarray]:
    """Fit on the complement of each test mask, predict the mask; return pooled predictions, the scored-row mask and per-fold scores.

    With a label mask, pairs set aside are not training positives (the fit sees them as unlabelled) and are left out
    of every score. With group_ids (the pathway-wise hold-outs, which select perturbations by seed gene), perturbations
    sharing a leakage group with a held-out one are left out of that fit as well. With cache_directory, the slow fits
    (random walk, TransE) are stored per fold under a key of all their inputs (result_cache), so a restarted run
    reuses the folds it finished."""
    training_outcomes = outcomes if label_mask is None else outcomes * label_mask
    predictions = np.zeros_like(outcomes)
    scored = np.zeros(outcomes.shape[0], dtype=bool)
    fold_of_row = np.full(outcomes.shape[0], -1, dtype=int)
    per_fold = []
    for fold_index, test in enumerate(test_masks):
        train = ~test if group_ids is None else np.array(training_mask_without_group_partners(test, group_ids))
        if test.sum() == 0 or train.sum() == 0:
            continue
        fit = partial(fit_and_predict, data, training_outcomes, train, test, model_name, restart_probability, normalized_adjacency, random_seed, label_mask)
        if cache_directory is not None and model_name in CACHED_BASELINES:
            predictions[test] = cached_result(cache_directory, fold_cache_key(data, training_outcomes, train, test, model_name, restart_probability,
                                                                              normalized_adjacency, random_seed, label_mask), fit)
        else:
            predictions[test] = fit()
        scored |= test
        fold_of_row[test] = fold_index
        test_label_mask = None if label_mask is None else label_mask[test]
        macro_auprc, macro_auroc = macro_scores(predictions[test], outcomes[test], test_label_mask) if test.sum() >= min_fold_size_for_macro else (float("nan"), float("nan"))
        per_fold.append({"fold": fold_index, "label": split_labels[fold_index] if split_labels else str(fold_index), "num_test": int(test.sum()), "num_positive_pairs": int(training_outcomes[test].sum()),
                         "num_group_partners_left_out": int((~test & ~train).sum()),
                         "macro_auprc": macro_auprc, "macro_auroc": macro_auroc,
                         "micro_auprc": micro_auprc(predictions[test], outcomes[test], test_label_mask) if test.sum() >= min_fold_size_for_macro else float("nan"),
                         "mean_reciprocal_rank": mean_reciprocal_rank(predictions[test], outcomes[test], test_label_mask), "hits_at_3": hits_at_k(predictions[test], outcomes[test], 3, test_label_mask)})
    return predictions, scored, per_fold, fold_of_row


def score(data, predictions: np.ndarray, outcomes: np.ndarray, rows: np.ndarray, per_fold: list[dict], num_bootstrap: int, fold_of_row: np.ndarray | None = None,
          label_mask: np.ndarray | None = None, coverage=None) -> dict:
    degrees = data.perturbation_degrees_for_strata[rows]
    scored_perturbation_ids = [perturbation_id for perturbation_id, keep in zip(data.perturbation_ids, rows) if keep] if coverage is not None else []
    fold_of_scored_row = fold_of_row[rows] if fold_of_row is not None else None
    predictions, outcomes = predictions[rows], outcomes[rows]
    mask = None if label_mask is None else label_mask[rows]
    normalised = rank_normalise_within_groups(predictions, fold_of_scored_row) if fold_of_scored_row is not None else None
    per_symptom = {}
    for symptom_index, symptom in enumerate(data.symptoms):
        if not scorable_symptom(outcomes, symptom_index, MINIMUM_POSITIVES_TO_SCORE, mask):
            continue
        labelled = scored_rows(len(outcomes), symptom_index, mask)
        positives = outcomes[labelled, symptom_index].sum()
        auprc = bootstrap_interval(lambda p, y, m=None: per_symptom_auprc(p, y, symptom_index, m), predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask)
        auroc = bootstrap_interval(lambda p, y, m=None: per_symptom_auroc(p, y, symptom_index, m), predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask)
        per_symptom[symptom] = {"positives": int(positives), "auprc": auprc.__dict__, "auroc": auroc.__dict__, "base_rate": float(outcomes[labelled, symptom_index].mean())}
        if normalised is not None:
            per_symptom[symptom]["auprc_rank_normalised"] = per_symptom_auprc(normalised, outcomes, symptom_index, mask)
            per_symptom[symptom]["auroc_rank_normalised"] = per_symptom_auroc(normalised, outcomes, symptom_index, mask)
            per_symptom[symptom]["auroc_stratified"] = stratified_auroc(predictions, outcomes, fold_of_scored_row, symptom_index, mask)
    fold_auprcs = [entry["macro_auprc"] for entry in per_fold if not np.isnan(entry["macro_auprc"])]
    fold_aurocs = [entry["macro_auroc"] for entry in per_fold if not np.isnan(entry["macro_auroc"])]
    fold_micro_auprcs = [entry["micro_auprc"] for entry in per_fold if not np.isnan(entry.get("micro_auprc", float("nan")))]
    return {
        "num_scored_perturbations": int(rows.sum()),
        "micro_auprc": bootstrap_interval(micro_auprc, predictions, outcomes, num_bootstrap=num_bootstrap, mask=mask).__dict__,
        "micro_auprc_rank_normalised": micro_auprc(normalised, outcomes, mask) if normalised is not None else float("nan"),
        "per_fold_micro_auprc_mean": float(np.mean(fold_micro_auprcs)) if fold_micro_auprcs else float("nan"),
        "per_fold_micro_auprc_sd": float(np.std(fold_micro_auprcs)) if fold_micro_auprcs else float("nan"),
        "per_symptom": per_symptom,
        "macro_auprc": float(np.mean([entry["auprc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auroc": float(np.mean([entry["auroc"]["point"] for entry in per_symptom.values()])) if per_symptom else float("nan"),
        "macro_auprc_rank_normalised": float(np.mean([entry["auprc_rank_normalised"] for entry in per_symptom.values()])) if per_symptom and normalised is not None else float("nan"),
        "macro_auroc_rank_normalised": float(np.mean([entry["auroc_rank_normalised"] for entry in per_symptom.values()])) if per_symptom and normalised is not None else float("nan"),
        "macro_auroc_stratified": float(np.nanmean([entry["auroc_stratified"] for entry in per_symptom.values()])) if per_symptom and normalised is not None else float("nan"),
        "per_fold_macro_auprc_mean": float(np.mean(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auprc_sd": float(np.std(fold_auprcs)) if fold_auprcs else float("nan"),
        "per_fold_macro_auroc_mean": float(np.mean(fold_aurocs)) if fold_aurocs else float("nan"),
        "per_fold_macro_auroc_sd": float(np.std(fold_aurocs)) if fold_aurocs else float("nan"),
        "per_fold_mrr_mean": float(np.mean([entry["mean_reciprocal_rank"] for entry in per_fold])) if per_fold else float("nan"),
        "per_fold_hits_at_3_mean": float(np.mean([entry["hits_at_3"] for entry in per_fold])) if per_fold else float("nan"),
        "per_fold": per_fold,
        "mean_reciprocal_rank": mean_reciprocal_rank(predictions, outcomes, mask),
        "hits_at_3": hits_at_k(predictions, outcomes, 3, mask),
        "macro_auprc_by_degree_bin": macro_auprc_by_degree_bin(predictions, outcomes, degrees, mask=mask),
        "macro_auprc_by_measurement_coverage": macro_auprc_by_measurement_coverage(predictions, outcomes, scored_perturbation_ids, coverage, mask),
    }


def macro_auprc_by_measurement_coverage(predictions: np.ndarray, outcomes: np.ndarray, scored_perturbation_ids: list[str], coverage, mask) -> dict[str, float]:
    """Macro AUPRC inside strata of how much binding measurement each drug has, or an empty result without the table."""
    if coverage is None:
        return {}
    stratum_of_row, labels = measurement_coverage_strata(scored_perturbation_ids, coverage)
    scores = macro_auprc_by_stratum(predictions, outcomes, stratum_of_row, mask=mask)
    return {label: scores.get(label, float("nan")) for label in labels}


def kept_positives(data) -> np.ndarray:
    """Positive pairs that are scored: all of them, or with a label selection only the kept ones."""
    return data.outcomes if data.label_mask is None else data.outcomes * data.label_mask


def pathway_wise_test_masks(data, curated_modules_path: Path, min_holdout_positives: int) -> tuple[list[np.ndarray], list[str]]:
    masks, module_ids = [], []
    for module_id, gene_symbols in read_curated_modules(curated_modules_path).items():
        module_nodes = {data.node_index[f"GENE:{symbol}"] for symbol in gene_symbols if f"GENE:{symbol}" in data.node_index}
        mask = np.array(perturbations_anchored_in_module(data.perturbation_seeds, module_nodes))
        if mask.sum() == 0 or kept_positives(data)[mask].sum() < min_holdout_positives:
            continue
        masks.append(mask)
        module_ids.append(module_id)
    return masks, module_ids


def subsystem_test_masks(data, min_holdout_positives: int) -> tuple[list[np.ndarray], list[str]]:
    if data.node_subsystem is None:
        return [], []
    primary = primary_subsystem_by_gene_node(data.node_subsystem, data.edge_source, data.edge_target, data.edge_relation, data.relation_types.index("catalyzed_by"))
    masks_by_subsystem = subsystem_holdout_masks(data.perturbation_seeds, primary, kept_positives(data).sum(axis=1), min_holdout_positives)
    return [np.array(mask) for mask in masks_by_subsystem.values()], list(masks_by_subsystem)


def holdout_section(title: str, description: str, entries: dict, control_entries: dict, labels: list[str], masks: list[np.ndarray], outcomes: np.ndarray) -> list[str]:
    header = ["| model | macro AUPRC (pooled, raw) | macro AUPRC (pooled, rank within hold-out) | macro AUROC (pooled, raw) | macro AUROC (rank within hold-out) | macro AUROC (stratified) | MRR (pooled) | hits@3 (pooled) | MRR (mean over hold-outs) | hits@3 (mean over hold-outs) | scored perturbations |", "|---|---|---|---|---|---|---|---|---|---|---|"]

    def row(name: str, entry: dict) -> str:
        return (f"| {name} | {entry['macro_auprc']:.3f} | {entry.get('macro_auprc_rank_normalised', float('nan')):.3f} | {entry['macro_auroc']:.3f} | {entry.get('macro_auroc_rank_normalised', float('nan')):.3f} | {entry.get('macro_auroc_stratified', float('nan')):.3f} | {entry['mean_reciprocal_rank']:.3f} | {entry['hits_at_3']:.3f} | "
                f"{entry['per_fold_mrr_mean']:.3f} | {entry['per_fold_hits_at_3_mean']:.3f} | {entry['num_scored_perturbations']} |")

    lines = [f"## {title}", "", description, " Raw pooled scores across hold-outs with different base rates carry the artifact of section 6.2 (review v0.4, finding 1); the rank-within-hold-out columns replace each score by its tie-averaged rank over (n + 1) inside its hold-out before pooling, and the stratified AUROC forms positive-negative pairs inside hold-outs only. The pathway-wise endpoint is defined on the rank-within-hold-out macro AUPRC. Perturbations sharing a leakage group with a held-out one are left out of that hold-out's fit.", "",
             "Held out (perturbations, scored positive pairs): " + "; ".join(f"{label} ({int(mask.sum())}, {int(outcomes[mask].sum())})" for label, mask in zip(labels, masks)) + ".", ""]
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


def time_split_evaluation(data, cutoff: date, model_name: str, restart_probability: float, normalized_adjacency, seed: int, require_publication_date: bool = False) -> dict | None:
    """Train on pairs dated on or before the cutoff; score the pairs that could still become positive.

    With require_publication_date, a pair counts as a new positive only when its date is a cited publication date;
    a pair whose only post-cutoff date is an HPO curation date (an upper bound on when the observation existed,
    review v0.4, finding 7) is excluded from scoring like an undated pair. Training positives are unchanged: a
    curation date on or before the cutoff proves the evidence existed by then.
    """
    if data.evidence_dates is None or (data.evidence_dates > 0).sum() == 0:
        return None
    if require_publication_date and data.evidence_date_is_publication is None:
        return None
    cutoff_ordinal = cutoff.toordinal()
    dated = data.evidence_dates > 0
    training_outcomes = ((data.outcomes > 0) & dated & (data.evidence_dates <= cutoff_ordinal)).astype(float)
    after_cutoff = (data.outcomes > 0) & dated & (data.evidence_dates > cutoff_ordinal)
    curation_dated_after_cutoff = after_cutoff & ~data.evidence_date_is_publication if require_publication_date else np.zeros_like(after_cutoff)
    new_positive = (after_cutoff & ~curation_dated_after_cutoff).astype(float)
    undated_positive = ((data.outcomes > 0) & ~dated) | curation_dated_after_cutoff
    scored_pairs = ~(training_outcomes > 0) & ~undated_positive
    all_rows = np.ones(len(data.perturbation_ids), dtype=bool)
    predictions = fit_and_predict(data, training_outcomes, all_rows, all_rows, model_name, restart_probability, normalized_adjacency, seed)
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
        "curation_dated_new_positives_excluded": int(curation_dated_after_cutoff.sum()), "require_publication_date": require_publication_date,
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
    parser.add_argument("--seed", type=int, default=0, help="seeds the fold assignment, the label permutation, the graph rewiring and the TransE initialisation and negative sampling")
    parser.add_argument("--measurement-coverage", type=Path, default=None,
                        help="coverage table of experiments/scope_off_target_binding.py --coverage-table; adds macro AUPRC inside strata of how "
                             "much binding measurement each drug has, the reading that says whether a result rides on the overstudied drugs")
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--restart-probability", type=float, default=0.3)
    parser.add_argument("--group-by", choices=["gene", "disease_cluster", "disease_cluster_and_targets"], default="gene")
    parser.add_argument("--label-grades", nargs="*", default=["A", "B"], help="evidence grades that count as positive labels; pass A B C to keep grade C rows as the version 0.3 ablation did")
    parser.add_argument("--label-selection", type=Path, default=None,
                        help="parquet of (perturbation_id, symptom, keep) from experiments/build_label_selection.py; positive pairs with keep False are masked out of fitting and scoring")
    parser.add_argument("--min-holdout-positives", type=int, default=10)
    parser.add_argument("--min-fold-size-for-macro", type=int, default=20)
    parser.add_argument("--skip-permutation-control", action="store_true")
    parser.add_argument("--rewiring-swaps-per-edge", type=int, default=50,
                        help="degree-preserving rewiring control for the random walk; 0 skips it. The default was 2 until 7 October 2026, "
                             "which left the slice graph under-mixed (docs/graph_content_null_results.md); 50 matches "
                             "experiments/run_rewiring_null_distribution.py. Documents generated before then say 2 in their rewiring heading")
    parser.add_argument("--rewiring-method", choices=["pair_sampling", "integer_draws", "walk_graph"], default="pair_sampling",
                        help="pair_sampling: degree_preserving_rewiring, O(edges) per swap, behind every rewiring result before 7 October 2026; "
                             "integer_draws: fast_degree_preserving_rewiring, the same swap rule at O(1) per swap (the full graph); "
                             "walk_graph: the walk's own undirected simple graph rewired so every node keeps its number of neighbours (the other two move apart "
                             "two stored edges joining one pair and give the walk a denser graph; negative_controls.rewire_walk_graph)")
    parser.add_argument("--lockbox", type=Path, default=None, help="lockbox file from experiments/draw_lockbox.py; its perturbations are removed before every split")
    parser.add_argument("--score-lockbox", action="store_true", help="with --lockbox: fit on the development set and score the lockbox once")
    parser.add_argument("--time-split-cutoff", type=date.fromisoformat, default=date(2015, 12, 31), help="monogenic time split: pairs dated on or before this day train")
    parser.add_argument("--with-kg-embedding", action="store_true", help="also run baseline B2 (TransE over graph plus training evidence triples)")
    parser.add_argument("--no-fold-cache", action="store_true",
                        help="recompute every random-walk and TransE fold instead of reusing those stored in OUTPUT_DIR/fold_cache under a key of all their inputs "
                             "(the store lets a run killed by a reclaimed container resume; a stored fold equals the recomputed one)")
    parser.add_argument("--with-type-popularity", action="store_true",
                        help="also run popularity and degree_popularity fitted per perturbation type (drug or gene), the graph-free reading of the type offset "
                             "(docs/perturbation_type_offset.md); off by default, so the confirmatory baseline set is unchanged unless the scorer is given these names")
    parser.add_argument("--kg-embedding-all-splits", action="store_true", help="run B2 on the pathway hold-outs too (many fits)")
    parser.add_argument("--metabolic-layer-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("runs/baselines"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/phase2_baselines.md"))
    arguments = parser.parse_args()
    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, metabolic_layer_only=arguments.metabolic_layer_only, group_by=arguments.group_by, label_grades=tuple(arguments.label_grades) if arguments.label_grades else None,
                                label_selection=arguments.label_selection)
    if data.label_selection_summary is not None:
        print(f"label selection: {data.label_selection_summary}")
    coverage = None if arguments.measurement_coverage is None else read_measurement_coverage(arguments.measurement_coverage)
    if arguments.score_lockbox and arguments.lockbox is None:
        raise SystemExit("--score-lockbox needs --lockbox")
    in_lockbox, lockbox_summary = None, None
    if arguments.lockbox is not None:
        in_lockbox = read_lockbox(arguments.lockbox, data, arguments.group_by, arguments.evidence_dir)
        lockbox_summary = {"path": str(arguments.lockbox), "num_lockbox_perturbations": int(in_lockbox.sum()),
                           "role": "scored" if arguments.score_lockbox else "removed before every split"}
        print(f"lockbox: {lockbox_summary}")
        if not arguments.score_lockbox:
            data = restrict_to_perturbations(data, ~in_lockbox)
            in_lockbox = None
    label_mask = data.label_mask  # None without --label-selection
    primary_split = "lockbox" if arguments.score_lockbox else "grouped"
    if arguments.score_lockbox:
        grouped_masks, grouped_labels = [in_lockbox], ["lockbox"]
        pathway_masks, pathway_module_ids, subsystem_masks, subsystem_labels = [], [], [], []
    else:
        fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.seed)
        fold_of_perturbation = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
        grouped_masks, grouped_labels = [fold_of_perturbation == fold for fold in range(arguments.num_folds)], [str(fold) for fold in range(arguments.num_folds)]
        pathway_masks, pathway_module_ids = pathway_wise_test_masks(data, arguments.curated_modules, arguments.min_holdout_positives)
        subsystem_masks, subsystem_labels = subsystem_test_masks(data, arguments.min_holdout_positives)
    normalized_adjacency = build_normalized_adjacency(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0])
    # rows move whole, so a pair keeps its label-mask entry; within the lockbox and the development set apart under --score-lockbox
    permutation_source_row = degree_stratified_row_permutation(data.perturbation_degrees_for_strata, random_seed=arguments.seed, partition=in_lockbox)
    permuted_outcomes = data.outcomes[permutation_source_row]
    permuted_label_mask = None if label_mask is None else label_mask[permutation_source_row]

    num_genes = sum(t == "gene" for t in data.perturbation_types)
    num_drugs = sum(t == "drug" for t in data.perturbation_types)
    results = {
        "num_perturbations": len(data.perturbation_ids), "num_genes": num_genes, "num_drugs": num_drugs, "num_symptoms": len(data.symptoms), "symptoms": data.symptoms,
        "group_by": arguments.group_by, "num_groups": len(set(data.group_ids)), "num_folds": 1 if arguments.score_lockbox else arguments.num_folds, "label_selection": data.label_selection_summary,
        "lockbox": lockbox_summary, "seed": arguments.seed,
        "kg_embedding_settings": KG_EMBEDDING_SETTINGS if arguments.with_kg_embedding else None, "type_popularity": arguments.with_type_popularity,
        "pathway_wise_modules": pathway_module_ids, "pathway_wise_holdout_sizes": [int(mask.sum()) for mask in pathway_masks],
        "pathway_wise_positives": [int(kept_positives(data)[mask].sum()) for mask in pathway_masks],
        "subsystem_wise_subsystems": subsystem_labels, "subsystem_wise_holdout_sizes": [int(mask.sum()) for mask in subsystem_masks],
        "subsystem_wise_positives": [int(kept_positives(data)[mask].sum()) for mask in subsystem_masks],
        "splits": {primary_split: {}, "pathway_wise": {}, "subsystem_wise": {}, f"{primary_split}_label_permutation": {}, "pathway_wise_label_permutation": {}, "subsystem_wise_label_permutation": {},
                   f"{primary_split}_rewired_graph": {}},
    }
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    np.save(arguments.output_dir / "permutation_source_rows.npy", permutation_source_row)  # the permuted labels are data.outcomes[these rows]
    fold_cache = None if arguments.no_fold_cache else arguments.output_dir / "fold_cache"
    if arguments.rewiring_swaps_per_edge > 0:
        original_edges = np.stack([data.edge_source, data.edge_target])
        if arguments.rewiring_method == "walk_graph":  # the walk's own undirected simple graph, every node keeping its number of neighbours
            rewired = rewire_walk_graph(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0],
                                        num_swaps_per_edge=arguments.rewiring_swaps_per_edge, random_seed=arguments.seed)
            real_pairs = rewire_walk_graph(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0], num_swaps_per_edge=0)
            results["rewiring"] = {"method": arguments.rewiring_method, "swaps_per_edge": arguments.rewiring_swaps_per_edge, "seed": arguments.seed,
                                   "num_walk_pairs": int(rewired.shape[1]),
                                   "share_of_walk_pairs_unchanged": float(len(set(map(tuple, real_pairs.T.tolist())) & set(map(tuple, rewired.T.tolist()))) / max(rewired.shape[1], 1))}
        else:
            rewire = fast_degree_preserving_rewiring if arguments.rewiring_method == "integer_draws" else degree_preserving_rewiring
            rewired = rewire(original_edges, data.edge_relation, num_swaps_per_edge=arguments.rewiring_swaps_per_edge, random_seed=arguments.seed)
            results["rewiring"] = {"method": arguments.rewiring_method, "swaps_per_edge": arguments.rewiring_swaps_per_edge, "seed": arguments.seed,
                                   "share_of_edges_unchanged": float((rewired == original_edges).all(axis=0).mean()),
                                   "duplicate_edges_before": duplicate_edge_count(original_edges, data.edge_relation), "duplicate_edges_after": duplicate_edge_count(rewired, data.edge_relation)}
        rewired_adjacency = build_normalized_adjacency(len(data.node_ids), rewired[0], rewired[1], np.where(data.is_currency)[0])
        predictions, rows, per_fold, fold_of_row = run_split(data, data.outcomes, grouped_masks, "random_walk_with_restart", arguments.restart_probability, rewired_adjacency, arguments.min_fold_size_for_macro,
                                                             grouped_labels, label_mask=label_mask, cache_directory=fold_cache)
        results["splits"][f"{primary_split}_rewired_graph"]["random_walk_with_restart"] = score(data, predictions, data.outcomes, rows, per_fold, arguments.num_bootstrap, fold_of_row, label_mask, coverage)
        np.save(arguments.output_dir / f"predictions_{primary_split}_rewired_graph_random_walk_with_restart.npy", predictions)
        entry = results["splits"][f"{primary_split}_rewired_graph"]["random_walk_with_restart"]
        print(f"{primary_split + '_rewired_graph':26s} {'random_walk_with_restart':26s} macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  macro AUROC {entry['macro_auroc']:.3f}")
    split_plan = [(primary_split, grouped_masks, grouped_labels), ("pathway_wise", pathway_masks, pathway_module_ids), ("subsystem_wise", subsystem_masks, subsystem_labels)]
    model_names = list(BASELINE_NAMES) + ([KG_EMBEDDING_NAME] if arguments.with_kg_embedding else []) + (list(TYPE_POPULARITY_NAMES) if arguments.with_type_popularity else [])
    for model_name in model_names:
        for split_name, masks, labels in split_plan:
            if not masks:
                continue
            if model_name == KG_EMBEDDING_NAME and split_name != primary_split and not arguments.kg_embedding_all_splits:
                continue
            partner_groups = None if split_name == primary_split else data.group_ids  # pathway-wise hold-outs are chosen by seed gene
            predictions, rows, per_fold, fold_of_row = run_split(data, data.outcomes, masks, model_name, arguments.restart_probability, normalized_adjacency, arguments.min_fold_size_for_macro, labels,
                                                                 label_mask=label_mask, random_seed=arguments.seed, group_ids=partner_groups, cache_directory=fold_cache)
            results["splits"][split_name][model_name] = score(data, predictions, data.outcomes, rows, per_fold, arguments.num_bootstrap, fold_of_row, label_mask, coverage)
            arguments.output_dir.mkdir(parents=True, exist_ok=True)
            np.save(arguments.output_dir / f"predictions_{split_name}_{model_name}.npy", predictions)  # pooled out-of-split predictions for paired comparisons
            np.save(arguments.output_dir / f"scored_rows_{split_name}.npy", rows)
            if not arguments.skip_permutation_control:
                predictions, rows, per_fold, fold_of_row = run_split(data, permuted_outcomes, masks, model_name, arguments.restart_probability, normalized_adjacency, arguments.min_fold_size_for_macro, labels,
                                                                     label_mask=permuted_label_mask, random_seed=arguments.seed, group_ids=partner_groups, cache_directory=fold_cache)
                results["splits"][f"{split_name}_label_permutation"][model_name] = score(data, predictions, permuted_outcomes, rows, per_fold, arguments.num_bootstrap, fold_of_row, permuted_label_mask, coverage)
                np.save(arguments.output_dir / f"predictions_{split_name}_label_permutation_{model_name}.npy", predictions)
        for split_name, entries in results["splits"].items():
            if model_name in entries:
                entry = entries[model_name]
                print(f"{split_name:26s} {model_name:26s} macro AUPRC {entry['macro_auprc']:.3f} (per fold {entry['per_fold_macro_auprc_mean']:.3f} ± {entry['per_fold_macro_auprc_sd']:.3f})  "
                      f"macro AUROC {entry['macro_auroc']:.3f} (per fold {entry['per_fold_macro_auroc_mean']:.3f}; rank within hold-out AUPRC {entry.get('macro_auprc_rank_normalised', float('nan')):.3f} AUROC {entry.get('macro_auroc_rank_normalised', float('nan')):.3f} stratified {entry.get('macro_auroc_stratified', float('nan')):.3f})  MRR {entry['mean_reciprocal_rank']:.3f}  hits@3 {entry['hits_at_3']:.3f}")
    results["time_split"] = {}
    if label_mask is not None or arguments.score_lockbox:
        model_names_for_time_split = []  # the time split scores new positive pairs; it is not defined under a label selection, and it would score the lockbox
        results["time_split_not_run"] = "not run with --score-lockbox" if arguments.score_lockbox else "not implemented with --label-selection"
    else:
        model_names_for_time_split = model_names
    for model_name in model_names_for_time_split:
        entry = time_split_evaluation(data, arguments.time_split_cutoff, model_name, arguments.restart_probability, normalized_adjacency, arguments.seed)
        if entry is not None:
            results["time_split"][model_name] = entry
            print(f"{'time_split':26s} {model_name:26s} macro AUPRC {entry['macro_auprc']:.3f} (permuted {entry['macro_auprc_permuted']:.3f})  macro AUROC {entry['macro_auroc']:.3f} (permuted {entry['macro_auroc_permuted']:.3f})  new positives {entry['new_positive_pairs']} of {entry['scored_pairs']} scored pairs")
    results["time_split_publication_dated"] = {}
    for model_name in model_names_for_time_split:
        entry = time_split_evaluation(data, arguments.time_split_cutoff, model_name, arguments.restart_probability, normalized_adjacency, arguments.seed, require_publication_date=True)
        if entry is not None:
            results["time_split_publication_dated"][model_name] = entry
            print(f"{'time_split_pubdated':26s} {model_name:26s} macro AUPRC {entry['macro_auprc']:.3f} (permuted {entry['macro_auprc_permuted']:.3f})  macro AUROC {entry['macro_auroc']:.3f} (permuted {entry['macro_auroc_permuted']:.3f})  new positives {entry['new_positive_pairs']} (curation-dated excluded {entry['curation_dated_new_positives_excluded']})")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    (arguments.output_dir / "results.json").write_text(json.dumps(results, indent=1))
    (arguments.output_dir / "perturbation_ids.json").write_text(json.dumps(data.perturbation_ids))

    lines = ["# Phase 2 baselines (generated by experiments/run_baselines.py)", "",
             f"{len(data.perturbation_ids)} perturbations ({num_genes} genes, {num_drugs} drugs), {len(data.symptoms)} symptoms, relation induces. "
             + (f"Lockbox: fitted once on the {int((~in_lockbox).sum())} development perturbations, scored on the {int(in_lockbox.sum())} lockbox perturbations of {arguments.lockbox} (leakage groups by {arguments.group_by}). " if arguments.score_lockbox else
              f"Grouped split: {arguments.num_folds}-fold cross-validation with leakage groups by {arguments.group_by} ({len(set(data.group_ids))} groups); out-of-fold predictions pooled. "
              + (f"The {lockbox_summary['num_lockbox_perturbations']} lockbox perturbations of {arguments.lockbox} were removed first. " if lockbox_summary else ""))
             + "Unobserved pairs count as negatives (positive-unlabelled convention). " + (f"Label selection {data.label_selection_summary['path']}: {data.label_selection_summary['masked_pairs']} of {data.label_selection_summary['positive_pairs']} positive pairs set aside, neither positive nor negative in fitting or scoring (positive counts in the hold-out tables include them). " if data.label_selection_summary else "") + "Pooled metrics are computed on the pooled out-of-fold predictions; per-fold values are the mean and standard deviation over folds of the macro metric, which is immune to the base-rate artifact that pulls pooled AUROC of a constant-per-fold predictor below 0.5.", ""]
    header = ["| model | macro AUPRC (pooled) | macro AUPRC (per fold) | macro AUROC (pooled) | macro AUROC (per fold) | MRR | hits@3 | scored perturbations |", "|---|---|---|---|---|---|---|---|"]
    lines += ["## Lockbox" if arguments.score_lockbox else "## Grouped perturbation-wise split", ""] + header + [summary_row(name, entry) for name, entry in results["splits"][primary_split].items()] + [""]
    degree_bins = list(next(iter(results["splits"][primary_split].values()))["macro_auprc_by_degree_bin"])
    lines += ["Macro AUPRC by perturbation degree tercile (pooled out-of-fold predictions; bins are [low degree, high degree] with the number of perturbations):", "",
              "| model | " + " | ".join(degree_bins) + " |", "|---|" + "---|" * len(degree_bins)]
    for name, entry in results["splits"][primary_split].items():
        lines.append(f"| {name} | " + " | ".join(f"{entry['macro_auprc_by_degree_bin'][b]:.3f}" for b in degree_bins) + " |")
    lines.append("")
    coverage_strata = list(next(iter(results["splits"][primary_split].values()))["macro_auprc_by_measurement_coverage"])
    if coverage_strata:
        lines += [f"Macro AUPRC inside strata of binding measurement coverage ({arguments.measurement_coverage}); a stratum where no symptom reaches "
                  f"{MINIMUM_POSITIVES_TO_SCORE} positives reads nan. A model whose advantage sits only in the most measured stratum is reading how well "
                  "the drug was studied:", "",
                  "| model | " + " | ".join(coverage_strata) + " |", "|---|" + "---|" * len(coverage_strata)]
        for name, entry in results["splits"][primary_split].items():
            lines.append(f"| {name} | " + " | ".join(f"{entry['macro_auprc_by_measurement_coverage'][stratum]:.3f}" for stratum in coverage_strata) + " |")
        lines.append("")
    if results["splits"]["pathway_wise"]:
        lines += holdout_section("Pathway-wise split, curated modules (each module of design section 3.2 held out in turn)",
                                 "Per-symptom AUPRC inside one held-out module is not meaningful (sets of 4 to 8 genes with homogeneous symptom profiles), so only pooled per-symptom metrics and per-hold-out ranking metrics are shown.",
                                 results["splits"]["pathway_wise"], results["splits"]["pathway_wise_label_permutation"], pathway_module_ids, pathway_masks, kept_positives(data))
    if results["splits"]["subsystem_wise"]:
        lines += holdout_section("Pathway-wise split, Human-GEM subsystems (every gene whose primary subsystem is the held-out one)",
                                 "Subsystems are the reconstruction's own pathway partition; they cover many more annotated genes than the curated modules and are the candidate definition of the pre-registered pathway-wise split.",
                                 results["splits"]["subsystem_wise"], results["splits"]["subsystem_wise_label_permutation"], subsystem_labels, subsystem_masks, kept_positives(data))
    split_title = "lockbox" if arguments.score_lockbox else "grouped split"
    if results["splits"][f"{primary_split}_label_permutation"]:
        lines += [f"## Negative control: labels permuted within degree strata ({split_title})", ""] + header + [summary_row(name, entry) for name, entry in results["splits"][f"{primary_split}_label_permutation"].items()] + [""]
    if results["splits"][f"{primary_split}_rewired_graph"]:
        lines += [f"## Negative control: degree-preserving rewiring of the graph ({arguments.rewiring_swaps_per_edge} swaps per edge, {arguments.rewiring_method}; {split_title}, random walk only)", ""] + header + [summary_row(name, entry) for name, entry in results["splits"][f"{primary_split}_rewired_graph"].items()] + [""]
    lines += [f"## Per-symptom AUPRC, {split_title}", "", "Point estimate with 95 percent bootstrap interval over perturbations; base rate in parentheses.", "",
              "| symptom | " + " | ".join(results["splits"][primary_split]) + " |", "|---|" + "---|" * len(results["splits"][primary_split])]
    for symptom in data.symptoms:
        cells = []
        for entry in results["splits"][primary_split].values():
            s = entry["per_symptom"].get(symptom)
            cells.append("n/a" if s is None else f"{s['auprc']['point']:.3f} [{s['auprc']['lower']:.3f}, {s['auprc']['upper']:.3f}] ({s['base_rate']:.3f})")
        lines.append(f"| {symptom} | " + " | ".join(cells) + " |")
    if results["time_split"]:
        first = next(iter(results["time_split"].values()))
        lines += ["", f"## Time split, monogenic pairs by availability date (cutoff {first['cutoff']})", "",
                  f"Training positives: {first['training_positive_pairs']} pairs dated on or before the cutoff; new positives: {first['new_positive_pairs']} pairs dated after it; undated (Orphanet-only) positives excluded: {first['undated_positive_pairs']}; scored pairs: {first['scored_pairs']} (every pair that is neither a training positive nor an undated positive). "
                  "A pair's date is the publication date of the PubMed reference cited by its OMIM annotation (docs/hpo_reference_publication_dates.json) or, for annotations citing only the OMIM entry, the HPO biocuration date, so this split measures prediction of links published later. Permuted: new-positive labels permuted among the scored pairs of each symptom.", "",
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
    if results.get("time_split_publication_dated"):
        first = next(iter(results["time_split_publication_dated"].values()))
        lines += ["", "### Time split with publication-dated new positives only", "",
                  f"Pairs whose only post-cutoff date is an HPO curation date are excluded from the new-positive set (review v0.4, finding 7): {first['curation_dated_new_positives_excluded']} excluded; new positives {first['new_positive_pairs']}; scored pairs {first['scored_pairs']}.", "",
                  "| model | macro AUPRC | macro AUPRC (permuted) | macro AUROC | macro AUROC (permuted) |", "|---|---|---|---|---|"]
        for name, entry in results["time_split_publication_dated"].items():
            lines.append(f"| {name} | {entry['macro_auprc']:.3f} | {entry['macro_auprc_permuted']:.3f} | {entry['macro_auroc']:.3f} | {entry['macro_auroc_permuted']:.3f} |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
