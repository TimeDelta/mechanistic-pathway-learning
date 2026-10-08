"""Is the random walk's score on the real graph outside the null distribution of degree-preserving rewirings?

The committed rewiring control (run_baselines.py --rewiring-swaps-per-edge) rewires once, at one seed, and reports a
single number. That cannot say whether the gap between the real graph and the rewired one is larger than the spread
across rewirings, so this repeats the rewiring over many seeds with everything else held fixed -- the same graph, the
same evidence, the same fold assignment, the same restart probability -- and places the real score against the
distribution of rewired scores. The reported p-value is the fraction of rewirings scoring at least as well as the
real graph, so it is bounded below by 1 / (number of seeds + 1).

Only the random walk is run, which is what makes many seeds affordable: the rewiring destroys which node joins which
while keeping every node's in- and out-degree per relation type, so a model reading only connectivity should be
unaffected and a model reading the wiring should fall.

--rewiring-method walk_graph (8 October 2026) rewires the walk's own graph, the undirected simple graph it builds, so
every node keeps the number of neighbours the walk sees (negative_controls.rewire_walk_graph); pair_sampling, the
default and the method of every earlier result, rewires the stored directed edges per relation, which moves apart two
stored edges joining one pair (binds in both directions, a reversible reaction's substrate and product edges) and so
gives the walk a denser graph than the real one.

Usage:
  OMP_NUM_THREADS=1 python experiments/run_rewiring_null_distribution.py --num-rewirings 20 --rewiring-method walk_graph
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.run_baselines import run_split, score
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.negative_controls import degree_preserving_rewiring, rewire_walk_graph
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds
from mechanistic_pathway_learning.models.baselines.random_walk_with_restart_baseline import build_normalized_adjacency


def per_fold_macro_auprc(data, test_masks, normalized_adjacency, restart_probability, minimum_fold_size, num_bootstrap):
    predictions, rows, per_fold, fold_of_row = run_split(data, data.outcomes, test_masks, "random_walk_with_restart",
                                                         restart_probability, normalized_adjacency, minimum_fold_size)
    entry = score(data, predictions, data.outcomes, rows, per_fold, num_bootstrap, fold_of_row)
    return entry["per_fold_macro_auprc_mean"], entry["macro_auprc"], entry["macro_auroc"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence"))
    parser.add_argument("--group-by", default="disease_cluster")
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--split-seed", type=int, default=0, help="fold assignment; held fixed across rewirings")
    parser.add_argument("--num-rewirings", type=int, default=20)
    parser.add_argument("--swaps-per-edge", type=int, default=50)
    parser.add_argument("--rewiring-method", choices=["pair_sampling", "walk_graph"], default="pair_sampling")
    parser.add_argument("--restart-probability", type=float, default=0.3)
    parser.add_argument("--min-fold-size-for-macro", type=int, default=10)
    parser.add_argument("--num-bootstrap", type=int, default=200)
    parser.add_argument("--output-json", type=Path, default=Path("runs/rewiring_check/null_distribution.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/rewiring_null_distribution.md"))
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=arguments.group_by)
    fold_by_perturbation = assign_grouped_folds(data.perturbation_ids, data.group_ids, arguments.num_folds, arguments.split_seed)
    fold_of_perturbation = np.array([fold_by_perturbation[p] for p in data.perturbation_ids])
    test_masks = [fold_of_perturbation == fold for fold in range(arguments.num_folds)]
    common = dict(test_masks=test_masks, restart_probability=arguments.restart_probability,
                  minimum_fold_size=arguments.min_fold_size_for_macro, num_bootstrap=arguments.num_bootstrap)

    real_adjacency = build_normalized_adjacency(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0])
    real_per_fold, real_pooled, real_auroc = per_fold_macro_auprc(data, normalized_adjacency=real_adjacency, **common)
    print(f"real graph: per fold {real_per_fold:.4f}, pooled {real_pooled:.4f}, AUROC {real_auroc:.4f}", flush=True)

    rewired = []
    for seed in range(arguments.num_rewirings):
        if arguments.rewiring_method == "walk_graph":
            edges = rewire_walk_graph(len(data.node_ids), data.edge_source, data.edge_target, np.where(data.is_currency)[0],
                                      num_swaps_per_edge=arguments.swaps_per_edge, random_seed=seed)
        else:
            edges = degree_preserving_rewiring(np.stack([data.edge_source, data.edge_target]), data.edge_relation,
                                               num_swaps_per_edge=arguments.swaps_per_edge, random_seed=seed)
        adjacency = build_normalized_adjacency(len(data.node_ids), edges[0], edges[1], np.where(data.is_currency)[0])
        per_fold, pooled, auroc = per_fold_macro_auprc(data, normalized_adjacency=adjacency, **common)
        rewired.append({"seed": seed, "per_fold_macro_auprc": per_fold, "pooled_macro_auprc": pooled, "macro_auroc": auroc})
        print(f"  rewiring seed {seed:2d}: per fold {per_fold:.4f}, pooled {pooled:.4f}, AUROC {auroc:.4f}", flush=True)

    draws = np.array([entry["per_fold_macro_auprc"] for entry in rewired])
    at_least_as_good = int((draws >= real_per_fold).sum())
    p_value = (at_least_as_good + 1) / (len(draws) + 1)
    summary = {"real": {"per_fold_macro_auprc": real_per_fold, "pooled_macro_auprc": real_pooled, "macro_auroc": real_auroc},
               "rewirings": rewired, "swaps_per_edge": arguments.swaps_per_edge, "rewiring_method": arguments.rewiring_method,
               "rewired_mean": float(draws.mean()), "rewired_sd": float(draws.std(ddof=1)) if len(draws) > 1 else 0.0,
               "rewired_min": float(draws.min()), "rewired_max": float(draws.max()),
               "rewirings_at_least_as_good_as_real": at_least_as_good, "p_value": p_value}
    arguments.output_json.parent.mkdir(parents=True, exist_ok=True)
    arguments.output_json.write_text(json.dumps(summary, indent=1) + "\n")

    lines = [f"# Random walk on the real graph against the null distribution of {len(draws)} degree-preserving rewirings",
             "", f"Generated by experiments/run_rewiring_null_distribution.py. Graph {arguments.graph_dir}, evidence "
             f"{arguments.evidence_dir}, {arguments.group_by} grouped split with {arguments.num_folds} folds held fixed "
             f"at split seed {arguments.split_seed}, {arguments.swaps_per_edge} degree-preserving double-edge swaps per "
             "edge, random walk with restart only. Rewiring method: "
             + ("walk_graph, the walk's undirected simple graph rewired so every node keeps its number of neighbours." if arguments.rewiring_method == "walk_graph"
                else "pair_sampling, the stored directed edges per relation (the walk's graph comes out denser than the real one).") + "", "",
             "| | per-fold macro AUPRC | pooled macro AUPRC | macro AUROC |", "|---|---|---|---|",
             f"| real graph | **{real_per_fold:.3f}** | {real_pooled:.3f} | {real_auroc:.3f} |",
             f"| rewirings, mean over {len(draws)} seeds | {draws.mean():.3f} ± {summary['rewired_sd']:.3f} | | |",
             f"| rewirings, best of {len(draws)} | {draws.max():.3f} | | |",
             f"| rewirings, worst of {len(draws)} | {draws.min():.3f} | | |", "",
             f"Rewirings scoring at least as well as the real graph: **{at_least_as_good} of {len(draws)}**, "
             f"so p = {p_value:.3f} (bounded below by {1 / (len(draws) + 1):.3f} at this number of seeds).", "",
             "| rewiring seed | per fold | pooled | AUROC |", "|---|---|---|---|"]
    lines += [f"| {e['seed']} | {e['per_fold_macro_auprc']:.3f} | {e['pooled_macro_auprc']:.3f} | {e['macro_auroc']:.3f} |" for e in rewired]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"\nreal {real_per_fold:.4f} against rewired {draws.mean():.4f} ± {summary['rewired_sd']:.4f}; "
          f"{at_least_as_good} of {len(draws)} at least as good, p = {p_value:.3f}")
    print(f"wrote {arguments.markdown_output} and {arguments.output_json}")


if __name__ == "__main__":
    main()
