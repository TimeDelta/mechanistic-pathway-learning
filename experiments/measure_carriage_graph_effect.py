"""What the plasma carriage edges change for the confirmatory reading, measured rather than asserted.

The confirmatory family trains on data/processed/graph_full_neuronal_split_binders (the user, 9 October 2026: "yes,
include binder edges in confirmatory graph"; docs/preregistration.md, amendment of 9 October 2026 (fourth)). That
variant is the split graph plus `binds` edges from the extracellular copy of each cargo metabolite to its binder
protein. The edges land on protein nodes, but `perturbation_degrees_for_strata` hops from a perturbed gene node over
`encodes` to its protein, so a binder gene's strata degree moves even though its own seed node is untouched. Whether
any perturbation thereby moves degree quintile decides whether the registered within-strata readings and the
within-stratum label permutation change, so it is measured here.

Reads both graph directories and writes its report. Trains nothing.

Usage:
  python experiments/measure_carriage_graph_effect.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import degree_strata


def added_edge_rows(source_graph: Path, variant_graph: Path) -> pd.DataFrame:
    """The edge rows the variant holds and the source graph does not."""
    source_edges = pd.read_parquet(source_graph / "edges.parquet")
    variant_edges = pd.read_parquet(variant_graph / "edges.parquet")
    merged = variant_edges.merge(source_edges, how="outer", indicator=True)
    return merged[merged["_merge"] == "left_only"].drop(columns="_merge")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source-graph", type=Path, default=Path("data/processed/graph_full_neuronal_split"))
    parser.add_argument("--variant-graph", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v2.parquet"))
    parser.add_argument("--group-by", default="disease_cluster_and_targets")
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v2.json"))
    parser.add_argument("--output", type=Path, default=Path("docs/plasma_binder_confirmatory_effect.md"))
    arguments = parser.parse_args()

    added = added_edge_rows(arguments.source_graph, arguments.variant_graph)
    nodes_touched = sorted(set(added["source_id"]) | set(added["target_id"]))

    without_carriage = load_experiment_data(arguments.source_graph, arguments.evidence_dir, group_by=arguments.group_by,
                                            label_selection=arguments.label_selection)
    with_carriage = load_experiment_data(arguments.variant_graph, arguments.evidence_dir, group_by=arguments.group_by,
                                         label_selection=arguments.label_selection)
    perturbation_ids = list(without_carriage.perturbation_ids)
    unchanged_inputs = {
        "perturbation ids": perturbation_ids == list(with_carriage.perturbation_ids),
        "symptoms": list(without_carriage.symptoms) == list(with_carriage.symptoms),
        "leakage group ids": list(without_carriage.group_ids) == list(with_carriage.group_ids),
        "outcome matrix": np.array_equal(np.asarray(without_carriage.outcomes), np.asarray(with_carriage.outcomes)),
        "label mask": np.array_equal(np.asarray(without_carriage.label_mask), np.asarray(with_carriage.label_mask)),
        "weights": np.array_equal(np.asarray(without_carriage.weights), np.asarray(with_carriage.weights)),
    }

    node_index_of_id = {node_id: index for index, node_id in enumerate(pd.read_parquet(arguments.variant_graph / "nodes.parquet")["node_id"])}
    touched_indices = {node_index_of_id[node_id] for node_id in nodes_touched}
    perturbations_seeded_on_touched_nodes = [identifier for identifier, seeds in zip(perturbation_ids, without_carriage.perturbation_seeds)
                                             if {int(seed) for seed in seeds} & touched_indices]

    degrees_without = np.asarray(without_carriage.perturbation_degrees_for_strata, dtype=float)
    degrees_with = np.asarray(with_carriage.perturbation_degrees_for_strata, dtype=float)
    strata_without, strata_with = np.asarray(degree_strata(degrees_without)), np.asarray(degree_strata(degrees_with))
    in_lockbox = np.array([identifier in set(json.loads(arguments.lockbox.read_text())["perturbation_ids"]) for identifier in perturbation_ids])

    changed_degree = np.flatnonzero(degrees_without != degrees_with)
    moved_stratum = np.flatnonzero(strata_without != strata_with)

    lines = ["# What the plasma carriage edges change for the confirmatory reading", "",
             f"Generated by `experiments/measure_carriage_graph_effect.py` from `{arguments.source_graph}` and "
             f"`{arguments.variant_graph}`. Do not edit by hand.", "",
             f"{len(added)} added edge rows, relations {sorted(added['relation_type'].unique())}, signs "
             f"{sorted(float(sign) for sign in added['sign'].unique())}, touching {len(nodes_touched)} nodes "
             f"({sum(not node.startswith('PROTEIN:') for node in nodes_touched)} cargo metabolites, "
             f"{sum(node.startswith('PROTEIN:') for node in nodes_touched)} binder proteins).", "",
             "## Inputs the edges leave identical", ""]
    for name, identical in unchanged_inputs.items():
        lines.append(f"- {name}: {'identical' if identical else 'CHANGED'}")
    lines += ["", f"- perturbations whose own seed nodes an added edge touches: "
                  f"{', '.join(perturbations_seeded_on_touched_nodes) if perturbations_seeded_on_touched_nodes else 'none'} "
                  f"(the carriage sits on protein nodes; gene knockouts seed gene nodes)", "",
              "## Strata degree, which hops over `encodes` and so does move", "",
              "| perturbation | in lockbox | degree without | degree with | stratum without | stratum with |",
              "|---|---|---|---|---|---|"]
    for index in changed_degree:
        lines.append(f"| {perturbation_ids[index]} | {'yes' if in_lockbox[index] else 'no'} | {degrees_without[index]:.0f} | "
                     f"{degrees_with[index]:.0f} | {strata_without[index]} | {strata_with[index]} |")
    lines += ["", f"Perturbations that move degree quintile: "
                  f"{', '.join(perturbation_ids[index] for index in moved_stratum) if len(moved_stratum) else 'none'}; "
                  f"of them in the lockbox: "
                  f"{', '.join(perturbation_ids[index] for index in moved_stratum if in_lockbox[index]) if any(in_lockbox[index] for index in moved_stratum) else 'none'}.", "",
              f"- quintile sizes over all {len(perturbation_ids)} perturbations: "
              f"{np.bincount(strata_without, minlength=5).tolist()} without, {np.bincount(strata_with, minlength=5).tolist()} with",
              f"- quintile sizes over the {int(in_lockbox.sum())} lockbox perturbations: "
              f"{np.bincount(strata_without[in_lockbox], minlength=5).tolist()} without, "
              f"{np.bincount(strata_with[in_lockbox], minlength=5).tolist()} with", "",
              "The two within-strata readings of H1 rank the scores inside these quintiles, and the secondary label "
              "permutation permutes rows inside them, so a quintile that moves changes both. The quintile of every "
              "lockbox perturbation is what those readings use."]
    arguments.output.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
