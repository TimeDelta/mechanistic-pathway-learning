"""Compare degree definitions on the split graph against the merged graph's seed degree, on the full graph and the slice.

The registered degree stratification is the merged graph's seed degree. A split graph gives a perturbed gene node only
the edges that are about the gene, so the question is which degree on the split graph reproduces it
(docs/gene_protein_split_results.md, docs/preregistration.md amendment of 9 October 2026). Variants per perturbation:
  seed                 sum of the degrees of the perturbation's seed nodes (what is registered on a merged graph)
  one_hop_all          seed nodes and every node one edge away
  one_hop_encodes      seed nodes and the nodes they reach over the encodes relation, less the encodes edges themselves
                       (what ExperimentData.perturbation_degrees_for_strata returns on a split graph)

    python experiments/compare_split_degree_definitions.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data


def terciles(values: np.ndarray, num_bins: int = 3) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    bin_of_row = np.empty(len(values), dtype=int)
    bin_of_row[order] = np.minimum(np.arange(len(values)) * num_bins // len(values), num_bins - 1)
    return bin_of_row


def every_neighbour_degrees(data) -> np.ndarray:
    """Summed degree of the seeds and of every node one edge away; kept here only as the comparison."""
    neighbours_of_node = [[] for _ in data.node_ids]
    for source, target in zip(data.edge_source, data.edge_target):
        neighbours_of_node[source].append(target)
        neighbours_of_node[target].append(source)
    degrees = []
    for seeds in data.perturbation_seeds:
        reached = {int(seed) for seed in seeds}
        for seed in seeds:
            reached.update(neighbours_of_node[int(seed)])
        degrees.append(float(data.node_degree[sorted(reached)].sum()) if reached else 0.0)
    return np.array(degrees)


def encodes_hop_degrees(data) -> np.ndarray:
    relation_index = list(data.relation_types).index("encodes")
    partners_of_node = [[] for _ in data.node_ids]
    for source, target, relation in zip(data.edge_source, data.edge_target, data.edge_relation):
        if relation != relation_index:
            continue
        partners_of_node[source].append(target)
        partners_of_node[target].append(source)
    degrees = []
    for seeds in data.perturbation_seeds:
        reached, encodes_edges = set(int(seed) for seed in seeds), 0
        for seed in seeds:
            partners = partners_of_node[int(seed)]
            encodes_edges += len(partners)
            reached.update(partners)
        # an encodes edge adds one to the gene's degree and one to the protein's, and the merged node had neither
        degrees.append(float(data.node_degree[sorted(reached)].sum() - 2 * encodes_edges))
    return np.array(degrees)


def compare(merged_dir: str, split_dir: str, evidence: str, selection: str | None, group_by: str) -> dict:
    merged = load_experiment_data(Path(merged_dir), Path(evidence), group_by=group_by, label_selection=Path(selection) if selection else None)
    split = load_experiment_data(Path(split_dir), Path(evidence), group_by=group_by, label_selection=Path(selection) if selection else None)
    assert list(merged.perturbation_ids) == list(split.perturbation_ids)
    reference = merged.perturbation_degrees
    variants = {"seed": split.perturbation_degrees, "one_hop_all": every_neighbour_degrees(split), "one_hop_encodes": encodes_hop_degrees(split)}
    return {"num_perturbations": len(reference),
            "merged_seed_degree_median": float(np.median(reference)),
            **{name: {"median": float(np.median(values)), "max": float(values.max()),
                      "spearman_against_merged_seed": round(float(spearmanr(values, reference).statistic), 4),
                      "same_tercile_share": round(float((terciles(values) == terciles(reference)).mean()), 4),
                      "exactly_equal_share": round(float((values == reference).mean()), 4)}
               for name, values in variants.items()}}


result = {
    "full": compare("data/processed/graph_full_neuronal", "data/processed/graph_full_neuronal_split", "data/processed/evidence_full_v2",
                    "data/processed/label_selection/better_v2_full_v2.parquet", "disease_cluster_and_targets"),
    "slice": compare("data/processed/graph", "data/processed/graph_split", "data/processed/evidence", None, "disease_cluster"),
}
print(json.dumps(result, indent=1))
output = Path("data/processed/off_target_scoping").parent / "split_degree_definitions.json"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=1) + "\n")
print(f"wrote {output}")
