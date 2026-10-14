"""Does hop-limited reached-node grouping of GENE perturbations balloon, and what does it catch that disease clusters
miss? (Open decision 9.36.)

The owner, 14 October 2026: "why not split on downstream nodes touched instead of direct target? i guess that's kind of
like the graph proximity metric you mentioned but i'm sure this could balloon very quickly so, if that happens, maybe
limiting it to a certain number of hops?", and then "for that ballooning hop idea, i was only talking about the gene
perturbations".

Two questions, because the answers go opposite ways. Is there a hop limit that both catches sharing and leaves a usable
split? And of the sharing it would catch, how much does the registered disease-cluster grouping already catch?

Usage:
  python -I experiments/measure_gene_hop_grouping.py
"""
import sys
from pathlib import Path
import numpy as np
import torch
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.models.forward_cone_divisor import forward_reachable_nodes

data = load_experiment_data(Path("data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes"),
                            Path("data/processed/evidence_full_v4"), group_by="disease_cluster_and_targets",
                            label_selection=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
num_nodes = len(data.node_ids)
edge_source = torch.as_tensor(np.asarray(data.edge_source), dtype=torch.long)
edge_target = torch.as_tensor(np.asarray(data.edge_target), dtype=torch.long)

gene_rows = [index for index, kind in enumerate(data.perturbation_types) if kind != "drug"]
gene_identifiers = [str(data.perturbation_ids[row]) for row in gene_rows]
disease_cluster_of_gene = {str(data.perturbation_ids[row]): str(data.group_ids[row]) for row in gene_rows}
print(f"gene perturbations: {len(gene_rows)}; graph nodes: {num_nodes}; edges: {len(edge_source)}")
print(f"distinct disease-cluster groups over those genes: {len(set(disease_cluster_of_gene.values()))}")

class UnionFind:
    def __init__(self, size): self.parent = list(range(size))
    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]; item = self.parent[item]
        return item
    def union(self, left, right):
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root: self.parent[right_root] = left_root

for maximum_hops in (1, 2, 3):
    genes_reaching_node = {}
    reached_counts = []
    for position, row in enumerate(gene_rows):
        seeds = torch.as_tensor(np.asarray(data.perturbation_seeds[row]), dtype=torch.long)
        reached = forward_reachable_nodes(edge_source, edge_target, seeds, num_nodes, maximum_steps=maximum_hops)
        reached_node_indices = torch.nonzero(reached).reshape(-1).tolist()
        reached_counts.append(len(reached_node_indices))
        for node_index in reached_node_indices:
            genes_reaching_node.setdefault(node_index, []).append(position)

    union_find = UnionFind(len(gene_rows))
    shared_pairs_total = 0
    shared_pairs_already_in_one_disease_cluster = 0
    for node_index, positions in genes_reaching_node.items():
        if len(positions) < 2:
            continue
        for other in positions[1:]:
            union_find.union(positions[0], other)
        if len(positions) <= 40:  # count pairs only where it is cheap; a hub node is reported separately
            for i in range(len(positions)):
                for j in range(i + 1, len(positions)):
                    shared_pairs_total += 1
                    if disease_cluster_of_gene[gene_identifiers[positions[i]]] == disease_cluster_of_gene[gene_identifiers[positions[j]]]:
                        shared_pairs_already_in_one_disease_cluster += 1

    component_of = {}
    for position in range(len(gene_rows)):
        component_of.setdefault(union_find.find(position), []).append(position)
    sizes = sorted((len(members) for members in component_of.values()), reverse=True)
    hub_nodes = sum(1 for positions in genes_reaching_node.values() if len(positions) > 40)
    print(f"\n--- k = {maximum_hops} hops ---")
    print(f"  reached nodes per gene: median {int(np.median(reached_counts))}, max {max(reached_counts)}")
    print(f"  groups: {len(sizes)}; largest {sizes[0]} ({sizes[0]/len(gene_rows):.1%} of genes); next {sizes[1:5]}")
    print(f"  genes in a group of one: {sum(1 for s in sizes if s == 1)}")
    print(f"  nodes reached by >40 genes (hubs, pair counting skipped): {hub_nodes}")
    if shared_pairs_total:
        print(f"  gene pairs sharing a reached node (non-hub nodes only): {shared_pairs_total}; "
              f"already in one disease cluster: {shared_pairs_already_in_one_disease_cluster} "
              f"({shared_pairs_already_in_one_disease_cluster/shared_pairs_total:.1%})")
    else:
        print("  no gene pair shares a non-hub reached node")
