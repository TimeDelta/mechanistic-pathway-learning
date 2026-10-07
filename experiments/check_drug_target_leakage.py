"""Does the grouped split separate a drug from the loss of function of its own target gene? (Writes docs/drug_target_leakage.md.)

Under the disease-cluster grouping a gene is held out with the genes that share a disease entry, and a drug with the drugs
of the same target set. Nothing keeps a drug and its target gene in one fold, so a model can learn "perturbing this node
causes this symptom" from the gene in training and be credited for it on the drug in test. This script counts, for each
grouping: the drugs with a labelled target gene in another fold, the kept drug positives that are also kept positives of
such a gene, how often a drug is positive for a symptom when one of its target genes is (against the drugs' base rate),
and the size of the largest leakage group and fold.

Usage:
  python experiments/check_drug_target_leakage.py --evidence-dir data/processed/evidence_full_v2 \
      --label-selection data/processed/label_selection/better_v1_full_v2.parquet
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import assign_grouped_folds

GROUPINGS = ("disease_cluster", "disease_cluster_and_targets")


def leakage_counts(data, num_folds: int, seed: int) -> dict:
    folds = assign_grouped_folds(data.perturbation_ids, data.group_ids, num_folds, seed)
    kept = data.outcomes * (data.label_mask if data.label_mask is not None else 1.0)
    gene_row = {p: i for i, (p, kind) in enumerate(zip(data.perturbation_ids, data.perturbation_types)) if kind == "gene"}
    gene_of_node = {data.node_ids[data.perturbation_seeds[i][0]]: p for p, i in gene_row.items() if len(data.perturbation_seeds[i])}
    labelled = data.label_mask if data.label_mask is not None else np.ones_like(kept, dtype=bool)
    counts = {"drugs": 0, "drugs_with_a_kept_positive": 0, "drugs_targeting_a_labelled_gene": 0, "drugs_with_a_target_gene_in_another_fold": 0,
              "kept_drug_positives": 0, "kept_drug_positives_shared_with_a_target_gene_in_another_fold": 0,
              "drug_symptom_pairs_whose_target_gene_is_positive": 0, "of_which_the_drug_is_positive": 0, "drug_pairs_labelled": 0}
    for i, (perturbation_id, kind) in enumerate(zip(data.perturbation_ids, data.perturbation_types)):
        if kind != "drug":
            continue
        counts["drugs"] += 1
        drug_positive = kept[i] > 0
        counts["drugs_with_a_kept_positive"] += int(drug_positive.any())
        counts["kept_drug_positives"] += int(drug_positive.sum())
        counts["drug_pairs_labelled"] += int(labelled[i].sum())
        target_genes = [gene_of_node[data.node_ids[s]] for s in data.perturbation_seeds[i] if data.node_ids[s] in gene_of_node]
        if not target_genes:
            continue
        counts["drugs_targeting_a_labelled_gene"] += 1
        in_other_fold = [gene for gene in target_genes if folds[gene] != folds[perturbation_id]]
        counts["drugs_with_a_target_gene_in_another_fold"] += int(bool(in_other_fold))
        target_positive = np.zeros(kept.shape[1], dtype=bool)
        for gene in target_genes:
            target_positive |= kept[gene_row[gene]] > 0
        other_fold_positive = np.zeros(kept.shape[1], dtype=bool)
        for gene in in_other_fold:
            other_fold_positive |= kept[gene_row[gene]] > 0
        counts["kept_drug_positives_shared_with_a_target_gene_in_another_fold"] += int((drug_positive & other_fold_positive).sum())
        counts["drug_symptom_pairs_whose_target_gene_is_positive"] += int((target_positive & labelled[i]).sum())
        counts["of_which_the_drug_is_positive"] += int((target_positive & labelled[i] & drug_positive).sum())
    group_sizes = pd.Series(data.group_ids).value_counts()
    fold_sizes = pd.Series([folds[p] for p in data.perturbation_ids]).value_counts().sort_index()
    counts.update(groups=int(len(group_sizes)), largest_group=int(group_sizes.iloc[0]), next_largest_groups=group_sizes.iloc[1:4].astype(int).tolist(),
                  fold_sizes=fold_sizes.astype(int).tolist(),
                  kept_positives_per_fold=[int(kept[[folds[p] == fold for p in data.perturbation_ids]].sum()) for fold in range(num_folds)])
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v1_full_v2.parquet"))
    parser.add_argument("--num-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/drug_target_leakage.md"))
    parser.add_argument("--json-output", type=Path, default=Path("runs/drug_target_leakage.json"))
    arguments = parser.parse_args()
    results = {}
    for grouping in GROUPINGS:
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by=grouping, label_selection=arguments.label_selection)
        results[grouping] = leakage_counts(data, arguments.num_folds, arguments.seed)
        print(grouping, json.dumps(results[grouping]))
    arguments.json_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.json_output.write_text(json.dumps(results, indent=1))
    rows = [("drugs", "drugs"), ("drugs with a kept positive", "drugs_with_a_kept_positive"), ("drugs targeting a labelled gene", "drugs_targeting_a_labelled_gene"),
            ("of these, with a target gene in another fold", "drugs_with_a_target_gene_in_another_fold"), ("kept drug positives", "kept_drug_positives"),
            ("kept drug positives that a target gene in another fold shares", "kept_drug_positives_shared_with_a_target_gene_in_another_fold"),
            ("leakage groups", "groups"), ("largest leakage group (perturbations)", "largest_group"), ("next largest groups", "next_largest_groups"),
            ("perturbations per fold", "fold_sizes"), ("kept positive pairs per fold", "kept_positives_per_fold")]
    first = results[GROUPINGS[0]]
    conditional = first["of_which_the_drug_is_positive"] / max(first["drug_symptom_pairs_whose_target_gene_is_positive"], 1)
    base_rate = first["kept_drug_positives"] / max(first["drug_pairs_labelled"], 1)
    lines = ["# Drug-target leakage in the grouped split (generated by experiments/check_drug_target_leakage.py)", "",
             f"Evidence {arguments.evidence_dir}, label selection {arguments.label_selection}, graph {arguments.graph_dir}; {arguments.num_folds} folds, seed {arguments.seed}. "
             "A drug's perturbation nodes are its targets; a target gene is labelled when its loss of function is a perturbation with labels. "
             "disease_cluster_and_targets holds a drug out with every labelled gene it targets and every drug sharing a target node "
             "(experiment_data.merge_drugs_with_their_targets).", "",
             "| count | " + " | ".join(GROUPINGS) + " |", "|---|" + "---|" * len(GROUPINGS)]
    for label, key in rows:
        lines.append(f"| {label} | " + " | ".join(str(results[grouping][key]) for grouping in GROUPINGS) + " |")
    lines += ["", f"How much a target gene's labels say about the drug: over the drug-symptom pairs where a labelled target gene of the drug is a kept positive, "
              f"the drug is a kept positive in {first['of_which_the_drug_is_positive']} of {first['drug_symptom_pairs_whose_target_gene_is_positive']} "
              f"({conditional:.2f}), against {base_rate:.2f} over all labelled drug-symptom pairs. Inhibiting a target and losing the gene need not "
              "give the same symptoms (a drug can be an agonist, and a gene loss is lifelong), so the share is a measure of what a label lookup "
              "across folds would earn, not of mechanism."]
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(f"wrote {arguments.markdown_output}")


if __name__ == "__main__":
    main()
