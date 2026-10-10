"""Counts of each evidence-table variant built for the label review.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/summarise_label_variants.py    (after build_label_variants.sh)
"""
import json, sys
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, ".")
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

root = Path("data/processed/label_review")
bridge = json.loads(Path("data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json").read_text())
reference_drugs = None
rows = []
for name in ("n_only_cap1", "any_atc_cap1", "n_only_cap3", "any_atc_cap3", "n_only_cap0", "any_atc_cap0"):
    evidence_directory, selection_path = root / f"evidence_{name}", root / f"selection_{name}.parquet"
    if not selection_path.exists():
        continue
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet")
    selection = pd.read_parquet(selection_path)
    drug_selection = selection[selection.perturbation_type == "drug"]
    kept = drug_selection[drug_selection.keep.astype(bool)]
    drug_records = records[records.perturbation_type == "drug"]
    label_of = dict(zip(drug_records.perturbation_id, drug_records.perturbation_label))
    drugs = set(drug_records.perturbation_id)
    data = load_experiment_data(Path("data/processed/graph_full_neuronal_split_binders"), evidence_directory, label_selection=selection_path, group_by="disease_cluster_and_targets")
    group_sizes = Counter(data.group_ids)
    is_drug = np.array(data.perturbation_types) == "drug"
    drug_groups = Counter(np.array(data.group_ids)[is_drug])
    largest_group, largest_size = group_sizes.most_common(1)[0]
    seeds_per_drug = [len(data.perturbation_seeds[index]) for index in np.flatnonzero(is_drug)]
    row = {"variant": name, "drugs": len(drugs), "drugs_loaded": int(is_drug.sum()), "perturbations_loaded": len(data.perturbation_ids), "drug_pairs": len(drug_selection), "kept_drug_pairs": len(kept),
           "drugs_with_a_kept_positive": kept.perturbation_id.nunique(), "kept_gene_pairs": int(selection[(selection.perturbation_type == "gene") & selection.keep.astype(bool)].shape[0]),
           "largest_group": largest_size, "drugs_in_the_largest_group": int(sum(1 for group in np.array(data.group_ids)[is_drug] if group == largest_group)),
           "leakage_groups": len(group_sizes), "groups_holding_a_drug": len(drug_groups), "median_target_nodes_per_drug": float(np.median(seeds_per_drug))}
    rows.append(row)
    print(json.dumps(row))
    print("   kept drug pairs by symptom:", kept.symptom.value_counts().to_dict())
    if name == "n_only_cap1":
        reference_drugs = drugs
        reference_kept = Counter(kept.perturbation_id)
    else:
        entering = sorted(drugs - reference_drugs, key=lambda drug: -Counter(kept.perturbation_id)[drug])
        kept_counts = Counter(kept.perturbation_id)
        print(f"   entering against n_only_cap1: {len(entering)} drugs, {sum(kept_counts[drug] for drug in entering)} kept pairs;",
              "most kept positives:", [(label_of.get(drug, drug), kept_counts[drug]) for drug in entering[:45]])
pd.DataFrame(rows).to_csv(root / "variant_counts.csv", index=False)
