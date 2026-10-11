"""How much of the Open Targets target-level curation the study's drug labels already say.

A curated triple is (target, direction, symptom): "inhibition of X causes symptom Y". A study drug says the same thing
when it acts on X in that direction (a mechanism edge of negative sign for inhibition, positive for activation) and
has Y as a kept positive label. The triples are those of experiments/scope_open_targets_symptom_labels.py, matched by
keyword and so an upper bound on a reviewed crosswalk.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/scope_open_targets_symptom_labels.py
  python experiments/label_review/open_targets_overlap.py
"""
import json
from collections import defaultdict

import pandas as pd

rows = pd.read_parquet("data/processed/off_target_scoping/open_targets_symptom_rows.parquet")
triples = rows.groupby(["target", "direction", "target_symptom"]).agg(
    animal_readout=("animal_readout", "all"), species_ambiguous=("species_ambiguous_readout", "all"),
    acute_only=("dosing", lambda values: set(values) == {"acute"}), references=("reference", lambda values: "; ".join(sorted(set(values))))).reset_index()
print("positive triples:", len(triples), "| by direction:", triples.direction.value_counts().to_dict(), "| targets:", triples.target.nunique(),
      "| rodent readout only:", int(triples.animal_readout.sum()), "| wording shared by rodent assay and clinic only:", int(triples.species_ambiguous.sum()),
      "| acute dosing only:", int(triples.acute_only.sum()))

for label, evidence_directory, selection_path in (
        ("current rule (154 drugs)", "data/processed/evidence_full_v3_parkinsonism", "data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"),
        ("up to three targets (222 drugs)", "data/processed/label_review/evidence_n_only_cap3", "data/processed/label_review/selection_n_only_cap3.parquet")):
    records = pd.read_parquet(f"{evidence_directory}/evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label", "perturbation_nodes"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")
    selection = pd.read_parquet(selection_path)
    selection = selection[selection.perturbation_type == "drug"]
    kept = defaultdict(set)
    listed = defaultdict(set)
    for row in selection.itertuples():
        listed[row.perturbation_id].add(row.symptom)
        if row.keep:
            kept[row.perturbation_id].add(row.symptom)
    drugs_on_target = defaultdict(list)  # (gene symbol, direction) -> drugs
    for row in drugs.itertuples():
        for node, sign, _ in json.loads(row.perturbation_nodes):
            if str(node).startswith("GENE:") and sign != 0:
                drugs_on_target[(str(node).split(":", 1)[1], "inhibition" if sign < 0 else "activation")].append(row.perturbation_id)
    genes_with_a_drug = {gene for gene, _ in drugs_on_target}
    said_by_a_kept_label = said_by_a_listed_label = target_has_drug_in_direction = target_has_any_drug = 0
    new_rows = []
    for triple in triples.itertuples():
        acting = set(drugs_on_target.get((triple.target, triple.direction), []))
        target_has_any_drug += triple.target in genes_with_a_drug
        target_has_drug_in_direction += bool(acting)
        if any(triple.target_symptom in kept[drug] for drug in acting):
            said_by_a_kept_label += 1
        elif any(triple.target_symptom in listed[drug] for drug in acting):
            said_by_a_listed_label += 1
        else:
            new_rows.append(triple)
    new = pd.DataFrame(new_rows)
    print(f"\n{label}: triples whose target a study drug acts on in the same direction: {target_has_drug_in_direction} (in either direction: {target_has_any_drug})")
    print(f"   already said by a kept drug label: {said_by_a_kept_label} | by a listed but masked drug label: {said_by_a_listed_label} | said by no drug label of the study: {len(new)}")
    print(f"   of those {len(new)}: target with no study drug in that direction {int((~new.apply(lambda row: (row.target, row.direction) in drugs_on_target, axis=1)).sum())}, "
          f"rodent readout only {int(new.animal_readout.sum())}, shared wording only {int(new.species_ambiguous.sum())}, acute dosing only {int(new.acute_only.sum())}, "
          f"activation {int((new.direction == 'activation').sum())}, inhibition {int((new.direction == 'inhibition').sum())}")
    human_new = new[~new.animal_readout & ~new.species_ambiguous]
    print(f"   not said by a drug label and not a rodent or shared-wording readout: {len(human_new)} over {human_new.target.nunique()} targets; by symptom {human_new.target_symptom.value_counts().to_dict()}")
    if label.startswith("up to three"):
        print("   targets of those:", sorted(human_new.target.unique()))
