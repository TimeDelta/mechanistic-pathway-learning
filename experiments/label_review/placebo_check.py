"""How many of the study's kept positive drug pairs have a placebo arm in SIDER, and in how many the drug arm's
frequency is no higher than the placebo arm's.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/placebo_check.py
"""
import pandas as pd

crosswalk = pd.read_csv("docs/symptom_crosswalk.csv")
symptom_of_term = {term.strip().lower(): row.target_symptom for row in crosswalk.itertuples() for term in str(row.meddra_preferred_terms).split(";") if term.strip() and term.strip() != "nan"}
frequency = pd.read_csv("data/raw/sider_4.1/meddra_freq.tsv.gz", sep="\t", header=None,
                        names=["stitch_flat", "stitch_stereo", "label_concept", "placebo", "frequency_text", "lower", "upper", "concept_type", "meddra_concept", "term"])
frequency = frequency[frequency.concept_type == "PT"].copy()
frequency["symptom"] = frequency.term.str.lower().map(symptom_of_term)
frequency = frequency.dropna(subset=["symptom"])
frequency["midpoint"] = (frequency.lower + frequency.upper) / 2
frequency["is_placebo"] = frequency.placebo.fillna("").astype(str).str.lower().eq("placebo")
selection = pd.read_parquet("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet")
drug_selection = selection[selection.perturbation_type == "drug"]
kept = drug_selection[drug_selection.keep.astype(bool)]
print("drug pairs in the selection:", len(drug_selection), "| kept:", len(kept), "| drugs:", drug_selection.perturbation_id.nunique())
study = frequency[frequency.stitch_flat.isin(set(drug_selection.perturbation_id))]
by_arm = study.groupby(["stitch_flat", "symptom", "is_placebo"]).midpoint.agg(["mean", "max", "count"]).unstack("is_placebo")
both = by_arm.dropna(subset=[("mean", True), ("mean", False)])
print("study drug-symptom pairs with a frequency row in SIDER:", len(by_arm), "| with a placebo-arm row as well:", len(both), "| drugs with one:", both.index.get_level_values(0).nunique())
kept_pairs = set(zip(kept.perturbation_id, kept.symptom))
all_pairs = set(zip(drug_selection.perturbation_id, drug_selection.symptom))
both = both.reset_index()
both.columns = ["perturbation_id", "symptom", "drug_mean", "placebo_mean", "drug_max", "placebo_max", "drug_rows", "placebo_rows"]
both["kept"] = [(p, s) in kept_pairs for p, s in zip(both.perturbation_id, both.symptom)]
both["in_selection"] = [(p, s) in all_pairs for p, s in zip(both.perturbation_id, both.symptom)]
both["not_above_placebo"] = both.drug_mean <= both.placebo_mean
both["ratio"] = both.drug_mean / both.placebo_mean.clip(lower=1e-9)
for label, subset in (("all pairs with both arms", both), ("kept positives with both arms", both[both.kept])):
    print(f"{label}: {len(subset)} | drug arm mean not above placebo mean: {int(subset.not_above_placebo.sum())} | drug mean under 1.5 times placebo: {int((subset.ratio < 1.5).sum())} | at least twice placebo: {int((subset.ratio >= 2).sum())}")
print(both[both.kept & both.not_above_placebo].sort_values("symptom").head(25).to_string(index=False))
print("kept positives by symptom with both arms and not above placebo:", both[both.kept & both.not_above_placebo].symptom.value_counts().to_dict())

# the same comparison inside one label at a time: a drug arm and a placebo arm of the same label come from the same trials
within = study.groupby(["stitch_flat", "symptom", "label_concept", "is_placebo"]).midpoint.mean().unstack("is_placebo").dropna()
within.columns = ["drug", "placebo"]
within = within.reset_index()
within["above"] = within.drug > within.placebo
pair_level = within.groupby(["stitch_flat", "symptom"]).agg(labels=("above", "size"), labels_above=("above", "sum"), largest_ratio=("drug", "max")).reset_index()
pair_level["kept"] = [(p, s) in kept_pairs for p, s in zip(pair_level.stitch_flat, pair_level.symptom)]
kept_level = pair_level[pair_level.kept]
print("\nwithin one label: label-level comparisons", len(within), "| pairs", len(pair_level), "| kept positive pairs", len(kept_level))
print("kept positive pairs in which no label has the drug arm above its own placebo arm:", int((kept_level.labels_above == 0).sum()),
      "| in which every label has it above:", int((kept_level.labels_above == kept_level.labels).sum()), "| mixed:", int(((kept_level.labels_above > 0) & (kept_level.labels_above < kept_level.labels)).sum()))
never = kept_level[kept_level.labels_above == 0]
print("by symptom, never above placebo:", never.symptom.value_counts().to_dict())
print("drugs involved:", never.stitch_flat.nunique(), "| kept positives of the study:", len(kept), "| share of kept positives with a within-label comparison:", round(len(kept_level) / len(kept), 3))
names = pd.read_parquet("data/processed/evidence_full_v3_parkinsonism/evidence_records.parquet", columns=["perturbation_id", "perturbation_label"]).drop_duplicates()
name_of = dict(zip(names.perturbation_id, names.perturbation_label))
print(sorted({(name_of.get(p, p), s) for p, s in zip(never.stitch_flat, never.symptom)})[:60])
