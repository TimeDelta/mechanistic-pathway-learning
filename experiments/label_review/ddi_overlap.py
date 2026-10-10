"""The curated drug-drug-event reference set (CRESCENDDI) against the study's drugs and symptoms.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/ddi_overlap.py data/raw/label_review/crescenddi/positive_controls.xlsx data/raw/label_review/crescenddi/negative_controls.xlsx
"""
import sys
import pandas as pd

positive = pd.read_excel(sys.argv[1])
negative = pd.read_excel(sys.argv[2])
print("positive controls:", len(positive), "columns:", list(positive.columns)[:10])
print("negative controls:", len(negative))
crosswalk = pd.read_csv("docs/symptom_crosswalk.csv")
symptom_of_term = {term.strip().lower(): row.target_symptom for row in crosswalk.itertuples() for term in str(row.meddra_preferred_terms).split(";") if term.strip() and term.strip() != "nan"}
def study_names(evidence_directory):
    records = pd.read_parquet(f"{evidence_directory}/evidence_records.parquet", columns=["perturbation_type", "perturbation_label"])
    return {name.lower() for name in records[records.perturbation_type == "drug"].perturbation_label.unique()}
name_sets = {"current rule (154 drugs)": study_names("data/processed/evidence_full_v3_parkinsonism"),
             "up to three targets (222)": study_names("data/processed/label_review/evidence_n_only_cap3"),
             "any ATC, up to three targets (924)": study_names("data/processed/label_review/evidence_any_atc_cap3")}
for table_name, table in (("positive", positive), ("negative", negative)):
    table = table.copy()
    table["first"], table["second"] = table.DRUG_1_CONCEPT_NAME.str.lower(), table.DRUG_2_CONCEPT_NAME.str.lower()
    table["symptom"] = table.EVENT_CONCEPT_NAME.str.lower().map(symptom_of_term)
    print(f"\n{table_name} controls: events {table.EVENT_CONCEPT_NAME.nunique()}, on a study symptom {int(table.symptom.notna().sum())} rows; events mapped: {sorted(table[table.symptom.notna()].EVENT_CONCEPT_NAME.unique())}")
    for label, names in name_sets.items():
        both = table[table["first"].isin(names) & table["second"].isin(names)]
        one = table[table["first"].isin(names) | table["second"].isin(names)]
        both_symptom = both[both.symptom.notna()]
        print(f"   {label}: both drugs in the study {len(both)} rows ({both[['first', 'second']].drop_duplicates().shape[0]} pairs); at least one {len(one)}; "
              f"both drugs and a study symptom {len(both_symptom)} rows ({both_symptom[['first', 'second']].drop_duplicates().shape[0]} pairs), by symptom {both_symptom.symptom.value_counts().to_dict()}")
    if table_name == "positive":
        print("   most frequent events overall:", table.EVENT_CONCEPT_NAME.value_counts().head(25).to_dict())
