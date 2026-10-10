"""Counts of HPO content the evidence table does not use: copy-number syndromes (DECIPHER entries), inheritance modes of
the diseases behind the study's gene labels, and the gain-of-function association types of Orphanet.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/hpo_unused.py
"""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    expand_symptom_terms, load_hpo_is_a_parents_from_obo, read_crosswalk_hpo_terms)

parents = load_hpo_is_a_parents_from_obo(Path("data/raw/hpo/hp.obo"))
crosswalk_terms = read_crosswalk_hpo_terms(Path("docs/symptom_crosswalk.csv"))
symptoms_of_term = expand_symptom_terms({symptom: roots for symptom, (roots, _) in crosswalk_terms.items()}, parents,
                                        {symptom: excluded for symptom, (_, excluded) in crosswalk_terms.items()})
hpoa = pd.read_csv("data/raw/hpo/phenotype.hpoa", sep="\t", comment="#", dtype=str).fillna("")
print("phenotype.hpoa rows:", len(hpoa), "| database prefixes:", hpoa.database_id.str.split(":").str[0].value_counts().to_dict())
decipher = hpoa[hpoa.database_id.str.startswith("DECIPHER")]
decipher_symptoms = decipher.assign(symptoms=decipher.hpo_id.map(lambda term: symptoms_of_term.get(term, [])))
decipher_symptoms = decipher_symptoms[decipher_symptoms.symptoms.map(len) > 0].explode("symptoms")
print("DECIPHER entries:", decipher.database_id.nunique(), "rows:", len(decipher), "| entries with a target-symptom term:", decipher_symptoms.database_id.nunique(),
      "| entry-symptom pairs:", len(decipher_symptoms.drop_duplicates(["database_id", "symptoms"])))
print(decipher_symptoms.drop_duplicates(["database_id", "symptoms"]).groupby("disease_name").symptoms.apply(lambda values: ", ".join(sorted(values))).to_string()[:3500])
print("by symptom:", decipher_symptoms.drop_duplicates(["database_id", "symptoms"]).symptoms.value_counts().to_dict())
