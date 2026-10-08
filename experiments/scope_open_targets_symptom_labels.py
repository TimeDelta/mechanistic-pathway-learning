"""Scope the Open Targets target-level safety curation as a source of symptom labels for gene perturbations: how many
of the 23 target symptoms it reaches, through how many targets, and with which sign (docs/off_target_scoping.md).

Source (downloaded into its own directory under data/raw, which git ignores): the Open Targets target safety
curation, adverse_effects.tsv (one row per curated effect: symptom, biologicalSystem, effect as direction and dosing,
efoId, ensemblId, target, pmid, ref) and secondary_pharmacology.json (the secondary pharmacology panel of
Brennan et al. 2024, JSON lines). Nothing is trained and no evidence table is changed.

The curation names a target, a direction (activation or inhibition) and a free-text effect, so a crosswalk is needed
to reach the study's symptom vocabulary. This script uses keyword patterns per symptom, prints the rows each pattern
catches and counts what a label set built this way would hold. Keyword matching is a scoping estimate, not the
crosswalk a label build would use: it reads English effect names, so it catches animal readouts ("decreased locomotor
activity") and misses what it has no word for, and a "decreased X" row is not a positive for X.

    python experiments/scope_open_targets_symptom_labels.py
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ADVERSE_EFFECTS = Path("data/raw/opentargets_curation/adverse_effects.tsv")
SECONDARY_PHARMACOLOGY = Path("data/raw/opentargets_curation/secondary_pharmacology.json")
# One pattern per target symptom of configs/evaluation.yaml. Written against the curation's own effect names.
SYMPTOM_PATTERNS = {
    "anxiety": r"\banxi",
    "apathy": r"\bapath",
    "catatonia": r"\bcataton|\bcatalep",
    "cognitive_impairment": r"cognit|\bmemory\b|\bamnes|\blearning\b|\bconfusion\b|\bdementia\b|\bdelirium\b",
    "compulsive_behavior": r"compulsi|stereotyp",
    "decreased_libido": r"\blibido\b|sexual dysfunction",
    "depressed_mood": r"\bdepress(ion|ed|ive)\b|dysphor",
    "disinhibition_or_impulsivity": r"disinhibit|impulsiv",
    "elevated_mood_or_mania": r"\bmania\b|\bmanic\b|\beuphor",
    "emotional_lability": r"emotional lability|mood swing",
    "fatigue": r"\bfatigue\b|\basthenia\b|\blethargy\b|\btired",
    "hyperactivity": r"hyperactiv|increased (locomotor|motor|spontaneous) activity|\bhyperkine",
    "increased_appetite": r"increased (appetite|food intake|body weight)|\bweight gain\b|hyperphag",
    "insomnia": r"\binsomnia\b|\bsleep disturb|decreased sleep",
    "irritability_or_aggression": r"irritab|aggress",
    "psychomotor_agitation": r"\bagitat|\brestless|akathisia",
    "psychomotor_retardation": r"decreased (locomotor|motor|spontaneous) activity|psychomotor retardation|hypokine|\bhypoactiv|\bbradykine",
    "psychosis": r"\bpsychos[ie]s\b|psychotic|hallucinat|\bdelusion|schizophren",
    "self_injury": r"self.injur|self.mutilat",
    "somnolence_or_hypersomnia": r"somnolen|sedat|\bdrowsi|hypersomn|\bsleepiness\b",
    "suicidality": r"suicid",
    "abnormal_dreams": r"\bdream|\bnightmare",
    "anhedonia": r"anhedon",
}
# A caught row whose text names the opposite change ("decreased anxiety"), both changes ("increased/decreased memory")
# or another thing ("respiratory depression") is not a positive for the target symptom; these patterns sort it out.
OPPOSITE_PATTERNS = {
    "anxiety": r"anxiolysis|anxiolytic|^decreased anxiety",
    "cognitive_impairment": r"^increased memory",
    "irritability_or_aggression": r"^decreased aggression",
}
AMBIGUOUS_PATTERN = r"increased/decreased|decreased/increased|^cognitive effects$"
NOT_THE_SYMPTOM_PATTERNS = {"depressed_mood": r"respiratory"}
# Readouts measured in rodents, not reported by patients.
ANIMAL_READOUT_PATTERN = r"locomotor|catalep|stereotyp"


def row_class(target_symptom: str, curated_symptom: str) -> str:
    text = curated_symptom.lower()
    if re.search(NOT_THE_SYMPTOM_PATTERNS.get(target_symptom, r"(?!)"), text):
        return "not the symptom"
    if re.search(AMBIGUOUS_PATTERN, text):
        return "both directions"
    if re.search(OPPOSITE_PATTERNS.get(target_symptom, r"(?!)"), text):
        return "opposite direction"
    return "positive"


def mapped_rows(adverse_effects: pd.DataFrame) -> pd.DataFrame:
    """One row per (target symptom, target, curated row) the patterns catch."""
    rows = []
    symptom_text = adverse_effects.symptom.astype(str).str.lower()
    for symptom, pattern in SYMPTOM_PATTERNS.items():
        caught = adverse_effects[symptom_text.str.contains(pattern, regex=True, na=False)]
        for record in caught.itertuples(index=False):
            direction, _, dosing = str(record.effect).partition("_")
            rows.append({"target_symptom": symptom, "curated_symptom": record.symptom, "row_class": row_class(symptom, str(record.symptom)),
                         "animal_readout": bool(re.search(ANIMAL_READOUT_PATTERN, str(record.symptom).lower())), "target": record.target,
                         "ensembl_id": record.ensemblId, "direction": direction, "dosing": dosing,
                         "biological_system": record.biologicalSystem, "reference": record.ref})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/off_target_scoping/open_targets_symptom_labels.json"))
    parser.add_argument("--rows-output", type=Path, default=Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet"))
    arguments = parser.parse_args()

    adverse_effects = pd.read_csv(ADVERSE_EFFECTS, sep="\t")
    secondary = [json.loads(line) for line in SECONDARY_PHARMACOLOGY.read_text().splitlines() if line.strip()]
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_type", "gene_symbol"])
    gene_symbols = set(nodes.gene_symbol[nodes.node_type == "gene"].dropna())

    rows = mapped_rows(adverse_effects)
    caught = rows.drop_duplicates(["target_symptom", "target", "direction", "row_class"])
    rows_by_class = {name: sorted(group.curated_symptom.str.lower().unique()) for name, group in rows[rows.row_class != "positive"].groupby("row_class")}
    rows = rows[rows.row_class == "positive"]
    pairs = rows.drop_duplicates(["target_symptom", "target", "direction"])
    both_directions = (pairs.groupby(["target_symptom", "target"]).direction.nunique() > 1).sum()
    summary = {
        "curation": {"rows": int(len(adverse_effects)), "targets": int(adverse_effects.target.nunique()),
                     "curated_symptoms": int(adverse_effects.symptom.nunique()),
                     "effects": adverse_effects.effect.value_counts().astype(int).to_dict(),
                     "references": adverse_effects.ref.value_counts().astype(int).to_dict()},
        "secondary_pharmacology": {"records": len(secondary), "events": len({record["event"] for record in secondary}),
                                   "targets": len({record["id"] for record in secondary}),
                                   "datasources": dict(Counter(record.get("datasource") for record in secondary))},
        "caught_by_class": caught.row_class.value_counts().astype(int).to_dict(), "curated_terms_set_aside": rows_by_class,
        "positive": {"rows": int(len(rows)), "symptom_target_direction_triples": int(len(pairs)),
                   "triples_from_animal_readouts": int(pairs.animal_readout.sum()),
                   "symptom_target_pairs": int(pairs.drop_duplicates(["target_symptom", "target"]).shape[0]),
                   "targets": int(pairs.target.nunique()), "target_symptoms_reached": int(pairs.target_symptom.nunique()),
                   "target_symptoms_total": len(SYMPTOM_PATTERNS),
                   "targets_that_are_graph_gene_nodes": int(pairs[pairs.target.isin(gene_symbols)].target.nunique()),
                   "pairs_curated_in_both_directions": int(both_directions),
                   "per_target_symptom": pairs.target_symptom.value_counts().astype(int).to_dict(),
                   "per_direction": pairs.direction.value_counts().astype(int).to_dict(),
                   "per_dosing": pairs.dosing.value_counts().astype(int).to_dict()},
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(summary, indent=2))
    rows.to_parquet(arguments.rows_output, index=False)
    print(json.dumps(summary, indent=2))
    print(f"wrote {arguments.output} and {arguments.rows_output}")


if __name__ == "__main__":
    main()
