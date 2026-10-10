"""Which ATC groups the drugs that enter without the nervous-system rule belong to, and what they add.

Part of the label source review of 10 October 2026 (docs/label_source_review.md). Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/atc_breakdown.py    (after build_label_variants.sh)
"""
import json
from collections import defaultdict
from pathlib import Path
import pandas as pd

root = Path("data/processed/label_review")
bridge = json.loads(Path("data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json").read_text())
entries = bridge.get("ingredients", bridge) if isinstance(bridge, dict) else bridge
entries = list(entries.values()) if isinstance(entries, dict) else entries
print("bridge entry keys:", sorted(entries[0].keys())[:20])
atc_of = defaultdict(set)
for entry in entries:
    for key in ("perturbation_id", "stitch_flat_id", "rxcui"):
        if entry.get(key):
            identifier = entry[key] if key != "rxcui" else f"RXCUI:{entry[key]}"
            atc_of[str(identifier)].update(entry.get("atc_codes") or [])
sider = pd.read_csv("data/raw/sider_4.1/drug_atc.tsv", sep="\t", header=None, names=["stitch_flat", "atc"])
for row in sider.itertuples():
    atc_of[row.stitch_flat].add(row.atc)
group_names = {"A": "alimentary tract and metabolism", "B": "blood", "C": "cardiovascular", "D": "dermatological", "G": "genito-urinary and sex hormones", "H": "systemic hormones",
               "J": "anti-infectives", "L": "antineoplastic and immunomodulating", "M": "musculo-skeletal", "N": "nervous system", "P": "antiparasitic", "R": "respiratory", "S": "sensory organs", "V": "various"}
current = set(pd.read_parquet(root / "selection_n_only_cap1.parquet").query("perturbation_type == 'drug'").perturbation_id)
selection = pd.read_parquet(root / "selection_any_atc_cap1.parquet").query("perturbation_type == 'drug'")
records = pd.read_parquet(root / "evidence_any_atc_cap1/evidence_records.parquet", columns=["perturbation_id", "perturbation_label"]).drop_duplicates()
label_of = dict(zip(records.perturbation_id, records.perturbation_label))
entering = selection[~selection.perturbation_id.isin(current)]
kept = entering[entering.keep.astype(bool)]
rows = []
unknown = 0
for letter, name in group_names.items():
    drugs = {drug for drug in entering.perturbation_id.unique() if any(code.startswith(letter) for code in atc_of.get(drug, ()))}
    kept_of_group = kept[kept.perturbation_id.isin(drugs)]
    top = kept_of_group.perturbation_id.value_counts().head(8)
    rows.append((letter, name, len(drugs), len(kept_of_group), kept_of_group.symptom.value_counts().head(4).to_dict(), [(label_of.get(drug, drug), int(count)) for drug, count in top.items()]))
no_code = [drug for drug in entering.perturbation_id.unique() if not atc_of.get(drug)]
print("entering drugs:", entering.perturbation_id.nunique(), "| kept pairs:", len(kept), "| entering drugs with no ATC code found:", len(no_code))
for row in sorted(rows, key=lambda row: -row[3]):
    print(row)
