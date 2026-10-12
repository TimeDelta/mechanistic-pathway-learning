"""Is there a measured cooperativity factor for the study's drugs, or only the binary alpha in {0, 1}?

Why this exists. For two drugs at one target, pure competition needs no free parameter: the Gaddum relation gives
the shift in occupancy as 1 + ([B]/K_B)/(1 + [A]/K_A), and the study has both affinities. For two drugs at
*different* sites on one target the allosteric ternary complex model applies, and it carries a cooperativity factor
alpha, which multiplies the orthosteric ligand's affinity when the modulator is bound: alpha below 1 is negative
cooperativity, above 1 positive, and 0 is the competitive limit. The off-target scoping document proposed alpha in
{0, 1}, which the owner rejected on 12 October 2026 and asked for "a better way that actually uses calculated or
measured values". This script measures whether such values exist.

Three sources are measured, each for the study's drugs and for the study's targets:
  1. ChEMBL activities whose standard_type is a cooperativity-style quantity (Cooperativity, Alpha, Fold shift,
     Shift, Dose ratio). Cached by --fetch; the cache is data/raw/label_review/chembl_cooperativity/.
  2. Guide to PHARMACOLOGY's "Receptor Site" column, a named binding site, which is the second site label beside
     ChEMBL's numeric site_id and crosses target records where site_id does not.
  3. Guide to PHARMACOLOGY's "Action comment" and "Assay Description", for cooperativity language.

Counting and quoting only: no label, table, graph or configuration of the study is changed.

Usage:
  python experiments/label_review/cooperativity_sources.py --fetch   # first run, writes the ChEMBL cache
  python experiments/label_review/cooperativity_sources.py
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.sider_chemical_matching import pubchem_cid_of_stitch_flat_id

CHEMBL_ACTIVITY_API = "https://www.ebi.ac.uk/chembl/api/data/activity"
# Every standard_type in ChEMBL that could carry a cooperativity factor or the dose ratio one is computed from.
# "Schild slope", "Fold-shift" and "alpha factor" were queried too and hold 0 activities in all of ChEMBL.
COOPERATIVITY_STANDARD_TYPES = ("Cooperativity", "Alpha", "Fold shift", "Shift", "Dose ratio")
CHEMBL_CACHE = Path("data/raw/label_review/chembl_cooperativity/cooperativity_activities.json")
GTOPDB_INTERACTIONS = Path("data/raw/gtopdb_2026/interactions.csv")
STATED_RANGE = re.compile(r"[Aa]ctivity range is (\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)")


def fetch_chembl_activities(cache_path: Path) -> None:
    """Every activity of the cooperativity-style standard types, cached whole so the counts are reproducible."""
    activities = []
    for standard_type in COOPERATIVITY_STANDARD_TYPES:
        offset = 0
        while True:
            query = urllib.parse.urlencode({"standard_type": standard_type, "format": "json", "limit": 1000, "offset": offset})
            request = urllib.request.Request(f"{CHEMBL_ACTIVITY_API}?{query}",
                                             headers={"Accept": "application/json", "User-Agent": "mechanistic-pathway-learning/0.1"})
            with urllib.request.urlopen(request, timeout=90) as response:
                page = json.load(response)
            for activity in page["activities"]:
                activity["_queried_standard_type"] = standard_type
                activities.append(activity)
            offset += 1000
            if offset >= page["page_meta"]["total_count"]:
                break
            time.sleep(1)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(activities, indent=1, sort_keys=True))


def study_drug_labels_by_chembl_id(evidence_directory: Path, chembl_directory: Path) -> dict[str, str]:
    """Every ChEMBL molecule id of a study drug -> the drug's label, through the PubChem mapping."""
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"][["perturbation_id", "perturbation_label"]].drop_duplicates()
    pubchem_to_chembl = json.loads((chembl_directory / "pubchem_to_chembl_with_parents.json").read_text())
    label_by_chembl_id = {}
    for drug in drugs.itertuples(index=False):
        for chembl_id in pubchem_to_chembl.get(str(pubchem_cid_of_stitch_flat_id(drug.perturbation_id)), []):
            label_by_chembl_id[chembl_id] = drug.perturbation_label
    return label_by_chembl_id


def stated_range_of(assay_description: str) -> tuple[float, float] | None:
    match = STATED_RANGE.search(str(assay_description))
    return (float(match.group(1)), float(match.group(2))) if match else None


def report_chembl(activities: pd.DataFrame, label_by_chembl_id: dict[str, str]) -> None:
    print(f"ChEMBL activities of a cooperativity-style type: {len(activities)} over "
          f"{activities.molecule_chembl_id.nunique()} molecules")
    print("  by type, target and organism:")
    for (standard_type, target_name, organism), group in activities.groupby(
            ["_queried_standard_type", "target_pref_name", "target_organism"]):
        print(f"    {standard_type}: {target_name} ({organism}), {len(group)} records")

    on_study_drugs = activities[activities.molecule_chembl_id.isin(label_by_chembl_id)]
    print(f"\n  on a study drug: {len(on_study_drugs)} records over "
          f"{on_study_drugs.molecule_chembl_id.map(label_by_chembl_id).nunique()} drugs")
    for row in on_study_drugs.itertuples(index=False):
        print(f"    {label_by_chembl_id[row.molecule_chembl_id]} | {row.target_pref_name} ({row.target_organism}) | "
              f"{row.standard_type} = {row.standard_value} {row.standard_units or ''}")
        print(f"        assay: {str(row.assay_description)[:200]}")

    cooperativity = activities[activities._queried_standard_type == "Cooperativity"]
    values = pd.to_numeric(cooperativity.standard_value, errors="coerce").dropna()
    print(f"\n  the Cooperativity records alone: {len(cooperativity)}, values from {values.min()} to {values.max()}, "
          f"median {values.median():.3f}")
    print(f"    above 1 (positive cooperativity): {int((values > 1).sum())}; "
          f"below 1 (negative): {int((values < 1).sum())}; at the competitive limit 0: {int((values == 0).sum())}")

    ranges = [stated_range_of(description) for description in cooperativity.assay_description]
    with_range = [(value, bounds) for value, bounds in zip(pd.to_numeric(cooperativity.standard_value, errors="coerce"), ranges)
                  if bounds is not None]
    inside = sum(1 for value, (low, high) in with_range if low <= value <= high)
    print(f"    records whose own assay description states a numeric range: {len(with_range)}; "
          f"value inside it: {inside}; outside it: {len(with_range) - inside}")


def report_gtopdb(interactions_path: Path) -> None:
    table = pd.read_csv(interactions_path, skiprows=1, low_memory=False)
    human = table[(table["Target Species"] == "Human") & table["Target Gene Symbol"].notna()]
    sited = human[human["Receptor Site"].notna()]
    print(f"\nGuide to PHARMACOLOGY, human rows: {len(human)}")
    print(f"  with a named Receptor Site: {len(sited)} over {sited['Ligand'].nunique()} ligands")
    for site, group in sited.groupby("Receptor Site"):
        plain_ligands = sorted({name.lower() for name in group["Ligand"] if "[" not in name and "&" not in name})
        print(f"    {site}: {len(group)} rows, {len(plain_ligands)} ligands that are not radioligands or coded compounds")
        print(f"      {plain_ligands}")

    with_cooperativity_language = human[human["Action comment"].fillna("").str.contains("cooperativ", case=False)]
    print(f"\n  rows whose Action comment says cooperativity: {len(with_cooperativity_language)}")
    for row in with_cooperativity_language.itertuples(index=False):
        print(f"    {row[human.columns.get_loc('Ligand')]} | {row[human.columns.get_loc('Target')]} | "
              f"\"{row[human.columns.get_loc('Action comment')]}\"")
    contexts = human["Ligand Context"].dropna().value_counts()
    print(f"  the Ligand Context column, which is topological and not an assay condition: {dict(contexts)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--fetch", action="store_true", help="refresh the ChEMBL cache before reporting")
    arguments = parser.parse_args()

    if arguments.fetch or not CHEMBL_CACHE.exists():
        fetch_chembl_activities(CHEMBL_CACHE)
    activities = pd.DataFrame(json.loads(CHEMBL_CACHE.read_text()))
    report_chembl(activities, study_drug_labels_by_chembl_id(arguments.evidence_dir, arguments.chembl_dir))
    report_gtopdb(GTOPDB_INTERACTIONS)


if __name__ == "__main__":
    main()
