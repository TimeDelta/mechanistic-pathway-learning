"""Build data/processed/label_review/drug_names.csv: per study drug, the names it is known by, for the name checks.

The file was first written by a research pass on 10 October 2026 and had no script; the tests that read it
(tests/test_drug_name_resolution.py, tests/test_label_review_changes.py) skipped on any machine that lacked it. This
script derives it from sources the rebuild already holds, so a new container can make it:

- `perturbation_id`, `perturbation_label`: the drug rows of the evidence table (--evidence-dir).
- `rxnorm_names`: the RxNorm ingredient names the OnSIDES identifier bridge unified with the drug
  (data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json, experiments/fetch_onsides_identifier_bridge.py),
  semicolon-joined; empty when no ingredient unified with it.
- `pubchem_title`: PubChem's title for the CID (PUG REST, property Title), cached under data/raw/pubchem so a rerun
  makes no request. The tests say why this is not an authoritative generic name (often an IUPAC name or a brand).
- `in_current`, `in_approved`, `in_approved_with_list`: whether the drug has a kept row in each of three label
  selections (experiments/label_review/effect_site_coverage.py names them); False throughout when a selection file
  is absent, and the summary says which were present.

Usage:
  python experiments/label_review/build_drug_names.py
  python experiments/label_review/build_drug_names.py --no-fetch   # titles from the cache only; uncached CIDs stay empty
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mechanistic_pathway_learning.evidence.load_drug_label_events import pubchem_cid_from_stitch_flat  # noqa: E402

PUBCHEM_TITLE_ADDRESS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/Title/JSON"
DEFAULT_SELECTIONS = {"in_current": Path("data/processed/label_selection/better_v2_full_v4.parquet"),
                      "in_approved": Path("data/processed/label_review/selection_approved.parquet"),
                      "in_approved_with_list": Path("data/processed/label_review/selection_approved_with_list.parquet")}
BATCH_SIZE = 50


def study_drugs(evidence_directory: Path) -> pd.DataFrame:
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")[["perturbation_id", "perturbation_label"]]
    return drugs.sort_values("perturbation_id").reset_index(drop=True)


def rxnorm_names_by_perturbation(bridge_document: dict) -> dict[str, str]:
    """perturbation_id -> the RxNorm ingredient names the bridge unified with it, sorted and semicolon-joined."""
    names: dict[str, set[str]] = {}
    for record in bridge_document.get("ingredients", {}).values():
        if record.get("unified_with_sider") and record.get("perturbation_id"):
            names.setdefault(record["perturbation_id"], set()).add(str(record.get("ingredient_name", "")).strip())
    return {perturbation_id: ";".join(sorted(name for name in collected if name)) for perturbation_id, collected in names.items()}


def fetch_titles(cids: list[int], retries: int = 4) -> dict[int, str]:
    """PubChem titles for the CIDs, one request per batch; a CID PubChem does not know is left out."""
    titles: dict[int, str] = {}
    for start in range(0, len(cids), BATCH_SIZE):
        batch = cids[start:start + BATCH_SIZE]
        address = PUBCHEM_TITLE_ADDRESS.format(cids=",".join(str(cid) for cid in batch))
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(urllib.request.Request(address, headers={"User-Agent": "perturbation-symptom-modules/drug-names"}), timeout=60) as response:
                    document = json.loads(response.read().decode("utf-8"))
                break
            except (urllib.error.URLError, ConnectionResetError, TimeoutError) as error:
                if attempt == retries - 1:
                    raise
                print(f"retry {attempt + 1} after {error}", flush=True)
                time.sleep(2.0 * (attempt + 1))
        for item in document.get("PropertyTable", {}).get("Properties", []):
            if "Title" in item:
                titles[int(item["CID"])] = str(item["Title"])
        time.sleep(0.25)  # PubChem asks for no more than five requests a second
    return titles


def pubchem_titles(cids: list[int], cache_path: Path, fetch: bool = True) -> dict[int, str]:
    """Titles from the cache, fetching the missing ones when allowed and extending the cache in place."""
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
    titles = {int(cid): title for cid, title in cache.items()}
    missing = sorted(cid for cid in set(cids) if cid not in titles)
    if missing and fetch:
        titles.update(fetch_titles(missing))
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps({str(cid): title for cid, title in sorted(titles.items())}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return {cid: titles[cid] for cid in cids if cid in titles}


def kept_drugs(selection_path: Path) -> set[str] | None:
    """The perturbation ids with a kept row in a label selection; None when the file is absent."""
    if not selection_path.exists():
        return None
    selection = pd.read_parquet(selection_path)
    kept = selection[selection.keep.astype(bool)] if "keep" in selection.columns else selection
    return set(kept.perturbation_id.astype(str))


def build_table(drugs: pd.DataFrame, rxnorm_names: dict[str, str], titles: dict[int, str], selections: dict[str, set[str] | None]) -> pd.DataFrame:
    table = drugs.copy()
    table["rxnorm_names"] = table.perturbation_id.map(rxnorm_names).fillna("")
    table["pubchem_title"] = [titles.get(pubchem_cid_from_stitch_flat(perturbation_id), "") for perturbation_id in table.perturbation_id]
    for column, kept in selections.items():
        table[column] = [bool(kept is not None and perturbation_id in kept) for perturbation_id in table.perturbation_id]
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--bridge", type=Path, default=Path("data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json"))
    parser.add_argument("--pubchem-cache", type=Path, default=Path("data/raw/pubchem/compound_titles.json"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/drug_names.csv"))
    parser.add_argument("--no-fetch", action="store_true", help="titles from the cache only")
    for column, path in DEFAULT_SELECTIONS.items():
        parser.add_argument(f"--{column.replace('_', '-')}", type=Path, default=path)
    arguments = parser.parse_args()
    drugs = study_drugs(arguments.evidence_dir)
    rxnorm_names = rxnorm_names_by_perturbation(json.loads(arguments.bridge.read_text(encoding="utf-8"))) if arguments.bridge.exists() else {}
    titles = pubchem_titles([pubchem_cid_from_stitch_flat(perturbation_id) for perturbation_id in drugs.perturbation_id], arguments.pubchem_cache, fetch=not arguments.no_fetch)
    selections = {column: kept_drugs(getattr(arguments, column)) for column in DEFAULT_SELECTIONS}
    table = build_table(drugs, rxnorm_names, titles, selections)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output, index=False)
    summary = {"drugs": int(len(table)), "with_rxnorm_names": int((table.rxnorm_names != "").sum()), "with_pubchem_title": int((table.pubchem_title != "").sum()),
               "bridge_present": arguments.bridge.exists(), "selections_present": {column: kept is not None for column, kept in selections.items()},
               "kept_per_selection": {column: int(table[column].sum()) for column in selections}}
    print(json.dumps(summary, indent=1))
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
