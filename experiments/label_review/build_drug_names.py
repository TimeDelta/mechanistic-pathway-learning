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
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from mechanistic_pathway_learning.evidence.load_drug_label_events import pubchem_cid_from_stitch_flat  # noqa: E402

PUBCHEM_TITLE_ADDRESS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/Title/JSON"
PUBCHEM_NAME_TITLE_ADDRESS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/Title/JSON"
DEFAULT_SELECTIONS = {"in_current": Path("data/processed/label_selection/better_v2_full_v4.parquet"),
                      "in_approved": Path("data/processed/label_review/selection_approved.parquet"),
                      "in_approved_with_list": Path("data/processed/label_review/selection_approved_with_list.parquet")}
BATCH_SIZE = 50


def study_drugs(evidence_directory: Path) -> pd.DataFrame:
    records = pd.read_parquet(evidence_directory / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = records[records.perturbation_type == "drug"].drop_duplicates("perturbation_id")[["perturbation_id", "perturbation_label"]]
    return drugs.sort_values("perturbation_id").reset_index(drop=True)


def rxnorm_names_by_perturbation(bridge_document: dict) -> dict[str, str]:
    """perturbation_id -> the RxNorm ingredient names the bridge unified with it (a SIDER drug), or the ingredient's own
    name for a drug that enters under its RxCUI (an OnSIDES-only drug, id "RXCUI:<rxcui>"); sorted, semicolon-joined."""
    names: dict[str, set[str]] = {}
    for ingredient_id, record in bridge_document.get("ingredients", {}).items():
        name = str(record.get("ingredient_name", "")).strip()
        if record.get("unified_with_sider") and record.get("perturbation_id"):
            names.setdefault(record["perturbation_id"], set()).add(name)
        if record.get("identifier_system") == "RXCUI":
            names.setdefault(f"RXCUI:{ingredient_id}", set()).add(name)
    return {perturbation_id: ";".join(sorted(name for name in collected if name)) for perturbation_id, collected in names.items()}


def pubchem_key(perturbation_id: str, perturbation_label: str) -> tuple[str, str]:
    """How PubChem is asked for the drug's title: ("cid", "<cid>") for a SIDER identifier, ("name", "<label>") for an RxCUI."""
    if perturbation_id.startswith("CID1"):
        return "cid", str(pubchem_cid_from_stitch_flat(perturbation_id))
    return "name", str(perturbation_label).strip().lower()


def read_json(address: str, retries: int = 4) -> dict:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(address, headers={"User-Agent": "perturbation-symptom-modules/drug-names"}), timeout=60) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code == 404:  # PubChem answers 404 for a name or CID it does not know
                return {}
            if attempt == retries - 1:
                raise
        except (urllib.error.URLError, ConnectionResetError, TimeoutError) as error:
            if attempt == retries - 1:
                raise
            print(f"retry {attempt + 1} after {error}", flush=True)
        time.sleep(2.0 * (attempt + 1))
    return {}


def fetch_titles(keys: list[tuple[str, str]]) -> dict[tuple[str, str], str]:
    """PubChem titles: the CIDs in batches of BATCH_SIZE, the names one request each; a key PubChem does not know is
    left out. PubChem asks for no more than five requests a second."""
    titles: dict[tuple[str, str], str] = {}
    cids = sorted({value for kind, value in keys if kind == "cid"}, key=int)
    for start in range(0, len(cids), BATCH_SIZE):
        batch = cids[start:start + BATCH_SIZE]
        document = read_json(PUBCHEM_TITLE_ADDRESS.format(cids=",".join(batch)))
        for item in document.get("PropertyTable", {}).get("Properties", []):
            if "Title" in item:
                titles[("cid", str(item["CID"]))] = str(item["Title"])
        time.sleep(0.25)
    for name in sorted({value for kind, value in keys if kind == "name"}):
        document = read_json(PUBCHEM_NAME_TITLE_ADDRESS.format(name=urllib.parse.quote(name)))
        properties = document.get("PropertyTable", {}).get("Properties", [])
        if properties and "Title" in properties[0]:
            titles[("name", name)] = str(properties[0]["Title"])
        time.sleep(0.25)
    return titles


def pubchem_titles(keys: list[tuple[str, str]], cache_path: Path, fetch: bool = True) -> dict[tuple[str, str], str]:
    """Titles from the cache (keyed "cid:<cid>" and "name:<name>"), fetching the missing ones when allowed and extending
    the cache in place."""
    cache = json.loads(cache_path.read_text(encoding="utf-8")) if cache_path.exists() else {}
    titles = {tuple(key.split(":", 1)): title for key, title in cache.items()}
    missing = sorted(set(keys) - set(titles))
    if missing and fetch:
        titles.update(fetch_titles(missing))
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps({f"{kind}:{value}": title for (kind, value), title in sorted(titles.items())}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return {key: titles[key] for key in keys if key in titles}


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
    table["pubchem_title"] = [titles.get(pubchem_key(perturbation_id, label), "") for perturbation_id, label in zip(table.perturbation_id, table.perturbation_label)]
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
    titles = pubchem_titles([pubchem_key(perturbation_id, label) for perturbation_id, label in zip(drugs.perturbation_id, drugs.perturbation_label)],
                            arguments.pubchem_cache, fetch=not arguments.no_fetch)
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
