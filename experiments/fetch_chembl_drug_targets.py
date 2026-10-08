"""Fetch and cache the ChEMBL side of evidence class E2 (design section 4.2).

Steps, each idempotent and cached under data/raw/chembl/:
  1. mechanisms.json                   every ChEMBL drug mechanism (molecule, target, action type)
  2. unichem_pubchem_to_chembl.json    SIDER PubChem CIDs -> ChEMBL molecule ids via UniChem
                                       (flat id first, stereo ids as fallback)
  3. molecule_parents.json             salt form -> parent molecule, because mechanisms are on parents
  4. pubchem_to_chembl_with_parents.json  the resolved mapping used downstream
  5. targets.json                      target records with gene symbols for every target hit

Usage:
  python experiments/fetch_chembl_drug_targets.py --sider-dir data/raw/sider_4.1 --chembl-dir data/raw/chembl
Rerun to resume; the time budget per invocation keeps it inside a cluster job's wall time.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import mechanisms_by_parent_and_molecule

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
UNICHEM_API = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNICHEM_PUBCHEM_SOURCE_ID = 22
USER_AGENT = "mechanistic-pathway-learning/0.1"


def get_json(url: str, retries: int = 3) -> dict:
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception:  # noqa: BLE001 - network retry
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"failed after retries: {url}")


def unichem_chembl_ids_for_pubchem_cid(pubchem_cid: str) -> list[str] | None:
    payload = json.dumps({"type": "sourceID", "compound": pubchem_cid, "sourceID": UNICHEM_PUBCHEM_SOURCE_ID}).encode()
    for attempt in range(3):
        try:
            request = urllib.request.Request(UNICHEM_API, data=payload, headers={"Content-Type": "application/json", "Accept": "application/json", "User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=30) as response:
                document = json.load(response)
            return sorted({source["compoundId"] for compound in document.get("compounds", []) for source in compound.get("sources", []) if source.get("shortName") == "chembl"})
        except Exception:  # noqa: BLE001 - network retry
            time.sleep(1 + attempt)
    return None


def load_or_empty(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def fetch_all_mechanisms(chembl_directory: Path) -> list[dict]:
    mechanisms_path = chembl_directory / "mechanisms.json"
    if mechanisms_path.exists():
        return json.loads(mechanisms_path.read_text())
    mechanisms: list[dict] = []
    offset = 0
    while True:
        page = get_json(f"{CHEMBL_API}/mechanism.json?limit=1000&offset={offset}")
        mechanisms.extend(page["mechanisms"])
        offset += 1000
        if offset >= page["page_meta"]["total_count"]:
            break
    mechanisms_path.write_text(json.dumps(mechanisms))
    return mechanisms


def sider_flat_to_stereo_ids(sider_directory: Path) -> dict[str, set[str]]:
    flat_to_stereo: dict[str, set[str]] = {}
    with gzip.open(sider_directory / "meddra_all_se.tsv.gz", "rt", encoding="utf-8") as side_effect_file:
        for row in csv.reader(side_effect_file, delimiter="\t"):
            flat_to_stereo.setdefault(row[0], set()).add(row[1])
    return flat_to_stereo


def map_sider_drugs_to_chembl(sider_directory: Path, chembl_directory: Path, time_budget_seconds: float) -> dict[str, list[str]]:
    mapping_path = chembl_directory / "unichem_pubchem_to_chembl.json"
    mapping = load_or_empty(mapping_path)
    flat_to_stereo = sider_flat_to_stereo_ids(sider_directory)
    started = time.time()
    jobs: list[tuple[str, str]] = []
    for flat_id, stereo_ids in flat_to_stereo.items():
        flat_cid = str(int(flat_id[4:]))
        if mapping.get(flat_cid):
            continue
        jobs.append((flat_cid, flat_cid))
        jobs.extend((flat_cid, str(int(stereo_id[4:]))) for stereo_id in stereo_ids)
    with ThreadPoolExecutor(max_workers=8) as pool:
        for (flat_cid, _), chembl_ids in zip(jobs, pool.map(lambda job: unichem_chembl_ids_for_pubchem_cid(job[1]), jobs)):
            if chembl_ids:
                mapping[flat_cid] = sorted(set(mapping.get(flat_cid, [])) | set(chembl_ids))
            elif chembl_ids is not None:
                mapping.setdefault(flat_cid, [])
            if time.time() - started > time_budget_seconds:
                break
    mapping_path.write_text(json.dumps(mapping))
    return mapping


def resolve_parent_molecules(mapping: dict[str, list[str]], mechanisms: list[dict], chembl_directory: Path) -> dict[str, list[str]]:
    parents_path = chembl_directory / "molecule_parents.json"
    parents = load_or_empty(parents_path)
    molecules_with_mechanisms = {mechanism["molecule_chembl_id"] for mechanism in mechanisms}
    to_check = sorted({chembl_id for chembl_ids in mapping.values() for chembl_id in chembl_ids if chembl_id not in molecules_with_mechanisms and chembl_id not in parents})
    for start in range(0, len(to_check), 50):
        batch = to_check[start : start + 50]
        page = get_json(f"{CHEMBL_API}/molecule.json?limit=50&only=molecule_chembl_id,molecule_hierarchy,pref_name,max_phase&molecule_chembl_id__in=" + ",".join(batch))
        for molecule in page.get("molecules", []):
            hierarchy = molecule.get("molecule_hierarchy") or {}
            parents[molecule["molecule_chembl_id"]] = {"parent": hierarchy.get("parent_chembl_id") or molecule["molecule_chembl_id"], "pref_name": molecule.get("pref_name"), "max_phase": molecule.get("max_phase")}
    parents_path.write_text(json.dumps(parents))
    resolved = {}
    for pubchem_cid, chembl_ids in mapping.items():
        expanded = set(chembl_ids)
        for chembl_id in chembl_ids:
            parent = parents.get(chembl_id, {}).get("parent")
            if parent:
                expanded.add(parent)
        resolved[pubchem_cid] = sorted(expanded)
    (chembl_directory / "pubchem_to_chembl_with_parents.json").write_text(json.dumps(resolved))
    return resolved


def fetch_targets(resolved: dict[str, list[str]], mechanisms: list[dict], chembl_directory: Path) -> dict[str, dict]:
    targets_path = chembl_directory / "targets.json"
    targets = load_or_empty(targets_path)
    mechanisms_by_molecule = mechanisms_by_parent_and_molecule(mechanisms)  # salt-form mechanisms count for their parent
    hit_targets = sorted({mechanism["target_chembl_id"] for chembl_ids in resolved.values() for chembl_id in chembl_ids for mechanism in mechanisms_by_molecule.get(chembl_id, []) if mechanism.get("target_chembl_id")})
    to_fetch = [target_id for target_id in hit_targets if target_id not in targets]
    for start in range(0, len(to_fetch), 50):
        batch = to_fetch[start : start + 50]
        page = get_json(f"{CHEMBL_API}/target.json?limit=50&only=target_chembl_id,pref_name,target_type,organism,target_components&target_chembl_id__in=" + ",".join(batch))
        for target in page.get("targets", []):
            gene_symbols = sorted({synonym["component_synonym"] for component in target.get("target_components", []) for synonym in component.get("target_component_synonyms", []) if synonym.get("syn_type") == "GENE_SYMBOL"})
            accessions = sorted({component.get("accession") for component in target.get("target_components", []) if component.get("accession")})
            targets[target["target_chembl_id"]] = {"pref_name": target.get("pref_name"), "target_type": target.get("target_type"), "organism": target.get("organism"), "gene_symbols": gene_symbols, "accessions": accessions}
    targets_path.write_text(json.dumps(targets))
    return targets


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"))
    parser.add_argument("--time-budget-seconds", type=float, default=240.0)
    arguments = parser.parse_args()
    arguments.chembl_dir.mkdir(parents=True, exist_ok=True)
    mechanisms = fetch_all_mechanisms(arguments.chembl_dir)
    mapping = map_sider_drugs_to_chembl(arguments.sider_dir, arguments.chembl_dir, arguments.time_budget_seconds)
    resolved = resolve_parent_molecules(mapping, mechanisms, arguments.chembl_dir)
    targets = fetch_targets(resolved, mechanisms, arguments.chembl_dir)
    status = get_json(f"{CHEMBL_API}/status.json")
    (arguments.chembl_dir / "chembl_release.json").write_text(json.dumps(status))
    print(f"mechanisms {len(mechanisms)}; SIDER drugs mapped {sum(1 for ids in mapping.values() if ids)}/{len(mapping)}; targets {len(targets)}; release {status.get('chembl_db_version')}")


if __name__ == "__main__":
    main()
