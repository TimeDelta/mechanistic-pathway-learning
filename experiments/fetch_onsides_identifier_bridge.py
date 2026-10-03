"""Build the identifier bridge for the OnSIDES later label slice (design section 4.2, evidence class E2).

Output: data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json, one record per RxNorm ingredient of
vocab_rxnorm_ingredient (1,866 in v3.1.1) with the route that resolved it to a ChEMBL parent molecule, the
identifiers seen on the way, its ATC classes, its unification with a SIDER drug and the version 1
pharmacological bar (single mechanism target, ATC N code). The decisions live in
mechanistic_pathway_learning/evidence/onsides_identifier_bridge.py; this script only fetches and joins.

Routes, in precedence order (the module's docstring has the reasons):
  (a) RxCUI -> UNII                 https://rxnav.nlm.nih.gov/REST/rxcui/<rxcui>/property.json?propName=UNII_CODE
      UNII  -> ChEMBL ids           POST https://www.ebi.ac.uk/unichem/api/v1/compounds {"type": "sourceID", "sourceID": 14}
      ChEMBL id -> parent           https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in=...
  (b) name -> ChEMBL pref_name      https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact=<name>
  ATC classes per RxCUI             https://rxnav.nlm.nih.gov/REST/rxclass/class/byRxcui.json?rxcui=<rxcui>&relaSource=ATC
Mechanism targets of parents missing from data/raw/chembl/targets.json are fetched with
fetch_chembl_drug_targets.fetch_targets, which extends that cache in place in its own format.

Every response is cached under data/raw/onsides/v3.1.1/http_cache/<host>/ by the shared
CachedRateLimitedClient (at most three requests per second per host, retries with backoff), so the run is
resumable and a rerun makes no request. Failures are recorded in the output under "errors" rather than
aborting the run.

Usage:
  OMP_NUM_THREADS=1 python experiments/fetch_onsides_identifier_bridge.py
  OMP_NUM_THREADS=1 python experiments/fetch_onsides_identifier_bridge.py --limit 20   # smoke run on the first ingredients
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))
sys.path.insert(0, str(REPOSITORY_ROOT / "experiments"))

from fetch_chembl_drug_targets import CHEMBL_API, fetch_all_mechanisms, fetch_targets  # noqa: E402
from mechanistic_pathway_learning.evidence.cached_http import USER_AGENT, CachedRateLimitedClient  # noqa: E402
from mechanistic_pathway_learning.evidence.load_drug_label_events import read_two_column_tsv  # noqa: E402
from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import (  # noqa: E402
    IDENTIFIER_SYSTEM_RXCUI,
    IngredientBridge,
    apply_bridge_precedence,
    bridge_counts,
    bridge_to_json_record,
    identifier_system_for_ingredient,
    mechanism_targets_for_parent,
    passes_nervous_system_rule,
    passes_single_target_rule,
    resolve_parent_from_name_route,
    resolve_parent_from_unii_route,
    rxnav_atc_classes_for_ingredient,
    rxnav_unii_codes,
    sider_parent_molecules,
    sider_pubchem_cids_by_parent,
    stitch_flat_id_from_pubchem_cid,
    unichem_chembl_ids,
    unify_with_sider,
)

RXNAV_API = "https://rxnav.nlm.nih.gov/REST"
UNICHEM_V1_COMPOUNDS = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
UNICHEM_LEGACY_ROUTE = "https://www.ebi.ac.uk/unichem/rest/src_compound_id/{compound}/{source_id}"
UNICHEM_FDA_SRS_SOURCE_ID = 14  # UniChem source "fdasrs": the FDA Substance Registration System, whose ids are UNIIs
MAX_REQUESTS_PER_SECOND = 3.0
HIERARCHY_BATCH_SIZE = 50
ONSIDES_RELEASE = "OnSIDES v3.1.1 (data release 3.1.1-20260422)"


class JsonPostCachedClient(CachedRateLimitedClient):
    """The shared client with a JSON content type on POST requests, which UniChem v1 requires (415 without it)."""

    def _request_with_backoff(self, url: str, post_body: bytes | None, accept: str) -> str:
        last_error: Exception | None = None
        for attempt in range(self.retries):
            self._wait_for_rate_limit()
            headers = {"Accept": accept, "User-Agent": USER_AGENT}
            if post_body is not None:
                headers["Content-Type"] = "application/json"
            request = urllib.request.Request(url, data=post_body, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    self.request_count += 1
                    return response.read().decode("utf-8")
            except urllib.error.HTTPError as error:
                last_error = error
                if error.code in (400, 404) and attempt == 0:
                    raise
            except Exception as error:  # noqa: BLE001 - network retry
                last_error = error
            time.sleep(min(60.0, 2.0 ** attempt))
        raise RuntimeError(f"request failed after {self.retries} attempts: {url}") from last_error


def load_ingredients(onsides_directory: Path, limit: int | None) -> pd.DataFrame:
    ingredients = pd.read_parquet(onsides_directory / "parquet" / "vocab_rxnorm_ingredient.parquet", columns=["rxnorm_id", "rxnorm_name"])
    ingredients = ingredients.drop_duplicates("rxnorm_id").sort_values("rxnorm_id").reset_index(drop=True)
    return ingredients.head(limit) if limit else ingredients


def probe_unichem_routes(unichem_client: JsonPostCachedClient, unii: str) -> dict:
    """Which UniChem REST form answers a UNII lookup: the legacy src_compound_id route and the v1 compounds POST, recorded for the pin."""
    outcomes: dict[str, dict] = {}
    legacy_url = UNICHEM_LEGACY_ROUTE.format(compound=unii, source_id=UNICHEM_FDA_SRS_SOURCE_ID)
    try:
        legacy_body = unichem_client.fetch_text(legacy_url)
        try:
            parsed = json.loads(legacy_body)
            outcomes["legacy_src_compound_id"] = {"url": legacy_url, "worked": bool(parsed), "note": "JSON body"}
        except json.JSONDecodeError:
            outcomes["legacy_src_compound_id"] = {"url": legacy_url, "worked": False, "note": "non-JSON body (HTML page)"}
    except urllib.error.HTTPError as error:
        outcomes["legacy_src_compound_id"] = {"url": legacy_url, "worked": False, "note": f"HTTP {error.code}"}
    except Exception as error:  # noqa: BLE001 - recorded, not fatal
        outcomes["legacy_src_compound_id"] = {"url": legacy_url, "worked": False, "note": str(error)[:200]}
    post_body = unichem_post_body(unii)
    try:
        document = unichem_client.fetch_json(UNICHEM_V1_COMPOUNDS, post_body)
        outcomes["v1_compounds_post"] = {"url": UNICHEM_V1_COMPOUNDS, "body": post_body.decode(), "worked": document.get("response") == "Success", "note": f"response {document.get('response')!r}, {len(document.get('compounds', []))} compounds"}
    except Exception as error:  # noqa: BLE001 - recorded, not fatal
        outcomes["v1_compounds_post"] = {"url": UNICHEM_V1_COMPOUNDS, "worked": False, "note": str(error)[:200]}
    outcomes["route_used"] = "v1_compounds_post" if outcomes["v1_compounds_post"]["worked"] else "none"
    return outcomes


def unichem_post_body(unii: str) -> bytes:
    return json.dumps({"type": "sourceID", "compound": unii, "sourceID": UNICHEM_FDA_SRS_SOURCE_ID}).encode("utf-8")


def fetch_unii_codes(rxnav_client: CachedRateLimitedClient, rxcuis: list[str], errors: list[dict]) -> dict[str, list[str]]:
    unii_by_rxcui: dict[str, list[str]] = {}
    for position, rxcui in enumerate(rxcuis, start=1):
        url = f"{RXNAV_API}/rxcui/{rxcui}/property.json?propName=UNII_CODE"
        try:
            unii_by_rxcui[rxcui] = rxnav_unii_codes(rxnav_client.fetch_json(url))
        except Exception as error:  # noqa: BLE001 - recorded, a rerun retries
            errors.append({"phase": "rxnav_unii", "ingredient_id": rxcui, "error": str(error)[:200]})
            unii_by_rxcui[rxcui] = []
        if position % 200 == 0:
            print(f"  RxNav UNII {position}/{len(rxcuis)} (requests {rxnav_client.request_count}, cache hits {rxnav_client.cache_hit_count})", flush=True)
    return unii_by_rxcui


def fetch_atc_classes(rxnav_client: CachedRateLimitedClient, rxcuis: list[str], errors: list[dict]) -> dict[str, list[str]]:
    atc_by_rxcui: dict[str, list[str]] = {}
    for position, rxcui in enumerate(rxcuis, start=1):
        url = f"{RXNAV_API}/rxclass/class/byRxcui.json?rxcui={rxcui}&relaSource=ATC"
        try:
            atc_by_rxcui[rxcui] = rxnav_atc_classes_for_ingredient(rxcui, rxnav_client.fetch_json(url))
        except Exception as error:  # noqa: BLE001 - recorded, a rerun retries
            errors.append({"phase": "rxnav_atc", "ingredient_id": rxcui, "error": str(error)[:200]})
            atc_by_rxcui[rxcui] = []
        if position % 200 == 0:
            print(f"  RxNav ATC {position}/{len(rxcuis)} (requests {rxnav_client.request_count}, cache hits {rxnav_client.cache_hit_count})", flush=True)
    return atc_by_rxcui


def fetch_unichem_chembl_ids(unichem_client: JsonPostCachedClient, uniis: list[str], errors: list[dict]) -> dict[str, list[str]]:
    chembl_ids_by_unii: dict[str, list[str]] = {}
    for position, unii in enumerate(uniis, start=1):
        try:
            chembl_ids_by_unii[unii] = unichem_chembl_ids(unichem_client.fetch_json(UNICHEM_V1_COMPOUNDS, unichem_post_body(unii)))
        except Exception as error:  # noqa: BLE001 - recorded, a rerun retries
            errors.append({"phase": "unichem", "unii": unii, "error": str(error)[:200]})
            chembl_ids_by_unii[unii] = []
        if position % 200 == 0:
            print(f"  UniChem {position}/{len(uniis)} (requests {unichem_client.request_count}, cache hits {unichem_client.cache_hit_count})", flush=True)
    return chembl_ids_by_unii


def fetch_molecule_hierarchy(chembl_client: CachedRateLimitedClient, chembl_ids: list[str], errors: list[dict]) -> dict[str, dict]:
    """ChEMBL id -> {"parent": parent id, "pref_name": name} in batches of 50 sorted ids, so the batch URLs are stable across reruns."""
    hierarchy: dict[str, dict] = {}
    ordered = sorted(set(chembl_ids), key=lambda chembl_id: int(chembl_id[6:]))
    for start in range(0, len(ordered), HIERARCHY_BATCH_SIZE):
        batch = ordered[start : start + HIERARCHY_BATCH_SIZE]
        url = f"{CHEMBL_API}/molecule.json?limit={HIERARCHY_BATCH_SIZE}&only=molecule_chembl_id,molecule_hierarchy,pref_name&molecule_chembl_id__in=" + ",".join(batch)
        try:
            page = chembl_client.fetch_json(url)
        except Exception as error:  # noqa: BLE001 - recorded, a rerun retries
            errors.append({"phase": "chembl_hierarchy", "chembl_ids": batch, "error": str(error)[:200]})
            continue
        for molecule in page.get("molecules", []):
            molecule_hierarchy = molecule.get("molecule_hierarchy") or {}
            hierarchy[molecule["molecule_chembl_id"]] = {"parent": molecule_hierarchy.get("parent_chembl_id") or molecule["molecule_chembl_id"], "pref_name": molecule.get("pref_name")}
    return hierarchy


def fetch_name_matches(chembl_client: CachedRateLimitedClient, names_by_ingredient: dict[str, str], errors: list[dict]) -> dict[str, list[dict]]:
    """ingredient id -> ChEMBL molecule records whose pref_name equals the ingredient name ignoring case."""
    molecules_by_ingredient: dict[str, list[dict]] = {}
    for position, (ingredient_id, name) in enumerate(names_by_ingredient.items(), start=1):
        url = f"{CHEMBL_API}/molecule.json?pref_name__iexact={urllib.parse.quote(name.strip())}&only=molecule_chembl_id,pref_name,molecule_hierarchy"
        try:
            molecules_by_ingredient[ingredient_id] = chembl_client.fetch_json(url).get("molecules", [])
        except Exception as error:  # noqa: BLE001 - recorded, a rerun retries
            errors.append({"phase": "chembl_name", "ingredient_id": ingredient_id, "error": str(error)[:200]})
            molecules_by_ingredient[ingredient_id] = []
        if position % 200 == 0:
            print(f"  ChEMBL names {position}/{len(names_by_ingredient)} (requests {chembl_client.request_count}, cache hits {chembl_client.cache_hit_count})", flush=True)
    return molecules_by_ingredient


def sider_name_bridge_comparison(ingredients: pd.DataFrame, bridges: dict[str, IngredientBridge], sider_drug_names: dict[str, str]) -> dict:
    """How the earlier name bridge (ingredient name equals a SIDER drug name, case-insensitive) relates to this bridge's unification."""
    cids_by_lowercase_name: dict[str, set[int]] = {}
    for stitch_flat_id, drug_name in sider_drug_names.items():
        cids_by_lowercase_name.setdefault(drug_name.strip().lower(), set()).add(int(stitch_flat_id[4:]))
    comparison = {"name_matched": 0, "name_matched_cids": 0, "name_matched_and_unified_same_cid": 0, "name_matched_and_unified_other_cid": 0, "name_matched_not_unified": 0, "name_matched_not_unified_without_parent": 0, "unified_without_name_match": 0, "neither": 0}
    matched_cids: set[int] = set()
    for ingredient_id, name in zip(ingredients.rxnorm_id, ingredients.rxnorm_name):
        bridge = bridges[ingredient_id]
        name_cids = cids_by_lowercase_name.get(str(name).strip().lower(), set())
        if name_cids:
            comparison["name_matched"] += 1
            matched_cids |= name_cids
            if bridge.unified_with_sider and bridge.sider_pubchem_cid in name_cids:
                comparison["name_matched_and_unified_same_cid"] += 1
            elif bridge.unified_with_sider:
                comparison["name_matched_and_unified_other_cid"] += 1
            else:
                comparison["name_matched_not_unified"] += 1
                if not bridge.chembl_parent:
                    comparison["name_matched_not_unified_without_parent"] += 1
        elif bridge.unified_with_sider:
            comparison["unified_without_name_match"] += 1
        else:
            comparison["neither"] += 1
    comparison["name_matched_cids"] = len(matched_cids)
    return comparison


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--onsides-dir", type=Path, default=REPOSITORY_ROOT / "data/raw/onsides/v3.1.1")
    parser.add_argument("--chembl-dir", type=Path, default=REPOSITORY_ROOT / "data/raw/chembl")
    parser.add_argument("--sider-dir", type=Path, default=REPOSITORY_ROOT / "data/raw/sider_4.1")
    parser.add_argument("--output", type=Path, default=None, help="defaults to <onsides-dir>/ingredient_identifier_bridge.json")
    parser.add_argument("--limit", type=int, default=None, help="only the first N ingredients by id, for a smoke run")
    parser.add_argument("--skip-target-fetch", action="store_true", help="do not extend data/raw/chembl/targets.json (offline check)")
    arguments = parser.parse_args()
    output_path = arguments.output or arguments.onsides_dir / "ingredient_identifier_bridge.json"
    started = time.time()
    errors: list[dict] = []

    cache_root = arguments.onsides_dir / "http_cache"
    rxnav_client = CachedRateLimitedClient(cache_root / "rxnav", MAX_REQUESTS_PER_SECOND)
    unichem_client = JsonPostCachedClient(cache_root / "unichem", MAX_REQUESTS_PER_SECOND)
    chembl_client = CachedRateLimitedClient(cache_root / "chembl", MAX_REQUESTS_PER_SECOND)

    ingredients = load_ingredients(arguments.onsides_dir, arguments.limit)
    identifier_systems = {ingredient_id: identifier_system_for_ingredient(ingredient_id) for ingredient_id in ingredients.rxnorm_id}
    rxcuis = [ingredient_id for ingredient_id, system in identifier_systems.items() if system == IDENTIFIER_SYSTEM_RXCUI]
    print(f"{len(ingredients)} ingredients: {len(rxcuis)} RxCUIs, {len(ingredients) - len(rxcuis)} OMOP extension ids", flush=True)

    # The ATC phase shares RxNav with the UNII phase, so it runs after it on the same limiter, in a thread next to the UniChem and ChEMBL phases.
    print("phase 1: RxNav UNII codes", flush=True)
    unii_by_rxcui = fetch_unii_codes(rxnav_client, rxcuis, errors)
    uniis = sorted({unii for unii_codes in unii_by_rxcui.values() for unii in unii_codes})
    print(f"  {sum(1 for codes in unii_by_rxcui.values() if codes)} RxCUIs with a UNII, {len(uniis)} distinct UNIIs", flush=True)

    with ThreadPoolExecutor(max_workers=1) as atc_pool:
        print("phase 2: RxNav ATC classes (in a thread) and UniChem FDA SRS lookups", flush=True)
        atc_future = atc_pool.submit(fetch_atc_classes, rxnav_client, rxcuis, errors)
        unichem_routes = probe_unichem_routes(unichem_client, uniis[0]) if uniis else {"route_used": "none"}
        print(f"  UniChem route: {unichem_routes['route_used']}", flush=True)
        chembl_ids_by_unii = fetch_unichem_chembl_ids(unichem_client, uniis, errors)
        chembl_ids_from_unichem = {rxcui: sorted({chembl_id for unii in unii_codes for chembl_id in chembl_ids_by_unii.get(unii, [])}, key=lambda chembl_id: int(chembl_id[6:])) for rxcui, unii_codes in unii_by_rxcui.items()}
        print(f"  {sum(1 for ids in chembl_ids_from_unichem.values() if ids)} RxCUIs with a ChEMBL id through UniChem", flush=True)

        print("phase 3: ChEMBL molecule hierarchy for the UniChem hits", flush=True)
        hierarchy = fetch_molecule_hierarchy(chembl_client, [chembl_id for ids in chembl_ids_from_unichem.values() for chembl_id in ids], errors)
        parent_by_chembl_id = {chembl_id: record["parent"] for chembl_id, record in hierarchy.items()}
        unii_routes = {ingredient_id: resolve_parent_from_unii_route(chembl_ids_from_unichem.get(ingredient_id, []), parent_by_chembl_id) for ingredient_id in ingredients.rxnorm_id}

        print("phase 4: ChEMBL pref_name matches for ingredients without a UniChem parent", flush=True)
        names_for_name_route = {ingredient_id: name for ingredient_id, name in zip(ingredients.rxnorm_id, ingredients.rxnorm_name) if unii_routes[ingredient_id][0] is None}
        molecules_by_ingredient = fetch_name_matches(chembl_client, names_for_name_route, errors)
        name_routes = {ingredient_id: resolve_parent_from_name_route(molecules_by_ingredient.get(ingredient_id, [])) for ingredient_id in ingredients.rxnorm_id}
        atc_by_rxcui = atc_future.result()

    print("phase 5: SIDER unification, ATC codes and mechanism targets", flush=True)
    mechanisms = fetch_all_mechanisms(arguments.chembl_dir)
    mechanisms_by_molecule: dict[str, list[dict]] = {}
    for mechanism in mechanisms:
        mechanisms_by_molecule.setdefault(mechanism["molecule_chembl_id"], []).append(mechanism)
    pubchem_to_chembl_with_parents = json.loads((arguments.chembl_dir / "pubchem_to_chembl_with_parents.json").read_text())
    molecule_parents = json.loads((arguments.chembl_dir / "molecule_parents.json").read_text()) if (arguments.chembl_dir / "molecule_parents.json").exists() else {}
    cids_by_parent = sider_pubchem_cids_by_parent(sider_parent_molecules(pubchem_to_chembl_with_parents, molecule_parents, mechanisms))
    sider_drug_names = {stitch_flat_id: names[0] for stitch_flat_id, names in read_two_column_tsv(arguments.sider_dir / "drug_names.tsv").items()}
    sider_atc_codes = read_two_column_tsv(arguments.sider_dir / "drug_atc.tsv")

    bridges: dict[str, IngredientBridge] = {}
    for ingredient_id, name in zip(ingredients.rxnorm_id, ingredients.rxnorm_name):
        bridge_method, chembl_parent, alternatives, pref_name = apply_bridge_precedence(identifier_systems[ingredient_id], unii_routes[ingredient_id], name_routes[ingredient_id])
        chembl_ids_seen = set(chembl_ids_from_unichem.get(ingredient_id, [])) | {molecule["molecule_chembl_id"] for molecule in molecules_by_ingredient.get(ingredient_id, [])}
        if chembl_parent:
            chembl_ids_seen.add(chembl_parent)
            pref_name = pref_name or hierarchy.get(chembl_parent, {}).get("pref_name")
        perturbation_id, sider_cid, collisions = unify_with_sider(ingredient_id, chembl_parent, cids_by_parent)
        stitch_flat_id = stitch_flat_id_from_pubchem_cid(sider_cid) if sider_cid is not None else None
        bridges[ingredient_id] = IngredientBridge(
            ingredient_id=ingredient_id, ingredient_name=str(name), identifier_system=identifier_systems[ingredient_id], bridge_method=bridge_method,
            unii_codes=unii_by_rxcui.get(ingredient_id, []), chembl_ids_seen=sorted(chembl_ids_seen, key=lambda chembl_id: int(chembl_id[6:])),
            chembl_parent=chembl_parent, chembl_parent_alternatives=alternatives, chembl_pref_name=pref_name,
            atc_codes_rxnav=atc_by_rxcui.get(ingredient_id, []), atc_codes_sider=sorted(set(sider_atc_codes.get(stitch_flat_id, []))) if stitch_flat_id else [],
            perturbation_id=perturbation_id, unified_with_sider=sider_cid is not None, sider_pubchem_cid=sider_cid,
            sider_drug_name=sider_drug_names.get(stitch_flat_id) if stitch_flat_id else None, sider_collision_pubchem_cids=collisions,
        )

    parents_for_targets = sorted({bridge.chembl_parent for bridge in bridges.values() if bridge.chembl_parent})
    targets_path = arguments.chembl_dir / "targets.json"
    targets = json.loads(targets_path.read_text()) if targets_path.exists() else {}
    target_ids_hit = {mechanism["target_chembl_id"] for parent in parents_for_targets for mechanism in mechanisms_by_molecule.get(parent, []) if mechanism.get("target_chembl_id")}
    missing_target_ids = sorted(target_ids_hit - set(targets))
    if missing_target_ids and not arguments.skip_target_fetch:
        print(f"  fetching {len(missing_target_ids)} targets missing from targets.json ({len(targets)} cached)", flush=True)
        targets = fetch_targets({"onsides_parents": parents_for_targets}, mechanisms, arguments.chembl_dir)
    elif missing_target_ids:
        errors.append({"phase": "chembl_targets", "error": f"{len(missing_target_ids)} targets not in targets.json and --skip-target-fetch given"})
    for bridge in bridges.values():
        drug_targets = mechanism_targets_for_parent(bridge.chembl_parent, mechanisms_by_molecule, targets)
        bridge.mechanism_target_ids = sorted(target.target_chembl_id for target in drug_targets)
        bridge.single_mechanism_target = passes_single_target_rule(drug_targets)
        bridge.nervous_system_atc = passes_nervous_system_rule(bridge)

    counts = bridge_counts(bridges)
    name_bridge_comparison = sider_name_bridge_comparison(ingredients, bridges, sider_drug_names)
    document = {
        "release": ONSIDES_RELEASE, "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "elapsed_seconds": round(time.time() - started, 1),
        "sources": {"rxnav": RXNAV_API, "unichem": UNICHEM_V1_COMPOUNDS, "unichem_fda_srs_source_id": UNICHEM_FDA_SRS_SOURCE_ID, "chembl": CHEMBL_API, "chembl_release": json.loads((arguments.chembl_dir / "chembl_release.json").read_text()).get("chembl_db_version") if (arguments.chembl_dir / "chembl_release.json").exists() else None},
        "unichem_routes": unichem_routes, "requests": {"rxnav": rxnav_client.request_count, "unichem": unichem_client.request_count, "chembl": chembl_client.request_count},
        "cache_hits": {"rxnav": rxnav_client.cache_hit_count, "unichem": unichem_client.cache_hit_count, "chembl": chembl_client.cache_hit_count},
        "counts": counts, "sider_name_bridge_comparison": name_bridge_comparison, "sider_collisions": {bridge.chembl_parent: [bridge.sider_pubchem_cid, *bridge.sider_collision_pubchem_cids] for bridge in bridges.values() if bridge.sider_collision_pubchem_cids},
        "errors": errors, "ingredients": {ingredient_id: bridge_to_json_record(bridge) for ingredient_id, bridge in bridges.items()},
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(".json.part")
    temporary_path.write_text(json.dumps(document, indent=1, ensure_ascii=False), encoding="utf-8")
    temporary_path.replace(output_path)
    print(json.dumps({"counts": counts, "name_bridge": name_bridge_comparison, "errors": len(errors), "elapsed_seconds": document["elapsed_seconds"]}, indent=1), flush=True)
    print(f"wrote {output_path}", flush=True)


if __name__ == "__main__":
    main()
