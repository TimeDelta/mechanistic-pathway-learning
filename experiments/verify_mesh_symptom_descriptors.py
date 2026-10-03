"""Verify the MeSH descriptors of docs/symptom_crosswalk.csv against the live vocabularies (design section 3.1).

For every descriptor in the crosswalk's mesh_descriptors column:
  1. the preferred label and tree numbers are read from the NLM MeSH RDF service (id.nlm.nih.gov/mesh/{id}.json),
  2. the PubTator3 autocomplete endpoint is queried with that label; an entity whose db_id is the descriptor shows
     that PubTator3 annotates mentions to this exact id. PubTator3 normalizes diseases to the MEDIC vocabulary (MeSH
     branches C and F03), so a behavior-branch (F01) descriptor such as Depression D003863 is a valid MeSH id that
     PubTator3 and CTD never use; the table records both facts apart (mesh_valid, in_pubtator3_vocabulary),
  3. the tree numbers are recorded next to the crosswalk's level flag: MeSH files signs and symptoms under
     C23.888 (Signs and Symptoms) or F01 (Behavior and Behavior Mechanisms), and diagnoses under F03 (Mental
     Disorders) or a C10 nervous-system disease branch, so the flag can be checked against the tree.

Every response is cached under data/raw/pubtator3/autocomplete and data/raw/pubtator3/mesh_lookup (gitignored), the
NCBI rate limit is respected, and the result is written to data/raw/pubtator3/mesh_descriptor_verification.json
plus a Markdown table on standard output for docs/literature_survey_spec.md. Only the identifiers and the
descriptor names of the target symptoms are recorded, not the MeSH entry-term lists.

Usage:
  OMP_NUM_THREADS=1 python experiments/verify_mesh_symptom_descriptors.py
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
from pathlib import Path

from mechanistic_pathway_learning.evidence.cached_http import CachedRateLimitedClient
from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import load_symptom_mesh_descriptors

PUBTATOR_AUTOCOMPLETE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/"
MESH_RDF_URL = "https://id.nlm.nih.gov/mesh/"
SYMPTOM_TREE_PREFIXES = ("C23.888", "F01", "C10.597")  # Signs and Symptoms, Behavior, Neurologic Manifestations
DIAGNOSIS_TREE_PREFIXES = ("F03", "C10.228", "C10.886")  # Mental Disorders, Central Nervous System Diseases, Sleep Wake Disorders


def mesh_descriptor_record(client: CachedRateLimitedClient, descriptor: str) -> dict:
    """Preferred label and tree numbers of one descriptor from the MeSH RDF service."""
    document = client.fetch_json(f"{MESH_RDF_URL}{descriptor}.json")
    label = document.get("label", {}).get("@value") if isinstance(document.get("label"), dict) else document.get("label")
    tree_numbers = document.get("treeNumber", [])
    if isinstance(tree_numbers, str):
        tree_numbers = [tree_numbers]
    return {"label": label, "tree_numbers": sorted(tree.rsplit("/", 1)[-1] for tree in tree_numbers), "active": document.get("active")}


def pubtator_autocomplete(client: CachedRateLimitedClient, query: str) -> list[dict]:
    return client.fetch_json(f"{PUBTATOR_AUTOCOMPLETE_URL}?query={urllib.parse.quote(query)}")


def tree_level_opinion(tree_numbers: list[str]) -> str:
    """What the tree numbers alone would say about the level, for comparison with the crosswalk flag."""
    symptom_branch = any(tree.startswith(SYMPTOM_TREE_PREFIXES) for tree in tree_numbers)
    diagnosis_branch = any(tree.startswith(DIAGNOSIS_TREE_PREFIXES) for tree in tree_numbers)
    if symptom_branch and not diagnosis_branch:
        return "symptom"
    if diagnosis_branch and not symptom_branch:
        return "diagnosis"
    if symptom_branch and diagnosis_branch:
        return "both branches"
    return "neither branch"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/pubtator3"))
    parser.add_argument("--output", type=Path, default=Path("data/raw/pubtator3/mesh_descriptor_verification.json"))
    arguments = parser.parse_args()

    autocomplete_client = CachedRateLimitedClient(arguments.cache_dir / "autocomplete")
    mesh_client = CachedRateLimitedClient(arguments.cache_dir / "mesh_lookup")
    verification: list[dict] = []
    for entry in load_symptom_mesh_descriptors(arguments.crosswalk):
        record = mesh_descriptor_record(mesh_client, entry.mesh_descriptor)
        candidates = pubtator_autocomplete(autocomplete_client, record["label"] or entry.mesh_descriptor)
        matching = [candidate for candidate in candidates if candidate.get("db_id") == entry.mesh_descriptor]
        verification.append(
            {
                "target_symptom": entry.target_symptom,
                "mesh_descriptor": entry.mesh_descriptor,
                "level": entry.level,
                "mesh_label": record["label"],
                "tree_numbers": record["tree_numbers"],
                "tree_level_opinion": tree_level_opinion(record["tree_numbers"]),
                "pubtator_entity_id": matching[0]["_id"] if matching else None,
                "pubtator_name": matching[0]["name"] if matching else None,
                "pubtator_match": matching[0]["match"] if matching else None,
                "mesh_valid": record["label"] is not None,
                "in_pubtator3_vocabulary": bool(matching),
            }
        )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(verification, indent=2))

    print("| Target symptom | MeSH id | Level | MeSH preferred label | Tree numbers | Tree opinion | MeSH valid | PubTator3 entity |")
    print("|---|---|---|---|---|---|---|---|")
    for item in verification:
        print(f"| {item['target_symptom']} | {item['mesh_descriptor']} | {item['level']} | {item['mesh_label']} | {'; '.join(item['tree_numbers'])} | {item['tree_level_opinion']} | {'yes' if item['mesh_valid'] else 'NO'} | {item['pubtator_entity_id'] or 'none (not in the PubTator3 disease vocabulary)'} |")
    invalid = [item for item in verification if not item["mesh_valid"]]
    absent = [item for item in verification if item["mesh_valid"] and not item["in_pubtator3_vocabulary"]]
    print(f"\n{len(verification)} descriptors: {len(invalid)} not valid MeSH ids, {len(absent)} valid but absent from the PubTator3 vocabulary; requests made: autocomplete {autocomplete_client.request_count}, mesh {mesh_client.request_count}")
    if invalid:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
