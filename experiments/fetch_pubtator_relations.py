"""Fetch the PubTator3 half of evidence class E3: bulk relations filtered to target-symptom descriptors and graph
perturbations (design sections 4.2, 4.3 and 6.1; docs/literature_survey_spec.md).

Steps, each cached under data/raw/pubtator3 or data/processed/literature so a rerun makes no request:
  1. Pin the bulk file relation2pubtator3.gz: HEAD gives Content-Length and Last-Modified; the file is downloaded when
     it is under --max-download-bytes (400 MB by default), otherwise it is streamed through gzip and only the rows
     with a target descriptor are kept. The pin (URL, size, date, SHA-256 when downloaded, row count) is written to
     relation2pubtator3.pin.json.
  2. Gene identifier map (gene_identifier_map.py): NCBI Gene id -> graph symbol.
  3. SIDER drug names -> MeSH chemical ids through the PubTator3 autocomplete endpoint, one cached query per drug
     name at three requests per second (sider_chemical_matching.py records the match method).
  4. One pass over the bulk file keeping every row with a target descriptor on either side (count per relation type
     and per descriptor), then the partner filters in memory.
  5. Chemical ids on kept rows that no SIDER name matched are looked up by MeSH label (id.nlm.nih.gov) and matched
     again by name; --mesh-label-lookup-limit caps that pass and the cap is logged with the number left unresolved.
     One SIDER drug is then kept per MeSH chemical id (sider_chemical_matching.collapse_matches_to_one_drug_per_mesh_id,
     preferring the exact-name match, then a STITCH id with a ChEMBL mapping in data/raw/chembl, then the lowest CID),
     so a paper is never counted under a brand, salt or code-name entry as well; the dropped ids go to the summary.
  6. With --fetch-publication-dates, every distinct PMID of the kept rows is dated through esummary (200 per
     request) so the time split can use evidence_date and publication_year.
  7. data/processed/literature/pubtator_relations_filtered.parquet in the evidence_reports.parquet schema plus the
     literature columns, with pubtator_relations_summary.json holding every count and every drop reason.

Literature reports are soft priors only (grades C to E): never evaluation positives; unobserved pairs stay
unlabelled.

Usage:
  OMP_NUM_THREADS=1 python experiments/fetch_pubtator_relations.py --fetch-publication-dates
"""
from __future__ import annotations

import argparse
import hashlib
import gzip
import json
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.cached_http import USER_AGENT, CachedRateLimitedClient
from mechanistic_pathway_learning.evidence.gene_identifier_map import build_gene_identifier_map, gene_map_summary, ncbi_gene_id_to_symbol
from mechanistic_pathway_learning.evidence.load_pubtator_relations import (
    disease_rows_from_stream,
    filter_relation_rows,
    iterate_bulk_relation_rows,
    parse_entity,
    relation_rows_to_reports,
    summarize_reports,
)
from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import descriptor_lookup, load_symptom_mesh_descriptors
from mechanistic_pathway_learning.evidence.pubmed_publication_dates import fetch_publication_dates
from mechanistic_pathway_learning.evidence.sider_chemical_matching import (
    ChemicalMatch,
    chemical_matches_by_mesh_id,
    collapse_matches_to_one_drug_per_mesh_id,
    drug_names_by_normalized_name,
    load_sider_drug_names,
    match_autocomplete_candidates,
    match_mesh_label,
    pubchem_cid_of_stitch_flat_id,
)

BULK_RELATION_URL = "https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/relation2pubtator3.gz"
BULK_FILE_NAME = "relation2pubtator3.gz"
PUBTATOR_AUTOCOMPLETE_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/entity/autocomplete/"
MESH_LABEL_LOOKUP_URL = "https://id.nlm.nih.gov/mesh/lookup/label"
DEFAULT_MAX_DOWNLOAD_BYTES = 400 * 1024 * 1024


def head_headers(url: str) -> dict[str, str]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return {key.lower(): value for key, value in response.headers.items()}


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pin_or_download_bulk_file(raw_directory: Path, url: str, max_download_bytes: int) -> tuple[Path | None, dict]:
    """The local bulk file when it exists or fits the download cap, else None (stream it); the pin so far."""
    raw_directory.mkdir(parents=True, exist_ok=True)
    target_path = raw_directory / BULK_FILE_NAME
    pin_path = raw_directory / "relation2pubtator3.pin.json"
    headers = head_headers(url)
    content_length = int(headers.get("content-length", "0") or 0)
    pin = {"url": url, "content_length": content_length, "last_modified": headers.get("last-modified"), "checked_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if pin_path.exists():
        pin.update({key: value for key, value in json.loads(pin_path.read_text()).items() if key in ("sha256", "downloaded_at", "bulk_rows_total")})
    if target_path.exists() and target_path.stat().st_size == content_length:
        pin["local_file"] = str(target_path)
        pin.setdefault("sha256", sha256_of_file(target_path))
        pin_path.write_text(json.dumps(pin, indent=2))
        return target_path, pin
    if content_length and content_length <= max_download_bytes:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        temporary_path = target_path.with_suffix(".gz.part")
        with urllib.request.urlopen(request, timeout=600) as response, open(temporary_path, "wb") as handle:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                handle.write(chunk)
        temporary_path.replace(target_path)
        pin.update({"local_file": str(target_path), "downloaded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "sha256": sha256_of_file(target_path)})
        pin_path.write_text(json.dumps(pin, indent=2))
        return target_path, pin
    pin["local_file"] = None
    pin["streamed_without_download"] = True
    pin_path.write_text(json.dumps(pin, indent=2))
    return None, pin


def stream_bulk_rows(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=600) as response:
        with gzip.open(response, "rt", encoding="utf-8", newline="") as handle:
            yield from iterate_bulk_relation_rows(handle)


def match_sider_drugs_by_autocomplete(client: CachedRateLimitedClient, drug_names: dict[str, str]) -> tuple[list[ChemicalMatch], Counter]:
    matches: list[ChemicalMatch] = []
    counts: Counter = Counter()
    for stitch_flat_id, drug_name in sorted(drug_names.items()):
        url = f"{PUBTATOR_AUTOCOMPLETE_URL}?query={urllib.parse.quote(drug_name)}"
        try:
            candidates = client.fetch_json(url)
        except Exception as error:  # noqa: BLE001 - a name the service rejects is recorded, not fatal
            counts["autocomplete_query_failed"] += 1
            counts[f"autocomplete_query_failed_example:{drug_name}"] = str(error)[:80]
            continue
        if not isinstance(candidates, list):
            counts["autocomplete_no_candidates"] += 1
            continue
        match = match_autocomplete_candidates(stitch_flat_id, drug_name, candidates)
        if match is None:
            counts["drugs_unmatched_by_autocomplete"] += 1
        else:
            matches.append(match)
            counts[f"drugs_matched:{match.match_method}"] += 1
    return matches, counts


def chemical_ids_on_rows(rows: list[tuple[str, str, str, str]]) -> Counter:
    """Chemical MeSH ids on the partner side of the kept disease rows, with their row counts."""
    counts: Counter = Counter()
    for _, _, first_field, second_field in rows:
        for field in (first_field, second_field):
            entity_type, identifiers = parse_entity(field)
            if entity_type == "Chemical":
                for identifier in identifiers:
                    counts[identifier] += 1
    return counts


def match_unresolved_chemicals_by_mesh_label(client: CachedRateLimitedClient, chemical_row_counts: Counter, already_matched: set[str], names_index: dict, lookup_limit: int | None) -> tuple[list[ChemicalMatch], Counter]:
    """MeSH label lookup for chemical ids on kept rows that no SIDER name matched, most frequent first; a cap is logged
    with the number of ids and rows it leaves unresolved."""
    unresolved = [(identifier, count) for identifier, count in chemical_row_counts.most_common() if identifier not in already_matched]
    counts: Counter = Counter({"chemical_ids_unresolved_after_autocomplete": len(unresolved), "chemical_rows_unresolved_after_autocomplete": sum(count for _, count in unresolved)})
    to_query = unresolved if lookup_limit is None else unresolved[:lookup_limit]
    if lookup_limit is not None and len(unresolved) > lookup_limit:
        counts["mesh_label_lookup_cap"] = lookup_limit
        counts["chemical_ids_dropped_by_lookup_cap"] = len(unresolved) - lookup_limit
        counts["chemical_rows_dropped_by_lookup_cap"] = sum(count for _, count in unresolved[lookup_limit:])
    matches: list[ChemicalMatch] = []
    for identifier, _ in to_query:
        url = f"{MESH_LABEL_LOOKUP_URL}?resource={urllib.parse.quote(f'http://id.nlm.nih.gov/mesh/{identifier}', safe='')}"
        try:
            labels = client.fetch_json(url)
        except Exception:  # noqa: BLE001 - unknown id at NLM, recorded
            counts["mesh_label_lookup_failed"] += 1
            continue
        counts["mesh_label_lookups"] += 1
        for label in labels if isinstance(labels, list) else []:
            hits = match_mesh_label(identifier, str(label), names_index)
            matches.extend(hits)
            if hits:
                counts["chemical_ids_matched_by_mesh_label"] += 1
    return matches, counts


def stitch_ids_with_chembl_mapping(chembl_directory: Path, stitch_ids: set[str]) -> set[str]:
    """The STITCH flat ids among stitch_ids whose PubChem CID has a non-empty ChEMBL list in pubchem_to_chembl_with_parents.json
    (the drugs the assembler can join to targets); empty when the file is absent."""
    mapping_path = Path(chembl_directory) / "pubchem_to_chembl_with_parents.json"
    if not mapping_path.exists():
        return set()
    mapping = json.loads(mapping_path.read_text())
    preferred: set[str] = set()
    for stitch_id in stitch_ids:
        entry = mapping.get(str(pubchem_cid_of_stitch_flat_id(stitch_id)))
        if entry:
            preferred.add(stitch_id)
    return preferred


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bulk-url", default=BULK_RELATION_URL)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/pubtator3"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/literature"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full"))
    parser.add_argument("--human-gem-genes", type=Path, default=Path("data/raw/Human-GEM/model/genes.tsv"))
    parser.add_argument("--hgnc-dir", type=Path, default=Path("data/raw/hgnc"))
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"), help="pubchem_to_chembl_with_parents.json names the STITCH ids preferred when a MeSH chemical matches several SIDER drugs")
    parser.add_argument("--max-download-bytes", type=int, default=DEFAULT_MAX_DOWNLOAD_BYTES)
    parser.add_argument("--mesh-label-lookup-limit", type=int, default=None, help="cap on MeSH label lookups for unresolved chemical ids; logged when it drops anything")
    parser.add_argument("--fetch-publication-dates", action="store_true")
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"run_started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    bulk_path, pin = pin_or_download_bulk_file(arguments.raw_dir, arguments.bulk_url, arguments.max_download_bytes)
    print(f"bulk file: {bulk_path or 'streamed'} ({pin['content_length']} bytes, last modified {pin['last_modified']})", flush=True)

    descriptors = descriptor_lookup(load_symptom_mesh_descriptors(arguments.crosswalk))
    summary["descriptors"] = {identifier: {"symptom": entry.target_symptom, "level": entry.level} for identifier, entry in descriptors.items()}

    gene_map_path = arguments.output_dir / "gene_identifier_map.parquet"
    gene_map = pd.read_parquet(gene_map_path) if gene_map_path.exists() else build_gene_identifier_map(arguments.graph_dir, arguments.human_gem_genes, arguments.hgnc_dir, gene_map_path)
    ncbi_to_symbol = ncbi_gene_id_to_symbol(gene_map)
    summary["gene_identifier_map"] = gene_map_summary(gene_map)
    print(f"gene map: {summary['gene_identifier_map']['with_ncbi_gene_id']} of {summary['gene_identifier_map']['gene_nodes']} gene nodes have an NCBI id", flush=True)

    drug_names = load_sider_drug_names(arguments.sider_dir)
    names_index = drug_names_by_normalized_name(drug_names)
    autocomplete_client = CachedRateLimitedClient(arguments.raw_dir / "autocomplete")
    autocomplete_matches, autocomplete_counts = match_sider_drugs_by_autocomplete(autocomplete_client, drug_names)
    summary["sider_drugs_total"] = len(drug_names)
    summary["autocomplete"] = dict(autocomplete_counts)
    summary["autocomplete_requests_made"] = autocomplete_client.request_count
    print(f"autocomplete: {len(autocomplete_matches)} of {len(drug_names)} SIDER drugs matched a MeSH chemical ({autocomplete_client.request_count} requests, {autocomplete_client.cache_hit_count} cache hits)", flush=True)

    rows_source = iterate_bulk_relation_rows(bulk_path) if bulk_path is not None else stream_bulk_rows(arguments.bulk_url)
    started = time.monotonic()
    disease_rows, stream_counts = disease_rows_from_stream(rows_source, descriptors)
    summary["bulk_stream"] = dict(stream_counts)
    summary["bulk_stream_seconds"] = round(time.monotonic() - started, 1)
    pin["bulk_rows_total"] = stream_counts["bulk_rows_total"]
    (arguments.raw_dir / "relation2pubtator3.pin.json").write_text(json.dumps(pin, indent=2))
    summary["bulk_pin"] = pin
    descriptor_row_counts: Counter = Counter()
    for _, _, first_field, second_field in disease_rows:
        for field in (first_field, second_field):
            entity_type, identifiers = parse_entity(field)
            if entity_type == "Disease":
                for identifier in identifiers:
                    if identifier in descriptors:
                        descriptor_row_counts[identifier] += 1
    summary["bulk_rows_by_descriptor"] = {identifier: int(descriptor_row_counts.get(identifier, 0)) for identifier in descriptors}
    print(f"bulk stream: {stream_counts['bulk_rows_total']} rows, {stream_counts['bulk_rows_with_target_descriptor']} with a target descriptor", flush=True)

    chemical_row_counts = chemical_ids_on_rows(disease_rows)
    matched_ids = {match.chemical_mesh_id for match in autocomplete_matches}
    mesh_client = CachedRateLimitedClient(arguments.raw_dir / "mesh_lookup")
    label_matches, label_counts = match_unresolved_chemicals_by_mesh_label(mesh_client, chemical_row_counts, matched_ids, names_index, arguments.mesh_label_lookup_limit)
    summary["chemical_ids_on_disease_rows"] = len(chemical_row_counts)
    summary["mesh_label_lookup"] = dict(label_counts)
    all_matches = autocomplete_matches + [match for match in label_matches if match.chemical_mesh_id not in matched_ids]
    preferred_stitch_ids = stitch_ids_with_chembl_mapping(arguments.chembl_dir, {match.stitch_flat_id for match in all_matches})
    kept_matches, collapsed_ids = collapse_matches_to_one_drug_per_mesh_id(all_matches, preferred_stitch_ids)
    matches_index = chemical_matches_by_mesh_id(kept_matches)
    mesh_map = pd.DataFrame([match.__dict__ for match in kept_matches])
    mesh_map.to_parquet(arguments.output_dir / "sider_chemical_mesh_map.parquet", index=False)
    summary["sider_drugs_matched_before_collapse"] = len({match.stitch_flat_id for match in all_matches})
    summary["sider_drugs_matched_total"] = len({match.stitch_flat_id for match in kept_matches})
    summary["mesh_ids_with_several_sider_drugs"] = len(collapsed_ids)
    summary["sider_drugs_collapsed_into_another"] = sorted({stitch_id for dropped in collapsed_ids.values() for stitch_id in dropped})
    summary["collapsed_sider_drugs_by_mesh_id"] = collapsed_ids
    summary["chemical_match_methods"] = dict(Counter(match.match_method for match in kept_matches))
    print(f"chemical matches: {len(all_matches)} before collapse, {len(kept_matches)} kept (one SIDER drug per MeSH id; {len(collapsed_ids)} ids had several)", flush=True)
    print(f"mesh label lookup: {label_counts.get('mesh_label_lookups', 0)} lookups, {label_counts.get('chemical_ids_matched_by_mesh_label', 0)} further ids matched", flush=True)

    kept_rows, filter_counts = filter_relation_rows(disease_rows, descriptors, ncbi_to_symbol, matches_index)
    summary["filter"] = dict(filter_counts)
    print(f"filter: {len(kept_rows)} rows kept", flush=True)

    publication_dates: dict[str, str] = {}
    if arguments.fetch_publication_dates:
        esummary_client = CachedRateLimitedClient(arguments.raw_dir / "esummary")
        publication_dates, date_counts = fetch_publication_dates(esummary_client, [row.pmid for row in kept_rows])
        summary["publication_dates"] = date_counts
        summary["esummary_requests_made"] = esummary_client.request_count
        (arguments.output_dir / "pubmed_publication_dates.json").write_text(json.dumps(publication_dates, sort_keys=True))
        print(f"publication dates: {date_counts}", flush=True)

    reports = relation_rows_to_reports(kept_rows, publication_dates)
    output_path = arguments.output_dir / "pubtator_relations_filtered.parquet"
    temporary_path = output_path.with_suffix(".parquet.part")
    reports.to_parquet(temporary_path, index=False)
    temporary_path.replace(output_path)
    summary["output"] = summarize_reports(reports)
    summary["run_finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (arguments.output_dir / "pubtator_relations_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["output"], indent=2))


if __name__ == "__main__":
    main()
