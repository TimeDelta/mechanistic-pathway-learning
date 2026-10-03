"""Fetch the CTD half of evidence class E3: curated chemical-disease statements on target-symptom descriptors for SIDER
drugs, one report per cited paper (design section 4.2; docs/literature_survey_spec.md).

Steps, cached under data/raw/ctd and data/processed/literature:
  1. Pin CTD_chemicals_diseases.tsv.gz (HEAD for Content-Length and Last-Modified; downloaded when absent; SHA-256,
     the "Report created" line of its header and the row count go to CTD_chemicals_diseases.pin.json).
  2. Stream the file, keep rows whose DiseaseID is a crosswalk descriptor and whose chemical is a SIDER drug by
     ChemicalName (case-insensitive) or, failing that, by a MeSH id the PubTator3 pass already matched
     (data/processed/literature/sider_chemical_mesh_map.parquet, when present).
  3. DirectEvidence marker/mechanism -> induces, therapeutic -> relieves, both -> two reports; inferred rows dropped
     with a count; one report per PubMed id.
  4. With --fetch-publication-dates the cited PMIDs are dated through esummary (shared cache with the PubTator3 pass).
  5. data/processed/literature/ctd_relations_filtered.parquet plus ctd_relations_summary.json.

Literature reports are soft priors only (grades C to E): never evaluation positives; unobserved pairs stay unlabelled.

Usage:
  OMP_NUM_THREADS=1 python experiments/fetch_ctd_chemical_disease.py --fetch-publication-dates
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.cached_http import USER_AGENT, CachedRateLimitedClient
from mechanistic_pathway_learning.evidence.load_ctd_relations import ctd_rows_to_reports, filter_ctd_rows, iterate_ctd_rows
from mechanistic_pathway_learning.evidence.load_pubtator_relations import summarize_reports
from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import descriptor_lookup, load_symptom_mesh_descriptors
from mechanistic_pathway_learning.evidence.pubmed_publication_dates import fetch_publication_dates
from mechanistic_pathway_learning.evidence.sider_chemical_matching import ChemicalMatch, chemical_matches_by_mesh_id, drug_names_by_normalized_name, load_sider_drug_names

CTD_CHEMICAL_DISEASE_URL = "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz"
CTD_FILE_NAME = "CTD_chemicals_diseases.tsv.gz"


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pin_or_download(raw_directory: Path, url: str) -> tuple[Path, dict]:
    raw_directory.mkdir(parents=True, exist_ok=True)
    target_path = raw_directory / CTD_FILE_NAME
    pin_path = raw_directory / "CTD_chemicals_diseases.pin.json"
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        headers = {key.lower(): value for key, value in response.headers.items()}
    content_length = int(headers.get("content-length", "0") or 0)
    pin = {"url": url, "content_length": content_length, "last_modified": headers.get("last-modified"), "checked_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    if pin_path.exists():
        pin.update({key: value for key, value in json.loads(pin_path.read_text()).items() if key in ("sha256", "downloaded_at", "report_created", "ctd_rows_total")})
    if not (target_path.exists() and target_path.stat().st_size == content_length):
        temporary_path = target_path.with_suffix(".gz.part")
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=600) as response, open(temporary_path, "wb") as handle:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                handle.write(chunk)
        temporary_path.replace(target_path)
        pin["downloaded_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    pin.setdefault("sha256", sha256_of_file(target_path))
    with gzip.open(target_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Report created:"):
                pin["report_created"] = line.partition(":")[2].strip()
                break
            if not line.startswith("#"):
                break
    pin["local_file"] = str(target_path)
    pin_path.write_text(json.dumps(pin, indent=2))
    return target_path, pin


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--url", default=CTD_CHEMICAL_DISEASE_URL)
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw/ctd"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/literature"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--esummary-cache-dir", type=Path, default=Path("data/raw/pubtator3/esummary"))
    parser.add_argument("--fetch-publication-dates", action="store_true")
    arguments = parser.parse_args()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    summary: dict = {"run_started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    ctd_path, pin = pin_or_download(arguments.raw_dir, arguments.url)
    print(f"CTD file: {ctd_path} ({pin['content_length']} bytes, last modified {pin['last_modified']}, report created {pin.get('report_created')})", flush=True)

    descriptors = descriptor_lookup(load_symptom_mesh_descriptors(arguments.crosswalk))
    drug_names = load_sider_drug_names(arguments.sider_dir)
    names_index = drug_names_by_normalized_name(drug_names)
    mesh_map_path = arguments.output_dir / "sider_chemical_mesh_map.parquet"
    matches_by_mesh_id: dict[str, list[ChemicalMatch]] = {}
    if mesh_map_path.exists():
        mesh_map = pd.read_parquet(mesh_map_path)
        matches_by_mesh_id = chemical_matches_by_mesh_id([ChemicalMatch(**row) for row in mesh_map.to_dict("records")])
    summary["sider_drugs_total"] = len(drug_names)
    summary["pubtator_mesh_map_rows"] = int(sum(len(matches) for matches in matches_by_mesh_id.values()))

    started = time.monotonic()
    kept_rows, counts = filter_ctd_rows(iterate_ctd_rows(ctd_path), descriptors, names_index, matches_by_mesh_id)
    summary["filter"] = dict(counts)
    summary["stream_seconds"] = round(time.monotonic() - started, 1)
    pin["ctd_rows_total"] = counts["ctd_rows_total"]
    (arguments.raw_dir / "CTD_chemicals_diseases.pin.json").write_text(json.dumps(pin, indent=2))
    summary["ctd_pin"] = pin
    print(f"filter: {counts['ctd_rows_total']} rows, {counts['ctd_rows_with_target_descriptor']} on target descriptors, {counts.get('reports_kept', 0)} reports kept", flush=True)

    publication_dates: dict[str, str] = {}
    if arguments.fetch_publication_dates:
        esummary_client = CachedRateLimitedClient(arguments.esummary_cache_dir)
        publication_dates, date_counts = fetch_publication_dates(esummary_client, [row.pmid for row in kept_rows if row.pmid])
        summary["publication_dates"] = date_counts
        summary["esummary_requests_made"] = esummary_client.request_count
        print(f"publication dates: {date_counts}", flush=True)

    reports = ctd_rows_to_reports(kept_rows, publication_dates)
    output_path = arguments.output_dir / "ctd_relations_filtered.parquet"
    temporary_path = output_path.with_suffix(".parquet.part")
    reports.to_parquet(temporary_path, index=False)
    temporary_path.replace(output_path)
    summary["output"] = summarize_reports(reports)
    summary["output"]["rows_by_direct_evidence"] = {key: int(value) for key, value in reports.direct_evidence.value_counts().sort_index().items()} if not reports.empty else {}
    summary["run_finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (arguments.output_dir / "ctd_relations_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["output"], indent=2))


if __name__ == "__main__":
    main()
