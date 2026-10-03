"""Fetch the documents behind the filtered literature relations: title, abstract, year, journal and species
annotations from the PubTator3 export, and publication types from E-utilities esummary; cached and resumable.

Usage:
  python experiments/fetch_pubtator_documents.py [--max-pmids 30000]
Writes data/processed/literature/documents.parquet. When the PMIDs exceed the cap, the PMIDs of cause, treat,
prevent and correlation relations are fetched first and associate relations fill the rest; the number left
unfetched is logged and written to the summary (no silent cap).
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.cached_http import CachedRateLimitedClient

PUBTATOR_EXPORT_URL = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/publications/export/biocjson"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PRIORITY_RELATION_TYPES = ("cause", "treat", "prevent", "positive_correlation", "negative_correlation", "positive_correlate", "negative_correlate")


def pmids_in_priority_order(relation_tables: list[pd.DataFrame]) -> tuple[list[str], list[str]]:
    """Distinct PMIDs with directional relations first, then the rest (associate relations and CTD rows without a type)."""
    priority: list[str] = []
    rest: list[str] = []
    seen: set[str] = set()
    for table in relation_tables:
        if table is None or "pmid" not in table.columns:
            continue
        for pmid, relation_type in zip(table.pmid.astype(str), table.get("relation_type", pd.Series([""] * len(table))).astype(str)):
            if pmid in seen or not pmid or pmid == "nan":
                continue
            seen.add(pmid)
            (priority if relation_type in PRIORITY_RELATION_TYPES else rest).append(pmid)
    return priority, rest


def parse_export(document: dict) -> list[dict]:
    rows: list[dict] = []
    for article in document.get("PubTator3", []):
        pmid = str(article.get("id") or article.get("_id", "").split("|")[0])
        title, abstract, year, journal = "", "", None, ""
        species_ids: set[str] = set()
        species_names: set[str] = set()
        for passage in article.get("passages", []):
            infons = passage.get("infons", {})
            if infons.get("type") == "title":
                title = passage.get("text", "")
                year = infons.get("year") or year
                journal = infons.get("journal", journal)
            elif infons.get("type") == "abstract":
                abstract = passage.get("text", "")
            for annotation in passage.get("annotations", []):
                annotation_infons = annotation.get("infons", {})
                if annotation_infons.get("type") == "Species":
                    identifier = str(annotation_infons.get("identifier") or annotation_infons.get("normalized_id") or "")
                    if identifier:
                        species_ids.add(identifier)
                        species_names.add(annotation_infons.get("name") or annotation.get("text", ""))
        rows.append({"pmid": pmid, "year": int(year) if year and str(year).isdigit() else None, "journal": journal, "title": title, "abstract": abstract,
                     "species_ids": ";".join(sorted(species_ids)), "species_names": ";".join(sorted(species_names))})
    return rows


def fetch_documents(client: CachedRateLimitedClient, pmids: list[str], batch_size: int = 100) -> list[dict]:
    rows: list[dict] = []
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start:start + batch_size]
        url = f"{PUBTATOR_EXPORT_URL}?pmids={','.join(batch)}"
        try:
            rows += parse_export(json.loads(client.fetch_text(url)))
        except Exception as error:  # a failed batch is logged, not fatal; the rerun retries it
            print(f"export batch at {start} failed: {error}")
    return rows


def fetch_publication_types(client: CachedRateLimitedClient, pmids: list[str], batch_size: int = 200) -> dict[str, list[str]]:
    types: dict[str, list[str]] = {}
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start:start + batch_size]
        body = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "id": ",".join(batch)}).encode("utf-8")
        try:
            result = json.loads(client.fetch_text(ESUMMARY_URL, post_body=body)).get("result", {})
        except Exception as error:
            print(f"esummary batch at {start} failed: {error}")
            continue
        for pmid in batch:
            summary = result.get(pmid)
            if isinstance(summary, dict) and "error" not in summary:
                types[pmid] = list(summary.get("pubtype", []))
    return types


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--literature-dir", type=Path, default=Path("data/processed/literature"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/pubtator3"))
    parser.add_argument("--max-pmids", type=int, default=30000)
    arguments = parser.parse_args()
    tables = [pd.read_parquet(path) for path in (arguments.literature_dir / "pubtator_relations_filtered.parquet", arguments.literature_dir / "ctd_relations_filtered.parquet") if path.exists()]
    priority, rest = pmids_in_priority_order(tables)
    selected = (priority + rest)[: arguments.max_pmids]
    unfetched = max(0, len(priority) + len(rest) - len(selected))
    print(f"pmids: {len(priority)} with directional relations, {len(rest)} other; fetching {len(selected)}, leaving {unfetched} unfetched")
    export_client = CachedRateLimitedClient(arguments.cache_dir / "documents")
    esummary_client = CachedRateLimitedClient(arguments.cache_dir / "esummary")
    documents = pd.DataFrame(fetch_documents(export_client, selected))
    publication_types = fetch_publication_types(esummary_client, selected)
    if len(documents):
        documents["publication_types"] = documents.pmid.map(lambda pmid: ";".join(publication_types.get(str(pmid), [])))
    documents.to_parquet(arguments.literature_dir / "documents.parquet", index=False)
    summary = {"pmids_with_directional_relations": len(priority), "pmids_other": len(rest), "pmids_fetched": len(selected), "pmids_unfetched": unfetched,
               "documents_returned": int(len(documents)), "documents_with_species": int((documents.species_ids != "").sum()) if len(documents) else 0,
               "documents_with_publication_types": int((documents.publication_types != "").sum()) if len(documents) else 0}
    (arguments.literature_dir / "documents_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
