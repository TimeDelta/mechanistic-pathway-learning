"""Fetch the documents behind the filtered literature relations: title, abstract, year, journal and species
annotations from the PubTator3 export, and publication types from E-utilities esummary; cached and resumable.

Usage:
  python experiments/fetch_pubtator_documents.py [--max-pmids 30000]
Writes data/processed/literature/documents.parquet. When the PMIDs exceed the cap, the PMIDs that carry any directional
relation in any table (PubTator3 cause, treat, prevent, stimulate, inhibit and the correlation types; CTD marker/mechanism
and therapeutic) are fetched first, in first-seen order, and associate-only PMIDs fill the rest; a PMID is judged on all
its rows, not on the first one seen. Documents already in documents.parquet are kept without a request, so a rerun after
a change of priority fetches only the newly selected PMIDs. The number left unfetched is logged and written to the
summary (no silent cap).
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
PRIORITY_RELATION_TYPES = ("cause", "treat", "prevent", "stimulate", "inhibit", "positive_correlation", "negative_correlation", "positive_correlate", "negative_correlate", "marker/mechanism", "therapeutic")


def pmids_in_priority_order(relation_tables: list[pd.DataFrame]) -> tuple[list[str], list[str]]:
    """Distinct PMIDs that carry a directional relation on any of their rows (any table) first, in first-seen order,
    then the rest (associate-only PMIDs). A PMID whose first row is an associate relation and whose later row is a
    cause relation is a priority PMID."""
    directional: set[str] = set()
    order: list[str] = []
    seen: set[str] = set()
    for table in relation_tables:
        if table is None or "pmid" not in table.columns:
            continue
        relation_types = table["relation_type"].astype(str) if "relation_type" in table.columns else pd.Series([""] * len(table), index=table.index)
        for pmid, relation_type in zip(table.pmid.astype(str), relation_types):
            if not pmid or pmid == "nan":
                continue
            if pmid not in seen:
                seen.add(pmid)
                order.append(pmid)
            if relation_type in PRIORITY_RELATION_TYPES:
                directional.add(pmid)
    priority = [pmid for pmid in order if pmid in directional]
    rest = [pmid for pmid in order if pmid not in directional]
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
    documents_path = arguments.literature_dir / "documents.parquet"
    existing = pd.read_parquet(documents_path) if documents_path.exists() else pd.DataFrame()
    existing_pmids = set(existing.pmid.astype(str)) if len(existing) else set()
    to_fetch = [pmid for pmid in selected if pmid not in existing_pmids]
    all_pmids = set(priority) | set(rest)
    print(f"pmids: {len(priority)} with directional relations, {len(rest)} other; selected {len(selected)} (cap {arguments.max_pmids}); {len(existing_pmids)} already fetched, fetching {len(to_fetch)} new")
    export_client = CachedRateLimitedClient(arguments.cache_dir / "documents")
    esummary_client = CachedRateLimitedClient(arguments.cache_dir / "esummary")
    new_documents = pd.DataFrame(fetch_documents(export_client, to_fetch))
    publication_types = fetch_publication_types(esummary_client, to_fetch)
    if len(new_documents):
        new_documents["publication_types"] = new_documents.pmid.map(lambda pmid: ";".join(publication_types.get(str(pmid), [])))
    documents = pd.concat([table for table in (existing, new_documents) if len(table)], ignore_index=True, sort=False) if (len(existing) or len(new_documents)) else new_documents
    if len(documents):
        documents = documents.drop_duplicates(subset="pmid", keep="first")
    documents.to_parquet(documents_path, index=False)
    fetched_pmids = set(documents.pmid.astype(str)) if len(documents) else set()
    unfetched = len(all_pmids - fetched_pmids)
    summary = {"pmids_with_directional_relations": len(priority), "pmids_other": len(rest), "pmids_selected_by_cap": len(selected), "max_pmids": arguments.max_pmids,
               "pmids_newly_fetched": len(to_fetch), "pmids_fetched": int(len(fetched_pmids)), "pmids_unfetched": unfetched,
               "directional_pmids_without_document": int(len(set(priority) - fetched_pmids)),
               "documents_returned": int(len(documents)), "documents_with_species": int((documents.species_ids != "").sum()) if len(documents) else 0,
               "documents_with_publication_types": int((documents.publication_types != "").sum()) if len(documents) else 0}
    (arguments.literature_dir / "documents_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
