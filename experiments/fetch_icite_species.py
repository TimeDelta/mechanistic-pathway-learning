"""Fetch, for every PubMed id behind a curated signalling or transcription edge, the species its MeSH terms record,
from NIH iCite (decision 9.63: "Curated edges for other species should not be considered labels unless there's
corroborating evidence from humans", the owner, 16 October 2026).

What is read. iCite's per-article human, animal and molecular/cellular scores, which NIH's Office of Portfolio
Analysis computes "for each article using Medical Subject Heading terms assigned by the National Library of Medicine
(NLM)": Human counts MeSH "M01 (Persons)" and the Humans descriptor, Animal "B01 (Eukaryota) not ... (Human)",
Molecular/Cellular the cell, microbe and chemical or cell-physiology branches (support.icite.nih.gov, "Human, Animal,
and Molecular/Cellular Scores", read 17 October 2026; Weber 2013, J Transl Med 11:126; Hutchins et al. 2019, PLoS
Biol 17(10):e3000416). A paper is human evidence when its human score is above 0. The scores are per paper, not per
interaction: a paper that shows an interaction in mice and mentions a human cohort counts as human evidence, which is
the coarseness of the route the owner accepted on 16 October 2026 (docs/open_decisions.md, 9.63).

Licence: the NIH Open Citation Collection that iCite serves is described by its makers as "a comprehensive, public
domain citation graph for biomedical research made freely available to the community" (Hutchins et al. 2019, PLoS
Biol 17(10):e3000385, PMC6786512, read 17 October 2026); the scores are NIH's derivation from NLM's MeSH indexing.

The API: https://icite.od.nih.gov/api/pubs?pmids=...&fl=..., no key and no e-mail; batches of BATCH_SIZE ids, since a
request naming 1,000 was refused as too large (HTTP 413); the requests are paced (PAUSE_SECONDS) and cached under
data/raw/icite, so a rerun asks only for what is missing.

Usage:
  PYTHONPATH=. python experiments/fetch_icite_species.py
  PYTHONPATH=. python experiments/fetch_icite_species.py --limit 2000   # a first batch only
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.edge_evidence import TEXT_MINED_REFERENCE_RESOURCES, pmids_of

API = "https://icite.od.nih.gov/api/pubs"
FIELDS = ("pmid", "human", "animal", "molecular_cellular", "year", "is_research_article")
BATCH_SIZE = 200  # the API takes up to 1,000 ids, but a GET of 1,000 was refused with HTTP 413 (17 October 2026)
PAUSE_SECONDS = 2.0
USER_AGENT = "perturbation-symptom-modules (research; https://github.com/TimeDelta/perturbation-symptom-modules)"
DEFAULT_CACHE = Path("data/raw/icite/icite_species.parquet")
DEFAULT_OMNIPATH = Path("data/raw/omnipath/omnipath_interactions.tsv")
DEFAULT_COLLECTRI = Path("data/raw/omnipath/collectri_interactions.tsv")


def curated_pubmed_ids(omnipath_path: Path, collectri_path: Path) -> list[str]:
    """Every PubMed id a curated signalling or transcription edge can carry: OmniPath's references outside the
    text-mined resource and CollecTRI's references, sorted."""
    pmids: set[str] = set()
    for references in pd.read_csv(omnipath_path, sep="\t", usecols=["references"], low_memory=False).references:
        kept, _ = pmids_of(references, TEXT_MINED_REFERENCE_RESOURCES)
        pmids |= kept
    for references in pd.read_csv(collectri_path, sep="\t", usecols=["references"], low_memory=False).references:
        kept, _ = pmids_of(references)
        pmids |= kept
    return sorted(pmids, key=int)


def fetch_batch(pmids: list[str]) -> list[dict]:
    query = urllib.parse.urlencode({"pmids": ",".join(pmids), "fl": ",".join(FIELDS), "limit": str(len(pmids))})
    request = urllib.request.Request(f"{API}?{query}", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=120) as response:
        body = json.load(response)
    return body.get("data", [])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--omnipath", type=Path, default=DEFAULT_OMNIPATH)
    parser.add_argument("--collectri", type=Path, default=DEFAULT_COLLECTRI)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--limit", type=int, default=None, help="ask for at most this many missing ids")
    parser.add_argument("--pause", type=float, default=PAUSE_SECONDS)
    arguments = parser.parse_args()
    wanted = curated_pubmed_ids(arguments.omnipath, arguments.collectri)
    cached = pd.read_parquet(arguments.cache) if arguments.cache.exists() else pd.DataFrame(columns=[*FIELDS, "found", "retrieved_at"])
    have = set(cached.pmid.astype(str)) if len(cached) else set()
    missing = [pmid for pmid in wanted if pmid not in have]
    if arguments.limit is not None:
        missing = missing[: arguments.limit]
    print(f"{len(wanted):,} PubMed ids behind curated edges; {len(have):,} cached; asking for {len(missing):,}", flush=True)
    arguments.cache.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for start in range(0, len(missing), BATCH_SIZE):
        batch = missing[start : start + BATCH_SIZE]
        retrieved_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        for attempt in range(3):
            try:
                returned = fetch_batch(batch)
                break
            except Exception as error:  # a transient failure is retried twice; then the batch is left for the next run
                print(f"batch at {start}: {error}; attempt {attempt + 1}", file=sys.stderr, flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            continue
        found = {str(item.get("pmid")): item for item in returned}
        for pmid in batch:
            item = found.get(pmid)
            rows.append({"pmid": pmid, "human": None if item is None else item.get("human"), "animal": None if item is None else item.get("animal"),
                         "molecular_cellular": None if item is None else item.get("molecular_cellular"), "year": None if item is None else item.get("year"),
                         "is_research_article": None if item is None else item.get("is_research_article"), "found": item is not None,
                         "retrieved_at": retrieved_at})
        if (start // BATCH_SIZE) % 10 == 9 or start + BATCH_SIZE >= len(missing):  # write as it goes, so a stopped run keeps what it read
            table = pd.concat([cached, pd.DataFrame(rows)], ignore_index=True) if len(cached) else pd.DataFrame(rows)
            table["pmid"] = table.pmid.astype(str)
            for column in ("human", "animal", "molecular_cellular"):
                table[column] = pd.to_numeric(table[column], errors="coerce")
            table.to_parquet(arguments.cache, index=False)
            print(f"{start + len(batch):,}/{len(missing):,} asked; {int(table.found.sum()):,} found in all", flush=True)
        time.sleep(arguments.pause)
    return 0


if __name__ == "__main__":
    sys.exit(main())
