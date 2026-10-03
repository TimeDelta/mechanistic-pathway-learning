"""Publication dates of PubMed identifiers through E-utilities esummary, cached and rate-limited (design section 6.1).

Why: the PubTator3 bulk relation file carries no date and the time split needs one per report. esummary returns a
pubdate such as "2022", "2022 Mar 29" or "2021 Jan-Feb"; it is parsed to an ISO date with the first day of the
month or year when the finer part is missing, so a report dated "2022" sorts before any cutoff inside 2022 and
after one at the end of 2021. Identifiers are sorted and sent in fixed batches so the cache key of a batch is
deterministic and a rerun makes no request.
"""
from __future__ import annotations

import json
import re
import urllib.parse
from collections.abc import Iterable

from mechanistic_pathway_learning.evidence.cached_http import CachedRateLimitedClient

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ESUMMARY_BATCH_SIZE = 200
MONTH_NUMBERS = {name: number for number, name in enumerate(("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), start=1)}
PUBDATE_PATTERN = re.compile(r"^(?P<year>\d{4})(?:\s+(?P<month>[A-Za-z]{3})[A-Za-z]*(?:-[A-Za-z]+)?)?(?:\s+(?P<day>\d{1,2}))?")


def parse_pubdate(pubdate: str | None) -> str | None:
    """ISO date from an esummary pubdate; None when no four-digit year leads the string."""
    if not pubdate:
        return None
    match = PUBDATE_PATTERN.match(pubdate.strip())
    if match is None:
        return None
    year = int(match.group("year"))
    month = MONTH_NUMBERS.get((match.group("month") or "").lower()[:3], 1)
    day = int(match.group("day") or 1)
    if not 1 <= day <= 31:
        day = 1
    return f"{year:04d}-{month:02d}-{day:02d}"


def fetch_publication_dates(client: CachedRateLimitedClient, pmids: Iterable[str], batch_size: int = ESUMMARY_BATCH_SIZE) -> tuple[dict[str, str], dict]:
    """PMID -> ISO date for every identifier esummary knows; the second value counts batches, hits and misses."""
    unique = sorted({str(pmid) for pmid in pmids if str(pmid).strip()}, key=lambda value: (len(value), value))
    dates: dict[str, str] = {}
    counts = {"pmids_requested": len(unique), "batches": 0, "pmids_dated": 0, "pmids_without_date": 0, "pmids_not_returned": 0}
    for start in range(0, len(unique), batch_size):
        batch = unique[start:start + batch_size]
        body = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "id": ",".join(batch)}).encode("utf-8")
        document = json.loads(client.fetch_text(ESUMMARY_URL, post_body=body))
        counts["batches"] += 1
        result = document.get("result", {})
        for pmid in batch:
            summary = result.get(pmid)
            if not isinstance(summary, dict) or "error" in summary:
                counts["pmids_not_returned"] += 1
                continue
            iso_date = parse_pubdate(summary.get("pubdate")) or parse_pubdate(summary.get("epubdate"))
            if iso_date is None:
                counts["pmids_without_date"] += 1
            else:
                dates[pmid] = iso_date
                counts["pmids_dated"] += 1
    return dates, counts
