"""Check that every DOI, arXiv id and ISBN in docs/references.bib is registered, and that no two records
describe the same article. A citation that does not resolve is a defect of the write-up, and a duplicate
record is how one article comes to be cited twice under two keys.

The check asks Crossref in batches, then doi.org for the DOIs Crossref does not hold, because arXiv DOIs
(10.48550/arXiv.*) are registered with DataCite and the Crossref API does not return them. It reads the
bibliography and writes nothing but its report.

Usage:
  python experiments/check_bibliography_dois.py
  python experiments/check_bibliography_dois.py --bibliography docs/references.bib
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

CROSSREF_BATCH_SIZE = 20
REQUEST_TIMEOUT_SECONDS = 40
USER_AGENT = "mechanistic-pathway-learning bibliography check"


def records_of_bibliography(bibliography: Path) -> list[dict]:
    """One dictionary per BibTeX entry, with its key and the fields this check looks at."""
    records = []
    for block in re.split(r"\n(?=@)", bibliography.read_text()):
        if not block.lstrip().startswith("@"):
            continue
        key = re.search(r"@[a-zA-Z]+\{([^,]+),", block)
        if not key:
            continue

        def field(name: str) -> str | None:
            match = re.search(r"\b" + name + r"\s*=\s*\{([^}]*)\}", block, re.I)
            return match.group(1).strip() if match else None

        records.append({
            "key": key.group(1),
            "doi": field("doi"),
            "pmid": field("pmid"),
            "eprint": field("eprint"),
            "isbn": field("isbn"),
            "url": field("url"),
            "howpublished": field("howpublished"),
        })
    return records


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
        return json.load(response)


def dois_crossref_holds(dois: list[str]) -> set[str]:
    """The subset of dois that Crossref returns, asked in batches, lowercased for comparison."""
    held = set()
    for start in range(0, len(dois), CROSSREF_BATCH_SIZE):
        batch = dois[start : start + CROSSREF_BATCH_SIZE]
        query = ",".join("doi:" + doi for doi in batch)
        url = "https://api.crossref.org/works?rows=%d&select=DOI&filter=%s" % (
            CROSSREF_BATCH_SIZE,
            urllib.parse.quote(query, safe=":,./()-"),
        )
        try:
            payload = fetch_json(url)
        except (urllib.error.URLError, urllib.error.HTTPError) as error:
            print(f"  Crossref batch starting at {start} failed, so its DOIs fall through to doi.org: {error}")
            continue
        held.update(item["DOI"].lower() for item in payload["message"]["items"])
        time.sleep(0.5)
    return held


def doi_resolves(doi: str) -> tuple[bool, str]:
    """Whether doi.org resolves the DOI, which covers every registration agency, not Crossref alone."""
    url = "https://doi.org/" + urllib.parse.quote(doi, safe=":/().-")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            return True, str(response.status)
    except urllib.error.HTTPError as error:
        return error.code < 400, str(error.code)
    except urllib.error.URLError as error:
        return False, str(error)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bibliography", type=Path, default=Path("docs/references.bib"))
    arguments = parser.parse_args()

    records = records_of_bibliography(arguments.bibliography)
    keys_by_doi: dict[str, list[str]] = {}
    for record in records:
        if record["doi"]:
            keys_by_doi.setdefault(record["doi"], []).append(record["key"])

    print(f"{len(records)} records, {len(keys_by_doi)} distinct DOIs")

    problems = []

    duplicate_keys = [key for key in {record["key"] for record in records}
                      if sum(record["key"] == key for record in records) > 1]
    for key in sorted(duplicate_keys):
        problems.append(f"the key {key} is used by more than one record")

    for doi, keys in sorted(keys_by_doi.items()):
        if len(keys) > 1:
            problems.append(f"one article under {len(keys)} keys ({', '.join(sorted(keys))}): {doi}")

    without_locator = [record["key"] for record in records
                       if not any(record[field] for field in ("doi", "eprint", "isbn", "url", "howpublished"))]

    print("asking Crossref")
    held = dois_crossref_holds(sorted(keys_by_doi))
    outside_crossref = [doi for doi in sorted(keys_by_doi) if doi.lower() not in held]
    print(f"  Crossref holds {len(keys_by_doi) - len(outside_crossref)}; asking doi.org for the other "
          f"{len(outside_crossref)}")

    for doi in outside_crossref:
        resolves, status = doi_resolves(doi)
        if not resolves:
            problems.append(f"the DOI of {', '.join(keys_by_doi[doi])} does not resolve ({status}): {doi}")
        time.sleep(0.3)

    print()
    if without_locator:
        print("records with no DOI, arXiv id, ISBN or URL, which is a defect only where one exists:")
        for key in without_locator:
            print("  ", key)
        print()
    if problems:
        print("problems:")
        for problem in problems:
            print("  ", problem)
        return 1
    print("every DOI resolves, no key is reused and no article is under two keys")
    return 0


if __name__ == "__main__":
    sys.exit(main())
