"""Fetch and pin the fraction unbound in plasma (fu) database of Al-Qassabi and colleagues, which is stage 2 of
docs/plasma_protein_binding.md: the free fraction per drug, with the binding protein named.

It is the one open source found that carries all three things orosomucoid needs and that no source already pinned
here holds: the fraction unbound in a reference population, the same drug's fraction unbound in a population whose
binder concentration has moved (hepatic or renal impairment, age, inflammatory disease, ethnicity) and, per drug,
whether albumin (HSA), orosomucoid (AAG) or both is the major binding protein. Every row cites a PMID or DOI.

  dataset: https://doi.org/10.48420/25243138.v1, CC BY 4.0, University of Manchester
  paper:   Al-Qassabi et al., J Pharm Sci 113:1664-1673, 2024, https://doi.org/10.1016/j.xphs.2024.02.024

The file is taken through the figshare API rather than from a fixed download link, so its md5 can be checked against
the one figshare publishes for it; a mismatch fails the fetch rather than pinning a corrupt file. Raw files are never
committed (data/raw is gitignored); the pin is recorded in docs/data_sources.md, the coverage in
docs/fraction_unbound_coverage.md.

Usage:
  python experiments/fetch_fraction_unbound_database.py [--refresh]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

from experiments.fetch_node_descriptor_sources import USER_AGENT, download_and_pin

FIGSHARE_ARTICLE_ID = 25243138
PINNED_FILE_NAME = "fraction_unbound_database.xlsx"
EXPECTED_MD5 = "e1458c9b73dcf2174e36429df2c41e66"  # as figshare published it on 9 October 2026; the fetch checks it


def figshare_record(article_id: int = FIGSHARE_ARTICLE_ID) -> dict:
    """The figshare metadata of the dataset: its files with their md5, its DOI and its licence."""
    url = f"https://api.figshare.com/v2/articles/{article_id}"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def md5_of_file(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    parser.add_argument("--refresh", action="store_true", help="fetch again even when a pinned file is present")
    arguments = parser.parse_args()

    record = figshare_record()
    spreadsheets = [entry for entry in record["files"] if entry["name"].endswith(".xlsx")]
    if len(spreadsheets) != 1:
        raise SystemExit(f"expected one spreadsheet in figshare {FIGSHARE_ARTICLE_ID}, found {len(spreadsheets)}")
    spreadsheet = spreadsheets[0]
    if spreadsheet["supplied_md5"] != EXPECTED_MD5:
        raise SystemExit(f"figshare now publishes md5 {spreadsheet['supplied_md5']} for {spreadsheet['name']}, not "
                         f"{EXPECTED_MD5}: the dataset was revised, so read it before changing the expected value")

    pin = download_and_pin(arguments.raw_root, "plasma_protein_binding", PINNED_FILE_NAME,
                           spreadsheet["download_url"], arguments.refresh)
    pinned_path = Path(pin["file"])
    pinned_md5 = md5_of_file(pinned_path)
    if pinned_md5 != EXPECTED_MD5:
        raise SystemExit(f"{pinned_path} has md5 {pinned_md5}, not the {EXPECTED_MD5} figshare publishes for it")

    sidecar = pinned_path.with_suffix(pinned_path.suffix + ".figshare.json")
    sidecar.write_text(json.dumps({"article_id": FIGSHARE_ARTICLE_ID, "doi": record["doi"], "title": record["title"],
                                   "license": record["license"]["name"], "published_date": record["published_date"],
                                   "citation": record["citation"], "file_name": spreadsheet["name"],
                                   "supplied_md5": spreadsheet["supplied_md5"]}, indent=2) + "\n")
    print(f"{pinned_path}: {pin['bytes']:,} bytes, sha256 {pin['sha256'][:16]}, md5 verified against figshare")
    print(f"  {record['doi']} ({record['license']['name']}), published {record['published_date']}")
    print(f"  figshare name: {spreadsheet['name']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
