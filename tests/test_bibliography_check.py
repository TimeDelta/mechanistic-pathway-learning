"""The bibliography checker's parser and the duplicate rules it enforces
(experiments/check_bibliography_dois.py). The network checks are not exercised here; the parser is, because
every finding the checker reports rests on it, and so is the live bibliography's shape."""
import re
from pathlib import Path

from experiments.check_bibliography_dois import records_of_bibliography

BIBLIOGRAPHY = Path("docs/references.bib")

SAMPLE = """% a comment line between records
@article{one_1999,
  author = {A. Author},
  title = {A title},
  pmid = {12345},
  doi = {10.1000/one}
}

@book{two_2000,
  author = {B. Author},
  title = {A book},
  isbn = {0-123-45678-9}
}

@inproceedings{three_2001,
  author = {C. Author},
  title = {A preprint},
  eprint = {1712.01312},
  archiveprefix = {arXiv}
}
"""


def test_the_parser_reads_the_key_and_the_locator_fields(tmp_path):
    bibliography = tmp_path / "sample.bib"
    bibliography.write_text(SAMPLE)
    records = {record["key"]: record for record in records_of_bibliography(bibliography)}
    assert set(records) == {"one_1999", "two_2000", "three_2001"}
    assert records["one_1999"]["doi"] == "10.1000/one"
    assert records["one_1999"]["pmid"] == "12345"
    assert records["two_2000"]["isbn"] == "0-123-45678-9"
    assert records["two_2000"]["doi"] is None, "a field a record does not carry reads as absent, not as empty"
    assert records["three_2001"]["eprint"] == "1712.01312"


def test_the_live_bibliography_has_no_duplicate_key_and_no_article_under_two_keys():
    """The two defects the checker exists for; both were in the file before 9 October 2026."""
    records = records_of_bibliography(BIBLIOGRAPHY)
    keys = [record["key"] for record in records]
    assert len(keys) == len(set(keys)), "a key is reused"
    keys_by_doi = {}
    for record in records:
        if record["doi"]:
            keys_by_doi.setdefault(record["doi"].lower(), []).append(record["key"])
    shared = {doi: keys for doi, keys in keys_by_doi.items() if len(keys) > 1}
    assert not shared, f"one article under several keys: {shared}"


def test_every_live_record_carries_a_locator_or_is_the_known_exception():
    """Kacser and Burns 1973 predates DOI registration, so its note names the 1995 reprint's DOI instead."""
    without_locator = {record["key"] for record in records_of_bibliography(BIBLIOGRAPHY)
                       if not any(record[field] for field in ("doi", "eprint", "isbn", "url", "howpublished"))}
    assert without_locator == {"ref32_kacser1973"}, without_locator


def test_every_live_doi_has_the_shape_of_a_doi():
    """A typed or truncated DOI fails here without a network call."""
    for record in records_of_bibliography(BIBLIOGRAPHY):
        if not record["doi"]:
            continue
        assert re.fullmatch(r"10\.\d{4,9}/\S+", record["doi"]), (record["key"], record["doi"])
