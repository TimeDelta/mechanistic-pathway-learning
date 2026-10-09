"""Cache the plasma-protein-binding sentences of each study drug's FDA label, with the label version that carries them.

Why this source. The carrier identity of a drug, which plasma protein holds it, reaches only 27 of this study's drugs
through the pinned fraction-unbound database (docs/fraction_unbound_coverage.md), and the plasma carriage edges of
docs/plasma_binder_graph.md need that identity per drug. A drug label states it when it matters clinically: methadone's
says "In plasma, methadone is predominantly bound to alpha-1-acid glycoprotein (85% to 90%)". Labels are a source
family this study already uses, because OnSIDES is built from them, and openFDA serves them with a `set_id` and an
`effective_time`, so the exact label version behind a sentence is recorded and the sentence can be quoted.

What this writes. One JSON per drug under --cache-dir, holding the query, the retrieval date, and per label its
set_id, effective_time and the sentences of its pharmacology sections that mention protein binding. Nothing is
interpreted here: which binder a sentence names is read by experiments/scope_label_plasma_binder_coverage.py, so a
change to the reading does not refetch. The cache is raw downloaded data, so it goes in its own directory and no code
is run from it.

Idempotent and resumable: a drug whose cache file exists is skipped unless --refetch. Rate limited to openFDA's
unauthenticated allowance.

Usage:
  python experiments/fetch_fda_label_protein_binding.py
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

LABEL_ENDPOINT = "https://api.fda.gov/drug/label.json"
PHARMACOLOGY_SECTIONS = ("clinical_pharmacology", "pharmacokinetics", "description", "dosage_and_administration")
SENTENCE_BREAK = re.compile(r"(?<=[.;])\s+")
BINDING_WORDS = ("protein binding", "protein-binding", "bound to", "binds to", "glycoprotein", "albumin")
SECONDS_BETWEEN_REQUESTS = 0.34  # openFDA allows 240 unauthenticated requests a minute
RETRY_WAITS = (2, 4, 8, 16)


def binding_sentences(text: str) -> list[str]:
    """The sentences of one label section that say something about plasma protein binding."""
    kept = []
    for sentence in SENTENCE_BREAK.split(text):
        lowered = sentence.lower()
        if any(word in lowered for word in BINDING_WORDS) and ("protein" in lowered or "albumin" in lowered):
            kept.append(" ".join(sentence.split()))
    return kept


def fetch_label(drug_label: str, limit: int) -> tuple[str, dict]:
    """openFDA's labels for one generic name, or an empty payload when it has none."""
    query = urllib.parse.quote(f'openfda.generic_name:"{drug_label}"')
    url = f"{LABEL_ENDPOINT}?search={query}&limit={limit}"
    for attempt, wait in enumerate((0,) + RETRY_WAITS):
        if wait:
            time.sleep(wait)
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return url, json.loads(response.read())
        except urllib.error.HTTPError as error:
            if error.code == 404:  # openFDA answers a search with no match with 404
                return url, {}
            if error.code not in (429, 500, 502, 503, 504) or attempt == len(RETRY_WAITS):
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == len(RETRY_WAITS):
                raise
    return url, {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/fda_labels"))
    parser.add_argument("--labels-per-drug", type=int, default=5)
    parser.add_argument("--refetch", action="store_true")
    arguments = parser.parse_args()

    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = (records[records.perturbation_type == "drug"][["perturbation_id", "perturbation_label"]]
             .drop_duplicates().sort_values("perturbation_label"))
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    fetched, skipped, without_a_label, failed = 0, 0, 0, []
    for row in drugs.itertuples(index=False):
        cache_path = arguments.cache_dir / f"{row.perturbation_id}.json"
        if cache_path.exists() and not arguments.refetch:
            skipped += 1
            continue
        try:
            url, payload = fetch_label(row.perturbation_label, arguments.labels_per_drug)
        except Exception as error:  # noqa: BLE001 - the drug is recorded as unfetched and the loop goes on
            failed.append((row.perturbation_label, str(error)))
            continue
        labels = [{"set_id": result.get("set_id", ""), "effective_time": result.get("effective_time", ""),
                   "sentences": sorted({sentence for section in PHARMACOLOGY_SECTIONS
                                        for text in result.get(section, []) for sentence in binding_sentences(text)})}
                  for result in payload.get("results", [])]
        cache_path.write_text(json.dumps({"perturbation_id": row.perturbation_id, "drug": row.perturbation_label,
                                          "query": url, "retrieved": date.today().isoformat(), "labels": labels}, indent=1))
        fetched += 1
        without_a_label += 0 if labels else 1
        time.sleep(SECONDS_BETWEEN_REQUESTS)
    print(f"drugs: {len(drugs)}; fetched {fetched}, already cached {skipped}, no openFDA label {without_a_label}")
    for drug, error in failed:
        print(f"unfetched: {drug}: {error}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
