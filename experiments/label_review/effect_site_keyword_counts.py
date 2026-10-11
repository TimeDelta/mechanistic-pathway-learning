"""PubMed keyword counts per study drug for three kinds of effect-site data: a receptor occupancy measured by PET or
SPECT, an effect-compartment equilibration rate, and an unbound brain concentration.

Part of the question of 10 October 2026 whether a prediction at the peak of a drug's concentration at its site of action
can be supplied with data (docs/label_source_review.md, section 10). A keyword count says that an abstract names the
drug together with the words of the query. It does not say that the article measured the quantity for that drug, so
every count here is an upper bound on the drugs with data, and the script writes the first PubMed identifiers of each
query so that a count can be checked by reading.

Reads data/processed/label_review/drug_names.csv (the study drugs under the current rule, with the three approved label
changes and with the draft list of admitted ATC groups). Writes the responses of the NCBI E-utilities under --cache-dir
(raw downloaded data; no code is run from it) and one row per drug to --output. Counting only: no label, table or
configuration of the study is changed.

Usage:
  python experiments/label_review/effect_site_keyword_counts.py
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
SECONDS_BETWEEN_REQUESTS = 0.4  # the E-utilities allow three requests a second without a key
LABELS_NOT_SEARCHED = {"delta", "lu", "implanon"}  # study labels that are not a drug name on their own; the RxNorm name is searched
EXTRA_NAMES = {"hydroxybutyrate": ["gamma-hydroxybutyrate", "sodium oxybate"], "valproate": ["valproic acid"], "acetaminophen": ["paracetamol"],
               "delta": ["tetrahydrocannabinol"], "dothiepin": ["dosulepin"], "meperidine": ["pethidine"], "levodopa": ["L-DOPA"]}
OCCUPANCY = '(occupancy[tiab]) AND ("positron emission tomography"[tiab] OR PET[tiab] OR SPECT[tiab] OR "single photon emission"[tiab])'
QUERIES = {
    "occupancy_by_pet_or_spect": OCCUPANCY,
    "occupancy_with_time_words": OCCUPANCY + ' AND ("time course"[tiab] OR "time-course"[tiab] OR "hours after"[tiab] OR "h after"[tiab] OR "hr after"[tiab] OR "half-life"[tiab] OR kinetic*[tiab] OR "duration of"[tiab])',
    "effect_compartment_rate": '("effect compartment"[tiab] OR "effect-compartment"[tiab] OR ke0[tiab] OR "effect site"[tiab] OR "effect-site"[tiab] OR biophase[tiab] OR "equilibration half-life"[tiab])',
    "unbound_brain_concentration": '("Kp,uu"[tiab] OR "unbound brain"[tiab] OR "brain unbound"[tiab] OR "unbound brain-to-plasma"[tiab] OR "brain interstitial fluid"[tiab] OR (microdialysis[tiab] AND brain[tiab] AND unbound[tiab]))',
}


def search_names(row) -> list[str]:
    names = [name for name in str(row.rxnorm_names).split(";") if name]
    label = str(row.perturbation_label)
    if label.lower() not in LABELS_NOT_SEARCHED and label.lower() not in {name.lower() for name in names}:
        names.append(label)
    names += EXTRA_NAMES.get(label.lower(), [])
    return sorted(set(names), key=str.lower)


def esearch(term: str) -> dict:
    url = ESEARCH + "?" + urllib.parse.urlencode({"db": "pubmed", "term": term, "retmode": "json", "retmax": 8, "sort": "relevance"})
    for wait in (0, 2, 5, 15, 40):
        if wait:
            time.sleep(wait)
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                result = json.loads(response.read())["esearchresult"]
            return {"term": term, "count": int(result["count"]), "pmids": result.get("idlist", [])}
        except (urllib.error.URLError, TimeoutError, KeyError, json.JSONDecodeError):
            continue
    raise RuntimeError(f"no answer from the E-utilities for {term!r}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--drug-names", type=Path, default=Path("data/processed/label_review/drug_names.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/label_review/pubmed_keyword_counts"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/effect_site_keyword_counts.csv"))
    arguments = parser.parse_args()
    drugs = pd.read_csv(arguments.drug_names).fillna("")
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for position, drug in enumerate(drugs.itertuples(index=False)):
        names = search_names(drug)
        name_clause = "(" + " OR ".join(f'"{name}"[tiab]' for name in names) + ")"
        cache_path = arguments.cache_dir / f"{drug.perturbation_id.replace(':', '_')}.json"
        cached = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        for query_name, query in QUERIES.items():
            if query_name not in cached or cached[query_name]["term"] != f"{name_clause} AND {query}":
                cached[query_name] = esearch(f"{name_clause} AND {query}")
                time.sleep(SECONDS_BETWEEN_REQUESTS)
        cache_path.write_text(json.dumps(cached, indent=1))
        rows.append({"perturbation_id": drug.perturbation_id, "drug": drug.perturbation_label, "searched_names": ";".join(names),
                     **{query_name: cached[query_name]["count"] for query_name in QUERIES},
                     **{f"{query_name}_first_pmids": ";".join(cached[query_name]["pmids"]) for query_name in QUERIES}})
        if position % 25 == 0:
            print(position, drug.perturbation_label, {query_name: cached[query_name]["count"] for query_name in QUERIES}, flush=True)
    table = pd.DataFrame(rows).merge(drugs[["perturbation_id", "in_current", "in_approved", "in_approved_with_list"]], on="perturbation_id")
    table.to_csv(arguments.output, index=False)
    for set_name in ("in_current", "in_approved", "in_approved_with_list"):
        of_set = table[table[set_name].astype(bool)]
        print(set_name, len(of_set), {query_name: int((of_set[query_name] > 0).sum()) for query_name in QUERIES})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
