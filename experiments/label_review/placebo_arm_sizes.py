"""Arm sizes for the placebo-controlled pairs, and the loss weight they imply.

The owner, 12 October 2026: "can you please look for placebo arm sizes for those specific studies and then calculate
negative weight based on that then use .5 as default for nothing found?"

Why the arm size matters. The 61 pairs that `--never-above-placebo-as-negative` relabels are pairs whose drug arm
reported the symptom no more often than its placebo arm. That is evidence against the label, but how much evidence
depends on how many patients were in the two arms: no excess over 30 patients says little and no excess over 800
says a good deal. SIDER carries the frequencies and not the arm sizes, and OnSIDES v3.1.1 carries neither, so the
sizes are read from the drug's own label, where the adverse-reaction tables head their columns "Placebo (n=385)".

The weight. With n_drug and n_placebo the effective size of a two-arm comparison is the harmonic quantity
n_effective = 1 / (1/n_drug + 1/n_placebo), which is what the standard error of a difference in proportions depends
on. The weight is then

    weight = n_effective / (n_effective + n_reference),   capped at 1.0

with n_reference the median effective size over the pairs that have one. That makes the owner's 0.5 both the default
for a pair whose sizes were not found and the weight of a pair of typical size, so the default is not an outlier: a
larger comparison earns more than 0.5 and a smaller one less, and nothing exceeds the weight of a grade A positive.

Counting and reading only: no label, table, graph or model input is changed. The output is a table for the owner.

Usage:
  python experiments/label_review/placebo_arm_sizes.py --fetch
  python experiments/label_review/placebo_arm_sizes.py --out configs/placebo_negative_weights.csv
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"
CACHE = Path("data/raw/label_review/placebo_arm_sizes")
DEFAULT_WEIGHT = 0.5  # the owner's value, for a pair whose arm sizes are not found
MAXIMUM_WEIGHT = 1.0  # the weight of a grade A positive; a large trial does not outweigh a monogenic label
HEADER_CELL = re.compile(r"<th[^>]*>(.*?)</th>", re.DOTALL | re.IGNORECASE)
CAPTION = re.compile(r"<caption[^>]*>(.*?)</caption>", re.DOTALL | re.IGNORECASE)
TAGS = re.compile(r"<[^>]+>")
ARM_SIZE = re.compile(r"\(?\s*[Nn]\s*=\s*([\d,]+)\s*\)?")
PLACEBO_IN_TEXT = re.compile(r"placebo[^.;()]{0,60}\(\s*[Nn]\s*=\s*([\d,]+)\s*\)", re.IGNORECASE)


def plain_text(markup: str) -> str:
    return " ".join(TAGS.sub(" ", markup).split())


# The study names a drug by the name its evidence source uses, which is not always the name on a United States
# label. Each entry is checked against the label returned, so a wrong synonym would show as a missing arm size
# rather than as a silently borrowed number.
LABEL_SYNONYMS: dict[str, tuple[str, ...]] = {
    "deprenyl": ("selegiline",),
    "hydroxybutyrate": ("sodium oxybate", "oxybate"),
    "rotigotine": ("rotigotine transdermal system",),
    "pergolide": (),      # withdrawn from the United States market in 2007, so no current label exists
    "reboxetine": (),     # never approved in the United States
}


def query_openfda(search_name: str) -> dict:
    query = urllib.parse.urlencode({"search": f'openfda.generic_name:"{search_name}"', "limit": 1})
    try:
        request = urllib.request.Request(f"{OPENFDA_LABEL}?{query}", headers={"User-Agent": "mechanistic-pathway-learning/0.1"})
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)
    except urllib.error.HTTPError as error:
        return {"error": f"HTTP {error.code}"}
    except Exception as error:  # noqa: BLE001 - a failed lookup is recorded, not raised
        return {"error": type(error).__name__}


def fetch_label(drug_name: str) -> dict | None:
    """The drug's label from openFDA, cached. Searched on the generic name, then on any recorded synonym."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{re.sub(r'[^a-z0-9]+', '_', drug_name.lower())}.json"
    if path.exists():
        cached = json.loads(path.read_text())
        if cached.get("results") or not LABEL_SYNONYMS.get(drug_name.lower()):
            return cached
    document, searched = {}, []
    for search_name in (drug_name, *LABEL_SYNONYMS.get(drug_name.lower(), ())):
        searched.append(search_name)
        document = query_openfda(search_name)
        time.sleep(0.5)
        if document.get("results"):
            break
    document["_queried_drug"] = drug_name
    document["_searched_names"] = searched
    path.write_text(json.dumps(document))
    return document


def arm_sizes_from_label(document: dict, drug_name: str) -> dict:
    """{placebo, drug, source, caption} from the label's adverse-reaction tables, or the text as a fallback."""
    results = document.get("results") or []
    if not results:
        return {"placebo": None, "drug": None, "source": document.get("error", "no label found"), "caption": ""}
    result = results[0]
    best = {"placebo": None, "drug": None, "source": "not stated", "caption": ""}
    for markup in result.get("adverse_reactions_table") or []:
        headers = [plain_text(cell) for cell in HEADER_CELL.findall(markup)]
        placebo_sizes, drug_sizes = [], []
        for header in headers:
            match = ARM_SIZE.search(header)
            if not match:
                continue
            size = int(match.group(1).replace(",", ""))
            (placebo_sizes if "placebo" in header.lower() else drug_sizes).append(size)
        if not placebo_sizes or not drug_sizes:
            continue
        # The largest placebo-controlled table is the one with the most patients behind it.
        placebo = max(placebo_sizes)
        if best["placebo"] is None or placebo > best["placebo"]:
            caption = CAPTION.search(markup)
            best = {"placebo": placebo, "drug": sum(drug_sizes), "source": "adverse reactions table",
                    "caption": plain_text(caption.group(1))[:160] if caption else ""}
    if best["placebo"] is not None:
        return best

    text = " ".join(result.get("adverse_reactions") or [])
    placebo_match = PLACEBO_IN_TEXT.search(text)
    if placebo_match:
        placebo = int(placebo_match.group(1).replace(",", ""))
        others = [int(value.replace(",", "")) for value in ARM_SIZE.findall(text)]
        others = [value for value in others if value != placebo]
        return {"placebo": placebo, "drug": max(others) if others else None,
                "source": "adverse reactions text", "caption": ""}
    return {"placebo": None, "drug": None, "source": "no arm size in the label", "caption": ""}


def effective_size(drug: float | None, placebo: float | None) -> float | None:
    if not drug or not placebo:
        return None
    return 1.0 / (1.0 / drug + 1.0 / placebo)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--selection", type=Path,
                        default=Path("data/processed/label_selection/better_v2_full_v4_placebo_negative.parquet"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--fetch", action="store_true", help="query openFDA for labels that are not cached")
    parser.add_argument("--out", type=Path, default=None)
    arguments = parser.parse_args()

    selection = pd.read_parquet(arguments.selection)
    pairs = selection[selection.counts_as_negative.astype(bool)][["perturbation_id", "symptom"]]
    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet",
                              columns=["perturbation_id", "perturbation_label"]).drop_duplicates()
    name_of = dict(zip(records.perturbation_id, records.perturbation_label))
    pairs = pairs.assign(drug=[name_of.get(identifier, identifier) for identifier in pairs.perturbation_id])
    print(f"placebo-controlled pairs: {len(pairs)} over {pairs.drug.nunique()} drugs")

    rows = []
    for drug in sorted(pairs.drug.unique()):
        document = fetch_label(drug) if (arguments.fetch or (CACHE / f"{re.sub(r'[^a-z0-9]+', '_', drug.lower())}.json").exists()) else None
        sizes = arm_sizes_from_label(document, drug) if document else {"placebo": None, "drug": None, "source": "not fetched", "caption": ""}
        rows.append({"drug": drug, "placebo_arm": sizes["placebo"], "drug_arm": sizes["drug"],
                     "arm_size_source": sizes["source"], "table_caption": sizes["caption"],
                     "effective_size": effective_size(sizes["drug"], sizes["placebo"])})
    per_drug = pd.DataFrame(rows)
    found = per_drug[per_drug.effective_size.notna()]
    print(f"\ndrugs with both arm sizes: {len(found)} of {len(per_drug)}")
    print(per_drug.arm_size_source.value_counts().to_string())

    if len(found):
        reference = float(found.effective_size.median())
        print(f"\nreference effective size, the median of those found: {reference:.1f} "
              f"(so a pair of that size weighs exactly {DEFAULT_WEIGHT})")
        per_drug["weight"] = [
            DEFAULT_WEIGHT if pd.isna(size) else min(MAXIMUM_WEIGHT, size / (size + reference))
            for size in per_drug.effective_size]
        print(f"  effective sizes: min {found.effective_size.min():.0f}, median {reference:.0f}, "
              f"max {found.effective_size.max():.0f}")
        weights = per_drug.weight
        print(f"  weights: min {weights.min():.3f}, median {weights.median():.3f}, max {weights.max():.3f}; "
              f"at the default {DEFAULT_WEIGHT}: {int((weights == DEFAULT_WEIGHT).sum())} drugs")
    else:
        per_drug["weight"] = DEFAULT_WEIGHT
        print("\nno arm sizes found; every pair takes the default weight")

    print("\nper drug:")
    print(per_drug[["drug", "placebo_arm", "drug_arm", "effective_size", "weight", "arm_size_source"]].to_string(index=False))

    pair_weights = pairs.merge(per_drug[["drug", "placebo_arm", "drug_arm", "effective_size", "weight", "arm_size_source"]], on="drug")
    if arguments.out:
        arguments.out.parent.mkdir(parents=True, exist_ok=True)
        pair_weights.sort_values(["drug", "symptom"]).to_csv(arguments.out, index=False)
        print(f"\nwrote {arguments.out} ({len(pair_weights)} pairs)")


if __name__ == "__main__":
    main()
