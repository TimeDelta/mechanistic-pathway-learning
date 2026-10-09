"""Read the carrier identity out of the cached FDA label sentences, and check it against the pinned measurements.

The plasma carriage edges of docs/plasma_binder_graph.md need to know which plasma protein holds a drug. The pinned
fraction-unbound database gives that for 27 of this study's drugs and names orosomucoid for 11 of them
(docs/fraction_unbound_coverage.md). This script asks how much further the drug labels reach, and whether they agree
with the database where both speak, which is the only check available on a new source of the same fact.

A drug is assigned a carrier only when its label sentences name exactly one of albumin and orosomucoid, outside a
negation. Sentences that name both are reported as "both" rather than resolved, and every assignment carries the
sentence and the label's set_id so a reader can check the claim rather than take it.

Reads the cache of experiments/fetch_fda_label_protein_binding.py and the pinned database. Writes
docs/label_plasma_binders.md. Adds no edge and changes no graph: whether these carriers become carriage edges is the
user's decision, and the confirmatory graph is not touched here.

Usage:
  python experiments/scope_label_plasma_binder_coverage.py
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

# The spellings a label uses for orosomucoid. "acid glycoprotein" covers "alpha-1-acid", "α 1 -acid" and the OCR
# variant "alpha-l-acid" with a letter l for the digit 1; some labels drop "acid" and write "α1-glycoprotein"; the
# acronyms are matched at a word boundary so they do not fire inside another word.
OROSOMUCOID_PATTERN = re.compile(r"acid glycoprotein|orosomucoid|(alpha|α)\s?-?\s?1?\s?-?\s?glycoprotein|\b(aag|agp)\b")
ALBUMIN_PATTERN = re.compile(r"albumin")
NEGATIONS = ("not bound", "not significantly bound", "no binding", "not appreciably bound", "independent of")
# a sentence saying the drug does not change another drug's binding is about the other drug's carrier, not this one's
ANOTHER_DRUG_S_BINDING = ("affect the binding of", "affected the binding of", "alter the binding of",
                          "altered the binding of", "displace the binding of")
OUTPUT_DOCUMENT = Path("docs/label_plasma_binders.md")
FRACTION_UNBOUND_SHEET_BINDERS = {"hsa": "albumin", "aag": "orosomucoid", "both": "both"}


def binders_named(sentence: str) -> set[str]:
    """Which plasma carriers one sentence names, empty when it names none, negates the binding or is about another drug."""
    lowered = sentence.lower()
    if any(negation in lowered for negation in NEGATIONS) or any(phrase in lowered for phrase in ANOTHER_DRUG_S_BINDING):
        return set()
    named = set()
    if OROSOMUCOID_PATTERN.search(lowered):
        named.add("orosomucoid")
    if ALBUMIN_PATTERN.search(lowered):
        named.add("albumin")
    return named


def carrier_of_drug(cached: dict) -> dict:
    """One drug's carrier reading from its cached label sentences."""
    drug = cached.get("drug", "")
    named_by_sentence = []
    for label in cached.get("labels", []):
        for sentence in label.get("sentences", []):
            named = binders_named(sentence)
            if named:
                names_the_drug = drug.lower() in sentence.lower()
                named_by_sentence.append((sorted(named), sentence, label.get("set_id", ""),
                                          label.get("effective_time", ""), names_the_drug))
    naming_the_drug = [entry for entry in named_by_sentence if entry[4]]
    # a sentence that names the drug is about this drug's binding; one that does not may be about a class or another
    # drug, so it is reported separately rather than mixed in
    used = naming_the_drug or named_by_sentence
    carriers = {binder for entry in used for binder in entry[0]}
    reading = {"drug": drug, "perturbation_id": cached.get("perturbation_id", ""),
               "labels": len(cached.get("labels", [])), "sentences_naming_a_carrier": len(named_by_sentence),
               "names_the_drug": bool(naming_the_drug),
               "carrier": "both" if len(carriers) > 1 else (next(iter(carriers)) if carriers else "")}
    if used:
        # the sentence the reading rests on: one that names the drug, then the most specific, then the fullest
        preferred = sorted(used, key=lambda item: (not item[4], len(item[0]), -len(item[1])))[0]
        reading.update({"sentence": preferred[1], "set_id": preferred[2], "effective_time": preferred[3]})
    return reading


def measured_carriers(database_path: Path) -> dict[str, str]:
    """The major binding protein of each drug in the pinned fraction-unbound database, lower-cased drug names."""
    sheets = pd.read_excel(database_path, sheet_name=None)
    carriers: dict[str, set[str]] = {}
    for name, frame in sheets.items():
        if not name[0].isdigit() or "Major Binding Protein" not in frame.columns:
            continue
        named = frame.rename(columns={"Major Binding Protein": "binder", "Drug": "drug"})
        for binder_cell, drug_cell in zip(named.binder, named.drug):
            binder = FRACTION_UNBOUND_SHEET_BINDERS.get(str(binder_cell).strip().lower())
            drug = str(drug_cell).strip().lower()
            if binder and drug and drug != "nan":
                carriers.setdefault(drug, set()).add(binder)
    return {drug: ("both" if len(binders) > 1 or "both" in binders else next(iter(binders)))
            for drug, binders in carriers.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/fda_labels"))
    parser.add_argument("--fraction-unbound-database", type=Path,
                        default=Path("data/raw/plasma_protein_binding/fraction_unbound_database.xlsx"))
    arguments = parser.parse_args()

    readings = [carrier_of_drug(json.loads(path.read_text())) for path in sorted(arguments.cache_dir.glob("*.json"))]
    table = pd.DataFrame(readings).sort_values("drug")
    measured = measured_carriers(arguments.fraction_unbound_database) if arguments.fraction_unbound_database.exists() else {}
    table["measured_carrier"] = table.drug.str.lower().map(measured).fillna("")
    with_a_carrier = table[table.carrier != ""]
    both_sources = with_a_carrier[with_a_carrier.measured_carrier != ""]
    agreeing = both_sources[both_sources.carrier == both_sources.measured_carrier]
    new_to_the_study = with_a_carrier[with_a_carrier.measured_carrier == ""]
    orosomucoid = with_a_carrier[with_a_carrier.carrier == "orosomucoid"]

    lines = [
        "# What the drug labels say about the plasma carrier, and how far it reaches",
        "",
        "Generated by `experiments/scope_label_plasma_binder_coverage.py`. Do not edit by hand.",
        "",
        "The user, 9 October 2026: \"i just want better coverage and i think it's possible to obtain in a valid way\". "
        "The pinned fraction-unbound database gives the carrier identity for 27 of this study's drugs. This is the "
        "second source of the same fact: the plasma-protein-binding sentences of each drug's FDA label, served by "
        "openFDA with the `set_id` and `effective_time` of the label version that carries them. A drug is assigned a "
        "carrier only when its sentences name exactly one of albumin and orosomucoid outside a negation.",
        "",
        f"- drugs with a cached label: **{len(table)}**",
        f"- drugs whose carrier sentence names the drug itself: **{int(with_a_carrier.names_the_drug.sum())}** of "
        f"{len(with_a_carrier)}; the rest rest on a sentence about a class or another drug and are marked below",
        f"- drugs whose label names a carrier: **{len(with_a_carrier)}** "
        f"({len(orosomucoid)} orosomucoid, {int((with_a_carrier.carrier == 'albumin').sum())} albumin, "
        f"{int((with_a_carrier.carrier == 'both').sum())} both)",
        f"- of those, drugs the fraction-unbound database does not hold: **{len(new_to_the_study)}**",
        f"- drugs both sources name: **{len(both_sources)}**, agreeing on **{len(agreeing)}**",
        "",
        "## Where the two sources disagree",
        "",
    ]
    disagreeing = both_sources[both_sources.carrier != both_sources.measured_carrier]
    if len(disagreeing):
        lines += ["| drug | label says | database says | the label sentence |", "| --- | --- | --- | --- |"]
        for row in disagreeing.itertuples(index=False):
            lines.append(f"| {row.drug} | {row.carrier} | {row.measured_carrier} | {getattr(row, 'sentence', '')[:240]} |")
    else:
        lines.append("Nowhere: every drug both sources name takes the same carrier from each.")
    lines += ["", "## The carrier each label names", "",
              "| drug | carrier | database | names the drug | label set_id | effective | sentence |",
              "| --- | --- | --- | --- | --- | --- | --- |"]
    for row in with_a_carrier.itertuples(index=False):
        lines.append(f"| {row.drug} | {row.carrier} | {row.measured_carrier or '-'} | "
                     f"{'yes' if row.names_the_drug else 'NO'} | `{getattr(row, 'set_id', '')[:8]}` | "
                     f"{getattr(row, 'effective_time', '')} | {getattr(row, 'sentence', '')[:240]} |")
    silent = table[table.carrier == ""]
    lines += ["", "## The drugs whose labels name no carrier", "",
              f"{len(silent)} drugs: " + ", ".join(sorted(silent.drug)) + ".", "",
              "A label states the carrier when it matters clinically, which is why the drugs it names are not a random "
              "sample of the study's drugs and why this source cannot be read as a measurement of how many drugs "
              "orosomucoid carries. It is a per-drug statement of carrier identity, which is what a carriage edge "
              "needs, and nothing here is a free fraction."]
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")
    print(f"cached drugs {len(table)}; carrier named for {len(with_a_carrier)} "
          f"({len(orosomucoid)} orosomucoid); new to the study {len(new_to_the_study)}; "
          f"both sources {len(both_sources)} agreeing {len(agreeing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
