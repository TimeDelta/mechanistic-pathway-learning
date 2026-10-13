#!/usr/bin/env python
"""Does a drug's own label state that it acts on the central nervous system? Read from openFDA, quoted, per drug.

Why this exists. The central-action rule of experiments/apply_central_action_rule.py has two measured halves, a
brain-expressed mechanism target and a measured Kp,uu,brain, and over the 924-drug pool those two admit metformin,
warfarin, omeprazole, etoposide and six NSAIDs while leaving carvedilol out (docs/label_source_review.md 29.2). The
reason is that neither half is a claim about the drug: almost every drug has a target expressed somewhere in the brain,
and having a Kp,uu,brain means the drug was used as a reference compound in a brain-penetration study. The half that was
missing is the one the ATC list carried by hand — the drug's established pharmacology — and this script supplies it as
evidence rather than as a group name, to the standard the six named drugs of configs/admitted_atc_prefixes_draft.csv
already met: the sentence of the label that states the central action, quoted, with the set identifier it was read from.

What is read, and what is deliberately not. Only the pharmacology sections: clinical pharmacology, mechanism of action,
pharmacodynamics, description. The adverse-reaction, warning and indication sections are NOT read, because a drug
admitted for the psychiatric events on its label would select the drug set on the outcome the study predicts, which is
the objection section 9.4 had to guard against by hand. A sentence about a metabolite, an excipient or another drug
class is not a statement about this drug, so every match is kept with its sentence for a reader to judge; this script
reports and does not decide.

What a verdict means. "states a central action" is a sentence in a pharmacology section that places the drug, its
action or its site in the central nervous system or across the blood-brain barrier. "states that it does not reach the brain" is the label saying so, which is positive
evidence for leaving the drug out. "only an adverse-event sentence" means the single match was an adverse-event phrase
that leaked into a pharmacology section, which is not a mechanism claim. "no such sentence" means the sections were read
and none was found, which is a finding. "no label section read" means openFDA returned no
pharmacology text for the name, which is not.

What this screen cannot do, measured on a control. Metoclopramide is already in
configs/admitted_atc_prefixes_draft.csv with the quoted sentence "The antiemetic properties of metoclopramide appear to
be a result of its antagonism of central and peripheral dopamine receptors", read by hand from openFDA set
0034a97a-fb05-418a-a234-de278907cd2a. This screen reads the labels openFDA returns first for a generic name, so it can
miss a sentence that sits in another label of the same drug: at a limit of 3 labels it returned "no such sentence" for
metoclopramide. The limit is 8 for that reason, and the screen is still a finder and not an oracle. A negative verdict
on a drug that matters should be checked against the drug's other labels before it decides anything.

Measurement only: no evidence table, configuration or model is changed. Responses are cached under
data/raw/label_review/openfda_central_action so a rerun reads no network.

Usage:
  python experiments/read_label_central_action.py --drugs propranolol metoprolol atenolol
  python experiments/read_label_central_action.py --from-rule-table configs/centrally_acting_drugs.csv --verdict admitted
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

LABEL_ENDPOINT = "https://api.fda.gov/drug/label.json"
CACHE_DIRECTORY = Path("data/raw/label_review/openfda_central_action")
OUTPUT_TABLE = Path("configs/label_central_action.csv")
OUTPUT_DOCUMENT = Path("docs/label_central_action.md")
USER_AGENT = "mechanistic-pathway-learning/0.1 (central action rule)"
# The pharmacology sections only. Adverse reactions, warnings and indications are not read; see the module docstring.
PHARMACOLOGY_FIELDS = ("mechanism_of_action", "clinical_pharmacology", "pharmacodynamics", "description")
# A sentence that places the drug, its action or its site inside the central nervous system or across the barrier.
CENTRAL_ACTION = re.compile(r"\bcentral nervous system\b|\bblood[- ]brain barrier\b|\bCNS\b|\bcentrally\b|\bcentral action\b"
                            r"|\bin the brain\b|\bbrain receptors?\b|\bcerebral\b|\bintracerebral\b", re.IGNORECASE)
# A sentence can carry a barrier or brain word and say the drug does NOT get there. The first screen written here
# counted four such sentences as positives - labetalol's "Only negligible amounts of the drug crossed the blood-brain
# barrier", nadolol's "only minimal amounts ... detected in the brain", sotalol's "crosses the blood brain barrier
# poorly" and montelukast's "minimal distribution across the blood-brain barrier" - because it looked for a negative
# word immediately before the verb and these put it elsewhere in the sentence. So the negative is matched anywhere in
# the sentence, by word rather than by position.
READS_AS_NOT_CENTRAL = re.compile(r"\bnegligible\b|\bminimal(?:ly)?\b|\bpoorly\b|\bdoes not\b|\bdo not\b|\bnot (?:readily|detect)"
                                  r"|\blower than\b|\blimited (?:distribution|penetration|entry)\b|\blittle or no\b|\bundetectable\b"
                                  r"|\bdoes-not\b|\bwithout (?:crossing|entering)\b", re.IGNORECASE)
# An adverse-event phrase that has leaked into a pharmacology section is not a statement of mechanism: ibuprofen's only
# match was "milder gastrointestinal ... and CNS side effects". Admitting a drug on such a sentence would be admitting
# it for its adverse events, which is what reading only the pharmacology sections was meant to prevent.
READS_AS_AN_ADVERSE_EVENT = re.compile(r"\bside effects?\b|\badverse\b|\btoxicit|\bdepression of\b.{0,24}\bfunction\b", re.IGNORECASE)
SENTENCE = re.compile(r"(?<=[.;])\s+(?=[A-Z(])")


def label_payload(drug_name: str, limit: int = 8) -> dict:
    """The openFDA label response for a generic name, cached on disk by name."""
    CACHE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIRECTORY / f"{re.sub(r'[^a-z0-9]+', '_', drug_name.lower())}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    url = f"{LABEL_ENDPOINT}?" + urllib.parse.urlencode({"search": f'openfda.generic_name:"{drug_name}"', "limit": limit})
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=60) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:  # 404 is openFDA's answer for a name it does not hold
        payload = {"results": [], "error": f"HTTP {error.code}"}
    except Exception as error:  # noqa: BLE001 - the reason is recorded rather than raised, so one name cannot stop a sweep
        payload = {"results": [], "error": f"{type(error).__name__}: {error}"}
    cache_path.write_text(json.dumps(payload, indent=1))
    time.sleep(0.4)  # openFDA asks for a modest rate without a key
    return payload


def central_action_sentences(payload: dict) -> tuple[list[dict], int]:
    """(the sentences that state a central action, how many pharmacology sections were read)."""
    found, sections_read = [], 0
    for result in payload.get("results", []) or []:
        set_id = result.get("set_id", "")
        for field in PHARMACOLOGY_FIELDS:
            text = " ".join(result.get(field, []) or [])
            if not text:
                continue
            sections_read += 1
            for sentence in SENTENCE.split(text):
                sentence = " ".join(sentence.split())
                if CENTRAL_ACTION.search(sentence) and len(sentence) < 600:
                    found.append({"section": field, "set_id": set_id, "sentence": sentence,
                                  "reads_as_not_central": bool(READS_AS_NOT_CENTRAL.search(sentence)),
                                  "reads_as_an_adverse_event": bool(READS_AS_AN_ADVERSE_EVENT.search(sentence))})
    return found, sections_read


def verdicts(drug_names: list[str]) -> pd.DataFrame:
    rows = []
    for drug_name in drug_names:
        payload = label_payload(drug_name)
        sentences, sections_read = central_action_sentences(payload)
        positive = [row for row in sentences if not row["reads_as_not_central"] and not row["reads_as_an_adverse_event"]]
        negative = [row for row in sentences if row["reads_as_not_central"]]
        adverse = [row for row in sentences if row["reads_as_an_adverse_event"] and not row["reads_as_not_central"]]
        # A negative sentence outranks a positive one. Nadolol showed why: its own sentence says "only minimal amounts of
        # nadolol were detected in the brain", while the class boilerplate every beta blocker's label carries ("the
        # mechanism of the antihypertensive effects ... has not been established; factors may include ... vasomotor
        # centers in the brain") reads positive. The drug-specific statement is the informative one, and a label that
        # says the drug does not reach the brain is evidence for leaving it out rather than a tie.
        verdict = ("no label section read" if not sections_read else
                   "states that it does not reach the brain" if negative else
                   "states a central action" if positive else
                   "only an adverse-event sentence" if adverse else "no such sentence")
        best = (negative or positive or adverse or [{}])[0]
        rows.append({"drug": drug_name, "verdict": verdict, "sections_read": sections_read,
                     "matching_sentences": len(positive), "sentences_reading_as_not_central": len(negative),
                     "adverse_event_sentences": len(adverse),
                     "section": best.get("section", ""), "set_id": best.get("set_id", ""), "quote": best.get("sentence", ""),
                     "error": payload.get("error", "")})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--drugs", nargs="*", default=[], help="generic names to read")
    parser.add_argument("--from-rule-table", type=Path, default=None, help="configs/centrally_acting_drugs.csv, whose drug column names the candidates")
    parser.add_argument("--verdict", nargs="*", default=None, help="only the rows of that table with these verdicts")
    parser.add_argument("--output-table", type=Path, default=OUTPUT_TABLE)
    parser.add_argument("--output-document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    drug_names = list(arguments.drugs)
    if arguments.from_rule_table is not None:
        table = pd.read_csv(arguments.from_rule_table)
        if arguments.verdict:
            table = table[table.verdict.isin(arguments.verdict)]
        drug_names += [name for name in table.drug.astype(str) if name not in drug_names]
    if not drug_names:
        parser.error("name some drugs with --drugs or --from-rule-table")

    table = verdicts(drug_names)
    arguments.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.output_table, index=False)
    print(table.verdict.value_counts().to_string())
    write_document(arguments.output_document, table)


def write_document(path: Path, table: pd.DataFrame) -> None:
    counts = table.verdict.value_counts()
    provenance = ("Read from openFDA by `experiments/read_label_central_action.py`, pharmacology sections only "
                  "(mechanism of action, clinical pharmacology, pharmacodynamics, description). The adverse-reaction, "
                  "warning and indication sections are deliberately not read: a drug admitted for the psychiatric events "
                  "on its label would select the drug set on the outcome the study predicts. A positive verdict is a "
                  "shortlist to read, not a decision: class boilerplate makes every beta blocker positive, and a "
                  "combination-product label can put a drug's only match on a sentence about the other ingredient.")
    lines = ["# What a drug's own label says about acting centrally", "", provenance, "",
             f"{len(table)} drugs read.", "", "| verdict | drugs |", "| --- | --- |"]
    lines += [f"| {verdict} | {count} |" for verdict, count in counts.items()]
    for verdict in ("states a central action", "states that it does not reach the brain", "only an adverse-event sentence",
                    "no such sentence", "no label section read"):
        rows = table[table.verdict == verdict]
        if not len(rows):
            continue
        lines += ["", f"## {verdict}", ""]
        if verdict in ("no such sentence", "no label section read", "only an adverse-event sentence"):
            lines += [", ".join(sorted(rows.drug)), ""]
            continue
        lines += ["| drug | section | quote | openFDA set |", "| --- | --- | --- | --- |"]
        for row in rows.sort_values("drug").itertuples(index=False):
            quote = row.quote.replace("|", "\\|")
            lines.append(f"| {row.drug} | {row.section} | {quote} | {row.set_id} |")
    lines.append("")
    path.write_text("\n".join(lines))
    print(f"wrote {path} and {OUTPUT_TABLE}")


if __name__ == "__main__":
    main()
