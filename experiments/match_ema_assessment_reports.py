"""Which of the study's drugs have a European public assessment report, and where it is.

Decision 9.28, the owner on 14 October 2026: "include the EPARs". An EPAR carries the clinical pharmacology a
modern US label states only in fragments across its PLR subsections, which is the gap part (a)'s matcher falls into,
and it often carries brain-distribution detail for drugs with no published Kp,uu,brain, which is part (c)'s gap.

The EMA medicines table is the entry point rather than guessing a brand name per drug: it lists every centrally
authorised product with its active substance, its marketing authorisation date and the address of its EPAR page.
Matching is on the active substance and the INN, lowercased, because the study keys drugs by generic name.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import pandas as pd

EMA_TABLE_HEADER_ROWS_TO_SKIP = 8

NEWER_DRUGS_WANTED = (
    "cariprazine", "brexpiprazole", "lumateperone", "bupropion", "agomelatine", "desvenlafaxine", "clobazam",
    "bromocriptine", "selegiline", "vortioxetine", "escitalopram", "topiramate", "guanfacine", "viloxazine",
    "levomilnacipran", "zuranolone", "brexanolone", "tianeptine", "flupentixol", "zuclopenthixol", "sertindole",
    "solriamfetol", "trimipramine", "mianserin", "dextromethorphan", "sodium oxybate", "lurasidone", "pitolisant",
    "istradefylline", "brivaracetam", "amisulpride", "moclobemide", "milnacipran",
)


def normalised_substance_names(cell: object) -> set[str]:
    """One EMA row can name several substances, separated by commas, semicolons or newlines."""
    if not isinstance(cell, str):
        return set()
    parts = re.split(r"[;,\n]", cell)
    return {part.strip().lower() for part in parts if part.strip()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ema-table", type=Path, default=Path("data/raw/ema/medicines-output-medicines-report_en.xlsx"))
    parser.add_argument("--label-verdicts", type=Path, default=Path("configs/label_central_action.csv"))
    parser.add_argument("--psychiatric-additions", type=Path, default=Path("configs/psychiatric_drug_additions_v4_draft.csv"))
    parser.add_argument("--output", type=Path, default=Path("configs/ema_assessment_reports.csv"))
    parser.add_argument("--summary", type=Path, default=Path("docs/ema_assessment_report_coverage.json"))
    arguments = parser.parse_args()

    ema = pd.read_excel(arguments.ema_table, skiprows=EMA_TABLE_HEADER_ROWS_TO_SKIP)
    ema = ema[ema.Category.astype(str).str.strip().str.lower() == "human"]
    inn_column = next(column for column in ema.columns if "non-proprietary" in str(column))
    rows_by_substance: dict[str, list] = {}
    for values in ema.itertuples(index=False):
        record = dict(zip(ema.columns, values))
        names = normalised_substance_names(record.get("Active substance"))
        names |= normalised_substance_names(record.get(inn_column))
        for name in names:
            rows_by_substance.setdefault(name, []).append(record)

    drugs_wanted: dict[str, str] = {}
    verdicts = pd.read_csv(arguments.label_verdicts)
    for row in verdicts.itertuples():
        if row.verdict in ("no such sentence", "no label section read", "only an adverse-event sentence"):
            drugs_wanted[str(row.drug).lower()] = f"item 31(a): {row.verdict}"
    # The psychiatric additions table is keyed by perturbation id and carries no drug name, so the names come from
    # this list instead: the modern drugs decision 1.3(b) found to lack an exposure chain, plus the additions of
    # decision 9.7 and the drugs whose only occupancy evidence is a PET study.
    for name in NEWER_DRUGS_WANTED:
        drugs_wanted.setdefault(name, "newer drug lacking a Kp,uu,brain or a central-exposure chain")

    matched_rows, unmatched = [], []
    for drug, reason in sorted(drugs_wanted.items()):
        candidates = rows_by_substance.get(drug, [])
        if not candidates:
            unmatched.append(drug)
            continue
        for candidate in candidates:
            matched_rows.append({
                "drug": drug,
                "why_wanted": reason,
                "ema_medicine_name": str(candidate.get("Name of medicine", "")),
                "ema_product_number": str(candidate.get("EMA product number", "")),
                "medicine_status": str(candidate.get("Medicine status", "")),
                "marketing_authorisation_date": str(candidate.get("Marketing authorisation date", "") or "")[:10],
                "atc_code": str(candidate.get("ATC code (human)", "")),
                "medicine_url": str(candidate.get("Medicine URL", "")),
            })

    frame = pd.DataFrame(matched_rows)
    if len(frame):
        frame = frame.sort_values(["why_wanted", "drug"])
        frame.to_csv(arguments.output, index=False)

    summary = {"drugs_wanted": len(drugs_wanted),
               "drugs_with_an_epar": int(frame.drug.nunique()) if len(frame) else 0,
               "epar_rows": len(frame),
               "drugs_with_no_centrally_authorised_product": len(unmatched),
               "unmatched_drugs": sorted(unmatched),
               "by_reason": (frame.drop_duplicates("drug").why_wanted.value_counts().to_dict() if len(frame) else {})}
    arguments.summary.write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
