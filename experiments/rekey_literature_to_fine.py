"""Re-key a literature report table from the coarse constructs to their namesake fine constructs (decision 9.58).

The literature sources (PubTator3, CTD) are mapped through MeSH descriptors, which the fine crosswalk keeps with each
coarse construct's namesake, so the fine-level literature table is the coarse one with its symptom column renamed:
no row is added, dropped or re-mapped. Used where the raw literature cache is not on disk and the release table is.

Usage:
  python experiments/rekey_literature_to_fine.py --input data/releases/v0.4/literature/literature_reports.parquet \\
      --output data/processed/literature_v5/literature_reports.parquet
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.symptom_constructs import COARSE_CROSSWALK_PATH, namesake_of_coarse, read_rows


def rekey(reports: pd.DataFrame, namesake: dict[str, str]) -> pd.DataFrame:
    unknown = sorted(set(reports.symptom) - set(namesake))
    if unknown:
        raise ValueError(f"symptoms with no coarse construct: {unknown}")
    renamed = reports.copy()
    renamed["coarse_symptom"] = renamed.symptom
    renamed["symptom"] = renamed.symptom.map(namesake)
    if "report_id" in renamed.columns:
        renamed["report_id"] = [report_id.replace(f"|{coarse}|", f"|{fine}|") if coarse != fine else report_id
                                for report_id, coarse, fine in zip(renamed.report_id, renamed.coarse_symptom, renamed.symptom)]
    return renamed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, default=COARSE_CROSSWALK_PATH)
    arguments = parser.parse_args()
    namesake = namesake_of_coarse(read_rows(arguments.crosswalk))
    reports = pd.read_parquet(arguments.input)
    renamed = rekey(reports, namesake)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    renamed.to_parquet(arguments.output, index=False)
    renamed_count = int((renamed.symptom != renamed.coarse_symptom).sum())
    summary = {"input": str(arguments.input), "rows": int(len(renamed)), "rows_renamed": renamed_count, "namesake_of_coarse": namesake}
    (arguments.output.parent / "rekey_summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {arguments.output}: {len(renamed)} rows, {renamed_count} renamed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
