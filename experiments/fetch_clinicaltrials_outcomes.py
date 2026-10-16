"""Fetch the outcome measures of the study drugs' trials and read placebo-controlled comparisons on named scales.

Decision 9.60 (the owner, 16 October 2026). The adverse-event cache (data/raw/clinicaltrials) holds no outcome module,
so this script queries ClinicalTrials.gov API v2 again for the same drugs with the OutcomeMeasuresModule and caches
each drug's studies under data/raw/clinicaltrials_outcomes/<drug>.json; a rerun reads the cache. The reading is in
mechanistic_pathway_learning/evidence/clinical_trial_outcomes.py and the scale-to-construct table is
configs/trial_outcome_scale_crosswalk.csv (fine constructs; --symptom-level coarse maps them to the coarse crosswalk's).

Outputs: <output-dir>/clinical_trial_outcome_reports.parquet, clinical_trial_outcome_summary.json and the generated
page docs/clinical_trial_outcome_evidence.md.

Usage:
  python experiments/fetch_clinicaltrials_outcomes.py --evidence-dir data/processed/evidence_full_v5
  python experiments/fetch_clinicaltrials_outcomes.py --evidence-dir data/processed/evidence_full_v4 --symptom-level coarse --output-dir data/processed/clinical_trial_outcomes_v4
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

from mechanistic_pathway_learning.documentation import provenance_stamp
from mechanistic_pathway_learning.evidence.clinical_trial_adverse_events import perturbations_from_reports
from mechanistic_pathway_learning.evidence.clinical_trial_outcomes import (
    DEFAULT_SCALE_CROSSWALK,
    NO_EFFECT_MARGIN,
    SOURCE_NAME,
    build_reports,
    read_scale_crosswalk,
    rows_to_dataframe,
)
from mechanistic_pathway_learning.evidence.symptom_constructs import FINE_CROSSWALK_PATH, coarse_of_fine, read_rows

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fetch_clinicaltrials_adverse_events import API, USER_AGENT, request_json  # noqa: E402

FIELDS = "NCTId,BriefTitle,OverallStatus,StudyType,DesignInfo,ResultsFirstPostDate,ArmsInterventionsModule,OutcomeMeasuresModule"


def fetch_studies(drug_name: str, pause: float, page_size: int = 50) -> tuple[list[dict], list[str]]:
    import time
    import urllib.parse

    studies, addresses, page_token = [], [], None
    while True:
        parameters = {"query.intr": drug_name, "filter.advanced": "AREA[HasResults]true", "fields": FIELDS, "pageSize": str(page_size), "format": "json"}
        if page_token:
            parameters["pageToken"] = page_token
        address = API + "?" + urllib.parse.urlencode(parameters)
        addresses.append(address)
        payload = request_json(address)
        studies.extend(payload.get("studies", []))
        page_token = payload.get("nextPageToken")
        time.sleep(pause)
        if not page_token:
            return studies, addresses


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v5"), help="the assembled table whose drugs are queried")
    parser.add_argument("--scale-crosswalk", type=Path, default=DEFAULT_SCALE_CROSSWALK)
    parser.add_argument("--symptom-level", choices=("fine", "coarse"), default="fine", help="the constructs the reports name; coarse maps each fine construct to its coarse one")
    parser.add_argument("--fine-crosswalk", type=Path, default=FINE_CROSSWALK_PATH)
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/clinicaltrials_outcomes"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/clinical_trial_outcomes"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/clinical_trial_outcome_evidence.md"))
    parser.add_argument("--no-effect-margin", type=float, default=NO_EFFECT_MARGIN, help="standardised difference inside which a construct is read as not moved; a choice (9.48's thresholds)")
    parser.add_argument("--pause", type=float, default=0.5)
    parser.add_argument("--limit-drugs", type=int, default=None)
    parser.add_argument("--refresh", action="store_true")
    arguments = parser.parse_args()

    perturbations = perturbations_from_reports(arguments.evidence_dir / "evidence_reports.parquet")
    rules = read_scale_crosswalk(arguments.scale_crosswalk)
    import mechanistic_pathway_learning.evidence.clinical_trial_outcomes as outcomes_module
    outcomes_module.NO_EFFECT_MARGIN = arguments.no_effect_margin
    to_coarse = coarse_of_fine(read_rows(arguments.fine_crosswalk)) if arguments.symptom_level == "coarse" else None
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    rows, per_drug = [], {}
    labels = sorted(perturbations)[: arguments.limit_drugs] if arguments.limit_drugs else sorted(perturbations)
    for index, label in enumerate(labels, start=1):
        cache_path = arguments.cache_dir / f"{label.replace('/', '_')}.json"
        if cache_path.exists() and not arguments.refresh:
            studies = json.loads(cache_path.read_text(encoding="utf-8"))["studies"]
        else:
            try:
                studies, addresses = fetch_studies(label, arguments.pause)
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ConnectionResetError) as error:
                print(f"{label}: request failed ({error}); skipped, not cached", file=sys.stderr)
                per_drug[label] = {"request_failed": 1}
                continue
            cache_path.write_text(json.dumps({"drug": label, "requested_at": datetime.now(timezone.utc).isoformat(timespec="seconds"), "user_agent": USER_AGENT,
                                              "addresses": addresses, "studies": studies}, indent=1), encoding="utf-8")
        drug_rows, counts = build_reports(studies, label, perturbations[label], rules, other_study_drugs=set(perturbations) - {label})
        if to_coarse is not None:
            for row in drug_rows:
                row["symptom"] = to_coarse.get(row["symptom"], row["symptom"])
        rows.extend(drug_rows)
        per_drug[label] = counts
        if index % 20 == 0:
            print(f"{index}/{len(labels)} drugs; {len(rows)} comparisons so far", flush=True)

    comparisons = rows_to_dataframe(rows)
    comparisons.to_parquet(arguments.output_dir / "clinical_trial_outcome_comparisons.parquet", index=False)  # every comparison, the undecided included
    table = comparisons[comparisons.report_value >= 0].reset_index(drop=True)  # the report table: a reader of report_value must never see -1
    table.to_parquet(arguments.output_dir / "clinical_trial_outcome_reports.parquet", index=False)
    totals: dict[str, int] = {}
    for counts in per_drug.values():
        for key, value in counts.items():
            if isinstance(value, int):
                totals[key] = totals.get(key, 0) + value
    summary = {"source": SOURCE_NAME, "symptom_level": arguments.symptom_level, "no_effect_margin": arguments.no_effect_margin, "drugs_queried": len(labels),
               "totals": totals, "per_drug": per_drug, "scale_crosswalk": str(arguments.scale_crosswalk)}
    (arguments.output_dir / "clinical_trial_outcome_summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")

    by_construct = table.groupby(["symptom", "relation", "report_value"]).size() if len(table) else None
    lines = [
        "# Trial outcome measures as evidence (generated by experiments/fetch_clinicaltrials_outcomes.py — do not edit by hand)", "",
        provenance_stamp("experiments/fetch_clinicaltrials_outcomes.py", ["ClinicalTrials.gov API v2 (data/raw/clinicaltrials_outcomes)", str(arguments.evidence_dir), str(arguments.scale_crosswalk)]), "",
        f"Source `{SOURCE_NAME}` at the {arguments.symptom_level} level: {len(labels)} study drugs queried; {totals.get('studies', 0):,} studies with posted results; "
        f"{totals.get('outcome_measures', 0):,} outcome measures, {totals.get('crosswalked', 0):,} on a crosswalked scale, {totals.get('placebo_controlled', 0):,} of those "
        f"with a drug arm and a placebo arm that could be told apart, {totals.get('compared', 0):,} read as a mean comparison: **{totals.get('relief', 0):,} relief**, "
        f"**{totals.get('worsening', 0):,} worsening**, **{totals.get('no_effect', 0):,} no effect** (standardised interval within {arguments.no_effect_margin} of 0) and "
        f"{totals.get('undecided', 0):,} undecided (in clinical_trial_outcome_comparisons.parquet with report_value -1, never in the report table).", "",
        "A relief is relation `relieves` with report_value 1, a worsening is `induces` 1, a no-effect is `relieves` 0. Class human_association, grade C: a soft prior, "
        "never a label. A total score that fuses several constructs (PANSS total, an ADHD total) matches no rule and is not read; the subscales are.", "",
    ]
    if by_construct is not None:
        lines += ["| construct | relation | value | comparisons |", "| --- | --- | --- | --- |"]
        for (symptom, relation, value), count in by_construct.items():
            lines.append(f"| {symptom} | {relation} | {value} | {count} |")
    lines.append("")
    arguments.markdown_output.parent.mkdir(parents=True, exist_ok=True)
    arguments.markdown_output.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {arguments.output_dir / 'clinical_trial_outcome_reports.parquet'} ({len(table)} rows of {len(comparisons)} comparisons) and {arguments.markdown_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
