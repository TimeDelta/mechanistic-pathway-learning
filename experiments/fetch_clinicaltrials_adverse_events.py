"""Fetch the results sections of placebo-controlled trials for the study drugs from ClinicalTrials.gov and write their
adverse-event arm comparisons as evidence reports (mechanistic_pathway_learning/evidence/clinical_trial_adverse_events.py).

Cached and resumable: each drug's studies are saved as data/raw/clinicaltrials/<label>.json with the request address
and time, and a drug with a cached file is not queried again (--refresh re-queries). Rate: one request every
--pause seconds (default 0.5) against the public API v2, which needs no key. A connection the egress proxy drops is
retried, since every request is idempotent.

Outputs (default --output-dir data/processed/clinical_trials):
  clinical_trial_reports.parquet   report rows in the evidence_reports schema, for --extra-reports of the assembler
  clinical_trial_comparisons.parquet   every arm comparison on a target symptom, whatever its verdict, with its own
                                   statistic (decision 9.48, option 3)
  clinical_trial_summary.json      counts per drug and in total
  docs/clinical_trial_arm_evidence.md (--markdown-output)   the generated summary, stamped

Usage:
  python experiments/fetch_clinicaltrials_adverse_events.py --evidence-dir data/processed/evidence_full_v4
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from mechanistic_pathway_learning.documentation import provenance_stamp
from mechanistic_pathway_learning.evidence.clinical_trial_adverse_events import (
    COMPARISON_COLUMNS,
    DEFAULT_MINIMUM_EXPOSURE_FOR_ABSENCE,
    SOURCE_NAME,
    build_comparisons,
    build_reports,
    comparisons_to_dataframe,
    perturbations_from_reports,
    read_lowest_level_to_preferred_term,
    rows_to_dataframe,
    write_summary,
)
from mechanistic_pathway_learning.evidence.load_drug_label_events import read_preferred_term_to_target_symptom

API = "https://clinicaltrials.gov/api/v2/studies"
FIELDS = "NCTId,BriefTitle,OverallStatus,StudyType,DesignInfo,ResultsFirstPostDate,ArmsInterventionsModule,AdverseEventsModule"
USER_AGENT = "perturbation-symptom-modules (research; contact via GitHub)"


def request_json(address: str) -> dict:
    request = urllib.request.Request(address, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code in (429, 500, 502, 503, 504) and attempt < 4:  # rate limited, or a transient server error
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, ConnectionResetError, TimeoutError):
            if attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            raise
    raise RuntimeError("unreachable")


def fetch_studies(drug_name: str, pause: float, page_size: int = 100) -> tuple[list[dict], list[str]]:
    """Every study with posted results whose interventions name the drug, through the API's paging."""
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
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"),
                        help="the assembled table whose drugs are queried; their identities and target nodes are reused")
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--sider-dir", type=Path, default=Path("data/raw/sider_4.1"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/clinicaltrials"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/clinical_trials"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/clinical_trial_arm_evidence.md"))
    parser.add_argument("--minimum-exposure-for-absence", type=int, default=DEFAULT_MINIMUM_EXPOSURE_FOR_ABSENCE)
    parser.add_argument("--pause", type=float, default=0.5)
    parser.add_argument("--limit-drugs", type=int, default=None, help="query only the first N drugs (for a dry run)")
    parser.add_argument("--refresh", action="store_true", help="re-query drugs that have a cached file")
    arguments = parser.parse_args()

    perturbations = perturbations_from_reports(arguments.evidence_dir / "evidence_reports.parquet")
    preferred_term_to_symptom = read_preferred_term_to_target_symptom(arguments.crosswalk)
    lowest_to_preferred = read_lowest_level_to_preferred_term(arguments.sider_dir)
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    rows, comparisons, per_drug = [], [], {}
    labels = sorted(perturbations)
    if arguments.limit_drugs:
        labels = labels[: arguments.limit_drugs]
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
            cache_path.write_text(json.dumps({"drug": label, "requested_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                                              "addresses": addresses, "studies": studies}, indent=1), encoding="utf-8")
        drug_rows, counts = build_reports(studies, label, perturbations[label], preferred_term_to_symptom, lowest_to_preferred,
                                          arguments.minimum_exposure_for_absence, other_study_drugs=set(perturbations) - {label})
        comparisons.extend(build_comparisons(studies, label, perturbations[label], preferred_term_to_symptom, lowest_to_preferred,
                                             arguments.minimum_exposure_for_absence, other_study_drugs=set(perturbations) - {label}))
        rows.extend(drug_rows)
        per_drug[label] = counts
        if index % 20 == 0:
            print(f"{index}/{len(labels)} drugs; {len(rows)} reports so far", flush=True)

    table = rows_to_dataframe(rows)
    table.to_parquet(arguments.output_dir / "clinical_trial_reports.parquet", index=False)
    # every arm comparison, whatever its verdict, with its own statistic (decision 9.48, option 3)
    comparisons_to_dataframe(comparisons, COMPARISON_COLUMNS).to_parquet(arguments.output_dir / "clinical_trial_comparisons.parquet", index=False)
    write_summary(arguments.output_dir / "clinical_trial_summary.json", per_drug, arguments.minimum_exposure_for_absence)

    totals: dict[str, int] = {}
    for counts in per_drug.values():
        for key, value in counts.items():
            if isinstance(value, int):
                totals[key] = totals.get(key, 0) + value
    positives = int((table.report_value == 1).sum()) if len(table) else 0
    negatives = int((table.report_value == 0).sum()) if len(table) else 0
    drugs_with_reports = int(table.perturbation_id.nunique()) if len(table) else 0
    lines = [
        "# Placebo-controlled trial arms as evidence (generated by experiments/fetch_clinicaltrials_adverse_events.py — do not edit by hand)", "",
        provenance_stamp("experiments/fetch_clinicaltrials_adverse_events.py", ["ClinicalTrials.gov API v2 (data/raw/clinicaltrials)", str(arguments.evidence_dir)]), "",
        f"Source `{SOURCE_NAME}`: {len(labels)} study drugs queried; {totals.get('studies', 0):,} studies with posted results named a drug among their interventions; "
        f"{totals.get('studies_with_usable_arms', 0):,} of them had a drug arm and a placebo arm that could be told apart, after "
        f"{totals.get('arms_excluded_for_another_study_drug', 0):,} arms naming a second study drug were left out of both sides (an arm's events reach one drug's "
        f"record only, so nothing crosses the held-out split through a combination); {totals.get('terms', 0):,} listed terms, "
        f"{totals.get('terms_on_a_symptom', 0):,} on a target symptom; **{positives:,} positive reports** (drug-arm rate above placebo) and "
        f"**{negatives:,} negative reports** (not above placebo, drug arm at least {arguments.minimum_exposure_for_absence} at risk) over {drugs_with_reports} drugs; "
        f"{totals.get('non_excess_below_exposure_floor', 0):,} non-excess terms fell below the exposure floor and produced no report.", "",
        "A term a trial does not list is below its frequency threshold in every arm, not absent, so it produces nothing. Arms are pooled across "
        "doses. The rows join the assembled table through `--extra-reports` and carry `rubric_placebo_controlled` 1 and the drug arm's size in "
        "`rubric_log_sample_size`; the reliability fit and the tier weights decide what a report is worth.", "",
    ]
    if len(table):
        by_symptom = table.groupby("symptom").report_value.agg(["count", "sum"])
        lines += ["| symptom | reports | positive | negative |", "| --- | --- | --- | --- |"]
        for symptom, row in by_symptom.sort_values("count", ascending=False).iterrows():
            lines.append(f"| {symptom} | {int(row['count'])} | {int(row['sum'])} | {int(row['count'] - row['sum'])} |")
        lines.append("")
    arguments.markdown_output.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {len(table)} reports for {drugs_with_reports} drugs to {arguments.output_dir}; {arguments.markdown_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
