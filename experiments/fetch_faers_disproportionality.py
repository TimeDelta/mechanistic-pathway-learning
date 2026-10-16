"""Fetch FAERS reaction counts for the study drugs from openFDA and write disproportionality evidence reports
(mechanistic_pathway_learning/evidence/faers_disproportionality.py).

Three kinds of query, all cached under data/raw/faers (resumable; --refresh re-queries):
  the reactions reported with each drug (count of patient.reaction.reactionmeddrapt.exact) and the number of reports
  naming the drug (openFDA's normalised generic name, falling back to the verbatim medicinal product);
  the number of reports carrying each reaction the crosswalk names, over all drugs (the background);
  the number of reports in the system.
openFDA allows 40 requests a minute and 1,000 a day without a key, and refuses a count limit above 100 without one;
--pause keeps under the first, and a reaction outside a drug's 100 most reported enters at the smallest listed count
as an upper bound. A connection the egress proxy drops is retried, since every request is idempotent.

Outputs (default --output-dir data/processed/faers): faers_reports.parquet, faers_summary.json, and the generated
docs/faers_disproportionality_evidence.md.

Usage:
  python experiments/fetch_faers_disproportionality.py --evidence-dir data/processed/evidence_full_v4
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
from mechanistic_pathway_learning.evidence.clinical_trial_adverse_events import perturbations_from_reports, rows_to_dataframe
from mechanistic_pathway_learning.evidence.faers_disproportionality import (
    DEFAULT_NEGATIVE_MINIMUM_DRUG_REPORTS,
    DEFAULT_SIGNAL_CHI_SQUARE,
    DEFAULT_SIGNAL_MINIMUM_REPORTS,
    DEFAULT_SIGNAL_PRR,
    SOURCE_NAME,
    build_reports,
    indicated_symptoms_by_drug,
    write_summary,
)
from mechanistic_pathway_learning.evidence.load_drug_label_events import read_preferred_term_to_target_symptom

API = "https://api.fda.gov/drug/event.json"
USER_AGENT = "perturbation-symptom-modules (research; contact via GitHub)"
COUNT_LIMIT_WITHOUT_KEY = 100


def query(parameters: dict, pause: float) -> dict:
    address = API + "?" + urllib.parse.urlencode(parameters)
    request = urllib.request.Request(address, headers={"User-Agent": USER_AGENT})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = json.loads(response.read().decode("utf-8"))
            time.sleep(pause)
            return payload
        except urllib.error.HTTPError as error:
            if error.code == 404:  # openFDA answers 404 for an empty result set
                time.sleep(pause)
                return {"meta": {"results": {"total": 0}}, "results": []}
            if error.code == 429 and attempt < 4:
                time.sleep(30 * (attempt + 1))
                continue
            if error.code >= 500 and attempt < 4:  # a transient server error (the rebuild of 16 October 2026 died on one HTTP 500)
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except (urllib.error.URLError, ConnectionResetError, TimeoutError):
            if attempt < 4:
                time.sleep(5 * (attempt + 1))
                continue
            raise
    raise RuntimeError("unreachable")


def drug_counts(label: str, pause: float) -> dict:
    """Reaction counts and the report total for a drug, by openFDA's generic name first, then the verbatim product name."""
    for field in ("patient.drug.openfda.generic_name.exact", "patient.drug.medicinalproduct.exact"):
        search = f'{field}:"{label.upper()}"'
        total = int(query({"search": search, "limit": "1"}, pause).get("meta", {}).get("results", {}).get("total", 0))
        if total > 0:
            counts = query({"search": search, "count": "patient.reaction.reactionmeddrapt.exact", "limit": str(COUNT_LIMIT_WITHOUT_KEY)}, pause)
            listed = {item["term"]: int(item["count"]) for item in counts.get("results", [])}
            return {"field": field, "search": search, "drug_total": total, "reaction_counts": listed,
                    "smallest_listed_count": min(listed.values()) if listed else 0, "listed_reactions": len(listed)}
    return {"field": None, "search": None, "drug_total": 0, "reaction_counts": {}, "smallest_listed_count": 0, "listed_reactions": 0}


def background(preferred_terms: list[str], pause: float) -> dict:
    all_reports = int(query({"limit": "1"}, pause).get("meta", {}).get("results", {}).get("total", 0))
    counts = {}
    for term in preferred_terms:
        counts[term.upper()] = int(query({"search": f'patient.reaction.reactionmeddrapt.exact:"{term.upper()}"', "limit": "1"}, pause)
                                   .get("meta", {}).get("results", {}).get("total", 0))
    return {"all_reports": all_reports, "reaction_totals": counts}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--crosswalk", type=Path, default=Path("docs/symptom_crosswalk.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/faers"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/faers"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/faers_disproportionality_evidence.md"))
    parser.add_argument("--signal-prr", type=float, default=DEFAULT_SIGNAL_PRR)
    parser.add_argument("--signal-chi-square", type=float, default=DEFAULT_SIGNAL_CHI_SQUARE)
    parser.add_argument("--signal-minimum-reports", type=int, default=DEFAULT_SIGNAL_MINIMUM_REPORTS)
    parser.add_argument("--negative-minimum-drug-reports", type=int, default=DEFAULT_NEGATIVE_MINIMUM_DRUG_REPORTS)
    parser.add_argument("--pause", type=float, default=1.6)
    parser.add_argument("--limit-drugs", type=int, default=None)
    parser.add_argument("--refresh", action="store_true")
    arguments = parser.parse_args()

    thresholds = {"signal_prr": arguments.signal_prr, "signal_chi_square": arguments.signal_chi_square,
                  "signal_minimum_reports": arguments.signal_minimum_reports, "negative_minimum_drug_reports": arguments.negative_minimum_drug_reports}
    perturbations = perturbations_from_reports(arguments.evidence_dir / "evidence_reports.parquet")
    indicated = indicated_symptoms_by_drug(arguments.evidence_dir / "evidence_reports.parquet")
    preferred_term_to_symptom = read_preferred_term_to_target_symptom(arguments.crosswalk)
    preferred_terms_upper = {term.upper() for term in preferred_term_to_symptom}
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    background_path = arguments.cache_dir / "_background.json"
    if background_path.exists() and not arguments.refresh:
        background_counts = json.loads(background_path.read_text(encoding="utf-8"))
    else:
        background_counts = background(sorted(preferred_term_to_symptom), arguments.pause)
        background_counts["requested_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        background_path.write_text(json.dumps(background_counts, indent=1), encoding="utf-8")

    rows, per_drug = [], {}
    labels = sorted(perturbations)
    if arguments.limit_drugs:
        labels = labels[: arguments.limit_drugs]
    for index, label in enumerate(labels, start=1):
        cache_path = arguments.cache_dir / f"{label.replace('/', '_')}.json"
        if cache_path.exists() and not arguments.refresh:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
        else:
            try:
                cached = drug_counts(label, arguments.pause)
            except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ConnectionResetError) as error:
                print(f"{label}: request failed ({error}); skipped, not cached", file=sys.stderr)
                per_drug[label] = {"request_failed": 1}
                continue
            cached["requested_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
            cache_path.write_text(json.dumps(cached, indent=1), encoding="utf-8")
        if not cached["drug_total"]:
            per_drug[label] = {"no_reports_naming_the_drug": 1}
            continue
        reaction_counts = {reaction: count for reaction, count in cached["reaction_counts"].items() if reaction in preferred_terms_upper}
        unlisted_bound = cached.get("smallest_listed_count") if cached.get("listed_reactions", 0) >= COUNT_LIMIT_WITHOUT_KEY else None
        drug_rows, counts = build_reports(label, perturbations[label], reaction_counts, cached["drug_total"], background_counts["reaction_totals"],
                                          background_counts["all_reports"], preferred_term_to_symptom, indicated.get(perturbations[label]["perturbation_id"], set()),
                                          cached.get("requested_at", ""), thresholds, unlisted_upper_bound=unlisted_bound)
        counts["drug_total"] = int(cached["drug_total"])
        counts["matched_by"] = cached.get("field") or ""
        rows.extend(drug_rows)
        per_drug[label] = counts
        if index % 20 == 0:
            print(f"{index}/{len(labels)} drugs; {len(rows)} reports so far", flush=True)

    table = rows_to_dataframe(rows)
    table.to_parquet(arguments.output_dir / "faers_reports.parquet", index=False)
    write_summary(arguments.output_dir / "faers_summary.json", per_drug, {"all_reports": background_counts["all_reports"]}, thresholds)
    signals = int((table.report_value == 1).sum()) if len(table) else 0
    absences = int((table.report_value == 0).sum()) if len(table) else 0
    drugs_with_reports = int(table.perturbation_id.nunique()) if len(table) else 0
    without = sum(1 for counts in per_drug.values() if "no_reports_naming_the_drug" in counts)
    confounded = sum(counts.get("indication_confounded", 0) for counts in per_drug.values())
    lines = [
        "# FAERS disproportionality as evidence (generated by experiments/fetch_faers_disproportionality.py — do not edit by hand)", "",
        provenance_stamp("experiments/fetch_faers_disproportionality.py", ["openFDA drug/event (data/raw/faers)", str(arguments.evidence_dir)]), "",
        f"Source `{SOURCE_NAME}`: {len(labels)} study drugs queried against {background_counts['all_reports']:,} reports in the system; {without} drugs had no report "
        f"naming them. **{signals:,} signals** (PRR >= {arguments.signal_prr:g}, chi-square >= {arguments.signal_chi_square:g}, at least "
        f"{arguments.signal_minimum_reports} reports) and **{absences:,} absences** (ROR upper bound below 1 with at least "
        f"{arguments.negative_minimum_drug_reports:,} reports naming the drug) over {drugs_with_reports} drugs; {confounded} reaction rows were set aside "
        "because the symptom is one the drug is indicated for, where disproportion measures the indication rather than the drug.", "",
        "Spontaneous reports carry no exposure denominator, so a signal is disproportion against every other drug and an absence is under-reporting "
        "against background, not a rate; one report names several drugs, so counts are correlated across co-prescribed drugs. The rows join the "
        "assembled table through `--extra-reports` with the number of reports naming the drug in `rubric_log_sample_size`; the tier and the "
        "reliability fit decide their weight.", "",
    ]
    if len(table):
        by_symptom = table.groupby("symptom").report_value.agg(["count", "sum"])
        lines += ["| symptom | reports | signals | absences |", "| --- | --- | --- | --- |"]
        for symptom, row in by_symptom.sort_values("count", ascending=False).iterrows():
            lines.append(f"| {symptom} | {int(row['count'])} | {int(row['sum'])} | {int(row['count'] - row['sum'])} |")
        lines.append("")
    arguments.markdown_output.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {len(table)} reports for {drugs_with_reports} drugs to {arguments.output_dir}; {arguments.markdown_output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
