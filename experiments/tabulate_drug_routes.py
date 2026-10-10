"""Tabulate the routes of administration of each study drug's FDA labels, against whether a plasma carrier is known.

Why. A carrier's sequestration edge (docs/drug_entry_nodes.md) says that more carrier means less free drug and so less
effect. Under the well-stirred model of hepatic clearance that holds for the unbound exposure of a drug only when the
dose bypasses the liver's first pass and the liver extracts most of what reaches it: Benet 2018 (J Clin Pharmacol
58:979, Table 1) gives "Orally dosed; hepatic elimination: No", "IV dosed; low ER: No", "IV dosed; high ER: Yes". So
whether the edge can mean what it says depends on the route and the extraction ratio of each drug, and this script
supplies the route half from the labels the study already reads. The extraction ratio is not in any source the
repository holds and is not tabulated here.

What it reads. openFDA's count of the `openfda.route` field over every label of a generic name, one request per drug,
cached under --cache-dir (raw downloaded data, in its own directory; no code is run from it). A count is over labels,
so a drug sold as a tablet by forty manufacturers and as an injection by two shows both routes with those counts, and a
combination product counts for each of its ingredients.

What it writes. --markdown-output, with the route classes per drug and their totals by carrier, and --table-output, a
CSV with one row per drug. The classes: `intravenous` when a label names the intravenous route; `other route past the
first pass` when a label names another route that reaches the circulation without passing the liver first
(ROUTES_PAST_THE_FIRST_PASS) and none names the intravenous route; `oral or local only` otherwise. A local route
(topical, ophthalmic, dental) is not counted as reaching the circulation.

Usage:
  python experiments/tabulate_drug_routes.py
"""
from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import pandas as pd

LABEL_ENDPOINT = "https://api.fda.gov/drug/label.json"
SECONDS_BETWEEN_REQUESTS = 0.34  # openFDA allows 240 unauthenticated requests a minute
RETRY_WAITS = (2, 4, 8, 16)
INTRAVENOUS_ROUTE = "INTRAVENOUS"
# routes that deliver the dose to the circulation without the liver's first pass, the intravenous route apart; the
# spinal and nerve-block routes are here because a local anaesthetic given by them is absorbed into the circulation
ROUTES_PAST_THE_FIRST_PASS = (
    "INTRAMUSCULAR", "SUBCUTANEOUS", "TRANSDERMAL", "SUBLINGUAL", "BUCCAL", "NASAL", "RESPIRATORY (INHALATION)",
    "EPIDURAL", "INTRATHECAL", "PERINEURAL", "INFILTRATION", "INTRACAUDAL", "INTRA-ARTERIAL", "INTRAOSSEOUS", "PARENTERAL",
    "INTRAVASCULAR", "INTRACAVERNOUS", "INTRASPINAL", "SUBARACHNOID",
)
ROUTE_CLASSES = ("intravenous", "other route past the first pass", "oral or local only", "no label route")


def fetch_route_counts(drug_label: str) -> tuple[str, list[dict]]:
    """openFDA's count of each route over the labels of one generic name; an empty list when it has no such label."""
    query = urllib.parse.quote(f'openfda.generic_name:"{drug_label}"')
    url = f"{LABEL_ENDPOINT}?search={query}&count=openfda.route.exact"
    for attempt, wait in enumerate((0,) + RETRY_WAITS):
        if wait:
            time.sleep(wait)
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return url, json.loads(response.read()).get("results", [])
        except urllib.error.HTTPError as error:
            if error.code == 404:  # openFDA answers a search with no match with 404
                return url, []
            if error.code not in (429, 500, 502, 503, 504) or attempt == len(RETRY_WAITS):
                raise
        except (urllib.error.URLError, TimeoutError):
            if attempt == len(RETRY_WAITS):
                raise
    return url, []


def route_class(routes: list[str]) -> str:
    if not routes:
        return "no label route"
    if INTRAVENOUS_ROUTE in routes:
        return "intravenous"
    if any(route in ROUTES_PAST_THE_FIRST_PASS for route in routes):
        return "other route past the first pass"
    return "oral or local only"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3_parkinsonism"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"))
    parser.add_argument("--carrier-table", type=Path, default=Path("configs/drug_plasma_carriers.csv"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/raw/fda_labels/routes"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/drug_administration_routes.md"))
    parser.add_argument("--table-output", type=Path, default=Path("configs/drug_administration_routes.csv"))
    parser.add_argument("--refetch", action="store_true")
    arguments = parser.parse_args()

    records = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet", columns=["perturbation_id", "perturbation_type", "perturbation_label"])
    drugs = (records[records.perturbation_type == "drug"][["perturbation_id", "perturbation_label"]]
             .drop_duplicates().sort_values("perturbation_label").reset_index(drop=True))
    arguments.cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for drug in drugs.itertuples(index=False):
        cache_path = arguments.cache_dir / f"{drug.perturbation_id}.json"
        if cache_path.exists() and not arguments.refetch:
            cached = json.loads(cache_path.read_text())
        else:
            url, counts = fetch_route_counts(drug.perturbation_label)
            cached = {"perturbation_id": drug.perturbation_id, "drug": drug.perturbation_label, "query": url,
                      "retrieved": date.today().isoformat(), "route_counts": counts}
            cache_path.write_text(json.dumps(cached, indent=1))
            time.sleep(SECONDS_BETWEEN_REQUESTS)
        route_counts = {entry["term"]: int(entry["count"]) for entry in cached["route_counts"]}
        ordered_routes = sorted(route_counts, key=lambda route: (-route_counts[route], route))
        rows.append({"perturbation_id": drug.perturbation_id, "drug": drug.perturbation_label, "route_class": route_class(ordered_routes),
                     "routes": "; ".join(f"{route} ({route_counts[route]})" for route in ordered_routes), "retrieved": cached["retrieved"]})
    table = pd.DataFrame(rows)

    carriers = pd.read_csv(arguments.carrier_table)
    carrier_of_drug = dict(zip(carriers.perturbation_id, carriers.carrier))
    table["carrier"] = table.perturbation_id.map(carrier_of_drug).fillna("")
    selection = pd.read_parquet(arguments.label_selection)
    kept_positives = selection[selection.keep.astype(bool)].groupby("perturbation_id").size()  # the positive pairs the selection keeps
    table["kept_positives"] = table.perturbation_id.map(kept_positives).fillna(0).astype(int)
    arguments.table_output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(arguments.table_output, index=False)

    def totals(subset: pd.DataFrame) -> str:
        return f"{len(subset)} | {int(subset.kept_positives.sum())}"

    lines = [
        "# Routes of administration of the study drugs, by whether a plasma carrier is known",
        "",
        "Generated by `experiments/tabulate_drug_routes.py`. Do not edit by hand.",
        "",
        f"Source: openFDA's count of the `openfda.route` field over every label of each drug's generic name, retrieved {table.retrieved.min()}"
        + ("" if table.retrieved.min() == table.retrieved.max() else f" to {table.retrieved.max()}") + ". A count is the number of labels that name the route, "
        "over every manufacturer and every combination product of the ingredient.",
        "",
        "Why this table exists. A carrier's sequestration edge states that more carrier lowers the drug's effect. For a drug cleared by the liver, the "
        "unbound exposure depends on the unbound fraction only when the dose bypasses the first pass and the liver extracts most of the drug "
        "(Benet 2018, Table 1; Hochman et al. 2015, Table 2; quotes in `docs/drug_entry_nodes.md`). This table gives the route half of that "
        "condition. **The extraction ratio is not tabulated: no source in this repository holds a clearance per drug.** The class of a drug is the "
        "most permissive route any of its labels names, so `intravenous` means an intravenous product exists, and the adverse-effect reports of this "
        "study do not say which product a patient took.",
        "",
        "## Totals",
        "",
        "| route class | drugs | kept positives | with a known carrier: drugs | kept positives |",
        "| --- | --- | --- | --- | --- |",
    ]
    with_carrier = table[table.carrier != ""]
    for name in ROUTE_CLASSES:
        of_class = table[table.route_class == name]
        if len(of_class):
            lines.append(f"| {name} | {totals(of_class)} | {totals(with_carrier[with_carrier.route_class == name])} |")
    lines += [f"| all | {totals(table)} | {totals(with_carrier)} |", "",
              "## The drugs with a known carrier", "", "| drug | carrier | route class | kept positives | routes (labels) |", "| --- | --- | --- | --- | --- |"]
    for row in with_carrier.sort_values(["route_class", "drug"]).itertuples(index=False):
        lines.append(f"| {row.drug} | {row.carrier} | {row.route_class} | {row.kept_positives} | {row.routes} |")
    lines += ["", "## The drugs without a known carrier", "", "| drug | route class | kept positives | routes (labels) |", "| --- | --- | --- | --- |"]
    for row in table[table.carrier == ""].sort_values(["route_class", "drug"]).itertuples(index=False):
        lines.append(f"| {row.drug} | {row.route_class} | {row.kept_positives} | {row.routes} |")
    arguments.markdown_output.write_text("\n".join(lines) + "\n")
    print(table.route_class.value_counts().to_string())
    print(f"wrote {arguments.markdown_output} and {arguments.table_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
