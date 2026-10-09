"""Check that adding parkinsonism to docs/symptom_crosswalk.csv changed only what a new symptom has to change.

The rebuild (runs/full/parkinsonism_evidence.sh) runs the commands recorded in docs/drug_targets_any_type.md, so the
first thing to establish is that nothing but the new symptom moved. Two columns of the evidence table are fitted or
derived across symptoms and so cannot stay fixed when a symptom is added:

  - `reliability_posterior`, the hierarchical fit over reports, which sees the new symptom's drug reports;
  - `disease_cluster_id`, the connected component of the disease graph, where a new symptom's annotations can join two
    components that no earlier symptom linked.

Every other column of every row that is not a parkinsonism row has to match the reference exactly, and so does the
whole label selection outside the new symptom. The two fitted columns are then reported with their size, because a
large move in either is a reason to look rather than to proceed.

The reference is evidence_full_v3. The salt-form and any-target-type change of 8 October
(docs/drug_targets_any_type.md: "The v2 tables rebuild byte-identically from the code before this change") is in the
code, so a rebuild today lands on v3 whatever the crosswalk says, and v2 is not reproducible from the current code.
The user's decision of 9 October 2026 ("just get rid of the v2 test at this point since it was never used anyways.
v3 is the usable one for now") is why v2 is neither rebuilt nor reported here.

Writes docs/parkinsonism_symptom.md. Reads no lockbox outcome. Idempotent.

Usage:
  python experiments/compare_parkinsonism_evidence.py
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

NEW_SYMPTOM = "parkinsonism"
# fitted or derived across symptoms, so a new symptom moves them by construction; reported rather than asserted
COLUMNS_A_NEW_SYMPTOM_MAY_MOVE = ("reliability_posterior", "disease_cluster_id")
ROW_ORDER = ["perturbation_id", "symptom", "source", "relation"]
OUTPUT_DOCUMENT = Path("docs/parkinsonism_symptom.md")


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def rows_outside_the_new_symptom(records: pd.DataFrame) -> pd.DataFrame:
    """The rows of a table that are not about the new symptom, in a fixed order, for a row-by-row comparison."""
    rows = records[records.symptom != NEW_SYMPTOM].copy()
    order = [column for column in ROW_ORDER if column in rows.columns] or list(rows.columns)
    return rows.sort_values(order, kind="mergesort").reset_index(drop=True)


def compare_evidence(reference: pd.DataFrame, new: pd.DataFrame) -> dict:
    """Column-by-column comparison of the rows outside the new symptom."""
    left, right = rows_outside_the_new_symptom(reference), rows_outside_the_new_symptom(new)
    if len(left) != len(right):
        return {"identical": False, "row_counts": (len(left), len(right)), "differing_columns": [],
                "posterior_shift": None, "cluster_changes": pd.DataFrame()}
    differing = [column for column in left.columns
                 if column not in COLUMNS_A_NEW_SYMPTOM_MAY_MOVE
                 and not left[column].astype(str).equals(right[column].astype(str))]
    posterior_shift = float((left.reliability_posterior - right.reliability_posterior).abs().max())
    changed = left.disease_cluster_id.astype(str) != right.disease_cluster_id.astype(str)
    cluster_changes = pd.DataFrame({"perturbation_id": left.perturbation_id[changed],
                                    "symptom": left.symptom[changed],
                                    "was": left.disease_cluster_id[changed],
                                    "now": right.disease_cluster_id[changed]})
    return {"identical": not differing, "row_counts": (len(left), len(right)), "differing_columns": differing,
            "posterior_shift": posterior_shift, "cluster_changes": cluster_changes}


def table_row(label: str, *cells: str) -> str:
    return f"| {label} | " + " | ".join(cells) + " |"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference-evidence", type=Path, default=Path("data/processed/evidence_full_v3"))
    parser.add_argument("--new-evidence", type=Path, default=Path("data/processed/evidence_full_v3_parkinsonism"))
    parser.add_argument("--reference-selection", type=Path,
                        default=Path("data/processed/label_selection/better_v2_full_v3.parquet"))
    parser.add_argument("--new-selection", type=Path,
                        default=Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"))
    arguments = parser.parse_args()

    reference = pd.read_parquet(arguments.reference_evidence / "evidence_records.parquet")
    new = pd.read_parquet(arguments.new_evidence / "evidence_records.parquet")
    evidence = compare_evidence(reference, new)

    reference_selection = pd.read_parquet(arguments.reference_selection)
    new_selection = pd.read_parquet(arguments.new_selection)
    selection_identical = rows_outside_the_new_symptom(reference_selection).equals(
        rows_outside_the_new_symptom(new_selection))
    sound = evidence["identical"] and selection_identical

    new_rows = new_selection[new_selection.symptom == NEW_SYMPTOM]
    kept = new_rows[new_rows.keep]
    by_type = kept.perturbation_type.value_counts().to_dict()
    kept_positives_by_type = ", ".join(f"{number} {name}" for name, number in sorted(by_type.items())) or "none"
    symptoms_above_the_macro_floor = int((reference_selection[reference_selection.keep]
                                          .symptom.value_counts() >= 5).sum())
    cluster_changes = evidence["cluster_changes"]
    selections = (reference_selection, new_selection)
    evidence_directories = (arguments.reference_evidence, arguments.new_evidence)
    selection_paths = (arguments.reference_selection, arguments.new_selection)

    lines = [
        "# Parkinsonism as the 24th symptom",
        "",
        "Generated by `experiments/compare_parkinsonism_evidence.py`. Do not edit by hand.",
        "",
        f"**Outside parkinsonism the rebuild {'reproduces' if sound else 'DOES NOT reproduce'} `evidence_full_v3` and "
        f"`better_v2_full_v3`.** Every evidence column a new symptom cannot move is "
        f"{'identical' if evidence['identical'] else 'NOT identical (' + ', '.join(evidence['differing_columns']) + ')'}"
        f" and the label selection outside the symptom is "
        f"{'identical' if selection_identical else 'NOT identical'}."
        + ("" if sound else " Something other than the symptom moved and nothing below should be read."),
        "",
        "The reference is v3, the table the current code produces. The salt-form and any-target-type change of",
        "8 October (`docs/drug_targets_any_type.md`) is in the code, so a rebuild today lands on v3 whatever the",
        "crosswalk says; v2 is not reproducible from the current code and is not reported here.",
        "",
        "| | v3, the current code | v3 with parkinsonism |",
        "| --- | --- | --- |",
    ]
    lines += [table_row(label, *cells) for label, cells in (
        ("evidence rows", (f"{len(reference):,}", f"{len(new):,}")),
        ("distinct symptoms", tuple(str(table.symptom.nunique()) for table in (reference, new))),
        ("selection rows", tuple(f"{len(table):,}" for table in selections)),
        ("kept positives", tuple(f"{int(table.keep.sum()):,}" for table in selections)),
        ("disease clusters", tuple(f"{table.disease_cluster_id.nunique():,}" for table in (reference, new))),
        ("evidence_records sha256 (first 16)",
         tuple(file_digest(directory / "evidence_records.parquet") for directory in evidence_directories)),
        ("label selection sha256 (first 16)", tuple(file_digest(path) for path in selection_paths)),
    )]
    lines += [
        "",
        "## What parkinsonism brings",
        "",
        f"- evidence rows: {int((new.symptom == NEW_SYMPTOM).sum()):,}",
        f"- selection rows: {len(new_rows):,}",
        f"- kept positives: {int(kept.keep.sum()):,} ({kept_positives_by_type})",
        f"- perturbations with a kept positive: {kept.perturbation_id.nunique():,}",
        "",
        f"The macro average takes a symptom with at least 5 kept positives. {symptoms_above_the_macro_floor} symptoms",
        "clear that floor in v3, so the symptom enters the macro average and moves its denominator by one.",
        "",
        "## What the new symptom moves outside itself",
        "",
        "Two columns are fitted or derived across symptoms, so a new symptom moves them by construction. Both are",
        "measured here rather than asserted away.",
        "",
        f"- `reliability_posterior` shifts by at most **{evidence['posterior_shift']:.2g}** over the "
        + f"{evidence['row_counts'][0]:,} rows outside parkinsonism. The posterior is fitted over reports, and the",
        "  symptom adds drug label reports, so every row's shrinkage moves a little. The label selection reads the",
        "  frequency and the grade rather than the posterior, which is why it is unchanged outside the symptom.",
        f"- `disease_cluster_id` changes on **{len(cluster_changes)}** of those rows. A cluster is a connected",
        "  component of the disease graph, so an annotation the new symptom brings can join two components that no",
        "  earlier symptom linked. The cluster is a leakage group, so each change moves those perturbations' fold.",
    ]
    if len(cluster_changes):
        lines += ["", "| perturbation | symptom | cluster was | cluster now |", "| --- | --- | --- | --- |"]
        for row in cluster_changes.itertuples(index=False):
            lines.append(f"| {row.perturbation_id} | {row.symptom} | `{row.was}` | `{row.now}` |")
    lines += [
        "",
        "## What has to follow before any confirmatory run",
        "",
        "`configs/lockbox_v2.json` pins the `evidence_records_sha256` and `label_selection_sha256` of the v2 tables,",
        "which the current code no longer produces, so it describes neither column above. A confirmatory run on this",
        "table needs a lockbox drawn from the same rule on it, and the development runs repeated. None of that is done",
        "here, and nothing in this script reads a lockbox outcome.",
    ]
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")
    print(f"evidence columns a new symptom cannot move, identical outside parkinsonism: {evidence['identical']}"
          + (f" (differing: {', '.join(evidence['differing_columns'])})" if evidence["differing_columns"] else ""))
    print(f"label selection identical outside parkinsonism: {selection_identical}")
    print(f"reliability_posterior maximum shift: {evidence['posterior_shift']:.3g}")
    print(f"rows whose disease cluster changed: {len(cluster_changes)}")
    return 0 if sound else 1


if __name__ == "__main__":
    raise SystemExit(main())
