"""How much of this study's scored signal rests on the drugs orosomucoid carries, and whether the binder genes can be
training examples in their own right (the user, 9 October 2026: "If orosomucoid carries 11 of 27 drugs by itself, we
NEED a physiologically and computationally viable way to deal with it ... If it just means dropping some portion of
its pairs, check which pairs and how important they are for psych symptoms. Knockouts for the orosomucoid binder
would also be fine as training examples.").

Two questions, answered from the repository's own inputs and nothing else:

  1. the stake. The pinned fraction-unbound database (docs/fraction_unbound_coverage.md) names orosomucoid as the
     major binding protein for 11 of this study's drug perturbations. This script counts the kept positive pairs
     those 11 hold, per drug and per symptom, against the whole label selection, so that "drop the orosomucoid
     drugs" has a price rather than a worry;
  2. the knockouts. A gene perturbation's labels come from HPO disease annotations, so a binder gene can be a
     training example only if HPO annotates it. This script reports the annotation count for every plasma binder,
     which is what decides whether ORM1 and ORM2 can enter at all.

Writes docs/orosomucoid_label_stake.md. Reads no outcome of any lockbox run. Idempotent.

Usage:
  python experiments/measure_orosomucoid_label_stake.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from experiments.measure_fraction_unbound_coverage import measurement_rows, study_drugs

HPO_ANNOTATIONS = Path("data/raw/hpo/genes_to_phenotype.txt")
LABEL_SELECTION = Path("data/processed/label_selection/better_v2_full_v2.parquet")
OUTPUT_DOCUMENT = Path("docs/orosomucoid_label_stake.md")
# every plasma binder the carriage table names (docs/plasma_protein_binding.md), orosomucoid's two genes first
PLASMA_BINDER_GENES = ["ORM1", "ORM2", "ALB", "SERPINA6", "SERPINA7", "TTR", "AFP", "SHBG", "RBP4", "APOA1", "APOB"]


def orosomucoid_carried_drugs(database: Path, evidence_directory: Path, lockbox: Path) -> pd.DataFrame:
    """The study's drug perturbations whose measured rows name orosomucoid as a major binding protein."""
    measurements = measurement_rows(database)
    drugs = study_drugs(evidence_directory, lockbox)
    carriers_of_drug = measurements.groupby("drug_lower").carrier.agg(set)
    carried = drugs[drugs.database_name.map(lambda name: "orosomucoid" in carriers_of_drug.get(name, set())
                                            or "both" in carriers_of_drug.get(name, set()))]
    return carried.reset_index(drop=True)


def label_stake(carried: pd.DataFrame, selection: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series, dict]:
    """Per drug and per symptom, the kept positive pairs the carried drugs hold, against the whole selection."""
    rows = selection[selection.perturbation_id.isin(set(carried.perturbation_id))]
    kept = rows[rows.keep]
    per_drug = (rows.groupby("perturbation_id")
                .agg(selection_rows=("keep", "size"), kept_positives=("keep", "sum"))
                .join(carried.set_index("perturbation_id")[["perturbation_label", "in_lockbox"]])
                .sort_values("kept_positives", ascending=False))
    per_drug["symptoms"] = [", ".join(sorted(rows.symptom[(rows.perturbation_id == identifier) & rows.keep])) or "none"
                            for identifier in per_drug.index]
    totals = {"kept_positives_of_the_carried_drugs": int(kept.keep.sum()),
              "kept_positives_in_the_study": int(selection.keep.sum()),
              "carried_drugs": int(len(carried)),
              "carried_drugs_with_no_kept_positive": int((per_drug.kept_positives == 0).sum()),
              "carried_drugs_in_the_lockbox": int(per_drug.in_lockbox.sum())}
    totals["share_of_kept_positives"] = totals["kept_positives_of_the_carried_drugs"] / totals["kept_positives_in_the_study"]
    return per_drug, kept.symptom.value_counts(), totals


def binder_annotation_counts(annotations_file: Path, selection: pd.DataFrame) -> pd.DataFrame:
    """HPO annotation rows per plasma binder gene, with whether the gene is already a perturbation of this study."""
    annotations = pd.read_csv(annotations_file, sep="\t")
    perturbations = set(selection.perturbation_id)
    rows = []
    for gene in PLASMA_BINDER_GENES:
        annotated = annotations[annotations.gene_symbol == gene]
        selected = selection[selection.perturbation_id == gene]
        rows.append({"gene": gene, "hpo_annotation_rows": len(annotated), "hpo_diseases": annotated.disease_id.nunique(),
                     "is_a_perturbation": gene in perturbations, "selection_rows": len(selected),
                     "kept_positives": int(selected.keep.sum())})
    return pd.DataFrame(rows)


def write_document(per_drug: pd.DataFrame, per_symptom: pd.Series, totals: dict, binders: pd.DataFrame) -> None:
    share = 100 * totals["share_of_kept_positives"]
    lines = [
        "# What the orosomucoid-carried drugs are worth to this study",
        "",
        "Generated by `experiments/measure_orosomucoid_label_stake.py`. Do not edit by hand.",
        "",
        f"The pinned fraction-unbound database names orosomucoid as a major binding protein for "
        f"{totals['carried_drugs']} of this study's drug perturbations ({totals['carried_drugs_in_the_lockbox']} in the "
        f"lockbox). Together they hold **{totals['kept_positives_of_the_carried_drugs']} of the study's "
        f"{totals['kept_positives_in_the_study']} kept positive pairs, {share:.1f} percent**, and "
        f"{totals['carried_drugs_with_no_kept_positive']} of them hold none at all.",
        "",
        "## Per drug",
        "",
        "| drug | in lockbox | selection rows | kept positives | symptoms |",
        "| --- | --- | --- | --- | --- |",
    ]
    for identifier, row in per_drug.iterrows():
        lines.append(f"| {row.perturbation_label} ({identifier}) | {'yes' if row.in_lockbox else 'no'} | "
                     f"{row.selection_rows} | {row.kept_positives} | {row.symptoms} |")
    lines += ["", "## Per symptom, over the carried drugs", "", "| symptom | kept positives |", "| --- | --- |"]
    for symptom, number in per_symptom.items():
        lines.append(f"| {symptom} | {number} |")
    lines += [
        "",
        "## Can a binder gene be a training example?",
        "",
        "A gene perturbation's labels come from HPO disease annotations, so a binder with no HPO annotation has no",
        "label and cannot be a training example however welcome it would be.",
        "",
        "| binder gene | HPO annotation rows | HPO diseases | already a perturbation | selection rows | kept positives |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in binders.itertuples(index=False):
        lines.append(f"| {row.gene} | {row.hpo_annotation_rows} | {row.hpo_diseases} | "
                     f"{'yes' if row.is_a_perturbation else 'no'} | {row.selection_rows} | {row.kept_positives} |")
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--database", type=Path, default=Path("data/raw/plasma_protein_binding/fraction_unbound_database.xlsx"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v2.json"))
    parser.add_argument("--label-selection", type=Path, default=LABEL_SELECTION)
    parser.add_argument("--hpo-annotations", type=Path, default=HPO_ANNOTATIONS)
    arguments = parser.parse_args()

    selection = pd.read_parquet(arguments.label_selection)
    carried = orosomucoid_carried_drugs(arguments.database, arguments.evidence_dir, arguments.lockbox)
    per_drug, per_symptom, totals = label_stake(carried, selection)
    binders = binder_annotation_counts(arguments.hpo_annotations, selection)
    write_document(per_drug, per_symptom, totals, binders)
    for key, value in totals.items():
        print(f"{key}: {value}")
    print(f"\nORM1 and ORM2 HPO annotation rows: "
          f"{binders.set_index('gene').loc['ORM1', 'hpo_annotation_rows']}, "
          f"{binders.set_index('gene').loc['ORM2', 'hpo_annotation_rows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
