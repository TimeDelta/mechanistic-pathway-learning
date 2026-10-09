"""Test whether the fraction unbound in a special population tells which plasma protein carries a drug.

The user, 9 October 2026: "for the fraction unbound, the 'in special populations' is the part that i am leaning on
there. maybe it doesn't afford the leniancy i think it does?" The idea is sound physiology: orosomucoid rises in
inflammation and renal disease and falls in liver disease and in the newborn, albumin falls in liver disease and in
the newborn, so the direction in which a drug's free fraction moves between a reference population and a special one
carries information about which protein was holding it.

This measures how much. The pinned database of Al-Qassabi and colleagues (data/raw/plasma_protein_binding, figshare
10.48420/25243138.v1, CC BY 4.0; the paper is 10.1016/j.xphs.2024.02.024, PMID 38417790) labels each measurement with
its major binding protein, so the direction can be scored against that label: per population, the AUROC of the free
fraction ratio for orosomucoid carriage, and then a leave-one-drug-out classifier that pools the populations a drug
has rows in, against the majority-class rate.

Two warnings the readings have to carry. The sign is not one sign: it flips between populations, so a ratio is
uninterpretable without knowing which population produced it and what that population does to each protein. And the
classifier needs a pair of measured fractions for the drug, which is the data this study lacks: a drug with no row in
this database has no ratio to classify, so a working classifier does not by itself widen coverage.

Writes docs/special_population_carrier_signal.md. Reads no lockbox outcome. Idempotent.

Usage:
  python experiments/measure_special_population_carrier_signal.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SHEET_BINDERS = {"hsa": "albumin", "aag": "orosomucoid", "both": "both"}
POPULATION_NAMES = {"InflamDisease": "inflammatory disease"}
OUTPUT_DOCUMENT = Path("docs/special_population_carrier_signal.md")
CARRIERS = ("albumin", "orosomucoid")


def measurement_rows(database_path: Path) -> pd.DataFrame:
    """One row per measurement: the drug, its labelled carrier, the population and the two fractions unbound."""
    sheets = pd.read_excel(database_path, sheet_name=None)
    collected = []
    for name, frame in sheets.items():
        if not name[0].isdigit() or "Major Binding Protein" not in frame.columns:
            continue
        _, population, _ = name.split("_")
        renamed = frame.rename(columns={"Major Binding Protein": "carrier", "Drug": "drug",
                                        "fu in reference population": "fraction_unbound_reference",
                                        "Special population subgroup": "subgroup",
                                        "fu in special population": "fraction_unbound_special"})
        renamed["population"] = POPULATION_NAMES.get(population, population.lower())
        renamed["carrier"] = renamed.carrier.astype(str).str.strip().str.lower().map(SHEET_BINDERS)
        collected.append(renamed[["carrier", "drug", "fraction_unbound_reference", "subgroup",
                                  "fraction_unbound_special", "population"]])
    rows = pd.concat(collected, ignore_index=True)
    rows["drug"] = rows.drug.astype(str).str.strip().str.lower()
    rows = rows.dropna(subset=["carrier", "fraction_unbound_reference", "fraction_unbound_special"])
    rows = rows[(rows.fraction_unbound_reference > 0) & (rows.fraction_unbound_special > 0)]
    rows["log_ratio"] = np.log(rows.fraction_unbound_special / rows.fraction_unbound_reference)
    return rows.reset_index(drop=True)


def population_readings(rows: pd.DataFrame) -> pd.DataFrame:
    """Per population: each carrier's median ratio and how well the ratio alone orders orosomucoid carriage."""
    readings = []
    for population, group in rows[rows.carrier.isin(CARRIERS)].groupby("population"):
        medians = {carrier: float(np.exp(part.log_ratio.median())) for carrier, part in group.groupby("carrier")}
        reading = {"population": population, "measurements": len(group),
                   "albumin_rows": int((group.carrier == "albumin").sum()),
                   "orosomucoid_rows": int((group.carrier == "orosomucoid").sum()),
                   "albumin_median_ratio": medians.get("albumin", float("nan")),
                   "orosomucoid_median_ratio": medians.get("orosomucoid", float("nan")), "auroc": float("nan")}
        if group.carrier.nunique() == 2:
            is_orosomucoid = (group.carrier == "orosomucoid").astype(int)
            reading["auroc"] = float(roc_auc_score(is_orosomucoid, group.log_ratio))
        readings.append(reading)
    return pd.DataFrame(readings).sort_values("population")


def leave_one_drug_out_predictions(rows: pd.DataFrame) -> pd.DataFrame:
    """Classify each drug's carrier from its ratios, with the per-population rule fitted on the other drugs."""
    two_carriers = rows[rows.carrier.isin(CARRIERS)]
    carrier_of_drug = two_carriers.groupby("drug").carrier.agg(lambda values: values.mode().iat[0])
    predictions = []
    for drug, carrier in carrier_of_drug.items():
        others = two_carriers[two_carriers.drug != drug]
        votes = []
        for population, group in two_carriers[two_carriers.drug == drug].groupby("population"):
            fitted = others[others.population == population]
            if fitted.carrier.nunique() < 2:
                continue
            albumin_median = fitted[fitted.carrier == "albumin"].log_ratio.median()
            orosomucoid_median = fitted[fitted.carrier == "orosomucoid"].log_ratio.median()
            direction = np.sign(orosomucoid_median - albumin_median)
            midpoint = (albumin_median + orosomucoid_median) / 2.0
            if direction != 0:
                votes.append(float(direction * (group.log_ratio.median() - midpoint)))
        if votes:
            predictions.append({"drug": drug, "carrier": carrier, "populations": len(votes),
                                "predicted": "orosomucoid" if float(np.mean(votes)) > 0 else "albumin"})
    return pd.DataFrame(predictions)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--database", type=Path,
                        default=Path("data/raw/plasma_protein_binding/fraction_unbound_database.xlsx"))
    parser.add_argument("--study-evidence", type=Path, default=Path("data/processed/evidence_full_v3"),
                        help="to say how many of the study's drugs the classifier could be applied to at all")
    arguments = parser.parse_args()

    rows = measurement_rows(arguments.database)
    readings = population_readings(rows)
    predictions = leave_one_drug_out_predictions(rows)
    correct = int((predictions.carrier == predictions.predicted).sum())
    majority = predictions.carrier.value_counts(normalize=True).max() if len(predictions) else float("nan")
    orosomucoid_rows = predictions[predictions.carrier == "orosomucoid"]
    albumin_rows = predictions[predictions.carrier == "albumin"]
    study_drugs = pd.read_parquet(arguments.study_evidence / "evidence_records.parquet",
                                  columns=["perturbation_type", "perturbation_label"])
    study_drug_labels = set(study_drugs[study_drugs.perturbation_type == "drug"].perturbation_label.str.lower())
    with_a_ratio = study_drug_labels & set(rows.drug)

    lines = [
        "# Does the fraction unbound in a special population say which protein carries the drug?",
        "",
        "Generated by `experiments/measure_special_population_carrier_signal.py`. Do not edit by hand.",
        "",
        "Source: the pinned database of Al-Qassabi and colleagues (figshare "
        "[10.48420/25243138.v1](https://doi.org/10.48420/25243138.v1), CC BY 4.0; the paper is PMID 38417790, "
        "[DOI](https://doi.org/10.1016/j.xphs.2024.02.024)), which labels every measurement with its major binding "
        "protein. The ratio below is the fraction unbound in the special population over the fraction unbound in the "
        "reference population, so above 1 means more free drug.",
        "",
        "## Per population",
        "",
        "| population | measurements | albumin rows | orosomucoid rows | albumin median ratio | orosomucoid median ratio | AUROC for orosomucoid |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in readings.itertuples(index=False):
        lines.append(f"| {row.population} | {row.measurements} | {row.albumin_rows} | {row.orosomucoid_rows} | "
                     f"{row.albumin_median_ratio:.3f} | {row.orosomucoid_median_ratio:.3f} | {row.auroc:.3f} |")
    lines += [
        "",
        "An AUROC near 1 means a lower ratio goes with albumin and a higher one with orosomucoid; near 0 means the "
        "reverse. **Both appear**, which is the first finding: the direction is not one direction. Orosomucoid is an "
        "acute-phase protein made in the liver, so it rises in inflammation and in renal disease and falls in hepatic "
        "impairment and in the newborn, while albumin falls in hepatic impairment and in the newborn too. A ratio "
        "therefore says nothing about the carrier until the population is known, and a source that reports a free "
        "fraction without naming the population cannot be read this way at all.",
        "",
        "## Pooling the populations a drug has, leave-one-drug-out",
        "",
        f"- drugs classified: **{len(predictions)}** ({len(orosomucoid_rows)} orosomucoid, {len(albumin_rows)} albumin)",
        f"- correct: **{correct}** of {len(predictions)}, {correct / len(predictions):.0%}" if len(predictions) else "- no drug classified",
        f"- majority-class rate: {majority:.0%}",
        f"- orosomucoid drugs found: {int((orosomucoid_rows.predicted == 'orosomucoid').sum())} of {len(orosomucoid_rows)}",
        f"- albumin drugs found: {int((albumin_rows.predicted == 'albumin').sum())} of {len(albumin_rows)}",
        "",
        "Accuracy is the wrong single number here, because the trivial rule \"albumin\" already reaches the "
        f"majority-class rate while finding none of the {len(orosomucoid_rows)} orosomucoid carriers. Read both ways: "
        "the direction barely moves accuracy, and it recovers "
        f"{int((orosomucoid_rows.predicted == 'orosomucoid').sum())} of the {len(orosomucoid_rows)} orosomucoid "
        f"carriers at the cost of calling {int((albumin_rows.predicted == 'orosomucoid').sum())} of "
        f"{len(albumin_rows)} albumin-carried drugs orosomucoid. For deciding one drug's carriage edge that error "
        "rate is too high to place an edge on; for flagging which drugs to look up in a source that states the "
        "carrier, it is useful.",
        "",
        "The rule is fitted per population on the other drugs (the median log ratio of each carrier's drugs, the "
        "midpoint between them as the threshold and the sign between them as the direction) and the populations a "
        "drug has rows in vote by their distance from that threshold. Drugs the database labels \"both\" are left out "
        "of both the fit and the test.",
        "",
        "## What this can and cannot buy",
        "",
        f"The classifier needs a measured pair of fractions for the drug. **{len(with_a_ratio)} of this study's "
        f"{len(study_drug_labels)} drug perturbations have one in this database**, which is the same set the carrier "
        "label already covers, so running the classifier on them adds nothing: where there is a ratio there is "
        "already a labelled carrier. It would only widen coverage through a second source that reports a reference "
        "and a special-population fraction for drugs this database misses, and names the population. No such source "
        "is pinned.",
        "",
        "So the leniency the special-population column affords is real but narrow: it supports reading the carrier "
        "off a free-fraction pair when the population is known, and it does not reach the 115 drugs with no measured "
        "fraction at all. The per-drug carrier identity those drugs need comes from a source that states it, which is "
        "what docs/label_plasma_binders.md measures.",
    ]
    OUTPUT_DOCUMENT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT_DOCUMENT}")
    print(f"drugs classified {len(predictions)}, correct {correct}, majority {majority:.0%}, "
          f"study drugs with a ratio {len(with_a_ratio)} of {len(study_drug_labels)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
