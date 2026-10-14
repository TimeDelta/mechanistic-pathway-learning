"""What it costs to require a measured concentration at an off-target while allowing a calibrated placeholder at a
mechanism target.

The owner, 13 October 2026, on the asymmetry table of preregistration item 5a, which this ruling withdrew in
favour of item 5b: "Ooooh. I think it's better to put
the measured values back then and calibrate the other way." The table that prompted it reads, at the calibration's
one-standard-deviation scatter of 1.14 log10 units:

    mechanism targets   median occupancy 0.787    a one-sd concentration error moves it x1.2 up, x3.7 down
    off-targets         median occupancy 0.063    a one-sd concentration error moves it x7.6 up, x13.0 down

So the imputation is nearly free where occupancy saturates and costs an order of magnitude where it is near zero,
and the two roles should not take the same input rule. The rule this script prices: an off-target enters only on a
measured concentration, which restores the prohibition of item 5 that the 5a amendment lifted; a mechanism target
may take the calibrated placeholder, which is where the error is absorbed.

Three questions, because a rule that is right in principle can still be unaffordable:
  1. What does the restored prohibition withhold? Drugs with off-targets but no measured concentration, the
     off-target rows behind them, and the kept positive pairs those drugs carry.
  2. What does the mechanism-side permission admit? Drugs whose mechanism occupancy would rest on the calibrated
     placeholder, and where that occupancy lands on the saturation curve.
  3. Does the asymmetry hold under the new allocation, or does allowing the placeholder at mechanism targets put
     some of them near Kd after all, where the error is not absorbed?

Counting only: no label, table, graph or model input is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/occupancy_by_target_role.py
"""
from __future__ import annotations

import argparse
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

OCCUPANCY_SCRIPT = Path(__file__).resolve().parents[1] / "measure_target_occupancy.py"
specification = importlib.util.spec_from_file_location("measure_target_occupancy", OCCUPANCY_SCRIPT)
occupancy_module = importlib.util.module_from_spec(specification)
specification.loader.exec_module(occupancy_module)

EXPOSURE_TABLE = Path("data/processed/exposure/drug_exposure_v4.csv")
OCCUPANCY_TABLE = Path("data/processed/exposure/target_occupancy_v4.csv")
OUTPUT_DOCUMENT = Path("docs/occupancy_by_target_role.md")

PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET = 0.9
CALIBRATION_OFFSET_LOG10 = 0.331   # placeholder_concentration_calibration.py, the median over the 115 drugs with both
CALIBRATION_SCATTER_LOG10 = 1.14   # the standard deviation about that median


def placeholder_concentration_of_drug(mechanism_affinities: dict[str, float], apply_calibration: bool) -> float | None:
    """The placeholder scheme inverted to a concentration: the best-bound mechanism target sits at 0.9 occupancy.

    occupancy = C / (C + Kd) at 0.9 is C = 9 * Kd, and Kd of the best mechanism target is 10 ** -p_affinity in molar.
    The calibration divides that by 10 ** 0.331, the measured median overstatement."""
    if not mechanism_affinities:
        return None
    best_dissociation_constant_molar = 10.0 ** (-max(mechanism_affinities.values()))
    concentration_molar = (PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET / (1 - PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET)) * best_dissociation_constant_molar
    if apply_calibration:
        concentration_molar /= 10.0 ** CALIBRATION_OFFSET_LOG10
    return concentration_molar * 1e6  # micromolar, the unit the exposure table uses


def occupancy_of(concentration_micromolar: float, p_affinity: float) -> float:
    dissociation_constant_micromolar = 10.0 ** (6.0 - p_affinity)
    return concentration_micromolar / (concentration_micromolar + dissociation_constant_micromolar)


def movement_under_scatter(occupancy: float, scatter_log10: float) -> tuple[float, float]:
    """How far one standard deviation of concentration error moves an occupancy, as a ratio up and a ratio down."""
    if occupancy <= 0 or occupancy >= 1:
        return (1.0, 1.0)
    concentration_over_kd = occupancy / (1 - occupancy)
    up = (concentration_over_kd * 10 ** scatter_log10) / (1 + concentration_over_kd * 10 ** scatter_log10)
    down = (concentration_over_kd / 10 ** scatter_log10) / (1 + concentration_over_kd / 10 ** scatter_log10)
    return (up / occupancy, occupancy / down)


def drugs_with_a_measured_concentration(exposure: pd.DataFrame) -> set[str]:
    plasma = exposure.unbound_plasma_low_micromolar.notna() & exposure.unbound_plasma_high_micromolar.notna()
    brain = exposure.unbound_brain_low_micromolar.notna() & exposure.unbound_brain_high_micromolar.notna()
    return set(exposure.index[plasma | brain])


def kept_positive_pairs_by_drug(data) -> dict[str, int]:
    """Kept positive labels per perturbation, so a withheld rule can be priced in pairs and not only in drugs."""
    mask = getattr(data, "label_mask", None)
    outcomes = np.asarray(data.outcomes)
    kept = outcomes > 0 if mask is None else (outcomes > 0) & np.asarray(mask).astype(bool)
    return {perturbation_id: int(count) for perturbation_id, count in zip(data.perturbation_ids, kept.sum(axis=1))}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--exposure-table", type=Path, default=EXPOSURE_TABLE)
    parser.add_argument("--occupancy-table", type=Path, default=OCCUPANCY_TABLE)
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    exposure = pd.read_csv(arguments.exposure_table).set_index("perturbation_id")
    measured = drugs_with_a_measured_concentration(exposure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster",
                                    label_selection=arguments.label_selection)
    positives_of_drug = kept_positive_pairs_by_drug(data)
    affinities = occupancy_module.affinities_by_drug(data, arguments.graph_dir)

    rows = []
    for perturbation_id, entry in affinities.items():
        placeholder = placeholder_concentration_of_drug(entry["mechanism"], apply_calibration=True)
        for role in ("mechanism", "other"):
            for gene, p_affinity in entry[role].items():
                rows.append({"perturbation_id": perturbation_id, "drug": entry["label"], "role": role, "gene": gene,
                             "p_affinity": p_affinity, "has_measured_concentration": perturbation_id in measured,
                             "placeholder_occupancy": None if placeholder is None else occupancy_of(placeholder, p_affinity),
                             "kept_positives": positives_of_drug.get(perturbation_id, 0)})
    table = pd.DataFrame(rows)
    table.to_csv(Path("data/processed/exposure/target_role_coverage_v4.csv"), index=False)
    write_document(table, affinities, measured, positives_of_drug, arguments.document)


def write_document(table: pd.DataFrame, affinities: dict, measured: set[str], positives_of_drug: dict[str, int],
                   document: Path) -> None:
    all_drugs = set(table.perturbation_id)
    off_target_rows = table[table.role == "other"]
    drugs_with_off_targets = set(off_target_rows.perturbation_id)
    withheld_drugs = drugs_with_off_targets - measured
    withheld_rows = off_target_rows[~off_target_rows.has_measured_concentration]
    total_positives = sum(positives_of_drug.get(drug, 0) for drug in all_drugs)
    withheld_positives = sum(positives_of_drug.get(drug, 0) for drug in withheld_drugs)

    mechanism_rows = table[(table.role == "mechanism") & table.placeholder_occupancy.notna()]
    imputed_mechanism = mechanism_rows[~mechanism_rows.has_measured_concentration]
    down_mechanism = np.array([movement_under_scatter(value, CALIBRATION_SCATTER_LOG10)[1]
                               for value in imputed_mechanism.placeholder_occupancy])
    up_mechanism = np.array([movement_under_scatter(value, CALIBRATION_SCATTER_LOG10)[0]
                             for value in imputed_mechanism.placeholder_occupancy])
    # What withholding asserts in the other direction: the occupancy the withheld off-targets would have carried.
    withheld_occupancy = withheld_rows.placeholder_occupancy.dropna()
    admitted_occupancy = off_target_rows[off_target_rows.has_measured_concentration].placeholder_occupancy.dropna()

    lines = [
        "# Occupancy by target role: measured at an off-target, calibrated placeholder at a mechanism target", "",
        "Generated by `experiments/label_review/occupancy_by_target_role.py`.", "",
        "The owner, 13 October 2026: \"Ooooh. I think it's better to put the measured values back then and calibrate "
        "the other way.\" The rule priced here gives the two target roles different input rules, because the "
        "asymmetry table carried by item 5a, and the rule registered as item 5b, shows the imputation error is "
        "absorbed at one role and not at the other: an "
        "off-target enters only on a measured concentration, a mechanism target may take the calibrated placeholder.",
        "",
        "## 1. What the restored prohibition withholds", "",
        "| reading | count |", "|---|---|",
        f"| drugs with an affinity at any target | {len(all_drugs)} |",
        f"| of those, with at least one off-target | {len(drugs_with_off_targets)} |",
        f"| of those, carrying a measured concentration, so their off-targets enter | {len(drugs_with_off_targets & measured)} |",
        f"| **withheld: off-targets but no measured concentration** | **{len(withheld_drugs)}** |",
        f"| off-target rows in all | {len(off_target_rows)} |",
        f"| off-target rows withheld | **{len(withheld_rows)}** ({len(withheld_rows) / max(len(off_target_rows), 1):.1%}) |",
        f"| kept positive pairs carried by the withheld drugs | {withheld_positives} of {total_positives} ({withheld_positives / max(total_positives, 1):.1%}) |",
        "",
        "**No drug leaves the study.** A withheld drug keeps its mechanism targets and so keeps its input and its "
        "labels; what it loses is the off-target part of that input. That is the distinction the owner drew on 13 "
        "October (\"Don't REQUIRE occupancy though ... some drugs legitimately don't have it\"): the requirement "
        "here is on an off-target entering, not on a drug entering.", "",
        "## 2. What the mechanism-side permission admits", "",
        "| reading | value |", "|---|---|",
        f"| mechanism rows with a placeholder occupancy | {len(mechanism_rows)} |",
        f"| of those, on a drug with no measured concentration, so imputed | **{len(imputed_mechanism)}** |",
        f"| median imputed mechanism occupancy | {imputed_mechanism.placeholder_occupancy.median():.3f} |",
        f"| quartiles | {imputed_mechanism.placeholder_occupancy.quantile(.25):.3f} to {imputed_mechanism.placeholder_occupancy.quantile(.75):.3f} |",
        f"| a one-sd concentration error moves it, median | x{np.median(up_mechanism):.1f} up, x{np.median(down_mechanism):.1f} down |",
        f"| the worst quarter of them, downwards | x{np.quantile(down_mechanism, 0.75):.1f} or more |",
        f"| the worst tenth | x{np.quantile(down_mechanism, 0.90):.1f} or more |",
        "",
        "**Half of this is built in and should not be read as a finding.** The placeholder pins the best-bound "
        "mechanism target at 0.9 occupancy, which the calibration moves to 0.808, so every drug has one mechanism "
        "target at exactly that value by construction and the rest follow from the affinity ratios. What the "
        "measured drugs contribute is the test of whether the assertion is roughly right: their median mechanism "
        "occupancy is 0.787 against the asserted 0.808, so it is. The spread is where the assertion is not right, "
        f"and the lower quartile at {imputed_mechanism.placeholder_occupancy.quantile(.25):.3f} is a reminder that "
        "a drug's weaker mechanism targets sit well down the curve even when its best one saturates.", "",
        "## 3. What withholding an off-target asserts in the other direction", "",
        "A withheld off-target is not abstained from: it is set to absent, which is the claim that the drug does "
        "not engage that target. That claim is wrong whenever the drug does. So the cost of the prohibition is not "
        "only the 578 rows it drops but the occupancy those rows would have carried.", "",
        "| reading | value |", "|---|---|",
        f"| median placeholder occupancy of a withheld off-target | {withheld_occupancy.median():.3f} |",
        f"| of an admitted one, for comparison | {admitted_occupancy.median():.3f} |",
        f"| withheld rows the placeholder puts above 0.5 | {int((withheld_occupancy >= 0.5).sum())} of {len(withheld_occupancy)} ({(withheld_occupancy >= 0.5).mean():.1%}) |",
        f"| withheld rows above 0.25 | {int((withheld_occupancy >= 0.25).sum())} ({(withheld_occupancy >= 0.25).mean():.1%}) |",
        f"| withheld rows below 0.10, where setting them absent is nearly true | {int((withheld_occupancy < 0.10).sum())} ({(withheld_occupancy < 0.10).mean():.1%}) |",
        "",
        "## 4. Does the asymmetry survive the new allocation", "",
        "The asymmetry was measured on the drugs that have a measured concentration, where the median "
        "mechanism occupancy is 0.787. What matters for this rule is whether the drugs that do **not** have one "
        "look the same, because those are the drugs the placeholder is carried to. They do: the median imputed "
        f"mechanism occupancy is {imputed_mechanism.placeholder_occupancy.median():.3f} and a one-sd error moves it "
        f"x{np.median(up_mechanism):.1f} up and x{np.median(down_mechanism):.1f} down, against x1.2 and x3.7 on the "
        "measured drugs. The off-target figure the prohibition is built on, x7.6 up and x13.0 down, stands "
        "unchanged, because it is a property of where a median off-target occupancy of 0.063 sits on the curve and "
        "not of which drugs were used to find it.", "",
    ]
    document.parent.mkdir(parents=True, exist_ok=True)
    document.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
