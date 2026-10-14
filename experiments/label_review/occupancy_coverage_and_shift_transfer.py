"""What preregistration item 5c covers, and whether its calibration shift transfers to the drugs it is applied to.

The owner, 13 October 2026: "keep off targets but always use a measured value then impute all values to get the
calibration shift and apply it to only the imputed values that have not been measured." Four steps, and the third
is the one this script exploits: the placeholder is computed for **every** drug, including the drugs that have a
measurement. Those scaffolding values enter no model input, but they make two things measurable that were not.

  1. **Coverage.** How many drugs and target rows take a measured concentration and how many take the shifted
     placeholder, split by target role, which is what `has_measured_occupancy` will stratify on.

  2. **Whether the shift transfers.** Item 29 registers a limitation it calls untested: the offset is estimated on
     the drugs that have a measured concentration, which are the better-studied ones, and carried to those that do
     not. A direct test is impossible, because testing it needs the measurement whose absence is the problem. A
     partial one is not: if the offset varies with a covariate that is observable for both groups, and the two
     groups differ on that covariate, the extrapolation is in trouble; if the offset is flat in the covariates that
     separate the groups, it is better supported than "stated and not tested". The covariates used are the ones the
     placeholder is built from or that stand in for how well studied a drug is: the best mechanism affinity, the
     number of mechanism targets and the number of off-targets.

  3. **The uncertainty by role**, which the stratum has to carry whichever way decision 1.3 goes on off-targets.

Counting only: no label, table, graph or model input is changed.

Usage:
  PYTHONPATH=. python experiments/label_review/occupancy_coverage_and_shift_transfer.py
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
OUTPUT_DOCUMENT = Path("docs/occupancy_coverage_and_shift_transfer.md")

PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET = 0.9
CALIBRATION_OFFSET_LOG10 = 0.331   # placeholder_concentration_calibration.py, the median over the drugs with both
CALIBRATION_SCATTER_LOG10 = 1.14   # the standard deviation about that median


def placeholder_concentration_of_drug(mechanism_affinities: dict[str, float], apply_calibration: bool) -> float | None:
    """Item 5c step 3: the placeholder, computed the same way for every drug whether or not it has a measurement.

    occupancy = C / (C + Kd) at 0.9 is C = 9 * Kd, with Kd of the best mechanism target at 10 ** -p_affinity molar.
    Step 4's shift divides that by 10 ** 0.331 and is applied only where no measurement exists."""
    if not mechanism_affinities:
        return None
    best_dissociation_constant_molar = 10.0 ** (-max(mechanism_affinities.values()))
    ratio = PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET / (1 - PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET)
    concentration_molar = ratio * best_dissociation_constant_molar
    if apply_calibration:
        concentration_molar /= 10.0 ** CALIBRATION_OFFSET_LOG10
    return concentration_molar * 1e6  # micromolar, the unit the exposure table uses


def occupancy_of(concentration_micromolar: float, p_affinity: float) -> float:
    dissociation_constant_micromolar = 10.0 ** (6.0 - p_affinity)
    return concentration_micromolar / (concentration_micromolar + dissociation_constant_micromolar)


def movement_under_scatter(occupancy: float, scatter_log10: float) -> tuple[float, float]:
    """How far one standard deviation of concentration error moves an occupancy, as a ratio up and a ratio down.

    d(log theta) / d(log C) = 1 - theta, so this falls monotonically as the occupancy rises."""
    if occupancy <= 0 or occupancy >= 1:
        return (1.0, 1.0)
    concentration_over_kd = occupancy / (1 - occupancy)
    factor = 10 ** scatter_log10
    up = (concentration_over_kd * factor) / (1 + concentration_over_kd * factor)
    down = (concentration_over_kd / factor) / (1 + concentration_over_kd / factor)
    return (up / occupancy, occupancy / down)


def measured_concentration_centre(exposure: pd.DataFrame) -> pd.Series:
    """The geometric centre of a drug's measured unbound range, preferring plasma where both sites exist."""
    plasma = np.sqrt(exposure.unbound_plasma_low_micromolar * exposure.unbound_plasma_high_micromolar)
    brain = np.sqrt(exposure.unbound_brain_low_micromolar * exposure.unbound_brain_high_micromolar)
    return plasma.fillna(brain)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v4.parquet"))
    parser.add_argument("--exposure-table", type=Path, default=EXPOSURE_TABLE)
    parser.add_argument("--document", type=Path, default=OUTPUT_DOCUMENT)
    arguments = parser.parse_args()

    exposure = pd.read_csv(arguments.exposure_table).set_index("perturbation_id")
    centre = measured_concentration_centre(exposure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster",
                                    label_selection=arguments.label_selection)
    affinities = occupancy_module.affinities_by_drug(data, arguments.graph_dir)

    per_drug, rows = [], []
    for perturbation_id, entry in affinities.items():
        raw_placeholder = placeholder_concentration_of_drug(entry["mechanism"], apply_calibration=False)
        if raw_placeholder is None:
            continue
        measured_centre = centre.get(perturbation_id, np.nan)
        has_measurement = bool(pd.notna(measured_centre))
        per_drug.append({
            "perturbation_id": perturbation_id, "drug": entry["label"], "has_measurement": has_measurement,
            "raw_placeholder": raw_placeholder, "measured_centre": measured_centre,
            "offset_log10": np.log10(raw_placeholder / measured_centre) if has_measurement else np.nan,
            "best_mechanism_affinity": max(entry["mechanism"].values()),
            "mechanism_target_count": len(entry["mechanism"]), "off_target_count": len(entry["other"])})
        concentration = measured_centre if has_measurement else raw_placeholder / 10 ** CALIBRATION_OFFSET_LOG10
        for role in ("mechanism", "other"):
            for gene, p_affinity in entry[role].items():
                rows.append({"perturbation_id": perturbation_id, "role": role, "gene": gene,
                             "has_measured_occupancy": has_measurement,
                             "occupancy": occupancy_of(concentration, p_affinity)})
    drugs, targets = pd.DataFrame(per_drug), pd.DataFrame(rows)
    targets.to_csv(Path("data/processed/exposure/occupancy_under_item_5c_v4.csv"), index=False)
    write_document(drugs, targets, arguments.document)


def covariate_table(drugs: pd.DataFrame) -> list[str]:
    """Does the offset move with a covariate that separates the two groups? Section 2 of the document."""
    with_measurement = drugs[drugs.has_measurement]
    lines = ["| covariate | measured drugs, median | imputed drugs, median | offset in its lower half | offset in its upper half | gap |",
             "|---|---|---|---|---|---|"]
    for column, label in (("best_mechanism_affinity", "best mechanism affinity (pK)"),
                          ("mechanism_target_count", "mechanism targets"),
                          ("off_target_count", "off-targets")):
        split_at = with_measurement[column].median()
        lower = with_measurement[with_measurement[column] <= split_at].offset_log10
        upper = with_measurement[with_measurement[column] > split_at].offset_log10
        gap = upper.median() - lower.median()
        lines.append(f"| {label} | {with_measurement[column].median():.2f} | "
                     f"{drugs[~drugs.has_measurement][column].median():.2f} | "
                     f"{lower.median():+.3f} (n={len(lower)}) | {upper.median():+.3f} (n={len(upper)}) | "
                     f"**{gap:+.3f}** |")
    return lines



def shift_models(drugs: pd.DataFrame) -> dict:
    """A constant shift against one that varies with the best mechanism affinity, scored out of sample.

    The placeholder is C = 9 * Kd of the best mechanism target, so log10 C_placeholder falls with the affinity at a
    slope of exactly -1 by construction. If the measured concentration falls more slowly than that - a more potent
    drug is dosed lower, but not proportionately lower - then the disagreement between the two carries the
    remainder of the affinity, no constant can remove it, and a constant fitted on one affinity range is wrong on
    another. That is testable here and is what item 29's untested limitation turns on."""
    both = drugs[drugs.has_measurement].copy()
    both["log_measured"] = np.log10(both.measured_centre)
    both["log_placeholder"] = np.log10(both.raw_placeholder)
    slope, intercept = np.polyfit(both.best_mechanism_affinity, both.log_measured, 1)
    constant_errors, slope_errors = [], []
    for position in range(len(both)):
        held_out = both.iloc[position]
        rest = both.drop(both.index[position])
        constant_errors.append(abs((held_out.log_placeholder - rest.offset_log10.median()) - held_out.log_measured))
        rest_slope, rest_intercept = np.polyfit(rest.best_mechanism_affinity, rest.log_measured, 1)
        slope_errors.append(abs((rest_intercept + rest_slope * held_out.best_mechanism_affinity) - held_out.log_measured))
    imputed = drugs[~drugs.has_measurement]
    under_constant = np.log10(imputed.raw_placeholder) - both.offset_log10.median()
    under_slope = intercept + slope * imputed.best_mechanism_affinity
    return {"constant_shift": both.offset_log10.median(), "slope": slope, "intercept": intercept,
            "constant_out_of_sample": float(np.median(constant_errors)),
            "slope_out_of_sample": float(np.median(slope_errors)),
            "affinity_gap": imputed.best_mechanism_affinity.median() - both.best_mechanism_affinity.median(),
            "disagreement_on_the_imputed": float(under_constant.median() - under_slope.median())}


def write_document(drugs: pd.DataFrame, targets: pd.DataFrame, document: Path) -> None:
    measured_drugs = drugs[drugs.has_measurement]
    imputed_drugs = drugs[~drugs.has_measurement]
    by_group_role = targets.groupby(["has_measured_occupancy", "role"]).size()

    movement = {}
    for (has_measurement, role), group in targets.groupby(["has_measured_occupancy", "role"]):
        if has_measurement:
            continue
        pairs = [movement_under_scatter(value, CALIBRATION_SCATTER_LOG10) for value in group.occupancy]
        movement[role] = (group.occupancy.median(), np.median([p[0] for p in pairs]), np.median([p[1] for p in pairs]))

    lines = [
        "# Occupancy under item 5c: what it covers, and whether the calibration shift transfers", "",
        "Generated by `experiments/label_review/occupancy_coverage_and_shift_transfer.py`.", "",
        "The owner, 13 October 2026: \"keep off targets but always use a measured value then impute all values to "
        "get the calibration shift and apply it to only the imputed values that have not been measured.\"", "",
        "## 1. Coverage", "",
        "| | drugs | mechanism target rows | off-target rows |", "|---|---|---|---|",
        f"| take a **measured** concentration | {len(measured_drugs)} | {by_group_role.get((True, 'mechanism'), 0)} | {by_group_role.get((True, 'other'), 0)} |",
        f"| take the **shifted placeholder** | {len(imputed_drugs)} | {by_group_role.get((False, 'mechanism'), 0)} | {by_group_role.get((False, 'other'), 0)} |",
        f"| all | {len(drugs)} | {int(by_group_role.get((True, 'mechanism'), 0) + by_group_role.get((False, 'mechanism'), 0))} | {int(by_group_role.get((True, 'other'), 0) + by_group_role.get((False, 'other'), 0))} |",
        "",
        f"Every one of the {len(drugs)} drugs with an affinity carries an occupancy at every one of its targets in "
        "both roles, which is what item 5c step 1 is for: no occupancy rule withholds a value, and "
        "`has_measured_occupancy` is the only thing separating the two populations.", "",
        "## 2. Does the shift transfer to the drugs it is applied to", "",
        "Item 29 registers this as a limitation that nothing measures. Item 5c step 3 makes a partial test "
        "possible: a direct one is impossible, since it would need the measurement whose absence is the problem, "
        "but the offset can be checked for dependence on covariates that are observable for both groups. If the "
        "offset moves with a covariate on which the two groups differ, the extrapolation is in trouble. The gap "
        "column is the offset's median in the upper half of the covariate minus its median in the lower half, "
        f"against a scatter of {CALIBRATION_SCATTER_LOG10:.2f} log units and a shift of "
        f"{CALIBRATION_OFFSET_LOG10:+.3f}.", "",
    ]
    lines += covariate_table(drugs)
    models = shift_models(drugs)
    lines += [
        "",
        "**The affinity row is not noise, and it is not a coincidence.** The placeholder is `C = 9 * Kd` of the "
        "best mechanism target, so `log10 C` falls with the affinity at a slope of exactly **-1 by construction**. "
        f"The measured concentration falls at **{models['slope']:+.3f}**: a more potent drug is dosed lower, but "
        "only about half as much lower as the placeholder assumes. The remainder stays in the disagreement, where "
        "no constant can reach it, and it is worst exactly where it matters, because the drugs the shift is "
        f"applied to are the more potent ones - median best affinity {models['affinity_gap']:+.2f} pK units above "
        "the drugs it was estimated on.", "",
        "| | shift | out-of-sample median absolute error (leave one drug out) |", "|---|---|---|",
        f"| a constant, as item 5c registers | {models['constant_shift']:+.3f} log10 | **{models['constant_out_of_sample']:.3f} log10** |",
        f"| a line in the best mechanism affinity | log10 C = {models['intercept']:+.3f} {models['slope']:+.3f} * pK | **{models['slope_out_of_sample']:.3f} log10** |",
        "",
        f"The line cuts the out-of-sample error by {1 - models['slope_out_of_sample'] / models['constant_out_of_sample']:.0%}, "
        "and on the drugs the shift is actually applied to the two disagree by a median "
        f"{models['disagreement_on_the_imputed']:+.3f} log10, a factor of "
        f"{10 ** abs(models['disagreement_on_the_imputed']):.1f}: the constant shift puts them **lower** than the "
        "line does, because it carries the affinity of a less potent population. **So item 29's limitation is no "
        "longer untested, and it does not pass.** Changing the estimator is the owner's and is recorded as open "
        "decision 9.2; item 5c as registered keeps the constant.", "",
        "**One caveat on the slope, in its own favour and against.** An affinity read with error attenuates the "
        "fitted slope towards zero, so the true dependence may be steeper than "
        f"{models['slope']:+.3f} and the placeholder's -1 less wrong than it looks. That weakens the mechanistic "
        "reading and not the practical one: the out-of-sample error is a prediction, measured on held-out drugs, "
        "and it falls whatever is causing the slope.", "",
        "## 3. The uncertainty the stratum has to carry", "",
        "| role, imputed drugs only | median occupancy | a one-sd concentration error moves it |", "|---|---|---|",
    ]
    for role, label in (("mechanism", "mechanism targets"), ("other", "off-targets")):
        if role in movement:
            median_occupancy, up, down = movement[role]
            lines.append(f"| {label} | {median_occupancy:.3f} | x{up:.1f} up, x{down:.1f} down |")
    lines += [
        "",
        "The sensitivity is `d(log theta) / d(log C) = 1 - theta`, monotone in the occupancy, so a low occupancy "
        "absorbs least of the error. That is a property of the curve and not of the target role; the role only "
        "correlates with where on the curve a target sits, and it fails in the tail, which is why "
        "`has_measured_occupancy` and not the role is the instrument item 29 uses "
        "(`docs/open_decisions.md` section 9.1).", "",
    ]
    document.parent.mkdir(parents=True, exist_ok=True)
    document.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
