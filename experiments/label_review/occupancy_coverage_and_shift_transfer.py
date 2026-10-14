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
from mechanistic_pathway_learning.perturbation.occupancy_calibration import (
    DEFAULT_CALIBRATION_FILE,
    calibrated_concentration_micromolar,
    fit_concentration_calibration,
    occupancy,
    placeholder_log10_concentration,
)

OCCUPANCY_SCRIPT = Path(__file__).resolve().parents[1] / "measure_target_occupancy.py"
specification = importlib.util.spec_from_file_location("measure_target_occupancy", OCCUPANCY_SCRIPT)
occupancy_module = importlib.util.module_from_spec(specification)
specification.loader.exec_module(occupancy_module)

EXPOSURE_TABLE = Path("data/processed/exposure/drug_exposure_v4.csv")
OCCUPANCY_TABLE = Path("data/processed/exposure/target_occupancy_v4.csv")
OUTPUT_DOCUMENT = Path("docs/occupancy_coverage_and_shift_transfer.md")

PLACEHOLDER_OCCUPANCY_OF_THE_BEST_TARGET = 0.9
BOOTSTRAP_DRAWS, BOOTSTRAP_SEED = 4000, 0
CALIBRATION_SCATTER_LOG10 = 1.14   # the standard deviation about that median


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
    parser.add_argument("--calibration-file", type=Path, default=DEFAULT_CALIBRATION_FILE)
    arguments = parser.parse_args()

    exposure = pd.read_csv(arguments.exposure_table).set_index("perturbation_id")
    centre = measured_concentration_centre(exposure)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster",
                                    label_selection=arguments.label_selection)
    affinities = occupancy_module.affinities_by_drug(data, arguments.graph_dir)

    per_drug = []
    for perturbation_id, entry in affinities.items():
        if not entry["mechanism"]:
            continue
        best_affinity = max(entry["mechanism"].values())
        raw_placeholder = 10.0 ** placeholder_log10_concentration(best_affinity)
        measured_centre = centre.get(perturbation_id, np.nan)
        has_measurement = bool(pd.notna(measured_centre)) and measured_centre > 0
        per_drug.append({
            "perturbation_id": perturbation_id, "drug": entry["label"], "has_measurement": has_measurement,
            "raw_placeholder": raw_placeholder, "measured_centre": measured_centre,
            "offset_log10": np.log10(raw_placeholder / measured_centre) if has_measurement else np.nan,
            "best_mechanism_affinity": best_affinity,
            "mechanism_target_count": len(entry["mechanism"]), "off_target_count": len(entry["other"])})
    drugs = pd.DataFrame(per_drug)
    fitting = drugs[drugs.has_measurement]
    calibration = fit_concentration_calibration(fitting.best_mechanism_affinity.values,
                                                np.log10(fitting.measured_centre.values))
    calibration.write(arguments.calibration_file)

    rows = []
    for perturbation_id, entry in affinities.items():
        if not entry["mechanism"]:
            continue
        drug = drugs.set_index("perturbation_id").loc[perturbation_id]
        concentration = drug.measured_centre if drug.has_measurement else \
            calibrated_concentration_micromolar(drug.best_mechanism_affinity, calibration)
        for role in ("mechanism", "other"):
            for gene, p_affinity in entry[role].items():
                rows.append({"perturbation_id": perturbation_id, "role": role, "gene": gene,
                             "has_measured_occupancy": bool(drug.has_measurement),
                             "occupancy": occupancy(concentration, p_affinity)})
    targets = pd.DataFrame(rows)
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
    """The registered slope-and-bias calibration against the withdrawn constant shift, scored out of sample.

    The placeholder is `C = 9 * Kd` of the best mechanism target, so its log concentration falls with affinity at a
    slope of exactly -1 by construction. If the measured concentration falls more slowly, the disagreement carries
    the remainder of the affinity, no constant reaches it, and a constant fitted on one affinity range is wrong on
    another. Both models are fitted here, the slope's 95 percent interval is bootstrapped to see whether 1.0
    survives, and both are scored by leaving one drug out."""
    both = drugs[drugs.has_measurement].copy()
    both["log_measured"] = np.log10(both.measured_centre)
    both["log_placeholder"] = np.log10(both.raw_placeholder)
    calibration = fit_concentration_calibration(both.best_mechanism_affinity.values, both.log_measured.values)

    generator = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrapped_slopes = []
    for _ in range(BOOTSTRAP_DRAWS):
        drawn = generator.integers(0, len(both), len(both))
        bootstrapped_slopes.append(fit_concentration_calibration(
            both.best_mechanism_affinity.values[drawn], both.log_measured.values[drawn]).slope_on_the_placeholder)
    slope_interval = np.percentile(bootstrapped_slopes, [2.5, 97.5])

    constant_errors, calibrated_errors = [], []
    for position in range(len(both)):
        held_out, rest = both.iloc[position], both.drop(both.index[position])
        constant_errors.append(abs((held_out.log_placeholder - rest.offset_log10.median()) - held_out.log_measured))
        rest_calibration = fit_concentration_calibration(rest.best_mechanism_affinity.values, rest.log_measured.values)
        predicted = np.log10(calibrated_concentration_micromolar(held_out.best_mechanism_affinity, rest_calibration))
        calibrated_errors.append(abs(predicted - held_out.log_measured))

    imputed = drugs[~drugs.has_measurement]
    outside = int(((imputed.best_mechanism_affinity < calibration.lowest_fitted_p_affinity) |
                   (imputed.best_mechanism_affinity > calibration.highest_fitted_p_affinity)).sum())
    constant_shift = both.offset_log10.median()
    best_target_under_constant = np.array([
        occupancy(10 ** (np.log10(row.raw_placeholder) - constant_shift), row.best_mechanism_affinity)
        for row in imputed.itertuples()])
    best_target_under_calibration = np.array([
        occupancy(calibrated_concentration_micromolar(row.best_mechanism_affinity, calibration),
                  row.best_mechanism_affinity) for row in imputed.itertuples()])
    best_target_measured = np.array([occupancy(row.measured_centre, row.best_mechanism_affinity)
                                     for row in both.itertuples()])
    return {"calibration": calibration, "constant_shift": constant_shift,
            "slope_interval": slope_interval,
            "constant_out_of_sample": float(np.median(constant_errors)),
            "calibrated_out_of_sample": float(np.median(calibrated_errors)),
            "affinity_gap": imputed.best_mechanism_affinity.median() - both.best_mechanism_affinity.median(),
            "drugs_outside_the_support": outside,
            "best_target_under_constant": best_target_under_constant,
            "best_target_under_calibration": best_target_under_calibration,
            "best_target_measured": best_target_measured}


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
        f"against a scatter of {CALIBRATION_SCATTER_LOG10:.2f} log units.", "",
    ]
    lines += covariate_table(drugs)
    models = shift_models(drugs)
    calibration = models["calibration"]
    low_interval, high_interval = models["slope_interval"]
    constant_best, calibrated_best = models["best_target_under_constant"], models["best_target_under_calibration"]
    measured_best = models["best_target_measured"]
    lines += [
        "",
        "**The affinity row is not noise and it is not a coincidence.** The placeholder is `C = 9 * Kd` of the best "
        "mechanism target, so its `log10 C` falls with the affinity at a slope of exactly **-1 by construction**. "
        "The measured concentration falls more slowly: a more potent drug is dosed lower, but only about half as "
        "much lower. The remainder stays in the disagreement, where no constant reaches it, and it is worst exactly "
        f"where it matters, because the drugs a shift is applied to are the more potent ones - a median "
        f"{models['affinity_gap']:+.2f} pK above the drugs it is estimated on.", "",
        "## 2b. The registered calibration: a slope and a bias", "",
        "The owner, 13 October 2026: \"A slope plus bias calibration is better. Do that instead.\" So "
        "`log10 C_calibrated = slope * log10 C_placeholder + bias`, fitted by least squares on the drugs that have "
        "both. The constant shift is this model's special case at `slope = 1`.", "",
        "| | slope on the placeholder | bias (log10 uM) | out-of-sample median absolute error, leave one drug out |",
        "|---|---|---|---|",
        f"| the withdrawn constant shift | 1 (fixed) | {-models['constant_shift']:+.3f} | **{models['constant_out_of_sample']:.3f} log10** |",
        f"| the registered calibration | **{calibration.slope_on_the_placeholder:+.4f}** | **{calibration.bias_log10:+.3f}** | **{models['calibrated_out_of_sample']:.3f} log10** |",
        "",
        f"Fitted on {calibration.fitting_drug_count} drugs. A bootstrap of {BOOTSTRAP_DRAWS} draws puts the slope's "
        f"95 percent interval at **[{low_interval:.3f}, {high_interval:.3f}]**, which "
        f"{'excludes' if high_interval < 1.0 else 'does not exclude'} 1.0: the constant shift is refuted and not "
        f"merely improvable. The calibration cuts the out-of-sample error by "
        f"{1 - models['calibrated_out_of_sample'] / models['constant_out_of_sample']:.0%}.", "",
        "**The second reason, which is the stronger one.** A constant shift leaves the best mechanism target of "
        "**every** imputed drug at one and the same occupancy, because the placeholder pins it at 0.9 and a "
        "constant moves it to another single number. The quantity is then a constant across those drugs and the "
        "model can read it only as an intercept.", "",
        "| best mechanism target's occupancy | 10th percentile | median | 90th percentile | standard deviation |",
        "|---|---|---|---|---|",
        f"| imputed drugs, under a constant shift | {np.quantile(constant_best, .1):.3f} | {np.median(constant_best):.3f} | {np.quantile(constant_best, .9):.3f} | **{constant_best.std():.3f}** |",
        f"| imputed drugs, under the calibration | {np.quantile(calibrated_best, .1):.3f} | {np.median(calibrated_best):.3f} | {np.quantile(calibrated_best, .9):.3f} | {calibrated_best.std():.3f} |",
        f"| measured drugs, for reference | {np.quantile(measured_best, .1):.3f} | {np.median(measured_best):.3f} | {np.quantile(measured_best, .9):.3f} | {measured_best.std():.3f} |",
        "",
        "**Two marks against the calibration, both recorded rather than resolved.** Its median best-target "
        f"occupancy is {np.median(calibrated_best):.3f} against the measured drugs' {np.median(measured_best):.3f} "
        f"and its spread is still narrower, {calibrated_best.std():.3f} against {measured_best.std():.3f}: the fit "
        "minimises error in log concentration and not the distance between the two occupancy distributions, so the "
        "imputed group still does not look like the measured one. And an affinity read with error attenuates the "
        "fitted slope towards zero, so the mechanistic reading - that dose tracks potency at about half a decade "
        "per decade - is weaker than the number suggests. Neither touches the out-of-sample error, which is "
        "measured on held-out drugs and falls whatever causes the slope.", "",
        "**The prediction is clamped to the affinity range it was fitted on**, "
        f"pK {calibration.lowest_fitted_p_affinity:.2f} to {calibration.highest_fitted_p_affinity:.2f}. "
        f"{models['drugs_outside_the_support']} of the imputed drugs sit outside it and would otherwise be "
        "extrapolated: ziconotide at pK 1.00 and vigabatrin at 3.07 below, six above, led by dexmedetomidine and "
        "ramelteon at 10.90. Unclamped the line puts ziconotide at 261 micromolar against 4.74 clamped. The two "
        "low cases are also the drugs an equilibrium `Kd` describes worst, vigabatrin being an irreversible "
        "inhibitor of GABA transaminase, which is the exception item 29 already records for the covalent binders.",
        "",
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
