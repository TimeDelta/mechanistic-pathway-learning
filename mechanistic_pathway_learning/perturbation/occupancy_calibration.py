"""Calibrating an imputed drug concentration against the measured ones: a slope and a bias in log10.

Preregistration items 5c and 29. A drug's target occupancy is `theta = C / (C + Kd)`, so the concentration `C` has
to come from somewhere. For 115 of the study's 183 drugs with a mechanism affinity it is measured. For the other 68
it is imputed from the placeholder scheme, which fixes the best-bound mechanism target at 0.9 occupancy and so
asserts `C = 9 * Kd` of that target.

**Why a constant shift was rejected** (the owner, 13 October 2026: "A slope plus bias calibration is better. Do
that instead"). Two measurements, both in `docs/occupancy_coverage_and_shift_transfer.md`:

  1. The placeholder's log concentration falls with the best mechanism affinity at a slope of exactly -1, because
     that is what `C = 9 * Kd` means. The measured concentration falls at about -0.53: a more potent drug is dosed
     lower, but only around half as much lower. A bootstrap over the 115 fitting drugs puts the slope on the
     placeholder at [0.42, 0.64], which excludes 1.0, so a shift that keeps the slope at 1 is refuted and not
     merely improvable.
  2. A constant shift leaves the best mechanism target of every imputed drug at one and the same occupancy, 0.838,
     with a standard deviation of zero across all 68 of them. The quantity carries no information for those drugs:
     it is a constant the model can only read as an intercept. The measured drugs' best targets spread over 0.135
     to 0.988. A slope restores a spread, 0.451 to 0.982.

So the calibration is `log10 C_calibrated = slope * log10 C_placeholder + bias`, fitted by least squares on the
drugs that have both. That is a two-parameter linear calibration in log space - a gain and an offset - and the
constant shift is its special case at `slope = 1`. Written against the placeholder rather than against the affinity
because the placeholder is the quantity being corrected, though the two forms are the same model: the placeholder's
log concentration is `log10 9 + 6 - pK`, so a line in one is a line in the other.

**What it is applied to, and what it is never applied to.** Only a drug with no measured concentration. A measured
value is never shifted (item 5c step 2). The imputed values of the drugs that have a measurement exist to fit this
calibration and enter no model input (item 5c step 3).

**The prediction is clamped to the affinity range it was fitted on.** A calibration line read outside its support
is extrapolation, and 8 of the 68 imputed drugs sit outside it: ziconotide at pK 1.00 and vigabatrin at 3.07 below,
and six above 10.40, led by dexmedetomidine and ramelteon at 10.90. Unclamped the line puts ziconotide at 261
micromolar; clamped it puts it at 4.74. The two low cases are also the drugs whose affinity an equilibrium Kd
describes worst - vigabatrin inhibits GABA transaminase irreversibly - which is the exception item 29 already
records for the covalent binders.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

DEFAULT_CALIBRATION_FILE = Path("configs/occupancy_calibration.json")
PLACEHOLDER_OCCUPANCY_OF_THE_BEST_MECHANISM_TARGET = 0.9


@dataclass(frozen=True)
class ConcentrationCalibration:
    """`log10 C_calibrated = slope * log10 C_placeholder + bias`, with the affinity support it was fitted on.

    slope_on_the_placeholder: 1.0 reproduces the withdrawn constant-shift scheme.
    bias_log10: the offset, in log10 micromolar.
    lowest_fitted_p_affinity / highest_fitted_p_affinity: the range of best mechanism affinities among the fitting
      drugs. A prediction outside it is clamped to the boundary rather than extrapolated.
    fitting_drug_count: how many drugs carried both a placeholder and a measured concentration.
    """

    slope_on_the_placeholder: float
    bias_log10: float
    lowest_fitted_p_affinity: float
    highest_fitted_p_affinity: float
    fitting_drug_count: int

    def write(self, path: Path = DEFAULT_CALIBRATION_FILE) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(asdict(self), indent=2) + "\n")

    @classmethod
    def read(cls, path: Path = DEFAULT_CALIBRATION_FILE) -> "ConcentrationCalibration":
        return cls(**json.loads(Path(path).read_text()))


def placeholder_log10_concentration(best_mechanism_p_affinity: float) -> float:
    """The placeholder scheme as a log10 concentration in micromolar.

    `theta = 0.9` means `C = 9 * Kd`, and `Kd` in micromolar is `10 ** (6 - pK)`."""
    occupancy = PLACEHOLDER_OCCUPANCY_OF_THE_BEST_MECHANISM_TARGET
    return float(np.log10(occupancy / (1 - occupancy)) + 6.0 - best_mechanism_p_affinity)


def fit_concentration_calibration(best_mechanism_p_affinity: np.ndarray,
                                  measured_log10_concentration: np.ndarray) -> ConcentrationCalibration:
    """Least squares of the measured log10 concentration on the placeholder's, over the drugs that have both.

    One row per drug, never one per target: the concentration is a property of the drug and a drug with many
    targets would otherwise set the fit (decision 9.3 of docs/open_decisions.md)."""
    affinities = np.asarray(best_mechanism_p_affinity, dtype=float)
    measured = np.asarray(measured_log10_concentration, dtype=float)
    if affinities.shape != measured.shape:
        raise ValueError(f"one measured concentration per drug is required, got {affinities.shape} and {measured.shape}")
    usable = np.isfinite(affinities) & np.isfinite(measured)
    if usable.sum() < 3:
        raise ValueError(f"a slope and a bias need at least three drugs with both, got {int(usable.sum())}")
    placeholder = np.array([placeholder_log10_concentration(value) for value in affinities[usable]])
    slope, bias = np.polyfit(placeholder, measured[usable], 1)
    return ConcentrationCalibration(slope_on_the_placeholder=float(slope), bias_log10=float(bias),
                                    lowest_fitted_p_affinity=float(affinities[usable].min()),
                                    highest_fitted_p_affinity=float(affinities[usable].max()),
                                    fitting_drug_count=int(usable.sum()))


def calibrated_concentration_micromolar(best_mechanism_p_affinity: float,
                                        calibration: ConcentrationCalibration) -> float:
    """The imputed concentration for a drug with no measured one, clamped to the fitted affinity support."""
    clamped_affinity = float(np.clip(best_mechanism_p_affinity,
                                     calibration.lowest_fitted_p_affinity,
                                     calibration.highest_fitted_p_affinity))
    log10_placeholder = placeholder_log10_concentration(clamped_affinity)
    return float(10.0 ** (calibration.slope_on_the_placeholder * log10_placeholder + calibration.bias_log10))


def concentration_for_drug(best_mechanism_p_affinity: float,
                           measured_concentration_micromolar: float | None,
                           calibration: ConcentrationCalibration) -> tuple[float, bool]:
    """Item 5c steps 2 and 4: a measurement always wins, and the calibration reaches only the drugs without one.

    Returns the concentration and whether it was measured, the second being the `has_measured_occupancy` flag that
    item 29 reports every result within."""
    if measured_concentration_micromolar is not None and np.isfinite(measured_concentration_micromolar) \
            and measured_concentration_micromolar > 0:
        return float(measured_concentration_micromolar), True
    return calibrated_concentration_micromolar(best_mechanism_p_affinity, calibration), False


def occupancy(concentration_micromolar: float, p_affinity: float) -> float:
    """`theta = C / (C + Kd)` with `Kd` from an affinity given as -log10 of a molar constant."""
    dissociation_constant_micromolar = 10.0 ** (6.0 - p_affinity)
    return float(concentration_micromolar / (concentration_micromolar + dissociation_constant_micromolar))
