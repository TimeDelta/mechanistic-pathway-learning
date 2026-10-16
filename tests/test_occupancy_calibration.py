"""Preregistration item 5c(ii): the imputed concentration is calibrated by a slope and a bias, clamped to support.

What is held here is the form and its guards, not the fitted numbers: those refit on the table the lockbox is drawn
from (decision 9.2), so a test pinning 0.5322 would pin the fit to the table in use.
"""
from __future__ import annotations


import numpy as np
import pytest

from mechanistic_pathway_learning.perturbation.occupancy_calibration import (
    DEFAULT_CALIBRATION_FILE,
    PLACEHOLDER_OCCUPANCY_OF_THE_BEST_MECHANISM_TARGET,
    ConcentrationCalibration,
    calibrated_concentration_micromolar,
    concentration_for_drug,
    fit_concentration_calibration,
    occupancy,
    placeholder_log10_concentration,
)


def calibration_at(slope: float, bias: float, lowest: float = 4.0, highest: float = 10.0) -> ConcentrationCalibration:
    return ConcentrationCalibration(slope_on_the_placeholder=slope, bias_log10=bias,
                                    lowest_fitted_p_affinity=lowest, highest_fitted_p_affinity=highest,
                                    fitting_drug_count=115)


def test_the_placeholder_puts_the_best_mechanism_target_at_nine_tenths_occupancy():
    for p_affinity in (5.0, 7.5, 9.0, 10.0):
        concentration = 10.0 ** placeholder_log10_concentration(p_affinity)
        assert occupancy(concentration, p_affinity) == pytest.approx(PLACEHOLDER_OCCUPANCY_OF_THE_BEST_MECHANISM_TARGET)


def test_slope_one_reproduces_the_withdrawn_constant_shift():
    """The constant shift is this model's special case, which is why the slope's interval is the test of it."""
    withdrawn_shift_log10 = 0.242
    calibrated = calibrated_concentration_micromolar(8.0, calibration_at(1.0, -withdrawn_shift_log10))
    assert calibrated == pytest.approx(10.0 ** placeholder_log10_concentration(8.0) / 10 ** withdrawn_shift_log10)


def test_a_constant_shift_leaves_the_best_target_at_one_occupancy_and_a_slope_does_not():
    """The second and stronger reason item 5c(ii) gives for the slope: at slope 1 the quantity is a constant."""
    affinities = [5.0, 7.0, 8.5, 10.0]
    under_constant = [occupancy(calibrated_concentration_micromolar(value, calibration_at(1.0, -0.242)), value)
                      for value in affinities]
    under_calibration = [occupancy(calibrated_concentration_micromolar(value, calibration_at(0.5322, -0.753)), value)
                         for value in affinities]
    assert np.std(under_constant) == pytest.approx(0.0, abs=1e-12)
    assert np.std(under_calibration) > 0.05
    assert sorted(under_calibration) == under_calibration, "a more potent drug should not end up less occupied"


def test_the_prediction_is_clamped_to_the_fitted_affinity_support():
    """Eight of the 68 imputed drugs sit outside it, ziconotide at pK 1.00 being the extreme."""
    calibration = calibration_at(0.5322, -0.753, lowest=4.27, highest=10.40)
    below = calibrated_concentration_micromolar(1.00, calibration)
    at_the_floor = calibrated_concentration_micromolar(4.27, calibration)
    above = calibrated_concentration_micromolar(12.0, calibration)
    at_the_ceiling = calibrated_concentration_micromolar(10.40, calibration)
    assert below == pytest.approx(at_the_floor)
    assert above == pytest.approx(at_the_ceiling)
    unclamped_at_one = 10 ** (0.5322 * placeholder_log10_concentration(1.00) - 0.753)
    assert unclamped_at_one > 50 * below, "the clamp has to be doing real work at pK 1.00"


def test_a_measurement_always_wins_and_is_never_calibrated():
    """Item 5c step 2. The flag returned is what item 29 reports every result within."""
    calibration = calibration_at(0.5322, -0.753)
    measured, was_measured = concentration_for_drug(8.0, 0.025, calibration)
    assert (measured, was_measured) == (0.025, True)
    imputed, was_measured = concentration_for_drug(8.0, None, calibration)
    assert was_measured is False
    assert imputed == pytest.approx(calibrated_concentration_micromolar(8.0, calibration))
    for absent in (None, float("nan"), 0.0, -1.0):
        assert concentration_for_drug(8.0, absent, calibration)[1] is False


def test_the_fit_recovers_a_planted_slope_and_bias():
    affinities = np.linspace(4.5, 10.0, 60)
    placeholder = np.array([placeholder_log10_concentration(value) for value in affinities])
    planted_slope, planted_bias = 0.53, -0.75
    fitted = fit_concentration_calibration(affinities, planted_slope * placeholder + planted_bias)
    assert fitted.slope_on_the_placeholder == pytest.approx(planted_slope)
    assert fitted.bias_log10 == pytest.approx(planted_bias)
    assert (fitted.lowest_fitted_p_affinity, fitted.highest_fitted_p_affinity) == (4.5, 10.0)
    assert fitted.fitting_drug_count == 60


def test_the_fit_refuses_mismatched_or_too_few_drugs():
    with pytest.raises(ValueError, match="one measured concentration per drug"):
        fit_concentration_calibration(np.array([8.0, 9.0]), np.array([-1.0]))
    with pytest.raises(ValueError, match="at least three drugs"):
        fit_concentration_calibration(np.array([8.0, 9.0]), np.array([-1.0, -1.2]))
    with pytest.raises(ValueError, match="at least three drugs"):
        fit_concentration_calibration(np.array([8.0, 9.0, np.nan]), np.array([-1.0, np.nan, -1.2]))


def test_the_registered_file_round_trips_and_carries_its_support():
    stored = ConcentrationCalibration.read(DEFAULT_CALIBRATION_FILE)
    assert stored.lowest_fitted_p_affinity < stored.highest_fitted_p_affinity
    assert stored.fitting_drug_count > 3
    assert 0.0 < stored.slope_on_the_placeholder < 1.0, "a fitted slope at or above 1 would reopen decision 9.2"
