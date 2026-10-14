"""Preregistration item 5c: a measured concentration always wins, every drug is imputed to estimate the shift, and
the shift is applied only where no measurement exists.

Item 5b, which this file was first written for, was withdrawn as a misreading on the day it was registered; what
survives it is the relation in the binding isotherm the stratum rests on, which is what this file holds. Nothing
here reads a data file: the counts items 5c and 29 register are recomputed by
`experiments/label_review/occupancy_coverage_and_shift_transfer.py` against the table the lockbox is drawn from,
and a test that pinned them would pin them to the table in use rather than to the derivation. The calibration
itself - the slope, the bias and the clamp of item 5c(ii) - is held in `tests/test_occupancy_calibration.py`.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "experiments" / "label_review" / "occupancy_coverage_and_shift_transfer.py"
specification = importlib.util.spec_from_file_location("occupancy_coverage_and_shift_transfer", SCRIPT)
shift_module = importlib.util.module_from_spec(specification)
specification.loader.exec_module(shift_module)

CALIBRATION_SCATTER_LOG10 = 1.14
MEDIAN_MECHANISM_OCCUPANCY = 0.787
MEDIAN_OFF_TARGET_OCCUPANCY = 0.063


def test_the_sensitivity_is_monotone_in_the_occupancy_and_not_a_u_shape():
    """d(log theta) / d(log C) = 1 - theta, so a lower occupancy always absorbs less of a concentration error.

    The prose corrected on 13 October 2026 had this as a U shape - "an occupancy near saturation or near zero
    absorbs and one near Kd does not" - which inverts what an imputed off-target costs, since an off-target sits
    near zero. The relation is what `has_measured_occupancy` is reported against."""
    occupancies = [0.01, 0.063, 0.1, 0.3, 0.5, 0.787, 0.9, 0.99]
    downward = [shift_module.movement_under_scatter(value, CALIBRATION_SCATTER_LOG10)[1] for value in occupancies]
    upward = [shift_module.movement_under_scatter(value, CALIBRATION_SCATTER_LOG10)[0] for value in occupancies]
    assert downward == sorted(downward, reverse=True), "a lower occupancy must move further, not less far"
    assert upward == sorted(upward, reverse=True)


def test_the_asymmetry_the_stratum_reports_reproduces_from_the_two_median_occupancies():
    """x1.2 up and x3.7 down at a mechanism target, x7.6 up and x13.0 down at an off-target."""
    mechanism_up, mechanism_down = shift_module.movement_under_scatter(MEDIAN_MECHANISM_OCCUPANCY, CALIBRATION_SCATTER_LOG10)
    off_target_up, off_target_down = shift_module.movement_under_scatter(MEDIAN_OFF_TARGET_OCCUPANCY, CALIBRATION_SCATTER_LOG10)
    assert mechanism_up == pytest.approx(1.2, abs=0.05)
    assert mechanism_down == pytest.approx(3.7, abs=0.05)
    assert off_target_up == pytest.approx(7.6, abs=0.05)
    assert off_target_down == pytest.approx(13.0, abs=0.05)


def test_the_one_sd_band_at_a_median_off_target_spans_absent_to_half_engaged():
    """What an imputed off-target occupancy is worth: the band runs 0.005 to 0.48. Item 5c keeps the value and
    reports it inside `has_measured_occupancy`; the band is why the stratum is not optional."""
    up_ratio, down_ratio = shift_module.movement_under_scatter(MEDIAN_OFF_TARGET_OCCUPANCY, CALIBRATION_SCATTER_LOG10)
    assert MEDIAN_OFF_TARGET_OCCUPANCY / down_ratio == pytest.approx(0.005, abs=0.001)
    assert MEDIAN_OFF_TARGET_OCCUPANCY * up_ratio == pytest.approx(0.48, abs=0.01)
