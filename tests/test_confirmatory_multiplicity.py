"""Amendment item 18: H1 first, then H2 and P1 together by Holm at one-sided 0.025 per model."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from score_confirmatory import confirmatory_decisions, holm_rejections

REGISTERED_ALPHA = 0.025
HOLM_LEVEL_OF_THE_SMALLER_OF_TWO = REGISTERED_ALPHA / 2


def test_holm_tests_the_smaller_p_value_at_half_the_level_and_the_larger_at_the_full_level():
    rejections = holm_rejections({"H2": 0.010, "P1": 0.020}, REGISTERED_ALPHA)
    assert rejections["H2"]["level"] == pytest.approx(HOLM_LEVEL_OF_THE_SMALLER_OF_TWO)
    assert rejections["P1"]["level"] == pytest.approx(REGISTERED_ALPHA)
    assert rejections["H2"]["rejected"] and rejections["P1"]["rejected"]


def test_holm_stops_at_the_first_failure_so_a_larger_p_value_below_alpha_is_not_rejected():
    rejections = holm_rejections({"H2": 0.020, "P1": 0.024}, REGISTERED_ALPHA)
    assert not rejections["H2"]["rejected"] and not rejections["P1"]["rejected"]


def test_h2_needs_half_the_level_unless_p1_clears_it_first():
    """The sentence item 18 states: 'H2 therefore needs p below 0.0125 unless P1 also passes'."""
    without_p1_passing = confirmatory_decisions(True, p_h1=0.001, p_h2=0.020, meets_minimum_differences=True,
                                                alpha=REGISTERED_ALPHA, p_p1=0.900, pathway_gates_passed=True)
    assert not without_p1_passing["H2_confirmed"]
    with_p1_passing = confirmatory_decisions(True, p_h1=0.001, p_h2=0.020, meets_minimum_differences=True,
                                             alpha=REGISTERED_ALPHA, p_p1=0.005, pathway_gates_passed=True)
    assert with_p1_passing["H2_confirmed"] and with_p1_passing["P1_confirmed"]
    assert with_p1_passing["second_stage_holm"]["H2"]["level"] == pytest.approx(REGISTERED_ALPHA)


def test_no_pathway_test_leaves_h2_alone_at_the_full_level():
    decisions = confirmatory_decisions(True, p_h1=0.001, p_h2=0.020, meets_minimum_differences=True, alpha=REGISTERED_ALPHA)
    assert decisions["second_stage_family"] == ["H2"]
    assert decisions["second_stage_holm"]["H2"]["level"] == pytest.approx(REGISTERED_ALPHA)
    assert decisions["H2_confirmed"] and "P1_confirmed" not in decisions


def test_a_failed_gate_blocks_the_p1_claim_without_changing_the_level_h2_faced():
    gates_passed = confirmatory_decisions(True, p_h1=0.001, p_h2=0.020, meets_minimum_differences=True,
                                          alpha=REGISTERED_ALPHA, p_p1=0.005, pathway_gates_passed=True)
    gates_failed = confirmatory_decisions(True, p_h1=0.001, p_h2=0.020, meets_minimum_differences=True,
                                          alpha=REGISTERED_ALPHA, p_p1=0.005, pathway_gates_passed=False)
    assert gates_passed["P1_confirmed"] and not gates_failed["P1_confirmed"]
    assert gates_failed["H2_confirmed"]
    assert gates_failed["second_stage_holm"]["H2"]["level"] == gates_passed["second_stage_holm"]["H2"]["level"]


def test_the_second_stage_is_not_entered_until_h1_is_confirmed():
    for blocked in (dict(complete=False, meets_minimum_differences=True), dict(complete=True, meets_minimum_differences=False)):
        decisions = confirmatory_decisions(p_h1=0.001, p_h2=0.001, alpha=REGISTERED_ALPHA, p_p1=0.001,
                                           pathway_gates_passed=True, **blocked)
        assert not decisions["H1_confirmed"] and not decisions["H2_confirmed"] and not decisions["P1_confirmed"]


def test_a_p_value_at_the_level_is_rejected_and_one_above_it_is_not():
    at_the_level = holm_rejections({"H2": HOLM_LEVEL_OF_THE_SMALLER_OF_TWO, "P1": REGISTERED_ALPHA}, REGISTERED_ALPHA)
    assert at_the_level["H2"]["rejected"] and at_the_level["P1"]["rejected"]
    above = holm_rejections({"H2": HOLM_LEVEL_OF_THE_SMALLER_OF_TWO + 1e-9, "P1": REGISTERED_ALPHA}, REGISTERED_ALPHA)
    assert not above["H2"]["rejected"] and not above["P1"]["rejected"]
