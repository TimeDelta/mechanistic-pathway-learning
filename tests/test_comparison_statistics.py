"""Each comparison's own statistic (decision 9.48, option 3; evidence/comparison_statistics.py): the shrunk
observed-to-expected ratio of FAERS read from its Gamma posterior, the log risk ratio of a trial arm with the Katz
interval, both in the one shape of STATISTIC_COLUMNS."""
from __future__ import annotations

import math

import pytest
from scipy import special, stats

from mechanistic_pathway_learning.evidence.comparison_statistics import (
    CENSORED_UPPER_BOUND,
    STATISTIC_COLUMNS,
    empty_statistic,
    expected_count,
    log_risk_ratio,
    shrunk_observed_to_expected,
)


def test_the_expected_count_is_drug_reports_times_reaction_reports_over_all_reports():
    assert expected_count(1000, 20000, 20_000_000) == 1.0
    with pytest.raises(ValueError):
        expected_count(1, 1, 0)


def test_the_shrunk_ratio_is_the_gamma_posterior_of_noren_and_colleagues():
    row = shrunk_observed_to_expected(3, 1.0)
    assert set(row) == set(STATISTIC_COLUMNS)
    assert abs(row["statistic_estimate"] - math.log(3.5 / 1.5)) < 1e-12, "log of (O + 1/2) / (E + 1/2)"
    lower, upper = stats.gamma.ppf([0.025, 0.975], a=3.5, scale=1 / 1.5)
    assert abs(row["statistic_lower_95"] - math.log(lower)) < 1e-12 and abs(row["statistic_upper_95"] - math.log(upper)) < 1e-12
    assert abs(row["statistic_probability_of_excess"] - stats.gamma.sf(1.0, a=3.5, scale=1 / 1.5)) < 1e-12
    assert abs(row["statistic_information"] - 1 / special.polygamma(1, 3.5)) < 1e-12
    assert row["statistic_censored"] == ""


def test_no_report_against_none_expected_is_no_disproportion_and_more_reports_say_more():
    nothing = shrunk_observed_to_expected(0, 0.0)
    assert nothing["statistic_estimate"] == 0.0, "the shrinkage pulls an empty comparison to a ratio of 1"
    many = shrunk_observed_to_expected(200, 100.0)
    few = shrunk_observed_to_expected(2, 1.0)
    assert many["statistic_probability_of_excess"] > 0.999 and 0.5 < few["statistic_probability_of_excess"] < 0.95
    assert many["statistic_information"] > 50 * few["statistic_information"] / 2.5, "information grows with the count, about O + 1/2"
    assert abs(many["statistic_information"] - 200.0) < 1.0
    under = shrunk_observed_to_expected(50, 100.0)
    assert under["statistic_upper_95"] < 0 and under["statistic_probability_of_excess"] < 0.001, "half the expected reports is an under-reporting"


def test_a_count_known_only_as_a_bound_is_marked_censored():
    assert shrunk_observed_to_expected(40, 80.0, censored=True)["statistic_censored"] == CENSORED_UPPER_BOUND
    with pytest.raises(ValueError):
        shrunk_observed_to_expected(-1, 1.0)


def test_the_trial_statistic_is_the_log_risk_ratio_with_the_katz_interval():
    row = log_risk_ratio(10, 100, 5, 100)
    standard_error = math.sqrt(1 / 10 - 1 / 100 + 1 / 5 - 1 / 100)
    assert abs(row["statistic_estimate"] - math.log(2.0)) < 1e-12
    assert abs(row["statistic_upper_95"] - row["statistic_estimate"] - 1.959963984540054 * standard_error) < 1e-12
    assert abs(row["statistic_probability_of_excess"] - stats.norm.cdf(math.log(2.0) / standard_error)) < 1e-12
    assert abs(row["statistic_information"] - 1 / standard_error ** 2) < 1e-12


def test_a_zero_cell_takes_the_haldane_anscombe_half_and_impossible_arms_give_nothing():
    row = log_risk_ratio(0, 100, 0, 100)
    assert row["statistic_estimate"] == 0.0 and abs(row["statistic_probability_of_excess"] - 0.5) < 1e-12
    assert abs(row["statistic_information"] - 1 / (2 / 0.5 - 2 / 101)) < 1e-12
    assert log_risk_ratio(5, 0, 1, 10) is None and log_risk_ratio(11, 10, 1, 10) is None
    assert log_risk_ratio(10, 10, 10, 10) is not None, "a zero cell (no unaffected) is corrected, not refused"


def test_a_row_without_a_statistic_has_every_column_and_no_values():
    row = empty_statistic()
    assert set(row) == set(STATISTIC_COLUMNS) and row["statistic_name"] == "" and math.isnan(row["statistic_estimate"])
