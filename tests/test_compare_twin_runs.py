import numpy as np
import pytest
from scipy import stats

from experiments.compare_twin_runs import argument_pairs, fold_t_interval, twin_pairs


def test_argument_pairs_reads_flags_with_and_without_values():
    assert argument_pairs(["--head", "sigmoid", "--cofactor-relations", "--graph-dir", "data/x"]) == {
        ("--head", "sigmoid"), ("--cofactor-relations", None), ("--graph-dir", "data/x")}


def test_twin_pairs_keeps_one_flag_differences_and_puts_the_change_on_a():
    configurations = {"base": ["--a", "1"], "plus": ["--a", "1", "--flag"], "other": ["--a", "2"], "two": ["--a", "2", "--flag"]}
    found = {(run, reference): change for run, reference, change in twin_pairs(configurations)}
    assert found == {("plus", "base"): "--flag", ("other", "base"): "--a 2", ("two", "other"): "--flag", ("two", "plus"): "--a 2"}


def test_fold_t_interval_uses_n_minus_one_degrees_of_freedom():
    mean, lower, upper = fold_t_interval(np.array([1.0, 2.0, 3.0]))
    half_width = stats.t.ppf(0.975, 2) / np.sqrt(3)
    assert mean == pytest.approx(2.0)
    assert (lower, upper) == pytest.approx((2.0 - half_width, 2.0 + half_width))


def test_fold_t_interval_is_wider_than_a_five_value_percentile_bootstrap():
    differences = np.array([0.011, 0.053, -0.003, -0.032, 0.033])
    _, lower, upper = fold_t_interval(differences)
    generator = np.random.default_rng(0)
    draws = [differences[generator.integers(0, 5, 5)].mean() for _ in range(4000)]
    assert upper - lower > 1.4 * (np.percentile(draws, 97.5) - np.percentile(draws, 2.5))
