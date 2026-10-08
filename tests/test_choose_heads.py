"""Tests for the head-choice rule (experiments/choose_heads.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

from choose_heads import choose, head_margin  # noqa: E402


def test_margin_is_the_weakest_h1_condition_above_its_threshold() -> None:
    means = {"macro": 0.06, "micro": 0.05, "macro_within_degree_strata": 0.03, "micro_within_degree_strata": 0.004}
    assert abs(head_margin(means, 0.041, 0.028) - 0.004) < 1e-12
    means["micro_within_degree_strata"] = 0.1
    assert abs(head_margin(means, 0.041, 0.028) - 0.019) < 1e-12  # macro 0.06 - 0.041


def test_the_larger_margin_wins_and_a_tie_goes_to_noisy_or() -> None:
    assert choose(0.010, 0.020) == "sigmoid"
    assert choose(0.020, 0.010) == "noisy_or"
    assert choose(0.0100, 0.0105) == "noisy_or"  # within 0.001
    assert choose(-0.05, -0.02) == "sigmoid"  # both below their thresholds: still the larger margin
