"""Tests for the post-hoc noisy-OR recomputation, the sufficiency test and the paired bootstrap comparison."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

from analyze_pathway_modules import noisy_or_probabilities, sufficiency_test  # noqa: E402

from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import paired_bootstrap_macro_difference  # noqa: E402


def test_noisy_or_recomputation_matches_closed_form_and_ablation_removes_one_module() -> None:
    activations = np.array([[1.0, 0.5], [0.0, 1.0]])
    links = np.array([[0.8, 0.0], [0.4, 0.6]])
    leaks = np.array([0.1, 0.0])
    full = noisy_or_probabilities(activations, links, leaks)
    expected_first = 1.0 - (1.0 - 0.1) * (1.0 - 0.8) * (1.0 - 0.5 * 0.4)
    assert abs(full[0, 0] - expected_first) < 1e-9
    without_module_0 = noisy_or_probabilities(activations, links, leaks, ablated_module=0)
    assert abs(without_module_0[0, 0] - (1.0 - 0.9 * (1.0 - 0.2))) < 1e-9
    assert np.allclose(without_module_0[:, 1], 1.0 - (1.0 - activations[:, 1] * 0.6))


def test_sufficiency_test_detects_the_module_that_carries_the_symptom() -> None:
    generator = np.random.default_rng(0)
    num = 60
    activations = np.zeros((num, 2))
    activations[:30, 0] = generator.uniform(0.6, 1.0, 30)  # module 0 drives the first half
    activations[30:, 1] = generator.uniform(0.6, 1.0, 30)  # module 1 drives the second half
    links = np.array([[0.9], [0.9]])
    leaks = np.array([0.02])
    outcomes = (activations.max(axis=1, keepdims=True) > 0.8).astype(float)  # the symptom follows the driving module's activation
    rows = sufficiency_test(activations, links, leaks, outcomes, ["symptom"], min_group_size=8, min_positives=3, num_bootstrap=100)
    assert {row["module"] for row in rows} == {"module_0", "module_1"}
    for row in rows:
        assert row["group_size"] == 30
        assert row["auprc_full"] > 0.99  # the full model ranks its own group perfectly
        assert row["auprc_ablate_own"] < row["auprc_full"] - 0.2  # without the dominant module only the leak remains: AUPRC falls to the base rate
        assert row["own_ablation_drop_ci"][0] > 0.0
        assert abs(row["auprc_ablate_others_mean"] - row["auprc_full"]) < 1e-9  # the other module contributes nothing to this group


def test_paired_bootstrap_difference_is_positive_for_a_better_model() -> None:
    generator = np.random.default_rng(1)
    outcomes = (generator.random((200, 3)) < 0.3).astype(float)
    better = outcomes * 0.6 + generator.random((200, 3)) * 0.4
    worse = generator.random((200, 3))
    result = paired_bootstrap_macro_difference(better, worse, outcomes, num_bootstrap=200)
    assert result["difference"] > 0.2 and result["lower"] > 0.0 and result["fraction_resamples_favoring_a"] > 0.99
    null = paired_bootstrap_macro_difference(worse, worse, outcomes, num_bootstrap=50)
    assert abs(null["difference"]) < 1e-12
