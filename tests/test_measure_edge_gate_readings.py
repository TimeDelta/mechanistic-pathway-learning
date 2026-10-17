"""The numbers decision 9.62 records: the closed forms of a hard concrete gate (Louizos et al. 2018, eqs. 12 and 13)
against sampling, the KL identity behind an evidence-weighted L0 cost, and the three binary readings of a (p, n) prior."""
from __future__ import annotations

import numpy as np

from experiments.measure_edge_gate_readings import (
    bernoulli_kl,
    binary_priors,
    hard_concrete_masses,
    log_alpha_for_open_probability,
    logit,
    sample_hard_concrete,
)


def test_the_closed_forms_match_sampled_gates():
    generator = np.random.default_rng(0)
    for prior in (0.55, 0.7, 0.9):
        log_alpha = log_alpha_for_open_probability(prior)
        masses = hard_concrete_masses(log_alpha)
        samples = sample_hard_concrete(log_alpha, 400_000, generator)
        assert abs(float((samples > 0).mean()) - prior) < 0.004, "eq. 12: the gate is open with the prior's probability"
        assert abs(float((samples == 0).mean()) - float(masses["closed"])) < 0.004
        assert abs(float((samples == 1).mean()) - float(masses["open"])) < 0.004
        assert abs(float(masses["closed"] + masses["fractional"] + masses["open"]) - 1.0) < 1e-12


def test_the_recorded_values_at_a_curated_edges_prior():
    masses = hard_concrete_masses(log_alpha_for_open_probability(0.9))
    assert round(float(masses["open"]), 2) == 0.27 and round(float(masses["fractional"]), 2) == 0.63 and round(float(masses["test_time_gate"]), 2) == 0.67


def test_an_l0_cost_of_the_prior_log_odds_less_the_entropy_is_the_kl_to_the_prior():
    for open_probability in (0.05, 0.3, 0.8, 0.99):
        for prior in (0.2, 0.55, 0.9):
            entropy = -(open_probability * np.log(open_probability) + (1 - open_probability) * np.log(1 - open_probability))
            l0_form = open_probability * np.log((1 - prior) / prior) - entropy - np.log(1 - prior)
            assert abs(float(bernoulli_kl(open_probability, prior)) - float(l0_form)) < 1e-12


def test_the_three_binary_readings_of_a_prior():
    readings = binary_priors(np.array([0.9, 0.9, 0.7]), np.array([1.0, 2.0, 0.5]))
    assert np.allclose(readings["p"], [0.9, 0.9, 0.7])
    assert np.allclose(readings["beta_mean"], [1.9 / 3.0, 2.8 / 4.0, 1.35 / 2.5])
    assert abs(float(readings["log_odds_times_n"][0]) - 0.9) < 1e-12, "n = 1 gives p"
    assert float(readings["log_odds_times_n"][1]) > 0.9 > float(readings["beta_mean"][1]), "more evidence raises the prior under the log-odds reading only"
    assert abs(float(logit(readings["log_odds_times_n"][2])) - 0.5 * float(logit(0.7))) < 1e-12
