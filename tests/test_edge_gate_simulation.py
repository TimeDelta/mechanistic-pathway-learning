"""The planted-truth simulation of the hard edge gates (experiments/simulate_edge_gate_recovery.py, decision 9.62), its
digest (experiments/summarise_edge_gate_simulations.py), the runner that reproduces it (scripts/run_edge_gate_simulations.sh)
and the trainer's agreement with the settings it tested (docs/edge_gate_simulation.md, section 5)."""
from __future__ import annotations

import json
import math
import os
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

from experiments.run_main_model import build_argument_parser as build_trainer_parser
from experiments.simulate_edge_gate_recovery import TIERS, build_argument_parser, make_world
from experiments.summarise_edge_gate_simulations import OUTPUT_PATH, RUNS, posterior_false, selected_scores

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DIGEST = REPOSITORY_ROOT / OUTPUT_PATH


def small_world(**changes):
    settings = dict(num_nodes=120, mean_in_degree=6.0, num_perturbations=60, num_symptoms=4, nodes_per_symptom=6, seed=3)
    settings.update(changes)
    return make_world(**settings)


def tier_masks(world):
    return {name: world.prior_p == p for name, p, _, _ in TIERS}


def test_a_world_is_the_same_world_from_the_same_seed():
    first, second = small_world(), small_world()
    assert np.array_equal(first.edges, second.edges)
    assert np.array_equal(first.real, second.real)
    assert torch.equal(first.labels, second.labels)
    other = small_world(seed=4)
    assert not (np.array_equal(first.edges, other.edges) and np.array_equal(first.real, other.real)), "another seed, another world"


def test_tiers_come_in_the_v5_proportions_and_the_truth_is_calibrated_to_the_priors():
    world = make_world(1000, 6.0, 20, 2, 6, seed=0)
    masks = tier_masks(world)
    for name, p, n, share in TIERS:
        prior = 1.0 / (1.0 + math.exp(-n * math.log(p / (1.0 - p))))
        assert abs(masks[name].mean() - share) < 0.02, name
        assert np.allclose(world.prior[masks[name]], prior), "pi = sigmoid(n logit p)"
        assert abs(world.real[masks[name]].mean() - prior) < 0.05, f"{name}: each edge is real with probability pi"


def test_a_miscalibrated_tier_changes_only_that_tier_of_the_truth():
    calibrated = small_world()
    miscalibrated = small_world(truth_real_probability={"text_mined_only": 0.1})
    masks = tier_masks(calibrated)
    assert np.array_equal(calibrated.edges, miscalibrated.edges)
    assert np.array_equal(calibrated.prior, miscalibrated.prior), "the fit still reads pi"
    assert np.array_equal(calibrated.real[~masks["text_mined_only"]], miscalibrated.real[~masks["text_mined_only"]])
    assert miscalibrated.real[masks["text_mined_only"]].mean() < calibrated.real[masks["text_mined_only"]].mean()
    with pytest.raises(ValueError):
        small_world(truth_real_probability={"no_such_tier": 0.5})


def test_label_sensitivity_hides_true_positives_and_adds_none():
    full = small_world(num_perturbations=200)
    half = small_world(num_perturbations=200, label_sensitivity=0.5)
    assert torch.equal(full.true_logits, half.true_logits)
    assert bool((half.labels <= full.labels).all()), "an observed positive is a true positive"
    kept = float(half.labels.sum() / full.labels.sum())
    assert 0.4 < kept < 0.6


def test_the_readout_is_scaled_to_logits_with_a_standard_deviation_of_two():
    world = small_world()
    assert float((world.true_logits - world.bias).std()) == pytest.approx(2.0, rel=1e-4)
    assert torch.all(world.bias == -1.0)


def test_the_posterior_that_a_closed_edge_is_false():
    for prior in (0.525, 0.604, 0.9):
        logit_prior = math.log(prior / (1.0 - prior))
        assert posterior_false(logit_prior, prior) == pytest.approx(0.5), "at the closing threshold E = logit pi the posterior is one half"
        assert posterior_false(logit_prior + 2.0, prior) == pytest.approx(1.0 / (1.0 + math.exp(-2.0)))
    assert posterior_false(1.0, 0.9) < posterior_false(1.0, 0.6), "the same evidence weighs less against a stronger prior"


def test_selected_scores_count_the_false_edges_among_those_kept():
    edges = [{"real": False, "evidence": 2.0}, {"real": True, "evidence": 1.5}, {"real": False, "evidence": 0.2}]
    assert selected_scores(edges, lambda edge: edge["evidence"] >= 1.0) == {"closed": 2, "false": 1, "precision": 0.5}
    assert selected_scores(edges, lambda edge: edge["evidence"] >= 5.0) == {"closed": 0, "false": 0, "precision": None}


def test_a_quick_run_with_the_curated_tier_taken_as_known_closes_no_curated_edge(tmp_path):
    output = tmp_path / "quick.json"
    subprocess.run([sys.executable, "experiments/simulate_edge_gate_recovery.py", "--quick", "--ungated-tiers", "curated", "--closure-evidence",
                    "--per-edge-detail", "--output", str(output)], cwd=REPOSITORY_ROOT, check=True, capture_output=True,
                   env={**os.environ, "PYTHONPATH": str(REPOSITORY_ROOT)})
    run = json.loads(output.read_text())
    assert run["settings"]["ungated_tiers"] == ["curated"]
    assert run["fit_by_tier"]["curated"]["closed"] == 0
    assert all(edge["tier"] != "curated" for edge in run["closed_edges"])
    for edge in run["closed_edges"]:
        if edge["evidence_nats"] is not None:
            assert 0.0 <= posterior_false(edge["evidence_nats"], edge["prior"]) <= 1.0
    loss = run["held_out_loss"]
    assert set(loss) == {"fitted_graph", "every_edge_open", "true_graph"}


def runner_arguments() -> dict[str, str]:
    listed = subprocess.run(["bash", "scripts/run_edge_gate_simulations.sh", "--list"], cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True).stdout
    return dict(line.split("\t", 1) for line in listed.splitlines() if line.strip())


def test_the_runner_defines_every_run_the_digest_reads_with_the_settings_it_recorded():
    digest = json.loads(DIGEST.read_text())
    arguments = runner_arguments()
    assert set(RUNS) <= set(arguments), f"runs the digest reads but the runner does not define: {sorted(set(RUNS) - set(arguments))}"
    parser = build_argument_parser()
    for name, run in digest["runs"].items():
        parsed = vars(parser.parse_args(shlex.split(arguments[name])))
        for key, recorded in run["settings"].items():
            assert parsed[key] == recorded, f"{name}: the runner gives {key} = {parsed[key]!r}, the run recorded {recorded!r}"


def test_the_trainer_takes_the_settings_the_simulations_tested():
    tested = json.loads(DIGEST.read_text())["runs"]["b3_s0"]["settings"]
    trainer = build_trainer_parser().parse_args([])
    assert trainer.edge_gate_binary_epochs == tested["binary_epochs"]
    assert trainer.edge_gate_binary_learning_rate == tested["binary_learning_rate"]
    assert trainer.edge_gate_binary_learning_rate_decay == tested.get("binary_learning_rate_decay", 1.0), "a run before the option kept the rate"
    assert trainer.edge_gate_surrogate_width == tested["surrogate_width"]
    assert trainer.edge_gate_binary_full_batch == tested["binary_full_batch"]
    assert trainer.edge_gate_beta == pytest.approx(tested["beta"])
    assert trainer.edge_gate_kind == "hard_concrete"
