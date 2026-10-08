"""--start-at-weighted-optimum (experiments/run_main_model.py): the constant it starts from minimises the training loss,
it sets the noisy-OR leak or the sigmoid output bias, and left off it changes neither the model nor the fingerprint."""
from argparse import Namespace

import numpy as np
import torch

from experiments.run_main_model import configuration_fingerprint, initialise_rates_at_weighted_optimum, weighted_constant_optimum
from mechanistic_pathway_learning.models.baselines.relational_gnn_sigmoid_baseline import RelationalGnnSigmoidHead
from mechanistic_pathway_learning.models.noisy_or_pathway_module_model import NoisyOrPathwayModuleHead
from mechanistic_pathway_learning.models.soft_constraint_losses import evidence_weighted_binary_cross_entropy


def small_data() -> tuple[Namespace, np.ndarray]:
    random_numbers = np.random.default_rng(0)
    outcomes = (random_numbers.random((40, 3)) < np.array([0.05, 0.2, 0.5])).astype(float)
    weights = np.where(outcomes > 0, random_numbers.uniform(0.3, 1.0, outcomes.shape), 0.0)
    frequencies = np.where(outcomes > 0, random_numbers.uniform(0.0, 1.0, outcomes.shape), np.nan)
    label_mask = random_numbers.random(outcomes.shape) > 0.1
    return Namespace(outcomes=outcomes, weights=weights, frequencies=frequencies), label_mask


def arguments_for(head: str, start: bool = True, from_frequency: bool = False) -> Namespace:
    return Namespace(head=head, start_at_weighted_optimum=start, positive_target=0.99, negative_weight=0.2,
                     positive_target_from_frequency=from_frequency, minimum_frequency_target=0.05)


def trainer_loss(data, rows, label_mask, constant: np.ndarray) -> float:
    """The trainer's loss (default targets and weights) for a constant prediction per symptom."""
    weights = np.where(data.outcomes[rows] > 0, np.maximum(data.weights[rows], 1e-3), 0.2) * label_mask[rows]
    probability = torch.as_tensor(np.broadcast_to(constant, data.outcomes[rows].shape).copy(), dtype=torch.float32)
    return float(evidence_weighted_binary_cross_entropy(probability, torch.as_tensor(data.outcomes[rows], dtype=torch.float32),
                                                        torch.as_tensor(weights, dtype=torch.float32), positive_target=0.99))


def test_the_constant_minimises_the_trainer_loss_per_symptom() -> None:
    data, label_mask = small_data()
    rows = np.arange(30)
    optimum = weighted_constant_optimum(data, rows, label_mask, arguments_for("sigmoid"))
    best = trainer_loss(data, rows, label_mask, optimum)
    for symptom in range(3):
        for step in (-0.02, 0.02):
            moved = optimum.copy()
            moved[symptom] += step
            assert trainer_loss(data, rows, label_mask, moved) > best
    raw_rate = (data.outcomes[rows] * label_mask[rows]).sum(0) / label_mask[rows].sum(0)
    assert np.all(optimum > raw_rate)  # negatives weigh 0.2, so the optimum sits above the raw rate


def test_frequency_targets_give_the_weighted_mean_of_the_frequencies() -> None:
    data, label_mask = small_data()
    rows = np.arange(40)
    optimum = weighted_constant_optimum(data, rows, label_mask, arguments_for("sigmoid", from_frequency=True))
    positive = data.outcomes > 0
    targets = np.where(positive, np.maximum(np.nan_to_num(data.frequencies, nan=0.99), 0.05), 0.0)
    weights = np.where(positive, 1.0, 0.2) * label_mask
    assert np.allclose(optimum, (weights * targets).sum(0) / weights.sum(0))


def test_the_flag_sets_the_noisy_or_leak_and_the_sigmoid_bias() -> None:
    data, label_mask = small_data()
    rows = np.arange(40)
    optimum = weighted_constant_optimum(data, rows, label_mask, arguments_for("noisy_or"))
    noisy_or_head = NoisyOrPathwayModuleHead(5, 4, 2, 3)
    initialise_rates_at_weighted_optimum(noisy_or_head, data, rows, label_mask, arguments_for("noisy_or"), "cpu")
    assert np.allclose(torch.sigmoid(noisy_or_head.symptom_leak_logit[0]).detach().numpy(), optimum, atol=1e-6)
    sigmoid_head = RelationalGnnSigmoidHead(4, 3)
    initialise_rates_at_weighted_optimum(sigmoid_head, data, rows, label_mask, arguments_for("sigmoid"), "cpu")
    assert np.allclose(torch.sigmoid(sigmoid_head.readout[-1].bias).detach().numpy(), optimum, atol=1e-6)


def test_off_the_flag_changes_nothing() -> None:
    data, label_mask = small_data()
    torch.manual_seed(0)
    head = RelationalGnnSigmoidHead(4, 3)
    before = {name: value.clone() for name, value in head.state_dict().items()}
    initialise_rates_at_weighted_optimum(head, data, np.arange(40), label_mask, arguments_for("sigmoid", start=False), "cpu")
    assert all(torch.equal(before[name], value) for name, value in head.state_dict().items())


def test_off_the_flag_leaves_the_fingerprint_as_it_was() -> None:
    """A run resumed from a checkpoint written before the flag existed must report no configuration change."""
    data, label_mask = small_data()
    data.edge_source, data.edge_target, data.edge_relation = np.array([0]), np.array([1]), np.array([0])
    data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes = [np.array([0])], [np.array([-1.0])], [np.array([1.0])]
    off = configuration_fingerprint(Namespace(head="sigmoid", start_at_weighted_optimum=False), data, label_mask)
    before_the_flag = configuration_fingerprint(Namespace(head="sigmoid"), data, label_mask)
    assert off == before_the_flag
    on = configuration_fingerprint(Namespace(head="sigmoid", start_at_weighted_optimum=True), data, label_mask)
    assert on["argument:start_at_weighted_optimum"] == "True"
