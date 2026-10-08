"""--equal-drug-shares and --normalise-drug-input (experiments/run_main_model.py): off, the loss weights and seeds are as
before and the fingerprint is unchanged; on, every drug holds the mean counts of positives and of negatives, gene rows keep
their weights and every drug's seed magnitudes sum to 1."""
from argparse import Namespace

import numpy as np

from experiments.run_main_model import configuration_fingerprint, loss_pair_weights, weighted_constant_optimum


def small_data() -> tuple[Namespace, np.ndarray]:
    random_numbers = np.random.default_rng(1)
    outcomes = (random_numbers.random((12, 4)) < 0.4).astype(float)
    weights = np.where(outcomes > 0, random_numbers.uniform(0.3, 1.0, outcomes.shape), 0.0)
    frequencies = np.where(outcomes > 0, random_numbers.uniform(0.0, 1.0, outcomes.shape), np.nan)
    label_mask = random_numbers.random(outcomes.shape) > 0.1
    perturbation_types = ["drug"] * 8 + ["gene"] * 4
    return Namespace(outcomes=outcomes, weights=weights, frequencies=frequencies, perturbation_types=perturbation_types), label_mask


def arguments_for(equal_drug_shares: bool, from_frequency: bool = False) -> Namespace:
    return Namespace(head="sigmoid", equal_drug_shares=equal_drug_shares, positive_target=0.99, negative_weight=0.2,
                     positive_target_from_frequency=from_frequency, minimum_frequency_target=0.05)


def test_off_the_flag_gives_the_trainer_weights_as_before() -> None:
    data, label_mask = small_data()
    before_the_flag = np.where(data.outcomes > 0, np.maximum(data.weights, 1e-3), 0.2) * label_mask
    assert np.allclose(loss_pair_weights(data, label_mask, arguments_for(False), np.arange(10)), before_the_flag)
    frequency_weights = np.where(data.outcomes > 0, 1.0, 0.2) * label_mask
    assert np.allclose(loss_pair_weights(data, label_mask, arguments_for(False, True), np.arange(10), frequency_targets=True), frequency_weights)


def test_on_the_flag_every_drug_holds_the_mean_counts_of_positives_and_negatives() -> None:
    data, label_mask = small_data()
    fit_indices = np.arange(10)  # drugs 0 to 7, genes 8 and 9
    before = loss_pair_weights(data, label_mask, arguments_for(False), fit_indices)
    after = loss_pair_weights(data, label_mask, arguments_for(True), fit_indices)
    positive = data.outcomes > 0
    negative_count = ((~positive) & (before > 0)).sum(axis=1)[:8]
    assert np.allclose(np.where(~positive, after, 0).sum(axis=1)[:8], 0.2 * negative_count.mean())
    positive_count = (positive & (before > 0)).sum(axis=1)[:8]
    with_positives = positive_count > 0
    scale = positive_count[with_positives].mean() / positive_count[with_positives]
    assert np.allclose(np.where(positive, after, 0)[:8][with_positives], np.where(positive, before, 0)[:8][with_positives] * scale[:, None])
    assert np.allclose(after[8:], before[8:])


def test_the_input_normalisation_rescales_drugs_to_one_and_leaves_genes() -> None:
    from experiments.run_main_model import drug_input_normalised

    magnitudes = [np.array([1.0, 1.0, 0.5]), np.array([1.0]), np.array([0.25, 0.25]), np.array([0.0])]
    normalised = drug_input_normalised(["drug", "gene", "drug", "drug"], magnitudes)
    assert np.allclose(normalised[0], [0.4, 0.4, 0.2])
    assert np.allclose(normalised[1], [1.0])
    assert np.allclose(normalised[2], [0.5, 0.5])
    assert np.allclose(normalised[3], [0.0])


def test_the_weighted_start_uses_the_same_weights() -> None:
    data, label_mask = small_data()
    fit_indices = np.arange(10)
    weights = loss_pair_weights(data, label_mask, arguments_for(True), fit_indices)[fit_indices]
    targets = np.where(data.outcomes[fit_indices] > 0, 0.99, 0.0)
    expected = (weights * targets).sum(0) / weights.sum(0)
    assert np.allclose(weighted_constant_optimum(data, fit_indices, label_mask, arguments_for(True)), expected)


def test_off_the_flag_leaves_the_fingerprint_as_it_was() -> None:
    data, label_mask = small_data()
    data.edge_source, data.edge_target, data.edge_relation = np.array([0]), np.array([1]), np.array([0])
    data.perturbation_seeds, data.perturbation_signs, data.perturbation_magnitudes = [np.array([0])], [np.array([-1.0])], [np.array([1.0])]
    off = configuration_fingerprint(Namespace(head="sigmoid", equal_drug_shares=False), data, label_mask)
    assert off == configuration_fingerprint(Namespace(head="sigmoid"), data, label_mask)
    on = configuration_fingerprint(Namespace(head="sigmoid", equal_drug_shares=True), data, label_mask)
    assert on["argument:equal_drug_shares"] == "True"
    off = configuration_fingerprint(Namespace(head="sigmoid", normalise_drug_input=False), data, label_mask)
    assert off == configuration_fingerprint(Namespace(head="sigmoid"), data, label_mask)
