"""Tests: the three controls for the input-size objection (docs/label_source_review.md, 13.1; the user's decision of
11 October 2026 to keep 1 per target and counter the confound in the evaluation). The input-size baseline sees a
perturbation's number of seeds and total input magnitude and nothing else; the strata report a metric by how many
targets a drug has; the permuted-target control gives every drug another drug's targets of the same number."""
from types import SimpleNamespace

import numpy as np
import pytest

from experiments.run_main_model import permuted_drug_targets, refuse_a_negative_loss_weight
from mechanistic_pathway_learning.evaluation.input_size_strata import GENE_PERTURBATION_STRATUM, input_size_strata, strata_of_scored_rows
from mechanistic_pathway_learning.models.baselines.input_size_baseline import InputSizeBaseline, input_size_features


def test_the_baseline_reads_input_size_and_learns_a_size_signal_only_when_there_is_one() -> None:
    generator = np.random.default_rng(0)
    counts = np.repeat([1, 2, 3], 40)
    seeds = [np.arange(count) for count in counts]
    magnitudes = [np.ones(count) for count in counts]
    features = input_size_features(seeds, magnitudes)
    np.testing.assert_array_equal(features[:, 0], counts)
    np.testing.assert_array_equal(features[:, 1], counts)  # the input rule as built: total magnitude is the target count

    # symptom 0 rises with the target count, symptom 1 does not
    with_size = (generator.random(len(counts)) < counts / 6.0).astype(float)
    without = (generator.random(len(counts)) < 0.5).astype(float)
    outcomes = np.stack([with_size, without], axis=1)
    model = InputSizeBaseline().fit(features, outcomes)
    predictions = model.predict(features)
    by_count = {count: predictions[counts == count, 0].mean() for count in (1, 2, 3)}
    assert by_count[1] < by_count[2] < by_count[3]  # the size signal is found
    assert predictions[counts == 1, 1].std() == pytest.approx(0.0, abs=1e-9)  # one prediction per input size, always
    spread_without_signal = max(predictions[counts == count, 1].mean() for count in (1, 2, 3)) - min(predictions[counts == count, 1].mean() for count in (1, 2, 3))
    assert spread_without_signal < 0.1  # no size signal: near the base rate whatever the count
    assert np.isfinite(predictions).all() and (predictions > 0).all() and (predictions < 1).all()


def test_the_baseline_is_well_posed_under_each_input_rule_and_respects_a_label_mask() -> None:
    counts = np.repeat([1, 2, 4], 20)
    seeds = [np.arange(count) for count in counts]
    outcomes = np.stack([(counts > 1).astype(float), np.zeros(len(counts))], axis=1)
    as_built = input_size_features(seeds, [np.ones(count) for count in counts])
    normalised = input_size_features(seeds, [np.full(count, 1.0 / count) for count in counts])  # --normalise-drug-input: total 1
    occupancy = input_size_features(seeds, [np.linspace(0.2, 0.9, count) for count in counts])
    for features in (as_built, normalised, occupancy):
        predictions = InputSizeBaseline().fit(features, outcomes).predict(features)
        assert np.isfinite(predictions).all()
    # a masked pair is neither a positive nor a negative, so a symptom masked away everywhere reads its intercept only
    mask = np.ones_like(outcomes, dtype=bool)
    mask[counts > 1, 0] = False
    masked = InputSizeBaseline().fit(as_built, outcomes, training_label_mask=mask).predict(as_built)
    assert masked[:, 0].max() < 0.5  # the only labelled rows of symptom 0 are its negatives
    with pytest.raises(RuntimeError, match="fit before predict"):
        InputSizeBaseline().predict(as_built)


def test_the_strata_name_each_target_count_and_keep_the_genes_apart() -> None:
    types = ["gene", "drug", "drug", "drug", "drug"]
    labels, order = input_size_strata(types, np.array([1, 1, 2, 3, 7]))
    assert labels == [GENE_PERTURBATION_STRATUM, "drug_targets_1", "drug_targets_2", "drug_targets_3", "drug_targets_4_or_more"]
    assert order == [GENE_PERTURBATION_STRATUM, "drug_targets_1", "drug_targets_2", "drug_targets_3", "drug_targets_4_or_more"]
    assert input_size_strata(types, np.array([1, 1, 2, 3, 7]))[1] == order  # order is reading order, not first-seen
    # a count that no perturbation has is left out of the order
    assert "drug_targets_2" not in input_size_strata(["drug", "drug"], np.array([1, 3]))[1]
    # on a graph with drug entry nodes a drug is seeded on its own node; the strata come from its targets, so they do not change
    data = SimpleNamespace(perturbation_types=types, perturbation_seeds=[np.array([0]), np.array([1]), np.array([2]), np.array([3]), np.array([4])],
                           perturbation_target_seeds=[np.array([0]), np.array([1]), np.array([2, 3]), np.array([4, 5, 6]), np.arange(7)])
    assert strata_of_scored_rows(data, np.array([True, True, True, True, True]))[0] == labels
    assert strata_of_scored_rows(data, np.array([1, 3]))[0] == ["drug_targets_1", "drug_targets_3"]


def test_the_permuted_target_control_keeps_every_input_size_and_moves_every_identity() -> None:
    types = ["gene"] + ["drug"] * 7
    seeds = [np.array([0])] + [np.array([1]), np.array([2]), np.array([3]), np.array([4, 5]), np.array([6, 7]), np.array([8, 9]), np.array([10, 11, 12])]
    signs = [np.array([-1.0])] + [np.full(len(group), 0.5) for group in seeds[1:]]
    magnitudes = [np.ones(len(group)) for group in seeds]
    permuted_seeds, permuted_signs, permuted_magnitudes, moved = permuted_drug_targets(types, seeds, signs, magnitudes, random_seed=3)

    assert [len(group) for group in permuted_seeds] == [len(group) for group in seeds]  # every input size kept
    assert [float(np.abs(group).sum()) for group in permuted_magnitudes] == [float(np.abs(group).sum()) for group in magnitudes]
    assert np.array_equal(permuted_seeds[0], seeds[0]) and np.array_equal(permuted_signs[0], signs[0])  # the gene is untouched
    assert moved == 6  # the three one-target drugs and the three two-target drugs; the lone three-target drug has no partner
    assert np.array_equal(permuted_seeds[7], seeds[7])
    for position in range(1, 7):
        assert not np.array_equal(permuted_seeds[position], seeds[position])  # no drug keeps its own targets
        assert any(np.array_equal(permuted_seeds[position], seeds[other]) for other in range(1, 8))  # every drug holds some drug's targets
    assert sorted(tuple(group) for group in permuted_seeds[1:]) == sorted(tuple(group) for group in seeds[1:])  # the target sets are a permutation
    assert permuted_drug_targets(types, seeds, signs, magnitudes, random_seed=3)[0][1].tolist() == permuted_seeds[1].tolist()  # the seed decides it
    assert np.array_equal(seeds[1], np.array([1]))  # the caller's lists are not modified


def test_the_permuted_target_control_is_refused_where_a_drug_is_seeded_on_its_own_node(tmp_path) -> None:
    from pathlib import Path

    import pandas as pd

    from experiments.run_main_model import build_argument_parser
    from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
    from tests.test_clearance_edges import variant_with_enzymes
    from tests.test_drug_entry_nodes import DRUG_SEEDS, evidence_row, write_graph

    _, _, new_nodes, new_edges, relations, entry_table, _ = variant_with_enzymes()
    write_graph(tmp_path / "graph", new_nodes, new_edges, relations, entry_table)
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples])
                  for perturbation, (_, triples) in DRUG_SEEDS.items()]).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    as_nodes = load_experiment_data(tmp_path / "graph", evidence_directory, drug_entry="nodes")
    as_targets = load_experiment_data(tmp_path / "graph", evidence_directory, drug_entry="targets")
    assert as_nodes.drug_entry == "nodes" and as_targets.drug_entry is None
    # seeded on its node a drug has one seed, so an exchange would say nothing about its targets; the trainer refuses it.
    # (The exception is a drug that is itself a graph compound: it keeps the compound's own seeds and gets no node.)
    on_its_own_node = [len(seeds) for perturbation_id, seeds in zip(as_nodes.perturbation_ids, as_nodes.perturbation_seeds) if perturbation_id.startswith("drug_")]
    assert on_its_own_node == [1, 1] and sorted(len(seeds) for seeds in as_targets.perturbation_seeds) == [1, 2, 2]
    parser = build_argument_parser()
    arguments = parser.parse_args(["--graph-dir", str(tmp_path / "graph"), "--evidence-dir", str(evidence_directory), "--permute-drug-targets"])
    assert arguments.permute_drug_targets and Path(arguments.graph_dir).exists()


def test_a_negative_loss_weight_is_refused_because_it_inverts_the_gradient() -> None:
    """The user wrote "-.5" for the placebo negative weight on 12 October 2026. A weight multiplies a per-pair
    cross-entropy, so a negative one makes the optimiser maximise the loss on those pairs and train the model
    towards the outcome the placebo-controlled comparison argues against. The run stops instead."""
    for negative_weight, placebo_negative_weight in ((-0.5, None), (0.2, -0.5), (-1.0, -1.0)):
        arguments = SimpleNamespace(negative_weight=negative_weight, placebo_negative_weight=placebo_negative_weight)
        with pytest.raises(ValueError, match="inverts the gradient"):
            refuse_a_negative_loss_weight(arguments)

    # zero is allowed, because it sets the pairs aside rather than reversing them, and an unset placebo weight is
    # the default: those pairs then take --negative-weight like any other negative
    refuse_a_negative_loss_weight(SimpleNamespace(negative_weight=0.0, placebo_negative_weight=0.0))
    refuse_a_negative_loss_weight(SimpleNamespace(negative_weight=0.2, placebo_negative_weight=None))
    refuse_a_negative_loss_weight(SimpleNamespace(negative_weight=0.2, placebo_negative_weight=0.5))


def test_the_clamp_report_names_a_labelled_positive_that_carries_no_gradient():
    """Amendment item 19 of 12 October 2026: the clamp stays off and every development run reports whether a labelled
    positive sits below it, because the evidence that raised the question was measured on a head the family no longer
    holds."""
    from types import SimpleNamespace

    from experiments.run_main_model import clamped_probability_report
    from mechanistic_pathway_learning.models.soft_constraint_losses import PROBABILITY_EPSILON

    outcomes = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    data = SimpleNamespace(outcomes=outcomes)
    indices = np.array([0, 1, 2])
    # row 0's positive is below the clamp, row 1's is not, and row 2's negative is above it
    predictions = np.array([[PROBABILITY_EPSILON / 10, 0.5], [0.4, 0.5], [1.0 - PROBABILITY_EPSILON / 10, 0.9]])
    report = clamped_probability_report(predictions, data, indices, None)
    assert report["below_the_clamp"] == 1 and report["above_the_clamp"] == 1
    assert report["labelled_positives_below_the_clamp"] == 1
    assert report["labelled_negatives_above_the_clamp"] == 1
    assert report["reopens_bce_in_log_space"] is True
    # a pair the label selection set aside is not a labelled positive, so it does not reopen the option
    masked = np.array([[False, True], [True, True], [True, True]])
    quiet = clamped_probability_report(predictions, data, indices, masked)
    assert quiet["below_the_clamp"] == 1 and quiet["labelled_positives_below_the_clamp"] == 0
    assert quiet["reopens_bce_in_log_space"] is False
