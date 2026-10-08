import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from mechanistic_pathway_learning.evaluation.experiment_data import ExperimentData, read_lockbox, restrict_to_perturbations
from mechanistic_pathway_learning.evaluation.negative_controls import (
    degree_preserving_rewiring, degree_stratified_row_permutation, duplicate_edge_count, fast_degree_preserving_rewiring,
    permute_symptom_labels_within_degree_strata)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
from draw_lockbox import draw_groups  # noqa: E402


def small_data(num_perturbations=6, label_selection_sha256="abc"):
    generator = np.random.default_rng(0)
    outcomes = (generator.random((num_perturbations, 2)) < 0.5).astype(float)
    return ExperimentData(
        node_ids=["a", "b", "c"], node_index={"a": 0, "b": 1, "c": 2}, node_types=np.array(["gene", "gene", "reaction"]), is_currency=np.zeros(3, dtype=bool),
        node_degree=np.array([1.0, 2.0, 3.0]), edge_source=np.array([0, 1]), edge_target=np.array([2, 2]), edge_relation=np.array([0, 0]), relation_types=["catalyzed_by"],
        symptoms=["anxiety", "fatigue"], perturbation_ids=[f"p{i}" for i in range(num_perturbations)], perturbation_labels=[f"P{i}" for i in range(num_perturbations)],
        perturbation_types=["gene"] * (num_perturbations - 1) + ["drug"], group_ids=["g0", "g0", "g1", "g2", "g2", "g3"][:num_perturbations],
        perturbation_seeds=[np.array([i % 3]) for i in range(num_perturbations)], perturbation_signs=[np.array([-1.0])] * num_perturbations,
        perturbation_magnitudes=[np.array([1.0])] * num_perturbations, outcomes=outcomes, weights=outcomes * 0.5, in_metabolic_layer=np.ones(num_perturbations, dtype=bool),
        frequencies=np.where(outcomes > 0, 0.4, np.nan), label_mask=np.ones_like(outcomes, dtype=bool),
        label_selection_summary={"sha256": label_selection_sha256},
    )


def test_restrict_keeps_the_rows_and_shares_the_graph():
    data = small_data()
    keep = np.array([True, False, True, True, False, True])
    restricted = restrict_to_perturbations(data, keep)
    assert restricted.perturbation_ids == ["p0", "p2", "p3", "p5"]
    assert restricted.group_ids == ["g0", "g1", "g2", "g3"]
    assert np.array_equal(restricted.outcomes, data.outcomes[keep]) and np.array_equal(restricted.weights, data.weights[keep])
    assert np.array_equal(restricted.label_mask, data.label_mask[keep])
    assert restricted.edge_source is data.edge_source
    assert np.array_equal(restricted.perturbation_degrees, data.perturbation_degrees[keep])


def write_lockbox(tmp_path, ids, sha256="abc"):
    path = tmp_path / "lockbox.json"
    path.write_text(json.dumps({"label_selection_sha256": sha256, "perturbation_ids": ids, "group_by": "disease_cluster_and_targets"}))
    return path


def test_read_lockbox_marks_members_and_refuses_mismatches(tmp_path):
    data = small_data()
    assert read_lockbox(write_lockbox(tmp_path, ["p2", "p5"]), data).tolist() == [False, False, True, False, False, True]
    with pytest.raises(ValueError, match="label selection"):
        read_lockbox(write_lockbox(tmp_path, ["p2"], sha256="other"), data)
    with pytest.raises(ValueError, match="not in the data"):
        read_lockbox(write_lockbox(tmp_path, ["p9"]), data)
    with pytest.raises(ValueError, match="both sides"):
        read_lockbox(write_lockbox(tmp_path, ["p3"]), data)  # p4 shares group g2


def test_draw_groups_reaches_the_share_per_stratum_and_skips_large_groups():
    group_ids = [f"g{i}" for i in range(40) for _ in range(1 + i % 3)]
    size_of_group = {f"g{i}": 1 + i % 3 for i in range(40)}
    stratum_of_group = {f"g{i}": "drug" if i < 10 else "other" for i in range(40)}
    chosen = draw_groups(group_ids, stratum_of_group, size_of_group, 0.2, seed=1, excluded_groups={"g2"})
    assert "g2" not in chosen
    for stratum in ("drug", "other"):
        total = sum(size for group, size in size_of_group.items() if stratum_of_group[group] == stratum)
        taken = sum(size_of_group[group] for group in chosen if stratum_of_group[group] == stratum)
        assert 0.2 * total <= taken < 0.2 * total + 3
    assert chosen == draw_groups(group_ids, stratum_of_group, size_of_group, 0.2, seed=1, excluded_groups={"g2"})


def test_row_permutation_reproduces_the_old_outcome_permutation_and_respects_the_partition():
    generator = np.random.default_rng(4)
    degrees = generator.integers(1, 50, size=300).astype(float)
    outcomes = generator.random((300, 4))
    old_generator = np.random.default_rng(7)
    edges = np.quantile(degrees, np.linspace(0, 1, 6))
    stratum = np.clip(np.searchsorted(edges, degrees, side="right") - 1, 0, 4)
    expected = outcomes.copy()
    for index in range(5):
        rows = np.where(stratum == index)[0]
        expected[rows] = outcomes[old_generator.permutation(rows)]
    assert np.array_equal(permute_symptom_labels_within_degree_strata(outcomes, degrees, random_seed=7), expected)
    partition = generator.random(300) < 0.2
    source = degree_stratified_row_permutation(degrees, random_seed=7, partition=partition)
    assert sorted(source.tolist()) == list(range(300))
    assert np.array_equal(partition[source], partition) and np.array_equal(stratum[source], stratum)


@pytest.mark.parametrize("rewire", [degree_preserving_rewiring, fast_degree_preserving_rewiring])
def test_rewiring_keeps_degrees_per_relation_and_makes_no_self_loop(rewire):
    generator = np.random.default_rng(3)
    edges = np.stack([generator.integers(0, 60, 400), generator.integers(60, 120, 400)])
    relation = generator.integers(0, 3, 400)
    rewired = rewire(edges, relation, num_swaps_per_edge=5, random_seed=0)
    for index in range(3):
        selected = relation == index
        assert np.array_equal(np.bincount(edges[0, selected], minlength=120), np.bincount(rewired[0, selected], minlength=120))
        assert np.array_equal(np.bincount(edges[1, selected], minlength=120), np.bincount(rewired[1, selected], minlength=120))
    assert not (rewired[0] == rewired[1]).any()
    assert (rewired != edges).any(axis=0).mean() > 0.9
    assert duplicate_edge_count(edges, relation) == len(edges[0]) - len({(s, t, r) for s, t, r in zip(edges[0], edges[1], relation)})


def test_fast_rewiring_adds_no_duplicate_edge():
    generator = np.random.default_rng(5)
    hub_heavy_sources = np.minimum(generator.geometric(0.2, 600), 30)  # a few hubs, so plain target swaps would repeat edges
    edges = np.stack([hub_heavy_sources, generator.integers(40, 80, 600)])
    relation = np.zeros(600, dtype=int)
    before = duplicate_edge_count(edges, relation)
    assert duplicate_edge_count(fast_degree_preserving_rewiring(edges, relation, num_swaps_per_edge=20, random_seed=1), relation) <= before



def test_scorer_average_precision_equals_sklearn_with_ties():
    from sklearn.metrics import average_precision_score

    from score_confirmatory import average_precision, macro_auprc, micro_auprc

    generator = np.random.default_rng(8)
    for _ in range(20):
        scores = np.round(generator.random(200), 1)  # many ties
        labels = (generator.random(200) < 0.2).astype(float)
        assert np.isclose(average_precision(scores, labels), average_precision_score(labels, scores))
    assert np.isnan(average_precision(np.array([0.1, 0.2]), np.array([0.0, 0.0])))
    predictions, outcomes = generator.random((50, 3)), (generator.random((50, 3)) < 0.3).astype(float)
    mask = np.ones_like(outcomes, dtype=bool)
    mask[0, 0] = False
    expected_macro = np.mean([average_precision_score(outcomes[mask[:, c], c], predictions[mask[:, c], c]) for c in (0, 2)])
    assert np.isclose(macro_auprc(predictions, outcomes, mask, [0, 2]), expected_macro)
    assert np.isclose(micro_auprc(predictions, outcomes, mask, [0, 2]), average_precision_score(outcomes[:, [0, 2]][mask[:, [0, 2]]], predictions[:, [0, 2]][mask[:, [0, 2]]]))


def test_holm_steps_down_and_stops_at_the_first_failure():
    from score_confirmatory import holm_rejections, one_sided_p

    rejected = holm_rejections({"a": 0.001, "b": 0.02, "c": 0.004, "d": 0.5}, alpha=0.025)
    assert rejected == {"a": True, "c": True, "b": False, "d": False}  # 0.001 <= 0.025/4, 0.004 <= 0.025/3, 0.02 > 0.025/2
    assert holm_rejections({"a": 0.01, "b": 0.0004}, alpha=0.001) == {"a": False, "b": True}
    assert holm_rejections({"a": 0.0004, "b": 0.0004}, alpha=0.001) == {"a": True, "b": True}
    assert one_sided_p(np.array([0.1, 0.2, -0.1, np.nan])) == (1 + 1) / (3 + 1)
