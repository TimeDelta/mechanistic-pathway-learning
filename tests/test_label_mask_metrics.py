import numpy as np

from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import (
    bootstrap_interval, expected_calibration_error, hits_at_k, mean_reciprocal_rank, paired_bootstrap_macro_difference,
    per_symptom_auprc, per_symptom_auroc, scorable_symptom)


def example(seed=0):
    generator = np.random.default_rng(seed)
    predictions = generator.random((40, 3))
    outcomes = (generator.random((40, 3)) < 0.3).astype(float)
    mask = np.ones_like(outcomes, dtype=bool)
    positive_cells = np.argwhere(outcomes > 0)
    for row, column in positive_cells[::3]:  # set aside every third positive
        mask[row, column] = False
    return predictions, outcomes, mask


def test_masked_per_symptom_metrics_equal_the_metrics_with_those_rows_removed():
    predictions, outcomes, mask = example()
    for symptom in range(3):
        kept = mask[:, symptom]
        assert per_symptom_auprc(predictions, outcomes, symptom, mask) == per_symptom_auprc(predictions[kept], outcomes[kept], symptom)
        assert per_symptom_auroc(predictions, outcomes, symptom, mask) == per_symptom_auroc(predictions[kept], outcomes[kept], symptom)


def test_a_masked_positive_is_not_counted_as_a_negative():
    predictions = np.array([[0.8], [0.9], [0.1], [0.2]])  # the set-aside positive ranks first
    outcomes = np.array([[1.0], [1.0], [0.0], [0.0]])
    mask = np.array([[True], [False], [True], [True]])
    relabelled = outcomes.copy()
    relabelled[1, 0] = 0.0
    assert per_symptom_auprc(predictions, outcomes, 0, mask) == 1.0
    assert per_symptom_auprc(predictions, relabelled, 0) < 1.0


def test_ranking_and_calibration_drop_masked_symptoms_of_a_perturbation():
    predictions = np.array([[0.9, 0.5, 0.1]])
    outcomes = np.array([[1.0, 0.0, 1.0]])
    mask = np.array([[False, True, True]])  # the top-ranked positive is set aside
    assert mean_reciprocal_rank(predictions, outcomes) == 1.0
    assert mean_reciprocal_rank(predictions, outcomes, mask) == 0.5
    assert hits_at_k(predictions, outcomes, 1, mask) == 0.0
    assert expected_calibration_error(predictions, outcomes, mask=mask) == expected_calibration_error(predictions[:, 1:], outcomes[:, 1:])


def test_scorable_symptom_counts_labelled_rows_only():
    outcomes = np.array([[1.0], [1.0], [0.0]])
    mask = np.array([[True], [True], [False]])
    assert scorable_symptom(outcomes, 0, 1) is True
    assert scorable_symptom(outcomes, 0, 1, mask) is False  # every labelled row is positive


def test_bootstrap_and_paired_difference_accept_the_mask_and_reduce_to_the_unmasked_form():
    predictions, outcomes, mask = example(1)
    unmasked = bootstrap_interval(lambda p, y, m=None: per_symptom_auprc(p, y, 0, m), predictions, outcomes, num_bootstrap=50)
    all_true = bootstrap_interval(lambda p, y, m=None: per_symptom_auprc(p, y, 0, m), predictions, outcomes, num_bootstrap=50, mask=np.ones_like(mask))
    assert (unmasked.point, unmasked.lower, unmasked.upper) == (all_true.point, all_true.lower, all_true.upper)
    masked = bootstrap_interval(lambda p, y, m=None: per_symptom_auprc(p, y, 0, m), predictions, outcomes, num_bootstrap=50, mask=mask)
    assert masked.point == per_symptom_auprc(predictions, outcomes, 0, mask)
    other = np.random.default_rng(2).random(predictions.shape)
    plain = paired_bootstrap_macro_difference(predictions, other, outcomes, num_bootstrap=50, minimum_positives=2)
    with_ones = paired_bootstrap_macro_difference(predictions, other, outcomes, num_bootstrap=50, minimum_positives=2, mask=np.ones_like(mask))
    assert plain == with_ones
    with_mask = paired_bootstrap_macro_difference(predictions, other, outcomes, num_bootstrap=50, minimum_positives=2, mask=mask)
    assert np.isfinite(with_mask["difference"])


def test_label_selection_masks_positive_pairs_in_the_loader(tmp_path):
    import json

    import pandas as pd

    from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

    graph = tmp_path / "graph"
    evidence = tmp_path / "evidence"
    graph.mkdir()
    evidence.mkdir()
    pd.DataFrame({"node_id": ["GENE:A", "GENE:B", "R1"], "node_type": ["gene", "gene", "reaction"], "degree": [1, 1, 2],
                  "is_currency": [False, False, False]}).to_parquet(graph / "nodes.parquet")
    pd.DataFrame({"source_id": ["GENE:A", "GENE:B"], "target_id": ["R1", "R1"], "relation_type": ["catalyzed_by"] * 2}).to_parquet(graph / "edges.parquet")
    (graph / "relation_types.json").write_text(json.dumps(["catalyzed_by"]))
    rows = []
    for gene, symptom in [("A", "anxiety"), ("A", "fatigue"), ("B", "anxiety")]:
        rows.append({"perturbation_id": gene, "perturbation_type": "gene", "perturbation_label": gene, "group_id": gene, "symptom": symptom,
                     "relation": "induces", "grade": "A", "weight": 1.0, "in_metabolic_layer": True,
                     "perturbation_nodes": json.dumps([[f"GENE:{gene}", -1.0, 1.0]])})
    pd.DataFrame(rows).to_parquet(evidence / "evidence_records.parquet")
    selection = tmp_path / "selection.parquet"
    pd.DataFrame({"perturbation_id": ["A", "A", "B"], "symptom": ["anxiety", "fatigue", "anxiety"], "keep": [True, False, True]}).to_parquet(selection)
    data = load_experiment_data(graph, evidence, label_selection=selection)
    fatigue = data.symptoms.index("fatigue")
    row_a = data.perturbation_ids.index("A")
    assert data.outcomes[row_a, fatigue] == 1.0  # still a positive in outcomes
    assert not data.label_mask[row_a, fatigue]
    assert data.label_mask.sum() == data.label_mask.size - 1
    assert data.label_selection_summary["masked_pairs"] == 1
    assert load_experiment_data(graph, evidence).label_mask is None



def test_mask_grades_masks_grade_c_only_pairs_and_leaves_positives(tmp_path):
    """--mask-grades C: a pair with grade C evidence only is neither positive nor negative; a pair with grade A and C stays a positive."""
    import json
    import subprocess
    import sys

    import pandas as pd

    from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data

    graph = tmp_path / "graph"
    evidence = tmp_path / "evidence"
    graph.mkdir()
    evidence.mkdir()
    pd.DataFrame({"node_id": ["GENE:A", "GENE:B", "R1"], "node_type": ["gene", "gene", "reaction"], "degree": [1, 1, 2],
                  "is_currency": [False, False, False]}).to_parquet(graph / "nodes.parquet")
    pd.DataFrame({"source_id": ["GENE:A", "GENE:B"], "target_id": ["R1", "R1"], "relation_type": ["catalyzed_by"] * 2}).to_parquet(graph / "edges.parquet")
    (graph / "relation_types.json").write_text(json.dumps(["catalyzed_by"]))
    rows = []
    for gene, symptom, grade in [("A", "anxiety", "A"), ("A", "anxiety", "C"), ("A", "fatigue", "C"), ("B", "anxiety", "A"), ("B", "apathy", "A")]:
        rows.append({"perturbation_id": gene, "perturbation_type": "gene", "perturbation_label": gene, "group_id": gene, "symptom": symptom,
                     "relation": "induces", "grade": grade, "weight": 1.0, "in_metabolic_layer": True, "label_frequency": 0.5, "source": "HPO",
                     "perturbation_nodes": json.dumps([[f"GENE:{gene}", -1.0, 1.0]])})
    pd.DataFrame(rows).to_parquet(evidence / "evidence_records.parquet")
    selection = tmp_path / "selection.parquet"
    subprocess.run([sys.executable, "experiments/build_label_selection.py", "--evidence-dir", str(evidence), "--mask-grades", "C", "--output", str(selection)],
                   check=True, capture_output=True)
    data = load_experiment_data(graph, evidence, symptoms=["anxiety", "apathy", "fatigue"], label_selection=selection)
    row_a, fatigue, anxiety = data.perturbation_ids.index("A"), data.symptoms.index("fatigue"), data.symptoms.index("anxiety")
    assert data.outcomes[row_a, fatigue] == 0 and not data.label_mask[row_a, fatigue]  # grade C only: masked
    assert data.outcomes[row_a, anxiety] == 1 and data.label_mask[row_a, anxiety]  # grade A and C: a kept positive
    assert data.label_mask.sum() == data.label_mask.size - 1
    assert data.label_selection_summary["masked_negative_pairs"] == 1 and data.label_selection_summary["masked_pairs"] == 0
    summary = json.loads(selection.with_suffix(".summary.json").read_text())
    assert summary["positive_pairs"] == 3 and summary["masked_negative_pairs"] == 1

def test_micro_auprc_pools_every_labelled_pair_and_skips_masked_ones():
    from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import micro_auprc
    outcomes = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    predictions = np.array([[0.9, 0.1], [0.2, 0.8], [0.3, 0.05]])
    assert micro_auprc(predictions, outcomes) == 1.0
    worse = predictions.copy()
    worse[2, 0] = 0.95  # a negative ranked first
    assert micro_auprc(worse, outcomes) < 1.0
    mask = np.ones_like(outcomes, dtype=bool)
    mask[2, 0] = False  # the same pair masked out is neither a positive nor a negative
    assert micro_auprc(worse, outcomes, mask) == 1.0
    assert np.isnan(micro_auprc(predictions, np.zeros_like(outcomes)))


def test_micro_auprc_counts_a_symptom_below_the_macro_cutoff():
    from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import micro_auprc, paired_bootstrap_macro_difference, paired_bootstrap_micro_difference
    generator = np.random.default_rng(0)
    outcomes = np.zeros((60, 2))
    outcomes[:20, 0] = 1.0
    outcomes[0, 1] = 1.0  # one positive: left out of the macro average, inside the micro average
    good = generator.random(outcomes.shape) * 0.5 + outcomes * 0.5
    rare_missed = good.copy()
    rare_missed[0, 1] = 0.0
    macro = paired_bootstrap_macro_difference(good, rare_missed, outcomes, num_bootstrap=50)
    micro = paired_bootstrap_micro_difference(good, rare_missed, outcomes, num_bootstrap=50)
    assert macro["difference"] == 0.0
    assert micro["difference"] > 0.0
    assert micro_auprc(good, outcomes) > micro_auprc(rare_missed, outcomes)


def test_drugs_join_the_groups_of_the_genes_they_target_and_of_drugs_sharing_a_target():
    from mechanistic_pathway_learning.evaluation.experiment_data import merge_drugs_with_their_targets
    node_ids = ["GENE:A", "GENE:B", "GENE:C", "GENE:D"]
    perturbation_ids = ["A", "B", "C", "drug1", "drug2", "drug3"]
    perturbation_types = ["gene", "gene", "gene", "drug", "drug", "drug"]
    group_ids = ["cluster:A", "cluster:B", "cluster:C", "T1", "T2", "T3"]
    seeds = [np.array([0]), np.array([1]), np.array([2]), np.array([0]), np.array([3]), np.array([3, 1])]
    merged = merge_drugs_with_their_targets(perturbation_ids, perturbation_types, group_ids, seeds, node_ids)
    group = dict(zip(perturbation_ids, merged))
    assert group["drug1"] == group["A"] == "cluster:A"  # drug1 targets gene A
    assert group["drug2"] == group["drug3"] == group["B"]  # drug2 and drug3 share node D; drug3 also targets gene B
    assert group["C"] == "cluster:C"  # untouched
    assert len(set(merged)) == 3


def test_popularity_base_rate_leaves_set_aside_pairs_out() -> None:
    from mechanistic_pathway_learning.models.baselines.popularity_baseline import PopularityBaseline
    outcomes = np.array([[1, 0], [1, 1], [0, 0], [0, 1]], dtype=float)
    mask = np.array([[True, True], [False, True], [True, True], [True, False]])  # two positives set aside
    unmasked = PopularityBaseline().fit(outcomes * mask)
    masked = PopularityBaseline().fit(outcomes * mask, training_label_mask=mask)
    assert np.allclose(unmasked.symptom_base_rate, [0.25, 0.25])  # a set-aside positive counted as a negative
    assert np.allclose(masked.symptom_base_rate, [1 / 3, 1 / 3])  # one kept positive over three labelled pairs


def test_bootstrap_leaves_out_resamples_without_a_positive() -> None:
    import numpy as np

    from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import bootstrap_interval, per_symptom_auprc

    generator = np.random.default_rng(1)
    outcomes = np.zeros((100, 1))
    outcomes[:5, 0] = 1.0  # five positives: about 0.6 percent of resamples hold none
    predictions = generator.random((100, 1)) + outcomes  # positives rank first
    assert np.isnan(per_symptom_auprc(predictions[5:], outcomes[5:], 0))
    interval = bootstrap_interval(lambda p, y: per_symptom_auprc(p, y, 0), predictions, outcomes, num_bootstrap=2000)
    assert interval.lower == 1.0  # every resample with a positive ranks it first; a resample without one is not a 0
