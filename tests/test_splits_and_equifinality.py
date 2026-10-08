"""Tests for splits, evidence grading, currency tagging and equifinality indices."""
import math
from datetime import date

from mechanistic_pathway_learning.evaluation.equifinality_independence_and_convergence import (
    active_modules_for_symptom,
    convergence_index,
    downstream_reachable_nodes,
    independence_index,
    jaccard_index,
)
from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
    assign_grouped_folds,
    pathway_wise_holdout_sets,
    time_split,
)
from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord, assign_evidence_grade, loss_weight_for_record
from mechanistic_pathway_learning.graph.tag_currency_metabolites import strip_compartment_suffix, tag_currency_metabolites


def test_grouped_folds_keep_shared_target_drugs_together() -> None:
    drugs = ["haloperidol", "chlorpromazine", "fluoxetine", "sertraline", "diazepam", "ketamine", "memantine"]
    targets = ["DRD2", "DRD2", "SLC6A4", "SLC6A4", "GABRA1", "GRIN1", "GRIN1"]
    fold_by_drug = assign_grouped_folds(drugs, targets, num_folds=3, random_seed=1)
    assert set(fold_by_drug) == set(drugs)
    assert fold_by_drug["haloperidol"] == fold_by_drug["chlorpromazine"]
    assert fold_by_drug["fluoxetine"] == fold_by_drug["sertraline"]
    assert fold_by_drug["ketamine"] == fold_by_drug["memantine"]
    assert max(fold_by_drug.values()) < 3


def test_pathway_wise_holdout_groups_by_module() -> None:
    held_out = pathway_wise_holdout_sets({"HMBS": "heme", "PPOX": "heme", "ATP7B": "copper"})
    assert sorted(held_out["heme"]) == ["HMBS", "PPOX"]
    assert held_out["copper"] == ["ATP7B"]


def test_time_split_respects_cutoff_and_drops_undated() -> None:
    train, test = time_split(["a", "b", "c"], [date(2014, 1, 1), date(2020, 6, 1), None], date(2015, 12, 31))
    assert train == ["a"] and test == ["b"]


def test_evidence_grades_and_weights() -> None:
    monogenic = EvidenceRecord("HMBS", "psychosis", "induces", "monogenic", has_omim_clinical_synopsis=True)
    assert assign_evidence_grade(monogenic) == "A" and loss_weight_for_record(monogenic) == 1.0
    weak_monogenic = EvidenceRecord("GENE", "anxiety", "induces", "monogenic", independent_case_series_count=1)
    assert assign_evidence_grade(weak_monogenic) == "C" and loss_weight_for_record(weak_monogenic) == 0.0
    drug = EvidenceRecord("ketamine", "psychosis", "induces", "pharmacological", cns_penetrant=True, has_dominant_target=True, label_event_frequency=0.05)
    assert assign_evidence_grade(drug) == "B" and abs(loss_weight_for_record(drug) - 0.3) < 1e-9
    literature = EvidenceRecord("KMO", "depressed_mood", "induces", "literature", predication_type="causes")
    assert assign_evidence_grade(literature) == "E" and loss_weight_for_record(literature) == 0.10


def test_currency_tagging_by_name_and_degree() -> None:
    names = {"atp_c": "ATP", "m1": "dopamine", "m2": "serotonin", "m3": "h2o"}
    degrees = {"atp_c": 900, "m1": 5, "m2": 4, "m3": 800}
    tagged = tag_currency_metabolites(names, degrees, degree_quantile_threshold=0.75)
    assert "atp_c" in tagged and "m3" in tagged
    assert "m1" not in tagged and "m2" not in tagged
    assert strip_compartment_suffix("atp[c]") == "atp" and strip_compartment_suffix("MAM01371c") == "MAM01371c"


def test_equifinality_indices() -> None:
    supports = {"heme": {"ALAS1", "HMBS", "PPOX"}, "copper": {"ATP7B", "CP"}, "overlap": {"HMBS", "PPOX", "ATP7B"}}
    assert jaccard_index(supports["heme"], supports["copper"]) == 0.0
    assert independence_index(supports, ["heme", "copper"]) == 1.0
    assert independence_index(supports, ["heme", "overlap"]) < 1.0
    assert math.isnan(independence_index(supports, ["heme"]))
    adjacency = {"HMBS": ["porphobilinogen"], "ATP7B": ["copper_ion"], "porphobilinogen": ["dopamine"], "copper_ion": ["dopamine"]}
    assert downstream_reachable_nodes(adjacency, {"HMBS"}, max_hops=1) == {"porphobilinogen"}
    assert convergence_index(supports, ["heme", "copper"], adjacency, max_hops=1) == 0.0
    assert convergence_index(supports, ["heme", "copper"], adjacency, max_hops=2) > 0.0
    assert active_modules_for_symptom({"heme": 0.8, "copper": 0.4}, link_threshold=0.5) == ["heme"]


def test_primary_subsystem_routes_catch_all_bins_to_pathways_and_holdout_masks() -> None:
    from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import (
        perturbations_anchored_in_module,
        primary_subsystem_by_gene_node,
        subsystem_holdout_masks,
    )

    # nodes 0-2 genes, 3-6 reactions; relation 2 is catalyzed_by
    node_subsystem = ["", "", "", "Transport reactions", "Urea cycle", "Urea cycle", "Porphyrin metabolism"]
    edge_source = [0, 0, 1, 2, 2]
    edge_target = [3, 4, 5, 6, 3]
    edge_relation = [2, 2, 2, 2, 2]
    primary = primary_subsystem_by_gene_node(node_subsystem, edge_source, edge_target, edge_relation, catalyzed_by_relation_index=2)
    assert primary[0] == "Urea cycle"  # one transport and one urea cycle reaction: the pathway wins over the catch-all
    assert primary[1] == "Urea cycle" and primary[2] == "Porphyrin metabolism"
    seeds = [[0], [1], [2], [0, 2]]
    positives = [2.0, 3.0, 5.0, 1.0]
    masks = subsystem_holdout_masks(seeds, primary, positives, min_holdout_positives=5)
    assert masks["Urea cycle"] == [True, True, False, True]  # a drug writing onto a urea cycle gene leaves with the pathway
    assert masks["Porphyrin metabolism"] == [False, False, True, True]
    assert perturbations_anchored_in_module(seeds, {1}) == [False, True, False, False]
    assert "Transport reactions" not in subsystem_holdout_masks(seeds, {0: "Transport reactions", 1: "Transport reactions", 2: "Transport reactions"}, positives, 1)


def test_pathway_wise_training_mask_leaves_out_group_partners() -> None:
    from mechanistic_pathway_learning.evaluation.perturbation_wise_and_pathway_wise_splits import training_mask_without_group_partners

    group_ids = ["cluster:A", "cluster:A", "cluster:B", "cluster:C", "cluster:B"]
    # perturbation 0 is held out; 1 shares its disease cluster and must not train; 2 to 4 train
    assert training_mask_without_group_partners([True, False, False, False, False], group_ids) == [False, False, True, True, True]
    # a hold-out of whole groups trains on the complement
    assert training_mask_without_group_partners([False, False, True, False, True], group_ids) == [True, True, False, True, False]


def test_baseline_pathway_hold_out_fits_without_group_partners() -> None:
    import importlib.util
    from pathlib import Path
    from types import SimpleNamespace

    import numpy as np

    spec = importlib.util.spec_from_file_location("run_baselines_module", Path(__file__).resolve().parents[1] / "experiments" / "run_baselines.py")
    run_baselines = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_baselines)
    # popularity predicts the training base rate, so a held-out row's prediction shows which rows were fitted
    outcomes = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
    data = SimpleNamespace(perturbation_degrees=np.ones(4), group_ids=["g0", "g0", "g1", "g2"])
    test = np.array([True, False, False, False])
    without_groups, _, per_fold, _ = run_baselines.run_split(data, outcomes, [test], "popularity", 0.3, None, min_fold_size_for_macro=100)
    with_groups, _, per_fold_grouped, _ = run_baselines.run_split(data, outcomes, [test], "popularity", 0.3, None, min_fold_size_for_macro=100, group_ids=data.group_ids)
    np.testing.assert_allclose(without_groups[0], [1 / 3, 2 / 3])  # rows 1 to 3 fitted, row 1 is the held-out row's partner
    np.testing.assert_allclose(with_groups[0], [0.0, 1.0])  # rows 2 and 3 only
    assert per_fold[0]["num_group_partners_left_out"] == 0 and per_fold_grouped[0]["num_group_partners_left_out"] == 1
