"""Soft priors from the evidence below the label bar (mechanistic_pathway_learning/evidence/soft_priors.py)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.soft_priors import (
    DEFAULT_TIERS_PATH,
    Tier,
    apply_soft_priors,
    eligible_cells,
    literature_popularity,
    permute_soft_priors,
    read_tiers,
    record_key,
    resolve_pairs,
    soft_evidence_from_records,
    soft_evidence_from_reports,
    soft_prior_matrices,
    summarise,
    tempered_count,
)

TIERS = {
    ("monogenic", "positive"): Tier("monogenic", "positive", 0.75, 0.5),
    ("human_association", "positive"): Tier("human_association", "positive", 0.70, 0.5),
    ("human_association", "negative"): Tier("human_association", "negative", 0.0, 0.5),
    ("animal", "positive"): Tier("animal", "positive", 0.60, 0.3),
    ("animal", "negative"): Tier("animal", "negative", 0.0, 0.3),
    ("literature", "positive"): Tier("literature", "positive", 0.55, 0.1),
}


def report(perturbation, symptom, evidence_class, source, value=1, relation="induces", pmid=None, references="", record="", report_id=None):
    return {"report_id": report_id or f"{source}|{perturbation}|{symptom}|{record or pmid}", "perturbation_id": perturbation, "symptom": symptom, "relation": relation,
            "evidence_class": evidence_class, "source": source, "report_value": value, "pmid": pmid, "references": references, "source_record_id": record}


def test_the_committed_tier_table_reads_and_names_every_soft_class():
    tiers = read_tiers(DEFAULT_TIERS_PATH)
    assert {key[0] for key in tiers} == {"monogenic", "human_association", "animal", "literature"}
    assert tiers[("monogenic", "negative")].unit_strength == 0.5, "an HPO NOT annotation is the strongest curated absence in the inputs"
    assert {key for key, tier in tiers.items() if tier.direction == "negative"} == {("monogenic", "negative"), ("human_association", "negative"), ("animal", "negative")}, \
        "every class that can state an absence has a tier for it; the literature sources state none"
    assert tiers[("human_association", "negative")].prior_probability == 0.0 and tiers[("human_association", "negative")].unit_strength == 0.5, "the registered placebo weight"
    assert tiers[("literature", "positive")].unit_strength < tiers[("animal", "positive")].unit_strength < tiers[("human_association", "positive")].unit_strength
    assert all(0.5 < tier.prior_probability < 1.0 for tier in tiers.values() if tier.direction == "positive"), "a positive prior pulls above a coin toss and below a label"


def test_a_tier_table_with_a_positive_absence_target_is_refused(tmp_path):
    path = tmp_path / "tiers.csv"
    path.write_text("evidence_class,direction,prior_probability,unit_strength,note\nanimal,negative,0.2,0.3,x\n")
    with pytest.raises(ValueError, match="trained toward 0"):
        read_tiers(path)
    path.write_text("evidence_class,direction,prior_probability,unit_strength,note\nanimal,sideways,0.2,0.3,x\n")
    with pytest.raises(ValueError, match="direction"):
        read_tiers(path)


def test_tempering_caps_the_count_and_takes_its_root():
    assert tempered_count(1, cap=2, exponent=0.5) == pytest.approx(1.0)
    assert tempered_count(2, cap=2, exponent=0.5) == pytest.approx(np.sqrt(2))
    assert tempered_count(50, cap=2, exponent=0.5) == pytest.approx(np.sqrt(2)), "the fiftieth record adds nothing"
    assert tempered_count(3, cap=3, exponent=1.0) == pytest.approx(3.0)
    assert list(tempered_count([1, 4], cap=2, exponent=0.0)) == [1.0, 1.0], "exponent 0 is an indicator"


def test_the_record_key_prefers_the_pmid_then_the_cited_reference_then_the_source_record():
    rows = pd.DataFrame([report("G", "s", "literature", "PubTator3-cause", pmid="123", references="PMID:999"),
                         report("G", "s", "monogenic", "HPO-OMIM", references="PMID:555;PMID:556", record="OMIM:1"),
                         report("D", "s", "human_association", "ClinicalTrials-arms", record="NCT1"),
                         report("D", "s", "human_association", "FAERS", report_id="FAERS|D|s")])
    keys = [record_key(row) for row in rows.itertuples(index=False)]
    assert keys == ["PMID:123", "PMID:555", "ClinicalTrials-arms|NCT1", "FAERS|D|s"]


def test_report_evidence_counts_distinct_records_per_direction_and_drops_co_mention_literature():
    reports = pd.DataFrame([
        report("D", "parkinsonism", "human_association", "ClinicalTrials-arms", 1, record="NCT1"),
        report("D", "parkinsonism", "human_association", "ClinicalTrials-arms", 1, record="NCT2"),
        report("D", "parkinsonism", "human_association", "ClinicalTrials-arms", 1, record="NCT2"),  # the same trial twice (two arms)
        report("D", "anxiety", "human_association", "ClinicalTrials-arms", 0, record="NCT1"),
        report("G", "anxiety", "literature", "PubTator3-cause", 1, pmid="1"),
        report("G", "anxiety", "literature", "PubTator3-cause", 1, pmid="2"),
        report("G", "anxiety", "literature", "PubTator3-associate", 1, pmid="3", relation="associated_with"),
        report("G", "anxiety", "literature", "PubTator3-treat", 1, pmid="4", relation="relieves"),
        report("G", "psychosis", "literature", "PubTator3-associate", 1, pmid="5"),  # an induces row from a co-mention source: not admitted
    ])
    evidence = soft_evidence_from_reports(reports, relation="induces")
    rows = {(row.perturbation_id, row.symptom, row.direction): (row.n_distinct, row.sources) for row in evidence.itertuples(index=False)}
    assert rows[("D", "parkinsonism", "positive")] == (2, "ClinicalTrials-arms")
    assert rows[("D", "anxiety", "negative")] == (1, "ClinicalTrials-arms")
    assert rows[("G", "anxiety", "positive")] == (2, "PubTator3-cause"), "two PMIDs; the co-mention and the relieves rows are not induces evidence"
    assert ("G", "psychosis", "positive") not in rows
    assert soft_evidence_from_reports(pd.DataFrame(), "induces").empty


def test_record_evidence_takes_the_soft_grades_with_their_reference_counts():
    records = pd.DataFrame([
        {"perturbation_id": "G1", "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "C", "positive_report_count": 3, "negative_report_count": 0,
         "distinct_pubmed_reference_count": 2, "source": "HPO genes_to_phenotype"},
        {"perturbation_id": "G2", "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "A", "positive_report_count": 1, "negative_report_count": 0,
         "distinct_pubmed_reference_count": 1, "source": "HPO"},
        {"perturbation_id": "G3", "symptom": "anxiety", "relation": "induces", "evidence_class": "monogenic", "grade": "C", "positive_report_count": 1, "negative_report_count": 1,
         "distinct_pubmed_reference_count": None, "source": "HPO"},
        {"perturbation_id": "G1", "symptom": "anxiety", "relation": "relieves", "evidence_class": "monogenic", "grade": "C", "positive_report_count": 1, "negative_report_count": 0,
         "distinct_pubmed_reference_count": 0, "source": "HPO"},
    ])
    evidence = soft_evidence_from_records(records, grades=("C",), relation="induces")
    rows = {(row.perturbation_id, row.direction): row.n_distinct for row in evidence.itertuples(index=False)}
    assert rows == {("G1", "positive"): 3, ("G3", "positive"): 1, ("G3", "negative"): 1}, "grade A is a label; the relieves row is another relation; G3 is contested"


def test_resolution_takes_the_highest_ranked_class_and_sets_contested_pairs_aside():
    evidence = pd.DataFrame([
        {"perturbation_id": "D", "symptom": "anxiety", "evidence_class": "literature", "direction": "positive", "n_distinct": 9, "sources": "PubTator3-cause"},
        {"perturbation_id": "D", "symptom": "anxiety", "evidence_class": "human_association", "direction": "positive", "n_distinct": 1, "sources": "ClinicalTrials-arms"},
        {"perturbation_id": "D", "symptom": "psychosis", "evidence_class": "human_association", "direction": "negative", "n_distinct": 2, "sources": "ClinicalTrials-arms"},
        {"perturbation_id": "D", "symptom": "psychosis", "evidence_class": "animal", "direction": "positive", "n_distinct": 1, "sources": "IMPC"},
        {"perturbation_id": "G", "symptom": "anxiety", "evidence_class": "animal", "direction": "negative", "n_distinct": 1, "sources": "IMPC"},
        {"perturbation_id": "G", "symptom": "insomnia", "evidence_class": "open_targets", "direction": "positive", "n_distinct": 1, "sources": "OT"},
    ])
    resolved = resolve_pairs(evidence, TIERS, tempering_exponent=0.5, count_cap=2)
    pairs = {(row.perturbation_id, row.symptom): row for row in resolved.pairs.itertuples(index=False)}
    anxiety = pairs[("D", "anxiety")]
    assert anxiety.evidence_class == "human_association" and anxiety.prior_probability == 0.70 and anxiety.strength == pytest.approx(0.5), "the trial decides, at one record"
    assert ("D", "psychosis") not in pairs and resolved.contested.iloc[0].to_dict() == {"perturbation_id": "D", "symptom": "psychosis", "classes": "animal;human_association", "directions": "negative;positive"}
    absence = pairs[("G", "anxiety")]
    assert absence.direction == "negative" and absence.prior_probability == 0.0 and absence.strength == pytest.approx(0.3)
    assert resolved.without_a_tier == {"open_targets": 1} and ("G", "insomnia") not in pairs


def test_the_same_pair_from_two_tables_sums_its_records_before_tempering():
    first = pd.DataFrame([{"perturbation_id": "G", "symptom": "anxiety", "evidence_class": "literature", "direction": "positive", "n_distinct": 1, "sources": "CTD-curated"}])
    second = pd.DataFrame([{"perturbation_id": "G", "symptom": "anxiety", "evidence_class": "literature", "direction": "positive", "n_distinct": 1, "sources": "PubTator3-cause"}])
    resolved = resolve_pairs(pd.concat([first, second]), TIERS, tempering_exponent=0.5, count_cap=2)
    row = resolved.pairs.iloc[0]
    assert row.n_distinct == 2 and row.strength == pytest.approx(0.1 * np.sqrt(2)) and row.sources == "CTD-curated;PubTator3-cause"


def small_universe():
    perturbations, symptoms = ["D", "G", "H"], ["anxiety", "psychosis", "insomnia"]
    outcomes = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    label_mask = np.array([[True, True, True], [True, False, True], [True, True, False]])
    masked_for_soft = np.array([[False, False, False], [False, True, False], [False, False, False]])  # G/psychosis set aside for grade C only
    evidenced_negative = np.array([[False, False, False], [False, False, False], [False, False, True]])
    return perturbations, symptoms, outcomes, label_mask, masked_for_soft, evidenced_negative


def test_eligible_cells_exclude_labels_evidenced_negatives_and_placebo_set_asides():
    _, _, outcomes, label_mask, masked_for_soft, evidenced_negative = small_universe()
    eligible = eligible_cells(outcomes, label_mask, evidenced_negative, masked_for_soft)
    assert eligible.tolist() == [[False, True, True], [True, True, True], [True, False, False]], "G/psychosis is eligible (set aside for grade C); H/insomnia is a placebo set-aside"
    assert eligible_cells(outcomes, None, None, None).tolist() == (outcomes == 0).tolist()


def test_matrices_align_to_the_universe_and_drop_pairs_outside_it():
    perturbations, symptoms, *_ = small_universe()
    resolved = resolve_pairs(pd.DataFrame([
        {"perturbation_id": "D", "symptom": "psychosis", "evidence_class": "human_association", "direction": "negative", "n_distinct": 1, "sources": "x"},
        {"perturbation_id": "G", "symptom": "psychosis", "evidence_class": "monogenic", "direction": "positive", "n_distinct": 2, "sources": "x"},
        {"perturbation_id": "ZZ", "symptom": "psychosis", "evidence_class": "monogenic", "direction": "positive", "n_distinct": 2, "sources": "x"},
        {"perturbation_id": "G", "symptom": "not_a_symptom", "evidence_class": "monogenic", "direction": "positive", "n_distinct": 2, "sources": "x"},
    ]), TIERS)
    probability, strength = soft_prior_matrices(resolved, perturbations, symptoms)
    assert np.isnan(probability[0, 0]) and probability[0, 1] == 0.0 and strength[0, 1] == pytest.approx(0.5)
    assert probability[1, 1] == 0.75 and strength[1, 1] == pytest.approx(0.5 * np.sqrt(2)) and (strength > 0).sum() == 2


def test_applying_priors_touches_only_eligible_cells_of_fit_rows():
    perturbations, symptoms, outcomes, label_mask, masked_for_soft, evidenced_negative = small_universe()
    weights = np.where(outcomes > 0, 1.0, 0.2) * label_mask
    probability = np.full(outcomes.shape, np.nan)
    strength = np.zeros(outcomes.shape)
    probability[0, 0], strength[0, 0] = 0.7, 0.5  # on a label: ignored
    probability[0, 1], strength[0, 1] = 0.0, 0.5  # an absence on an unlabelled cell: weight rises from 0.2 to 0.5
    probability[1, 1], strength[1, 1] = 0.75, 0.7  # on the grade-C set-aside cell: target 0.75, weight 0.7
    probability[2, 0], strength[2, 0] = 0.6, 0.3  # row 2 is not a fit row
    probability[2, 2], strength[2, 2] = 0.6, 0.3  # a placebo set-aside: ineligible
    eligible = eligible_cells(outcomes, label_mask, evidenced_negative, masked_for_soft)
    new_outcomes, new_weights, summary = apply_soft_priors(outcomes, weights, probability, strength, eligible, fit_rows=np.array([0, 1]), global_strength=1.0)
    assert new_outcomes[0, 0] == 1.0 and new_weights[0, 0] == 1.0, "the label stands"
    assert new_outcomes[0, 1] == 0.0 and new_weights[0, 1] == pytest.approx(0.5)
    assert new_outcomes[1, 1] == 0.75 and new_weights[1, 1] == pytest.approx(0.7)
    assert new_outcomes[2, 0] == 0.0 and new_weights[2, 0] == pytest.approx(0.2), "row 2 is scored, not fitted"
    assert new_weights[2, 2] == 0.0
    assert summary["positive_cells"] == 1 and summary["absence_cells"] == 1 and summary["priors_outside_eligible_cells"] == 1
    assert outcomes[1, 1] == 0.0 and weights[0, 1] == 0.2, "the inputs are not modified"
    halved_outcomes, halved_weights, _ = apply_soft_priors(outcomes, weights, probability, strength, eligible, np.array([0, 1]), global_strength=0.5)
    assert halved_weights[1, 1] == pytest.approx(0.35) and halved_weights[0, 1] == pytest.approx(0.25) and halved_outcomes[1, 1] == 0.75, "the global strength scales weights, not targets"


def test_global_strength_zero_leaves_the_matrices_as_they_were():
    _, _, outcomes, label_mask, masked_for_soft, evidenced_negative = small_universe()
    weights = np.where(outcomes > 0, 1.0, 0.2)
    probability, strength = np.full(outcomes.shape, 0.7), np.full(outcomes.shape, 0.5)
    eligible = eligible_cells(outcomes, label_mask, evidenced_negative, masked_for_soft)
    new_outcomes, new_weights, summary = apply_soft_priors(outcomes, weights, probability, strength, eligible, np.arange(3), global_strength=0.0)
    assert np.array_equal(new_outcomes, outcomes) and np.array_equal(new_weights, weights) and summary["positive_cells"] == 0


def test_the_permutation_keeps_every_row_and_column_count_and_lands_on_eligible_cells():
    generator = np.random.default_rng(1)
    rows, columns = 40, 12
    outcomes = (generator.random((rows, columns)) < 0.08).astype(float)
    eligible = outcomes == 0
    strength = np.where((generator.random((rows, columns)) < 0.15) & eligible, 0.5, 0.0)
    probability = np.where(strength > 0, 0.7, np.nan)
    strength[strength > 0] = generator.choice([0.3, 0.5, 0.7], size=int((strength > 0).sum()))
    fit_rows = np.arange(0, 30)
    new_probability, new_strength, record = permute_soft_priors(probability, strength, eligible, fit_rows, random_seed=3)
    assert record["swaps_made"] > 0
    assert np.array_equal((new_strength > 0).sum(axis=1), (strength > 0).sum(axis=1)), "each perturbation keeps its number of soft pairs"
    assert np.array_equal((new_strength > 0).sum(axis=0), (strength > 0).sum(axis=0)), "each symptom keeps its number of soft pairs"
    assert not ((new_strength > 0) & ~eligible).any(), "no prior lands on a labelled cell"
    assert np.array_equal(new_strength[30:], strength[30:]), "rows outside the fit keep their priors"
    assert sorted(new_strength[new_strength > 0].tolist()) == sorted(strength[strength > 0].tolist()), "the strengths travel with the priors"
    assert np.isnan(new_probability[new_strength == 0]).all() and (new_probability[new_strength > 0] == 0.7).all()
    assert (new_strength[:30] != strength[:30]).any(), "the layout moved"
    again = permute_soft_priors(probability, strength, eligible, fit_rows, random_seed=3)
    assert again[2]["layout_sha256"] == record["layout_sha256"], "the same seed gives the same layout"
    assert permute_soft_priors(probability, strength, eligible, fit_rows, random_seed=4)[2]["layout_sha256"] != record["layout_sha256"]


def test_literature_popularity_counts_distinct_pmids_per_perturbation():
    literature = pd.DataFrame({"perturbation_id": ["G", "G", "G", "D"], "pmid": ["1", "1", "2", "3"], "symptom": ["a", "b", "a", "a"]})
    assert literature_popularity(["D", "G", "Z"], literature).tolist() == [1.0, 2.0, 0.0]
    assert literature_popularity(["D"], pd.DataFrame()).tolist() == [0.0]


def test_the_summary_counts_pairs_by_class_and_where_they_fell():
    perturbations, symptoms, outcomes, label_mask, masked_for_soft, evidenced_negative = small_universe()
    resolved = resolve_pairs(pd.DataFrame([
        {"perturbation_id": "D", "symptom": "anxiety", "evidence_class": "animal", "direction": "positive", "n_distinct": 2, "sources": "IMPC"},  # on a label
        {"perturbation_id": "D", "symptom": "psychosis", "evidence_class": "human_association", "direction": "negative", "n_distinct": 1, "sources": "x"},
        {"perturbation_id": "G", "symptom": "insomnia", "evidence_class": "literature", "direction": "positive", "n_distinct": 1, "sources": "x"},
        {"perturbation_id": "G", "symptom": "insomnia", "evidence_class": "animal", "direction": "negative", "n_distinct": 1, "sources": "x"},
    ]), TIERS)
    probability, strength = soft_prior_matrices(resolved, perturbations, symptoms)
    summary = summarise(resolved, probability, strength, eligible_cells(outcomes, label_mask, evidenced_negative, masked_for_soft))
    assert summary["resolved_pairs"] == 2 and summary["contested_pairs"] == 1 and summary["pairs_on_eligible_cells"] == 1 and summary["pairs_on_labelled_or_set_aside_cells"] == 1
    assert summary["by_class_and_direction"]["animal:positive"] == {"pairs": 1, "mean_strength": pytest.approx(0.3 * np.sqrt(2)), "n_distinct_at_least_two": 1}


def test_trial_arm_rows_enter_only_where_their_interval_supports_them():
    def trial(symptom, value, record, lower, upper):
        row = report("D", symptom, "human_association", "ClinicalTrials-arms", value, record=record)
        row.update({"risk_difference_lower_95": lower, "risk_difference_upper_95": upper})
        return row

    reports = pd.DataFrame([
        trial("anxiety", 1, "NCT1", 0.02, 0.10),       # an excess the interval supports
        trial("anxiety", 1, "NCT2", -0.01, 0.05),      # 3 of 100 against 2 of 100: an excess the interval does not support
        trial("insomnia", 0, "NCT3", -0.04, 0.03),     # an absence with the excess bounded below five in a hundred
        trial("insomnia", 0, "NCT4", -0.10, 0.12),     # an absence that leaves a large excess open
        report("G", "anxiety", "animal", "IMPC", 1, record="MGI:1"),  # another source passes as it is
    ])
    counts = {}
    evidence = soft_evidence_from_reports(reports, "induces", counts)
    rows = {(row.perturbation_id, row.symptom, row.direction): row.n_distinct for row in evidence.itertuples(index=False)}
    assert rows == {("D", "anxiety", "positive"): 1, ("D", "insomnia", "negative"): 1, ("G", "anxiety", "positive"): 1}
    assert counts == {"trial_arm_excesses_admitted": 1, "trial_arm_excesses_not_supported_by_the_interval": 1, "trial_arm_absences_admitted": 1,
                      "trial_arm_absences_with_an_excess_left_open": 1}
    without_interval = reports.drop(columns=["risk_difference_lower_95", "risk_difference_upper_95"])
    counts = {}
    evidence = soft_evidence_from_reports(without_interval, "induces", counts)
    assert len(evidence) == 3 and counts == {"trial_arm_rows_without_an_interval": 4}, "a table without the interval passes as it is and says so"
