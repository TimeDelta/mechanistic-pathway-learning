"""Trial outcome measures as a therapeutic-direction source (evidence/clinical_trial_outcomes.py, decision 9.60)."""
from __future__ import annotations

import math

import pytest

from mechanistic_pathway_learning.evidence.clinical_trial_outcomes import (
    EVIDENCE_CLASS,
    NO_EFFECT_MARGIN,
    OUTCOME_COLUMNS,
    SOURCE_NAME,
    ScaleRule,
    build_reports,
    compare_outcome,
    first_measurements,
    group_counts,
    read_scale_crosswalk,
    rows_to_dataframe,
    rule_for,
    select_outcome_groups,
    sponsor_comparison,
    standard_deviation,
    verdict_of,
    welch_interval,
)
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

import re

MADRS = ScaleRule(re.compile("MADRS", re.IGNORECASE), "MADRS", "depressed_mood", True, "")
MMSE = ScaleRule(re.compile("MMSE", re.IGNORECASE), "MMSE", "cognitive_impairment", False, "")


def outcome(title="Change From Baseline in MADRS Total Score", param_type="MEAN", dispersion="Standard Deviation", groups=None, values=None, spreads=None,
            counts=None, analyses=None, limits=None):
    groups = groups or [("OG000", "Placebo"), ("OG001", "Drugx 20 mg")]
    values = values or ["-8.0", "-12.0"] * len(groups)
    spreads = spreads or ["10.0"] * len(groups)
    counts = counts or ["100"] * len(groups)
    measurements = []
    for index, (group_id, _) in enumerate(groups):
        measurement = {"groupId": group_id, "value": values[index], "spread": spreads[index] if spreads[index] is not None else None}
        if limits:
            measurement["lowerLimit"], measurement["upperLimit"] = limits[index]
        measurements.append(measurement)
    return {"type": "PRIMARY", "title": title, "paramType": param_type, "dispersionType": dispersion, "unitOfMeasure": "score on a scale",
            "groups": [{"id": group_id, "title": group_title, "description": ""} for group_id, group_title in groups],
            "denoms": [{"units": "Participants", "counts": [{"groupId": group_id, "value": count} for (group_id, _), count in zip(groups, counts)]}],
            "classes": [{"categories": [{"measurements": measurements}]}], "analyses": analyses or []}


def study(outcomes, nct="NCT00000001", title="A trial of drugx"):
    return {"protocolSection": {"identificationModule": {"nctId": nct, "briefTitle": title}, "statusModule": {"overallStatus": "COMPLETED"},
                                "designModule": {"designInfo": {"allocation": "RANDOMIZED", "maskingInfo": {"masking": "DOUBLE"}}}},
            "resultsSection": {"outcomeMeasuresModule": {"outcomeMeasures": outcomes}}}


def test_group_selection_names_the_drug_against_placebo_and_leaves_combinations_out():
    measure = outcome(groups=[("OG000", "Placebo"), ("OG001", "Drugx 20 mg"), ("OG002", "Drugx + Othery"), ("OG003", "Othery alone")])
    selection = select_outcome_groups(measure, "drugx", other_study_drugs=["othery", "drugx"])
    assert selection.drug_groups == ("OG001",) and selection.placebo_groups == ("OG000",) and selection.excluded == ("OG002",)
    assert select_outcome_groups(outcome(groups=[("OG000", "Drugx plus placebo background"), ("OG001", "Drugx")]), "drugx").placebo_groups == ()


def test_counts_measurements_and_dispersions_are_read_from_the_module_shapes():
    measure = outcome(dispersion="Standard Error")
    assert group_counts(measure) == {"OG000": 100, "OG001": 100}
    assert first_measurements(measure)["OG001"]["value"] == "-12.0"
    assert standard_deviation({"spread": "2.0"}, "Standard Deviation", 100) == 2.0
    assert standard_deviation({"spread": "1.0"}, "Standard Error", 100) == pytest.approx(10.0)
    assert standard_deviation({"lowerLimit": "-14", "upperLimit": "-10"}, "95% Confidence Interval", 100) == pytest.approx(2 / 1.959964 * 10)
    assert standard_deviation({"spread": "1.0"}, "Inter-Quartile Range", 100) is None
    assert group_counts({"denoms": [{"units": "Eyes", "counts": []}]}) == {}


def test_the_welch_interval_and_the_sponsor_analysis_are_oriented_as_drug_minus_placebo():
    difference, lower, upper = welch_interval(-12.0, 10.0, 100, -8.0, 10.0, 100)
    assert difference == -4.0 and lower == pytest.approx(-4 - 1.959964 * math.sqrt(2)) and upper == pytest.approx(-4 + 1.959964 * math.sqrt(2))
    analysis = {"groupIds": ["OG000", "OG001"], "paramType": "Mean Difference (Final Values)", "paramValue": "4.0", "ciLowerLimit": "1.0", "ciUpperLimit": "7.0", "ciPctValue": "95"}
    assert sponsor_comparison({"analyses": [analysis]}, "OG001", "OG000") == (-4.0, -7.0, -1.0), "the sponsor states placebo minus drug; flipped"
    assert sponsor_comparison({"analyses": [analysis]}, "OG000", "OG001") == (4.0, 1.0, 7.0)
    assert sponsor_comparison({"analyses": [{**analysis, "groupIds": ["OG000", "OG002"]}]}, "OG001", "OG000") is None
    assert sponsor_comparison({"analyses": [{**analysis, "paramType": "Hazard Ratio"}]}, "OG001", "OG000") is None
    assert sponsor_comparison({"analyses": [{**analysis, "ciPctValue": "90"}]}, "OG001", "OG000") is None


def test_verdicts_follow_the_interval_and_the_scale_direction():
    assert verdict_of((-0.4, -0.6, -0.2), (-4.0, -6.0, -2.0), higher_is_worse=True) == "relief"
    assert verdict_of((0.4, 0.2, 0.6), (4.0, 2.0, 6.0), higher_is_worse=True) == "worsening"
    assert verdict_of((0.4, 0.2, 0.6), (4.0, 2.0, 6.0), higher_is_worse=False) == "relief", "a higher score is better on MMSE"
    assert verdict_of((0.0, -0.1, 0.15), (0.0, -1.0, 1.5), higher_is_worse=True) == "no_effect"
    assert verdict_of((0.0, -0.3, 0.1), (0.0, -3.0, 1.0), higher_is_worse=True) == "undecided", f"the lower bound is beyond {NO_EFFECT_MARGIN}"
    assert verdict_of(None, (0.0, -1.0, 1.0), higher_is_worse=True) == "undecided", "no standard deviation, so no no-effect reading"


def test_a_comparison_needs_one_drug_arm_one_placebo_arm_and_a_mean():
    selection = select_outcome_groups(outcome(), "drugx")
    comparison = compare_outcome(outcome(), selection, MADRS)
    assert comparison["verdict"] == "relief" and comparison["comparison_source"] == "arm means, Welch interval"
    assert comparison["standardised_difference"] == pytest.approx(-0.4) and comparison["drug_n"] == 100
    assert compare_outcome(outcome(param_type="MEDIAN"), selection, MADRS) is None
    two_doses = outcome(groups=[("OG000", "Placebo"), ("OG001", "Drugx 10 mg"), ("OG002", "Drugx 20 mg")], values=["-8", "-10", "-12"], spreads=["10", "10", "10"], counts=["100", "100", "100"])
    assert compare_outcome(two_doses, select_outcome_groups(two_doses, "drugx"), MADRS) is None, "which dose arm to read is not this module's choice"
    no_spread = outcome(spreads=[None, None])
    assert compare_outcome(no_spread, selection, MADRS) is None
    sponsor = outcome(spreads=[None, None], analyses=[{"groupIds": ["OG001", "OG000"], "paramType": "Mean Difference (Net)", "paramValue": "-4", "ciLowerLimit": "-6", "ciUpperLimit": "-2"}])
    comparison = compare_outcome(sponsor, selection, MADRS)
    assert comparison["comparison_source"] == "sponsor analysis" and comparison["verdict"] == "relief" and comparison["standardised_difference"] is None


def test_reports_carry_the_relation_value_and_outcome_columns_and_make_a_table():
    rules = [MADRS, MMSE]
    studies = [study([outcome(), outcome(title="Change in MMSE", values=["0.5", "0.7"], spreads=["3", "3"], counts=["200", "200"]),
                      outcome(title="PANSS Total"), outcome(title="MADRS at week 2", groups=[("OG000", "Drugx"), ("OG001", "Othery")])])]
    rows, counts = build_reports(studies, "drugx", {"perturbation_id": "CID1", "perturbation_label": "drugx", "targets": "T1", "perturbation_nodes": "[]"}, rules, other_study_drugs=["othery"])
    assert counts == {"studies": 1, "outcome_measures": 4, "crosswalked": 3, "placebo_controlled": 2, "compared": 2, "relief": 1, "worsening": 0, "no_effect": 0, "undecided": 1}
    by_symptom = {row["symptom"]: row for row in rows}
    assert by_symptom["depressed_mood"]["relation"] == "relieves" and by_symptom["depressed_mood"]["report_value"] == 1
    assert by_symptom["cognitive_impairment"]["report_value"] == -1 and by_symptom["cognitive_impairment"]["verdict"] == "undecided"
    assert by_symptom["depressed_mood"]["source"] == SOURCE_NAME and by_symptom["depressed_mood"]["evidence_class"] == EVIDENCE_CLASS
    assert by_symptom["depressed_mood"]["rubric_placebo_controlled"] == 1.0 and by_symptom["depressed_mood"]["references"] == "NCT:NCT00000001"
    assert by_symptom["depressed_mood"]["report_id"] == f"{SOURCE_NAME}|CID1|NCT00000001|MADRS|depressed_mood|relieves|1"
    table = rows_to_dataframe(rows)
    assert list(table.columns) == [*REPORT_COLUMNS, *OUTCOME_COLUMNS] and str(table.drug_n.dtype) == "Int64" and table.report_value.tolist() == [1, -1]
    assert rows_to_dataframe([]).shape[0] == 0


def test_the_committed_scale_crosswalk_reads_and_orders_its_rules():
    rules = read_scale_crosswalk()
    assert len(rules) >= 25
    assert rule_for("Change From Baseline in ADHD-RS-IV Hyperactivity/Impulsivity Subscale Score", rules).construct == "hyperactivity"
    assert rule_for("Change in PANSS Positive Subscale", rules).construct == "psychosis" and rule_for("PANSS Total Score", rules) is None
    assert rule_for("Mean Change in CAARS ADHD Symptoms Total Score", rules) is None, "a total that fuses inattention with hyperactivity is not read"
    facit, fatigue = rule_for("FACIT-Fatigue Scale", rules), rule_for("Fatigue Severity Scale", rules)
    assert facit.higher_is_worse is False and fatigue.higher_is_worse is True
    assert rule_for("Changes in Sexual Functioning Questionnaire", rules).higher_is_worse is False and rule_for("Arizona Sexual Experience Scale", rules).higher_is_worse is True
    assert rule_for("Mini-Mental State Examination", rules).higher_is_worse is False and rule_for("ADAS-Cog", rules).higher_is_worse is True
