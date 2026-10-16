"""Placebo-controlled trial arms as evidence reports, positive and negative
(mechanistic_pathway_learning/evidence/clinical_trial_adverse_events.py; the owner, 15 October 2026: "it doesn't have
to be positive only if we find a way to incorporate negatives")."""
from __future__ import annotations

import gzip
import json
import math

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.clinical_trial_adverse_events import (
    ARM_COLUMNS,
    EVIDENCE_CODE,
    SOURCE_NAME,
    build_reports,
    name_pattern,
    perturbations_from_reports,
    read_lowest_level_to_preferred_term,
    newcombe_risk_difference_interval,
    rows_to_dataframe,
    select_arms,
    study_descriptors,
    symptom_of_term,
    term_comparisons,
)
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS

PERTURBATION = {"perturbation_id": "CID100000001", "perturbation_label": "examplone", "perturbation_nodes": "[[\"GENE:HTR2A\", -1.0, 1.0]]", "targets": "CHEMBL1:ANTAGONIST"}
TERMS = {"insomnia": "insomnia", "somnolence": "somnolence_or_hypersomnia", "anxiety": "anxiety"}


def study(nct: str, groups: list[tuple[str, str]], other: list[dict], serious: list[dict] | None = None, threshold: str = "5",
          allocation: str = "RANDOMIZED", masking: str = "DOUBLE", status: str = "COMPLETED", posted: str = "2021-03-04") -> dict:
    return {
        "protocolSection": {
            "identificationModule": {"nctId": nct, "briefTitle": f"Trial {nct}"},
            "statusModule": {"overallStatus": status, "resultsFirstPostDateStruct": {"date": posted}},
            "designModule": {"studyType": "INTERVENTIONAL", "designInfo": {"allocation": allocation, "maskingInfo": {"masking": masking}}},
        },
        "resultsSection": {"adverseEventsModule": {"frequencyThreshold": threshold, "eventGroups": [{"id": identifier, "title": title} for identifier, title in groups],
                                                   "otherEvents": other, "seriousEvents": serious or []}},
    }


def event(term: str, stats: dict[str, tuple[int, int]], organ: str = "Psychiatric disorders") -> dict:
    return {"term": term, "organSystem": organ, "sourceVocabulary": "MedDRA 20.0",
            "stats": [{"groupId": group, "numAffected": affected, "numAtRisk": at_risk} for group, (affected, at_risk) in stats.items()]}


def test_the_drug_name_matches_as_a_whole_word():
    pattern = name_pattern("Lu")
    assert pattern.search("lu 10 mg") and not pattern.search("flu shot") and not pattern.search("lumateperone")
    assert name_pattern("doxepin").search("doxepin 6 mg") and not name_pattern("doxepin").search("doxepinone")


def test_arms_are_the_groups_naming_the_drug_against_the_groups_naming_placebo():
    trial = study("NCT1", [("EG0", "Examplone 10 mg"), ("EG1", "Examplone 20 mg"), ("EG2", "Placebo")], [])
    selection = select_arms(trial, "examplone")
    assert selection.drug_group_ids == ("EG0", "EG1") and selection.placebo_group_ids == ("EG2",) and selection.usable


def test_a_group_naming_both_the_drug_and_placebo_belongs_to_neither_side():
    """An adjunctive trial of another drug on top of the study drug: both arms name it, neither is its arm."""
    trial = study("NCT2", [("EG0", "Otherdrug + examplone"), ("EG1", "Placebo + examplone")], [])
    selection = select_arms(trial, "examplone")
    assert not selection.usable and selection.drug_group_ids == ("EG0",) and selection.placebo_group_ids == ()


def test_an_arm_naming_a_second_study_drug_is_left_out_of_both_sides():
    """The held-out split groups drugs by target; an arm naming two study drugs would carry one arm's events into
    two records, one of which may be held out while the other trains. It is excluded, on either side."""
    trial = study("NCT2b", [("EG0", "Examplone + partnerol"), ("EG1", "Examplone alone"), ("EG2", "Placebo"), ("EG3", "Placebo + partnerol")], [])
    selection = select_arms(trial, "examplone", other_study_drugs={"partnerol", "examplone", "unrelatedin"})
    assert selection.drug_group_ids == ("EG1",) and selection.placebo_group_ids == ("EG2",)
    assert selection.excluded_for_another_study_drug == ("Examplone + partnerol", "Placebo + partnerol")
    without_the_rule = select_arms(trial, "examplone")
    assert without_the_rule.drug_group_ids == ("EG0", "EG1"), "the rule is what removes the combination arm"


def test_term_comparisons_pool_the_dose_arms_and_skip_terms_missing_an_arm():
    trial = study("NCT3", [("EG0", "Examplone 10 mg"), ("EG1", "Examplone 20 mg"), ("EG2", "Placebo")],
                  [event("Insomnia", {"EG0": (3, 50), "EG1": (5, 50), "EG2": (2, 100)}), event("Headache", {"EG0": (1, 50)})],
                  serious=[event("Anxiety", {"EG0": (0, 50), "EG1": (1, 50), "EG2": (0, 100)})])
    comparisons = term_comparisons(trial, select_arms(trial, "examplone"))
    assert [(c.term, c.seriousness, c.drug_affected, c.drug_at_risk, c.placebo_affected, c.placebo_at_risk) for c in comparisons] == [
        ("Insomnia", "other", 8, 100, 2, 100), ("Anxiety", "serious", 1, 100, 0, 100)]
    assert comparisons[0].drug_rate == pytest.approx(0.08) and comparisons[0].placebo_rate == pytest.approx(0.02)


def test_terms_map_through_the_preferred_terms_and_through_the_lowest_level_terms(tmp_path):
    sider = tmp_path / "sider"
    sider.mkdir()
    with gzip.open(sider / "meddra_all_se.tsv.gz", "wt", encoding="utf-8") as se_file:
        se_file.write("CID1\tCID2\tC0001\tPT\tC0001\tInsomnia\n")
        se_file.write("CID1\tCID2\tC0002\tLLT\tC0001\tSleeplessness\n")
        se_file.write("CID1\tCID2\tC0003\tLLT\tC0999\tOrphan lowest term\n")
    lowest = read_lowest_level_to_preferred_term(sider)
    assert lowest == {"sleeplessness": "insomnia"}
    assert symptom_of_term("INSOMNIA", TERMS, lowest) == "insomnia"
    assert symptom_of_term("Sleeplessness", TERMS, lowest) == "insomnia"
    assert symptom_of_term("Headache", TERMS, lowest) is None
    assert read_lowest_level_to_preferred_term(tmp_path / "absent") == {}


def test_reports_are_positive_above_placebo_and_negative_only_above_the_exposure_floor():
    trial = study("NCT4", [("EG0", "Examplone"), ("EG1", "Placebo")],
                  [event("Insomnia", {"EG0": (10, 200), "EG1": (4, 200)}),       # excess: positive
                   event("Somnolence", {"EG0": (5, 200), "EG1": (6, 200)}),      # no excess, exposure 200: negative
                   event("Anxiety", {"EG0": (0, 40), "EG1": (1, 40)})])          # no excess, exposure 40: no report
    rows, counts = build_reports([trial], "examplone", PERTURBATION, TERMS, {}, minimum_exposure_for_absence=100)
    assert [(row["symptom"], row["report_value"]) for row in rows] == [("insomnia", 1), ("somnolence_or_hypersomnia", 0)]
    assert counts == {"studies": 1, "studies_with_usable_arms": 1, "arms_excluded_for_another_study_drug": 0, "terms": 3, "terms_on_a_symptom": 3,
                      "reports_positive": 1, "reports_negative": 1, "non_excess_below_exposure_floor": 1}


def test_a_report_carries_the_arms_the_threshold_and_the_design_in_its_text():
    trial = study("NCT5", [("EG0", "Examplone 5 mg"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (9, 120), "EG1": (3, 118)})],
                  threshold="3", allocation="RANDOMIZED", masking="QUADRUPLE", status="COMPLETED")
    rows, _ = build_reports([trial], "examplone", PERTURBATION, TERMS, {})
    row = rows[0]
    assert row["source"] == SOURCE_NAME and row["evidence_code"] == EVIDENCE_CODE and row["placebo_flag"] is True
    assert row["frequency"] == pytest.approx(9 / 120) and row["frequency_denominator"] == 120
    assert "9/120" in row["model_description"] and "3/118" in row["model_description"] and "CHEMBL1:ANTAGONIST" in row["model_description"]
    assert "threshold 3 percent" in row["limitations"] and "QUADRUPLE" in row["limitations"] and "Examplone 5 mg" in row["limitations"]
    assert row["references"] == "NCT:NCT5" and row["evidence_date"] == "2021-03-04" and row["evidence_date_source"] == "release"
    assert row["rubric_placebo_controlled"] == 1.0 and row["rubric_log_sample_size"] == pytest.approx(math.log1p(120))
    assert row["perturbation_id"] == "CID100000001" and row["perturbation_nodes"] == PERTURBATION["perturbation_nodes"]
    assert row["report_id"].startswith(f"{SOURCE_NAME}|CID100000001|NCT5|insomnia|other|")


def test_an_unlisted_term_and_a_study_without_usable_arms_produce_nothing():
    trial = study("NCT6", [("EG0", "Otherdrug"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (9, 120), "EG1": (3, 118)})])
    rows, counts = build_reports([trial], "examplone", PERTURBATION, TERMS, {})
    assert rows == [] and counts["studies"] == 1 and counts["studies_with_usable_arms"] == 0


def test_a_combination_arm_produces_no_report_for_either_drug():
    trial = study("NCT6b", [("EG0", "Examplone + partnerol"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (9, 120), "EG1": (3, 118)})])
    for drug in ("examplone", "partnerol"):
        rows, counts = build_reports([trial], drug, PERTURBATION, TERMS, {}, other_study_drugs={"examplone", "partnerol"})
        assert rows == [] and counts["arms_excluded_for_another_study_drug"] == 1, drug


def test_study_descriptors_read_the_threshold_and_fall_back_to_five_percent():
    trial = study("NCT7", [], [], threshold="")
    descriptors = study_descriptors(trial)
    assert descriptors["frequency_threshold"] == 5.0 and descriptors["allocation"] == "RANDOMIZED" and descriptors["results_first_posted"] == "2021-03-04"


def test_rows_become_a_table_with_the_report_schema_and_dtypes():
    trial = study("NCT8", [("EG0", "Examplone"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (10, 200), "EG1": (4, 200)})])
    rows, _ = build_reports([trial], "examplone", PERTURBATION, TERMS, {})
    table = rows_to_dataframe(rows)
    assert list(table.columns) == [*REPORT_COLUMNS, *ARM_COLUMNS]
    assert str(table.report_value.dtype) == "int64" and str(table.frequency_denominator.dtype) == "Int64" and table.placebo_flag.dtype == bool
    assert str(table.drug_at_risk.dtype) == "Int64" and str(table.risk_difference_lower_95.dtype) == "float64"
    assert table.drug_affected.tolist() == [10] and table.placebo_at_risk.tolist() == [200] and table.risk_difference.iloc[0] == pytest.approx(0.03)
    assert rows_to_dataframe([]).empty and list(rows_to_dataframe([]).columns) == [*REPORT_COLUMNS, *ARM_COLUMNS]


def test_the_newcombe_interval_has_width_at_zero_events_and_excludes_zero_for_a_clear_excess():
    difference, lower, upper = newcombe_risk_difference_interval(0, 100, 0, 100)
    assert difference == 0.0 and lower < 0 < upper and upper == pytest.approx(0.0370, abs=0.001), "0 of 100 against 0 of 100 leaves about 3.7 percent open either way"
    difference, lower, upper = newcombe_risk_difference_interval(30, 200, 10, 200)
    assert difference == pytest.approx(0.10) and lower > 0, "30 of 200 against 10 of 200 is an excess the interval supports"
    difference, lower, upper = newcombe_risk_difference_interval(3, 100, 2, 100)
    assert lower < 0 < upper, "3 of 100 against 2 of 100 is an excess the interval does not support"
    rows, _ = build_reports([study("NCT9", [("EG0", "Examplone"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (3, 100), "EG1": (2, 100)})])], "examplone", PERTURBATION, TERMS, {})
    assert rows[0]["report_value"] == 1 and rows[0]["risk_difference_lower_95"] < 0, "the row is positive by the source's rule and says so with its interval"


def test_perturbations_are_read_from_the_assembled_report_table(tmp_path):
    reports = pd.DataFrame({"perturbation_id": ["CID1", "CID1", "HTR2A"], "perturbation_type": ["drug", "drug", "gene"],
                            "perturbation_label": ["Examplone", "Examplone", "HTR2A"], "perturbation_nodes": ["[[\"GENE:X\", -1.0, 1.0]]"] * 3,
                            "model_description": ["human; drug label; ChEMBL targets CHEMBL1:ANTAGONIST; more", "x", "human loss-of-function"]})
    path = tmp_path / "evidence_reports.parquet"
    reports.to_parquet(path, index=False)
    identities = perturbations_from_reports(path)
    assert set(identities) == {"examplone"}
    assert identities["examplone"]["targets"] == "CHEMBL1:ANTAGONIST" and identities["examplone"]["perturbation_id"] == "CID1"


def test_the_fetch_script_builds_the_same_rows_from_a_cached_file(tmp_path, monkeypatch):
    """The script's cache is the raw record; a run with cached files makes no request."""
    import experiments.fetch_clinicaltrials_adverse_events as script

    reports = pd.DataFrame({"perturbation_id": ["CID1"], "perturbation_type": ["drug"], "perturbation_label": ["examplone"],
                            "perturbation_nodes": ["[]"], "model_description": ["human; drug label; ChEMBL targets CHEMBL1:ANTAGONIST"]})
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    reports.to_parquet(evidence_dir / "evidence_reports.parquet", index=False)
    crosswalk = tmp_path / "crosswalk.csv"
    crosswalk.write_text("target_symptom,meddra_preferred_terms\ninsomnia,Insomnia\n")
    cache = tmp_path / "cache"
    cache.mkdir()
    trial = study("NCT9", [("EG0", "Examplone"), ("EG1", "Placebo")], [event("Insomnia", {"EG0": (10, 200), "EG1": (4, 200)})])
    (cache / "examplone.json").write_text(json.dumps({"drug": "examplone", "studies": [trial]}))
    monkeypatch.setattr(script, "fetch_studies", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("no request expected")))
    monkeypatch.setattr("sys.argv", ["x", "--evidence-dir", str(evidence_dir), "--crosswalk", str(crosswalk), "--sider-dir", str(tmp_path / "nosider"),
                                     "--cache-dir", str(cache), "--output-dir", str(tmp_path / "out"), "--markdown-output", str(tmp_path / "doc.md")])
    assert script.main() == 0
    table = pd.read_parquet(tmp_path / "out" / "clinical_trial_reports.parquet")
    assert len(table) == 1 and table.report_value.iloc[0] == 1 and table.symptom.iloc[0] == "insomnia"
    assert "Generated by `experiments/fetch_clinicaltrials_adverse_events.py`" in (tmp_path / "doc.md").read_text()
    summary = json.loads((tmp_path / "out" / "clinical_trial_summary.json").read_text())
    assert summary["totals"]["reports_positive"] == 1
