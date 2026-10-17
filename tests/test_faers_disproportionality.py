"""FAERS disproportionality as evidence, signals and measured absences
(mechanistic_pathway_learning/evidence/faers_disproportionality.py)."""
from __future__ import annotations

import json
import math

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.comparison_statistics import STATISTIC_COLUMNS
from mechanistic_pathway_learning.evidence.evidence_reports import REPORT_COLUMNS
from mechanistic_pathway_learning.evidence.faers_disproportionality import COMPARISON_COLUMNS as FAERS_COMPARISON_COLUMNS
from mechanistic_pathway_learning.evidence.faers_disproportionality import (
    EVIDENCE_CODE,
    SOURCE_NAME,
    Disproportionality,
    build_comparisons,
    build_reports,
    disproportionalities,
    indicated_symptoms_by_drug,
)

PERTURBATION = {"perturbation_id": "CID100000001", "perturbation_label": "examplone", "perturbation_nodes": "[]", "targets": "CHEMBL1:ANTAGONIST"}
TERMS = {"insomnia": "insomnia", "somnolence": "somnolence_or_hypersomnia", "depression": "depressed_mood"}  # crosswalk keys are lowercase; FAERS names are uppercase


def test_the_two_by_two_table_and_the_prr_follow_evans_2001():
    """a=20 of 1,000 drug reports; the reaction in 2,000 of 1,000,000 reports overall: PRR = (20/1000)/(1980/999000)."""
    item = Disproportionality("INSOMNIA", 20, 1000, 2000, 1_000_000)
    a, b, c, d = item.table
    assert (a, b, c, d) == (20, 980, 1980, 997020)
    assert item.prr == pytest.approx((20 / 1000) / (1980 / 999000))
    assert item.ror == pytest.approx((20.5 * 997020.5) / (980.5 * 1980.5))
    low, high = item.ror_interval
    assert low < item.ror < high
    assert item.chi_square > 4


def test_a_signal_needs_the_three_evans_conditions_together():
    strong = Disproportionality("INSOMNIA", 20, 1000, 2000, 1_000_000)
    assert strong.classify() == 1
    too_few = Disproportionality("INSOMNIA", 2, 100, 2000, 1_000_000)  # PRR 10 but two reports
    assert too_few.classify() is None
    weak = Disproportionality("INSOMNIA", 3, 1000, 2000, 1_000_000)  # three reports, PRR 1.5
    assert weak.classify() is None


def test_an_absence_needs_an_upper_bound_below_one_and_enough_reports_naming_the_drug():
    under_reported = Disproportionality("INSOMNIA", 5, 50_000, 20_000, 1_000_000)  # proportion 1e-4 against 2e-2
    assert under_reported.ror_interval[1] < 1.0 and under_reported.classify() == 0
    too_little_exposure = Disproportionality("INSOMNIA", 0, 500, 20_000, 1_000_000)
    assert too_little_exposure.classify() is None
    assert too_little_exposure.classify(negative_minimum_drug_reports=100) == 0


def test_empty_margins_give_no_statistic_rather_than_a_division_error():
    assert Disproportionality("X", 0, 0, 0, 10).prr is None
    assert Disproportionality("X", 3, 10, 3, 1000).prr == math.inf  # reported with the drug and never without it
    assert Disproportionality("X", 0, 10, 0, 1000).prr is None


def test_disproportionalities_cover_reactions_reported_with_the_drug_or_in_the_background():
    items = disproportionalities({"INSOMNIA": 5}, 100, {"INSOMNIA": 500, "SOMNOLENCE": 300}, 100_000)
    assert [(item.reaction, item.with_drug, item.reaction_total) for item in items] == [("INSOMNIA", 5, 500), ("SOMNOLENCE", 0, 300)]


def test_an_unlisted_reaction_enters_at_the_upper_bound_and_can_only_be_an_absence():
    """openFDA lists a drug's 100 most reported reactions; a crosswalk reaction outside them has a count no larger
    than the smallest listed, which can make an under-reporting more conservative and never a signal."""
    items = disproportionalities({"NAUSEA": 5000}, 50_000, {"NAUSEA": 300_000, "INSOMNIA": 20_000}, 2_000_000, unlisted_upper_bound=5000)
    bounded = [item for item in items if item.reaction == "INSOMNIA"][0]
    assert bounded.bounded and bounded.with_drug == 5000
    assert bounded.classify() is None, "at the bound the PRR is 10, but a bound can never be a signal"
    small_bound = disproportionalities({"NAUSEA": 20}, 50_000, {"INSOMNIA": 20_000}, 2_000_000, unlisted_upper_bound=20)[0]
    assert small_bound.classify() == 0, "at most 20 of 50,000 against 1 percent background is an under-reporting"
    rows, _ = build_reports("examplone", PERTURBATION, {"NAUSEA": 20}, 50_000, {"INSOMNIA": 20_000}, 2_000_000, TERMS, set(), "2026-10-16", unlisted_upper_bound=20)
    assert rows and "upper bound" in rows[0]["limitations"] and "at most 20 of 50000" in rows[0]["model_description"]


def test_reports_carry_signals_and_absences_and_set_aside_the_indications():
    reactions = {"INSOMNIA": 400, "SOMNOLENCE": 2, "DEPRESSION": 900}
    background = {"INSOMNIA": 20_000, "SOMNOLENCE": 30_000, "DEPRESSION": 40_000}
    rows, counts = build_reports("examplone", PERTURBATION, reactions, 10_000, background, 2_000_000, TERMS, {"depressed_mood"}, "2026-10-16T00:00:00+00:00")
    by_symptom = {row["symptom"]: row for row in rows}
    assert by_symptom["insomnia"]["report_value"] == 1 and by_symptom["somnolence_or_hypersomnia"]["report_value"] == 0
    assert "depressed_mood" not in by_symptom and counts["indication_confounded"] == 1
    assert counts["signals"] == 1 and counts["absences"] == 1
    row = by_symptom["insomnia"]
    assert row["source"] == SOURCE_NAME and row["evidence_code"] == EVIDENCE_CODE and row["placebo_flag"] is False
    assert "PRR" in row["model_description"] and "400 of 10000 reports" in row["model_description"]
    assert row["rubric_log_sample_size"] == pytest.approx(math.log1p(10_000)) and row["evidence_date"] == "2026-10-16"
    assert "no exposure denominator" in row["limitations"] and "co-medication" in row["limitations"]
    assert set(row) == set(REPORT_COLUMNS) | set(STATISTIC_COLUMNS)
    # the comparison's own statistic beside the verdict (decision 9.48, option 3): 400 observed against 100 expected
    assert row["statistic_name"] == "log_shrunk_observed_to_expected" and row["statistic_estimate"] == pytest.approx(math.log(400.5 / 100.5))
    assert row["statistic_probability_of_excess"] > 0.999 and row["statistic_censored"] == ""


def test_every_comparison_is_kept_with_its_statistic_and_an_unlisted_count_only_when_its_bound_says_something():
    """build_comparisons (decision 9.48, option 3): the neither rows build_reports drops are kept with their statistic;
    an indication is still set aside; an unlisted reaction enters, censored, only when its bound is below the count
    expected, since a bound above it says nothing about the drug."""
    reactions = {"INSOMNIA": 100}
    background = {"INSOMNIA": 20_000, "SOMNOLENCE": 30_000, "DEPRESSION": 40_000}
    rows = build_comparisons("examplone", PERTURBATION, reactions, 10_000, background, 2_000_000, TERMS, {"depressed_mood"}, unlisted_upper_bound=50)
    by_symptom = {row["symptom"]: row for row in rows}
    assert set(by_symptom) == {"insomnia", "somnolence_or_hypersomnia"}, "depression is the indication; nothing else is in the crosswalk"
    assert by_symptom["insomnia"]["verdict"] is None and by_symptom["insomnia"]["expected"] == pytest.approx(100.0)
    assert set(by_symptom["insomnia"]) == set(FAERS_COMPARISON_COLUMNS)
    censored = by_symptom["somnolence_or_hypersomnia"]  # at most 50 against 150 expected
    assert censored["statistic_censored"] == "upper_bound" and censored["statistic_upper_95"] < 0
    above = build_comparisons("examplone", PERTURBATION, reactions, 10_000, background, 2_000_000, TERMS, set(), unlisted_upper_bound=500)
    assert {row["symptom"] for row in above} == {"insomnia"}, "a bound of 500 against 150 and 200 expected says nothing"


def test_a_reaction_that_is_neither_signal_nor_absence_is_counted_and_not_reported():
    # 100 of 10,000 against 20,000 of 2,000,000: PRR 1.0, neither a signal nor an under-reporting
    rows, counts = build_reports("examplone", PERTURBATION, {"INSOMNIA": 100}, 10_000, {"INSOMNIA": 20_000}, 2_000_000, TERMS, set(), "2026-10-16")
    assert rows == [] and counts["neither"] == 1


def test_indicated_symptoms_come_from_relieves_reports(tmp_path):
    reports = pd.DataFrame({"perturbation_id": ["CID1", "CID1", "CID2"], "symptom": ["depressed_mood", "insomnia", "anxiety"],
                            "relation": ["relieves", "induces", "relieves"]})
    path = tmp_path / "evidence_reports.parquet"
    reports.to_parquet(path, index=False)
    assert indicated_symptoms_by_drug(path) == {"CID1": {"depressed_mood"}, "CID2": {"anxiety"}}


def test_the_fetch_script_runs_from_cached_counts_without_a_request(tmp_path, monkeypatch):
    import experiments.fetch_faers_disproportionality as script

    reports = pd.DataFrame({"perturbation_id": ["CID1", "CID1"], "perturbation_type": ["drug", "drug"], "perturbation_label": ["examplone", "examplone"],
                            "perturbation_nodes": ["[]", "[]"], "model_description": ["human; drug label; ChEMBL targets CHEMBL1:ANTAGONIST", "x"],
                            "symptom": ["depressed_mood", "insomnia"], "relation": ["relieves", "induces"]})
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    reports.to_parquet(evidence_dir / "evidence_reports.parquet", index=False)
    crosswalk = tmp_path / "crosswalk.csv"
    crosswalk.write_text("target_symptom,meddra_preferred_terms\ninsomnia,Insomnia\ndepressed_mood,Depression\n")
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "_background.json").write_text(json.dumps({"all_reports": 2_000_000, "reaction_totals": {"INSOMNIA": 20_000, "DEPRESSION": 40_000}}))
    (cache / "examplone.json").write_text(json.dumps({"field": "patient.drug.openfda.generic_name.exact", "drug_total": 10_000,
                                                      "reaction_counts": {"INSOMNIA": 400, "DEPRESSION": 900, "NAUSEA": 50}, "requested_at": "2026-10-16T00:00:00+00:00"}))
    monkeypatch.setattr(script, "query", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("no request expected")))
    monkeypatch.setattr("sys.argv", ["x", "--evidence-dir", str(evidence_dir), "--crosswalk", str(crosswalk), "--cache-dir", str(cache),
                                     "--output-dir", str(tmp_path / "out"), "--markdown-output", str(tmp_path / "doc.md")])
    assert script.main() == 0
    table = pd.read_parquet(tmp_path / "out" / "faers_reports.parquet")
    assert list(table.symptom) == ["insomnia"] and int(table.report_value.iloc[0]) == 1
    summary = json.loads((tmp_path / "out" / "faers_summary.json").read_text())
    assert summary["per_drug"]["examplone"]["indication_confounded"] == 1
    assert "Generated by `experiments/fetch_faers_disproportionality.py`" in (tmp_path / "doc.md").read_text()
