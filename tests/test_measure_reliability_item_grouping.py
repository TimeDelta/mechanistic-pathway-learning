"""The item-grouping measurement of decision 9.55 on a toy report table (experiments/measure_reliability_item_grouping.py)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.measure_reliability_item_grouping import (  # noqa: E402
    GROUPINGS,
    coverage_family_labels,
    fit_under_grouping,
    item_groups_for,
    render,
)
from mechanistic_pathway_learning.evidence.evidence_reports import RubricWeightDefaults  # noqa: E402


def toy_tables():
    """Six items: three gene pairs the HPO sources agree on, three mouse-only pairs IMPC tested, one of them a call."""
    observations = pd.DataFrame({
        "perturbation_id": ["G1", "G1", "G2", "G3", "G3", "G4"], "symptom": ["anxiety", "insomnia", "anxiety", "anxiety", "insomnia", "anxiety"],
        "relation": ["induces"] * 6, "perturbation_type": ["gene"] * 6, "evidence_class": ["monogenic", "monogenic", "monogenic", "animal", "animal", "animal"],
        "grade": ["A", "A", "A", "D", "D", "D"],
    })
    rows = []
    for perturbation, symptom in [("G1", "anxiety"), ("G1", "insomnia"), ("G2", "anxiety")]:
        for source in ("HPO-OMIM", "HPO-Orphanet"):
            rows.append({"perturbation_id": perturbation, "symptom": symptom, "relation": "induces", "source": source, "report_value": 1, "source_record_id": f"{source}:{perturbation}", "rubric_weight": 1.0})
    for perturbation, symptom, value in [("G3", "anxiety", 0), ("G3", "insomnia", 0), ("G4", "anxiety", 1)]:
        rows.append({"perturbation_id": perturbation, "symptom": symptom, "relation": "induces", "source": "IMPC-knockout", "report_value": value, "source_record_id": f"MGI:{perturbation}", "rubric_weight": 1.0})
    return observations, pd.DataFrame(rows)


def test_coverage_family_labels_join_the_families_of_the_covering_sources():
    coverage = np.array([[1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 0, 0]])
    assert coverage_family_labels(coverage, ["HPO-OMIM", "SIDER-label", "IMPC-knockout"]) == ["genetics+label", "impc", "genetics+impc", "none"]


def test_each_grouping_partitions_the_items_as_named():
    observations, reports = toy_tables()
    coverage = np.array([[1, 1, 0]] * 3 + [[0, 0, 1]] * 3)
    sources = ["HPO-OMIM", "HPO-Orphanet", "IMPC-knockout"]
    assert item_groups_for("evidence_class", observations, coverage, sources) == ["monogenic"] * 3 + ["animal"] * 3
    assert item_groups_for("perturbation_type", observations, coverage, sources) == ["gene"] * 6
    assert item_groups_for("coverage_family", observations, coverage, sources) == ["genetics"] * 3 + ["impc"] * 3
    assert item_groups_for("pooled", observations, coverage, sources) == ["all"] * 6


def test_the_fits_report_their_groups_sources_and_posteriors_and_the_table_renders():
    observations, reports = toy_tables()
    results = [fit_under_grouping(grouping, observations, reports, RubricWeightDefaults()) for grouping in GROUPINGS]
    by_name = {result["grouping"]: result for result in results}
    assert by_name["evidence_class"]["groups"] == 2 and by_name["perturbation_type"]["groups"] == 1 and by_name["coverage_family"]["groups"] == 2 and by_name["pooled"]["groups"] == 1
    for result in results:
        assert set(result["per_source"]) == {"HPO-OMIM", "HPO-Orphanet", "IMPC-knockout"}
        assert set(result["posterior_by_grade"]) == {"A", "D"} and result["posterior_by_grade"]["A"]["items"] == 3
        assert sum(entry["items"] for entry in result["prevalence_by_group"].values()) == 6
    assert by_name["evidence_class"]["posterior_by_grade"]["A"]["median"] > 0.9, "two agreeing HPO sources make the label pairs near-certain"
    text = render(results)
    assert "IMPC-knockout" in text and "posterior median by grade" in text and "prevalence by group" in text
    for grouping in GROUPINGS:
        assert grouping in text
