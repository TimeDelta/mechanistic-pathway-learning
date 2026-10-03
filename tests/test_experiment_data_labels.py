"""Tests: load_experiment_data labels a pair positive only when it has a positive report, and keeps the older table layout as it was."""
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data


def write_toy_graph(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"node_id": ["GENE:A", "GENE:B", "R1"], "node_type": ["gene", "gene", "reaction"], "gene_symbol": ["A", "B", None],
                  "is_currency": [False, False, False], "degree": [1.0, 1.0, 2.0]}).to_parquet(directory / "nodes.parquet", index=False)
    pd.DataFrame({"source_id": ["GENE:A", "GENE:B"], "target_id": ["R1", "R1"], "relation_type": ["catalyzed_by", "catalyzed_by"]}).to_parquet(directory / "edges.parquet", index=False)
    (directory / "relation_types.json").write_text(json.dumps(["catalyzed_by"]))


def evidence_row(gene: str, symptom: str, weight: float, **extra) -> dict:
    return {"perturbation_id": gene, "perturbation_type": "gene", "perturbation_label": gene, "group_id": gene, "disease_cluster_id": "cluster:" + gene,
            "symptom": symptom, "relation": "induces", "evidence_class": "monogenic", "grade": extra.pop("grade", "A"), "weight": weight,
            "label_frequency": None, "in_metabolic_layer": True, "evidence_date": None, "perturbation_nodes": json.dumps([[f"GENE:{gene}", -1.0, 1.0]]), **extra}


def test_experiment_data_never_labels_an_all_negative_pair_positive(tmp_path: Path) -> None:
    write_toy_graph(tmp_path / "graph")
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([
        evidence_row("A", "psychosis", 1.0, positive_report_count=2, negative_report_count=0, report_count=2),
        evidence_row("A", "anxiety", 0.0, positive_report_count=0, negative_report_count=1, report_count=1),
        evidence_row("B", "anxiety", 0.0, positive_report_count=1, negative_report_count=1, report_count=2),
    ]).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    data = load_experiment_data(tmp_path / "graph", evidence_directory)
    assert data.symptoms == ["anxiety", "psychosis"] and data.perturbation_ids == ["A", "B"]
    assert data.outcomes.tolist() == [[0.0, 1.0], [1.0, 0.0]]  # A anxiety has no positive report: outcome 0, not a positive label
    assert data.weights.tolist() == [[0.0, 1.0], [0.0, 0.0]]  # B anxiety: a conflicting pair with weight 0 is still a positive claim


def test_experiment_data_keeps_zero_weight_positives_without_report_counts(tmp_path: Path) -> None:
    write_toy_graph(tmp_path / "graph")
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row("A", "psychosis", 1.0), evidence_row("B", "psychosis", 0.0, grade="C")]).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    data = load_experiment_data(tmp_path / "graph", evidence_directory, label_grades=None)
    assert data.outcomes.tolist() == [[1.0], [1.0]] and data.weights.tolist() == [[1.0], [0.0]]  # the version 0.3 layout: a grade C row is a test positive only when every grade is requested
    default = load_experiment_data(tmp_path / "graph", evidence_directory)
    assert default.perturbation_ids == ["A"] and default.outcomes.tolist() == [[1.0]]  # by default grade C (human association) rows are not labels


def test_label_grades_filter_keeps_grade_b_drug_rows(tmp_path: Path) -> None:
    write_toy_graph(tmp_path / "graph")
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    pd.DataFrame([evidence_row("A", "psychosis", 1.0, grade="A"), evidence_row("B", "psychosis", 0.6, grade="B"), evidence_row("B", "anxiety", 0.0, grade="C")]).to_parquet(evidence_directory / "evidence_records.parquet", index=False)
    data = load_experiment_data(tmp_path / "graph", evidence_directory)
    assert data.symptoms == ["psychosis"] and data.outcomes.tolist() == [[1.0], [1.0]] and data.weights.tolist() == [[1.0], [0.6]]
