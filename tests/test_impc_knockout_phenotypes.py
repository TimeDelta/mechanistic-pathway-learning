"""IMPC knockout phenotypes as cross-species evidence, positive and negative
(mechanistic_pathway_learning/evidence/impc_knockout_phenotypes.py)."""
from __future__ import annotations

import csv
import gzip
import json
import math

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.impc_knockout_phenotypes import (
    DEFAULT_CROSSWALK_PATH,
    EVIDENCE_CLASS,
    EVIDENCE_CODE_NOT_SIGNIFICANT,
    EVIDENCE_CODE_SIGNIFICANT,
    SOURCE_NAME,
    VALIDITY_READ_BY_DEFAULT,
    VALIDITY_READ_WITH_WEAK,
    CrosswalkRow,
    build_reports,
    collect_assertions,
    collect_tests,
    mouse_to_human_symbols,
    procedure_code,
    read_crosswalk,
    read_statistical_results,
)

CROSSWALK = [
    CrosswalkRow("hyperactivity", "induces", ("MP:0001399",), ("OFD",), "strong", "open field"),
    CrosswalkRow("anxiety", "induces", ("MP:0001363", "MP:0002797"), ("OFD", "LDT"), "moderate", "thigmotaxis"),
    CrosswalkRow("anxiety", "relieves", ("MP:0001364",), (), "moderate", "reverse"),
    CrosswalkRow("parkinsonism", "induces", ("MP:0000745",), (), "weak", "tremors"),
    CrosswalkRow("psychosis", "induces", (), (), "endophenotype_only", "never read"),
]
MOUSE_TO_HUMAN = {"MGI:1": "DRD2", "MGI:2": "HTR2A", "MGI:3": "ORPHAN"}
NODES = {"DRD2": "GENE:DRD2", "HTR2A": "GENE:HTR2A"}


def assertion(term: str, name: str, p: float = 1e-5, zygosity: str = "homozygote") -> dict:
    return {"p_value": p, "procedure_stable_id": "IMPC_OFD_001", "procedure_name": "Open Field", "parameter_name": "Distance travelled", "zygosity": zygosity,
            "effect_size": 1.2, "percentage_change": "+40%", "marker_symbol": "Drd2", "allele_symbol": "Drd2<tm1>", "center": "JAX", "mp_term_name": name}


def procedure_record(parameters: int, significant: set[str] = frozenset(), mutants: int = 16, p: float = 0.4):
    from mechanistic_pathway_learning.evidence.impc_knockout_phenotypes import GeneProcedureTests
    record = GeneProcedureTests(parameters_tested={f"P{index}" for index in range(parameters)}, significant_terms=set(significant), smallest_p_value=p,
                                mutant_count=mutants, control_count=400, procedure_names={"Open Field"}, centers={"JAX"})
    return record


def test_procedure_codes_are_the_family_across_centres_and_legacy_numbers_are_unread():
    assert procedure_code("IMPC_OFD_001") == "OFD" and procedure_code("JAX_OFD_001") == "OFD" and procedure_code("HRWLLA_FEA_001") == "FEA"
    assert procedure_code("ESLIM_008_001") is None and procedure_code("M-G-P_007_001") is None and procedure_code("odd") is None


def test_the_committed_crosswalk_reads_and_names_every_symptom_once_per_relation():
    rows = read_crosswalk(DEFAULT_CROSSWALK_PATH)
    pairs = [(row.symptom, row.relation) for row in rows]
    assert len(pairs) == len(set(pairs))
    symptoms = {row.symptom for row in rows}
    crosswalk_symptoms = {row["target_symptom"] for row in csv.DictReader(open("docs/symptom_crosswalk.csv", encoding="utf-8"))}
    assert symptoms == crosswalk_symptoms, "every study symptom has a row saying what its mouse proxy is, or that it has none"
    assert all(row.validity in ("strong", "moderate", "weak", "endophenotype_only", "none") for row in rows)
    for row in rows:
        if row.validity in ("endophenotype_only", "none"):
            assert not row.mp_term_ids and not row.negative_procedure_codes, row.symptom
        for term in row.mp_term_ids:
            assert term.startswith("MP:") and len(term) == 10, term


def test_mouse_genes_map_to_human_genes_one_to_one_only(tmp_path):
    hgnc = tmp_path / "hgnc.txt"
    hgnc.write_text("hgnc_id\tsymbol\tmgd_id\nHGNC:1\tDRD2\tMGI:1\nHGNC:2\tHTR2A\tMGI:2\nHGNC:3\tDUP1\tMGI:9\nHGNC:4\tDUP2\tMGI:9\nHGNC:5\tTWO\tMGI:7|MGI:8\nHGNC:6\tNONE\t\n")
    assert mouse_to_human_symbols(hgnc) == {"MGI:1": "DRD2", "MGI:2": "HTR2A"}


def test_a_mapped_assertion_is_a_positive_with_its_class_and_provenance():
    assertions = {"MGI:1": {"MP:0001399": assertion("MP:0001399", "hyperactivity")}}
    rows, counts = build_reports(assertions, {}, MOUSE_TO_HUMAN, NODES, CROSSWALK)
    assert counts["positives"] == 1 and counts["negatives"] == 0
    row = rows[0]
    assert row["perturbation_id"] == "DRD2" and row["perturbation_type"] == "gene" and row["symptom"] == "hyperactivity" and row["report_value"] == 1
    assert row["evidence_class"] == EVIDENCE_CLASS == "animal" and row["source"] == SOURCE_NAME and row["evidence_code"] == EVIDENCE_CODE_SIGNIFICANT
    assert json.loads(row["perturbation_nodes"]) == [["GENE:DRD2", -1.0, 1.0]] and row["source_record_id"] == "MGI:1"
    assert "homozygote" in row["model_description"] and "Drd2<tm1>" in row["source_record_label"]
    assert row["rubric_cross_species"] == 1.0 and row["rubric_proxy_validity"] == 1.0 and "strong validity" in row["limitations"]


def test_a_term_mapped_to_the_reverse_direction_is_a_relieves_report():
    assertions = {"MGI:2": {"MP:0001364": assertion("MP:0001364", "decreased anxiety-related response")}}
    rows, _ = build_reports(assertions, {}, MOUSE_TO_HUMAN, NODES, CROSSWALK)
    assert [(row["symptom"], row["relation"], row["report_value"]) for row in rows] == [("anxiety", "relieves", 1)]


def test_a_tested_knockout_with_no_call_is_a_negative_at_its_power_and_a_hit_blocks_it():
    tests = {("MGI:1", "OFD", "homozygote"): procedure_record(20, significant={"MP:0001399"}, mutants=16),   # hyperactive: no hyperactivity negative, an anxiety negative
             ("MGI:2", "OFD", "homozygote"): procedure_record(20, mutants=14),                                # nothing: negatives for both
             ("MGI:2", "LDT", "heterozygote"): procedure_record(6, mutants=8)}                                 # light-dark: anxiety only
    assertions = {"MGI:1": {"MP:0001399": assertion("MP:0001399", "hyperactivity")}}
    rows, counts = build_reports(assertions, tests, MOUSE_TO_HUMAN, NODES, CROSSWALK)
    negatives = sorted((row["perturbation_id"], row["symptom"], row["report_id"].rsplit("|", 1)[1]) for row in rows if row["report_value"] == 0)
    assert negatives == [("DRD2", "anxiety", "homozygote"), ("HTR2A", "anxiety", "heterozygote"), ("HTR2A", "anxiety", "homozygote"), ("HTR2A", "hyperactivity", "homozygote")]
    assert counts["negatives_blocked_by_a_hit"] == 1 and counts["positives"] == 1
    drd2_anxiety = [row for row in rows if row["perturbation_id"] == "DRD2" and row["report_value"] == 0][0]
    assert drd2_anxiety["evidence_code"] == EVIDENCE_CODE_NOT_SIGNIFICANT and drd2_anxiety["rubric_log_sample_size"] == pytest.approx(math.log1p(16))
    assert "20 parameters" in drd2_anxiety["model_description"] and "mutants 16" in drd2_anxiety["model_description"]
    assert "absence at the power of 16 mutants" in drd2_anxiety["limitations"]


def test_a_hit_for_the_symptom_in_either_direction_blocks_its_negative():
    tests = {("MGI:2", "OFD", "homozygote"): procedure_record(20, significant={"MP:0001364"})}  # decreased anxiety: a relieves hit
    rows, counts = build_reports({}, tests, MOUSE_TO_HUMAN, NODES, CROSSWALK)
    assert [row["symptom"] for row in rows if row["report_value"] == 0] == ["hyperactivity"] and counts["negatives_blocked_by_a_hit"] == 1


def test_weak_proxies_are_read_only_when_asked_and_endophenotypes_never():
    assertions = {"MGI:1": {"MP:0000745": assertion("MP:0000745", "tremors")}}
    rows, _ = build_reports(assertions, {}, MOUSE_TO_HUMAN, NODES, CROSSWALK, VALIDITY_READ_BY_DEFAULT)
    assert rows == []
    rows, _ = build_reports(assertions, {}, MOUSE_TO_HUMAN, NODES, CROSSWALK, VALIDITY_READ_WITH_WEAK)
    assert [row["symptom"] for row in rows] == ["parkinsonism"] and rows[0]["rubric_proxy_validity"] == 0.3
    assert "psychosis" not in {row["symptom"] for row in rows}


def test_a_weak_hit_still_blocks_a_negative_for_its_symptom():
    """Reading weak proxies as positives is optional; ignoring a weak hit when calling an absence is not."""
    crosswalk = CROSSWALK + [CrosswalkRow("parkinsonism", "induces", ("MP:0001406",), ("OFD",), "weak", "gait")]
    tests = {("MGI:1", "OFD", "homozygote"): procedure_record(20, significant={"MP:0001406"})}
    rows, counts = build_reports({}, tests, MOUSE_TO_HUMAN, NODES, crosswalk, VALIDITY_READ_BY_DEFAULT)
    assert "parkinsonism" not in {row["symptom"] for row in rows}
    assert counts["negatives_blocked_by_a_hit"] == 0, "the weak row's procedure is not read at all under the default validities"


def test_genes_without_an_orthologue_or_a_node_produce_nothing_and_are_counted():
    assertions = {"MGI:3": {"MP:0001399": assertion("MP:0001399", "hyperactivity")}, "MGI:99": {"MP:0001399": assertion("MP:0001399", "hyperactivity")}}
    rows, counts = build_reports(assertions, {}, MOUSE_TO_HUMAN, NODES, CROSSWALK)
    assert rows == [] and counts["human_gene_not_in_graph"] == 1 and counts["no_one_to_one_human_gene"] == 1


def test_the_statistical_results_stream_keeps_only_the_read_families_and_collects_per_zygosity(tmp_path):
    header = ["marker_accession_id", "procedure_stable_id", "procedure_name", "parameter_stable_id", "zygosity", "status", "significant", "p_value",
              "genotype_effect_p_value", "male_mutant_count", "female_mutant_count", "male_control_count", "female_control_count", "mp_term_id", "phenotyping_center"]
    rows = [["MGI:1", "IMPC_OFD_001", "Open Field", "IMPC_OFD_010_001", "homozygote", "Successful", "false", "0.3", "", "7", "8", "200", "210", "", "JAX"],
            ["MGI:1", "IMPC_OFD_001", "Open Field", "IMPC_OFD_011_001", "homozygote", "Successful", "true", "1e-6", "", "7", "8", "200", "210", "MP:0001399", "JAX"],
            ["MGI:1", "IMPC_OFD_001", "Open Field", "IMPC_OFD_012_001", "homozygote", "Failed", "false", "", "", "0", "0", "0", "0", "", "JAX"],
            ["MGI:1", "IMPC_CSD_003", "SHIRPA", "IMPC_CSD_001_001", "homozygote", "Successful", "false", "0.9", "", "7", "8", "200", "210", "", "JAX"],
            ["MGI:1", "IMPC_OFD_001", "Open Field", "IMPC_OFD_010_001", "heterozygote", "Successful", "false", "0.5", "", "4", "4", "200", "210", "", "TCP"]]
    path = tmp_path / "statistical-results-ALL.csv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as results_file:
        writer = csv.writer(results_file)
        writer.writerow(header)
        writer.writerows(rows)
    streamed = list(read_statistical_results(path, {"OFD"}))
    assert len(streamed) == 4, "the SHIRPA row is not of a read family"
    tests = collect_tests(streamed)
    homozygote = tests[("MGI:1", "OFD", "homozygote")]
    assert homozygote.parameters_tested == {"IMPC_OFD_010_001", "IMPC_OFD_011_001"} and homozygote.significant_terms == {"MP:0001399"}
    assert homozygote.smallest_p_value == pytest.approx(1e-6) and homozygote.mutant_count == 15 and homozygote.control_count == 410
    assert tests[("MGI:1", "OFD", "heterozygote")].mutant_count == 8 and tests[("MGI:1", "OFD", "heterozygote")].centers == {"TCP"}


def test_assertions_keep_the_strongest_per_gene_and_term(tmp_path):
    header = ["marker_accession_id", "marker_symbol", "allele_symbol", "phenotyping_center", "zygosity", "procedure_stable_id", "procedure_name", "parameter_name",
              "mp_term_id", "mp_term_name", "p_value", "effect_size", "percentage_change"]
    rows = [["MGI:1", "Drd2", "a", "JAX", "homozygote", "IMPC_OFD_001", "Open Field", "Distance", "MP:0001399", "hyperactivity", "1e-3", "0.8", "+20%"],
            ["MGI:1", "Drd2", "b", "TCP", "heterozygote", "IMPC_OFD_001", "Open Field", "Distance", "MP:0001399", "hyperactivity", "1e-8", "1.5", "+60%"],
            ["MGI:1", "Drd2", "a", "JAX", "homozygote", "IMPC_GRS_001", "Grip", "Grip", "", "", "1e-3", "", ""]]
    path = tmp_path / "genotype-phenotype-assertions-ALL.csv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as assertions_file:
        writer = csv.writer(assertions_file)
        writer.writerow(header)
        writer.writerows(rows)
    assertions = collect_assertions(path)
    assert set(assertions) == {"MGI:1"} and set(assertions["MGI:1"]) == {"MP:0001399"}
    assert assertions["MGI:1"]["MP:0001399"]["allele_symbol"] == "b" and assertions["MGI:1"]["MP:0001399"]["p_value"] == pytest.approx(1e-8)


def test_the_build_script_writes_the_table_the_summary_and_a_stamped_document(tmp_path, monkeypatch):
    import experiments.build_impc_reports as script

    impc = tmp_path / "impc"
    impc.mkdir()
    with gzip.open(impc / "genotype-phenotype-assertions-ALL.csv.gz", "wt", encoding="utf-8", newline="") as assertions_file:
        writer = csv.writer(assertions_file)
        writer.writerow(["marker_accession_id", "marker_symbol", "allele_symbol", "phenotyping_center", "zygosity", "procedure_stable_id", "procedure_name",
                         "parameter_name", "mp_term_id", "mp_term_name", "p_value", "effect_size", "percentage_change"])
        writer.writerow(["MGI:1", "Drd2", "a", "JAX", "homozygote", "IMPC_OFD_001", "Open Field", "Distance", "MP:0001399", "hyperactivity", "1e-6", "1.2", "+40%"])
    with gzip.open(impc / "statistical-results-ALL.csv.gz", "wt", encoding="utf-8", newline="") as results_file:
        writer = csv.writer(results_file)
        writer.writerow(["marker_accession_id", "procedure_stable_id", "procedure_name", "parameter_stable_id", "zygosity", "status", "significant", "p_value",
                         "genotype_effect_p_value", "male_mutant_count", "female_mutant_count", "male_control_count", "female_control_count", "mp_term_id", "phenotyping_center"])
        writer.writerow(["MGI:2", "IMPC_OFD_001", "Open Field", "IMPC_OFD_010_001", "homozygote", "Successful", "false", "0.3", "", "7", "8", "200", "210", "", "JAX"])
    hgnc = tmp_path / "hgnc.txt"
    hgnc.write_text("hgnc_id\tsymbol\tmgd_id\nHGNC:1\tDRD2\tMGI:1\nHGNC:2\tHTR2A\tMGI:2\n")
    graph = tmp_path / "graph"
    graph.mkdir()
    pd.DataFrame({"node_id": ["GENE:DRD2", "GENE:HTR2A", "MAM1c"], "node_type": ["gene", "gene", "metabolite"], "gene_symbol": ["DRD2", "HTR2A", None]}).to_parquet(graph / "nodes.parquet", index=False)
    monkeypatch.setattr("sys.argv", ["x", "--impc-dir", str(impc), "--hgnc", str(hgnc), "--graph-dir", str(graph), "--output-dir", str(tmp_path / "out"),
                                     "--markdown-output", str(tmp_path / "doc.md")])
    assert script.main() == 0
    table = pd.read_parquet(tmp_path / "out" / "impc_reports.parquet")
    assert sorted(zip(table.perturbation_id, table.symptom, table.report_value)) == [("DRD2", "hyperactivity", 1), ("HTR2A", "anxiety", 0), ("HTR2A", "hyperactivity", 0),
                                                                                       ("HTR2A", "psychomotor_retardation", 0)]
    assert set(table.evidence_class) == {"animal"} and "rubric_proxy_validity" in table.columns and str(table.report_value.dtype) == "int64"
    summary = json.loads((tmp_path / "out" / "impc_summary.json").read_text())
    assert summary["counts"]["positives"] == 1 and summary["counts"]["negatives"] == 3 and summary["release"] == "24.0"
    assert "Generated by `experiments/build_impc_reports.py`" in (tmp_path / "doc.md").read_text()
