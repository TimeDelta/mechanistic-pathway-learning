"""Tests for the report table, its aggregation to observations and the report-level reliability fit (synthetic fixtures only)."""
import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.assemble_evidence_table import (
    aggregate_reports_to_observations,
    disease_cluster_ids_from_reports,
    fit_report_reliability,
    summarize_reports,
)
from mechanistic_pathway_learning.evidence.assign_evidence_grades import assign_evidence_grade, loss_weight_for_record
from mechanistic_pathway_learning.evidence.evidence_reliability_model import (
    fit_dawid_skene,
    fit_weighted_dawid_skene,
    report_count_matrices,
)
from mechanistic_pathway_learning.evidence.evidence_reports import (
    REPORT_COLUMNS,
    EvidenceReport,
    RubricWeightDefaults,
    load_rubric_weight_defaults,
    reports_to_dataframe,
    rubric_weight,
    rubric_weights_for_table,
)
from mechanistic_pathway_learning.evidence.load_drug_label_events import drug_label_reports, load_sider_events
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    load_hpo_annotation_rows,
    load_hpo_is_a_parents_from_obo,
    monogenic_evidence_records,
    monogenic_evidence_reports,
    parse_genes_to_phenotype,
    read_crosswalk_hpo_terms,
)
from tests.test_drug_label_and_target_loaders import write_synthetic_sider
from tests.test_evidence_reliability_and_appraisal import simulate_reports
from tests.test_monogenic_loader_frequency_and_exclusions import write_inputs

HPOA_HEADER = "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration"


def write_hpoa(path: Path, rows: list[str]) -> Path:
    path.write_text("\n".join(["#description: test", HPOA_HEADER, *rows]) + "\n")
    return path


def monogenic_fixture(tmp_path: Path, hpoa_rows: list[str] | None, genes: set[str], publication_dates=None):
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    excluded = {symptom: terms[1] for symptom, terms in crosswalk.items()}
    hpoa_rows_by_key = load_hpo_annotation_rows(write_hpoa(tmp_path / "phenotype.hpoa", hpoa_rows)) if hpoa_rows is not None else None
    rows = parse_genes_to_phenotype(annotations_path)
    reports = monogenic_evidence_reports(rows, roots, parents, genes, excluded, hpoa_rows_by_key, publication_dates)
    for report in reports:
        report.perturbation_nodes = json.dumps([[f"GENE:{report.perturbation_id}", -1.0, 1.0]])
    return reports, rows, roots, parents, excluded


def test_not_qualified_and_excluded_annotations_become_negative_reports(tmp_path: Path) -> None:
    # OTC's disease has a NOT row and an Excluded-frequency row on Dementia (a cognitive_impairment term)
    reports, *_ = monogenic_fixture(tmp_path, [
        "OMIM:311250\tOTC deficiency\tNOT\tHP:0000726\tPMID:1\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
        "OMIM:311250\tOTC deficiency\t\tHP:0000726\tOMIM:311250\tTAS\t\tHP:0040285\t\t\tP\tHPO:a[2013-03-04]",
    ], {"OTC"})
    assert [report.report_value for report in reports] == [0, 0]
    not_report, excluded_report = reports
    assert not_report.frequency is None and excluded_report.frequency == 0.0
    assert not_report.source == "HPO-OMIM" and not_report.evidence_class == "monogenic" and not_report.relation == "induces"
    assert not_report.model_description == "human loss-of-function; OMIM:311250 OTC deficiency"
    assert not_report.perturbation_nodes == json.dumps([["GENE:OTC", -1.0, 1.0]])


def test_evidence_code_references_and_date_join_from_phenotype_hpoa(tmp_path: Path) -> None:
    publication_dates = {"10": date(1999, 5, 1), "20": date(2004, 1, 1)}
    reports, *_ = monogenic_fixture(tmp_path, [
        "OMIM:176000\tAcute intermittent porphyria\t\tHP:0000709\tPMID:10;PMID:20\tPCS\tHP:0003593\t2/7\tMALE\t\tP\tHPO:a[2012-03-04]",
        "OMIM:176000\tAcute intermittent porphyria\t\tHP:0000709\tOMIM:176000\tTAS\t\tHP:0040283\t\t\tP\tHPO:b[2010-05-06]",
        "ORPHA:79276\tAcute intermittent porphyria\t\tHP:0000709\tORPHA:79276\tTAS\t\tHP:0040283\t\t\tP\tORPHA:orphadata[2026-09-02]",
    ], {"HMBS"}, publication_dates)
    by_code = {(report.source, report.evidence_code): report for report in reports}
    pcs, tas, orphanet = by_code[("HPO-OMIM", "PCS")], by_code[("HPO-OMIM", "TAS")], by_code[("HPO-Orphanet", "TAS")]
    assert pcs.pubmed_reference_count == 2 and tas.pubmed_reference_count == 0
    assert pcs.rubric_curated_synopsis == 0.0 and tas.rubric_curated_synopsis == 1.0 and orphanet.rubric_curated_synopsis == 1.0
    assert pcs.evidence_date == "1999-05-01" and tas.evidence_date == "2010-05-06" and orphanet.evidence_date is None
    assert pcs.frequency == pytest.approx(2 / 7) and pcs.frequency_denominator == 7 and pcs.onset == "HP:0003593" and pcs.sex == "MALE"
    assert pcs.references == "PMID:10;PMID:20" and pcs.rubric_evidence_code_pcs == 1.0 and tas.rubric_evidence_code_tas == 1.0
    assert pcs.source_record_label == "Acute intermittent porphyria" and pcs.source_term_label == "Psychosis"
    assert "n = 7 patients" in pcs.limitations and "undated" in orphanet.limitations and "no PubMed reference" in tas.limitations
    assert [report.report_id.rsplit("|", 1)[1] for report in reports] == ["1", "2", "1"]  # ordinal within the (disease, term) group


def test_unjoined_rows_fall_back_to_genes_to_phenotype(tmp_path: Path) -> None:
    reports, *_ = monogenic_fixture(tmp_path, None, {"HMBS", "OTC", "SHAREDGENE"})
    assert reports and all(report.evidence_code == "unjoined" for report in reports)
    by_pair_and_disease = {(report.perturbation_id, report.symptom, report.source_record_id): report for report in reports}
    assert by_pair_and_disease[("OTC", "cognitive_impairment", "OMIM:311250")].report_value == 1
    assert by_pair_and_disease[("OTC", "cognitive_impairment", "OMIM:311250")].frequency == pytest.approx(0.6)
    assert by_pair_and_disease[("SHAREDGENE", "cognitive_impairment", "ORPHA:2")].report_value == 0  # genes_to_phenotype frequency HP:0040285
    assert by_pair_and_disease[("HMBS", "psychosis", "OMIM:176000")].frequency is None and by_pair_and_disease[("HMBS", "psychosis", "OMIM:176000")].evidence_date is None
    assert all(report.report_id.endswith("|0") and report.source_record_label == report.source_record_id for report in reports)


def test_default_weighting_reproduces_grade_weights(tmp_path: Path) -> None:
    reports, rows, roots, parents, excluded = monogenic_fixture(tmp_path, None, {"HMBS", "OTC"})
    table = reports_to_dataframe(reports)
    table["rubric_weight"] = rubric_weights_for_table(table, RubricWeightDefaults())
    observations = aggregate_reports_to_observations(table, "curated_synopsis_provenance", disease_cluster_ids_from_reports(table), {"GENE:HMBS"})
    records = {(record.perturbation_identifier, record.symptom_identifier): record for record in monogenic_evidence_records(rows, roots, parents, {"HMBS", "OTC"}, excluded)}
    assert len(observations) == len(records) == 2
    for row in observations.itertuples():
        record = records[(row.perturbation_id, row.symptom)]
        assert row.grade == assign_evidence_grade(record) and row.weight == loss_weight_for_record(record)
        assert (row.label_frequency == record.max_annotation_frequency) or (pd.isna(row.label_frequency) and record.max_annotation_frequency is None)
    by_gene = observations.set_index("perturbation_id")
    assert by_gene.loc["HMBS", "weight"] == pytest.approx(0.68) and by_gene.loc["OTC", "weight"] == 1.0
    assert bool(by_gene.loc["HMBS", "in_metabolic_layer"]) and not bool(by_gene.loc["OTC", "in_metabolic_layer"])
    assert by_gene.loc["HMBS", "omim_entry_count"] == 1 and by_gene.loc["HMBS", "orpha_entry_count"] == 1 and by_gene.loc["HMBS", "report_count"] == 2
    # SIDER: one event per preferred term; the insomnia frequency is the mean of the treatment midpoints with the placebo row excluded
    sider_directory = tmp_path / "sider"
    write_synthetic_sider(sider_directory)
    events = [event for event in load_sider_events(sider_directory, Path("docs/symptom_crosswalk.csv")) if event.target_symptom]
    drug_reports = [report for event in events for report in drug_label_reports(event, "human; drug label; ChEMBL targets CHEMBL2093872:POSITIVE ALLOSTERIC MODULATOR", json.dumps([["GENE:GABRA1", 0.5, 1.0]]))]
    drug_table = reports_to_dataframe(drug_reports)
    drug_table["rubric_weight"] = rubric_weights_for_table(drug_table, RubricWeightDefaults())
    drug_observations = aggregate_reports_to_observations(drug_table, "curated_synopsis_provenance", {}, set(), {"CID100003016": "CHEMBL2093872:POSITIVE ALLOSTERIC MODULATOR"})
    by_pair = drug_observations.set_index(["symptom", "relation"])
    insomnia = by_pair.loc[("insomnia", "induces")]
    assert insomnia.label_frequency == pytest.approx(0.055) and insomnia.weight == pytest.approx(0.6 * 0.55) and insomnia.grade == "B"
    assert insomnia.group_id == "CHEMBL2093872" and insomnia.source == "SIDER 4.1; targets CHEMBL2093872:POSITIVE ALLOSTERIC MODULATOR"
    assert by_pair.loc[("psychosis", "induces")].weight == 0.6 and pd.isna(by_pair.loc[("psychosis", "induces")].label_frequency)
    assert by_pair.loc[("anxiety", "relieves")].grade == "B" and drug_table[drug_table.relation == "relieves"].source.unique().tolist() == ["SIDER-indication"]
    insomnia_report = drug_table[drug_table.symptom == "insomnia"].iloc[0]
    assert insomnia_report.evidence_code == "placebo_controlled_frequency" and bool(insomnia_report.placebo_flag) and insomnia_report.source_term_id == "C0917801"


def test_all_negative_pair_is_kept_with_zero_weight(tmp_path: Path) -> None:
    reports, *_ = monogenic_fixture(tmp_path, [
        "OMIM:311250\tOTC deficiency\tNOT\tHP:0000726\tPMID:1\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
    ], {"OTC"})
    table = reports_to_dataframe(reports)
    table["rubric_weight"] = rubric_weights_for_table(table, RubricWeightDefaults())
    observations = aggregate_reports_to_observations(table, "curated_synopsis_provenance", {}, set())
    row = observations.iloc[0]
    assert (row.perturbation_id, row.symptom) == ("OTC", "cognitive_impairment")
    assert row.positive_report_count == 0 and row.negative_report_count == 1 and row.weight == 0.0 and row.grade == "C"
    assert row.disease_identifiers == "OMIM:311250" and row.omim_entry_count == 1 and pd.isna(row.label_frequency)


def synthetic_report_table(rows: list[tuple[str, str, str, int, float]]) -> pd.DataFrame:
    """(perturbation_id, symptom, source, report_value, rubric_weight) -> a minimal report table."""
    table = pd.DataFrame(rows, columns=["perturbation_id", "symptom", "source", "report_value", "rubric_weight"])
    table["relation"] = "induces"
    return table


def test_conflicting_pair_posterior_lies_strictly_between() -> None:
    sources = ["A", "B", "C"]
    rng = np.random.default_rng(0)
    background = [(f"p{i}", "s", source, int(rng.random() < 0.7), 1.0) for i in range(60) for source in sources]
    conflicting = background + [("x", "s", "A", 1, 1.0), ("x", "s", "B", 0, 1.0)]
    all_positive = background + [("x", "s", "A", 1, 1.0), ("x", "s", "B", 1, 1.0)]
    all_negative = background + [("x", "s", "A", 0, 1.0), ("x", "s", "B", 0, 1.0)]
    items = [(f"p{i}", "s", "induces") for i in range(60)] + [("x", "s", "induces")]
    posteriors = []
    for rows in (conflicting, all_positive, all_negative):
        positive, negative, coverage = report_count_matrices(synthetic_report_table(rows), items, sources)
        assert coverage[-1].tolist() == [1.0, 1.0, 0.0]  # C never reports on x, so it does not cover it
        posteriors.append(fit_weighted_dawid_skene(positive, negative, sources, items).posterior[-1])
    posterior_conflicting, posterior_positive, posterior_negative = posteriors
    assert posterior_negative < posterior_conflicting < posterior_positive


def test_weighted_fit_reduces_to_dawid_skene_on_unit_weights() -> None:
    _, reports, coverage = simulate_reports(3000, [0.95, 0.80, 0.55], [0.95, 0.85, 0.55], prevalence=0.3)
    parameters, posterior = fit_dawid_skene(reports, coverage, ["monogenic", "label", "comention"])
    fit = fit_weighted_dawid_skene(reports * coverage, (1 - reports) * coverage, ["monogenic", "label", "comention"])
    assert np.allclose(fit.posterior, posterior) and np.allclose(fit.sensitivity, parameters.sensitivity) and np.allclose(fit.specificity, parameters.specificity)
    assert fit.prevalence == pytest.approx(parameters.prevalence) and fit.log_likelihood_trace == parameters.log_likelihood_trace
    assert not fit.weakly_identified and fit.num_iterations == len(fit.log_likelihood_trace)


def make_report(evidence_code: str, denominator: int | None = None) -> EvidenceReport:
    return EvidenceReport("id", "G", "gene", "G", "s", "induces", "monogenic", "HPO-OMIM", 1, "OMIM:1", "d", "HP:1", "t", evidence_code, "m",
                          0.5, denominator, False, None, None, "", 0, None, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def test_rubric_weights_order_evidence_codes_and_sample_sizes(tmp_path: Path) -> None:
    defaults = RubricWeightDefaults()
    assert rubric_weight(make_report("PCS"), defaults) > rubric_weight(make_report("TAS"), defaults) > rubric_weight(make_report("IEA"), defaults)
    assert rubric_weight(make_report("PCS", 1), defaults) < rubric_weight(make_report("PCS", 20), defaults) == 1.0
    assert rubric_weight(make_report("PCS", 1), defaults) == pytest.approx(0.5 + 0.5 * np.log1p(1) / np.log1p(20))
    assert rubric_weight(make_report("placebo_controlled_frequency"), defaults) > rubric_weight(make_report("label_text"), defaults)
    assert rubric_weight(make_report("unjoined"), defaults) == 0.7 == rubric_weight(make_report("ZZZ"), defaults)
    overrides = tmp_path / "rubric.json"
    overrides.write_text(json.dumps({"evidence_code_factors": {"TAS": 0.3}, "minimum_rubric_weight": 0.2}))
    changed = load_rubric_weight_defaults(overrides)
    assert changed.evidence_code_factors["TAS"] == 0.3 and changed.evidence_code_factors["PCS"] == 1.0 and changed.minimum_rubric_weight == 0.2
    assert changed.implicit_negative_weight == 1.0 and rubric_weight(make_report("TAS"), changed) == pytest.approx(0.3)
    (tmp_path / "bad.json").write_text(json.dumps({"evidence_code_factor": {"TAS": 0.3}}))
    with pytest.raises(ValueError):
        load_rubric_weight_defaults(tmp_path / "bad.json")
    table = reports_to_dataframe([make_report("PCS", 5), make_report("IEA")])
    assert list(table.columns) == list(REPORT_COLUMNS)
    assert np.allclose(rubric_weights_for_table(table, defaults), [rubric_weight(make_report("PCS", 5), defaults), 0.6])


def test_weakly_identified_flag() -> None:
    two_sources = synthetic_report_table([(f"p{i}", "s", source, 1, 1.0) for i in range(5) for source in ("A", "B")])
    items = [(f"p{i}", "s", "induces") for i in range(5)]
    positive, negative, _ = report_count_matrices(two_sources, items, ["A", "B"])
    assert fit_weighted_dawid_skene(positive, negative, ["A", "B"], items).weakly_identified
    three_on_one = synthetic_report_table([(f"p{i}", "s", source, 1, 1.0) for i in range(5) for source in ("A", "B")] + [("p0", "s", "C", 1, 1.0)])
    positive, negative, _ = report_count_matrices(three_on_one, items, ["A", "B", "C"])
    assert not fit_weighted_dawid_skene(positive, negative, ["A", "B", "C"], items).weakly_identified


def test_implicit_negative_only_inside_coverage() -> None:
    table = synthetic_report_table([("g1", "s1", "A", 1, 0.8), ("g1", "s2", "B", 1, 0.5), ("g2", "s1", "A", 1, 0.8), ("g1", "s1", "A", 0, 0.6)])
    items = [("g1", "s1", "induces"), ("g1", "s2", "induces"), ("g2", "s1", "induces")]
    positive, negative, coverage = report_count_matrices(table, items, ["A", "B"], implicit_negative_weight=0.4)
    assert positive.tolist() == [[0.8, 0.0], [0.0, 0.5], [0.8, 0.0]]
    assert negative.tolist() == [[0.6, 0.4], [0.4, 0.0], [0.0, 0.0]]  # g2 has no B report: B contributes nothing; g1 s2 is silent in A: implicit negative
    assert coverage.tolist() == [[1.0, 1.0], [1.0, 1.0], [1.0, 0.0]]


def test_report_ids_are_unique_and_deterministic(tmp_path: Path) -> None:
    hpoa_rows = [
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
        "OMIM:176000\tAIP\t\tHP:0000709\tOMIM:176000\tTAS\t\t\t\t\tP\tHPO:b[2010-05-06]",
        "ORPHA:79276\tAIP\t\tHP:0000709\tORPHA:79276\tTAS\t\t\t\t\tP\tORPHA:orphadata[2026-09-02]",
        "OMIM:1\tX\t\tHP:0000709\tPMID:2\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
    ]
    first, *_ = monogenic_fixture(tmp_path / "a", hpoa_rows, {"HMBS", "OTC", "NOTINGRAPH"})
    second, *_ = monogenic_fixture(tmp_path / "b", hpoa_rows, {"HMBS", "OTC", "NOTINGRAPH"})
    identifiers = [report.report_id for report in first]
    assert len(identifiers) == len(set(identifiers)) == 5 and identifiers == [report.report_id for report in second]
    assert identifiers[0] == "HPO-OMIM|HMBS|OMIM:176000|HP:0000709|psychosis|1"
    sider_directory = tmp_path / "sider"
    write_synthetic_sider(sider_directory)
    drug_identifiers = [report.report_id for event in load_sider_events(sider_directory, Path("docs/symptom_crosswalk.csv")) for report in drug_label_reports(event, "m", "[]")]
    assert len(drug_identifiers) == len(set(drug_identifiers)) == 3


def test_fit_report_reliability_and_summary_over_a_small_table(tmp_path: Path) -> None:
    reports, *_ = monogenic_fixture(tmp_path, [
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
        "ORPHA:79276\tAIP\t\tHP:0000709\tORPHA:79276\tTAS\t\t\t\t\tP\tORPHA:orphadata[2026-09-02]",
        "OMIM:311250\tOTC\tNOT\tHP:0000726\tPMID:1\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
        "OMIM:311250\tOTC\t\tHP:0000726\tPMID:3\tPCS\t\t1/4\t\t\tP\tHPO:a[2012-03-04]",
    ], {"HMBS", "OTC"})
    table = reports_to_dataframe(reports)
    defaults = RubricWeightDefaults()
    table["rubric_weight"] = rubric_weights_for_table(table, defaults)
    observations = aggregate_reports_to_observations(table, "curated_synopsis_provenance", {}, set())
    fit = fit_report_reliability(table, observations, defaults)
    assert fit.source_names == ["HPO-OMIM", "HPO-Orphanet"] and fit.weakly_identified and len(fit.posterior) == len(observations) == 2
    assert fit.item_keys == [("HMBS", "psychosis", "induces"), ("OTC", "cognitive_impairment", "induces")]
    assert fit.per_source_summary["HPO-OMIM"]["items_covered"] == 2 and fit.per_source_summary["HPO-Orphanet"]["items_covered"] == 1
    assert fit.per_source_summary["HPO-Orphanet"]["implicit_negative_cells"] == 0 and fit.per_source_summary["HPO-OMIM"]["implicit_negative_cells"] == 0
    summary = summarize_reports(table, fit, defaults, observations, reports_dropped_unknown_provenance=0)
    assert summary["reports"]["rows"] == 4 and summary["reports"]["pairs_conflicting"] == 1 and summary["reports"]["pairs_all_negative"] == 0
    assert summary["reports"]["by_source_and_value"] == {"HPO-OMIM|0": 1, "HPO-OMIM|1": 2, "HPO-Orphanet|1": 1}
    assert summary["reliability"]["weakly_identified"] and summary["reliability"]["posterior_quantiles"]["conflicting"]["pairs"] == 1
    assert summary["reliability"]["rubric_weight_defaults"]["implicit_negative_weight"] == 1.0
    json.dumps(summary)  # serializable
