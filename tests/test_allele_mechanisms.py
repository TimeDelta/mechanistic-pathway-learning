"""The allele as the perturbation (evidence/allele_mechanisms.py, decision 9.54 (d)): a typed mechanism carried across
one disease's Orphanet and OMIM records and down to its narrower forms but never through a group, a direct type kept
over an inherited one, conflicts left untyped; each gene's alleles seeded by their class and grouped under the gene."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.allele_mechanisms import (
    GAIN_OF_FUNCTION,
    HYPERMORPH,
    LOSS_OF_FUNCTION,
    NEOMORPH,
    UNTYPED,
    DiseaseMechanism,
    allele_of,
    disease_mechanisms,
    gene_of_perturbation,
    load_orphanet_omim_mappings,
    read_neomorphic_alleles,
)
from mechanistic_pathway_learning.evidence.assemble_evidence_table import (
    aggregate_reports_to_observations,
    assign_alleles,
    disease_cluster_ids_from_reports,
    summarise_alleles,
)
from mechanistic_pathway_learning.evidence.evidence_reports import EvidenceReport, reports_to_dataframe

GAIN = "Disease-causing germline mutation(s) (gain of function) in"
LOSS = "Disease-causing germline mutation(s) (loss of function) in"
UNDIRECTED = "Disease-causing germline mutation(s) in"


def test_a_typed_mechanism_crosses_exact_and_narrower_mappings_and_never_a_group():
    types = {("G", "ORPHA:1"): GAIN, ("G", "ORPHA:2"): UNDIRECTED, ("G", "ORPHA:3"): LOSS, ("G", "ORPHA:9"): UNDIRECTED}
    mappings = {"ORPHA:1": [("OMIM:10", "E")], "ORPHA:2": [("OMIM:10", "NTBT")], "ORPHA:9": [("OMIM:20", "BTNT"), ("OMIM:10", "BTNT")]}
    pairs = [("G", "ORPHA:1"), ("G", "OMIM:10"), ("G", "ORPHA:2"), ("G", "ORPHA:3"), ("G", "ORPHA:9"), ("G", "OMIM:20"), ("H", "OMIM:10")]
    result = disease_mechanisms(types, mappings, pairs)
    assert result[("G", "ORPHA:1")] == DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed")
    assert result[("G", "OMIM:10")] == DiseaseMechanism(GAIN_OF_FUNCTION, "same disease as ORPHA:1"), "one disease, two records"
    assert result[("G", "ORPHA:2")].mechanism == GAIN_OF_FUNCTION and "narrower than OMIM:10" in result[("G", "ORPHA:2")].basis
    assert result[("G", "ORPHA:3")] == DiseaseMechanism(LOSS_OF_FUNCTION, "orphanet typed")
    assert result[("G", "ORPHA:9")].mechanism == UNTYPED and result[("G", "OMIM:20")].mechanism == UNTYPED, "a broader code is a group and passes nothing"
    assert result[("H", "OMIM:10")].mechanism == UNTYPED, "a mechanism belongs to the gene it was typed for"


def test_a_direct_type_wins_and_two_inherited_mechanisms_that_disagree_leave_the_pair_untyped():
    types = {("H", "ORPHA:4"): GAIN, ("H", "ORPHA:5"): LOSS, ("H", "ORPHA:6"): LOSS, ("H", "ORPHA:7"): GAIN}
    mappings = {"ORPHA:4": [("OMIM:30", "E")], "ORPHA:5": [("OMIM:30", "E")], "ORPHA:7": [("OMIM:40", "E")], "ORPHA:6": [("OMIM:40", "NTBT")]}
    result = disease_mechanisms(types, mappings, [("H", "OMIM:30"), ("H", "ORPHA:6")])
    assert result[("H", "OMIM:30")].mechanism == UNTYPED and result[("H", "OMIM:30")].basis.startswith("conflict")
    assert result[("H", "ORPHA:6")] == DiseaseMechanism(LOSS_OF_FUNCTION, "orphanet typed"), "Orphanet's own type for the pair is kept"


def test_product_one_mappings_are_read_with_their_relation_codes(tmp_path):
    xml = """<?xml version="1.0"?><JDBOR><DisorderList count="1"><Disorder id="1"><OrphaCode>399</OrphaCode>
    <ExternalReferenceList count="2"><ExternalReference id="1"><Source>OMIM</Source><Reference>143100</Reference>
    <DisorderMappingRelation id="2"><Name lang="en">E (Exact mapping: the two concepts are equivalent)</Name></DisorderMappingRelation></ExternalReference>
    <ExternalReference id="3"><Source>ICD-10</Source><Reference>G10</Reference></ExternalReference></ExternalReferenceList></Disorder></DisorderList></JDBOR>"""
    path = tmp_path / "product1.xml"
    path.write_text(xml)
    assert load_orphanet_omim_mappings(path) == {"ORPHA:399": [("OMIM:143100", "E")]}


def test_the_committed_neomorphic_alleles_cite_every_edge():
    alleles = read_neomorphic_alleles()
    huntingtin = alleles["HTT"]
    assert huntingtin.allele_id == "HTT|CAG_repeat_expansion" and huntingtin.node_id == "ALLELE:HTT_CAG_repeat_expansion"
    assert {edge.target_symbol for edge in huntingtin.edges} == {"CREBBP", "TP53", "PPARGC1A", "REST"}
    assert all(edge.pmid.isdigit() and len(edge.quote) > 40 and edge.magnitude == 1.0 for edge in huntingtin.edges)
    assert {edge.target_symbol: edge.sign for edge in huntingtin.edges}["REST"] == 1.0, "the wild-type function lost: REST accumulates in the nucleus"


@pytest.mark.parametrize("row,message", [({"effect": "blocks"}, "effect"), ({"pmid": ""}, "PubMed"), ({"target_kind": "complex"}, "target_kind"),
                                         ({"allele_id": "CAG_repeat_expansion"}, "must start with its gene"), ({"magnitude": "2"}, "magnitude")])
def test_a_malformed_allele_row_is_refused(tmp_path, row, message):
    base = {"allele_id": "HTT|CAG_repeat_expansion", "gene": "HTT", "label": "HTT (CAG)", "target_symbol": "CREBBP", "target_kind": "protein",
            "effect": "decreases", "magnitude": "1", "pmid": "11264541", "quote": "Expanded polyglutamine repeats specifically interfere with CBP-activated gene transcription"}
    path = tmp_path / "alleles.csv"
    pd.DataFrame([{**base, **row}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match=message):
        read_neomorphic_alleles(path)


def test_each_mechanism_class_is_its_own_perturbation_with_its_own_seed():
    neomorphs = read_neomorphic_alleles()
    assert allele_of("HTT", DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed"), neomorphs).allele_class == NEOMORPH
    hypermorph = allele_of("SCN8A", DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed"), neomorphs)
    assert (hypermorph.perturbation_id, hypermorph.allele_class, hypermorph.seed_sign) == ("SCN8A|gain_of_function", HYPERMORPH, 1.0)
    for mechanism in (LOSS_OF_FUNCTION, UNTYPED):
        allele = allele_of("HTT", DiseaseMechanism(mechanism, "x"), neomorphs)
        assert (allele.perturbation_id, allele.seed_sign) == ("HTT", -1.0)
    assert gene_of_perturbation("HTT|CAG_repeat_expansion") == "HTT" == gene_of_perturbation("HTT")


def monogenic_report(gene: str, disease: str, symptom: str, index: int) -> EvidenceReport:
    return EvidenceReport(report_id=f"{gene}|{disease}|{symptom}|{index}", perturbation_id=gene, perturbation_type="gene", perturbation_label=gene, symptom=symptom,
                          relation="induces", evidence_class="monogenic", source="HPO-OMIM" if disease.startswith("OMIM") else "HPO-ORPHA", report_value=1,
                          source_record_id=disease, source_record_label=disease, source_term_id="HP:0000001", source_term_label="term", evidence_code="PCS",
                          model_description="human loss-of-function; HPO", frequency=None, frequency_denominator=None, placebo_flag=False, onset=None, sex=None,
                          references="", pubmed_reference_count=0, evidence_date=None, rubric_log_sample_size=0.0, rubric_frequency_known=0.0, rubric_evidence_code_pcs=1.0,
                          rubric_evidence_code_tas=0.0, rubric_evidence_code_iea=0.0, rubric_placebo_controlled=0.0, rubric_curated_synopsis=0.0,
                          association_type=GAIN if disease == "ORPHA:399" else UNDIRECTED, rubric_causal_association=1.0)


def test_alleles_are_seeded_by_class_and_grouped_and_clustered_under_their_gene():
    reports = [monogenic_report("HTT", "ORPHA:399", "irritability", 0), monogenic_report("HTT", "OMIM:143100", "chorea_like", 1),
               monogenic_report("HTT", "OMIM:617435", "intellectual_disability", 2), monogenic_report("SCN8A", "ORPHA:442835", "anxiety", 3),
               monogenic_report("SCN8A", "OMIM:614558", "anxiety", 4)]
    types = {("HTT", "ORPHA:399"): GAIN, ("SCN8A", "ORPHA:442835"): GAIN}
    mappings = {"ORPHA:399": [("OMIM:143100", "E")]}
    mechanisms = disease_mechanisms(types, mappings, {(report.perturbation_id, report.source_record_id) for report in reports})
    nodes = {symbol: f"GENE:{symbol}" for symbol in ("HTT", "SCN8A", "CREBBP", "TP53", "PPARGC1A", "REST")}
    assigned = assign_alleles(reports, mechanisms, read_neomorphic_alleles(), nodes)
    by_id = {report.report_id: report for report in reports}
    assert by_id["HTT|ORPHA:399|irritability|0"].perturbation_id == by_id["HTT|OMIM:143100|chorea_like|1"].perturbation_id == "HTT|CAG_repeat_expansion"
    assert json.loads(by_id["HTT|ORPHA:399|irritability|0"].perturbation_nodes) == [["GENE:CREBBP", -1.0, 1.0], ["GENE:TP53", -1.0, 1.0], ["GENE:PPARGC1A", -1.0, 1.0], ["GENE:REST", 1.0, 1.0]]
    assert by_id["HTT|OMIM:617435|intellectual_disability|2"].perturbation_id == "HTT" and json.loads(by_id["HTT|OMIM:617435|intellectual_disability|2"].perturbation_nodes) == [["GENE:HTT", -1.0, 1.0]]
    assert by_id["SCN8A|ORPHA:442835|anxiety|3"].perturbation_id == "SCN8A|gain_of_function" and "gain-of-function" in by_id["SCN8A|ORPHA:442835|anxiety|3"].model_description
    assert by_id["SCN8A|OMIM:614558|anxiety|4"].perturbation_id == "SCN8A", "an OMIM record no Orphanet disorder maps exactly keeps the default"
    summary = summarise_alleles(assigned, reports)
    assert summary["alleles_by_class"] == {HYPERMORPH: 1, NEOMORPH: 1, UNTYPED: 2} and summary["genes_with_several_alleles_list"] == ["HTT", "SCN8A"]

    table = reports_to_dataframe(reports)
    table["allele_class"] = table.report_id.map(lambda report_id: assigned[report_id][1])
    clusters = disease_cluster_ids_from_reports(table)
    assert set(clusters) == {"HTT", "SCN8A"}, "clusters are of genes, not of alleles"
    observations = aggregate_reports_to_observations(table, "curated_synopsis_provenance", clusters, set())
    assert set(observations.group_id) == {"HTT", "SCN8A"} and set(observations.disease_cluster_id) == {"cluster:HTT", "cluster:SCN8A"}
    assert dict(zip(observations.perturbation_id, observations.allele_class))["HTT|CAG_repeat_expansion"] == NEOMORPH


def test_the_two_seeding_rules_are_not_combined(tmp_path: Path):
    from mechanistic_pathway_learning.evidence.assemble_evidence_table import assemble

    with pytest.raises(ValueError, match="replaces seed_gain_of_function_genes_positive"):
        assemble(tmp_path, tmp_path, tmp_path, tmp_path, None, None, allele_perturbations=True, seed_gain_of_function_genes_positive=True)


def test_literature_on_a_gene_with_several_alleles_is_set_aside_and_knockouts_stay(tmp_path):
    """A literature statement names a gene, not an allele; a knockout is the loss-of-function allele itself."""
    from argparse import Namespace

    from experiments.run_main_model import build_soft_priors, genes_with_several_alleles
    from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
    from tests.test_soft_priors_in_the_trainer import REPOSITORY_ROOT, write_toy_dataset

    assert genes_with_several_alleles(["HTT", "HTT|CAG_repeat_expansion", "SCN8A|gain_of_function", "A+B"]) == {"HTT"}
    paths = write_toy_dataset(tmp_path, num_genes=12)
    records = pd.read_parquet(paths["evidence"] / "evidence_records.parquet")
    allele = records[records.perturbation_id == "G1"].assign(perturbation_id="G1|gain_of_function", perturbation_label="G1 (gain of function)",
                                                              perturbation_nodes=json.dumps([["GENE:G1", 1.0, 1.0]]))
    pd.concat([records, allele], ignore_index=True).to_parquet(paths["evidence"] / "evidence_records.parquet", index=False)
    data = load_experiment_data(paths["graph"], paths["evidence"])
    assert {"G1", "G1|gain_of_function"} <= set(data.perturbation_ids)
    arguments = Namespace(soft_prior_strength=1.0, soft_prior_tiers=REPOSITORY_ROOT / "configs/soft_prior_tiers.csv", evidence_dir=paths["evidence"],
                          soft_prior_record_grades=["C"], soft_prior_reports=[paths["reports"]], soft_prior_tempering=0.5, soft_prior_count_cap=2, soft_prior_probability="tier")
    probability, strength, summary, _, _ = build_soft_priors(data, arguments)
    reports = pd.read_parquet(paths["reports"])
    literature_on_g1 = int(((reports.perturbation_id == "G1") & (reports.evidence_class == "literature")).sum())
    assert literature_on_g1 > 0 and summary["literature_reports_set_aside_for_genes_with_several_alleles"] == literature_on_g1
    records.to_parquet(paths["evidence"] / "evidence_records.parquet", index=False)  # the same table without the allele
    plain = load_experiment_data(paths["graph"], paths["evidence"])
    _, plain_strength, plain_summary, _, _ = build_soft_priors(plain, arguments)
    assert "literature_reports_set_aside_for_genes_with_several_alleles" not in plain_summary, "a table without alleles keeps the summary it had"
    g1, psychosis = data.perturbation_ids.index("G1"), data.symptoms.index("psychosis")
    plain_g1 = plain.perturbation_ids.index("G1")
    literature_only = [column for column in range(len(data.symptoms)) if plain_strength[plain_g1, column] > 0 and strength[g1, column] == 0]
    assert literature_only == [psychosis], "the one prior G1 loses is its literature relation on psychosis"
    others = [column for column in range(len(data.symptoms)) if column != psychosis]
    assert np.array_equal(strength[g1, others], plain_strength[plain_g1, others]), "the knockout's and the records' priors stay with the gene's default allele"
