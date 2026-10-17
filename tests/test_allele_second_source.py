"""Gene2Phenotype as the second source of a disease's mechanism (evidence/allele_mechanisms.py, decision 9.54 (d)): read
only where Orphanet types nothing, carried to an Orphanet disorder through exact and narrower mappings and never a
group, typed only at definitive, strong or moderate confidence, a dominant-negative product read as a loss of function,
G2P's own disagreements left untyped, and a disagreement with Orphanet recorded without overruling it."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.allele_mechanisms import (
    DEFAULT_GENE2PHENOTYPE_PATH,
    GAIN_OF_FUNCTION,
    GENE2PHENOTYPE_COLUMNS,
    HYPERMORPH,
    LOSS_OF_FUNCTION,
    UNTYPED,
    DiseaseMechanism,
    disease_mechanisms,
    first_source_unless_absent,
    gene2phenotype_mechanisms,
    load_gene2phenotype,
    read_neomorphic_alleles,
)
from mechanistic_pathway_learning.evidence.assemble_evidence_table import assign_alleles, summarise_alleles
from tests.test_allele_mechanisms import GAIN, monogenic_report

HEADER = ["g2p id", "gene symbol", "gene mim", "hgnc id", "previous gene symbols", "disease name", "disease mim", "disease MONDO", "allelic requirement",
          "cross cutting modifier", "confidence", "variant consequence", "variant types", "molecular mechanism", "molecular mechanism support",
          "molecular mechanism categorisation", "molecular mechanism evidence", "phenotypes", "publications", "additional mined publications", "panel",
          "comments", "date of last review", "review"]


def g2p_row(g2p_id: str, gene: str, mim: str, mechanism: str, confidence: str = "definitive", previous: str = "", allelic: str = "monoallelic_autosomal") -> dict:
    row = {column: "" for column in HEADER}
    row.update({"g2p id": g2p_id, "gene symbol": gene, "previous gene symbols": previous, "disease mim": mim, "confidence": confidence,
                "molecular mechanism": mechanism, "allelic requirement": allelic})
    return row


def write_g2p(tmp_path: Path, rows: list[dict]) -> Path:
    path = tmp_path / "G2P_all_panels.csv"
    pd.DataFrame(rows, columns=HEADER).to_csv(path, index=False)
    return path


def test_records_are_read_and_those_without_an_omim_disease_are_counted(tmp_path):
    path = write_g2p(tmp_path, [g2p_row("G2P1", "A", "100", "loss of function", previous="OLD1; OLD2"), g2p_row("G2P2", "B", "", "gain of function")])
    records, counts = load_gene2phenotype(path)
    assert [(record.g2p_id, record.disease, record.previous_symbols) for record in records] == [("G2P1", "OMIM:100", ("OLD1", "OLD2"))]
    assert counts == {"records_with_an_omim_disease": 1, "records_without_an_omim_disease": 1}
    pd.DataFrame([{"g2p id": "x"}]).to_csv(tmp_path / "other.csv", index=False)
    with pytest.raises(ValueError, match="lacks the G2P columns"):
        load_gene2phenotype(tmp_path / "other.csv")


def test_the_mechanism_reaches_orphanet_through_exact_and_narrower_mappings_only(tmp_path):
    records, _ = load_gene2phenotype(write_g2p(tmp_path, [g2p_row("G2P1", "A", "100", "gain of function"), g2p_row("G2P2", "A", "200", "dominant negative", "strong")]))
    mappings = {"ORPHA:1": [("OMIM:100", "E")], "ORPHA:2": [("OMIM:100", "NTBT")], "ORPHA:3": [("OMIM:100", "BTNT")], "ORPHA:4": [("OMIM:200", "E")]}
    pairs = [("A", "OMIM:100"), ("A", "ORPHA:1"), ("A", "ORPHA:2"), ("A", "ORPHA:3"), ("A", "ORPHA:4"), ("B", "OMIM:100")]
    result = gene2phenotype_mechanisms(records, mappings, pairs)
    assert result[("A", "OMIM:100")] == DiseaseMechanism(GAIN_OF_FUNCTION, "gene2phenotype G2P1 (gain of function, definitive)")
    assert result[("A", "ORPHA:1")].mechanism == GAIN_OF_FUNCTION and "for OMIM:100, the same disease" in result[("A", "ORPHA:1")].basis
    assert result[("A", "ORPHA:2")].mechanism == GAIN_OF_FUNCTION and "narrower form" in result[("A", "ORPHA:2")].basis
    assert ("A", "ORPHA:3") not in result, "a broader Orphanet code is a group and receives nothing"
    assert result[("A", "ORPHA:4")].mechanism == LOSS_OF_FUNCTION and "dominant negative" in result[("A", "ORPHA:4")].basis, "a dominant-negative product lowers the gene's net activity"
    assert ("B", "OMIM:100") not in result, "a record speaks for its own gene"


def test_only_curated_confidence_types_a_pair_and_undetermined_types_nothing(tmp_path):
    records, _ = load_gene2phenotype(write_g2p(tmp_path, [
        g2p_row("G2P1", "A", "100", "gain of function", "limited"), g2p_row("G2P2", "A", "101", "gain of function", "disputed"),
        g2p_row("G2P3", "A", "102", "gain of function", "moderate"), g2p_row("G2P4", "A", "103", "undetermined non-loss-of-function")]))
    result = gene2phenotype_mechanisms(records, {}, [("A", "OMIM:100"), ("A", "OMIM:101"), ("A", "OMIM:102"), ("A", "OMIM:103")])
    assert ("A", "OMIM:100") not in result and ("A", "OMIM:101") not in result
    assert result[("A", "OMIM:102")].mechanism == GAIN_OF_FUNCTION
    assert result[("A", "OMIM:103")] == DiseaseMechanism(UNTYPED, "untyped; gene2phenotype G2P4 (undetermined non-loss-of-function, definitive)")


def test_disagreeing_records_leave_the_pair_untyped_and_a_previous_symbol_matches(tmp_path):
    records, _ = load_gene2phenotype(write_g2p(tmp_path, [
        g2p_row("G2P1", "A", "100", "loss of function", allelic="biallelic_autosomal"), g2p_row("G2P2", "A", "100", "gain of function"),
        g2p_row("G2P3", "NEW", "300", "loss of function", previous="OLD; A")]))
    result = gene2phenotype_mechanisms(records, {}, [("A", "OMIM:100"), ("OLD", "OMIM:300"), ("A", "OMIM:300")])
    assert result[("A", "OMIM:100")].mechanism == UNTYPED and result[("A", "OMIM:100")].basis.startswith("conflict: ")
    assert result[("OLD", "OMIM:300")].mechanism == LOSS_OF_FUNCTION, "a gene renamed since the HPO release is still found"
    assert ("A", "OMIM:300") not in result, "a previous symbol that is another G2P gene's current symbol is not used"


def test_the_second_source_fills_an_absence_and_never_overrules_the_first():
    first = {("A", "OMIM:1"): DiseaseMechanism(UNTYPED, "untyped"), ("A", "OMIM:2"): DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed"),
             ("A", "OMIM:3"): DiseaseMechanism(UNTYPED, "conflict: gain_of_function (x); loss_of_function (y)"), ("A", "OMIM:4"): DiseaseMechanism(LOSS_OF_FUNCTION, "orphanet typed")}
    second = {("A", "OMIM:1"): DiseaseMechanism(GAIN_OF_FUNCTION, "gene2phenotype G2P1"), ("A", "OMIM:2"): DiseaseMechanism(LOSS_OF_FUNCTION, "gene2phenotype G2P2"),
              ("A", "OMIM:3"): DiseaseMechanism(LOSS_OF_FUNCTION, "gene2phenotype G2P3"), ("A", "OMIM:4"): DiseaseMechanism(LOSS_OF_FUNCTION, "gene2phenotype G2P4")}
    pairs = list(first) + [("A", "OMIM:5")]
    result = first_source_unless_absent(first, second, pairs)
    assert result[("A", "OMIM:1")] == DiseaseMechanism(GAIN_OF_FUNCTION, "gene2phenotype G2P1")
    assert result[("A", "OMIM:2")] == DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed; disputed by gene2phenotype G2P2")
    assert result[("A", "OMIM:3")].mechanism == UNTYPED and result[("A", "OMIM:3")].basis.endswith("; gene2phenotype G2P3"), "a conflict in the first source stays one"
    assert result[("A", "OMIM:4")] == DiseaseMechanism(LOSS_OF_FUNCTION, "orphanet typed"), "agreement adds nothing to the basis"
    assert result[("A", "OMIM:5")] == DiseaseMechanism(UNTYPED, "untyped")


def test_a_gain_of_function_only_g2p_states_moves_the_reports_to_the_gain_allele_and_is_counted(tmp_path):
    reports = [monogenic_report("SCN8A", "OMIM:614558", "anxiety", 0), monogenic_report("SCN8A", "ORPHA:442835", "anxiety", 1),
               monogenic_report("KCNC1", "OMIM:616187", "anxiety", 2)]
    pairs = {(report.perturbation_id, report.source_record_id) for report in reports}
    orphanet = disease_mechanisms({("SCN8A", "ORPHA:442835"): GAIN, ("KCNC1", "ORPHA:5"): "Disease-causing germline mutation(s) (loss of function) in"},
                                  {"ORPHA:5": [("OMIM:616187", "E")]}, pairs)
    records, _ = load_gene2phenotype(write_g2p(tmp_path, [g2p_row("G2P10", "SCN8A", "614558", "gain of function"), g2p_row("G2P11", "KCNC1", "616187", "gain of function")]))
    mechanisms = first_source_unless_absent(orphanet, gene2phenotype_mechanisms(records, {}, pairs), pairs)
    assigned = assign_alleles(reports, mechanisms, read_neomorphic_alleles(), {"SCN8A": "GENE:SCN8A", "KCNC1": "GENE:KCNC1"})
    by_id = {report.report_id: report for report in reports}
    assert by_id["SCN8A|OMIM:614558|anxiety|0"].perturbation_id == "SCN8A|gain_of_function" == by_id["SCN8A|ORPHA:442835|anxiety|1"].perturbation_id
    assert json.loads(by_id["SCN8A|OMIM:614558|anxiety|0"].perturbation_nodes) == [["GENE:SCN8A", 1.0, 1.0]]
    assert by_id["KCNC1|OMIM:616187|anxiety|2"].perturbation_id == "KCNC1", "Orphanet's loss of function stands against G2P"
    summary = summarise_alleles(assigned, reports)
    assert summary["pairs_typed_by_the_second_source_only"] == 1 and summary["reports_by_mechanism_basis"]["gene2phenotype"] == 1
    assert summary["pairs_disputed_by_the_second_source"] == ["KCNC1 OMIM:616187: same disease as ORPHA:5; disputed by gene2phenotype G2P11 (gain of function, definitive)"]
    assert summary["alleles_by_class"] == {HYPERMORPH: 1, LOSS_OF_FUNCTION: 1}


@pytest.mark.skipif(not DEFAULT_GENE2PHENOTYPE_PATH.exists(), reason="the G2P download is fetched by scripts/fetch_raw_sources.sh")
def test_the_pinned_download_has_the_columns_and_mechanisms_the_reader_expects():
    records, counts = load_gene2phenotype(DEFAULT_GENE2PHENOTYPE_PATH)
    assert set(GENE2PHENOTYPE_COLUMNS) <= set(pd.read_csv(DEFAULT_GENE2PHENOTYPE_PATH, nrows=0).columns)
    assert counts["records_with_an_omim_disease"] > 2000
    assert {record.mechanism for record in records} <= {"loss of function", "dominant negative", "gain of function", "undetermined", "undetermined non-loss-of-function"}


def test_the_committed_adjudications_cite_every_row_and_settle_the_disputes_read_on_16_october():
    from mechanistic_pathway_learning.evidence.allele_mechanisms import read_mechanism_adjudications

    adjudications = read_mechanism_adjudications()
    assert adjudications[("CDC73", "ORPHA:99880")].mechanism == LOSS_OF_FUNCTION, "a tumour suppressor: inactivating germline mutations"
    assert adjudications[("PRPS1", "ORPHA:411536")].mechanism == GAIN_OF_FUNCTION, "catalytic superactivity: the normal enzyme overexpressed"
    assert {pair: mechanism.mechanism for pair, mechanism in adjudications.items() if pair[0] in ("ARHGAP31", "KCNC1", "SLC40A1")} == {
        ("ARHGAP31", "ORPHA:974"): GAIN_OF_FUNCTION, ("ARHGAP31", "OMIM:100300"): GAIN_OF_FUNCTION, ("KCNC1", "OMIM:616187"): LOSS_OF_FUNCTION,
        ("SLC40A1", "OMIM:606069"): GAIN_OF_FUNCTION}
    table = pd.read_csv("configs/allele_mechanism_adjudications.csv", dtype=str)
    assert table.pmid.str.fullmatch(r"\d+").all() and (table.quote.str.len() > 40).all() and table.note.str.contains("Gene2Phenotype").all()


@pytest.mark.parametrize("row,message", [({"mechanism": "neomorph"}, "mechanism must be one of"), ({"pmid": "x"}, "PubMed id"), ({"quote": " "}, "PubMed id and a quote"),
                                         ({"disease_id": "MONDO:1"}, "not an OMIM or Orphanet")])
def test_a_malformed_adjudication_is_refused(tmp_path, row, message):
    from mechanistic_pathway_learning.evidence.allele_mechanisms import read_mechanism_adjudications

    base = {"gene": "A", "disease_id": "OMIM:1", "mechanism": "loss_of_function", "pmid": "1", "quote": "a sentence from the paper", "note": ""}
    pd.DataFrame([{**base, **row}]).to_csv(tmp_path / "adjudications.csv", index=False)
    with pytest.raises(ValueError, match=message):
        read_mechanism_adjudications(tmp_path / "adjudications.csv")
    pd.DataFrame([base, base]).to_csv(tmp_path / "twice.csv", index=False)
    with pytest.raises(ValueError, match="adjudicated twice"):
        read_mechanism_adjudications(tmp_path / "twice.csv")


def test_an_adjudication_replaces_both_sources_and_keeps_their_view_in_the_basis():
    from mechanistic_pathway_learning.evidence.allele_mechanisms import apply_adjudications

    mechanisms = {("CDC73", "ORPHA:99880"): DiseaseMechanism(GAIN_OF_FUNCTION, "orphanet typed; disputed by gene2phenotype G2P01833"),
                  ("B", "OMIM:2"): DiseaseMechanism(UNTYPED, "untyped")}
    adjudicated = apply_adjudications(mechanisms, {("CDC73", "ORPHA:99880"): DiseaseMechanism(LOSS_OF_FUNCTION, "adjudicated from PMID 12434154"),
                                                   ("C", "OMIM:3"): DiseaseMechanism(GAIN_OF_FUNCTION, "adjudicated from PMID 1")})
    assert adjudicated[("CDC73", "ORPHA:99880")] == DiseaseMechanism(LOSS_OF_FUNCTION, "adjudicated from PMID 12434154 (the sources: orphanet typed; disputed by gene2phenotype G2P01833)")
    assert adjudicated[("B", "OMIM:2")] == mechanisms[("B", "OMIM:2")] and ("C", "OMIM:3") not in adjudicated, "a pair without a report is not added"
    reports = [monogenic_report("CDC73", "ORPHA:99880", "fatigue", 0)]
    assigned = assign_alleles(reports, adjudicated, {}, {"CDC73": "GENE:CDC73"})
    assert reports[0].perturbation_id == "CDC73" and json.loads(reports[0].perturbation_nodes) == [["GENE:CDC73", -1.0, 1.0]]
    summary = summarise_alleles(assigned, reports)
    assert summary["pairs_adjudicated"] == ["CDC73 ORPHA:99880: loss_of_function, adjudicated from PMID 12434154 (the sources: orphanet typed; disputed by gene2phenotype G2P01833)"]
    assert summary["pairs_disputed_by_the_second_source"] == [] and summary["reports_by_mechanism_basis"] == {"adjudicated": 1}
