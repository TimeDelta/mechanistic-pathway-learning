"""Tests for the three label changes of 10 October 2026 (docs/label_source_review.md): the gain-of-function seed sign,
the written list of admitted ATC groups and the mask on drug pairs never more frequent on the drug than on placebo.
Each is off by default, so the tables built before it are reproduced without a flag."""
import gzip
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.assemble_evidence_table import seed_monogenic_reports
from mechanistic_pathway_learning.evidence.gene_disease_association_types import (
    GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE,
    genes_with_gain_of_function_associations_only,
)
from mechanistic_pathway_learning.evidence.load_drug_label_events import (
    NERVOUS_SYSTEM_ATC_PREFIXES,
    has_admitted_atc_code,
    is_nervous_system_atc,
    load_sider_frequency_rows,
    read_admitted_atc_prefixes,
)
from mechanistic_pathway_learning.evidence.load_onsides_label_events import aggregate_label_statements, onsides_reports
from mechanistic_pathway_learning.evidence.placebo_arm_comparison import PlaceboArmComparison, compare_drug_arm_with_placebo_arm, pairs_never_above_placebo
from tests.test_evidence_reports_and_aggregation import monogenic_fixture
from tests.test_onsides_label_events import CROSSWALK, MECHANISMS, NODES, TARGETS, sertraline_bridge, toy_tables

LOSS_OF_FUNCTION_TYPE = "Disease-causing germline mutation(s) (loss of function) in"
HPOA_ROWS = ["OMIM:311250\tOTC deficiency\t\tHP:0000726\tPMID:2\tPCS\t\t5/10\t\t\tP\tHPO:a[2014-03-04]"]


def otc_reports(tmp_path: Path, association_type: str, causal: float):
    reports, *_ = monogenic_fixture(tmp_path, HPOA_ROWS, {"OTC"})
    assert reports
    for report in reports:
        report.perturbation_nodes, report.association_type, report.rubric_causal_association = "", association_type, causal
    return reports


def test_a_gene_is_gain_of_function_only_when_every_causal_association_says_so() -> None:
    gain = GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE
    assert genes_with_gain_of_function_associations_only({"A": {gain}, "B": {gain, "MENDELIAN"}, "C": {LOSS_OF_FUNCTION_TYPE}, "D": set()}) == {"A"}


def test_every_gene_is_seeded_as_a_loss_of_function_without_the_flag(tmp_path: Path) -> None:
    reports = otc_reports(tmp_path, GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE, 1.0)
    assert seed_monogenic_reports(reports, {"OTC": "GENE:OTC"}) == set()
    assert all(report.perturbation_nodes == json.dumps([["GENE:OTC", -1.0, 1.0]]) for report in reports)
    assert all(report.model_description.startswith("human loss-of-function") for report in reports)


def test_the_flag_seeds_a_gain_of_function_gene_with_sign_plus_one(tmp_path: Path) -> None:
    reports = otc_reports(tmp_path, GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE, 1.0)
    assert seed_monogenic_reports(reports, {"OTC": "GENE:OTC"}, seed_gain_of_function_genes_positive=True) == {"OTC"}
    assert all(report.perturbation_nodes == json.dumps([["GENE:OTC", 1.0, 1.0]]) for report in reports)
    assert all(report.model_description.startswith("human gain-of-function") for report in reports)


def test_the_flag_leaves_a_gene_with_another_causal_association_or_a_non_causal_gain_as_it_was(tmp_path: Path) -> None:
    mixed = otc_reports(tmp_path, GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE, 1.0) + otc_reports(tmp_path, "MENDELIAN", 1.0)
    assert seed_monogenic_reports(mixed, {"OTC": "GENE:OTC"}, seed_gain_of_function_genes_positive=True) == set()
    assert all(report.perturbation_nodes == json.dumps([["GENE:OTC", -1.0, 1.0]]) for report in mixed)
    # a candidate-gene report states no mechanism: with it the gene is still gain of function only, without a causal report it is not
    with_candidate = otc_reports(tmp_path, GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE, 1.0) + otc_reports(tmp_path, "Candidate gene tested in", 0.0)
    assert seed_monogenic_reports(with_candidate, {"OTC": "GENE:OTC"}, seed_gain_of_function_genes_positive=True) == {"OTC"}
    candidate_only = otc_reports(tmp_path, "Candidate gene tested in", 0.0)
    assert seed_monogenic_reports(candidate_only, {"OTC": "GENE:OTC"}, seed_gain_of_function_genes_positive=True) == set()


def write_prefix_list(path: Path, prefixes: list[str]) -> Path:
    path.write_text("atc_prefix,name,kind,reason,source\n" + "".join(f"{prefix},name,group,reason,source\n" for prefix in prefixes))
    return path


def test_the_default_admitted_list_is_the_nervous_system_rule() -> None:
    for codes in (["N05BA01"], ["C02AC02"], ["M03AB", "N01AX"], []):
        assert has_admitted_atc_code(codes) == is_nervous_system_atc(codes)
    assert NERVOUS_SYSTEM_ATC_PREFIXES == ("N",)


def test_a_written_list_admits_groups_and_single_substances(tmp_path: Path) -> None:
    prefixes = read_admitted_atc_prefixes(write_prefix_list(tmp_path / "list.csv", ["N", "H", "C02A", "A03FA01"]))
    assert prefixes == ("N", "H", "C02A", "A03FA01")
    assert has_admitted_atc_code(["C02AC02"], prefixes) and has_admitted_atc_code(["H02AB06"], prefixes) and has_admitted_atc_code(["A03FA01"], prefixes)
    assert not has_admitted_atc_code(["C02CA01"], prefixes)  # a peripherally acting antiadrenergic agent
    assert not has_admitted_atc_code(["A03FA03"], prefixes)  # another substance of the named drug's subgroup
    assert not has_admitted_atc_code(["A03FA"], prefixes)  # the subgroup code alone does not name the substance
    assert not has_admitted_atc_code(["N05BA01"], ("H",))  # the list replaces group N, it is not added to it


def test_a_written_list_is_refused_when_empty_repeated_or_nested(tmp_path: Path) -> None:
    for prefixes in ([], ["N", "N"], ["N", "N05"]):
        with pytest.raises(ValueError):
            read_admitted_atc_prefixes(write_prefix_list(tmp_path / "list.csv", prefixes))


def test_onsides_ingredients_qualify_by_the_same_written_list() -> None:
    statements = aggregate_label_statements(toy_tables(), CROSSWALK)
    cardiac = {"36437": sertraline_bridge(atc_codes_rxnav=["C02AC02"])}
    rejected, rejected_counts = onsides_reports(statements, cardiac, MECHANISMS, TARGETS, NODES)
    assert not rejected.qualifies.any() and set(rejected_counts["disqualified"]) == {"atc_not_nervous_system"}
    admitted, _ = onsides_reports(statements, cardiac, MECHANISMS, TARGETS, NODES, admitted_atc_prefixes=("N", "C02A"))
    assert admitted[admitted.single_ingredient_label_count > 0].qualifies.all()
    other_list, other_counts = onsides_reports(statements, cardiac, MECHANISMS, TARGETS, NODES, admitted_atc_prefixes=("N", "H"))
    assert not other_list.qualifies.any() and set(other_counts["disqualified"]) == {"atc_not_admitted"}


FREQUENCY_ROWS = [
    # drug A, insomnia: one label, drug arm below its placebo arm
    "CID100000001\tCID000000001\tLABEL1\t\t4%\t0.04\t0.04\tPT\tC0917801\tInsomnia",
    "CID100000001\tCID000000001\tLABEL1\tplacebo\t6%\t0.06\t0.06\tPT\tC0917801\tInsomnia",
    # drug A, anxiety: below placebo in one label, above it in another
    "CID100000001\tCID000000001\tLABEL1\t\t2%\t0.02\t0.02\tPT\tC0003467\tAnxiety",
    "CID100000001\tCID000000001\tLABEL1\tplacebo\t3%\t0.03\t0.03\tPT\tC0003467\tAnxiety",
    "CID100000001\tCID000000001\tLABEL2\t\t5%\t0.05\t0.05\tPT\tC0003467\tAnxiety",
    "CID100000001\tCID000000001\tLABEL2\tplacebo\t3%\t0.03\t0.03\tPT\tC0003467\tAnxiety",
    # drug B, insomnia: a drug arm in one label and a placebo arm in another, so no entry gives both
    "CID100000002\tCID000000002\tLABEL3\t\t4%\t0.04\t0.04\tPT\tC0917801\tInsomnia",
    "CID100000002\tCID000000002\tLABEL4\tplacebo\t9%\t0.09\t0.09\tPT\tC0917801\tInsomnia",
    # drug B, anxiety: equal arms count as not above
    "CID100000002\tCID000000002\tLABEL3\t\t3%\t0.03\t0.03\tPT\tC0003467\tAnxiety",
    "CID100000002\tCID000000002\tLABEL3\tplacebo\t3%\t0.03\t0.03\tPT\tC0003467\tAnxiety",
    # a term outside the crosswalk
    "CID100000002\tCID000000002\tLABEL3\t\t1%\t0.01\t0.01\tPT\tC0015230\tRash",
    "CID100000002\tCID000000002\tLABEL3\tplacebo\t5%\t0.05\t0.05\tPT\tC0015230\tRash",
]


def write_frequency_table(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / "meddra_freq.tsv.gz", "wt", encoding="utf-8") as frequency_file:
        frequency_file.write("\n".join(FREQUENCY_ROWS) + "\n")
    return directory


def test_arms_are_compared_within_one_label_and_one_term(tmp_path: Path) -> None:
    comparisons = compare_drug_arm_with_placebo_arm(load_sider_frequency_rows(write_frequency_table(tmp_path / "sider")), {"insomnia": "insomnia", "anxiety": "anxiety"})
    assert comparisons == {
        ("CID100000001", "insomnia"): PlaceboArmComparison(1, 0),
        ("CID100000001", "anxiety"): PlaceboArmComparison(2, 1),
        ("CID100000002", "anxiety"): PlaceboArmComparison(1, 0),
    }
    assert pairs_never_above_placebo(comparisons) == {("CID100000001", "insomnia"), ("CID100000002", "anxiety")}


def test_two_terms_of_one_symptom_are_not_pooled(tmp_path: Path) -> None:
    # one term above its placebo arm keeps the pair, however far the other term's placebo arm is above its drug arm
    directory = tmp_path / "sider"
    directory.mkdir()
    with gzip.open(directory / "meddra_freq.tsv.gz", "wt", encoding="utf-8") as frequency_file:
        frequency_file.write("\n".join([
            "CID100000001\tCID000000001\tLABEL1\t\t10%\t0.10\t0.10\tPT\tC1\tInsomnia",
            "CID100000001\tCID000000001\tLABEL1\tplacebo\t8%\t0.08\t0.08\tPT\tC1\tInsomnia",
            "CID100000001\tCID000000001\tLABEL1\t\t1%\t0.01\t0.01\tPT\tC2\tMiddle insomnia",
            "CID100000001\tCID000000001\tLABEL1\tplacebo\t9%\t0.09\t0.09\tPT\tC2\tMiddle insomnia",
        ]) + "\n")
    comparisons = compare_drug_arm_with_placebo_arm(load_sider_frequency_rows(directory), {"insomnia": "insomnia", "middle insomnia": "insomnia"})
    assert comparisons == {("CID100000001", "insomnia"): PlaceboArmComparison(2, 1)}


def evidence_rows() -> pd.DataFrame:
    def row(perturbation_id: str, perturbation_type: str, symptom: str, frequency: float) -> dict:
        return {"perturbation_id": perturbation_id, "perturbation_type": perturbation_type, "perturbation_label": perturbation_id, "symptom": symptom, "relation": "induces",
                "grade": "A" if perturbation_type == "gene" else "B", "label_frequency": frequency, "source": "HPO-OMIM" if perturbation_type == "gene" else "SIDER-label",
                "positive_report_count": 1}
    return pd.DataFrame([
        row("CID100000001", "drug", "insomnia", 0.04), row("CID100000001", "drug", "anxiety", 0.035), row("CID100000002", "drug", "insomnia", 0.04),
        row("CID100000002", "drug", "anxiety", 0.03), row("CID100000002", "drug", "fatigue", 0.001), row("GENEA", "gene", "insomnia", 0.5),
    ])


def test_the_selection_sets_aside_kept_drug_pairs_never_above_placebo_only_with_the_flag(tmp_path: Path) -> None:
    evidence_directory = tmp_path / "evidence"
    evidence_directory.mkdir()
    evidence_rows().to_parquet(evidence_directory / "evidence_records.parquet")
    sider_directory = write_frequency_table(tmp_path / "sider")
    crosswalk = tmp_path / "crosswalk.csv"
    crosswalk.write_text("target_symptom,meddra_preferred_terms\ninsomnia,Insomnia\nanxiety,Anxiety\nfatigue,Fatigue\n")

    def build(name: str, *flags: str) -> tuple[pd.DataFrame, dict]:
        output = tmp_path / f"{name}.parquet"
        subprocess.run([sys.executable, "experiments/build_label_selection.py", "--evidence-dir", str(evidence_directory), "--output", str(output), *flags], check=True, capture_output=True)
        return pd.read_parquet(output).set_index(["perturbation_id", "symptom"]), json.loads(output.with_suffix(".summary.json").read_text())

    plain, plain_summary = build("plain")
    assert plain.keep.sum() == 5 and "mask_never_above_placebo" not in plain_summary
    masked, masked_summary = build("masked", "--mask-never-above-placebo", "--sider-dir", str(sider_directory), "--crosswalk", str(crosswalk))
    set_aside = masked[plain.keep & ~masked.keep]
    assert sorted(set_aside.index) == [("CID100000001", "insomnia"), ("CID100000002", "anxiety")]
    assert set_aside.reason.str.startswith("drug pair never more frequent on the drug than on placebo").all()
    assert masked.loc[("CID100000001", "anxiety"), "keep"]  # above placebo in one of its two labels
    assert masked.loc[("CID100000002", "insomnia"), "keep"]  # no entry with both arms: not judged
    assert masked.loc[("GENEA", "insomnia"), "keep"]  # a gene pair is never judged
    assert masked.loc[("CID100000002", "fatigue"), "reason"] == plain.loc[("CID100000002", "fatigue"), "reason"]  # already set aside for its frequency
    assert masked_summary["mask_never_above_placebo"] == {"drug_pairs_with_a_label_entry_giving_both_arms": 3, "of_them_never_above_placebo": 2,
                                                          "kept_pairs_set_aside": 2, "kept_pairs_set_aside_by_symptom": {"anxiety": 1, "insomnia": 1}}
