"""experiments/measure_genebass_label_agreement.py: the monogenic labels against Genebass's pLoF burden results."""
import gzip
import json

import pandas as pd
import pytest

from experiments.measure_genebass_label_agreement import (EXCLUDED_ITEMS, ITEM_TO_CONSTRUCT, LOWER_IS_THE_SYMPTOM, agreement,
                                                          analysis_id_of_address, labelled_genes, mapped_analyses,
                                                          nominal_hits_in_symptom_direction, pooled, read_saved_json,
                                                          requirement_of_gene, saved_pages)


def gene_table(rows):
    return pd.DataFrame(rows, columns=["gene_symbol", "BETA_Burden", "Pvalue_Burden"])


def test_a_hit_needs_a_nominal_p_and_the_symptom_direction():
    table = gene_table([("A", 0.2, 0.01), ("B", -0.2, 0.01), ("C", 0.2, 0.2), ("D", None, 0.01)])
    assert nominal_hits_in_symptom_direction(table, lower_is_the_symptom=False) == {"A": True, "B": False, "C": False}


def test_for_fluid_intelligence_a_lower_score_is_the_symptom():
    table = gene_table([("A", 0.2, 0.01), ("B", -0.2, 0.01)])
    assert nominal_hits_in_symptom_direction(table, lower_is_the_symptom=True) == {"A": False, "B": True}
    assert "20016" in LOWER_IS_THE_SYMPTOM and ITEM_TO_CONSTRUCT["20016"][0] == "cognitive_impairment"


def test_the_threshold_is_strict():
    table = gene_table([("A", 0.2, 0.05)])
    assert nominal_hits_in_symptom_direction(table, lower_is_the_symptom=False) == {"A": False}


def test_conditional_and_two_directional_items_are_not_mapped():
    for code in EXCLUDED_ITEMS:
        assert code not in ITEM_TO_CONSTRUCT, f"{code} is excluded and must not be mapped"
    assert {"20534", "20535"} <= set(EXCLUDED_ITEMS), "the worst-episode sleep items have a depressed denominator"


def test_only_both_sexes_analyses_of_mapped_items_are_kept():
    phenotypes = [
        {"phenocode": "20510", "pheno_sex": "both_sexes", "analysis_id": "continuous-20510-both_sexes--", "trait_type": "continuous", "n_cases": 10, "n_controls": None},
        {"phenocode": "20510", "pheno_sex": "females", "analysis_id": "continuous-20510-females--", "trait_type": "continuous", "n_cases": 6, "n_controls": None},
        {"phenocode": "20534", "pheno_sex": "both_sexes", "analysis_id": "categorical-20534-both_sexes--", "trait_type": "categorical", "n_cases": 5, "n_controls": 7},
        {"phenocode": "120", "pheno_sex": "both_sexes", "analysis_id": "categorical-120-both_sexes-1-", "trait_type": "categorical", "n_cases": 5, "n_controls": 7},
    ]
    analyses = mapped_analyses(phenotypes)
    assert analyses["analysis_id"].tolist() == ["continuous-20510-both_sexes--"]
    assert analyses["construct"].tolist() == ["depressed_mood"]


def test_an_analysis_id_is_read_from_a_gene_manhattan_address_only():
    assert analysis_id_of_address("https://main.genebass.org/api/analysis/continuous-20510-both_sexes--/gene-manhattan?burdenSet=pLoF") == "continuous-20510-both_sexes--"
    assert analysis_id_of_address("https://main.genebass.org/api/phenotypes") is None


def test_saved_pages_are_found_by_their_metadata_and_gzip_is_undone(tmp_path):
    payload = [{"gene_symbol": "A", "BETA_Burden": 0.1, "Pvalue_Burden": 0.01}]
    (tmp_path / "abc.txt").write_bytes(gzip.compress(json.dumps(payload).encode()))
    (tmp_path / "abc.txt.json").write_text(json.dumps({"address": "https://main.genebass.org/api/x", "read": True}))
    (tmp_path / "refused.txt.json").write_text(json.dumps({"address": "https://example.org/y", "read": False}))
    (tmp_path / "plain.txt").write_text(json.dumps({"ok": 1}))
    (tmp_path / "plain.txt.json").write_text(json.dumps({"address": "https://main.genebass.org/api/z", "read": True}))
    pages = saved_pages(tmp_path)
    assert set(pages) == {"https://main.genebass.org/api/x", "https://main.genebass.org/api/z"}, "a refused page is not a saved page"
    assert read_saved_json(pages["https://main.genebass.org/api/x"]) == payload
    assert read_saved_json(pages["https://main.genebass.org/api/z"]) == {"ok": 1}


def test_allelic_requirement_counts_autosomal_records_that_are_not_disputed():
    g2p = pd.DataFrame({
        "gene symbol": ["MONO", "BI", "BOTH", "BOTH", "DISPUTED", "XLINKED"],
        "allelic requirement": ["monoallelic_autosomal", "biallelic_autosomal", "monoallelic_autosomal", "biallelic_autosomal",
                                "monoallelic_autosomal", "monoallelic_X_hemizygous"],
        "confidence": ["definitive", "strong", "moderate", "limited", "disputed", "definitive"],
    })
    assert requirement_of_gene(g2p) == {"MONO": "monoallelic", "BI": "biallelic", "BOTH": "both"}


def test_labels_are_gene_perturbations_that_induce_the_construct():
    evidence = pd.DataFrame({
        "perturbation_id": ["A", "B", "CID1", "C"],
        "perturbation_type": ["gene", "gene", "drug", "gene"],
        "relation": ["induces", "induces", "induces", "relieves"],
        "symptom": ["anxiety", "insomnia", "anxiety", "anxiety"],
    })
    assert labelled_genes(evidence) == {"anxiety": {"A"}, "insomnia": {"B"}}


def test_agreement_counts_hits_against_the_items_background():
    analyses = pd.DataFrame([{"analysis_id": "continuous-20506-both_sexes--", "phenocode": "20506", "construct": "anxiety", "cases": 100}])
    results = {"continuous-20506-both_sexes--": gene_table([("A", 0.3, 0.001), ("B", 0.3, 0.6), ("C", -0.3, 0.001), ("D", 0.1, 0.9)])}
    labels = {"anxiety": {"A", "B", "UNTESTED"}}
    rows = agreement(analyses, results, labels, lambda gene: "all")
    row = rows.iloc[0]
    assert (row.labelled, row.tested, row.hits) == (3, 2, 1)
    assert row.background_rate == pytest.approx(0.25), "one hit in the symptom direction among the four genes tested"
    assert row.hit_genes == "A"


def test_agreement_skips_an_item_with_no_saved_results_and_splits_by_group():
    analyses = pd.DataFrame([
        {"analysis_id": "continuous-20506-both_sexes--", "phenocode": "20506", "construct": "anxiety", "cases": 100},
        {"analysis_id": "continuous-20510-both_sexes--", "phenocode": "20510", "construct": "depressed_mood", "cases": 100},
    ])
    results = {"continuous-20506-both_sexes--": gene_table([("A", 0.3, 0.001), ("B", 0.3, 0.6)])}
    rows = agreement(analyses, results, {"anxiety": {"A", "B"}, "depressed_mood": {"A"}}, lambda gene: "group " + gene)
    assert rows["analysis_id"].unique().tolist() == ["continuous-20506-both_sexes--"]
    assert dict(zip(rows["group"], rows["hits"])) == {"group A": 1, "group B": 0}


def test_pooling_weights_each_item_by_its_own_background():
    rows = pd.DataFrame([
        {"group": "g", "tested": 10, "hits": 2, "background_rate": 0.1},
        {"group": "g", "tested": 20, "hits": 1, "background_rate": 0.05},
    ])
    result = pooled(rows).iloc[0]
    assert (result.pairs_tested, result.hits) == (30, 3)
    assert result.expected == pytest.approx(2.0)
    assert result.ratio == pytest.approx(1.5)
    assert 0 < result.binomial_p < 1
