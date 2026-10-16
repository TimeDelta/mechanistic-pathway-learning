"""Per-edge evidence and the edge priors computed from it (graph/edge_evidence.py, decision 9.57)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.edge_evidence import (
    INTACT_ABSENCE_PSEUDO_COUNT,
    EdgeTier,
    accession_of,
    edge_evidence_table,
    edge_priors,
    family_of,
    members,
    permute_edge_priors,
    pmids_of,
    read_edge_tiers,
    signalling_references,
    summarise,
    symbol_of,
    tempered,
    transcription_references,
)


def test_families_symbols_and_accessions_are_read_from_the_edge_itself():
    assert family_of("OmniPath") == "signalling" and family_of("OmniPath small_molecule_protein") == "small_molecule"
    assert family_of("CollecTRI") == "transcription" and family_of("BRENDA CC BY 4.0 (3 references, accession stated by BRENDA)") == "enzyme_regulation"
    assert family_of("curated oxidative modification (PMID 14585973)") == "curated" and family_of("Human-GEM") == "other"
    assert symbol_of("PROTEIN:CALCA:P01258") == "CALCA" and symbol_of("GENE:CHEBI:15996") == "CHEBI:15996" and symbol_of("GENE:JUN") == "JUN"
    assert symbol_of("MAM01570e") == "MAM01570e"
    assert accession_of("PROTEIN:CALCA:P01258") == "P01258" and accession_of("PROTEIN:A0A087WWS6") == "A0A087WWS6" and accession_of("PROTEIN:JUN") is None
    assert members("CALM1_CALM2") == ["CALM1", "CALM2"] and members("COMPLEX:P05412_P01100-2") == ["P05412", "P01100"]


def test_pubmed_ids_are_split_by_resource_and_deduplicated():
    kept, excluded = pmids_of("SIGNOR:1;HPRD:1;ProtMapper:2;ProtMapper:1;KEA:x;bad", ("ProtMapper",))
    assert kept == {"1"} and excluded == {"2", "1"}
    assert pmids_of(float("nan")) == (set(), set())


def omnipath_rows():
    return pd.DataFrame({
        "source": ["P1", "P3", "COMPLEX:P5_P6", "P7"], "target": ["P2", "P4", "P8", "P9"],
        "source_genesymbol": ["A", "C", "E_F", "G"], "target_genesymbol": ["B", "D", "H", "I"],
        "sources": ["SIGNOR", "ProtMapper", "CellChatDB", "SIGNOR;ProtMapper"],
        "references": ["SIGNOR:11;SIGNOR:12;HPRD:11", "ProtMapper:21", None, "ProtMapper:31"],
    })


def collectri_rows():
    return pd.DataFrame({
        "source": ["Q1", "Q2"], "target": ["Q3", "Q4"], "source_genesymbol": ["TF1", "TF2"], "target_genesymbol": ["T1", "T2"],
        "sources": ["CollecTRI;TRRUST_CollecTRI;ExTRI_CollecTRI", "CollecTRI;ExTRI_CollecTRI;ExTRI2_CollecTRI2"],
        "references": ["CollecTRI:41;TRRUST:41;CollecTRI:42", "CollecTRI:51"],
    })


def test_reference_tables_count_distinct_pmids_expand_complexes_and_flag_text_mining():
    signalling = signalling_references(omnipath_rows()).set_index(["source_symbol", "target_symbol"])
    assert signalling.loc[("A", "B")].curated_pmids == 2 and not signalling.loc[("A", "B")].text_mined_only
    assert signalling.loc[("C", "D")].curated_pmids == 0 and signalling.loc[("C", "D")].text_mined_only and signalling.loc[("C", "D")].text_mined_pmids == 1
    assert ("E", "H") in signalling.index and ("F", "H") in signalling.index and ("P5", "P8") in signalling.index, "complex members, by symbol and by accession"
    assert signalling.loc[("E", "H")].curated_pmids == 0 and not signalling.loc[("E", "H")].text_mined_only
    transcription = transcription_references(collectri_rows()).set_index(["source_symbol", "target_symbol"])
    assert transcription.loc[("TF1", "T1")].curated_pmids == 2 and not transcription.loc[("TF1", "T1")].text_mined_only
    assert transcription.loc[("TF2", "T2")].text_mined_only, "every supporting resource but the aggregators is an ExTRI layer"


def graph_edges():
    return pd.DataFrame({
        "source_id": ["PROTEIN:A", "PROTEIN:B", "PROTEIN:C", "PROTEIN:E", "PROTEIN:Z", "PROTEIN:TF1", "PROTEIN:TF2", "CHEBI1c", "MAM02041c", "MAM1c", "PROTEIN:B"],
        "target_id": ["PROTEIN:B", "PROTEIN:A", "PROTEIN:D", "PROTEIN:H", "PROTEIN:Y", "GENE:T1", "GENE:T2", "PROTEIN:AMY2A", "PROTEIN:KEAP1", "MAR1", "PROTEIN:A"],
        "relation_type": ["activates", "binds", "activates", "binds", "activates", "regulates_transcription_of", "regulates_transcription_of", "inhibits",
                          "oxidatively_modifies", "substrate_of", "activates"],
        "evidence_source": ["OmniPath", "OmniPath", "OmniPath", "OmniPath", "OmniPath", "CollecTRI", "CollecTRI", "BRENDA CC BY 4.0 (5 references, accession stated by BRENDA)",
                            "curated oxidative modification (PMID 14585973)", "Human-GEM", "OmniPath"],
    })


def test_every_graph_edge_gets_a_family_kind_record_count_and_contest_flag():
    references = pd.concat([signalling_references(omnipath_rows()), transcription_references(collectri_rows())], ignore_index=True)
    evidence = edge_evidence_table(graph_edges(), references, {frozenset(("A", "B"))})
    rows = list(evidence.itertuples(index=False))
    assert (rows[0].family, rows[0].kind, rows[0].records, rows[0].contested) == ("signalling", "curated", 2, True)
    assert (rows[1].kind, rows[1].records) == ("curated", 2), "an undirected binding takes the evidence of the reverse raw row"
    assert rows[10].kind == "unmatched", "a directed edge never takes the evidence of its reverse, which is another claim"
    assert (rows[2].kind, rows[2].records) == ("text_mined_only", 1) and rows[3].kind == "unreferenced" and rows[4].kind == "unmatched"
    assert (rows[5].kind, rows[5].records) == ("curated", 2) and (rows[6].kind, rows[6].records) == ("text_mined_only", 1)
    assert (rows[7].family, rows[7].kind, rows[7].records) == ("enzyme_regulation", "in_vitro", 5)
    assert (rows[8].family, rows[8].kind, rows[8].records) == ("curated", "curated", 1) and (rows[9].family, rows[9].kind) == ("other", "none")
    assert list(evidence.index) == list(graph_edges().index)


TIERS = {("signalling", "curated"): EdgeTier("signalling", "curated", 0.9, 1.0), ("signalling", "text_mined_only"): EdgeTier("signalling", "text_mined_only", 0.6, 0.5),
         ("signalling", "unreferenced"): EdgeTier("signalling", "unreferenced", 0.7, 0.5), ("enzyme_regulation", "in_vitro"): EdgeTier("enzyme_regulation", "in_vitro", 0.6, 0.5)}


def test_priors_are_tiered_tempered_capped_and_updated_by_a_curated_absence():
    evidence = pd.DataFrame({"family": ["signalling", "signalling", "signalling", "enzyme_regulation", "other", "signalling"],
                             "kind": ["curated", "curated", "unreferenced", "in_vitro", "none", "unmatched"], "records": [1, 9, 0, 5, 0, 0],
                             "contested": [False, True, False, False, False, False]})
    gated, probability, strength = edge_priors(evidence, TIERS, global_strength=2.0)
    assert gated.tolist() == [True, True, True, True, False, False]
    assert probability[0] == pytest.approx(0.9) and strength[0] == pytest.approx(2.0 * 1.0 * 1.0)
    capped = 2.0 * 1.0 * float(tempered(2))
    alpha, beta = 0.9 * capped, 0.1 * capped + 2.0 * INTACT_ABSENCE_PSEUDO_COUNT
    assert strength[1] == pytest.approx(alpha + beta) and probability[1] == pytest.approx(alpha / (alpha + beta)) and probability[1] < 0.9
    assert probability[2] == pytest.approx(0.7) and strength[2] == pytest.approx(2.0 * 0.5), "an unreferenced tier still gives its unit strength once"
    assert strength[3] == pytest.approx(2.0 * 0.5 * float(tempered(5))) and np.isnan(probability[4]) and strength[4] == 0.0 and not gated[5]
    summary = summarise(evidence, gated, probability, strength)
    assert summary["gated"] == 4 and summary["contested"] == 1 and summary["by_family_and_kind"]["other:none"]["mean_prior"] is None


def test_the_permuted_control_shuffles_priors_within_a_relation_and_keeps_the_rest():
    gated = np.array([True, True, True, False, True, True])
    probability = np.array([0.9, 0.6, 0.7, np.nan, 0.85, 0.55])
    strength = np.array([1.0, 0.5, 0.5, 0.0, 1.2, 0.6])
    relation = np.array([0, 0, 0, 0, 1, 1])
    moved_probability, moved_strength = permute_edge_priors(gated, probability, strength, relation, random_seed=3)
    assert sorted(moved_probability[:3]) == sorted(probability[:3]) and sorted(moved_probability[4:]) == sorted(probability[4:])
    assert np.isnan(moved_probability[3]) and moved_strength[3] == 0.0
    pairs = {(p, s) for p, s in zip(probability[gated], strength[gated])}
    assert {(p, s) for p, s in zip(moved_probability[gated], moved_strength[gated])} == pairs, "a prior's probability and strength travel together"
    again = permute_edge_priors(gated, probability, strength, relation, random_seed=3)
    assert np.array_equal(again[0], moved_probability, equal_nan=True)


def test_the_committed_tier_table_is_valid_and_ordered_by_evidence_type():
    tiers = read_edge_tiers()
    assert tiers[("signalling", "curated")].prior_probability > tiers[("signalling", "unreferenced")].prior_probability > tiers[("signalling", "text_mined_only")].prior_probability
    assert tiers[("transcription", "curated")].prior_probability > tiers[("transcription", "text_mined_only")].prior_probability
    assert all(0.0 < tier.prior_probability < 1.0 for tier in tiers.values())
    assert ("other", "none") not in tiers and ("signalling", "unmatched") not in tiers, "ungated families have no row"
