"""Tests for the literature class E3 retrieval on small synthetic fixtures: relation-type mapping, descriptor-level flags,
gene identifier precedence, CTD evidence mapping, chemical name matching and the report schema. No network."""
from __future__ import annotations

import io
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from mechanistic_pathway_learning.evidence.gene_identifier_map import map_gene_nodes, ncbi_gene_id_to_symbol
from mechanistic_pathway_learning.evidence.load_ctd_relations import (
    DIRECT_EVIDENCE_TO_RELATION,
    ctd_rows_to_reports,
    expand_ctd_row,
    filter_ctd_rows,
    iterate_ctd_rows,
)
from mechanistic_pathway_learning.evidence.load_pubtator_relations import (
    EVIDENCE_REPORT_COLUMNS,
    LITERATURE_REPORT_COLUMNS,
    RELATION_TYPE_TO_RELATION,
    disease_rows_from_stream,
    filter_relation_rows,
    iterate_bulk_relation_rows,
    map_relation_type,
    parse_entity,
    relation_rows_to_reports,
    summarize_reports,
)
from mechanistic_pathway_learning.evidence.mesh_symptom_descriptors import (
    SymptomMeshDescriptor,
    descriptor_lookup,
    load_symptom_mesh_descriptors,
    parse_descriptor_columns,
)
from mechanistic_pathway_learning.evidence.pubmed_publication_dates import parse_pubdate
from mechanistic_pathway_learning.evidence.sider_chemical_matching import (
    ChemicalMatch,
    chemical_matches_by_mesh_id,
    drug_names_by_normalized_name,
    match_autocomplete_candidates,
    match_ctd_chemical_name,
    match_mesh_label,
)

FIXTURE_DESCRIPTORS = descriptor_lookup(
    [
        SymptomMeshDescriptor("anxiety", "D001007", "symptom"),
        SymptomMeshDescriptor("anxiety", "D001008", "diagnosis"),
        SymptomMeshDescriptor("fatigue", "D005221", "symptom"),
    ]
)
FLUOXETINE = ChemicalMatch("CID100003386", "fluoxetine", "D005473", "Fluoxetine", "autocomplete_name")


# relation-type mapping


def test_relation_types_map_to_induces_relieves_and_associated_with() -> None:
    assert map_relation_type("cause") == ("induces", 1)
    assert map_relation_type("positive_correlate") == ("induces", 1)
    assert map_relation_type("positive_correlation") == ("induces", 1)
    assert map_relation_type("treat") == ("relieves", 1)
    assert map_relation_type("prevent") == ("relieves", 1)
    assert map_relation_type("negative_correlate") == ("relieves", 1)
    assert map_relation_type("Associate") == ("associated_with", 1)  # kept as its own undirected relation, never merged into induces
    assert map_relation_type("association") == ("associated_with", 1)
    for chemical_gene_type in ("stimulate", "inhibit"):
        assert map_relation_type(chemical_gene_type) is None
    assert map_relation_type("cotreat") is None
    assert set(RELATION_TYPE_TO_RELATION.values()) == {"induces", "relieves", "associated_with"}


def test_parse_entity_strips_mesh_prefix_and_splits_several_gene_ids() -> None:
    assert parse_entity("Disease|MESH:D001007") == ("Disease", ["D001007"])
    assert parse_entity("Chemical|MESH:C035444") == ("Chemical", ["C035444"])
    assert parse_entity("Gene|1234;5678") == ("Gene", ["1234", "5678"])
    assert parse_entity("Species|9606") == ("Species", ["9606"])
    assert parse_entity("malformed") == ("malformed", [])


def test_filter_keeps_graph_genes_and_sider_drugs_and_counts_every_drop() -> None:
    rows = [
        ("11", "cause", "Gene|8813", "Disease|MESH:D001007"),  # kept: graph gene, symptom-level
        ("12", "treat", "Chemical|MESH:D005473", "Disease|MESH:D001008"),  # kept: SIDER drug, diagnosis-level
        ("13", "associate", "Disease|MESH:D001007", "Gene|8813"),  # kept: disease first, associated_with
        ("14", "stimulate", "Chemical|MESH:D005473", "Disease|MESH:D001007"),  # dropped: chemical-gene type
        ("15", "inhibit", "Gene|8813", "Disease|MESH:D001007"),  # dropped: chemical-gene type
        ("16", "cotreat", "Chemical|MESH:D005473", "Disease|MESH:D001007"),  # dropped: not mapped
        ("17", "cause", "Gene|99", "Disease|MESH:D001007"),  # dropped: gene not in graph
        ("18", "treat", "Chemical|MESH:D000001", "Disease|MESH:D001007"),  # dropped: chemical not SIDER
        ("19", "associate", "Disease|MESH:D001007", "Disease|MESH:D005221"),  # dropped: both sides targets
        ("20", "associate", "Disease|MESH:D001007", "Species|9606"),  # dropped: partner type
        ("21", "cause", "Gene|8813;77", "Disease|MESH:D005221"),  # kept once: one of two ids in graph
    ]
    kept, counts = filter_relation_rows(rows, FIXTURE_DESCRIPTORS, {8813: "DPM1"}, chemical_matches_by_mesh_id([FLUOXETINE]))
    assert [row.pmid for row in kept] == ["11", "12", "13", "21"]
    assert [row.relation for row in kept] == ["induces", "relieves", "associated_with", "induces"]
    assert [row.entity_role for row in kept] == ["disease_second", "disease_second", "disease_first", "disease_second"]
    assert kept[1].perturbation_id == "CID100003386" and kept[1].chemical_match_method == "autocomplete_name"
    assert kept[0].perturbation_ncbi_gene_id == 8813 and kept[0].perturbation_id == "DPM1"
    assert counts["dropped_chemical_gene_relation_type:stimulate"] == 1
    assert counts["dropped_chemical_gene_relation_type:inhibit"] == 1
    assert counts["dropped_relation_type_not_mapped:cotreat"] == 1
    assert counts["dropped_gene_not_in_graph"] == 2  # row 17 and the second id of row 21
    assert counts["dropped_chemical_not_matched_to_sider"] == 1
    assert counts["dropped_both_sides_target_descriptors"] == 1
    assert counts["dropped_partner_type:Species"] == 1
    assert counts["rows_with_several_partner_identifiers"] == 1
    assert counts["rows_kept"] == 4


def test_bulk_stream_reader_and_first_pass_keep_only_target_descriptor_rows() -> None:
    text = "1\tcause\tGene|8813\tDisease|MESH:D001007\n2\ttreat\tChemical|MESH:D005473\tDisease|MESH:D009999\nshort\n"
    rows = list(iterate_bulk_relation_rows(io.StringIO(text)))
    assert len(rows) == 2
    kept, counts = disease_rows_from_stream(rows, FIXTURE_DESCRIPTORS)
    assert len(kept) == 1 and counts["bulk_rows_total"] == 2 and counts["bulk_rows_with_target_descriptor"] == 1
    assert counts["bulk_rows_by_type:treat"] == 1


# descriptor-level flags


def test_descriptor_columns_parse_aligned_lists_and_reject_malformed_ones() -> None:
    parsed = parse_descriptor_columns("anxiety", "D001007;D001008", "symptom;diagnosis")
    assert [(entry.mesh_descriptor, entry.level) for entry in parsed] == [("D001007", "symptom"), ("D001008", "diagnosis")]
    assert parse_descriptor_columns("psychomotor_retardation", "", "") == []
    with pytest.raises(ValueError):
        parse_descriptor_columns("anxiety", "D001007;D001008", "symptom")
    with pytest.raises(ValueError):
        parse_descriptor_columns("anxiety", "HP:0000739", "symptom")
    with pytest.raises(ValueError):
        parse_descriptor_columns("anxiety", "D001007", "disorder")


def test_crosswalk_loader_reads_levels_and_rejects_a_descriptor_under_two_symptoms(tmp_path: Path) -> None:
    crosswalk = tmp_path / "crosswalk.csv"
    crosswalk.write_text("target_symptom,notes,mesh_descriptors,mesh_descriptor_level\nanxiety,,D001007;D001008,symptom;diagnosis\nfatigue,,D005221,symptom\n")
    entries = load_symptom_mesh_descriptors(crosswalk)
    assert [(entry.target_symptom, entry.mesh_descriptor, entry.level) for entry in entries] == [("anxiety", "D001007", "symptom"), ("anxiety", "D001008", "diagnosis"), ("fatigue", "D005221", "symptom")]
    crosswalk.write_text("target_symptom,notes,mesh_descriptors,mesh_descriptor_level\nanxiety,,D001007,symptom\nfatigue,,D001007,symptom\n")
    with pytest.raises(ValueError):
        load_symptom_mesh_descriptors(crosswalk)


def test_descriptor_level_travels_into_the_report_and_its_limitations() -> None:
    rows = [("1", "treat", "Chemical|MESH:D005473", "Disease|MESH:D001008"), ("1", "treat", "Chemical|MESH:D005473", "Disease|MESH:D001007")]
    kept, _ = filter_relation_rows(rows, FIXTURE_DESCRIPTORS, {}, chemical_matches_by_mesh_id([FLUOXETINE]))
    reports = relation_rows_to_reports(kept)
    by_descriptor = reports.set_index("mesh_descriptor")
    assert by_descriptor.loc["D001008", "mesh_descriptor_level"] == "diagnosis"
    assert "diagnosis-level descriptor" in by_descriptor.loc["D001008", "limitations"]
    assert by_descriptor.loc["D001007", "mesh_descriptor_level"] == "symptom"
    assert "diagnosis-level descriptor" not in by_descriptor.loc["D001007", "limitations"]
    assert all("soft prior only" in text for text in reports.limitations)


def test_committed_crosswalk_has_the_two_mesh_columns_with_symptom_level_descriptors() -> None:
    entries = load_symptom_mesh_descriptors(Path(__file__).resolve().parents[1] / "docs" / "symptom_crosswalk.csv")
    by_symptom: dict[str, list[SymptomMeshDescriptor]] = {}
    for entry in entries:
        by_symptom.setdefault(entry.target_symptom, []).append(entry)
    assert ("D001007", "symptom") in [(entry.mesh_descriptor, entry.level) for entry in by_symptom["anxiety"]]
    assert ("D001008", "diagnosis") in [(entry.mesh_descriptor, entry.level) for entry in by_symptom["anxiety"]]
    assert "psychomotor_retardation" not in by_symptom
    for symptom, symptom_entries in by_symptom.items():
        assert any(entry.level == "symptom" for entry in symptom_entries), symptom


# gene identifier precedence


def _human_gem_fixture() -> pd.DataFrame:
    return pd.DataFrame({"ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3"], "gene_symbol": ["DPM1", "OLDNAME", "NOENTREZ"], "ncbi_gene_id": pd.array([8813, 2519, None], dtype="Int64")})


def _hgnc_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "symbol": ["DPM1", "FUCA2", "GBA1", "NOENTREZ", "AMBIG1", "AMBIG2"],
            "ncbi_gene_id": pd.array([999, 2519, 2629, 42, 1, 2], dtype="Int64"),
            "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG4", "ENSG3", "ENSG5", "ENSG6"],
            "prev_symbols": [[], ["OLDNAME"], ["GBA"], [], ["SHARED"], ["SHARED"]],
            "alias_symbols": [[], [], ["GCB"], [], ["ALIASX"], ["ALIASX"]],
        }
    )


def test_gene_id_mapping_precedence_human_gem_then_hgnc_symbol_previous_alias() -> None:
    nodes = pd.DataFrame(
        {
            "gene_symbol": ["DPM1", "OLDNAME", "FUCA2", "GBA", "GCB", "SHARED", "ALIASX", "NOENTREZ", "UNKNOWN"],
            "ensembl_gene_id": ["ENSG1", None, None, None, None, None, None, "ENSG3", None],
        }
    )
    mapped = map_gene_nodes(nodes, _human_gem_fixture(), _hgnc_fixture()).set_index("gene_symbol")
    assert mapped.loc["DPM1", "ncbi_gene_id"] == 8813 and mapped.loc["DPM1", "mapping_source"] == "human_gem_ensembl"  # Human-GEM beats HGNC's 999
    assert mapped.loc["OLDNAME", "ncbi_gene_id"] == 2519 and mapped.loc["OLDNAME", "mapping_source"] == "human_gem_symbol"
    assert mapped.loc["OLDNAME", "ensembl_gene_id"] == "ENSG2"  # Ensembl id filled from the matched row
    assert mapped.loc["FUCA2", "mapping_source"] == "hgnc_symbol"
    assert mapped.loc["GBA", "ncbi_gene_id"] == 2629 and mapped.loc["GBA", "mapping_source"] == "hgnc_previous_symbol" and mapped.loc["GBA", "mapped_symbol"] == "GBA1"
    assert mapped.loc["GCB", "mapping_source"] == "hgnc_alias_symbol"
    assert mapped.loc["SHARED", "mapping_source"] == "unmapped"  # previous symbol of two genes: ambiguous
    assert mapped.loc["ALIASX", "mapping_source"] == "unmapped"
    assert mapped.loc["NOENTREZ", "ncbi_gene_id"] == 42 and mapped.loc["NOENTREZ", "mapping_source"] == "hgnc_symbol"  # Human-GEM row without an id falls through
    assert mapped.loc["UNKNOWN", "mapping_source"] == "unmapped" and pd.isna(mapped.loc["UNKNOWN", "ncbi_gene_id"])
    lookup = ncbi_gene_id_to_symbol(mapped.reset_index())
    assert lookup[8813] == "DPM1" and lookup[2629] == "GBA" and 999 not in lookup


# chemical name matching


def test_autocomplete_match_prefers_exact_name_then_matched_synonym_and_skips_other_biotypes() -> None:
    candidates = [
        {"_id": "@GENE_FLUOXETINE", "biotype": "gene", "db": "ncbi_gene", "db_id": "1", "name": "fluoxetine", "match": "Matched on name <m>fluoxetine</m>"},
        {"_id": "@CHEMICAL_olanzapine_fluoxetine_combination", "biotype": "chemical", "db": "ncbi_mesh", "db_id": "C492572", "name": "olanzapine-fluoxetine combination", "match": "Multiple matches"},
        {"_id": "@CHEMICAL_Fluoxetine", "biotype": "chemical", "db": "ncbi_mesh", "db_id": "D005473", "name": "Fluoxetine", "match": "Matched on name <m>Fluoxetine</m>"},
    ]
    match = match_autocomplete_candidates("CID100003386", "fluoxetine", candidates)
    assert match == FLUOXETINE
    synonym_candidates = [{"biotype": "chemical", "db": "ncbi_mesh", "db_id": "D005680", "name": "gamma-Aminobutyric Acid", "match": "Matched on synonyms <m>GABA</m>"}]
    synonym_match = match_autocomplete_candidates("CID100000119", "GABA", synonym_candidates)
    assert synonym_match is not None and synonym_match.match_method == "autocomplete_synonym" and synonym_match.chemical_mesh_id == "D005680"
    assert match_autocomplete_candidates("CID100000119", "gamma-aminobutyric", synonym_candidates) is None


def test_mesh_label_and_ctd_name_matching_are_case_insensitive_and_record_their_method() -> None:
    names_index = drug_names_by_normalized_name({"CID100003386": "fluoxetine", "CID100002118": "Amitriptyline"})
    label_hits = match_mesh_label("MESH:D005473", "FLUOXETINE", names_index)
    assert [(hit.stitch_flat_id, hit.chemical_mesh_id, hit.match_method) for hit in label_hits] == [("CID100003386", "D005473", "mesh_label_lookup")]
    ctd_hits = match_ctd_chemical_name("D000639", "amitriptyline", names_index)
    assert [(hit.stitch_flat_id, hit.match_method) for hit in ctd_hits] == [("CID100002118", "ctd_chemical_name")]
    assert match_ctd_chemical_name("D000001", "unknown drug", names_index) == []


# CTD evidence mapping


CTD_HEADER = "# Fields:\n# ChemicalName\tChemicalID\tCasRN\tDiseaseName\tDiseaseID\tDirectEvidence\tInferenceGeneSymbol\tInferenceScore\tOmimIDs\tPubMedIDs\n#\n"


def test_ctd_direct_evidence_maps_marker_to_induces_and_therapeutic_to_relieves_one_report_per_paper() -> None:
    assert DIRECT_EVIDENCE_TO_RELATION == {"marker/mechanism": "induces", "therapeutic": "relieves"}
    text = CTD_HEADER + "\n".join(
        [
            "Fluoxetine\tD005473\t54910-89-3\tAnxiety\tMESH:D001007\ttherapeutic\t\t\t\t111|222",
            "Fluoxetine\tD005473\t54910-89-3\tAnxiety Disorders\tMESH:D001008\tmarker/mechanism|therapeutic\t\t\t\t333",
            "Fluoxetine\tD005473\t54910-89-3\tFatigue\tMESH:D005221\tmarker/mechanism\t\t\t\t",
            "Fluoxetine\tD005473\t54910-89-3\tAnxiety\tMESH:D001007\t\tSLC6A4\t4.1\t\t444",
            "Unknownium\tC999999\t\tAnxiety\tMESH:D001007\ttherapeutic\t\t\t\t555",
            "Fluoxetine\tD005473\t54910-89-3\tBreast Neoplasms\tMESH:D001943\ttherapeutic\t\t\t\t666",
        ]
    ) + "\n"
    names_index = drug_names_by_normalized_name({"CID100003386": "fluoxetine"})
    kept, counts = filter_ctd_rows(iterate_ctd_rows(io.StringIO(text)), FIXTURE_DESCRIPTORS, names_index)
    assert counts["ctd_rows_total"] == 6
    assert counts["dropped_disease_not_target_descriptor"] == 1
    assert counts["dropped_chemical_not_matched_to_sider"] == 1
    assert counts["dropped_inferred_no_direct_evidence"] == 1
    assert counts["direct_evidence_rows_without_pubmed_id"] == 1
    assert counts["rows_with_both_direct_evidence_values"] == 2
    assert [(row.relation, row.pmid, row.symptom) for row in kept] == [
        ("relieves", "111", "anxiety"),
        ("relieves", "222", "anxiety"),
        ("induces", "333", "anxiety"),
        ("relieves", "333", "anxiety"),
        ("induces", None, "fatigue"),
    ]
    assert all(row.chemical_match_method == "ctd_chemical_name" for row in kept)
    reports = ctd_rows_to_reports(kept, {"111": "2001-05-01"})
    assert reports.report_id.is_unique
    assert reports.loc[0, "evidence_date"] == "2001-05-01" and reports.loc[0, "publication_year"] == 2001
    assert reports.loc[1, "evidence_date"] is None and "undated" in reports.loc[1, "limitations"]
    assert reports.loc[4, "references"] == "" and reports.loc[4, "pubmed_reference_count"] == 0 and "no PubMed reference" in reports.loc[4, "limitations"]
    assert set(reports.evidence_code) == {"ctd_direct_evidence_therapeutic", "ctd_direct_evidence_marker_mechanism"}
    assert reports.source.eq("CTD").all() and reports.evidence_class.eq("literature").all() and reports.perturbation_type.eq("drug").all()


def test_ctd_rows_fall_back_to_the_pubtator_mesh_map_when_the_name_differs() -> None:
    text = CTD_HEADER + "Fluoxetine Hydrochloride\tD005473\t\tAnxiety\tMESH:D001007\ttherapeutic\t\t\t\t111\n"
    names_index = drug_names_by_normalized_name({"CID100003386": "fluoxetine"})
    kept, _ = filter_ctd_rows(iterate_ctd_rows(io.StringIO(text)), FIXTURE_DESCRIPTORS, names_index, chemical_matches_by_mesh_id([FLUOXETINE]))
    assert len(kept) == 1 and kept[0].chemical_match_method == "autocomplete_name"


def test_expand_ctd_row_drops_unknown_direct_evidence_values_with_a_count() -> None:
    counts: Counter = Counter()
    row = {"ChemicalName": "x", "ChemicalID": "D1", "CasRN": "", "DiseaseName": "Anxiety", "DiseaseID": "MESH:D001007", "DirectEvidence": "unknown_value", "InferenceGeneSymbol": "", "InferenceScore": "", "OmimIDs": "", "PubMedIDs": "1"}
    assert expand_ctd_row(row, FIXTURE_DESCRIPTORS["D001007"], [FLUOXETINE], counts) == []
    assert counts["dropped_direct_evidence_not_mapped:unknown_value"] == 1


# report schema


def test_pubtator_reports_follow_the_evidence_reports_schema_and_are_soft_priors() -> None:
    rows = [("7", "cause", "Gene|8813", "Disease|MESH:D001007"), ("7", "associate", "Gene|8813", "Disease|MESH:D001007")]
    kept, _ = filter_relation_rows(rows, FIXTURE_DESCRIPTORS, {8813: "DPM1"}, {})
    reports = relation_rows_to_reports(kept, {"7": "2015 Mar"} and {"7": parse_pubdate("2015 Mar")})
    assert tuple(reports.columns) == LITERATURE_REPORT_COLUMNS
    assert tuple(reports.columns[: len(EVIDENCE_REPORT_COLUMNS)]) == EVIDENCE_REPORT_COLUMNS
    assert reports.report_id.tolist() == ["PubTator3|PMID:7|MESH:D001007|anxiety|1", "PubTator3|PMID:7|MESH:D001007|anxiety|2"]
    assert reports.evidence_class.eq("literature").all() and reports.report_value.eq(1).all()
    assert reports.evidence_date.tolist() == ["2015-03-01", "2015-03-01"] and reports.publication_year.tolist() == [2015, 2015]
    assert json.loads(reports.loc[0, "perturbation_nodes"]) == [["GENE:DPM1", -1.0, 1.0]]
    assert reports.loc[1, "relation"] == "associated_with" and "undirected association" in reports.loc[1, "limitations"]
    assert (reports[["rubric_evidence_code_pcs", "rubric_evidence_code_tas", "rubric_evidence_code_iea", "rubric_placebo_controlled", "rubric_curated_synopsis"]] == 0.0).all().all()
    summary = summarize_reports(reports)
    assert summary["rows"] == 2 and summary["distinct_pmids"] == 1 and summary["rows_by_relation"] == {"associated_with": 1, "induces": 1}


def test_report_columns_match_the_shared_report_dataclass_when_it_is_present() -> None:
    evidence_reports = pytest.importorskip("mechanistic_pathway_learning.evidence.evidence_reports")
    assert tuple(evidence_reports.REPORT_COLUMNS) == EVIDENCE_REPORT_COLUMNS


def test_pubdate_parsing_defaults_missing_month_and_day_to_one() -> None:
    assert parse_pubdate("2022 Mar 29") == "2022-03-29"
    assert parse_pubdate("2021 Jan-Feb") == "2021-01-01"
    assert parse_pubdate("2020") == "2020-01-01"
    assert parse_pubdate("2019 Winter") == "2019-01-01"
    assert parse_pubdate(None) is None and parse_pubdate("n.d.") is None
