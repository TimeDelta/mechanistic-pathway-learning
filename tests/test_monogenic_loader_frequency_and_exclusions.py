"""Offline tests for the HPO loader: term expansion with exclusions, frequency parsing, provenance, grade policies and disease clusters."""
from pathlib import Path

from mechanistic_pathway_learning.evidence.assemble_evidence_table import disease_cluster_ids
from mechanistic_pathway_learning.evidence.assign_evidence_grades import assign_evidence_grade, loss_weight_for_record
from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    expand_symptom_terms,
    load_hpo_is_a_parents_from_obo,
    load_hpo_term_names_from_obo,
    monogenic_evidence_records,
    parse_frequency_denominator,
    parse_frequency_qualifier,
    parse_genes_to_phenotype,
    read_crosswalk_hpo_terms,
)

SYNTHETIC_OBO = """format-version: 1.2

[Term]
id: HP:0100543
name: Cognitive impairment

[Term]
id: HP:0007018
name: Attention deficit hyperactivity disorder
is_a: HP:0100543 ! Cognitive impairment

[Term]
id: HP:0000736
name: Short attention span
is_a: HP:0007018 ! Attention deficit hyperactivity disorder
is_a: HP:0100543 ! Cognitive impairment

[Term]
id: HP:0000726
name: Dementia
is_a: HP:0100543 ! Cognitive impairment

[Term]
id: HP:0000709
name: Psychosis
"""

SYNTHETIC_ANNOTATIONS = "\n".join([
    "ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id",
    "1\tHMBS\tHP:0000709\tPsychosis\t-\tOMIM:176000",
    "1\tHMBS\tHP:0000709\tPsychosis\tHP:0040283\tORPHA:79276",
    "2\tOTC\tHP:0000726\tDementia\t3/5\tOMIM:311250",
    "3\tADHDGENE\tHP:0007018\tADHD\tHP:0040282\tORPHA:1",
    "4\tSHAREDGENE\tHP:0000736\tShort attention span\t25%\tORPHA:1",
    "4\tSHAREDGENE\tHP:0000726\tDementia\tHP:0040285\tORPHA:2",
    "5\tNOTINGRAPH\tHP:0000709\tPsychosis\t-\tOMIM:1",
]) + "\n"

CROSSWALK = "\n".join([
    "target_symptom,hpo_ids,hpo_names,excluded_hpo_ids,meddra_preferred_terms",
    "cognitive_impairment,HP:0100543,Cognitive impairment,HP:0007018,Memory impairment",
    "psychosis,HP:0000709,Psychosis,,Psychotic disorder",
]) + "\n"


def write_inputs(directory: Path) -> tuple[Path, Path, Path]:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "hp.obo").write_text(SYNTHETIC_OBO)
    (directory / "genes_to_phenotype.txt").write_text(SYNTHETIC_ANNOTATIONS)
    (directory / "crosswalk.csv").write_text(CROSSWALK)
    return directory / "hp.obo", directory / "genes_to_phenotype.txt", directory / "crosswalk.csv"


def test_frequency_parsing() -> None:
    assert parse_frequency_qualifier("HP:0040281") == 0.895
    assert abs(parse_frequency_qualifier("3/5") - 0.6) < 1e-9
    assert abs(parse_frequency_qualifier("25%") - 0.25) < 1e-9
    assert parse_frequency_qualifier("-") is None and parse_frequency_qualifier("") is None and parse_frequency_qualifier("0/0") is None
    assert parse_frequency_denominator("3/5") == 5 and parse_frequency_denominator("HP:0040282") is None


def test_exclusion_removes_subtree_but_keeps_terms_reachable_elsewhere(tmp_path: Path) -> None:
    obo_path, _, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    names = load_hpo_term_names_from_obo(obo_path)
    assert names["HP:0000726"] == "Dementia"
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    assert crosswalk["cognitive_impairment"] == (["HP:0100543"], ["HP:0007018"])
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    excluded = {symptom: terms[1] for symptom, terms in crosswalk.items()}
    expanded = expand_symptom_terms(roots, parents, excluded)
    assert "HP:0007018" not in expanded  # excluded diagnosis-level term
    assert "HP:0000736" not in expanded  # its descendant goes with it even though it also hangs under the root
    assert expanded["HP:0000726"] == ["cognitive_impairment"] and expanded["HP:0000709"] == ["psychosis"]
    without_exclusion = expand_symptom_terms(roots, parents)
    assert "HP:0007018" in without_exclusion and "HP:0000736" in without_exclusion


def test_records_carry_frequency_provenance_and_drop_excluded_rows(tmp_path: Path) -> None:
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    excluded = {symptom: terms[1] for symptom, terms in crosswalk.items()}
    records = monogenic_evidence_records(parse_genes_to_phenotype(annotations_path), roots, parents, {"HMBS", "OTC", "ADHDGENE", "SHAREDGENE"}, excluded)
    by_pair = {(record.perturbation_identifier, record.symptom_identifier): record for record in records}
    hmbs = by_pair[("HMBS", "psychosis")]
    assert hmbs.omim_entry_count == 1 and hmbs.orpha_entry_count == 1 and hmbs.annotation_row_count == 2
    assert abs(hmbs.max_annotation_frequency - 0.17) < 1e-9 and hmbs.has_omim_clinical_synopsis
    otc = by_pair[("OTC", "cognitive_impairment")]
    assert abs(otc.max_annotation_frequency - 0.6) < 1e-9 and otc.annotation_patient_count == 5 and otc.orpha_entry_count == 0
    assert ("ADHDGENE", "cognitive_impairment") not in by_pair  # only annotated through the excluded term
    assert ("SHAREDGENE", "cognitive_impairment") not in by_pair  # Short attention span excluded; Dementia row is an Excluded qualifier
    assert ("NOTINGRAPH", "psychosis") not in by_pair


def test_grade_policies_and_frequency_scaled_weights(tmp_path: Path) -> None:
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    records = {(r.perturbation_identifier, r.symptom_identifier): r for r in monogenic_evidence_records(parse_genes_to_phenotype(annotations_path), roots, parents, {"HMBS", "OTC"})}
    hmbs, otc = records[("HMBS", "psychosis")], records[("OTC", "cognitive_impairment")]
    # provenance policy: both curated entries are grade A; weights follow the largest frequency
    assert assign_evidence_grade(hmbs) == "A" and assign_evidence_grade(otc) == "A"
    assert abs(loss_weight_for_record(hmbs) - 0.68) < 1e-9  # Occasional midpoint 0.17 * 4
    assert loss_weight_for_record(otc) == 1.0  # 3/5 = 0.6 -> capped at full weight
    assert loss_weight_for_record(otc, scale_by_frequency=False) == 1.0
    # version 0.3 proxy: OMIM entry or two entries -> A, else C with zero weight
    assert assign_evidence_grade(otc, grade_a_policy="two_distinct_disease_entries") == "A"
    orpha_only = type(otc)("GENE", "psychosis", "induces", "monogenic", orpha_entry_count=1, independent_case_series_count=1)
    assert assign_evidence_grade(orpha_only) == "A" and assign_evidence_grade(orpha_only, grade_a_policy="two_distinct_disease_entries") == "C"
    assert loss_weight_for_record(orpha_only, grade_a_policy="two_distinct_disease_entries") == 0.0


def test_disease_clusters_join_genes_sharing_a_disease(tmp_path: Path) -> None:
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    records = monogenic_evidence_records(parse_genes_to_phenotype(annotations_path), roots, parents, {"HMBS", "OTC", "ADHDGENE", "SHAREDGENE"})
    clusters = disease_cluster_ids(records)
    assert clusters["ADHDGENE"] == clusters["SHAREDGENE"]  # both annotated to ORPHA:1
    assert clusters["HMBS"] != clusters["OTC"] and clusters["HMBS"] != clusters["ADHDGENE"]


def test_annotation_dates_use_omim_rows_only_and_the_earliest_curation(tmp_path: Path) -> None:
    from datetime import date

    from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import load_hpo_annotation_dates

    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text("\n".join([
        "#description: test",
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration",
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1\tPCS\t\t\t\t\tP\tHPO:skoehler[2012-03-04];HPO:probinson[2021-06-21]",
        "OMIM:176000\tAIP\t\tHP:0000709\tOMIM:176000\tIEA\t\t\t\t\tP\tHPO:iea[2009-02-17]",
        "OMIM:176000\tAIP\tNOT\tHP:0000726\tPMID:2\tPCS\t\t\t\t\tP\tHPO:skoehler[2010-01-01]",
        "ORPHA:79276\tAIP\t\tHP:0000709\tORPHA:79276\tTAS\t\t\t\t\tP\tORPHA:orphadata[2026-09-02]",
    ]) + "\n")
    dates = load_hpo_annotation_dates(hpoa)
    assert dates[("OMIM:176000", "HP:0000709")] == date(2009, 2, 17)  # earliest curation across rows and curators
    assert ("OMIM:176000", "HP:0000726") not in dates  # NOT-qualified rows are skipped
    assert ("ORPHA:79276", "HP:0000709") not in dates  # Orphanet rows carry the import date and are excluded
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    records = {r.perturbation_identifier: r for r in monogenic_evidence_records(parse_genes_to_phenotype(annotations_path), roots, parents, {"HMBS", "OTC"}, annotation_dates=dates)}
    assert records["HMBS"].evidence_available_date == date(2009, 2, 17)
    assert records["OTC"].evidence_available_date is None


def test_publication_dates_take_precedence_over_curation_dates(tmp_path: Path) -> None:
    import json
    from datetime import date

    from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import load_hpo_annotation_dates, load_reference_publication_dates

    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text("\n".join([
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration",
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:123;PMID:456\tPCS\t\t\t\t\tP\tHPO:skoehler[2012-03-04]",
        "OMIM:176000\tAIP\t\tHP:0000726\tOMIM:176000\tTAS\t\t\t\t\tP\tHPO:skoehler[2010-05-06]",
    ]) + "\n")
    lookup = tmp_path / "dates.json"
    lookup.write_text(json.dumps({"publication_dates_by_pmid": {"123": {"year": 1998, "month": 7, "day": None}, "456": {"year": 2005, "month": None}}}))
    publication_dates = load_reference_publication_dates(lookup)
    assert publication_dates["123"] == date(1998, 7, 1) and publication_dates["456"] == date(2005, 1, 1)
    dates = load_hpo_annotation_dates(hpoa, publication_dates_by_pmid=publication_dates)
    assert dates[("OMIM:176000", "HP:0000709")] == date(1998, 7, 1)  # earliest cited publication
    assert dates[("OMIM:176000", "HP:0000726")] == date(2010, 5, 6)  # no PubMed reference: curation date stands


def test_reference_counts_are_descriptive_and_do_not_change_weights(tmp_path: Path) -> None:
    from mechanistic_pathway_learning.evidence.assign_evidence_grades import loss_weight_for_record
    from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import load_hpo_annotation_references

    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text("\n".join([
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration",
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1;OMIM:176000\tPCS\t\t\t\t\tP\tHPO:a[2012-03-04]",
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1\tPCS\t\t\t\t\tP\tHPO:b[2013-03-04]",
        "OMIM:176000\tAIP\tNOT\tHP:0000709\tPMID:99\tPCS\t\t\t\t\tP\tHPO:a[2013-03-04]",
        "ORPHA:79276\tAIP\t\tHP:0000709\tORPHA:79276\tTAS\t\t\t\t\tP\tORPHA:orphadata[2026-09-02]",
    ]) + "\n")
    references = load_hpo_annotation_references(hpoa)
    assert references[("OMIM:176000", "HP:0000709")] == {"PMID:1", "OMIM:176000"}  # NOT-qualified row skipped, duplicates collapsed
    assert references[("ORPHA:79276", "HP:0000709")] == {"ORPHA:79276"}
    obo_path, annotations_path, crosswalk_path = write_inputs(tmp_path)
    parents = load_hpo_is_a_parents_from_obo(obo_path)
    crosswalk = read_crosswalk_hpo_terms(crosswalk_path)
    roots = {symptom: terms[0] for symptom, terms in crosswalk.items()}
    rows = parse_genes_to_phenotype(annotations_path)
    with_references = {r.perturbation_identifier: r for r in monogenic_evidence_records(rows, roots, parents, {"HMBS", "OTC"}, annotation_references=references)}
    without_references = {r.perturbation_identifier: r for r in monogenic_evidence_records(rows, roots, parents, {"HMBS", "OTC"})}
    assert with_references["HMBS"].distinct_reference_count == 3 and with_references["HMBS"].distinct_pubmed_reference_count == 1  # PMID:1, OMIM:176000, ORPHA:79276 across its two disease entries
    assert with_references["OTC"].distinct_reference_count == 0  # no phenotype.hpoa row for its disease
    for gene in ("HMBS", "OTC"):
        assert loss_weight_for_record(with_references[gene]) == loss_weight_for_record(without_references[gene])


def test_annotation_rows_reader_matches_dates_and_references(tmp_path: Path) -> None:
    from datetime import date

    from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
        hpoa_row_availability_date,
        load_hpo_annotation_dates,
        load_hpo_annotation_references,
        load_hpo_annotation_rows,
    )

    hpoa = tmp_path / "phenotype.hpoa"
    hpoa.write_text("\n".join([
        "#description: test",
        "database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration",
        "OMIM:176000\tAIP\t\tHP:0000709\tPMID:1;OMIM:176000\tPCS\t\t3/5\t\t\tP\tHPO:skoehler[2012-03-04];HPO:probinson[2021-06-21]",
        "OMIM:176000\tAIP\t\tHP:0000709\tOMIM:176000\tIEA\t\t\t\t\tP\tHPO:iea[2009-02-17]",
        "OMIM:176000\tAIP\tNOT\tHP:0000726\tPMID:2\tPCS\t\t\t\t\tP\tHPO:skoehler[2010-01-01]",
        "ORPHA:79276\tAIP\t\tHP:0000709\tORPHA:79276\tTAS\t\t\t\t\tP\tORPHA:orphadata[2026-09-02]",
    ]) + "\n")
    rows_by_key = load_hpo_annotation_rows(hpoa)
    omim_rows = rows_by_key[("OMIM:176000", "HP:0000709")]
    assert [row.ordinal for row in omim_rows] == [1, 2] and omim_rows[0].references == ["PMID:1", "OMIM:176000"] and omim_rows[0].frequency == "3/5"
    assert omim_rows[0].biocuration_dates == [date(2012, 3, 4), date(2021, 6, 21)] and omim_rows[0].evidence == "PCS" and omim_rows[0].disease_name == "AIP"
    assert rows_by_key[("OMIM:176000", "HP:0000726")][0].qualifier == "NOT" and rows_by_key[("OMIM:176000", "HP:0000726")][0].is_negated
    assert hpoa_row_availability_date(omim_rows[0], publication_dates_by_pmid={"1": date(1990, 1, 1)}) == date(1990, 1, 1)
    assert hpoa_row_availability_date(rows_by_key[("ORPHA:79276", "HP:0000709")][0]) is None
    # the two older readers, now wrappers over the row reader, give the dictionaries they always gave
    assert load_hpo_annotation_dates(hpoa) == {("OMIM:176000", "HP:0000709"): date(2009, 2, 17)}
    assert load_hpo_annotation_dates(hpoa, publication_dates_by_pmid={"1": date(1990, 1, 1)}) == {("OMIM:176000", "HP:0000709"): date(1990, 1, 1)}
    assert load_hpo_annotation_references(hpoa) == {("OMIM:176000", "HP:0000709"): {"PMID:1", "OMIM:176000"}, ("ORPHA:79276", "HP:0000709"): {"ORPHA:79276"}}
