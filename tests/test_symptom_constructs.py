"""The fine symptom constructs derived from the coarse crosswalk (evidence/symptom_constructs.py, decision 9.58)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.build_fine_symptom_crosswalk import complete_addition, render  # noqa: E402
from mechanistic_pathway_learning.evidence.symptom_constructs import (  # noqa: E402
    FINE_COLUMNS,
    IMPC_COLUMNS,
    HpoVocabulary,
    assign_term,
    build_fine_crosswalk,
    coarse_of_fine,
    fine_constructs_of,
    fine_impc_rows,
    phrase,
    read_rows,
    read_term_overrides,
    slug,
    stems,
    write_rows,
)

OBO = """format-version: 1.2

[Term]
id: HP:0000001
name: All

[Term]
id: HP:0000709
name: Psychosis
synonym: "Psychotic episodes" EXACT []
is_a: HP:0000001 ! All

[Term]
id: HP:0000738
name: Hallucinations
synonym: "Hallucination" EXACT layperson []
is_a: HP:0000709 ! Psychosis

[Term]
id: HP:0000746
name: Delusion
synonym: "Delusions" EXACT layperson []
is_a: HP:0000709 ! Psychosis

[Term]
id: HP:0008765
name: Visual hallucinations
synonym: "Hallucination, visual" EXACT []
is_a: HP:0000738 ! Hallucinations

[Term]
id: HP:0000713
name: Agitation
is_a: HP:0000001 ! All

[Term]
id: HP:0000711
name: Restlessness
synonym: "Motor restlessness" EXACT []
is_a: HP:0000001 ! All

[Term]
id: HP:0001249
name: Intellectual disability
synonym: "Low intelligence" RELATED []
is_a: HP:0000001 ! All

[Term]
id: HP:9999999
name: Obsolete thing
is_obsolete: true
"""


def coarse_rows():
    return [
        {"target_symptom": "psychosis", "hpo_ids": "HP:0000709;HP:0000738;HP:0000746", "hpo_names": "Psychosis;Hallucinations;Delusion", "excluded_hpo_ids": "",
         "meddra_preferred_terms": "Psychotic disorder;Hallucination;Delusion;Hallucination, visual;Acute psychosis", "questionnaire_item": "", "rdoc_domain": "perception",
         "identifier_status": "verified", "notes": "four roots", "mesh_descriptors": "D011618", "mesh_descriptor_level": "diagnosis"},
        {"target_symptom": "psychomotor_agitation", "hpo_ids": "HP:0000713;HP:0000711", "hpo_names": "Agitation;Restlessness", "excluded_hpo_ids": "",
         "meddra_preferred_terms": "Agitation;Restlessness;Akathisia", "questionnaire_item": "PHQ-9 item 8", "rdoc_domain": "arousal", "identifier_status": "verified",
         "notes": "", "mesh_descriptors": "", "mesh_descriptor_level": ""},
        {"target_symptom": "anxiety", "hpo_ids": "HP:0000739", "hpo_names": "Anxiety", "excluded_hpo_ids": "", "meddra_preferred_terms": "Anxiety;Nervousness",
         "questionnaire_item": "GAD-7", "rdoc_domain": "negative_valence", "identifier_status": "verified", "notes": "one root", "mesh_descriptors": "D001007",
         "mesh_descriptor_level": "symptom"},
        {"target_symptom": "catatonia", "hpo_ids": "", "hpo_names": "", "excluded_hpo_ids": "", "meddra_preferred_terms": "Catatonia", "questionnaire_item": "",
         "rdoc_domain": "", "identifier_status": "", "notes": "no HPO root", "mesh_descriptors": "", "mesh_descriptor_level": ""},
    ]


@pytest.fixture
def vocabulary(tmp_path):
    obo = tmp_path / "hp.obo"
    obo.write_text(OBO, encoding="utf-8")
    return HpoVocabulary.from_obo(obo)


def test_slugs_stems_and_phrases_normalise_spelling_plurals_and_word_order():
    assert slug("Pill-rolling tremor") == "pill_rolling_tremor" and slug("  Aggressive behavior ") == "aggressive_behavior"
    assert stems("Suicidal behaviour") == {"suicid"} and stems("Suicidal behaviour", keep_stop_words=True) == {"suicid", "behavior"}
    assert phrase("Hallucination, visual") == phrase("Visual hallucinations")


def test_a_fused_row_splits_at_its_roots_and_a_single_root_row_keeps_its_name():
    psychosis, agitation, anxiety, catatonia = (fine_constructs_of(row) for row in coarse_rows())
    assert [construct.name for construct in psychosis] == ["psychosis", "hallucinations", "delusion"]
    assert [construct.is_namesake for construct in psychosis] == [True, False, False], "the root whose name slugs to the coarse name is the namesake"
    assert [construct.name for construct in agitation] == ["agitation", "restlessness"] and agitation[0].is_namesake, "no root slugs to the coarse name, so the first is"
    assert [(construct.name, construct.hpo_id) for construct in anxiety] == [("anxiety", "HP:0000739")]
    assert [(construct.name, construct.hpo_id) for construct in catatonia] == [("catatonia", "")]
    with pytest.raises(ValueError, match="hpo_names"):
        fine_constructs_of({"target_symptom": "x", "hpo_ids": "HP:1;HP:2", "hpo_names": "One"})


def test_the_vocabulary_reads_names_synonyms_and_descendants_and_drops_obsolete_terms(vocabulary):
    assert vocabulary.names_by_term["HP:0000738"] == "Hallucinations" and "HP:9999999" not in vocabulary.parents_by_term
    phrases = vocabulary.phrases_of("HP:0000738", [])
    assert phrase("Hallucination, visual") in phrases and phrase("Visual hallucinations") in phrases and phrase("Delusion") not in phrases
    assert phrase("Visual hallucinations") not in vocabulary.phrases_of("HP:0000738", ["HP:0008765"]), "an excluded descendant leaves the construct's vocabulary"
    assert "episod" in vocabulary.root_stems_of("HP:0000709")
    assert vocabulary.ancestors_of("HP:0008765") == {"HP:0000738", "HP:0000709", "HP:0000001"} and vocabulary.ancestors_of("HP:0000001") == set()


def test_terms_are_assigned_by_override_then_vocabulary_then_stems_then_namesake(vocabulary):
    constructs = fine_constructs_of(coarse_rows()[0])
    assert assign_term("Hallucination, visual", constructs, {}, vocabulary) == ("hallucinations", "hpo vocabulary"), "a descendant's synonym, in another word order"
    assert assign_term("Delusion", constructs, {}, vocabulary) == ("delusion", "hpo vocabulary")
    assert assign_term("Acute psychosis", constructs, {}, vocabulary) == ("psychosis", "stem match")
    assert assign_term("Psychotic disorder", constructs, {}, vocabulary) == ("psychosis", "stem match"), "through the root's synonym 'Psychotic episodes'"
    assert assign_term("Unrelated term", constructs, {}, vocabulary) == ("psychosis", "namesake")
    assert assign_term("Psychotic disorder", constructs, {"psychotic disorder": "delusion"}, vocabulary) == ("delusion", "override")
    with pytest.raises(ValueError, match="not a construct"):
        assign_term("Psychotic disorder", constructs, {"psychotic disorder": "anxiety"}, vocabulary)
    assert assign_term("Anxiety", fine_constructs_of(coarse_rows()[2]), {}, vocabulary) == ("anxiety", "single construct")


def test_the_fine_crosswalk_carries_terms_apart_keeps_mesh_with_the_namesake_and_adds_constructs(vocabulary):
    additions = [complete_addition({"target_symptom": "intellectual_disability", "hpo_ids": "HP:0001249", "rdoc_domain": "cognitive_systems", "mesh_descriptors": "D008607",
                                    "mesh_descriptor_level": "diagnosis", "notes": "added"}, vocabulary, ["Intellectual disability", "Low intelligence", "Anxiety"])]
    rows, assignments = build_fine_crosswalk(coarse_rows(), additions, {"akathisia": "restlessness"}, vocabulary)
    by_name = {row["target_symptom"]: row for row in rows}
    assert list(by_name) == ["psychosis", "hallucinations", "delusion", "agitation", "restlessness", "anxiety", "catatonia", "intellectual_disability"]
    assert by_name["hallucinations"]["meddra_preferred_terms"] == "Hallucination;Hallucination, visual" and by_name["delusion"]["meddra_preferred_terms"] == "Delusion"
    assert by_name["psychosis"]["meddra_preferred_terms"] == "Psychotic disorder;Acute psychosis"
    assert by_name["psychosis"]["mesh_descriptors"] == "D011618" and by_name["hallucinations"]["mesh_descriptors"] == "", "MeSH stays with the namesake"
    assert by_name["psychosis"]["excluded_hpo_ids"] == "HP:0000738;HP:0000746", "the parent root excludes the sibling roots that descend from it"
    assert by_name["hallucinations"]["excluded_hpo_ids"] == ""
    assert by_name["agitation"]["questionnaire_item"] == "PHQ-9 item 8" and by_name["restlessness"]["questionnaire_item"] == ""
    assert by_name["restlessness"]["meddra_preferred_terms"] == "Restlessness;Akathisia", "the override moved Akathisia"
    assert by_name["anxiety"]["notes"] == "one root" and by_name["hallucinations"]["notes"].startswith("derived from psychosis")
    assert by_name["intellectual_disability"]["hpo_names"] == "Intellectual disability"
    assert by_name["intellectual_disability"]["meddra_preferred_terms"] == "Intellectual disability;Low intelligence", "the preferred terms whose words the ontology uses"
    assert coarse_of_fine(rows)["delusion"] == "psychosis" and coarse_of_fine(rows)["intellectual_disability"] == "intellectual_disability"
    assert {assignment["how"] for assignment in assignments} == {"hpo vocabulary", "override", "stem match", "single construct"}
    with pytest.raises(ValueError, match="named twice"):
        build_fine_crosswalk(coarse_rows() + [coarse_rows()[2]], [], {}, vocabulary)


def test_the_impc_crosswalk_is_rekeyed_to_the_fine_constructs():
    impc = [
        {"target_symptom": "psychosis", "relation": "induces", "mp_term_ids": "", "mp_term_names": "", "negative_procedure_codes": "", "validity": "endophenotype_only", "note": "never read"},
        {"target_symptom": "psychomotor_agitation", "relation": "induces", "mp_term_ids": "MP:0001399", "mp_term_names": "hyperactivity;restlessness", "negative_procedure_codes": "OFD", "validity": "moderate", "note": ""},
        {"target_symptom": "anxiety", "relation": "induces", "mp_term_ids": "MP:0001363", "mp_term_names": "increased anxiety-related response", "negative_procedure_codes": "OFD", "validity": "strong", "note": ""},
    ]
    rows = fine_impc_rows(impc, coarse_rows(), [{"target_symptom": "intellectual_disability", "hpo_ids": "HP:0001249"}])
    by_key = {(row["target_symptom"], row["relation"]): row for row in rows}
    assert by_key[("psychosis", "induces")]["validity"] == "endophenotype_only" and by_key[("hallucinations", "induces")]["validity"] == "none"
    assert by_key[("restlessness", "induces")]["mp_term_ids"] == "MP:0001399", "the MP term names match restlessness, not agitation"
    assert by_key[("agitation", "induces")]["validity"] == "none" and "restlessness" in by_key[("agitation", "induces")]["note"]
    assert by_key[("anxiety", "induces")]["validity"] == "strong"
    assert by_key[("intellectual_disability", "induces")]["validity"] == "none"
    with pytest.raises(ValueError, match="does not hold"):
        fine_impc_rows([{**impc[0], "target_symptom": "unknown"}], coarse_rows(), [])


def test_rows_round_trip_through_the_writer_and_the_page_lists_every_assignment(tmp_path, vocabulary):
    rows, assignments = build_fine_crosswalk(coarse_rows(), [], {}, vocabulary)
    impc_rows = fine_impc_rows([], coarse_rows(), [])
    write_rows(tmp_path / "fine.csv", rows, FINE_COLUMNS)
    write_rows(tmp_path / "impc.csv", impc_rows, IMPC_COLUMNS)
    assert read_rows(tmp_path / "fine.csv") == [{column: row.get(column, "") for column in FINE_COLUMNS} for row in rows]
    assert read_rows(tmp_path / "impc.csv") == []
    with open(tmp_path / "overrides.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["term", "fine_construct", "note"])
        writer.writerow(["Akathisia", "restlessness", "motor"])
    assert read_term_overrides(tmp_path / "overrides.csv") == {"akathisia": "restlessness"} and read_term_overrides(tmp_path / "missing.csv") == {}
    page = render("stamp", rows, assignments, impc_rows)
    assert "| hallucinations | psychosis | HP:0000738 Hallucinations |" in page and "| psychosis | Hallucination, visual | hallucinations | hpo vocabulary |" in page
    assert "7 fine constructs from 4 coarse ones" in page and "Counts:" in page


def test_the_committed_fine_crosswalk_is_a_fresh_build_of_the_committed_inputs():
    """The derived files are generated: a change to the coarse crosswalk, the additions or the overrides without a rerun fails here."""
    import subprocess
    root = Path(__file__).resolve().parents[1]
    if not (root / "data/raw/hpo/hp.obo").exists() or not (root / "data/raw/sider_4.1/meddra_all_se.tsv.gz").exists():
        pytest.skip("hp.obo and the SIDER vocabulary are not fetched")
    result = subprocess.run([sys.executable, "experiments/build_fine_symptom_crosswalk.py", "--check"], cwd=root, capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin"})
    assert result.returncode == 0, result.stdout + result.stderr
    fine = read_rows(root / "docs/symptom_crosswalk_fine.csv")
    coarse = read_rows(root / "docs/symptom_crosswalk.csv")
    assert len(fine) == sum(max(1, len([part for part in row["hpo_ids"].split(";") if part])) for row in coarse) + 1
    assert {row["coarse_symptom"] for row in fine} == {row["target_symptom"] for row in coarse} | {"intellectual_disability"}


def test_the_literature_table_is_rekeyed_to_the_namesakes_without_a_row_moving():
    import pandas as pd
    from mechanistic_pathway_learning.evidence.symptom_constructs import namesake_of_coarse
    from experiments.rekey_literature_to_fine import rekey
    namesake = namesake_of_coarse(coarse_rows())
    assert namesake == {"psychosis": "psychosis", "psychomotor_agitation": "agitation", "anxiety": "anxiety", "catatonia": "catatonia"}
    reports = pd.DataFrame({"report_id": ["PubTator3|PMID:1|MESH:D011595|psychomotor_agitation|1", "CTD|PMID:2|MESH:D001007|anxiety|1"],
                            "symptom": ["psychomotor_agitation", "anxiety"], "report_value": [1, 1]})
    renamed = rekey(reports, namesake)
    assert list(renamed.symptom) == ["agitation", "anxiety"] and list(renamed.coarse_symptom) == ["psychomotor_agitation", "anxiety"]
    assert renamed.report_id[0] == "PubTator3|PMID:1|MESH:D011595|agitation|1" and renamed.report_id[1] == reports.report_id[1]
    with pytest.raises(ValueError, match="no coarse construct"):
        rekey(pd.DataFrame({"symptom": ["unknown"], "report_value": [1]}), namesake)
