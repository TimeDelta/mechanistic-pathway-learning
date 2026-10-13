"""A drug's targets per subunit from the Guide to PHARMACOLOGY, in place of ChEMBL's complex targets (the owner,
13 October 2026: "agreed on both 1 and 2")."""
import pandas as pd
import pytest

from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import (
    GTOPDB_ACTION_SIGN,
    MODULATES_THE_ENDOGENOUS_LIGAND,
    ligand_id_of_drug,
    one_gene_per_row,
    pubchem_cid_from_perturbation_id,
    sign_of_action,
    targets_per_drug,
)
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import ACTION_TYPE_SIGN

VERSION_LINE = "# GtoPdb Version: 2026.3 - published: 2026-09-16\n"


def write_gtopdb(directory, interaction_rows, mapping_rows):
    directory.mkdir(parents=True, exist_ok=True)
    for name, rows in (("interactions.csv", interaction_rows), ("ligand_id_mapping.csv", mapping_rows)):
        path = directory / name
        path.write_text(VERSION_LINE + pd.DataFrame(rows).to_csv(index=False))
    return directory


def test_a_stitch_flat_identifier_gives_its_pubchem_cid():
    assert pubchem_cid_from_perturbation_id("CID100000444") == 444  # bupropion
    assert pubchem_cid_from_perturbation_id("CID100000896") == 896  # melatonin
    assert pubchem_cid_from_perturbation_id("RXCUI:10239") is None
    assert pubchem_cid_from_perturbation_id("CID10000044") is None  # not eight digits


def test_every_action_sign_agrees_with_the_chembl_table_where_both_name_the_action():
    """The two routes must put a drug on one scale, or a drug's input would depend on which source described it."""
    shared = {"AGONIST", "PARTIAL AGONIST", "ANTAGONIST", "INVERSE AGONIST", "BINDING AGENT"}
    for action in shared & set(GTOPDB_ACTION_SIGN) & set(ACTION_TYPE_SIGN):
        assert GTOPDB_ACTION_SIGN[action] == ACTION_TYPE_SIGN[action], action
    assert sign_of_action("Full agonist") == 1.0
    assert sign_of_action("Partial agonist") == 0.5 == GTOPDB_ACTION_SIGN["POSITIVE"]
    assert sign_of_action("Negative") == -0.5
    assert sign_of_action("Pore blocker") == sign_of_action("Antagonist") == -1.0
    assert sign_of_action("Binding") == 0.0
    assert sign_of_action(None) is None and sign_of_action(float("nan")) is None
    assert sign_of_action("Nonsense") is None  # reported, never taken as zero


def test_the_modulating_actions_are_the_ones_that_need_an_endogenous_ligand():
    assert set(MODULATES_THE_ENDOGENOUS_LIGAND) <= set(GTOPDB_ACTION_SIGN)
    for action in MODULATES_THE_ENDOGENOUS_LIGAND:
        assert abs(GTOPDB_ACTION_SIGN[action]) == 0.5  # a half-magnitude today, which is what item 24 replaces


def test_a_drug_is_matched_by_its_identifier_first_and_its_name_second(tmp_path):
    directory = write_gtopdb(tmp_path / "gtopdb", [], [
        {"Ligand id": 11, "Name": "bupropion", "PubChem CID": "444"},
        {"Ligand id": 22, "Name": "oxprenolol", "PubChem CID": ""},
        {"Ligand id": 33, "Name": "Bupropion the second", "PubChem CID": "999"},
    ])
    ligand_of_drug, route = ligand_id_of_drug({"CID100000444": 444}, {"CID100000444": "bupropion", "RXCUI:1": "Oxprenolol",
                                                                      "RXCUI:2": "not a ligand"}, directory)
    assert ligand_of_drug == {"CID100000444": 11, "RXCUI:1": 22}
    assert route == {"CID100000444": "pubchem_cid", "RXCUI:1": "name"}  # the identifier wins where both could match


def test_only_human_rows_naming_one_gene_with_a_known_action_can_replace_a_complex(tmp_path):
    directory = write_gtopdb(tmp_path / "gtopdb", [
        {"Ligand ID": 11, "Target Species": "Human", "Target Gene Symbol": "GABRA1", "Action": "Positive",
         "Type": "Allosteric modulator", "Primary Target": "t"},
        {"Ligand ID": 11, "Target Species": "Human", "Target Gene Symbol": "GABRA2", "Action": "Positive",
         "Type": "Allosteric modulator", "Primary Target": "f"},
        {"Ligand ID": 11, "Target Species": "Rat", "Target Gene Symbol": "GABRA3", "Action": "Positive",
         "Type": "Allosteric modulator", "Primary Target": "t"},
        {"Ligand ID": 11, "Target Species": "Human", "Target Gene Symbol": "GABRA1|GABRB2|GABRG2", "Action": "Positive",
         "Type": "Allosteric modulator", "Primary Target": "t"},
        {"Ligand ID": 11, "Target Species": "Human", "Target Gene Symbol": "GABRB3", "Action": "Nonsense",
         "Type": "Agonist", "Primary Target": "t"},
        {"Ligand ID": 11, "Target Species": "Human", "Target Gene Symbol": "", "Action": "Antagonist",
         "Type": "Antagonist", "Primary Target": "t"},
    ], [{"Ligand id": 11, "Name": "a drug", "PubChem CID": "444"}])

    per_drug = targets_per_drug({"CID100000444": 11}, gtopdb_directory=directory)
    rows = per_drug["CID100000444"]
    assert len(rows) == 5  # the rat row is gone, the rest are kept so a caller can count them
    usable = one_gene_per_row(rows)
    assert sorted(usable.gene_symbol) == ["GABRA1", "GABRA2"]  # the subunit list is not expanded, which is the point
    assert set(usable.sign) == {0.5} and bool(usable.modulates_the_endogenous_ligand.all())
    primary = one_gene_per_row(targets_per_drug({"CID100000444": 11}, gtopdb_directory=directory, primary_only=True)["CID100000444"])
    assert sorted(primary.gene_symbol) == ["GABRA1"]


def test_a_drug_with_no_human_row_is_absent_rather_than_empty(tmp_path):
    directory = write_gtopdb(tmp_path / "gtopdb", [
        {"Ligand ID": 11, "Target Species": "Rat", "Target Gene Symbol": "GABRA1", "Action": "Positive",
         "Type": "Allosteric modulator", "Primary Target": "t"},
    ], [{"Ligand id": 11, "Name": "a drug", "PubChem CID": "444"}])
    assert targets_per_drug({"CID100000444": 11}, gtopdb_directory=directory) == {}


def test_the_release_on_disk_names_one_gene_for_nearly_every_human_row():
    """A guard on the real file, which is what the re-derivation rests on: if a release changed the shape of
    Target Gene Symbol, the complex expansion would come back without anything else failing."""
    from mechanistic_pathway_learning.perturbation.gtopdb_drug_targets import GTOPDB_DIRECTORY, read_interactions
    if not (GTOPDB_DIRECTORY / "interactions.csv").exists():
        pytest.skip("the GtoPdb release is not on disk")
    human = read_interactions()
    named = human[human["Target Gene Symbol"].notna()]
    one_gene = ~named["Target Gene Symbol"].astype(str).str.contains(r"[|,;]", regex=True, na=False)
    assert float(one_gene.mean()) > 0.95, "a GtoPdb release has started naming subunit lists where it named one gene"
    assert set(human.Action.dropna().astype(str).str.upper()) <= set(GTOPDB_ACTION_SIGN), "a new Action has no sign"
