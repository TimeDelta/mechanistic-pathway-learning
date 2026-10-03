"""Offline tests of the OnSIDES identifier bridge on small fixtures: precedence, unification with collisions, OMOP handling."""
import json
from pathlib import Path

from mechanistic_pathway_learning.evidence.onsides_identifier_bridge import (
    BRIDGE_METHOD_NAME_PREF_NAME,
    BRIDGE_METHOD_NONE,
    BRIDGE_METHOD_UNII_UNICHEM,
    IngredientBridge,
    apply_bridge_precedence,
    bridge_counts,
    bridge_to_json_record,
    fallback_perturbation_id,
    identifier_system_for_ingredient,
    load_ingredient_identifier_bridge,
    mechanism_targets_for_parent,
    passes_nervous_system_rule,
    passes_single_target_rule,
    perturbation_id_by_ingredient,
    resolve_parent_from_name_route,
    resolve_parent_from_unii_route,
    rxnav_atc_classes_for_ingredient,
    rxnav_unii_codes,
    sider_parent_molecules,
    sider_pubchem_cids_by_parent,
    stitch_flat_id_from_pubchem_cid,
    unichem_chembl_ids,
    unify_with_sider,
)

# A salt (CHEMBL1200789) whose parent is CHEMBL1201317; CHEMBL703 and CHEMBL455 are their own parents.
PARENT_BY_CHEMBL_ID = {"CHEMBL1200789": "CHEMBL1201317", "CHEMBL1201317": "CHEMBL1201317", "CHEMBL703": "CHEMBL703", "CHEMBL455": "CHEMBL455"}
NAME_ROUTE_MOLECULES = [{"molecule_chembl_id": "CHEMBL1200789", "pref_name": "TRICLOFOS SODIUM", "molecule_hierarchy": {"parent_chembl_id": "CHEMBL1201317"}}]


def test_precedence_takes_the_unii_route_before_the_name_route() -> None:
    unii_route = resolve_parent_from_unii_route(["CHEMBL1200789"], PARENT_BY_CHEMBL_ID)
    assert unii_route == ("CHEMBL1201317", [])
    name_route = resolve_parent_from_name_route([{"molecule_chembl_id": "CHEMBL703", "pref_name": "SUXAMETHONIUM", "molecule_hierarchy": None}])
    assert name_route == ("CHEMBL703", [], "SUXAMETHONIUM")
    method, parent, alternatives, pref_name = apply_bridge_precedence("RXCUI", unii_route, name_route)
    assert (method, parent, alternatives) == (BRIDGE_METHOD_UNII_UNICHEM, "CHEMBL1201317", [])
    assert pref_name == "SUXAMETHONIUM"  # the name the name route saw is kept for display even when route (a) wins


def test_precedence_falls_back_to_the_name_route_and_then_to_none() -> None:
    empty_unii_route = resolve_parent_from_unii_route([], PARENT_BY_CHEMBL_ID)
    assert empty_unii_route == (None, [])
    method, parent, alternatives, _ = apply_bridge_precedence("RXCUI", empty_unii_route, resolve_parent_from_name_route(NAME_ROUTE_MOLECULES))
    assert (method, parent, alternatives) == (BRIDGE_METHOD_NAME_PREF_NAME, "CHEMBL1201317", [])
    method, parent, alternatives, pref_name = apply_bridge_precedence("RXCUI", empty_unii_route, resolve_parent_from_name_route([]))
    assert (method, parent, alternatives, pref_name) == (BRIDGE_METHOD_NONE, None, [], None)


def test_several_parents_pick_the_lowest_id_and_keep_the_rest_as_alternatives() -> None:
    parent, alternatives = resolve_parent_from_unii_route(["CHEMBL703", "CHEMBL455", "CHEMBL1200789"], PARENT_BY_CHEMBL_ID)
    assert parent == "CHEMBL455" and alternatives == ["CHEMBL703", "CHEMBL1201317"]
    unknown_parent, _ = resolve_parent_from_unii_route(["CHEMBL99999"], PARENT_BY_CHEMBL_ID)
    assert unknown_parent == "CHEMBL99999"  # an id without a hierarchy record is its own parent


def test_omop_ids_skip_the_unii_route_and_get_an_omop_perturbation_id() -> None:
    assert identifier_system_for_ingredient("OMOP997977") == "OMOP"
    assert identifier_system_for_ingredient("83367") == "RXCUI"
    assert fallback_perturbation_id("OMOP997977") == "OMOP:997977"
    assert fallback_perturbation_id("83367") == "RXCUI:83367"
    unii_route = ("CHEMBL703", [])  # a caller passing a UNII result for an OMOP id must not have it honoured
    method, parent, _, _ = apply_bridge_precedence("OMOP", unii_route, resolve_parent_from_name_route(NAME_ROUTE_MOLECULES))
    assert (method, parent) == (BRIDGE_METHOD_NAME_PREF_NAME, "CHEMBL1201317")
    method, parent, _, _ = apply_bridge_precedence("OMOP", unii_route, resolve_parent_from_name_route([]))
    assert (method, parent) == (BRIDGE_METHOD_NONE, None)


def test_unification_with_sider_picks_the_lowest_cid_and_records_collisions() -> None:
    pubchem_to_chembl_with_parents = {"5314": ["CHEMBL703"], "5320": ["CHEMBL455"], "90000": ["CHEMBL1200789", "CHEMBL1201317"], "80000": ["CHEMBL1201317"], "70000": []}
    molecule_parents = {"CHEMBL1200789": {"parent": "CHEMBL1201317", "pref_name": "TRICLOFOS SODIUM", "max_phase": None}}
    mechanisms = [{"molecule_chembl_id": "CHEMBL703", "parent_molecule_chembl_id": "CHEMBL703", "target_chembl_id": "CHEMBL1"}]
    parents_by_cid = sider_parent_molecules(pubchem_to_chembl_with_parents, molecule_parents, mechanisms)
    assert parents_by_cid == {"5314": {"CHEMBL703"}, "5320": {"CHEMBL455"}, "90000": {"CHEMBL1201317"}, "80000": {"CHEMBL1201317"}}
    cids_by_parent = sider_pubchem_cids_by_parent(parents_by_cid)
    assert cids_by_parent["CHEMBL1201317"] == [80000, 90000]
    assert unify_with_sider("10154", "CHEMBL703", cids_by_parent) == ("CID100005314", 5314, [])
    assert unify_with_sider("10156", "CHEMBL1201317", cids_by_parent) == ("CID100080000", 80000, [90000])
    assert unify_with_sider("10156", "CHEMBL999", cids_by_parent) == ("RXCUI:10156", None, [])
    assert unify_with_sider("OMOP997977", None, cids_by_parent) == ("OMOP:997977", None, [])
    assert stitch_flat_id_from_pubchem_cid(5314) == "CID100005314"


def test_response_parsers_handle_empty_and_combination_entries() -> None:
    assert rxnav_unii_codes({}) == []
    assert rxnav_unii_codes({"propConceptGroup": {"propConcept": [{"propCategory": "CODES", "propName": "UNII_CODE", "propValue": "A0JWA85V8F"}]}}) == ["A0JWA85V8F"]
    rxclass_document = {"rxclassDrugInfoList": {"rxclassDrugInfo": [
        {"minConcept": {"rxcui": "1422085", "tty": "MIN"}, "rxclassMinConceptItem": {"classId": "C10BA"}},
        {"minConcept": {"rxcui": "83367", "tty": "IN"}, "rxclassMinConceptItem": {"classId": "C10AA"}},
        {"minConcept": {"rxcui": "83367", "tty": "IN"}, "rxclassMinConceptItem": {"classId": "N05AA"}},
    ]}}
    assert rxnav_atc_classes_for_ingredient("83367", rxclass_document) == ["C10AA", "N05AA"]  # the combination's code is dropped
    assert rxnav_atc_classes_for_ingredient("83367", {}) == []
    assert unichem_chembl_ids({"compounds": [], "response": "Not found"}) == []
    unichem_document = {"compounds": [{"sources": [{"shortName": "chembl", "compoundId": "CHEMBL1487"}, {"shortName": "fdasrs", "compoundId": "A0JWA85V8F"}, {"shortName": "chembl", "compoundId": "CHEMBL12"}]}], "response": "Success"}
    assert unichem_chembl_ids(unichem_document) == ["CHEMBL12", "CHEMBL1487"]


def test_pharmacological_bar_reuses_the_single_target_and_atc_n_rules() -> None:
    mechanisms_by_molecule = {
        "CHEMBL703": [{"molecule_chembl_id": "CHEMBL703", "target_chembl_id": "CHEMBL1", "action_type": "ANTAGONIST"}, {"molecule_chembl_id": "CHEMBL703", "target_chembl_id": "CHEMBL1", "action_type": "ANTAGONIST"}],
        "CHEMBL455": [{"molecule_chembl_id": "CHEMBL455", "target_chembl_id": "CHEMBL2", "action_type": "INHIBITOR"}, {"molecule_chembl_id": "CHEMBL455", "target_chembl_id": "CHEMBL3", "action_type": "INHIBITOR"}],
    }
    targets = {"CHEMBL1": {"organism": "Homo sapiens", "gene_symbols": ["CHRM1"]}, "CHEMBL2": {"organism": "Homo sapiens", "gene_symbols": ["DHPS"]}, "CHEMBL3": {"organism": "Escherichia coli", "gene_symbols": []}}
    single = mechanism_targets_for_parent("CHEMBL703", mechanisms_by_molecule, targets)
    assert [target.target_chembl_id for target in single] == ["CHEMBL1"] and passes_single_target_rule(single)
    human_only = mechanism_targets_for_parent("CHEMBL455", mechanisms_by_molecule, targets)
    assert [target.target_chembl_id for target in human_only] == ["CHEMBL2"]  # the E. coli target is dropped
    assert mechanism_targets_for_parent("CHEMBL455", mechanisms_by_molecule, targets, human_only=False).__len__() == 2
    assert mechanism_targets_for_parent(None, mechanisms_by_molecule, targets) == []
    assert not passes_single_target_rule([])
    bridge = IngredientBridge("10154", "succinylcholine", "RXCUI", BRIDGE_METHOD_UNII_UNICHEM, atc_codes_rxnav=["M03AB"], atc_codes_sider=["N01AX"])
    assert bridge.atc_codes == ["M03AB", "N01AX"] and passes_nervous_system_rule(bridge)
    assert not passes_nervous_system_rule(IngredientBridge("1", "x", "RXCUI", BRIDGE_METHOD_NONE, atc_codes_rxnav=["M03AB"]))


def test_bridge_json_roundtrip_and_counts(tmp_path: Path) -> None:
    unified = IngredientBridge("10154", "succinylcholine", "RXCUI", BRIDGE_METHOD_UNII_UNICHEM, unii_codes=["J2R869A8YF"], chembl_ids_seen=["CHEMBL703"], chembl_parent="CHEMBL703",
                               atc_codes_rxnav=["M03AB"], atc_codes_sider=["N01AX"], perturbation_id="CID100005314", unified_with_sider=True, sider_pubchem_cid=5314, sider_drug_name="succinylcholine",
                               sider_collision_pubchem_cids=[90000], mechanism_target_ids=["CHEMBL1"], single_mechanism_target=True, nervous_system_atc=True)
    onsides_only = IngredientBridge("OMOP997977", "Influenza vaccine", "OMOP", BRIDGE_METHOD_NONE, perturbation_id="OMOP:997977")
    bridges = {"10154": unified, "OMOP997977": onsides_only}
    path = tmp_path / "ingredient_identifier_bridge.json"
    path.write_text(json.dumps({"release": "fixture", "ingredients": {key: bridge_to_json_record(bridge) for key, bridge in bridges.items()}}), encoding="utf-8")
    loaded = load_ingredient_identifier_bridge(path)
    assert loaded == bridges
    assert perturbation_id_by_ingredient(loaded) == {"10154": "CID100005314", "OMOP997977": "OMOP:997977"}
    counts = bridge_counts(loaded)
    assert counts["ingredients"] == 2 and counts["method_unii_unichem"] == 1 and counts["method_none"] == 1 and counts["identifier_system_omop"] == 1
    assert counts["unified_with_sider"] == 1 and counts["onsides_only"] == 1 and counts["unified_with_collision"] == 1
    assert counts["single_mechanism_target"] == 1 and counts["with_atc_n_code"] == 1 and counts["single_target_and_atc_n"] == 1 and counts["single_target_and_atc_n_onsides_only"] == 0
