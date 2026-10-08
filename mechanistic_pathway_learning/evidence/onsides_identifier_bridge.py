"""Identifier bridge for the OnSIDES later label slice (design section 4.2, evidence class E2).

OnSIDES keys its label events by RxNorm ingredient (RxCUI, or an OMOP extension id for ingredients RxNorm
lacks) while SIDER keys its events by STITCH flat identifier (a PubChem CID) and the graph's drug
perturbations are built from ChEMBL mechanism targets on parent molecules. The two label slices only
report on the same perturbation when both resolve to the same ChEMBL parent, so this module holds the
decisions of that resolution as pure functions; experiments/fetch_onsides_identifier_bridge.py does the
network work and writes data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json, which
load_ingredient_identifier_bridge reads back.

Precedence per ingredient (the first route that yields a parent wins and is recorded as bridge_method):
  (a) unii_unichem   RxCUI -> UNII (RxNav property UNII_CODE) -> ChEMBL ids (UniChem, FDA SRS source)
                     -> parent molecule (ChEMBL molecule_hierarchy)
  (b) name_pref_name ingredient name -> ChEMBL pref_name, exact match ignoring case -> parent
  (c) none
An OMOP extension id has no RxNav record and goes straight to (b). When a route yields several distinct
parents the lowest-numbered one is taken and the rest are kept as alternatives, so the choice is
deterministic and visible.

Unification with SIDER (recorded as unification_method): an ingredient whose parent is among the ChEMBL
parents of a SIDER drug's STITCH flat identifier takes that identifier as its perturbation_id, the same
string SIDER reports carry in the report table (docs/evidence_reports_spec.md), so the evidence table sees
one perturbation with two label slices; otherwise an ingredient whose name equals a SIDER drug name
(case-insensitive) takes that drug's identifier, which covers SIDER drugs whose own ChEMBL mapping
(pubchem_to_chembl_with_parents.json) is empty or resolves to a stereo or racemate variant of the
ingredient's parent. The join is parent-to-flat-CID, not parent-to-parent: a STITCH flat identifier drops
stereochemistry, so one SIDER CID can carry several ChEMBL parents (enantiomers, epimers, deuterated
forms and in a few cases distinct drugs such as betamethasone and dexamethasone). The parents of the
chosen CID are recorded as sider_cid_chembl_parents and sider_cid_merges_distinct_parents flags a CID
with more than one, so a unified perturbation that pools stereoisomers stays visible. Several SIDER CIDs
can share a parent (salt forms listed as separate SIDER entries); the lowest PubChem CID is taken and the
others are recorded as a collision. An ingredient without a SIDER partner keeps "RXCUI:<id>" or
"OMOP:<id>" as its perturbation_id and is OnSIDES-only.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

from mechanistic_pathway_learning.evidence.load_drug_label_events import is_nervous_system_atc
from mechanistic_pathway_learning.perturbation.map_drug_targets_to_graph_nodes import (  # noqa: F401 (mechanisms_by_parent_and_molecule is re-exported)
    DrugTarget,
    has_dominant_target,
    mechanism_targets,
    mechanisms_by_parent_and_molecule,
)

BRIDGE_METHOD_UNII_UNICHEM = "unii_unichem"
BRIDGE_METHOD_NAME_PREF_NAME = "name_pref_name"
BRIDGE_METHOD_NONE = "none"
BRIDGE_METHODS: tuple[str, ...] = (BRIDGE_METHOD_UNII_UNICHEM, BRIDGE_METHOD_NAME_PREF_NAME, BRIDGE_METHOD_NONE)

IDENTIFIER_SYSTEM_RXCUI = "RXCUI"
IDENTIFIER_SYSTEM_OMOP = "OMOP"
OMOP_IDENTIFIER_PREFIX = "OMOP"

UNIFICATION_METHOD_CHEMBL_PARENT = "chembl_parent"
UNIFICATION_METHOD_SIDER_DRUG_NAME = "sider_drug_name"
UNIFICATION_METHOD_NONE = "none"
UNIFICATION_METHODS: tuple[str, ...] = (UNIFICATION_METHOD_CHEMBL_PARENT, UNIFICATION_METHOD_SIDER_DRUG_NAME, UNIFICATION_METHOD_NONE)

RXCLASS_INGREDIENT_TERM_TYPES = ("IN", "PIN")  # RxNorm ingredient and precise ingredient; MIN is a multi-ingredient concept


@dataclass
class IngredientBridge:
    """Everything the bridge recorded for one OnSIDES ingredient; one JSON object per ingredient id."""

    ingredient_id: str  # the vocab_rxnorm_ingredient.rxnorm_id string as shipped ("83367" or "OMOP997977")
    ingredient_name: str
    identifier_system: str  # RXCUI or OMOP
    bridge_method: str  # one of BRIDGE_METHODS
    unii_codes: list[str] = field(default_factory=list)
    chembl_ids_seen: list[str] = field(default_factory=list)  # every ChEMBL id any route returned
    chembl_parent: str | None = None
    chembl_parent_alternatives: list[str] = field(default_factory=list)  # other parents the winning route returned
    chembl_pref_name: str | None = None
    atc_codes_rxnav: list[str] = field(default_factory=list)  # ATC classes of the ingredient itself from RxNav
    atc_codes_sider: list[str] = field(default_factory=list)  # drug_atc.tsv codes of the unified SIDER drug
    perturbation_id: str = ""
    unified_with_sider: bool = False
    unification_method: str = UNIFICATION_METHOD_NONE  # one of UNIFICATION_METHODS
    sider_pubchem_cid: int | None = None
    sider_drug_name: str | None = None
    sider_collision_pubchem_cids: list[int] = field(default_factory=list)
    sider_cid_chembl_parents: list[str] = field(default_factory=list)  # every ChEMBL parent behind the chosen (stereo-flattened) SIDER CID
    sider_cid_merges_distinct_parents: bool = False  # the chosen CID carries more than one ChEMBL parent
    mechanism_target_ids: list[str] = field(default_factory=list)
    single_mechanism_target: bool = False
    nervous_system_atc: bool = False

    @property
    def atc_codes(self) -> list[str]:
        """Union of the RxNav and SIDER codes, sorted, for the ATC N test and for display."""
        return sorted(set(self.atc_codes_rxnav) | set(self.atc_codes_sider))


def identifier_system_for_ingredient(ingredient_id: str) -> str:
    """RXCUI for an all-digit id, OMOP for an OMOP extension id ("OMOP997977"); anything else is an error."""
    if ingredient_id.isdigit():
        return IDENTIFIER_SYSTEM_RXCUI
    if ingredient_id.startswith(OMOP_IDENTIFIER_PREFIX) and ingredient_id[len(OMOP_IDENTIFIER_PREFIX):].isdigit():
        return IDENTIFIER_SYSTEM_OMOP
    raise ValueError(f"unrecognized OnSIDES ingredient id: {ingredient_id!r}")


def fallback_perturbation_id(ingredient_id: str) -> str:
    """The perturbation id of an ingredient with no SIDER partner: "RXCUI:83367" or "OMOP:997977"."""
    if identifier_system_for_ingredient(ingredient_id) == IDENTIFIER_SYSTEM_OMOP:
        return f"{IDENTIFIER_SYSTEM_OMOP}:{ingredient_id[len(OMOP_IDENTIFIER_PREFIX):]}"
    return f"{IDENTIFIER_SYSTEM_RXCUI}:{ingredient_id}"


def stitch_flat_id_from_pubchem_cid(pubchem_cid: int) -> str:
    """The SIDER STITCH flat identifier of a PubChem CID: "CID1" and the CID zero-padded to eight digits."""
    return f"CID1{int(pubchem_cid):08d}"


def chembl_numeric_id(chembl_id: str) -> int:
    return int(chembl_id[len("CHEMBL"):])


def lowest_chembl_id(chembl_ids: set[str] | list[str]) -> str | None:
    """The lowest-numbered id, the deterministic choice when a route returns several parents; None when empty."""
    return min(chembl_ids, key=chembl_numeric_id) if chembl_ids else None


def parents_of_chembl_ids(chembl_ids: list[str], parent_by_chembl_id: dict[str, str]) -> set[str]:
    """Parent molecule of each id from the fetched molecule_hierarchy; an id without a record is its own parent."""
    return {parent_by_chembl_id.get(chembl_id) or chembl_id for chembl_id in chembl_ids}


def resolve_parent_from_unii_route(chembl_ids_from_unichem: list[str], parent_by_chembl_id: dict[str, str]) -> tuple[str | None, list[str]]:
    """Route (a): the parent of the UniChem hits, lowest id first, with the other parents as alternatives."""
    parents = parents_of_chembl_ids(chembl_ids_from_unichem, parent_by_chembl_id)
    winner = lowest_chembl_id(parents)
    return winner, sorted(parents - {winner}, key=chembl_numeric_id) if winner else []


def resolve_parent_from_name_route(molecules: list[dict]) -> tuple[str | None, list[str], str | None]:
    """Route (b): the parent of the ChEMBL molecules whose pref_name matched the ingredient name.

    molecules are ChEMBL molecule records with molecule_chembl_id, pref_name and molecule_hierarchy (the
    hierarchy is null for some records; such a molecule is then its own parent).
    """
    parents: set[str] = set()
    pref_name: str | None = None
    for molecule in molecules:
        hierarchy = molecule.get("molecule_hierarchy") or {}
        parents.add(hierarchy.get("parent_chembl_id") or molecule["molecule_chembl_id"])
        pref_name = pref_name or molecule.get("pref_name")
    winner = lowest_chembl_id(parents)
    return winner, sorted(parents - {winner}, key=chembl_numeric_id) if winner else [], pref_name


def apply_bridge_precedence(identifier_system: str, unii_route: tuple[str | None, list[str]], name_route: tuple[str | None, list[str], str | None]) -> tuple[str, str | None, list[str], str | None]:
    """(bridge_method, chembl_parent, alternatives, pref_name) under the precedence (a) then (b) then none.

    An OMOP id never takes route (a) even when a UNII route result is passed, since RxNav knows no OMOP ids;
    passing one is a caller mistake that would otherwise go unnoticed.
    """
    unii_parent, unii_alternatives = unii_route
    name_parent, name_alternatives, pref_name = name_route
    if identifier_system == IDENTIFIER_SYSTEM_RXCUI and unii_parent is not None:
        return BRIDGE_METHOD_UNII_UNICHEM, unii_parent, unii_alternatives, pref_name
    if name_parent is not None:
        return BRIDGE_METHOD_NAME_PREF_NAME, name_parent, name_alternatives, pref_name
    return BRIDGE_METHOD_NONE, None, [], pref_name


def sider_parent_molecules(pubchem_to_chembl_with_parents: dict[str, list[str]], molecule_parents: dict[str, dict], mechanisms: list[dict]) -> dict[str, set[str]]:
    """PubChem CID (string, as the cache keys it) -> the ChEMBL parent molecules of its UniChem hits.

    pubchem_to_chembl_with_parents.json lists each hit together with its parent; molecule_parents.json
    carries the hierarchy for hits without a mechanism and the mechanism table carries
    parent_molecule_chembl_id for the rest. A hit found in neither is taken as its own parent, which is
    what resolve_parent_molecules assumes for molecules that carry mechanisms.
    """
    parent_by_chembl_id: dict[str, str] = {chembl_id: record.get("parent") or chembl_id for chembl_id, record in molecule_parents.items()}
    for mechanism in mechanisms:
        molecule_id = mechanism.get("molecule_chembl_id")
        if molecule_id and molecule_id not in parent_by_chembl_id:
            parent_by_chembl_id[molecule_id] = mechanism.get("parent_molecule_chembl_id") or molecule_id
    return {pubchem_cid: parents_of_chembl_ids(chembl_ids, parent_by_chembl_id) for pubchem_cid, chembl_ids in pubchem_to_chembl_with_parents.items() if chembl_ids}


def sider_pubchem_cids_by_parent(sider_parents_by_cid: dict[str, set[str]]) -> dict[str, list[int]]:
    """ChEMBL parent -> SIDER PubChem CIDs sharing it, ascending, so collisions are visible before unification."""
    cids_by_parent: dict[str, set[int]] = {}
    for pubchem_cid, parents in sider_parents_by_cid.items():
        for parent in parents:
            cids_by_parent.setdefault(parent, set()).add(int(pubchem_cid))
    return {parent: sorted(cids) for parent, cids in cids_by_parent.items()}


def sider_drug_cids_by_lowercase_name(sider_drug_names: dict[str, str]) -> dict[str, list[int]]:
    """lowercase SIDER drug name -> PubChem CIDs carrying it (drug_names.tsv keyed by STITCH flat id), ascending."""
    cids_by_name: dict[str, set[int]] = {}
    for stitch_flat_id, drug_name in sider_drug_names.items():
        cids_by_name.setdefault(drug_name.strip().lower(), set()).add(int(stitch_flat_id[4:]))
    return {name: sorted(cids) for name, cids in cids_by_name.items()}


def unify_with_sider(ingredient_id: str, chembl_parent: str | None, cids_by_parent: dict[str, list[int]], ingredient_name: str = "", cids_by_lowercase_name: dict[str, list[int]] | None = None) -> tuple[str, int | None, list[int], str]:
    """(perturbation_id, SIDER PubChem CID or None, colliding CIDs, unification_method) for one ingredient.

    The parent route comes first: the CIDs whose ChEMBL parents include chembl_parent. When it yields nothing the
    name route takes the CIDs whose SIDER drug name equals the ingredient name ignoring case. Either way the
    perturbation_id is the STITCH flat identifier of the lowest CID, which is the perturbation_id SIDER reports
    carry, and the other CIDs of the route are returned as the collision list.
    """
    candidate_cids = cids_by_parent.get(chembl_parent, []) if chembl_parent else []
    unification_method = UNIFICATION_METHOD_CHEMBL_PARENT
    if not candidate_cids and cids_by_lowercase_name and ingredient_name.strip():
        candidate_cids = cids_by_lowercase_name.get(ingredient_name.strip().lower(), [])
        unification_method = UNIFICATION_METHOD_SIDER_DRUG_NAME
    if not candidate_cids:
        return fallback_perturbation_id(ingredient_id), None, [], UNIFICATION_METHOD_NONE
    chosen_cid = candidate_cids[0]
    return stitch_flat_id_from_pubchem_cid(chosen_cid), chosen_cid, list(candidate_cids[1:]), unification_method


def mechanism_targets_for_parent(chembl_parent: str | None, mechanisms_by_molecule: dict[str, list[dict]], targets: dict[str, dict], human_only: bool = True, missing_target_ids: set[str] | None = None) -> list[DrugTarget]:
    """Mechanism targets of a parent molecule, one entry per ChEMBL target, as drug_targets_for_pubchem_cid builds them for SIDER.

    A target id without a record in targets is kept (its organism is unknown, so the human filter cannot drop it
    and its gene symbols are empty) and added to missing_target_ids when a set is passed, so callers can report
    a missing record rather than a missing graph gene.
    """
    if chembl_parent is None:
        return []
    return mechanism_targets([chembl_parent], mechanisms_by_molecule, targets, human_only=human_only, missing_target_ids=missing_target_ids)


def passes_single_target_rule(drug_targets: list[DrugTarget]) -> bool:
    """The version 1 pharmacological bar on targets, the same rule SIDER drugs pass (has_dominant_target with one target)."""
    return has_dominant_target(drug_targets, max_targets=1)


def passes_nervous_system_rule(bridge: IngredientBridge) -> bool:
    """The ATC N proxy for central action on the union of the RxNav and SIDER codes (is_nervous_system_atc)."""
    return is_nervous_system_atc(bridge.atc_codes)


def rxnav_atc_classes_for_ingredient(rxcui: str, rxclass_document: dict) -> list[str]:
    """ATC classes RxNav attributes to the ingredient or its precise ingredients, sorted; multi-ingredient concepts are dropped.

    The byRxcui response lists every concept the ingredient belongs to: the ingredient itself (minConcept tty IN),
    its salts and esters (tty PIN, for example olmesartan medoxomil under olmesartan, which is where RxClass
    attaches the ATC class of many drugs) and combination products (tty MIN) with their own ATC codes. IN and
    PIN entries describe the ingredient; MIN entries do not. An entry without a term type is kept only when its
    minConcept.rxcui is the queried RxCUI.
    """
    entries = (rxclass_document.get("rxclassDrugInfoList") or {}).get("rxclassDrugInfo") or []
    classes: set[str] = set()
    for entry in entries:
        concept = entry.get("minConcept") or {}
        term_type = concept.get("tty")
        describes_ingredient = term_type in RXCLASS_INGREDIENT_TERM_TYPES if term_type else str(concept.get("rxcui")) == str(rxcui)
        class_id = (entry.get("rxclassMinConceptItem") or {}).get("classId")
        if describes_ingredient and class_id:
            classes.add(class_id)
    return sorted(classes)


def rxnav_unii_codes(property_document: dict) -> list[str]:
    """UNII codes in a RxNav property.json?propName=UNII_CODE response, sorted; [] for the {} of an unknown or UNII-less RxCUI."""
    concepts = (property_document.get("propConceptGroup") or {}).get("propConcept") or []
    return sorted({concept["propValue"] for concept in concepts if concept.get("propName") == "UNII_CODE" and concept.get("propValue")})


def unichem_chembl_ids(unichem_document: dict) -> list[str]:
    """ChEMBL ids among the sources of every compound in a UniChem v1 compounds response, sorted; [] for "Not found"."""
    return sorted({source["compoundId"] for compound in unichem_document.get("compounds", []) for source in compound.get("sources", []) if source.get("shortName") == "chembl" and source.get("compoundId")}, key=chembl_numeric_id)


def bridge_counts(bridges: dict[str, IngredientBridge]) -> dict[str, int]:
    """The counts docs/onsides_label_slice_spec.md reports for the bridge, computed from the bridge table alone."""
    values = list(bridges.values())
    counts = {"ingredients": len(values)}
    for method in BRIDGE_METHODS:
        counts[f"method_{method}"] = sum(1 for bridge in values if bridge.bridge_method == method)
    counts["identifier_system_omop"] = sum(1 for bridge in values if bridge.identifier_system == IDENTIFIER_SYSTEM_OMOP)
    counts["with_unii"] = sum(1 for bridge in values if bridge.unii_codes)
    counts["with_chembl_parent"] = sum(1 for bridge in values if bridge.chembl_parent)
    counts["distinct_chembl_parents"] = len({bridge.chembl_parent for bridge in values if bridge.chembl_parent})
    counts["with_parent_alternatives"] = sum(1 for bridge in values if bridge.chembl_parent_alternatives)
    counts["unified_with_sider"] = sum(1 for bridge in values if bridge.unified_with_sider)
    counts["distinct_sider_drugs_unified"] = len({bridge.sider_pubchem_cid for bridge in values if bridge.unified_with_sider})
    for method in UNIFICATION_METHODS:
        counts[f"unification_{method}"] = sum(1 for bridge in values if bridge.unification_method == method)
    counts["unified_with_collision"] = sum(1 for bridge in values if bridge.sider_collision_pubchem_cids)
    counts["unified_cid_merges_distinct_parents"] = sum(1 for bridge in values if bridge.sider_cid_merges_distinct_parents)
    counts["onsides_only"] = sum(1 for bridge in values if not bridge.unified_with_sider)
    counts["onsides_only_with_chembl_parent"] = sum(1 for bridge in values if not bridge.unified_with_sider and bridge.chembl_parent)
    counts["with_any_mechanism_target"] = sum(1 for bridge in values if bridge.mechanism_target_ids)
    counts["single_mechanism_target"] = sum(1 for bridge in values if bridge.single_mechanism_target)
    counts["with_atc_code"] = sum(1 for bridge in values if bridge.atc_codes)
    counts["with_atc_n_code"] = sum(1 for bridge in values if bridge.nervous_system_atc)
    counts["single_target_and_atc_n"] = sum(1 for bridge in values if bridge.single_mechanism_target and bridge.nervous_system_atc)
    counts["single_target_and_atc_n_onsides_only"] = sum(1 for bridge in values if bridge.single_mechanism_target and bridge.nervous_system_atc and not bridge.unified_with_sider)
    return counts


def bridge_to_json_record(bridge: IngredientBridge) -> dict:
    record = asdict(bridge)
    record["atc_codes"] = bridge.atc_codes
    return record


def bridge_from_json_record(record: dict) -> IngredientBridge:
    known_fields = {field_name for field_name in IngredientBridge.__dataclass_fields__}
    return IngredientBridge(**{key: value for key, value in record.items() if key in known_fields})


def load_ingredient_identifier_bridge(path: Path) -> dict[str, IngredientBridge]:
    """ingredient id -> IngredientBridge from ingredient_identifier_bridge.json (the "ingredients" object of the file)."""
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    return {ingredient_id: bridge_from_json_record(record) for ingredient_id, record in document["ingredients"].items()}


def perturbation_id_by_ingredient(bridges: dict[str, IngredientBridge]) -> dict[str, str]:
    """ingredient id -> perturbation_id, the join key the OnSIDES event loader needs."""
    return {ingredient_id: bridge.perturbation_id for ingredient_id, bridge in bridges.items()}
