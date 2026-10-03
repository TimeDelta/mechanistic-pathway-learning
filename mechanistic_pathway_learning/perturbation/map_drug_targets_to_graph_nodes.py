"""Drug perturbations from ChEMBL mechanisms (design section 4.2, E2).

A drug becomes a perturbation of the graph through its mechanism targets: an
inhibitor, antagonist or blocker is a negative sign on the target protein nodes,
an agonist or activator a positive sign. Targets that are enzymes or transporters
reach the metabolic layer directly; receptor targets reach it only through the
signaling and transcription layers.

Inputs are the caches written by experiments/fetch_chembl_drug_targets.py.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

ACTION_TYPE_SIGN: dict[str, float] = {
    "INHIBITOR": -1.0,
    "ANTAGONIST": -1.0,
    "BLOCKER": -1.0,
    "NEGATIVE ALLOSTERIC MODULATOR": -0.5,
    "INVERSE AGONIST": -1.0,
    "CHELATING AGENT": -1.0,
    "SEQUESTERING AGENT": -1.0,
    "DEGRADER": -1.0,
    "AGONIST": 1.0,
    "PARTIAL AGONIST": 0.5,
    "ACTIVATOR": 1.0,
    "POSITIVE ALLOSTERIC MODULATOR": 0.5,
    "OPENER": 1.0,
    "RELEASING AGENT": 1.0,
    "STABILISER": 0.5,
    "SUBSTRATE": 0.0,
    "MODULATOR": 0.0,
    "OTHER": 0.0,
}


@dataclass
class DrugTarget:
    target_chembl_id: str
    action_type: str | None
    gene_symbols: list[str] = field(default_factory=list)
    target_type: str | None = None
    pref_name: str | None = None

    @property
    def sign(self) -> float:
        return ACTION_TYPE_SIGN.get((self.action_type or "OTHER").upper(), 0.0)


@dataclass
class ChemblCaches:
    pubchem_to_chembl: dict[str, list[str]]
    mechanisms_by_molecule: dict[str, list[dict]]
    targets: dict[str, dict]


def load_chembl_caches(chembl_directory: Path) -> ChemblCaches:
    pubchem_to_chembl = json.loads((chembl_directory / "pubchem_to_chembl_with_parents.json").read_text())
    mechanisms_by_molecule: dict[str, list[dict]] = {}
    for mechanism in json.loads((chembl_directory / "mechanisms.json").read_text()):
        mechanisms_by_molecule.setdefault(mechanism["molecule_chembl_id"], []).append(mechanism)
    targets = json.loads((chembl_directory / "targets.json").read_text())
    return ChemblCaches(pubchem_to_chembl, mechanisms_by_molecule, targets)


def drug_targets_for_pubchem_cid(pubchem_cid: int, caches: ChemblCaches, human_only: bool = True) -> list[DrugTarget]:
    """Mechanism targets for a SIDER drug, one entry per ChEMBL target (families count as one)."""
    targets_by_id: dict[str, DrugTarget] = {}
    for chembl_id in caches.pubchem_to_chembl.get(str(pubchem_cid), []):
        for mechanism in caches.mechanisms_by_molecule.get(chembl_id, []):
            target_id = mechanism.get("target_chembl_id")
            if not target_id:
                continue
            record = caches.targets.get(target_id, {})
            if human_only and record.get("organism") not in (None, "Homo sapiens"):
                continue
            targets_by_id[target_id] = DrugTarget(target_id, mechanism.get("action_type"), record.get("gene_symbols", []), record.get("target_type"), record.get("pref_name"))
    return list(targets_by_id.values())


def has_dominant_target(drug_targets: list[DrugTarget], max_targets: int = 1) -> bool:
    """Version 1 dominant-target rule: at most max_targets distinct ChEMBL targets with a mechanism.

    Affinity margins are not yet used; a ChEMBL protein-family target (for example the GABA-A
    receptor) counts as one target even though it maps to many gene nodes.
    """
    return 0 < len(drug_targets) <= max_targets


def perturbation_from_drug_targets(drug_targets: list[DrugTarget], node_index_by_gene_symbol: dict[str, int]) -> list[tuple[int, float, float]]:
    """(node index, sign, magnitude) triples; magnitude is 1 / number of gene nodes under each target."""
    triples: list[tuple[int, float, float]] = []
    for target in drug_targets:
        mapped = [node_index_by_gene_symbol[symbol] for symbol in target.gene_symbols if symbol in node_index_by_gene_symbol]
        for node_index in mapped:
            triples.append((node_index, target.sign, 1.0 / len(mapped)))
    return triples
