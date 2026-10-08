"""Drug perturbations from ChEMBL mechanisms (design section 4.2, E2).

A drug becomes a perturbation of the graph through its mechanism targets: an
inhibitor, antagonist or blocker is a negative sign on the target nodes, an
agonist or activator a positive sign. Targets that are enzymes or transporters
reach the metabolic layer directly; receptor targets reach it only through the
signaling and transcription layers.

Every ChEMBL target type is placed where the graph has a node for it
(GraphNodeLookup.nodes_of_target, docs/drug_targets_any_type.md):
- protein targets (single proteins, complexes, families): the gene nodes of their gene symbols;
- nucleic-acid targets that name a gene (the mRNA or pre-mRNA of TTR, SMN2, APOB): that gene's node, by Ensembl id;
- metals, small molecules, lipids and oligosaccharides: the Human-GEM metabolites of
  configs/non_protein_drug_targets.csv, every compartment of each;
- a drug with no mechanism target under any of its ChEMBL forms whose active compound is itself a graph metabolite
  (GABA, tryptophan; configs/drugs_acting_as_graph_compounds.csv): that metabolite, raised (EXOGENOUS SUPPLY). Lithium is
  not listed: ChEMBL records IMPA1 and GSK3 inhibition on lithium carbonate and lithium citrate, so it falls under the
  single-target rule like any drug with protein targets (check_graph_compounds_have_no_mechanism enforces this).
Targets in another organism stay out (human_only), and DNA, RNA, antibodies and unnamed classes have no node.

Inputs are the caches written by experiments/fetch_chembl_drug_targets.py.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

NON_PROTEIN_TARGETS_PATH = Path("configs/non_protein_drug_targets.csv")
DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH = Path("configs/drugs_acting_as_graph_compounds.csv")
EXOGENOUS_SUPPLY = "EXOGENOUS SUPPLY"  # action of a drug that is itself the graph compound it raises
GRAPH_COMPOUND_TARGET_PREFIX = "COMPOUND:"  # target id of that action: COMPOUND:<Human-GEM base metabolite ids>

ACTION_TYPE_SIGN: dict[str, float] = {
    "INHIBITOR": -1.0,
    "ANTAGONIST": -1.0,
    "BLOCKER": -1.0,
    "NEGATIVE ALLOSTERIC MODULATOR": -0.5,
    "INVERSE AGONIST": -1.0,
    "CHELATING AGENT": -1.0,
    "SEQUESTERING AGENT": -1.0,
    "DEGRADER": -1.0,
    "RNAI INHIBITOR": -1.0,  # the target is a gene's mRNA, so its product falls
    "ANTISENSE INHIBITOR": -1.0,
    "HYDROLYTIC ENZYME": -1.0,  # an enzyme given as a drug (enzyme replacement) removes the substrate it targets
    "OXIDATIVE ENZYME": -1.0,
    "REDUCING AGENT": -1.0,  # the target molecule is reduced away (cysteamine on cystine)
    "AGONIST": 1.0,
    "PARTIAL AGONIST": 0.5,
    "ACTIVATOR": 1.0,
    "POSITIVE ALLOSTERIC MODULATOR": 0.5,
    "POSITIVE MODULATOR": 0.5,  # SMN2 pre-mRNA splicing modifiers raise full-length SMN
    "OPENER": 1.0,
    "RELEASING AGENT": 1.0,
    "STABILISER": 0.5,
    EXOGENOUS_SUPPLY: 1.0,
    "SUBSTRATE": 0.0,
    "MODULATOR": 0.0,
    "BINDING AGENT": 0.0,
    "CROSS-LINKING AGENT": 0.0,
    "OTHER": 0.0,
}


@dataclass
class DrugTarget:
    target_chembl_id: str
    action_type: str | None
    gene_symbols: list[str] = field(default_factory=list)
    target_type: str | None = None
    pref_name: str | None = None
    accessions: list[str] = field(default_factory=list)  # UniProt accessions, or the Ensembl gene id of a nucleic-acid target

    @property
    def sign(self) -> float:
        return ACTION_TYPE_SIGN.get((self.action_type or "OTHER").upper(), 0.0)


@dataclass
class ChemblCaches:
    pubchem_to_chembl: dict[str, list[str]]
    mechanisms_by_molecule: dict[str, list[dict]]
    targets: dict[str, dict]


def mechanisms_by_parent_and_molecule(mechanisms: list[dict]) -> dict[str, list[dict]]:
    """ChEMBL molecule id -> its mechanisms, keyed by molecule_chembl_id and by parent_molecule_chembl_id.

    ChEMBL records many mechanisms on a salt or prodrug form with the parent in parent_molecule_chembl_id (1,626 of the
    7,561 rows of ChEMBL 37; amitriptyline's are on its hydrochloride CHEMBL1200964, parent CHEMBL629), so a table keyed
    by molecule_chembl_id alone misses them for the parent. The SIDER route, the OnSIDES bridge, the target fetch and
    the report builder all key the table through this function so they apply one rule.
    """
    by_molecule: dict[str, list[dict]] = {}
    for mechanism in mechanisms:
        for key in sorted({mechanism.get("molecule_chembl_id"), mechanism.get("parent_molecule_chembl_id")} - {None, ""}):
            by_molecule.setdefault(key, []).append(mechanism)
    return by_molecule


def load_chembl_caches(chembl_directory: Path) -> ChemblCaches:
    pubchem_to_chembl = json.loads((chembl_directory / "pubchem_to_chembl_with_parents.json").read_text())
    mechanisms_by_molecule = mechanisms_by_parent_and_molecule(json.loads((chembl_directory / "mechanisms.json").read_text()))
    targets = json.loads((chembl_directory / "targets.json").read_text())
    return ChemblCaches(pubchem_to_chembl, mechanisms_by_molecule, targets)


def mechanism_targets(chembl_ids: list[str], mechanisms_by_molecule: dict[str, list[dict]], targets: dict[str, dict], human_only: bool = True,
                      missing_target_ids: set[str] | None = None) -> list[DrugTarget]:
    """Mechanism targets of the given ChEMBL molecules, one entry per ChEMBL target (families count as one).

    A target id without a record in targets is kept (its organism is unknown, so the human filter cannot drop it, and
    it has no gene symbols) and added to missing_target_ids when a set is passed.
    """
    targets_by_id: dict[str, DrugTarget] = {}
    for chembl_id in chembl_ids:
        for mechanism in mechanisms_by_molecule.get(chembl_id, []):
            target_id = mechanism.get("target_chembl_id")
            if not target_id:
                continue
            record = targets.get(target_id)
            if record is None:
                if missing_target_ids is not None:
                    missing_target_ids.add(target_id)
                record = {}
            if human_only and record.get("organism") not in (None, "Homo sapiens"):
                continue
            targets_by_id[target_id] = DrugTarget(target_id, mechanism.get("action_type"), record.get("gene_symbols", []), record.get("target_type"),
                                                  record.get("pref_name"), record.get("accessions", []))
    return list(targets_by_id.values())


def drug_targets_for_pubchem_cid(pubchem_cid: int, caches: ChemblCaches, human_only: bool = True, extra_chembl_ids: list[str] | None = None) -> list[DrugTarget]:
    """Mechanism targets for a SIDER drug. extra_chembl_ids adds molecules found another way (the ChEMBL parent the
    OnSIDES bridge unified with the drug when UniChem returned none for its PubChem id)."""
    chembl_ids = list(caches.pubchem_to_chembl.get(str(pubchem_cid), [])) + [chembl_id for chembl_id in extra_chembl_ids or [] if chembl_id]
    return mechanism_targets(chembl_ids, caches.mechanisms_by_molecule, caches.targets, human_only=human_only)


def has_dominant_target(drug_targets: list[DrugTarget], max_targets: int = 1) -> bool:
    """Version 1 dominant-target rule: at most max_targets distinct ChEMBL targets with a mechanism.

    Affinity margins are not yet used; a ChEMBL protein-family target (for example the GABA-A
    receptor) counts as one target even though it maps to many gene nodes.
    """
    return 0 < len(drug_targets) <= max_targets


def load_drugs_acting_as_graph_compounds(path: Path = DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH) -> tuple[dict[int, list[str]], dict[str, list[str]]]:
    """(PubChem CID -> base metabolite ids, ChEMBL id -> base metabolite ids) of configs/drugs_acting_as_graph_compounds.csv."""
    if not Path(path).exists():
        return {}, {}
    table = pd.read_csv(path, dtype=str).fillna("")
    by_pubchem_cid: dict[int, list[str]] = {}
    by_chembl_id: dict[str, list[str]] = {}
    for row in table.itertuples(index=False):
        base_ids = [base_id for base_id in row.base_metabolite_ids.split(";") if base_id]
        if row.pubchem_cid:
            by_pubchem_cid[int(row.pubchem_cid)] = base_ids
        for chembl_id in (identifier for identifier in row.chembl_ids.split(";") if identifier):
            by_chembl_id[chembl_id] = base_ids
    return by_pubchem_cid, by_chembl_id


def check_graph_compounds_have_no_mechanism(compounds_by_chembl_id: dict[str, list[str]], mechanisms_by_molecule: dict[str, list[dict]]) -> None:
    """Refuse a graph-compound entry whose ChEMBL form carries a mechanism: such a drug has protein (or other) targets and
    follows the mechanism rule, so one drug cannot enter through one identifier and be excluded through another."""
    with_mechanisms = sorted(chembl_id for chembl_id in compounds_by_chembl_id if mechanisms_by_molecule.get(chembl_id))
    if with_mechanisms:
        raise ValueError(f"{DRUGS_ACTING_AS_GRAPH_COMPOUNDS_PATH} lists ChEMBL forms that carry a mechanism: {with_mechanisms}")


def graph_compound_target(base_metabolite_ids: list[str]) -> DrugTarget:
    """The one target of a drug that is itself a graph compound: its metabolite nodes, raised."""
    return DrugTarget(GRAPH_COMPOUND_TARGET_PREFIX + ";".join(base_metabolite_ids), EXOGENOUS_SUPPLY, target_type="GRAPH COMPOUND")


def graph_compound_targets(drug_targets: list[DrugTarget], pubchem_cid: int | None, chembl_ids: list[str],
                           compounds_by_pubchem_cid: dict[int, list[str]], compounds_by_chembl_id: dict[str, list[str]]) -> list[DrugTarget]:
    """drug_targets unchanged when ChEMBL gives the drug any mechanism target; otherwise the graph-compound target of
    configs/drugs_acting_as_graph_compounds.csv when the drug is listed there, else none."""
    if drug_targets:
        return drug_targets
    base_ids = compounds_by_pubchem_cid.get(pubchem_cid) if pubchem_cid is not None else None
    for chembl_id in chembl_ids:
        base_ids = base_ids or compounds_by_chembl_id.get(chembl_id)
    return [graph_compound_target(base_ids)] if base_ids else []


@dataclass
class GraphNodeLookup:
    """Graph nodes a drug target stands for, for every ChEMBL target type (module docstring)."""

    node_by_symbol: dict[str, str]
    node_by_ensembl_gene: dict[str, str]
    metabolite_nodes_by_base_id: dict[str, list[str]]
    metabolites_by_target_id: dict[str, list[str]] = field(default_factory=dict)  # configs/non_protein_drug_targets.csv

    @classmethod
    def from_nodes(cls, nodes: pd.DataFrame, non_protein_targets_path: Path | None = NON_PROTEIN_TARGETS_PATH) -> "GraphNodeLookup":
        genes = nodes[nodes.node_type == "gene"]
        node_by_symbol = {symbol: node_id for symbol, node_id in zip(genes.gene_symbol, genes.node_id) if isinstance(symbol, str) and symbol}
        node_by_ensembl_gene = {}
        if "ensembl_gene_id" in genes.columns:
            node_by_ensembl_gene = {ensembl: node_id for ensembl, node_id in zip(genes.ensembl_gene_id, genes.node_id) if isinstance(ensembl, str) and ensembl}
        metabolites = nodes[nodes.node_type == "metabolite"]
        metabolite_nodes_by_base_id: dict[str, list[str]] = {}
        if "base_metabolite_id" in metabolites.columns:
            for base_id, node_id in zip(metabolites.base_metabolite_id, metabolites.node_id):
                if isinstance(base_id, str):
                    metabolite_nodes_by_base_id.setdefault(base_id, []).append(node_id)
        metabolites_by_target_id: dict[str, list[str]] = {}
        if non_protein_targets_path is not None and Path(non_protein_targets_path).exists():
            table = pd.read_csv(non_protein_targets_path, dtype=str).fillna("")
            metabolites_by_target_id = {row.target_chembl_id: [base_id for base_id in row.base_metabolite_ids.split(";") if base_id]
                                        for row in table.itertuples(index=False)}
        return cls(node_by_symbol, node_by_ensembl_gene, {base_id: sorted(node_ids) for base_id, node_ids in metabolite_nodes_by_base_id.items()}, metabolites_by_target_id)

    def nodes_of_target(self, target: DrugTarget) -> list[str]:
        gene_nodes = [self.node_by_symbol[symbol] for symbol in target.gene_symbols if symbol in self.node_by_symbol]
        if gene_nodes:
            return gene_nodes
        if target.target_type == "NUCLEIC-ACID":
            return [self.node_by_ensembl_gene[accession] for accession in target.accessions if accession in self.node_by_ensembl_gene]
        if target.target_chembl_id.startswith(GRAPH_COMPOUND_TARGET_PREFIX):
            base_ids = target.target_chembl_id.removeprefix(GRAPH_COMPOUND_TARGET_PREFIX).split(";")
        else:
            base_ids = self.metabolites_by_target_id.get(target.target_chembl_id, [])
        return [node_id for base_id in base_ids for node_id in self.metabolite_nodes_by_base_id.get(base_id, [])]

    def perturbation_nodes(self, drug_targets: list[DrugTarget]) -> list[list]:
        """[node id, sign, magnitude] triples; magnitude is 1 / number of graph nodes under each target."""
        triples: list[list] = []
        for target in drug_targets:
            mapped = self.nodes_of_target(target)
            triples.extend([node_id, target.sign, 1.0 / len(mapped)] for node_id in mapped)
        return triples
