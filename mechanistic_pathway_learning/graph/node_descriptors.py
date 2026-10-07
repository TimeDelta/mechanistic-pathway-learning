"""Fixed node descriptors from known chemistry and sequence, one block per node type (design section 5.2).

They enter the encoders as extra input columns: each node type fills its own block and leaves the others at zero, so the
encoder's existing learned input layer (feature_projection of the message passing encoder, output_gate of the
linear-response encoder) acts as one linear map per node type, shared by all nodes of that type. Nothing is learned per
node, and a held-out gene still has its descriptors.

Blocks:
  - metabolite: interpretable physicochemical properties of the structure, no statistical reduction: log molecular
    weight, Crippen logP, net charge in the model (SBML fbc:charge, the form Human-GEM uses near pH 7.3), log
    topological polar surface area, log counts of hydrogen-bond donors, acceptors, rotatable bonds and rings, a flag for
    having a complete structure and one for a partial structure (of 4,165 base metabolites, 778 have no SMILES, mostly
    generic pools, and 422 a SMILES with R-group atoms '*', mostly acyl chains of acyl-CoAs, whose properties would
    describe a fragment; both get the properties' mean, and net charge, known for all, stays);
  - reaction: the EC class (first digit, multi-hot over classes 1 to 7) of the reaction's annotated EC numbers, and a
    flag for having one;
  - gene: protein descriptors from protein_descriptors.py (ESM-2 embeddings reduced by reduced-rank regression onto
    function annotations), passed in as a table.
Continuous columns are standardised over the nodes of their type that have a value; missing values become 0 (the
mean) with the flag at 0.
"""
from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

METABOLITE_CONTINUOUS_COLUMNS = ["log_molecular_weight", "logp", "net_charge", "log_polar_surface_area", "log_hydrogen_bond_donors",
                                 "log_hydrogen_bond_acceptors", "log_rotatable_bonds", "log_rings"]
METABOLITE_FLAG_COLUMNS = ["has_structure", "partial_structure"]
NUMBER_OF_EC_CLASSES = 7
SPECIES_PATTERN = re.compile(r'<species [^>]*id="M_(MAM\d+)[a-z]+"[^>]*fbc:charge="(-?\d+)"')
REACTION_PATTERN = re.compile(r'<reaction [^>]*id="R_([^"]+)"')
EC_PATTERN = re.compile(r"ec-code/(\d)\.")


def metabolite_charges_from_sbml(sbml_text: str) -> dict[str, int]:
    """Net charge per Human-GEM base metabolite (the same in every compartment)."""
    return {base_id: int(charge) for base_id, charge in SPECIES_PATTERN.findall(sbml_text)}


def physicochemical_properties(smiles: str) -> dict[str, float] | None:
    """RDKit properties of one structure, or None when RDKit cannot parse it."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    return {
        "log_molecular_weight": math.log(max(Descriptors.MolWt(molecule), 1.0)),
        "logp": Crippen.MolLogP(molecule),
        "log_polar_surface_area": math.log1p(rdMolDescriptors.CalcTPSA(molecule)),
        "log_hydrogen_bond_donors": math.log1p(Lipinski.NumHDonors(molecule)),
        "log_hydrogen_bond_acceptors": math.log1p(Lipinski.NumHAcceptors(molecule)),
        "log_rotatable_bonds": math.log1p(rdMolDescriptors.CalcNumRotatableBonds(molecule)),
        "log_rings": math.log1p(rdMolDescriptors.CalcNumRings(molecule)),
    }


def metabolite_descriptors(metabolites_table: pd.DataFrame, sbml_text: str) -> pd.DataFrame:
    """One row per base metabolite (metsNoComp): METABOLITE_CONTINUOUS_COLUMNS unstandardised (NaN where unknown) and
    the flags of METABOLITE_FLAG_COLUMNS."""
    charges = metabolite_charges_from_sbml(sbml_text)
    rows = {}
    for base_id, smiles in metabolites_table.drop_duplicates("metsNoComp")[["metsNoComp", "metSmiles"]].itertuples(index=False):
        has_smiles = isinstance(smiles, str) and bool(smiles)
        is_partial = has_smiles and "*" in smiles
        properties = physicochemical_properties(smiles) if has_smiles and not is_partial else None
        rows[base_id] = {**(properties or {}), "net_charge": charges.get(base_id, np.nan), "has_structure": float(properties is not None),
                         "partial_structure": float(is_partial)}
    return pd.DataFrame.from_dict(rows, orient="index").reindex(columns=METABOLITE_CONTINUOUS_COLUMNS + METABOLITE_FLAG_COLUMNS)


def reaction_enzyme_classes(sbml_text: str) -> pd.DataFrame:
    """One row per reaction id: ec_class_1 ... ec_class_7 (multi-hot) and has_ec."""
    matches = list(REACTION_PATTERN.finditer(sbml_text))
    rows = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(sbml_text)
        block_end = sbml_text.find("</reaction>", match.start(), end)
        classes = {int(digit) for digit in EC_PATTERN.findall(sbml_text, match.start(), block_end if block_end > 0 else end)}
        row = {f"ec_class_{number}": float(number in classes) for number in range(1, NUMBER_OF_EC_CLASSES + 1)}
        row["has_ec"] = float(bool(classes))
        rows[match.group(1)] = row
    return pd.DataFrame.from_dict(rows, orient="index")


def standardise(values: pd.DataFrame) -> pd.DataFrame:
    """Z-score each column over its non-missing values; missing values become 0."""
    mean, standard_deviation = values.mean(), values.std(ddof=0).replace(0.0, 1.0)
    return ((values - mean) / standard_deviation).fillna(0.0)


def assemble_node_descriptor_table(nodes: pd.DataFrame, metabolite_table: pd.DataFrame | None = None, reaction_table: pd.DataFrame | None = None,
                                   protein_table: pd.DataFrame | None = None) -> pd.DataFrame:
    """node_id-indexed table with one block per node type (columns prefixed metabolite_, reaction_, protein_), zero
    outside the block's own type. protein_table is indexed by gene symbol."""
    blocks = []
    node_index = pd.Index(nodes.node_id, name="node_id")
    if metabolite_table is not None:
        is_metabolite = (nodes.node_type == "metabolite").to_numpy()
        per_node = metabolite_table.reindex(nodes.base_metabolite_id.where(is_metabolite)).set_axis(node_index)
        continuous = standardise(per_node.loc[is_metabolite, METABOLITE_CONTINUOUS_COLUMNS])
        block = pd.DataFrame(0.0, index=node_index, columns=METABOLITE_CONTINUOUS_COLUMNS + METABOLITE_FLAG_COLUMNS)
        block.loc[is_metabolite, METABOLITE_CONTINUOUS_COLUMNS] = continuous.to_numpy()
        block.loc[is_metabolite, METABOLITE_FLAG_COLUMNS] = per_node.loc[is_metabolite, METABOLITE_FLAG_COLUMNS].fillna(0.0).to_numpy()
        blocks.append(block.add_prefix("metabolite_"))
    if reaction_table is not None:
        is_reaction = (nodes.node_type == "reaction").to_numpy()
        block = pd.DataFrame(0.0, index=node_index, columns=reaction_table.columns)
        block.loc[is_reaction] = reaction_table.reindex(nodes.node_id[is_reaction]).fillna(0.0).to_numpy()
        blocks.append(block.add_prefix("reaction_"))
    if protein_table is not None:
        is_gene = (nodes.node_type == "gene").to_numpy()
        per_gene = protein_table.reindex(nodes.gene_symbol[is_gene])
        block = pd.DataFrame(0.0, index=node_index, columns=list(protein_table.columns) + ["has_protein_descriptors"])
        block.loc[is_gene, list(protein_table.columns)] = per_gene.fillna(0.0).to_numpy()
        block.loc[is_gene, "has_protein_descriptors"] = per_gene.notna().all(axis=1).astype(float).to_numpy()
        blocks.append(block.add_prefix("protein_"))
    return pd.concat(blocks, axis=1) if blocks else pd.DataFrame(index=node_index)
