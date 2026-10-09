"""Descriptors for the Reactome protein entities (complexes and entity sets): the entity's own annotation vector,
projected onto the same directions the protein descriptors use.

Why not the mean of the members' descriptors (the user, 9 October 2026: "I like the separate annotation vector much
better for the complexes"): averaging a kinase, a scaffold and a regulatory subunit gives a vector resembling none of
them, and over many members it tends to the population mean, so the result mostly encodes complex size.

What is built instead. The protein descriptors are the ESM-2 embeddings projected onto the k directions of annotation
space that a sequence predicts best (graph/protein_descriptors.py: Z = standardised X B V_k). A complex has no
sequence, but it has an annotation vector of its own, which the same V_k turns into the same coordinates. The vector is
built from the members' annotations by what each block means:

  - union (maximum) for Pfam families, EC numbers and GO molecular function: a complex carries every domain and every
    activity its subunits bring;
  - intersection (minimum) for GO cellular component and UniProt subcellular location: a complex exists only where its
    subunits sit together;
  - mean over the members for a DefinedSet or a CandidateSet, which is a disjunction ("one of these") rather than a
    machine, so no term is certain and the expected member is the honest reading. Reactome's 212-member entity is such
    a set, not a 212-protein complex.

A block whose annotation can be missing (GO function, GO component, location) is aggregated over the members annotated
in it, and the entity counts as unannotated in the block when none of its members is; EC numbers and Pfam families are
listed when they apply, so a missing one is a negative and every member counts.

The entity's vector is standardised with the statistics of the reviewed human proteome, not of the entities, so the
entities land in the proteins' units. The projection of a protein's *observed* annotations is not the same quantity as
its descriptor, which is the projection of the annotations *predicted from its sequence*; the build script reports the
correlation between the two over the proteome so the comparability is a measurement rather than an assumption.
"""
from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd

UNION_BLOCKS = ("ec_level_2", "go_function", "pfam")
INTERSECTION_BLOCKS = ("go_component", "location")
BLOCKS_WITH_MISSING_ANNOTATION = ("go_function", "go_component", "location")
SET_KINDS = frozenset({"DefinedSet", "CandidateSet", "OpenSet"})
MEMBER_OF_RELATION = "member_of"
PROTEIN_ENTITY_TYPE = "protein_entity"


def member_genes_of_entity(edges: pd.DataFrame, entity_node_ids) -> dict[str, list[str]]:
    """Protein entity node -> the gene symbols of its members, from the member_of edges the Reactome import writes."""
    wanted = set(entity_node_ids)
    members: dict[str, set[str]] = defaultdict(set)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        if relation == MEMBER_OF_RELATION and target in wanted and str(source).startswith("GENE:"):
            members[target].add(str(source)[5:])
    return {entity: sorted(genes) for entity, genes in members.items()}


def entity_annotation_rows(member_keys_of_entity: dict[str, list[str]], kind_of_entity: dict[str, str],
                           targets: pd.DataFrame, observed: pd.DataFrame, block_columns: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The entities' own 0/1 annotation rows and whether each entity is annotated in each block.

    One row per entity with at least one member whose annotations are known; a set kind takes the mean over its
    members, a complex the union or the intersection of its members by block (see the module docstring). The member
    keys index targets and observed: the caller passes a gene-level matrix (a gene's reviewed entries collapsed by
    union) so that a gene with several entries counts once."""
    entities = [entity for entity, accessions in member_keys_of_entity.items() if accessions]
    values = pd.DataFrame(0.0, index=pd.Index(entities, name="node_id"), columns=targets.columns)
    annotated = pd.DataFrame(True, index=values.index, columns=list(block_columns))
    target_values = targets.to_numpy(dtype=np.float64)
    row_of_member = {key: row for row, key in enumerate(targets.index)}
    for entity in entities:
        rows = [row_of_member[key] for key in member_keys_of_entity[entity] if key in row_of_member]
        is_a_set = kind_of_entity.get(entity) in SET_KINDS
        for block, columns in block_columns.items():
            if not columns or not rows:
                annotated.loc[entity, block] = bool(rows) and block not in BLOCKS_WITH_MISSING_ANNOTATION
                continue
            member_rows = rows
            if block in BLOCKS_WITH_MISSING_ANNOTATION:
                member_rows = [row for row in rows if observed[block].to_numpy()[row]]
                if not member_rows:
                    annotated.loc[entity, block] = False
                    continue
            positions = [targets.columns.get_loc(column) for column in columns]
            block_values = target_values[np.ix_(member_rows, positions)]
            if is_a_set:
                aggregated = block_values.mean(axis=0)
            elif block in INTERSECTION_BLOCKS:
                aggregated = block_values.min(axis=0)
            else:
                aggregated = block_values.max(axis=0)
            values.loc[entity, columns] = aggregated
    return values, annotated


def standardise_with_proteome_statistics(entity_values: pd.DataFrame, entity_annotated: pd.DataFrame, protein_targets: pd.DataFrame,
                                         protein_observed: pd.DataFrame, block_columns: dict[str, list[str]]) -> np.ndarray:
    """The entities' annotation rows standardised exactly as graph/protein_descriptors.weighted_standardised_targets
    standardises the proteins', but with the proteome's own block means and deviations, so the two live in one space.
    An entity unannotated in a block gets 0 there, the proteome's mean, as an unannotated protein does."""
    values = entity_values.to_numpy(dtype=np.float64).copy()
    protein_values = protein_targets.to_numpy(dtype=np.float64)
    for block, columns in block_columns.items():
        if not columns:
            continue
        positions = [entity_values.columns.get_loc(column) for column in columns]
        protein_rows = protein_observed[block].to_numpy()
        block_protein_values = protein_values[np.ix_(protein_rows, positions)]
        mean, deviation = block_protein_values.mean(axis=0), block_protein_values.std(axis=0)
        rows = entity_annotated[block].to_numpy()
        standardised = np.zeros((len(entity_values), len(positions)))
        standardised[rows] = (values[np.ix_(rows, positions)] - mean) / np.where(deviation > 0, deviation, 1.0)
        values[:, positions] = standardised / np.sqrt(len(positions))
    return values
