"""Shared loading of the graph and the evidence table into arrays for every experiment.

Produces the perturbation-by-symptom outcome matrix for one relation, the node
indices each perturbation seeds, leakage group ids for the grouped splits and the
edge arrays for the encoder. Unobserved pairs are outcome 0 (unlabelled, not
confirmed negative) in evaluation, which is the standard positive-unlabelled
convention and is stated as such in the design (section 5.4).
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class ExperimentData:
    node_ids: list[str]
    node_index: dict[str, int]
    node_types: np.ndarray
    is_currency: np.ndarray
    node_degree: np.ndarray
    edge_source: np.ndarray
    edge_target: np.ndarray
    edge_relation: np.ndarray
    relation_types: list[str]
    symptoms: list[str]
    perturbation_ids: list[str]
    perturbation_labels: list[str]
    perturbation_types: list[str]
    group_ids: list[str]
    perturbation_seeds: list[np.ndarray]  # node indices per perturbation
    perturbation_signs: list[np.ndarray]
    perturbation_magnitudes: list[np.ndarray]
    outcomes: np.ndarray  # [num_perturbations, num_symptoms]
    weights: np.ndarray  # [num_perturbations, num_symptoms]
    in_metabolic_layer: np.ndarray
    node_subsystem: np.ndarray | None = None  # reconstruction subsystem per node ("" for non-reactions or when absent)

    @property
    def perturbation_degrees(self) -> np.ndarray:
        return np.array([self.node_degree[seeds].sum() if len(seeds) else 0.0 for seeds in self.perturbation_seeds])


GROUPING_COLUMNS = {"gene": "group_id", "disease_cluster": "disease_cluster_id"}


def load_experiment_data(graph_directory: Path, evidence_directory: Path, relation: str = "induces", symptoms: list[str] | None = None, metabolic_layer_only: bool = False, group_by: str = "gene") -> ExperimentData:
    """group_by selects the leakage group for the grouped split: "gene" (the gene itself; drugs by dominant
    target) or "disease_cluster" (genes sharing a disease entry in HPO are held out together)."""
    if group_by not in GROUPING_COLUMNS:
        raise ValueError(f"group_by must be one of {sorted(GROUPING_COLUMNS)}")
    group_column = GROUPING_COLUMNS[group_by]
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    relation_types = json.loads((graph_directory / "relation_types.json").read_text())
    node_ids = list(nodes.node_id)
    node_index = {node_id: index for index, node_id in enumerate(node_ids)}
    relation_index = {name: index for index, name in enumerate(relation_types)}
    evidence = pd.read_parquet(evidence_directory / "evidence_records.parquet")
    evidence = evidence[evidence.relation == relation]
    if group_column not in evidence.columns:  # evidence tables written before disease clusters existed
        evidence = evidence.assign(**{group_column: evidence.group_id})
    if metabolic_layer_only:
        evidence = evidence[evidence.in_metabolic_layer == True]  # noqa: E712
    if symptoms is None:
        symptoms = sorted(evidence.symptom.unique())
    symptom_index = {symptom: index for index, symptom in enumerate(symptoms)}
    evidence = evidence[evidence.symptom.isin(symptom_index)]
    perturbation_ids = sorted(evidence.perturbation_id.unique())
    perturbation_position = {perturbation_id: index for index, perturbation_id in enumerate(perturbation_ids)}
    outcomes = np.zeros((len(perturbation_ids), len(symptoms)))
    weights = np.zeros_like(outcomes)
    labels, types, groups, seeds, signs, magnitudes, metabolic = {}, {}, {}, {}, {}, {}, {}
    for row in evidence.itertuples(index=False):
        position = perturbation_position[row.perturbation_id]
        outcomes[position, symptom_index[row.symptom]] = 1.0
        weights[position, symptom_index[row.symptom]] = max(weights[position, symptom_index[row.symptom]], float(row.weight))
        labels[row.perturbation_id] = row.perturbation_label
        types[row.perturbation_id] = row.perturbation_type
        groups[row.perturbation_id] = getattr(row, group_column)
        metabolic[row.perturbation_id] = bool(row.in_metabolic_layer)
        if row.perturbation_id not in seeds:
            triples = json.loads(row.perturbation_nodes)
            seeds[row.perturbation_id] = np.array([node_index[node_id] for node_id, _, _ in triples if node_id in node_index], dtype=int)
            signs[row.perturbation_id] = np.array([sign for node_id, sign, _ in triples if node_id in node_index])
            magnitudes[row.perturbation_id] = np.array([magnitude for node_id, _, magnitude in triples if node_id in node_index])
    return ExperimentData(
        node_ids=node_ids,
        node_index=node_index,
        node_types=nodes.node_type.to_numpy(),
        is_currency=nodes.is_currency.fillna(False).to_numpy().astype(bool),
        node_degree=nodes.degree.to_numpy().astype(float),
        edge_source=edges.source_id.map(node_index).to_numpy(),
        edge_target=edges.target_id.map(node_index).to_numpy(),
        edge_relation=edges.relation_type.map(relation_index).to_numpy(),
        relation_types=relation_types,
        symptoms=symptoms,
        perturbation_ids=perturbation_ids,
        perturbation_labels=[labels[p] for p in perturbation_ids],
        perturbation_types=[types[p] for p in perturbation_ids],
        group_ids=[groups[p] for p in perturbation_ids],
        perturbation_seeds=[seeds[p] for p in perturbation_ids],
        perturbation_signs=[signs[p] for p in perturbation_ids],
        perturbation_magnitudes=[magnitudes[p] for p in perturbation_ids],
        outcomes=outcomes,
        weights=weights,
        in_metabolic_layer=np.array([metabolic[p] for p in perturbation_ids]),
        node_subsystem=nodes.subsystem.fillna("").to_numpy().astype(str) if "subsystem" in nodes.columns else None,
    )
