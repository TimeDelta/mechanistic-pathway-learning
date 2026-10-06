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
from datetime import date
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
    frequencies: np.ndarray | None = None  # [num_perturbations, num_symptoms] reported frequency of a positive pair, NaN when unknown or negative
    evidence_dates: np.ndarray | None = None  # [num_perturbations, num_symptoms] proleptic Gregorian ordinal of the earliest dated evidence behind a positive pair, 0 when undated or negative
    evidence_date_is_publication: np.ndarray | None = None  # [num_perturbations, num_symptoms] True when the earliest date is a cited publication date rather than an HPO curation date (review v0.4, finding 7)
    node_compartment: np.ndarray | None = None  # compartment string per node ("" when none; "c;m" for a transport reaction)
    node_is_transport: np.ndarray | None = None
    node_is_reversible: np.ndarray | None = None
    node_brain_expression: np.ndarray | None = None  # log1p of the largest GTEx brain median TPM per node (genes, and reactions through their genes); 0 when unknown
    node_brain_expressed: np.ndarray | None = None  # median TPM at least the build threshold in one brain tissue; False when unknown

    def structural_node_features(self) -> np.ndarray:
        """Fixed per-node features with no node identity: one-hot type, multi-hot compartment, log degree, currency, transport and reversibility flags, brain expression (log TPM and expressed flag).

        These are what the inductive encoder variant reads instead of a learned embedding per node, so a
        model built on them can only use graph structure and node kinds (design section 5.2, version 0.4 ablation).
        """
        node_types = sorted(set(self.node_types.tolist()))
        compartments = sorted({part for value in (self.node_compartment if self.node_compartment is not None else []) for part in str(value).split(";") if part})
        num_nodes = len(self.node_ids)
        features = np.zeros((num_nodes, len(node_types) + len(compartments) + 6), dtype=np.float32)
        for index, node_type in enumerate(self.node_types):
            features[index, node_types.index(node_type)] = 1.0
        if self.node_compartment is not None:
            for index, value in enumerate(self.node_compartment):
                for part in str(value).split(";"):
                    if part:
                        features[index, len(node_types) + compartments.index(part)] = 1.0
        offset = len(node_types) + len(compartments)
        features[:, offset] = np.log1p(self.node_degree)
        features[:, offset + 1] = self.is_currency.astype(np.float32)
        if self.node_is_transport is not None:
            features[:, offset + 2] = self.node_is_transport.astype(np.float32)
        if self.node_is_reversible is not None:
            features[:, offset + 3] = self.node_is_reversible.astype(np.float32)
        if self.node_brain_expression is not None:
            features[:, offset + 4] = self.node_brain_expression.astype(np.float32)
        if self.node_brain_expressed is not None:
            features[:, offset + 5] = self.node_brain_expressed.astype(np.float32)
        return features

    @property
    def perturbation_degrees(self) -> np.ndarray:
        return np.array([self.node_degree[seeds].sum() if len(seeds) else 0.0 for seeds in self.perturbation_seeds])


GROUPING_COLUMNS = {"gene": "group_id", "disease_cluster": "disease_cluster_id"}


DEFAULT_LABEL_GRADES: tuple[str, ...] = ("A", "B")  # grades that count as positive labels; grade C (human association) and lower are soft evidence, not labels
SOFT_PRIOR_ONLY_GRADES: tuple[str, ...] = ("D", "E")  # literature grades: soft priors by design (section 4.2), never labels and never evaluation positives


def load_experiment_data(graph_directory: Path, evidence_directory: Path, relation: str = "induces", symptoms: list[str] | None = None, metabolic_layer_only: bool = False, group_by: str = "gene", label_grades: tuple[str, ...] | None = DEFAULT_LABEL_GRADES) -> ExperimentData:
    """group_by selects the leakage group for the grouped split: "gene" (the gene itself; drugs by dominant
    target) or "disease_cluster" (genes sharing a disease entry in HPO are held out together). label_grades restricts
    the rows that become positive labels (default A and B; None keeps every grade, which the version 0.3 ablation
    table of grade C weight-0 positives relies on)."""
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
    if "positive_report_count" in evidence.columns:  # tables written before reports existed have no such column and every row is a positive claim
        evidence = evidence[evidence.positive_report_count > 0]
    if label_grades is not None and set(label_grades) & set(SOFT_PRIOR_ONLY_GRADES):
        raise ValueError(f"label_grades {sorted(label_grades)} include a literature grade; grades {SOFT_PRIOR_ONLY_GRADES} are soft priors and never labels (design section 4.2)")
    if label_grades is not None and "grade" in evidence.columns:
        evidence = evidence[evidence.grade.isin(label_grades)]
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
    frequencies = np.full_like(outcomes, np.nan)
    evidence_dates = np.zeros(outcomes.shape, dtype=np.int64)
    has_frequency_column = "label_frequency" in evidence.columns
    has_date_column = "evidence_date" in evidence.columns
    has_date_source_column = "evidence_date_source" in evidence.columns
    evidence_date_is_publication = np.zeros(outcomes.shape, dtype=bool)
    labels, types, groups, seeds, signs, magnitudes, metabolic = {}, {}, {}, {}, {}, {}, {}
    for row in evidence.itertuples(index=False):
        position = perturbation_position[row.perturbation_id]
        outcomes[position, symptom_index[row.symptom]] = 1.0
        weights[position, symptom_index[row.symptom]] = max(weights[position, symptom_index[row.symptom]], float(row.weight))
        if has_frequency_column and row.label_frequency is not None and not (isinstance(row.label_frequency, float) and np.isnan(row.label_frequency)):
            frequencies[position, symptom_index[row.symptom]] = np.nanmax([frequencies[position, symptom_index[row.symptom]], float(row.label_frequency)])
        if has_date_column and isinstance(row.evidence_date, str) and row.evidence_date:
            ordinal = date.fromisoformat(row.evidence_date).toordinal()
            current = evidence_dates[position, symptom_index[row.symptom]]
            evidence_dates[position, symptom_index[row.symptom]] = ordinal if current == 0 else min(current, ordinal)
            if has_date_source_column and (current == 0 or ordinal <= current):
                evidence_date_is_publication[position, symptom_index[row.symptom]] = str(row.evidence_date_source) == "publication"
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
        frequencies=frequencies,
        evidence_dates=evidence_dates,
        evidence_date_is_publication=evidence_date_is_publication if has_date_source_column else None,
        node_compartment=nodes.compartment.fillna("").to_numpy().astype(str) if "compartment" in nodes.columns else None,
        node_is_transport=(nodes.is_transport == True).to_numpy() if "is_transport" in nodes.columns else None,  # noqa: E712 - NaN rows become False
        node_is_reversible=(nodes.reversible == True).to_numpy() if "reversible" in nodes.columns else None,  # noqa: E712
        node_brain_expression=np.log1p(pd.to_numeric(nodes.brain_median_tpm_max, errors="coerce").fillna(0.0).to_numpy(dtype=float)) if "brain_median_tpm_max" in nodes.columns else None,
        node_brain_expressed=(nodes.brain_expressed == True).to_numpy() if "brain_expressed" in nodes.columns else None,  # noqa: E712
    )
