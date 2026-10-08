"""Shared loading of the graph and the evidence table into arrays for every experiment.

Produces the perturbation-by-symptom outcome matrix for one relation, the node
indices each perturbation seeds, leakage group ids for the grouped splits and the
edge arrays for the encoder. Unobserved pairs are outcome 0 (unlabelled, not
confirmed negative) in evaluation, which is the standard positive-unlabelled
convention and is stated as such in the design (section 5.4).
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import warnings
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
    edge_sign: np.ndarray | None = None  # +1 activating, -1 inhibiting or repressing, 0 unsigned (binds); +1 where the graph has no sign column
    node_base_metabolite_id: np.ndarray | None = None  # metabolite id without compartment ("" for other node types)
    node_display_name: np.ndarray | None = None
    label_mask: np.ndarray | None = None  # [num_perturbations, num_symptoms] True where a pair is labelled; False for a positive pair a label selection set aside (neither positive nor negative); None when no selection was given
    label_selection_summary: dict | None = None
    node_descriptor_table: pd.DataFrame | None = None  # set by a rewired run whose reaction expression follows the rewiring; None: read --node-descriptors
    cell_class_weight_table: pd.DataFrame | None = None  # the same for --cell-class-weights

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


GROUPING_COLUMNS = {"gene": "group_id", "disease_cluster": "disease_cluster_id", "disease_cluster_and_targets": "disease_cluster_id"}


PER_PERTURBATION_LISTS = ("perturbation_ids", "perturbation_labels", "perturbation_types", "group_ids", "perturbation_seeds", "perturbation_signs", "perturbation_magnitudes")
PER_PERTURBATION_ARRAYS = ("outcomes", "weights", "in_metabolic_layer", "frequencies", "evidence_dates", "evidence_date_is_publication", "label_mask")


def restrict_to_perturbations(data: "ExperimentData", keep: np.ndarray) -> "ExperimentData":
    """A copy of data holding only the perturbations where keep is True (graph fields shared, not copied). Used to take the
    lockbox out of every pilot run, baseline and score, so nothing fitted or scored before confirmation sees it."""
    keep = np.asarray(keep, dtype=bool)
    if keep.shape != (len(data.perturbation_ids),):
        raise ValueError("keep must have one entry per perturbation")
    positions = np.flatnonzero(keep)
    changes = {name: [getattr(data, name)[i] for i in positions] for name in PER_PERTURBATION_LISTS}
    for name in PER_PERTURBATION_ARRAYS:
        value = getattr(data, name)
        changes[name] = None if value is None else value[positions]
    if data.label_selection_summary is not None and changes["label_mask"] is not None:
        # the summary's counts are over all perturbations; these are over the rows kept (development, without the lockbox)
        positive, kept = changes["outcomes"] > 0, np.asarray(changes["label_mask"], dtype=bool)
        changes["label_selection_summary"] = {**data.label_selection_summary, "after_restriction": {
            "perturbations": int(len(positions)), "positive_pairs": int(positive.sum()), "masked_pairs": int((positive & ~kept).sum()), "masked_negative_pairs": int((~positive & ~kept).sum()),
            "kept_positive_pairs": int((positive & kept).sum())}}
    return dataclasses.replace(data, **changes)


def read_lockbox(lockbox_path: Path, data: "ExperimentData", group_by: str | None = None, evidence_directory: Path | None = None) -> np.ndarray:
    """Boolean per perturbation of data: in the lockbox of lockbox_path (experiments/draw_lockbox.py). Refuses a lockbox drawn
    on another label selection or another leakage grouping, and one naming perturbations the data does not hold. With
    group_by and evidence_directory (the arguments the data was loaded with) it also refuses data loaded with another
    grouping, which a finer grouping would otherwise pass (its groups never straddle the lockbox, but its development
    folds would split drug-target groups), or from another evidence table."""
    lockbox = json.loads(Path(lockbox_path).read_text())
    if group_by is not None and group_by != lockbox["group_by"]:
        raise ValueError(f"{lockbox_path} was drawn with group_by {lockbox['group_by']!r}; the data was loaded with {group_by!r}")
    if evidence_directory is not None and lockbox.get("evidence_records_sha256"):
        evidence_sha256 = hashlib.sha256((Path(evidence_directory) / "evidence_records.parquet").read_bytes()).hexdigest()
        if evidence_sha256 != lockbox["evidence_records_sha256"]:
            raise ValueError(f"{lockbox_path} was drawn on evidence {lockbox['evidence_records_sha256'][:8]}, not on {evidence_directory} ({evidence_sha256[:8]})")
    selection_sha256 = (data.label_selection_summary or {}).get("sha256")
    if lockbox["label_selection_sha256"] != selection_sha256:
        raise ValueError(f"{lockbox_path} was drawn on label selection {lockbox['label_selection_sha256'][:8]}, not on the one loaded ({str(selection_sha256)[:8]})")
    held_out = set(lockbox["perturbation_ids"])
    missing = held_out - set(data.perturbation_ids)
    if missing:
        raise ValueError(f"{lockbox_path}: {len(missing)} lockbox perturbations are not in the data, e.g. {sorted(missing)[:3]}")
    in_lockbox = np.array([perturbation_id in held_out for perturbation_id in data.perturbation_ids])
    groups_inside = {group for group, inside in zip(data.group_ids, in_lockbox) if inside}
    straddling = [group for group, inside in zip(data.group_ids, in_lockbox) if not inside and group in groups_inside]
    if straddling:
        raise ValueError(f"{lockbox_path}: leakage group {straddling[0]} has members on both sides; load the data with group_by {lockbox['group_by']!r}")
    return in_lockbox


def merge_drugs_with_their_targets(perturbation_ids: list[str], perturbation_types: list[str], group_ids: list[str],
                                   perturbation_seeds: list, node_ids: list[str]) -> list[str]:
    """Leakage groups in which a drug is held out together with every labelled gene it targets and every drug that shares
    a target node with it (group_by "disease_cluster_and_targets").

    Under the disease-cluster grouping a drug's group is its target set, so a drug and the loss of function of its target
    gene can sit in different folds: on evidence_full with better_v1, 71 of the 72 drugs that target a labelled gene had a
    target gene in another fold, and 59 of the 340 kept drug positives were also positives of such a gene (docs/drug_target_leakage.md).
    A model then learns "perturbing this node causes this symptom" from the gene and is credited for it on the drug.
    Groups are the connected components of: perturbation - its group (disease cluster), drug - each gene perturbation
    whose node it targets, drug - each target node. Each component is named by the smallest disease cluster of a gene in it, or the smallest
    group id when it holds drugs only."""
    parent: dict[str, str] = {}

    def find(item: str) -> str:
        parent.setdefault(item, item)
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(first: str, second: str) -> None:
        parent[find(first)] = find(second)

    gene_perturbation_of_node = {node_ids[seeds[0]]: perturbation_id for perturbation_id, kind, seeds in zip(perturbation_ids, perturbation_types, perturbation_seeds)
                                 if kind == "gene" and len(seeds)}
    for perturbation_id, group_id in zip(perturbation_ids, group_ids):
        union("perturbation:" + perturbation_id, "group:" + group_id)
    for perturbation_id, kind, seeds in zip(perturbation_ids, perturbation_types, perturbation_seeds):
        if kind != "drug":
            continue
        for seed in seeds:
            union("perturbation:" + perturbation_id, "node:" + node_ids[seed])
            target_gene = gene_perturbation_of_node.get(node_ids[seed])
            if target_gene is not None:
                union("perturbation:" + perturbation_id, "perturbation:" + target_gene)
    name_of_root: dict[str, tuple[bool, str]] = {}
    for group_id, kind in zip(group_ids, perturbation_types):
        root = find("group:" + group_id)
        candidate = (kind != "gene", group_id)  # a gene's disease cluster names the component when it has one
        if root not in name_of_root or candidate < name_of_root[root]:
            name_of_root[root] = candidate
    return [name_of_root[find("perturbation:" + perturbation_id)][1] for perturbation_id in perturbation_ids]


GENE_TO_PROTEIN_FILE = "gene_to_protein.parquet"  # written by experiments/build_gene_protein_split.py; absent from merged graphs


def read_protein_of_gene(graph_directory: Path) -> dict[str, list[str]]:
    """gene node id -> the protein node ids a drug acting on the gene seeds, on a split graph (docs/gene_protein_split.md):
    its one protein node, or for a gene with several the ones drug targets name (the drug_target column); empty for a
    merged graph."""
    path = Path(graph_directory) / GENE_TO_PROTEIN_FILE
    if not path.exists():
        return {}
    table = pd.read_parquet(path)
    if "drug_target" in table.columns:
        table = table[table.drug_target]
    return table.groupby("gene_node_id", sort=False).protein_node_id.agg(list).to_dict()


def drug_seeds_on_proteins(triples: list, protein_of_gene: dict[str, list[str]]) -> list:
    """A drug acts on proteins: each target gene node is replaced by its protein nodes, with the target's sign and
    magnitude. When two target genes encode one shared protein node, the first triple is kept, so the node is seeded
    once."""
    seen, result = set(), []
    for node_id, sign, magnitude in triples:
        for protein_node_id in protein_of_gene.get(node_id, [node_id]):
            if protein_node_id not in seen:
                seen.add(protein_node_id)
                result.append([protein_node_id, sign, magnitude])
    return result


DEFAULT_LABEL_GRADES: tuple[str, ...] = ("A", "B")  # grades that count as positive labels; grade C (human association) and lower are soft evidence, not labels
SOFT_PRIOR_ONLY_GRADES: tuple[str, ...] = ("D", "E")  # literature grades: soft priors by design (section 4.2), never labels and never evaluation positives


def load_experiment_data(graph_directory: Path, evidence_directory: Path, relation: str = "induces", symptoms: list[str] | None = None, metabolic_layer_only: bool = False, group_by: str = "gene", label_grades: tuple[str, ...] | None = DEFAULT_LABEL_GRADES,
                         label_selection: Path | None = None) -> ExperimentData:
    """group_by selects the leakage group for the grouped split: "gene" (the gene itself; drugs by dominant
    target), "disease_cluster" (genes sharing a disease entry in HPO are held out together) or "disease_cluster_and_targets"
    (disease clusters with each drug joined to the genes it targets and to the drugs sharing a target: merge_drugs_with_their_targets). label_grades restricts
    the rows that become positive labels (default A and B; None keeps every grade, which the version 0.3 ablation
    table of grade C weight-0 positives relies on).

    label_selection is a parquet of (perturbation_id, symptom, keep) written by experiments/build_label_selection.py.
    A positive pair with keep False stays a positive in outcomes but is masked out (label_mask False): the trainer gives
    it zero weight and the metrics leave it out, so it is neither a positive nor a negative. Keeping outcome 1 means a
    code path that ignores the mask behaves as it did before the selection, never as if the pair were negative. A row
    with masks_a_negative True (build_label_selection.py --mask-grades) masks a pair that is not a positive (grade C
    evidence only) in the same way."""
    if group_by not in GROUPING_COLUMNS:
        raise ValueError(f"group_by must be one of {sorted(GROUPING_COLUMNS)}")
    group_column = GROUPING_COLUMNS[group_by]
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    relation_types = json.loads((graph_directory / "relation_types.json").read_text())
    protein_of_gene = read_protein_of_gene(graph_directory)
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
    group_seeds = {}  # the seeds the leakage groups are built from: the evidence's own gene nodes, also on a split graph
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
            triples = [triple for triple in json.loads(row.perturbation_nodes) if triple[0] in node_index]
            group_seeds[row.perturbation_id] = np.array([node_index[node_id] for node_id, _, _ in triples], dtype=int)
            if protein_of_gene and row.perturbation_type == "drug":  # split graph: a knockout seeds its gene, a drug its targets' proteins
                triples = drug_seeds_on_proteins(triples, protein_of_gene)
            seeds[row.perturbation_id] = np.array([node_index[node_id] for node_id, _, _ in triples], dtype=int)
            signs[row.perturbation_id] = np.array([sign for _, sign, _ in triples])
            magnitudes[row.perturbation_id] = np.array([magnitude for _, _, magnitude in triples])
    label_mask, label_selection_summary = None, None
    if label_selection is not None:
        selection = pd.read_parquet(label_selection)
        label_mask = np.ones_like(outcomes, dtype=bool)
        perturbation_row = {perturbation_id: index for index, perturbation_id in enumerate(perturbation_ids)}
        masks_a_negative = selection.masks_a_negative.astype(bool) if "masks_a_negative" in selection.columns else pd.Series(False, index=selection.index)
        set_aside = selection[~selection.keep.astype(bool) & ~masks_a_negative]
        masked, not_positive, masked_negatives = 0, 0, 0
        for row in set_aside.itertuples(index=False):
            position, column = perturbation_row.get(row.perturbation_id), symptom_index.get(row.symptom)
            if position is None or column is None or outcomes[position, column] == 0:
                not_positive += 1
                continue
            label_mask[position, column] = False
            masked += 1
        for row in selection[masks_a_negative].itertuples(index=False):  # --mask-grades rows: a pair with grade C evidence only
            position, column = perturbation_row.get(row.perturbation_id), symptom_index.get(row.symptom)
            if position is not None and column is not None and outcomes[position, column] == 0:
                label_mask[position, column] = False
                masked_negatives += 1
        selection = selection[~masks_a_negative]
        selected_pairs = set(zip(selection.perturbation_id, selection.symptom))
        positive_rows, positive_columns = np.nonzero(outcomes)
        without_a_row = sum((perturbation_ids[r], symptoms[c]) not in selected_pairs for r, c in zip(positive_rows.tolist(), positive_columns.tolist()))
        if without_a_row:  # kept, as before, but the selection did not judge them
            warnings.warn(f"{without_a_row} positive pairs have no row in {label_selection} and are kept unjudged")
        label_selection_summary = {"path": str(label_selection), "sha256": hashlib.sha256(Path(label_selection).read_bytes()).hexdigest(),
                                   "positive_pairs": int(outcomes.sum()), "masked_pairs": masked,
                                   "selection_rows_not_matching_a_positive": not_positive, "kept_positive_pairs": int((outcomes * label_mask).sum()),
                                   "positive_pairs_without_a_selection_row": int(without_a_row)}
        if masks_a_negative.any():
            label_selection_summary["masked_negative_pairs"] = masked_negatives
    group_ids = [groups[p] for p in perturbation_ids]
    if group_by == "disease_cluster_and_targets":
        group_ids = merge_drugs_with_their_targets(perturbation_ids, [types[p] for p in perturbation_ids], group_ids, [group_seeds[p] for p in perturbation_ids], node_ids)
    return ExperimentData(
        node_ids=node_ids,
        node_index=node_index,
        node_types=nodes.node_type.to_numpy(),
        is_currency=nodes.is_currency.fillna(False).to_numpy().astype(bool),
        node_degree=nodes.degree.to_numpy().astype(float),
        edge_source=edges.source_id.map(node_index).to_numpy(),
        edge_target=edges.target_id.map(node_index).to_numpy(),
        edge_relation=edges.relation_type.map(relation_index).to_numpy(),
        edge_sign=(edges["sign"].fillna(1.0) if "sign" in edges.columns else pd.Series(1.0, index=edges.index)).to_numpy(dtype=float),
        relation_types=relation_types,
        symptoms=symptoms,
        perturbation_ids=perturbation_ids,
        perturbation_labels=[labels[p] for p in perturbation_ids],
        perturbation_types=[types[p] for p in perturbation_ids],
        group_ids=group_ids,
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
        node_base_metabolite_id=nodes.base_metabolite_id.fillna("").to_numpy().astype(str) if "base_metabolite_id" in nodes.columns else None,
        node_display_name=nodes.display_name.fillna("").to_numpy().astype(str) if "display_name" in nodes.columns else None,
        label_mask=label_mask,
        label_selection_summary=label_selection_summary,
    )
