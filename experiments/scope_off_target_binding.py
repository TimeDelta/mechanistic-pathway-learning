"""Scope off-target binding for the drugs of an evidence table: how many measured off-targets each drug has under an
occupancy rule, how many carry an action (sign), how many are graph nodes, and how far the leakage groups merge when
every kept off-target joins its drug's group. Nothing is trained and no input is changed (docs/off_target_scoping.md).

Sources (downloaded into their own directories under data/raw, which git ignores):
  - DrugCentral drug.target.interaction.tsv.gz (release 2021_09_01, the latest flat file; newer releases ship only as
    a database dump): one row per drug, target and measurement, ACT_VALUE as -log10 molar, ACTION_TYPE mostly on
    mechanism rows, MOA = 1 for the mechanism of action.
  - IUPHAR/BPS Guide to PHARMACOLOGY interactions.csv and ligand_id_mapping.csv (release 2026.3): affinities with the
    ligand's type and action (full, partial or inverse agonist, antagonist, inhibitor, allosteric modulator, ...).

Per drug and target gene the affinity is the maximum over every measurement of that pair in both sources (the same
drug at the same target, never another ligand's value; the user's decision of 9 October 2026). A pair with more than
one value holds two curators' aggregates of different papers rather than one assay run twice (DrugCentral keeps one
value per pair 96.5 percent of the time, GtoPdb already medians its own papers), and by Cheng-Prusoff an IC50-derived
pK is biased downward against a Ki, so the larger value is the one closer to the dissociation constant. Occupancy assumes the therapeutic free concentration
occupies 90 percent of the drug's best-bound mechanism target, C = 9 K_primary, and Langmuir binding at one site:
occupancy = C / (C + K_target) = 1 / (1 + 10 ** (pK_primary - pK_target) / 9).

    python experiments/scope_off_target_binding.py
"""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, merge_drugs_with_their_targets

DRUGCENTRAL_INTERACTIONS = Path("data/raw/drugcentral_dti_2021_09_01/drug.target.interaction.tsv.gz")
GTOPDB_INTERACTIONS = Path("data/raw/gtopdb_2026/interactions.csv")
GTOPDB_LIGAND_MAPPING = Path("data/raw/gtopdb_2026/ligand_id_mapping.csv")
PRIMARY_OCCUPANCY = 0.9
OCCUPANCY_THRESHOLDS = (0.1, 0.25, 0.5)
AFFINITY_TYPES = {"Ki", "Kd", "IC50", "EC50", "Kb", "A2", "AC50"}
NEGATIVE_ACTIONS = re.compile(r"inhibit|antagon|block|inverse|negative|gating", re.IGNORECASE)
POSITIVE_ACTIONS = re.compile(r"(?<!ant)agonist|activat|opener|positive|potentiat|releasing", re.IGNORECASE)  # "antagonist" holds "agonist"


def occupancy(primary_p_affinity: float, target_p_affinity: float) -> float:
    """Fraction of the target bound when the primary target is PRIMARY_OCCUPANCY bound (one-site Langmuir binding)."""
    concentration_over_primary_k = PRIMARY_OCCUPANCY / (1 - PRIMARY_OCCUPANCY)
    return 1.0 / (1.0 + 10 ** (primary_p_affinity - target_p_affinity) / concentration_over_primary_k)


def action_sign(action: str) -> tuple[int | None, bool]:
    """(-1, +1 or None when unknown, is partial agonist)."""
    if not isinstance(action, str) or not action.strip():
        return None, False
    partial = "partial" in action.lower()
    if re.search("inverse", action, re.IGNORECASE):
        return -1, False
    if NEGATIVE_ACTIONS.search(action) and not POSITIVE_ACTIONS.search(action):
        return -1, partial
    if POSITIVE_ACTIONS.search(action) and not NEGATIVE_ACTIONS.search(action):
        return 1, partial
    return None, partial


def normalised_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def drugcentral_measurements(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    table = table[(table.ORGANISM == "Homo sapiens") & table.GENE.notna()]
    rows = []
    for record in table.itertuples(index=False):
        for gene in str(record.GENE).split("|"):
            rows.append({"drug_name": normalised_name(record.DRUG_NAME), "drugcentral_id": int(record.STRUCT_ID), "gene": gene.strip(),
                         "p_affinity": record.ACT_VALUE if record.ACT_TYPE in AFFINITY_TYPES else np.nan,
                         "action": record.ACTION_TYPE if isinstance(record.ACTION_TYPE, str) else "", "mechanism": record.MOA == 1,
                         "source": "DrugCentral:" + str(record.ACT_SOURCE)})
    return pd.DataFrame(rows)


def gtopdb_measurements(interactions_path: Path, mapping_path: Path) -> pd.DataFrame:
    table = pd.read_csv(interactions_path, skiprows=1, low_memory=False)
    table = table[(table["Target Species"] == "Human") & table["Target Gene Symbol"].notna()]
    mapping = pd.read_csv(mapping_path, skiprows=1, low_memory=False).set_index("Ligand id")
    median = pd.to_numeric(table["Affinity Median"], errors="coerce")
    middle = (pd.to_numeric(table["Affinity High"], errors="coerce") + pd.to_numeric(table["Affinity Low"], errors="coerce")) / 2
    p_affinity = median.fillna(middle)
    rows = []
    for record, value in zip(table.itertuples(index=False), p_affinity):
        ligand_id = record[table.columns.get_loc("Ligand ID")]
        mapped = mapping.loc[ligand_id] if ligand_id in mapping.index else None
        action = " ".join(str(part) for part in (record[table.columns.get_loc("Type")], record[table.columns.get_loc("Action")]) if isinstance(part, str))
        for gene in str(record[table.columns.get_loc("Target Gene Symbol")]).split("|"):
            rows.append({"drug_name": normalised_name(record[table.columns.get_loc("Ligand")]),
                         "pubchem_cid": None if mapped is None or pd.isna(mapped["PubChem CID"]) else int(mapped["PubChem CID"]),
                         "drugcentral_id": None if mapped is None or pd.isna(mapped["Drug Central ID"]) else int(mapped["Drug Central ID"]),
                         "inn": "" if mapped is None or pd.isna(mapped["INN"]) else normalised_name(mapped["INN"]),
                         "gene": gene.strip(), "p_affinity": value, "action": action,
                         "mechanism": str(record[table.columns.get_loc("Primary Target")]).lower() == "true", "source": "GtoPdb"})
    return pd.DataFrame(rows)


def measurement_coverage(per_drug: list[dict]) -> dict:
    """How unevenly the drugs were measured: human genes with an affinity per matched drug, and the drugs with rows
    from the two broad panels in DrugCentral (DrugMatrix, a fixed panel run on every compound it holds, and the NIMH
    Psychoactive Drug Screening Program)."""
    genes_with_affinity = np.array([entry["genes_with_affinity"] for entry in per_drug])
    return {"human_genes_with_an_affinity_per_drug": {"median": float(np.median(genes_with_affinity)),
                                                        "lower_quartile": float(np.quantile(genes_with_affinity, 0.25)),
                                                        "upper_quartile": float(np.quantile(genes_with_affinity, 0.75)),
                                                        "max": int(genes_with_affinity.max())},
            "drugs_with_drug_matrix_rows": sum("DrugCentral:DRUG MATRIX" in entry["measurement_sources"] for entry in per_drug),
            "drugs_with_pdsp_rows": sum("DrugCentral:PDSP" in entry["measurement_sources"] for entry in per_drug)}


def measurement_coverage_table(per_drug: list[dict], data, thresholds: tuple[float, ...]) -> pd.DataFrame:
    """One row per perturbation saying how well its drug was measured, for stratified scoring.

    A drug that has been screened against a whole panel shows more off-targets than one measured at its mechanism
    target alone, whatever its pharmacology, so any reading that uses measured off-targets has to be reported inside
    strata of how much measurement the drug has. A gene perturbation has no drug measurement and is its own stratum; a
    drug with no row in either source is "unmeasured".
    """
    by_position = {entry["position"]: entry for entry in per_drug}
    rows = []
    for position, (perturbation_id, kind, label) in enumerate(zip(data.perturbation_ids, data.perturbation_types, data.perturbation_labels)):
        entry = by_position.get(position)
        row = {"perturbation_id": perturbation_id, "perturbation_type": kind, "label": label,
               "genes_with_an_affinity": -1 if entry is None else entry["genes_with_affinity"],
               "genes_without_an_affinity": -1 if entry is None else entry["genes_without_affinity"],
               "primary_p_affinity": None if entry is None else entry["primary_p_affinity"],
               "measurement_sources": "" if entry is None else "|".join(sorted(entry["measurement_sources"])),
               "has_drug_matrix_rows": False if entry is None else "DrugCentral:DRUG MATRIX" in entry["measurement_sources"],
               "has_pdsp_rows": False if entry is None else "DrugCentral:PDSP" in entry["measurement_sources"]}
        for threshold in thresholds:
            kept = None if entry is None else entry["off_targets"][entry["off_targets"].occupancy >= threshold]
            row[f"off_targets_at_{threshold}"] = -1 if kept is None else int(len(kept))
            row[f"off_targets_in_the_graph_at_{threshold}"] = -1 if kept is None else int(kept.graph_node.notna().sum())
        rows.append(row)
    return pd.DataFrame(rows)


def largest_group(groups: list[str], kept_positive: np.ndarray, perturbation_types: list[str]) -> dict:
    """The largest leakage group by perturbations (as in docs/drug_targets_any_type.md) and the group holding the most
    kept positive pairs, with the share each takes."""
    positives, drugs, members = Counter(), Counter(), Counter()
    for group, count, kind in zip(groups, kept_positive, perturbation_types):
        positives[group] += int(count)
        drugs[group] += kind == "drug"
        members[group] += 1
    by_size, size = members.most_common(1)[0]
    by_positives, most_positives = positives.most_common(1)[0]
    total_drugs = sum(kind == "drug" for kind in perturbation_types)
    by_drugs, most_drugs = drugs.most_common(1)[0]
    return {"leakage_groups": len(members), "drug_holding_groups": sum(count > 0 for count in drugs.values()),
            "largest_by_drugs": {"drugs": int(most_drugs), "share_of_drugs": round(most_drugs / total_drugs, 3)},
            "largest_by_perturbations": {"perturbations": int(size), "share_of_perturbations": round(size / len(groups), 3), "drugs": int(drugs[by_size]),
                                         "share_of_drugs": round(drugs[by_size] / total_drugs, 3)},
            "most_kept_positives": {"perturbations": int(members[by_positives]), "drugs": int(drugs[by_positives]),
                                    "share_of_kept_positive_pairs": round(most_positives / max(int(kept_positive.sum()), 1), 3)}}


def join_off_target_genes_only(mechanism_groups: list[str], perturbation_ids: list[str], extra_nodes_of_drug: dict[int, list[str]],
                               gene_perturbation_of_node: dict[str, int]) -> list[str]:
    """Groups of mechanism targets (merge_drugs_with_their_targets), then each drug joined to the labelled gene
    perturbation of each of its off-target genes; two drugs sharing only an off-target are not joined."""
    parent = {group: group for group in mechanism_groups}

    def find(group: str) -> str:
        while parent[group] != group:
            parent[group] = parent[parent[group]]
            group = parent[group]
        return group

    for drug_position, extra_nodes in extra_nodes_of_drug.items():
        for node in extra_nodes:
            gene_position = gene_perturbation_of_node.get(node)
            if gene_position is not None:
                parent[find(mechanism_groups[drug_position])] = find(mechanism_groups[gene_position])
    return [find(group) for group in mechanism_groups]


def join_drugs_to_target_genes_only(disease_cluster_groups: list[str], target_nodes_of_drug: dict[int, list[str]],
                                    gene_perturbation_of_node: dict[str, int]) -> list[str]:
    """Disease-cluster groups (a drug's own group is its target set) with each drug joined to the labelled gene
    perturbation of each gene it binds; drugs that share a target but no labelled gene stay apart."""
    parent = {group: group for group in disease_cluster_groups}

    def find(group: str) -> str:
        while parent[group] != group:
            parent[group] = parent[parent[group]]
            group = parent[group]
        return group

    for drug_position, nodes in target_nodes_of_drug.items():
        for node in nodes:
            gene_position = gene_perturbation_of_node.get(node)
            if gene_position is not None:
                parent[find(disease_cluster_groups[drug_position])] = find(disease_cluster_groups[gene_position])
    return [find(group) for group in disease_cluster_groups]


def drug_keys(perturbation_id: str, label: str) -> tuple[int | None, str]:
    match = re.fullmatch(r"CID[01](\d+)", perturbation_id)
    return (int(match.group(1)) if match else None), normalised_name(label)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v3_max3_targets"))
    parser.add_argument("--label-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_max3_targets.parquet"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/off_target_scoping/summary.json"))
    parser.add_argument("--coverage-table", type=Path, default=Path("data/processed/off_target_scoping/measurement_coverage.parquet"),
                        help="per-perturbation measurement coverage, read by the stratified scoring")
    arguments = parser.parse_args()

    data = load_experiment_data(arguments.graph_dir, arguments.evidence_dir, group_by="disease_cluster", label_selection=arguments.label_selection)
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet", columns=["node_id", "node_type", "gene_symbol"])
    gene_node_of_symbol = dict(zip(nodes.gene_symbol[nodes.node_type == "gene"], nodes.node_id[nodes.node_type == "gene"]))
    node_index = {node_id: index for index, node_id in enumerate(data.node_ids)}
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))

    drugcentral = drugcentral_measurements(DRUGCENTRAL_INTERACTIONS)
    gtopdb = gtopdb_measurements(GTOPDB_INTERACTIONS, GTOPDB_LIGAND_MAPPING)
    drugcentral_ids_of_name = drugcentral.groupby("drug_name").drugcentral_id.first().to_dict()

    per_drug, unmatched = [], []
    for position, (perturbation_id, kind, label, seeds) in enumerate(zip(data.perturbation_ids, data.perturbation_types, data.perturbation_labels, data.perturbation_seeds)):
        if kind != "drug":
            continue
        pubchem_cid, name = drug_keys(perturbation_id, label)
        from_gtopdb = gtopdb[((gtopdb.pubchem_cid == pubchem_cid) if pubchem_cid else False) | (gtopdb.drug_name == name) | (gtopdb.inn == name)]
        drugcentral_ids = set(from_gtopdb.drugcentral_id.dropna().astype(int)) | ({drugcentral_ids_of_name[name]} if name in drugcentral_ids_of_name else set())
        from_drugcentral = drugcentral[drugcentral.drugcentral_id.isin(drugcentral_ids)]
        measurements = pd.concat([from_drugcentral, from_gtopdb], ignore_index=True)
        if measurements.empty:
            unmatched.append(label)
            continue
        mechanism_genes = {symbol_of_node.get(data.node_ids[seed]) for seed in seeds} - {None}
        by_gene = measurements.groupby("gene").agg(p_affinity=("p_affinity", "max"), measurements=("p_affinity", "count"),
                                                   actions=("action", lambda values: sorted({value for value in values if value})),
                                                   listed_as_mechanism=("mechanism", "any"))
        with_affinity = by_gene[by_gene.p_affinity.notna()]
        primary_candidates = with_affinity[with_affinity.index.isin(mechanism_genes)]
        primary_source = "ChEMBL mechanism target"
        if primary_candidates.empty:
            primary_candidates, primary_source = with_affinity[with_affinity.listed_as_mechanism], "source's own mechanism flag"
        primary = float(primary_candidates.p_affinity.max()) if not primary_candidates.empty else None
        off_targets = with_affinity[~with_affinity.index.isin(mechanism_genes)].copy()
        off_targets["occupancy"] = [occupancy(primary, value) if primary is not None else np.nan for value in off_targets.p_affinity]
        signs = [action_sign(" ".join(actions)) for actions in off_targets.actions]
        off_targets["sign"] = [sign for sign, _ in signs]
        off_targets["partial_agonist"] = [partial for _, partial in signs]
        off_targets["graph_node"] = [gene_node_of_symbol.get(gene) for gene in off_targets.index]
        per_drug.append({"position": position, "perturbation_id": perturbation_id, "label": label, "mechanism_genes": sorted(mechanism_genes),
                         "primary_p_affinity": primary, "primary_from": primary_source if primary is not None else None,
                         "genes_with_affinity": int(len(with_affinity)), "genes_without_affinity": int(len(by_gene) - len(with_affinity)),
                         "measurement_sources": set(measurements.source), "off_targets": off_targets})

    summary = {"evidence_dir": str(arguments.evidence_dir), "drugs": int(sum(kind == "drug" for kind in data.perturbation_types)),
               "drugs_matched": len(per_drug), "drugs_unmatched": sorted(unmatched),
               "drugs_with_primary_affinity": sum(entry["primary_p_affinity"] is not None for entry in per_drug),
               "primary_from": dict(Counter(entry["primary_from"] for entry in per_drug)),
               "sources": {"DrugCentral": str(DRUGCENTRAL_INTERACTIONS), "GtoPdb": str(GTOPDB_INTERACTIONS)},
               "measurement_coverage": measurement_coverage(per_drug), "by_threshold": {}}
    base_groups = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids, data.perturbation_seeds, data.node_ids)
    kept_positive = ((data.outcomes > 0) & (data.label_mask > 0 if data.label_mask is not None else True)).sum(axis=1)
    gene_perturbation_of_node = {data.node_ids[seeds[0]]: position for position, (kind, seeds) in enumerate(zip(data.perturbation_types, data.perturbation_seeds))
                                 if kind == "gene" and len(seeds)}
    for threshold in OCCUPANCY_THRESHOLDS:
        counts, signed, partial, in_graph, total_pairs = [], 0, 0, 0, 0
        seeds_with_off_targets = list(data.perturbation_seeds)
        extra_nodes_of_drug = {}
        for entry in per_drug:
            kept = entry["off_targets"][entry["off_targets"].occupancy >= threshold]
            counts.append(len(kept))
            total_pairs += len(kept)
            signed += int(kept.sign.notna().sum())
            partial += int(kept.partial_agonist.sum())
            in_graph += int(kept.graph_node.notna().sum())
            extra_nodes_of_drug[entry["position"]] = [node for node in kept.graph_node.dropna() if node in node_index]
            extra = [node_index[node] for node in extra_nodes_of_drug[entry["position"]]]
            seeds_with_off_targets[entry["position"]] = np.concatenate([np.asarray(data.perturbation_seeds[entry["position"]], dtype=int), np.asarray(extra, dtype=int)])
        every_off_target_joins = merge_drugs_with_their_targets(data.perturbation_ids, data.perturbation_types, data.group_ids, seeds_with_off_targets, data.node_ids)
        genes_only = join_off_target_genes_only(base_groups, data.perturbation_ids, extra_nodes_of_drug, gene_perturbation_of_node)
        every_target_node = {position: [data.node_ids[seed] for seed in seeds_with_off_targets[position]] for position in extra_nodes_of_drug}
        drug_to_gene_only = join_drugs_to_target_genes_only(list(data.group_ids), every_target_node, gene_perturbation_of_node)
        drugs_joined_to_a_gene = sum(any(node in gene_perturbation_of_node for node in nodes_) for nodes_ in extra_nodes_of_drug.values())
        summary["by_threshold"][str(threshold)] = {
            "off_target_pairs": total_pairs, "with_action": signed, "partial_agonist": partial, "gene_node_in_graph": in_graph,
            "off_targets_per_drug_median": float(np.median(counts)), "off_targets_per_drug_max": int(max(counts)),
            "drugs_with_any": int(sum(count > 0 for count in counts)), "drugs_with_an_off_target_that_is_a_labelled_gene": int(drugs_joined_to_a_gene),
            "every_off_target_joins": largest_group(every_off_target_joins, kept_positive, data.perturbation_types),
            "off_target_joins_only_its_gene_perturbation": largest_group(genes_only, kept_positive, data.perturbation_types),
            "drugs_joined_only_to_labelled_genes_they_bind": largest_group(drug_to_gene_only, kept_positive, data.perturbation_types)}
    summary["without_off_targets"] = largest_group(base_groups, kept_positive, data.perturbation_types)
    mechanism_nodes = {position: [data.node_ids[seed] for seed in data.perturbation_seeds[position]]
                       for position, kind in enumerate(data.perturbation_types) if kind == "drug"}
    summary["without_off_targets_drugs_joined_only_to_labelled_genes"] = largest_group(
        join_drugs_to_target_genes_only(list(data.group_ids), mechanism_nodes, gene_perturbation_of_node), kept_positive, data.perturbation_types)
    examples = {}
    for name in ("amitriptyline", "sertraline", "haloperidol", "lithium", "methylphenidate"):
        entry = next((entry for entry in per_drug if entry["label"].lower() == name), None)
        if entry is not None:
            top = entry["off_targets"].sort_values("occupancy", ascending=False).head(6)
            examples[name] = {"mechanism_genes": entry["mechanism_genes"], "primary_p_affinity": entry["primary_p_affinity"],
                              "top_off_targets": [{"gene": gene, "p_affinity": round(float(row.p_affinity), 2), "occupancy": round(float(row.occupancy), 2) if not np.isnan(row.occupancy) else None,
                                                   "sign": None if pd.isna(row.sign) else int(row.sign), "actions": row.actions} for gene, row in top.iterrows()]}
    summary["examples"] = examples
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    coverage = measurement_coverage_table(per_drug, data, OCCUPANCY_THRESHOLDS)
    arguments.coverage_table.parent.mkdir(parents=True, exist_ok=True)
    coverage.to_parquet(arguments.coverage_table, index=False)
    summary["coverage_table"] = {"path": str(arguments.coverage_table), "rows": int(len(coverage)),
                                 "drug_rows": int((coverage.perturbation_type == "drug").sum()),
                                 "drug_rows_unmeasured": int(((coverage.perturbation_type == "drug") & (coverage.genes_with_an_affinity < 0)).sum())}
    arguments.output.write_text(json.dumps(summary, indent=2, default=lambda value: None if value is None or (isinstance(value, float) and np.isnan(value)) else str(value)))
    print(json.dumps({key: value for key, value in summary.items() if key != "examples"}, indent=2, default=str))
    print(json.dumps(examples, indent=1, default=str)[:4000])


if __name__ == "__main__":
    main()
