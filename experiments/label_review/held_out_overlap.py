"""Which side of the held-out set would each added example fall on? (the user's rule of 10 October 2026)

The rule: an added example whose mechanism overlaps a held-out perturbation is held out with it, and one that does not
trains. Applied here to two kinds of addition, against the membership of configs/lockbox_v2.json:

  drugs             the drugs a wider drug rule admits (--new-evidence against --current-evidence);
  target statements the Open Targets statements "inhibition (or activation) of X causes symptom Y" that no study drug
                    label already makes (experiments/scope_open_targets_symptom_labels.py).

Overlap is tested on nodes, one step and not chained:

  1. the held-out set is carried to the current table by its leakage groups, as the registered rule has it (every member
     of a group that holds a lockbox perturbation is held out);
  2. an added drug is held out when it acts on a node that a carried held-out perturbation acts on or is;
  3. a target statement is held out when its target is such a node, or a node of a drug held out under 2;
  4. statements that share their direction, symptom and references and whose targets are subunits or members of one
     receptor (one ChEMBL target record lists both genes, or the gene symbols share their stem) were written as one
     sentence about that receptor, so if one of them is held out all of them are.

The registered rule chains instead: it joins every perturbation that shares a node with another. The script reports what
chaining would do to the wider table, because that is why it is not used here: one added drug that acts on a held-out
node and on a development node would pull the development node's whole group into the held-out set.

What step 2 leaves open, and the script counts: a held-out added drug with a second target shares that target with the
development drugs that act on it.

A graph distance is reported for reference only. On the confirmatory graph most development perturbations already sit
one or two edges from a held-out seed, so distance does not separate anything.

Only the membership of the lockbox is read (its perturbation ids), never a prediction or a score. Counting only: no label,
table or configuration of the study is changed.

Usage:
  python experiments/label_review/held_out_overlap.py
  python experiments/label_review/held_out_overlap.py --new-evidence DIR --new-selection FILE --name NAME
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from mechanistic_pathway_learning.evaluation.experiment_data import ExperimentData, load_experiment_data, merge_drugs_with_their_targets  # noqa: E402

STATEMENT_ROWS = Path("data/processed/off_target_scoping/open_targets_symptom_rows.parquet")


def seed_node_ids(data: ExperimentData, position: int) -> set[str]:
    return {data.node_ids[seed] for seed in data.perturbation_seeds[position]}


def carried_held_out(current: ExperimentData, lockbox_ids: set[str]) -> tuple[set[str], set[str]]:
    """(perturbation ids, seed node ids) of the held-out set on the current table: every member of a leakage group that holds a lockbox perturbation."""
    groups = merge_drugs_with_their_targets(current.perturbation_ids, current.perturbation_types, current.group_ids, current.perturbation_seeds, current.node_ids)
    held_groups = {group for perturbation_id, group in zip(current.perturbation_ids, groups) if perturbation_id in lockbox_ids}
    held = {perturbation_id for perturbation_id, group in zip(current.perturbation_ids, groups) if group in held_groups}
    nodes = set().union(*[seed_node_ids(current, position) for position, perturbation_id in enumerate(current.perturbation_ids) if perturbation_id in held])
    return held, nodes


def distance_to_nearest(source_nodes: set[int], neighbours: list[np.ndarray], blocked: np.ndarray) -> np.ndarray:
    """Edges to the nearest of source_nodes for every node, -1 when unreachable; a blocked node is never passed through."""
    distance = np.full(len(neighbours), -1, dtype=int)
    queue = deque()
    for node in source_nodes:
        distance[node] = 0
        queue.append(node)
    while queue:
        node = queue.popleft()
        if blocked[node] and distance[node] > 0:
            continue
        for neighbour in neighbours[node]:
            if distance[neighbour] < 0:
                distance[neighbour] = distance[node] + 1
                queue.append(neighbour)
    return distance


def symbol_stem(gene_symbol: str) -> str:
    """The letters of a gene symbol before its first digit (SCN for SCN1B, CHRM for CHRM2); "" when shorter than three letters."""
    stem = ""
    for character in gene_symbol:
        if character.isdigit():
            break
        stem += character
    return stem if len(stem) >= 3 else ""


def statement_families(statements: pd.DataFrame, chembl_targets: dict) -> list[str]:
    """One family name per statement. Two statements are one family when they share direction, symptom and references and
    their targets are subunits or members of one receptor: both gene symbols listed by one ChEMBL target record of several
    genes (a complex, a complex group or a protein family), or the same symbol stem. Chained, so A-B and B-C make A, B, C."""
    multi_gene_records = [set(record.get("gene_symbols") or []) for record in chembl_targets.values() if len(record.get("gene_symbols") or []) > 1]
    parent = list(range(len(statements)))

    def find(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    rows = list(statements.itertuples(index=False))
    for first in range(len(rows)):
        for second in range(first + 1, len(rows)):
            a, b = rows[first], rows[second]
            if (a.direction, a.target_symptom, a.references) != (b.direction, b.target_symptom, b.references) or a.target == b.target:
                continue
            same_stem = symbol_stem(a.target) != "" and symbol_stem(a.target) == symbol_stem(b.target)
            if same_stem or any(a.target in genes and b.target in genes for genes in multi_gene_records):
                parent[find(first)] = find(second)
    members: dict[int, list[str]] = defaultdict(list)
    for position, row in enumerate(rows):
        members[find(position)].append(row.target)
    return [f"{rows[position].direction}|{rows[position].target_symptom}|{'+'.join(sorted(set(members[find(position)])))}" for position in range(len(rows))]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lockbox", type=Path, default=Path("configs/lockbox_v2.json"))
    parser.add_argument("--grouping-graph", type=Path, default=Path("data/processed/graph_full_neuronal"), help="the graph the lockbox was drawn on")
    parser.add_argument("--distance-graph", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders"), help="the confirmatory graph")
    parser.add_argument("--current-evidence", type=Path, default=Path("data/processed/evidence_full_v3_parkinsonism"))
    parser.add_argument("--current-selection", type=Path, default=Path("data/processed/label_selection/better_v2_full_v3_parkinsonism.parquet"))
    parser.add_argument("--new-evidence", type=Path, default=Path("data/processed/label_review/evidence_n_only_cap3"))
    parser.add_argument("--new-selection", type=Path, default=Path("data/processed/label_review/selection_n_only_cap3.parquet"))
    parser.add_argument("--name", default="three_targets", help="names the output files")
    parser.add_argument("--chembl-dir", type=Path, default=Path("data/raw/chembl"), help="targets.json names the genes of each receptor complex and family")
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/label_review"))
    arguments = parser.parse_args()

    lockbox_ids = set(json.loads(arguments.lockbox.read_text())["perturbation_ids"])
    current = load_experiment_data(arguments.grouping_graph, arguments.current_evidence, group_by="disease_cluster", label_selection=arguments.current_selection)
    new = load_experiment_data(arguments.grouping_graph, arguments.new_evidence, group_by="disease_cluster", label_selection=arguments.new_selection)
    carried, carried_nodes = carried_held_out(current, lockbox_ids)
    current_types = dict(zip(current.perturbation_ids, current.perturbation_types))
    print(f"held-out set carried to the current table: {len(carried)} perturbations, {sum(current_types[p] == 'drug' for p in carried)} of its {current.perturbation_types.count('drug')} drugs "
          f"({len(lockbox_ids & set(current.perturbation_ids))} of the {len(lockbox_ids)} lockbox perturbations are in the table)")

    # ---- added drugs
    kept = (new.outcomes * new.label_mask)
    kept_per_perturbation = kept.sum(axis=1)
    in_current = set(current.perturbation_ids)
    missing_from_new = sorted(in_current - set(new.perturbation_ids))
    if missing_from_new:
        print(f"note: {len(missing_from_new)} perturbations of the current table are not in the new one, e.g. {missing_from_new[:5]}")
    added = [position for position, perturbation_id in enumerate(new.perturbation_ids) if perturbation_id not in in_current]
    added_held = [position for position in added if seed_node_ids(new, position) & carried_nodes]
    added_development = [position for position in added if position not in set(added_held)]
    held_positions = set(added_held) | {position for position, perturbation_id in enumerate(new.perturbation_ids) if perturbation_id in carried}
    print(f"\nadded perturbations: {len(added)} ({sum(new.perturbation_types[p] == 'drug' for p in added)} drugs), {int(sum(kept_per_perturbation[p] for p in added))} kept pairs")
    print(f"   acting on a node of a held-out perturbation, so held out: {len(added_held)}, {int(sum(kept_per_perturbation[p] for p in added_held))} kept pairs")
    print(f"   not, so they train: {len(added_development)}, {int(sum(kept_per_perturbation[p] for p in added_development))} kept pairs")
    print("   held out:", ", ".join(f"{new.perturbation_labels[p]} ({int(kept_per_perturbation[p])})" for p in sorted(added_held, key=lambda p: -kept_per_perturbation[p])))
    print("   train:", ", ".join(f"{new.perturbation_labels[p]} ({int(kept_per_perturbation[p])})" for p in sorted(added_development, key=lambda p: -kept_per_perturbation[p])))
    other_nodes = set().union(*[seed_node_ids(new, position) for position in added_held]) - carried_nodes if added_held else set()
    sharing = [position for position in range(len(new.perturbation_ids)) if position not in held_positions and seed_node_ids(new, position) & other_nodes]
    print(f"   other targets of the held-out added drugs: {len(other_nodes)} nodes {sorted(node.replace('GENE:', '') for node in other_nodes)}")
    print(f"   development perturbations that share one of those targets: {len(sharing)} ({sum(new.perturbation_types[p] == 'drug' for p in sharing)} drugs), {int(sum(kept_per_perturbation[p] for p in sharing))} kept pairs: "
          + ", ".join(sorted(new.perturbation_labels[p] for p in sharing)))
    held_out_drugs = [p for p in held_positions if new.perturbation_types[p] == "drug"]
    all_drugs = [p for p, kind in enumerate(new.perturbation_types) if kind == "drug"]
    print(f"   under this rule: {len(held_positions)} of {len(new.perturbation_ids)} perturbations held out; {len(held_out_drugs)} of {len(all_drugs)} drugs; "
          f"{int(kept_per_perturbation[sorted(held_positions)].sum())} of {int(kept_per_perturbation.sum())} kept pairs; {int(kept_per_perturbation[held_out_drugs].sum())} of {int(kept_per_perturbation[all_drugs].sum())} kept drug pairs")
    held_mask = np.zeros(len(new.perturbation_ids), dtype=bool)
    held_mask[sorted(held_positions)] = True
    print("   kept pairs by symptom (train / held out):", {symptom: (int(kept[~held_mask, column].sum()), int(kept[held_mask, column].sum())) for column, symptom in enumerate(new.symptoms)})

    chained = merge_drugs_with_their_targets(new.perturbation_ids, new.perturbation_types, new.group_ids, new.perturbation_seeds, new.node_ids)
    chained_held_groups = {group for perturbation_id, group in zip(new.perturbation_ids, chained) if perturbation_id in carried}
    pulled_in = [position for position, (perturbation_id, group) in enumerate(zip(new.perturbation_ids, chained)) if group in chained_held_groups and perturbation_id not in carried]
    print(f"   chaining instead (the registered group rule on the new table) would hold out {len(pulled_in)} more perturbations, "
          f"{sum(new.perturbation_types[p] == 'drug' for p in pulled_in)} of them drugs, with {int(kept_per_perturbation[pulled_in].sum())} kept pairs")

    # ---- target statements
    rows = pd.read_parquet(STATEMENT_ROWS)
    statements = rows.groupby(["target", "direction", "target_symptom"]).agg(
        animal_readout=("animal_readout", "all"), species_ambiguous=("species_ambiguous_readout", "all"),
        acute_only=("dosing", lambda values: set(values) == {"acute"}), references=("reference", lambda values: "; ".join(sorted(set(values))))).reset_index()
    statements["target_node"] = "GENE:" + statements.target
    selection = pd.read_parquet(arguments.new_selection)
    selection = selection[~selection.masks_a_negative.astype(bool)] if "masks_a_negative" in selection.columns else selection
    listed: dict[str, set[str]] = defaultdict(set)
    for row in selection.itertuples(index=False):
        listed[row.perturbation_id].add(row.symptom)
    drugs_on_node_in_direction: dict[tuple[str, str], list[int]] = defaultdict(list)
    perturbations_on_node: dict[str, list[int]] = defaultdict(list)
    gene_perturbation_at_node: dict[str, int] = {}
    for position, (kind, seeds, signs) in enumerate(zip(new.perturbation_types, new.perturbation_seeds, new.perturbation_signs)):
        for seed, sign in zip(seeds, signs):
            node_id = new.node_ids[seed]
            perturbations_on_node[node_id].append(position)
            if kind == "gene":
                gene_perturbation_at_node[node_id] = position
            elif sign != 0:
                drugs_on_node_in_direction[(node_id, "inhibition" if sign < 0 else "activation")].append(position)
    acting = [drugs_on_node_in_direction.get((row.target_node, row.direction), []) for row in statements.itertuples(index=False)]
    statements["made_by_a_study_drug_label"] = [any(row.target_symptom in listed[new.perturbation_ids[position]] for position in positions) for row, positions in zip(statements.itertuples(index=False), acting)]
    statements["study_drugs_in_the_statement_direction"] = [len(positions) for positions in acting]
    not_made = statements[~statements.made_by_a_study_drug_label]
    human = not_made[~not_made.animal_readout & ~not_made.species_ambiguous]
    usable = human[human.target_node.isin(new.node_index)].reset_index(drop=True)
    print(f"\ntarget statements: {len(statements)} | not made by a study drug label: {len(not_made)} | of them in human wording: {len(human)} | of them with the target in the graph: {len(usable)} over {usable.target.nunique()} targets"
          + (f" (not in the graph: {sorted(set(human.target) - set(usable.target))})" if len(usable) < len(human) else ""))
    held_nodes = carried_nodes | set().union(*[seed_node_ids(new, position) for position in added_held]) if added_held else set(carried_nodes)
    usable["held_out_by_node"] = usable.target_node.isin(held_nodes)
    family_key = pd.Series(statement_families(usable, json.loads((arguments.chembl_dir / "targets.json").read_text())), index=usable.index)
    usable["family"] = family_key
    usable["family_size"] = family_key.map(family_key.value_counts())
    usable["held_out"] = usable.held_out_by_node | family_key.map(usable.groupby(family_key).held_out_by_node.any())
    usable["target_is_a_gene_perturbation"] = usable.target_node.isin(gene_perturbation_at_node)
    usable["gene_perturbation_lists_the_symptom"] = [row.target_node in gene_perturbation_at_node and row.target_symptom in listed[new.perturbation_ids[gene_perturbation_at_node[row.target_node]]] for row in usable.itertuples(index=False)]
    usable["development_perturbations_on_the_target"] = [sum(position not in held_positions for position in perturbations_on_node.get(node_id, [])) for node_id in usable.target_node]

    nodes, edges = pd.read_parquet(arguments.distance_graph / "nodes.parquet"), pd.read_parquet(arguments.distance_graph / "edges.parquet")
    index_of = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    source, target = edges.source_id.map(index_of).to_numpy(), edges.target_id.map(index_of).to_numpy()
    starts_of_edges, ends_of_edges = np.concatenate([source, target]), np.concatenate([target, source])
    order = np.argsort(starts_of_edges, kind="stable")
    boundaries = np.searchsorted(starts_of_edges[order], np.arange(len(nodes) + 1))
    ends_in_start_order = ends_of_edges[order]
    neighbours = [ends_in_start_order[boundaries[index]: boundaries[index + 1]] for index in range(len(nodes))]
    blocked = nodes.is_currency.astype("boolean").fillna(False).to_numpy(dtype=bool)
    proteins_of_gene: dict[str, list[str]] = defaultdict(list)
    for row in pd.read_parquet(arguments.distance_graph / "gene_to_protein.parquet").itertuples(index=False):
        proteins_of_gene[row.gene_node_id].append(row.protein_node_id)

    def graph_nodes(node_id: str) -> set[int]:
        return {index_of[candidate] for candidate in [node_id, *proteins_of_gene.get(node_id, [])] if candidate in index_of}

    distance = distance_to_nearest(set().union(*[graph_nodes(node_id) for node_id in held_nodes]), neighbours, blocked)

    def nearest(node_ids: set[str]) -> int:
        reachable = [int(distance[index]) for node_id in node_ids for index in graph_nodes(node_id) if distance[index] >= 0]
        return min(reachable) if reachable else -1

    usable["edges_to_nearest_held_out_seed"] = [nearest({node_id}) for node_id in usable.target_node]
    development_distance = Counter(nearest(seed_node_ids(new, position)) for position in range(len(new.perturbation_ids)) if position not in held_positions and len(new.perturbation_seeds[position]))
    training, held = usable[~usable.held_out], usable[usable.held_out]
    print(f"   held out: {len(held)} statements over {held.target.nunique()} targets ({int(held.held_out_by_node.sum())} by their own target, {int((~held.held_out_by_node).sum())} with their family)")
    print(f"   train: {len(training)} statements over {training.target.nunique()} targets; activation {int((training.direction == 'activation').sum())}, inhibition {int((training.direction == 'inhibition').sum())}; acute dosing only {int(training.acute_only.sum())}")
    print("   train, by symptom:", training.target_symptom.value_counts().to_dict())
    print("   held out, by symptom:", held.target_symptom.value_counts().to_dict())
    print("   families of more than one target (family, held out):", [(family, bool(group.held_out.any())) for family, group in usable.groupby(family_key) if len(group) > 1])
    print(f"   independent statements (families): train {training.family.nunique()}, held out {held.family.nunique()}")
    print("   edges to the nearest held-out seed, statements that train:", dict(sorted(Counter(training.edges_to_nearest_held_out_seed).items())),
          "| development perturbations of the study:", dict(sorted(development_distance.items())))
    contradicted = usable[usable.study_drugs_in_the_statement_direction > 0]
    print(f"   statements whose target a study drug acts on in that direction without listing the symptom: {len(contradicted)} "
          f"({int((~contradicted.held_out).sum())} of them train): {[(row.target, row.direction, row.target_symptom) for row in contradicted.itertuples(index=False)]}")
    lowered_twice = training[(training.direction == "inhibition") & training.target_is_a_gene_perturbation]
    print(f"   inhibition statements that train and whose target is also a gene perturbation of the study: {len(lowered_twice)}; the gene's own labels list the symptom in {int(lowered_twice.gene_perturbation_lists_the_symptom.sum())}")
    print(f"   statements that train on a target no development perturbation acts on or is: {int((training.development_perturbations_on_the_target == 0).sum())}")
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    usable.drop(columns=["animal_readout", "species_ambiguous", "made_by_a_study_drug_label"]).to_csv(arguments.output_dir / f"target_statements_by_side_{arguments.name}.csv", index=False)
    pd.DataFrame({"perturbation_id": [new.perturbation_ids[p] for p in added], "label": [new.perturbation_labels[p] for p in added], "type": [new.perturbation_types[p] for p in added],
                  "kept_pairs": [int(kept_per_perturbation[p]) for p in added], "held_out": [p in set(added_held) for p in added]}).to_csv(arguments.output_dir / f"added_perturbations_by_side_{arguments.name}.csv", index=False)
    print(usable[["target", "direction", "target_symptom", "held_out", "held_out_by_node", "family_size", "study_drugs_in_the_statement_direction", "edges_to_nearest_held_out_seed"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
