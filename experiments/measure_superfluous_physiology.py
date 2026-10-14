"""How much of the physiology graph cannot carry a response, how much of it is duplicated, and what the content that
cannot carry a response costs the nodes that can through the in-degree divisor.

The owner, 14 October 2026: "Do you think there's superfluous physiology in the graph?" This script answers it with
counts. "Superfluous" is read three ways, each measured separately:

  1. **Unreachable.** A node that no perturbation seed can reach holds a field of exactly zero in every example, so it
     contributes nothing to any prediction. An edge leaving such a node always carries a message of zero.
  2. **Diluting.** The registered normalisation is item 28's `total_in_degree`: the message into a node is divided by
     the sum of the weights of every edge entering it, over all relations together
     (`linear_response_encoder.signed_stacked_adjacency`, `normalisation == "total_in_degree"`). An edge that always
     carries zero still counts in that divisor, so unreachable content divides down the real signal at every node it
     feeds. This is the sense in which superfluous content is harmful rather than merely inert.
  3. **Duplicated.** Two nodes that receive from exactly the same sources under exactly the same relations hold the
     same field in every example, so the second adds parameters and no information.

The adjacency measured here is the one the encoder builds, not the stored edge table:

  - the derived `depletes_substrate` relation, reaction -> metabolite, one per `substrate_of` edge;
  - the derived `product_opposes` relation, metabolite -> reaction, one per `product_of` edge of a reversible
    reaction (preregistration item 27);
  - edges leaving a currency metabolite are dropped, so a currency node receives and does not send, and its edges
    are in nobody's divisor;
  - a drug's `targets` edges leave the normalised adjacency and are added on their own, so they carry a response and
    are **not** in any divisor (`entry_relation_names`); every other relation is in the divisor.

The readout side. The registered symptom head reads the field at every graph node: the module support is a learned
hard-concrete gate over all N nodes and the pool is a sum over them
(`noisy_or_pathway_module_model.NoisyOrPathwayModuleModel.module_activation`), and the output gate is a sigmoid of
node features that zeroes nothing structurally. So under the registered head there is no readout node set to run a
backward cone from: every node that holds a field is read. The backward cone is reported from two stricter readout
sets instead, as a sensitivity reading, and each is labelled as the counterfactual it is.

Nothing here changes a graph, a label, a model or a model input. It reads two parquet files and prints counts.

Usage:
  PYTHONPATH=. python experiments/measure_superfluous_physiology.py
  PYTHONPATH=. python experiments/measure_superfluous_physiology.py --graph-dir data/processed/graph_full_neuronal
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

SUBSTRATE_RELATION = "substrate_of"
PRODUCT_RELATION = "product_of"
DEPLETES_SUBSTRATE_RELATION = "depletes_substrate"
PRODUCT_OPPOSES_RELATION = "product_opposes"
ENTRY_RELATIONS = ("targets",)  # graph/drug_entry_nodes.py ENTRY_RELATIONS; out of the divisor, still propagating
DRUG_NODE_TYPE = "drug"
PATHOLOGICAL_NODE = "VM_c"  # 405 changes_membrane_potential edges, measured median field 2.9e-15; reported apart
REGISTERED_LABEL_GRADES = ("A", "B")
REGISTERED_RELATION = "induces"
WORST_CASES_REPORTED = 15


def read_graph(graph_directory: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_parquet(graph_directory / "nodes.parquet")
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    return nodes, edges


def effective_edges(nodes: pd.DataFrame, stored_edges: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """The edge table the encoder propagates on: the stored edges, the two derived relations, currency sources dropped.

    `in_divisor` says whether the edge counts in its destination's `total_in_degree`; a drug's entry edges do not.
    Returns the table and the decomposition of how its row count is reached, so the arithmetic can be checked: the
    derived counts are the ones built **before** the currency drop, and the drop is counted over all three blocks
    together. Reading a post-drop derived count beside the drop total double-counts the derived edges that were
    themselves dropped.
    """
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    is_currency_node = nodes.is_currency.fillna(False).to_numpy().astype(bool)
    is_reversible_reaction = nodes.reversible.fillna(False).to_numpy().astype(bool)

    source_position = stored_edges.source_id.map(position_of_node).to_numpy()
    target_position = stored_edges.target_id.map(position_of_node).to_numpy()
    relation = stored_edges.relation_type.to_numpy()
    sign = stored_edges.sign.to_numpy(dtype=np.float64)

    is_substrate_edge = relation == SUBSTRATE_RELATION
    is_product_edge = relation == PRODUCT_RELATION
    product_of_reversible = is_product_edge & is_reversible_reaction[source_position]

    derived_source = np.concatenate([target_position[is_substrate_edge], target_position[product_of_reversible]])
    derived_target = np.concatenate([source_position[is_substrate_edge], source_position[product_of_reversible]])
    derived_relation = np.concatenate([
        np.full(int(is_substrate_edge.sum()), DEPLETES_SUBSTRATE_RELATION, dtype=object),
        np.full(int(product_of_reversible.sum()), PRODUCT_OPPOSES_RELATION, dtype=object)])
    derived_sign = -np.ones(len(derived_source))

    all_source = np.concatenate([source_position, derived_source])
    all_target = np.concatenate([target_position, derived_target])
    all_relation = np.concatenate([relation.astype(object), derived_relation])
    all_sign = np.concatenate([sign, derived_sign])
    is_derived = np.concatenate([np.zeros(len(source_position), dtype=bool), np.ones(len(derived_source), dtype=bool)])

    table = pd.DataFrame({"source": all_source, "target": all_target, "relation": all_relation,
                          "sign": all_sign, "is_derived": is_derived})
    table["in_divisor"] = ~table.relation.isin(ENTRY_RELATIONS)
    dropped_for_currency_source = is_currency_node[table.source.to_numpy()]
    decomposition = {
        "stored_edges": len(stored_edges),
        "derived_depletes_substrate_before_the_drop": int(is_substrate_edge.sum()),
        "derived_product_opposes_before_the_drop": int(product_of_reversible.sum()),
        "edges_before_the_currency_drop": len(table),
        "dropped_for_a_currency_source": int(dropped_for_currency_source.sum()),
    }
    table = table[~dropped_for_currency_source].reset_index(drop=True)
    decomposition["propagated_edges"] = len(table)
    return table, decomposition


def registered_perturbation_seeds(nodes: pd.DataFrame, graph_directory: Path, evidence_directory: Path) -> dict:
    """The node positions each registered perturbation is seeded on, under the registered `nodes` drug-entry mode.

    A gene knockout seeds its gene node; a drug seeds its own drug node and its `targets` edges carry the sign and
    magnitude on (`evaluation/experiment_data.load_experiment_data`). The `targets` mode, which seeds a drug on its
    target proteins instead, is returned beside it because the two give different forward cones.
    """
    position_of_node = {node_id: index for index, node_id in enumerate(nodes.node_id)}
    evidence = pd.read_parquet(evidence_directory / "evidence_records.parquet")
    evidence = evidence[evidence.relation == REGISTERED_RELATION]
    if "positive_report_count" in evidence.columns:
        evidence = evidence[evidence.positive_report_count > 0]
    if "grade" in evidence.columns:
        evidence = evidence[evidence.grade.isin(REGISTERED_LABEL_GRADES)]
    evidence = evidence.drop_duplicates(subset=["perturbation_id"])

    entry_path = graph_directory / "drug_entry_nodes.parquet"
    entry_table = pd.read_parquet(entry_path) if entry_path.exists() else None
    drug_node_of_perturbation = {} if entry_table is None else dict(zip(entry_table.perturbation_id, entry_table.drug_node_id))

    gene_to_protein_path = graph_directory / "gene_to_protein.parquet"
    protein_of_gene: dict[str, list[str]] = {}
    if gene_to_protein_path.exists():
        gene_to_protein = pd.read_parquet(gene_to_protein_path)
        if "drug_target" in gene_to_protein.columns:
            gene_to_protein = gene_to_protein[gene_to_protein.drug_target]
        protein_of_gene = gene_to_protein.groupby("gene_node_id", sort=False).protein_node_id.agg(list).to_dict()

    seeds_under_entry_nodes: set[int] = set()
    seeds_under_targets: set[int] = set()
    gene_perturbations, drug_perturbations, drugs_without_an_entry_node = 0, 0, 0
    for row in evidence.itertuples(index=False):
        named_nodes = [node_id for node_id, _, _ in json.loads(row.perturbation_nodes) if node_id in position_of_node]
        if row.perturbation_type == "drug":
            drug_perturbations += 1
            on_targets = []
            for node_id in named_nodes:
                on_targets += [protein for protein in protein_of_gene.get(node_id, [node_id]) if protein in position_of_node]
            seeds_under_targets.update(position_of_node[node_id] for node_id in on_targets)
            drug_node = drug_node_of_perturbation.get(row.perturbation_id)
            if drug_node in position_of_node:
                seeds_under_entry_nodes.add(position_of_node[drug_node])
            else:
                drugs_without_an_entry_node += 1
                seeds_under_entry_nodes.update(position_of_node[node_id] for node_id in on_targets)
        else:
            gene_perturbations += 1
            positions = {position_of_node[node_id] for node_id in named_nodes}
            seeds_under_entry_nodes.update(positions)
            seeds_under_targets.update(positions)
    return {"under_entry_nodes": np.array(sorted(seeds_under_entry_nodes), dtype=np.int64),
            "under_targets": np.array(sorted(seeds_under_targets), dtype=np.int64),
            "gene_perturbations": gene_perturbations, "drug_perturbations": drug_perturbations,
            "drugs_without_an_entry_node": drugs_without_an_entry_node}


def reachable_from(num_nodes: int, source: np.ndarray, target: np.ndarray, start: np.ndarray) -> np.ndarray:
    """Boolean mask of the nodes a directed walk from `start` can reach, start included."""
    forward = sparse.coo_matrix((np.ones(len(source), dtype=bool), (target, source)), shape=(num_nodes, num_nodes)).tocsr()
    reached = np.zeros(num_nodes, dtype=bool)
    reached[start] = True
    frontier = reached.copy()
    while frontier.any():
        newly = (forward @ frontier.astype(np.int8)).astype(bool) & ~reached
        reached |= newly
        frontier = newly
    return reached


def identical_neighbourhood_groups(table: pd.DataFrame, num_nodes: int, incoming_only: bool) -> dict:
    """Nodes grouped by the multiset of their edges, so a group's members are interchangeable.

    With incoming_only the key is the set of (relation, source) pairs entering the node: two nodes with the same key
    receive the same message from the same places under the same relations, so they hold the same field in every
    example and the readout cannot tell them apart. With incoming_only off the outgoing edges join the key, which is
    the stricter reading of a duplicate node.
    """
    keys: list[list[tuple]] = [[] for _ in range(num_nodes)]
    for source, target, relation in zip(table.source.to_numpy(), table.target.to_numpy(), table.relation.to_numpy()):
        keys[target].append(("in", relation, int(source)))
        if not incoming_only:
            keys[source].append(("out", relation, int(target)))
    grouped: dict[tuple, list[int]] = {}
    for node_position, key in enumerate(keys):
        if not key:
            continue  # a node with no edge at all is counted as unreachable, not as a duplicate
        grouped.setdefault(tuple(sorted(key)), []).append(node_position)
    groups = [members for members in grouped.values() if len(members) > 1]
    return {"groups": groups, "nodes_in_a_group": sum(len(members) for members in groups),
            "redundant_nodes": sum(len(members) - 1 for members in groups),
            "largest_group": max((len(members) for members in groups), default=0)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v4"))
    parser.add_argument("--curated-modules", type=Path, default=Path("docs/curated_pathway_modules.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/label_review/superfluous_physiology.json"))
    arguments = parser.parse_args()

    nodes, stored_edges = read_graph(arguments.graph_dir)
    num_nodes = len(nodes)
    table, decomposition = effective_edges(nodes, stored_edges)
    node_id = nodes.node_id.to_numpy()
    node_type = nodes.node_type.to_numpy()
    display_name = nodes.display_name.to_numpy()
    print(f"graph {arguments.graph_dir}")
    print(f"  stored: {num_nodes:,} nodes, {len(stored_edges):,} edges")
    print(f"  propagated: {len(table):,} edges after the two derived relations and the currency-source drop "
          f"({int((~table.in_divisor).sum()):,} of them outside every divisor, the drug entry edges)")
    print(f"  how that count is reached: {decomposition['stored_edges']:,} stored "
          f"+ {decomposition['derived_depletes_substrate_before_the_drop']:,} derived depletes_substrate "
          f"+ {decomposition['derived_product_opposes_before_the_drop']:,} derived product_opposes "
          f"= {decomposition['edges_before_the_currency_drop']:,}, less "
          f"{decomposition['dropped_for_a_currency_source']:,} whose source is a currency metabolite "
          f"= {decomposition['propagated_edges']:,}")

    seeds = registered_perturbation_seeds(nodes, arguments.graph_dir, arguments.evidence_dir)
    print(f"  registered perturbations: {seeds['gene_perturbations']:,} gene, {seeds['drug_perturbations']:,} drug "
          f"({seeds['drugs_without_an_entry_node']} drugs with no entry-node row, seeded on their targets)")
    print(f"  seed nodes: {len(seeds['under_entry_nodes']):,} under the registered drug-entry mode, "
          f"{len(seeds['under_targets']):,} under --drug-entry targets")

    source = table.source.to_numpy()
    target = table.target.to_numpy()
    forward_cone = reachable_from(num_nodes, source, target, seeds["under_entry_nodes"])
    forward_cone_on_targets = reachable_from(num_nodes, source, target, seeds["under_targets"])

    # 1. Reachability
    print("\n=== 1. Reachability ===")
    print(f"nodes a seed can reach: {int(forward_cone.sum()):,} of {num_nodes:,} ({forward_cone.mean():.3%})")
    print(f"nodes no seed can reach: {int((~forward_cone).sum()):,} ({(~forward_cone).mean():.3%})")
    unreachable_by_type = pd.Series(node_type[~forward_cone]).value_counts()
    for kind, count in unreachable_by_type.items():
        of_kind = int((node_type == kind).sum())
        print(f"    {kind:<20} {count:>6,} of {of_kind:,} ({count / of_kind:.1%})")
    print(f"nodes with no propagated edge at all: {int((~np.isin(np.arange(num_nodes), np.union1d(source, target))).sum()):,}")
    print(f"under --drug-entry targets the forward cone is {int(forward_cone_on_targets.sum()):,} nodes")

    edge_source_reachable = forward_cone[source]
    print(f"edges that can carry a message (source reachable): {int(edge_source_reachable.sum()):,} of {len(table):,} "
          f"({edge_source_reachable.mean():.3%})")
    print(f"edges that can never carry a message: {int((~edge_source_reachable).sum()):,} "
          f"({(~edge_source_reachable).mean():.3%})")
    dead_edges_in_divisor = (~edge_source_reachable) & table.in_divisor.to_numpy()
    print(f"    of those, in a destination's divisor: {int(dead_edges_in_divisor.sum()):,}")

    # What the unreachable content is, and what it costs the nodes that can carry a field.
    edge_into_a_reachable_node = np.zeros(num_nodes, dtype=bool)
    np.logical_or.at(edge_into_a_reachable_node, source, forward_cone[target])
    print(f"unreachable nodes with at least one edge into a reachable node (so they dilute rather than sit idle): "
          f"{int(((~forward_cone) & edge_into_a_reachable_node).sum()):,} of {int((~forward_cone).sum()):,}")
    stored_in_degree = np.zeros(num_nodes, dtype=np.int64)
    np.add.at(stored_in_degree, target, 1)
    for kind in ("gene", "protein"):
        of_kind = (~forward_cone) & (node_type == kind)
        print(f"unreachable {kind} nodes with no incoming propagated edge at all: "
              f"{int((stored_in_degree[of_kind] == 0).sum()):,} of {int(of_kind.sum()):,}")
    print("relations that can enter a gene node (the gene-protein split leaves only the transcriptional one):")
    into_gene = table[node_type[table.target.to_numpy()] == "gene"]
    print("    " + str(into_gene.relation.value_counts().to_dict()))
    unreachable_protein = np.flatnonzero((~forward_cone) & (node_type == "protein"))
    leaving_unreachable_proteins = table[np.isin(table.source.to_numpy(), unreachable_protein)]
    print(f"edges leaving the {len(unreachable_protein):,} unreachable proteins, by relation:")
    print("    " + str(leaving_unreachable_proteins.relation.value_counts().to_dict()))
    unreachable_gene_symbols = set(nodes.gene_symbol.fillna("").to_numpy()[(~forward_cone) & (node_type == "gene")]) - {""}
    every_grade = pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet")
    print(f"unreachable gene symbols that appear in the evidence table at some grade or relation the registered "
          f"filter drops: {len(unreachable_gene_symbols & set(every_grade.perturbation_id.astype(str))):,} of "
          f"{len(unreachable_gene_symbols):,}")

    # The backward cone. Under the registered head every node is read, so two stricter readouts are reported instead.
    readouts: dict[str, np.ndarray] = {}
    readouts["every node (the registered symptom head)"] = np.ones(num_nodes, dtype=bool)
    laboratory_readout = (node_type == "metabolite")
    readouts["every metabolite (the laboratory readout's node kind)"] = laboratory_readout
    if arguments.curated_modules.exists():
        curated = pd.read_csv(arguments.curated_modules)
        curated_genes = {gene for genes in curated.genes.fillna("") for gene in str(genes).split(";") if gene}
        gene_symbol = nodes.gene_symbol.fillna("").to_numpy()
        curated_module_nodes = np.isin(gene_symbol, sorted(curated_genes)) & np.isin(node_type, ["gene", "protein"])
        readouts[f"the {len(curated)} curated pathway modules' genes and proteins"] = curated_module_nodes
    backward_cones = {}
    print("\nbackward cones, by readout set:")
    for name, readout_mask in readouts.items():
        backward = reachable_from(num_nodes, target, source, np.flatnonzero(readout_mask))
        backward_cones[name] = backward
        in_neither = (~forward_cone) & (~backward)
        edges_in_neither = (~forward_cone[source]) & (~backward[target])
        print(f"  {name}: readout {int(readout_mask.sum()):,} nodes; backward cone {int(backward.sum()):,} nodes "
              f"({backward.mean():.1%}); in neither cone {int(in_neither.sum()):,} nodes and "
              f"{int(edges_in_neither.sum()):,} edges")
        on_a_route = forward_cone & backward
        print(f"      on a seed-to-readout route: {int(on_a_route.sum()):,} nodes ({on_a_route.mean():.1%})")

    # 2. Dilution
    print("\n=== 2. Dilution of the total_in_degree divisor ===")
    in_divisor = table.in_divisor.to_numpy()
    divisor_total = np.zeros(num_nodes, dtype=np.int64)
    np.add.at(divisor_total, target[in_divisor], 1)
    divisor_dead = np.zeros(num_nodes, dtype=np.int64)
    np.add.at(divisor_dead, target[in_divisor & ~edge_source_reachable], 1)
    has_divisor = divisor_total > 0
    on_route = forward_cone & has_divisor
    on_route_without_pathological = on_route & (node_id != PATHOLOGICAL_NODE)
    dead_share = np.zeros(num_nodes)
    dead_share[has_divisor] = divisor_dead[has_divisor] / divisor_total[has_divisor]

    measured = np.flatnonzero(on_route_without_pathological)
    print(f"nodes with a divisor and a field (VM_c excluded): {len(measured):,}")
    print(f"total divisor mass: {int(divisor_total[measured].sum()):,} edges, of which "
          f"{int(divisor_dead[measured].sum()):,} can never carry anything "
          f"({divisor_dead[measured].sum() / max(divisor_total[measured].sum(), 1):.3%})")
    print(f"dead share of the divisor: median {np.median(dead_share[measured]):.4f}, "
          f"mean {dead_share[measured].mean():.4f}, "
          f"90th percentile {np.percentile(dead_share[measured], 90):.4f}, "
          f"99th percentile {np.percentile(dead_share[measured], 99):.4f}")
    print(f"nodes whose divisor is at least half dead: {int((dead_share[measured] >= 0.5).sum()):,}; "
          f"entirely dead: {int((dead_share[measured] >= 1.0).sum()):,}")
    entirely_dead_divisor = measured[divisor_dead[measured] >= divisor_total[measured]]
    print("    the entirely dead ones, named (each carries a field only through its own seed input or an "
          "out-of-divisor entry edge):")
    for position in entirely_dead_divisor:
        print(f"      {node_id[position]:<22} {node_type[position]:<10} divisor {divisor_total[position]}, "
              f"dead {divisor_dead[position]}")
    print(f"divisor over all relations (the factor a single live edge is divided by): "
          f"median {np.median(divisor_total[measured]):.0f}, "
          f"90th percentile {np.percentile(divisor_total[measured], 90):.0f}, "
          f"largest {divisor_total[measured].max():,}")

    print("\nworst by dead edges in the divisor (VM_c excluded):")
    by_dead_count = measured[np.argsort(-divisor_dead[measured])][:WORST_CASES_REPORTED]
    for position in by_dead_count:
        print(f"    {node_id[position]:<34} {display_name[position][:40]:<42} {node_type[position]:<15} "
              f"dead {divisor_dead[position]:>6,} of {divisor_total[position]:>6,} ({dead_share[position]:.1%})")
    print("\nworst by dead share among nodes with a divisor of at least 20 (VM_c excluded):")
    sizeable = measured[divisor_total[measured] >= 20]
    by_share = sizeable[np.lexsort((-divisor_total[sizeable], -dead_share[sizeable]))][:WORST_CASES_REPORTED]
    for position in by_share:
        print(f"    {node_id[position]:<34} {display_name[position][:40]:<42} {node_type[position]:<15} "
              f"dead {divisor_dead[position]:>6,} of {divisor_total[position]:>6,} ({dead_share[position]:.1%})")
    print("\nlargest divisors, dead or not (VM_c excluded):")
    by_divisor = measured[np.argsort(-divisor_total[measured])][:WORST_CASES_REPORTED]
    for position in by_divisor:
        print(f"    {node_id[position]:<34} {display_name[position][:40]:<42} {node_type[position]:<15} "
              f"divisor {divisor_total[position]:>6,}, dead {divisor_dead[position]:>6,} ({dead_share[position]:.1%})")

    live_divisor = divisor_total - divisor_dead
    gain_from_dropping_dead = np.ones(num_nodes)
    with np.errstate(divide="ignore", invalid="ignore"):
        gain_from_dropping_dead[live_divisor > 0] = (divisor_total[live_divisor > 0] / live_divisor[live_divisor > 0])
    with_live = measured[live_divisor[measured] > 0]
    print(f"\nthe factor a live message would gain if the dead edges left the divisor: "
          f"median {np.median(gain_from_dropping_dead[with_live]):.3f}, "
          f"90th percentile {np.percentile(gain_from_dropping_dead[with_live], 90):.3f}, "
          f"99th percentile {np.percentile(gain_from_dropping_dead[with_live], 99):.3f}, "
          f"largest {gain_from_dropping_dead[with_live].max():.3f}")

    # Where the dilution actually lands. A dead edge into a node that itself holds no field dilutes nothing, so the
    # destinations are split: only the ones that DO carry a field are harmed. Counting the others as diluted would
    # name the harmless case (dead content feeding dead content) as the harmful one.
    print("\nwhat the dead edges of the two worst relations dilute, split by whether the destination holds a field:")
    dilution_by_relation = {}
    for relation_name in ("catalyzed_by", "encodes"):
        dead_rows = table[(table.relation == relation_name) & (~edge_source_reachable) & in_divisor]
        if not len(dead_rows):
            continue
        destinations = np.unique(dead_rows.target.to_numpy())
        destinations_with_a_field = destinations[forward_cone[destinations]]
        destinations_without_a_field = destinations[~forward_cone[destinations]]
        if not len(destinations_with_a_field):
            continue
        hit_share = divisor_dead[destinations_with_a_field] / divisor_total[destinations_with_a_field]
        hit_gain = divisor_total[destinations_with_a_field] / np.maximum(
            divisor_total[destinations_with_a_field] - divisor_dead[destinations_with_a_field], 1)
        kinds = pd.Series(node_type[destinations_with_a_field]).value_counts().to_dict()
        print(f"  {relation_name}: {len(dead_rows):,} dead edges reach {len(destinations):,} destinations")
        print(f"      {len(destinations_with_a_field):,} of them hold a field and are diluted {kinds}")
        print(f"      {len(destinations_without_a_field):,} hold no field, so these dilute nothing "
              f"{pd.Series(node_type[destinations_without_a_field]).value_counts().to_dict()}")
        print(f"      dead share of the diluted nodes' divisor: median {np.median(hit_share):.3f}, "
              f"90th percentile {np.percentile(hit_share, 90):.3f}")
        print(f"      the factor a live message there would gain if they left: median {np.median(hit_gain):.3f}, "
              f"90th percentile {np.percentile(hit_gain, 90):.3f}, largest {hit_gain.max():.3f}")
        dilution_by_relation[relation_name] = {
            "dead_edges": len(dead_rows),
            "destinations": len(destinations),
            "destinations_holding_a_field": len(destinations_with_a_field),
            "destinations_holding_no_field": len(destinations_without_a_field),
            "destination_kinds_holding_a_field": {str(kind): int(count) for kind, count in kinds.items()},
            "dead_share_median": float(np.median(hit_share)),
            "dead_share_p90": float(np.percentile(hit_share, 90)),
            "gain_median": float(np.median(hit_gain)),
            "gain_p90": float(np.percentile(hit_gain, 90)),
            "gain_largest": float(hit_gain.max()),
        }

    print("\nwhere the dead divisor mass comes from, by the node type of the edge's source:")
    dead_in_divisor_rows = table[(~edge_source_reachable) & in_divisor]
    for kind, count in pd.Series(node_type[dead_in_divisor_rows.source.to_numpy()]).value_counts().items():
        print(f"    {kind:<20} {count:>7,}")
    print("and by relation:")
    for relation, count in dead_in_divisor_rows.relation.value_counts().items():
        print(f"    {relation:<30} {count:>7,}")

    is_currency_node = nodes.is_currency.fillna(False).to_numpy().astype(bool)
    print(f"\ndivisor mass at the {int(is_currency_node.sum())} currency metabolites: "
          f"{int(divisor_total[is_currency_node].sum()):,} of {int(divisor_total.sum()):,} "
          f"({divisor_total[is_currency_node].sum() / divisor_total.sum():.1%}); "
          f"their edges send nothing on, and the field at them is still read by the pooled head")

    pathological = np.flatnonzero(node_id == PATHOLOGICAL_NODE)
    for position in pathological:
        by_relation = table[(table.target.to_numpy() == position) & in_divisor].relation.value_counts().to_dict()
        print(f"\n{PATHOLOGICAL_NODE}, reported apart: divisor {divisor_total[position]:,}, "
              f"dead {divisor_dead[position]:,} ({dead_share[position]:.1%}), reachable {bool(forward_cone[position])}; "
              f"by relation {by_relation}")

    # 3. Duplication
    print("\n=== 3. Duplication ===")
    is_metabolite = node_type == "metabolite"
    base = nodes.base_metabolite_id.fillna("").to_numpy()
    compartment_groups = pd.Series(base[is_metabolite & (base != "")]).value_counts()
    multi = compartment_groups[compartment_groups > 1]
    print(f"metabolite species with more than one compartment copy: {len(multi):,} of {len(compartment_groups):,} species; "
          f"{int(multi.sum()):,} nodes, {int((multi - 1).sum()):,} of them copies beyond the first")
    print("    largest: " + ", ".join(f"{name} x{count}" for name, count in multi.head(6).items()))
    print(f"    copies per species: median {int(compartment_groups.median())}, largest {int(compartment_groups.max())}")
    propagated_degree = np.zeros(num_nodes, dtype=np.int64)
    np.add.at(propagated_degree, source, 1)
    np.add.at(propagated_degree, target, 1)
    copy_node_positions = np.flatnonzero(is_metabolite & np.isin(base, multi.index.to_numpy()))
    print(f"    of the {len(copy_node_positions):,} copy nodes, with a total propagated degree of 2 or less: "
          f"{int((propagated_degree[copy_node_positions] <= 2).sum()):,}")
    subsystem_of_node = nodes.subsystem.fillna("").to_numpy()
    is_reaction_node = node_type == "reaction"
    metabolites_touching_only_exchange = []
    for position in np.flatnonzero(is_metabolite):
        partners = np.concatenate([target[source == position], source[target == position]])
        reaction_partners = partners[is_reaction_node[partners]]
        if len(reaction_partners) and len(reaction_partners) == len(partners) \
                and set(subsystem_of_node[reaction_partners]) == {"Exchange/demand reactions"}:
            metabolites_touching_only_exchange.append(position)
    metabolites_touching_only_exchange = np.array(metabolites_touching_only_exchange, dtype=np.int64)
    print(f"    metabolite nodes whose only partners are Exchange/demand reactions: "
          f"{len(metabolites_touching_only_exchange)}; reachable "
          f"{int(forward_cone[metabolites_touching_only_exchange].sum())}; compartment suffixes "
          f"{sorted({str(identifier)[-1] for identifier in node_id[metabolites_touching_only_exchange]})}")

    pool_prefixes = ("NEURON_", "REDOX_")
    pool_nodes = np.array([any(str(identifier).startswith(prefix) for prefix in pool_prefixes) for identifier in node_id])
    print(f"\npool nodes ({int(pool_nodes.sum())}): " + ", ".join(node_id[pool_nodes]))
    for position in np.flatnonzero(pool_nodes):
        species = base[position]
        twins = [node_id[other] for other in np.flatnonzero(is_metabolite & (base == species)) if other != position] if species else []
        print(f"    {node_id[position]:<20} base {species or '(none)':<12} divisor {divisor_total[position]:>5,} "
              f"reachable {bool(forward_cone[position])!s:<5} plain nodes of the same species: {twins}")

    drug_positions = np.flatnonzero(node_type == DRUG_NODE_TYPE)
    targets_of_drug = {}
    entry_rows = table[table.relation.isin(ENTRY_RELATIONS)]
    for drug, target_position, sign in zip(entry_rows.source.to_numpy(), entry_rows.target.to_numpy(), entry_rows.sign.to_numpy()):
        targets_of_drug.setdefault(int(drug), set()).add((int(target_position), float(sign)))
    drug_groups: dict[tuple, list[int]] = {}
    for drug in drug_positions:
        drug_groups.setdefault(tuple(sorted(targets_of_drug.get(int(drug), ()))), []).append(int(drug))
    shared = [members for members in drug_groups.values() if len(members) > 1]
    print(f"\ndrugs sharing an identical target set with sign: {sum(len(members) for members in shared)} of "
          f"{len(drug_positions)} drugs in {len(shared)} groups, largest {max((len(members) for members in shared), default=0)}")
    for members in sorted(shared, key=len, reverse=True)[:4]:
        print(f"    size {len(members):>3}: " + ", ".join(str(display_name[position]) for position in members[:6])
              + (" ..." if len(members) > 6 else ""))

    for incoming_only, label in ((True, "identical incoming neighbourhood (field-identical)"),
                                 (False, "identical incoming and outgoing neighbourhood")):
        result = identical_neighbourhood_groups(table, num_nodes, incoming_only=incoming_only)
        groups = result["groups"]
        reachable_groups = [members for members in groups if forward_cone[members].any()]
        print(f"\n{label}: {len(groups):,} groups, {result['nodes_in_a_group']:,} nodes, "
              f"{result['redundant_nodes']:,} redundant, largest group {result['largest_group']:,}")
        print(f"    groups with any reachable member: {len(reachable_groups):,} "
              f"({sum(len(members) - 1 for members in reachable_groups):,} redundant nodes)")
        largest = sorted(groups, key=len, reverse=True)[:5]
        for members in largest:
            kinds = pd.Series(node_type[members]).value_counts().to_dict()
            incoming_to_the_first = table[table.target.to_numpy() == members[0]]
            feeders = sorted({f"{node_id[int(feeding)]} [{relation}]" for feeding, relation
                              in zip(incoming_to_the_first.source.to_numpy(), incoming_to_the_first.relation.to_numpy())})
            print(f"      size {len(members):>4}: {kinds}; fed by {feeders if len(feeders) <= 3 else feeders[:3] + ['...']}"
                  f"; e.g. " + ", ".join(f"{node_id[position]} ({str(display_name[position])[:24]})" for position in members[:3]))
        if incoming_only:
            base_metabolite = nodes.base_metabolite_id.fillna("").to_numpy()
            one_species_groups = [members for members in groups
                                  if len({base_metabolite[position] for position in members}) == 1
                                  and base_metabolite[members[0]] != ""
                                  and all(node_type[position] == "metabolite" for position in members)]
            print(f"    field-identical groups that are compartment copies of ONE metabolite species: "
                  f"{len(one_species_groups)}")
            for members in one_species_groups:
                print(f"      base {base_metabolite[members[0]]}: " + ", ".join(node_id[position] for position in members)
                      + f"  ({str(display_name[members[0]])[:60]})")

    # 4. Dead relations
    print("\n=== 4. Relations, by what their edges can carry ===")
    relation_summary = []
    for relation, group in table.groupby("relation"):
        group_source = group.source.to_numpy()
        group_target = group.target.to_numpy()
        live = int(forward_cone[group_source].sum())
        into_unreachable = int((~forward_cone[group_target]).sum())
        into_currency = int(nodes.is_currency.fillna(False).to_numpy().astype(bool)[group_target].sum())
        into_entry_node = int((node_type[group_target] == DRUG_NODE_TYPE).sum())
        signed = int((group.sign.to_numpy() != 0).sum())
        relation_summary.append({"relation": relation, "edges": len(group), "live_edges": live,
                                 "edges_into_an_unreachable_node": into_unreachable,
                                 "edges_into_a_currency_node": into_currency,
                                 "edges_into_an_entry_node": into_entry_node,
                                 "in_divisor": bool(group.in_divisor.all()), "signed_edges": signed})
    relation_summary.sort(key=lambda row: -row["edges"])
    print(f"{'relation':<30} {'edges':>8} {'live':>8} {'live %':>8} {'into dead':>10} {'into currency':>14} {'into drug':>10}")
    for row in relation_summary:
        print(f"{row['relation']:<30} {row['edges']:>8,} {row['live_edges']:>8,} "
              f"{row['live_edges'] / row['edges']:>7.1%} {row['edges_into_an_unreachable_node']:>10,} "
              f"{row['edges_into_a_currency_node']:>14,} {row['edges_into_an_entry_node']:>10,}")
    dead_relations = [row["relation"] for row in relation_summary if row["live_edges"] == 0]
    print(f"relations with no live edge at all: {dead_relations or 'none'}")
    # sequesters is live by the structural reading and inert by the per-perturbation one: an entry node that is not
    # the seeded drug is held at zero, so a carrier -> drug edge can carry something in one perturbation only.
    for relation_name in ("sequesters", "clears"):
        carrier_rows = table[table.relation == relation_name]
        if not len(carrier_rows):
            continue
        carriers = sorted(set(node_id[carrier_rows.source.to_numpy()]))
        destinations = np.unique(carrier_rows.target.to_numpy())
        print(f"{relation_name}: {len(carrier_rows)} edges from {len(carriers)} carriers {carriers} into "
              f"{len(destinations)} nodes "
              f"{pd.Series(node_type[destinations]).value_counts().to_dict()}")
    small = [row for row in relation_summary if row["edges"] <= 100]
    print("relations with 100 edges or fewer: " + ", ".join(f"{row['relation']} ({row['edges']})" for row in small))

    # 5. Metabolic subsystems
    print("\n=== 5. Human-GEM subsystems ===")
    is_reaction = node_type == "reaction"
    subsystem = nodes.subsystem.fillna("").to_numpy()
    touching_edges = np.zeros(num_nodes, dtype=np.int64)
    np.add.at(touching_edges, source, 1)
    np.add.at(touching_edges, target, 1)
    rows = []
    for name in sorted(set(subsystem[is_reaction]) - {""}):
        members = np.flatnonzero(is_reaction & (subsystem == name))
        rows.append({"subsystem": name, "reactions": len(members),
                     "reactions_reachable": int(forward_cone[members].sum()),
                     "edges_touching": int(touching_edges[members].sum()),
                     "divisor_mass": int(divisor_total[members].sum()),
                     "dead_divisor_mass": int(divisor_dead[members].sum())})
    strict_readout_name = next((name for name in backward_cones if name.startswith("the ")), None)
    strict_backward = backward_cones[strict_readout_name] if strict_readout_name else np.ones(num_nodes, dtype=bool)
    for row, name in zip(rows, sorted(set(subsystem[is_reaction]) - {""})):
        members = np.flatnonzero(is_reaction & (subsystem == name))
        row["reactions_on_a_route_to_the_strict_readout"] = int((forward_cone[members] & strict_backward[members]).sum())
    rows.sort(key=lambda row: (row["reactions_reachable"] / row["reactions"], -row["edges_touching"]))
    unreached = [row for row in rows if row["reactions_reachable"] == 0]
    print(f"subsystems with {len(rows)} named: {len(unreached)} have no reaction any seed can reach")
    for row in unreached:
        print(f"    {row['subsystem']:<55} {row['reactions']:>6,} reactions, {row['edges_touching']:>8,} edges")
    print("\nsubsystems with the lowest reachable share:")
    for row in rows[:WORST_CASES_REPORTED]:
        print(f"    {row['subsystem']:<55} {row['reactions_reachable']:>6,} of {row['reactions']:>6,} reachable, "
              f"{row['edges_touching']:>8,} edges")
    off_route = [row for row in rows if row["reactions_on_a_route_to_the_strict_readout"] == 0]
    print(f"\nsubsystems with no reaction on a route to {strict_readout_name}: {len(off_route)}")
    for row in sorted(off_route, key=lambda item: -item["edges_touching"])[:WORST_CASES_REPORTED]:
        print(f"    {row['subsystem']:<55} {row['reactions']:>6,} reactions, {row['edges_touching']:>8,} edges, "
              f"divisor mass {row['divisor_mass']:>7,}")
    print("\nthe four subsystems that are modelling devices rather than physiology:")
    for name in ("Transport reactions", "Exchange/demand reactions", "Isolated", "Pool reactions"):
        picked = next((row for row in rows if row["subsystem"] == name), None)
        if picked is None:
            continue
        not_on_a_route = picked["reactions"] - picked["reactions_on_a_route_to_the_strict_readout"]
        print(f"    {name:<30} {picked['reactions']:>6,} reactions, {picked['reactions_reachable']:>6,} reachable, "
              f"{picked['reactions_on_a_route_to_the_strict_readout']:>6,} on a strict route "
              f"({not_on_a_route:>5,} not), {picked['edges_touching']:>8,} edges, "
              f"divisor mass {picked['divisor_mass']:>7,}")

    print("\nsubsystems by edges touching them:")
    for row in sorted(rows, key=lambda item: -item["edges_touching"])[:WORST_CASES_REPORTED]:
        print(f"    {row['subsystem']:<55} {row['reactions_reachable']:>6,} of {row['reactions']:>6,} reachable, "
              f"{row['edges_touching']:>8,} edges")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps({
        "graph_dir": str(arguments.graph_dir),
        "stored_nodes": int(num_nodes), "stored_edges": len(stored_edges),
        "propagated_edges": len(table),
        "propagated_edge_decomposition": decomposition,
        "dilution_by_relation": dilution_by_relation,
        "entirely_dead_divisor_nodes": [str(node_id[position]) for position in entirely_dead_divisor],
        "seed_nodes_under_entry_nodes": len(seeds["under_entry_nodes"]),
        "seed_nodes_under_targets": len(seeds["under_targets"]),
        "gene_perturbations": seeds["gene_perturbations"], "drug_perturbations": seeds["drug_perturbations"],
        "nodes_reachable": int(forward_cone.sum()), "nodes_unreachable": int((~forward_cone).sum()),
        "unreachable_by_node_type": {str(kind): int(count) for kind, count in unreachable_by_type.items()},
        "live_edges": int(edge_source_reachable.sum()), "dead_edges": int((~edge_source_reachable).sum()),
        "dead_edges_in_a_divisor": int(dead_edges_in_divisor.sum()),
        "divisor_mass_total": int(divisor_total[measured].sum()),
        "divisor_mass_dead": int(divisor_dead[measured].sum()),
        "dead_share_median": float(np.median(dead_share[measured])),
        "dead_share_mean": float(dead_share[measured].mean()),
        "dead_share_p90": float(np.percentile(dead_share[measured], 90)),
        "nodes_with_a_half_dead_divisor": int((dead_share[measured] >= 0.5).sum()),
        "vm_c": {"divisor": int(divisor_total[pathological[0]]), "dead": int(divisor_dead[pathological[0]])} if len(pathological) else None,
        "compartment_copy_groups": len(multi), "compartment_copy_nodes": int(multi.sum()),
        "relations": relation_summary,
        "subsystems": rows,
    }, indent=1, sort_keys=True))
    print(f"\nwrote {arguments.output}")


if __name__ == "__main__":
    main()
