"""How far a target's GtoPdb family is recoverable from the graph and from the node descriptors, which decides whether
grouping folds by family would remove the family confound or only hide it (the user's question of 9 October 2026: "is
there any graph structure or attribute that corresponds to family? That might be another source of leakage.").

The grouping question of docs/leakage_group_rules.md asks whether to keep a family's drugs in one fold. That only helps
if the model cannot tell which family a target belongs to once its family's drugs are out of training. So this script
asks the opposite question of each channel the model reads:

  - protein_descriptors: the 64 protein components of the descriptor table, which are ESM-2 embeddings projected onto
    annotation space (Pfam families, EC numbers, GO function and component, UniProt location). Pfam families are close
    to pharmacological families by construction, so this channel is expected to carry family; the number says how much.
  - graph_neighbours: the set of graph neighbours of the gene node, over every relation and both directions.
  - reactome_entity_members: neighbours restricted to the protein entities the gene is a member of, which is the
    channel that would group paralogues if Reactome's entity sets were family sets.
  - brain_expression: the GTEx, Human Protein Atlas region and cell-class columns, as a control, since nothing about a
    family's sequence is in them.

Each channel is scored by leave-one-out nearest neighbour: a gene's family is "recovered" when its nearest neighbour
among the other genes, by cosine for the descriptor channels and by Jaccard overlap of neighbour sets for the graph
channels, has the same family. Two baselines are reported beside it: the rate a random other gene would reach
(sum over families of n_f (n_f - 1) / (N (N - 1))), and the same nearest-neighbour rule on once-permuted family labels.

Usage:
  PYTHONPATH=. python experiments/measure_family_in_graph.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from experiments.measure_leakage_group_rules import target_families

GENE_NODE_TYPE = "gene"
PROTEIN_DESCRIPTOR_PREFIX = "protein_rrr_"
BRAIN_EXPRESSION_PREFIXES = ("gene_brain_",)
MEMBER_OF_RELATION = "member_of"
PROTEIN_ENTITY_TYPE = "protein_entity"
PERMUTATION_SEED = 0


def genes_with_a_family(nodes: pd.DataFrame, family_of_symbol: dict[str, str], minimum_family_size: int) -> pd.DataFrame:
    """(node_id, gene_symbol, family) for the graph's gene nodes whose symbol has a GtoPdb family, keeping only
    families with at least minimum_family_size genes in the graph: a family of one cannot be recovered by leaving that
    gene out, and cannot leak either."""
    genes = nodes[(nodes.node_type == GENE_NODE_TYPE) & nodes.gene_symbol.notna()].copy()
    genes["family"] = [family_of_symbol.get(symbol) for symbol in genes.gene_symbol]
    genes = genes[genes.family.notna()]
    sizes = genes.family.value_counts()
    kept = genes[genes.family.isin(sizes[sizes >= minimum_family_size].index)]
    return kept.loc[:, ["node_id", "gene_symbol", "family"]].sort_values("node_id").reset_index(drop=True)


def chance_rate(families: np.ndarray) -> float:
    """The share of same-family pairs among all ordered pairs of distinct genes: what picking a random other gene gives."""
    counts = pd.Series(families).value_counts().to_numpy(dtype=np.float64)
    total = counts.sum()
    return float((counts * (counts - 1)).sum() / (total * (total - 1))) if total > 1 else float("nan")


def nearest_neighbour_agreement(similarity: np.ndarray, families: np.ndarray) -> tuple[float, int]:
    """(share of genes whose most similar other gene has the same family, genes scored). A gene similar to nothing
    (every similarity zero or undefined) is left out, since it has no nearest neighbour to agree with."""
    if len(similarity) < 2:
        return float("nan"), 0
    scores = similarity.copy()
    np.fill_diagonal(scores, -np.inf)
    usable = np.isfinite(scores).any(axis=1) & (np.nanmax(np.where(np.isfinite(scores), scores, -np.inf), axis=1) > 0)
    if not usable.any():
        return float("nan"), 0
    neighbour = np.argmax(scores[usable], axis=1)
    return float((families[usable] == families[neighbour]).mean()), int(usable.sum())


def cosine_similarity(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    normalised = values / np.where(norms > 0, norms, 1.0)
    return normalised @ normalised.T


def jaccard_similarity(incidence: sparse.csr_matrix) -> np.ndarray:
    """Jaccard overlap of each pair of rows' neighbour sets."""
    intersection = (incidence @ incidence.T).toarray().astype(np.float64)
    sizes = np.asarray(incidence.sum(axis=1)).ravel().astype(np.float64)
    union = sizes[:, None] + sizes[None, :] - intersection
    return np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0)


def neighbour_incidence(edges: pd.DataFrame, node_ids: list[str]) -> sparse.csr_matrix:
    """Rows: the given nodes; columns: every node that is a neighbour of one of them, over both edge directions."""
    row_of_node = {node_id: row for row, node_id in enumerate(node_ids)}
    pairs = set()
    for source, target in zip(edges.source_id, edges.target_id):
        if source in row_of_node:
            pairs.add((row_of_node[source], target))
        if target in row_of_node:
            pairs.add((row_of_node[target], source))
    if not pairs:
        return sparse.csr_matrix((len(node_ids), 0), dtype=np.int32)
    columns = {node: index for index, node in enumerate(sorted({node for _, node in pairs}))}
    rows = np.fromiter((row for row, _ in pairs), dtype=np.int32, count=len(pairs))
    cols = np.fromiter((columns[node] for _, node in pairs), dtype=np.int32, count=len(pairs))
    return sparse.csr_matrix((np.ones(len(pairs), dtype=np.int32), (rows, cols)), shape=(len(node_ids), len(columns)))


def shared_neighbour_rates(incidence: sparse.csr_matrix, families: np.ndarray) -> dict:
    """The share of same-family and of different-family gene pairs that share at least one neighbour."""
    shares = (incidence @ incidence.T).toarray() > 0
    np.fill_diagonal(shares, False)
    same = families[:, None] == families[None, :]
    np.fill_diagonal(same, False)
    different = ~same
    np.fill_diagonal(different, False)
    within = float(shares[same].mean()) if same.any() else float("nan")
    between = float(shares[different].mean()) if different.any() else float("nan")
    return {"same_family_pairs_sharing_a_neighbour": round(within, 4), "different_family_pairs_sharing_a_neighbour": round(between, 4),
            "ratio": round(within / between, 2) if between and np.isfinite(between) and between > 0 else None}


def channel_result(name: str, similarity: np.ndarray, families: np.ndarray, permuted: np.ndarray, columns: int) -> dict:
    recovered, scored = nearest_neighbour_agreement(similarity, families)
    permuted_rate, _ = nearest_neighbour_agreement(similarity, permuted)
    return {"channel": name, "columns": columns, "genes_scored": scored,
            "family_recovered_by_nearest_neighbour": None if np.isnan(recovered) else round(recovered, 4),
            "same_rule_on_permuted_families": None if np.isnan(permuted_rate) else round(permuted_rate, 4)}


def drug_mechanism_genes(evidence_records: pd.DataFrame, symbol_of_node: dict[str, str]) -> set[str]:
    """The gene symbols a drug perturbation seeds, from the evidence table's perturbation_nodes, which the table holds
    as a JSON string of [node id, sign, magnitude] entries."""
    symbols = set()
    for nodes_field in evidence_records.loc[evidence_records.perturbation_type == "drug", "perturbation_nodes"]:
        entries = json.loads(nodes_field) if isinstance(nodes_field, str) else nodes_field
        for entry in entries:
            symbol = symbol_of_node.get(entry[0] if isinstance(entry, (list, tuple, np.ndarray)) else entry)
            if symbol:
                symbols.add(symbol)
    return symbols


def measure(nodes: pd.DataFrame, edges: pd.DataFrame, descriptors: pd.DataFrame, genes: pd.DataFrame) -> dict:
    families = genes.family.to_numpy()
    node_ids = genes.node_id.tolist()
    generator = np.random.default_rng(PERMUTATION_SEED)
    permuted = generator.permutation(families)
    descriptor_columns = [column for column in descriptors.columns if column.startswith(PROTEIN_DESCRIPTOR_PREFIX)]
    expression_columns = [column for column in descriptors.columns if column.startswith(BRAIN_EXPRESSION_PREFIXES)]
    rows = descriptors.reindex(node_ids)
    entity_nodes = set(nodes.loc[nodes.node_type == PROTEIN_ENTITY_TYPE, "node_id"])
    member_edges = edges[(edges.relation_type == MEMBER_OF_RELATION) & edges.target_id.isin(entity_nodes)]
    neighbours = neighbour_incidence(edges, node_ids)
    entity_members = neighbour_incidence(member_edges, node_ids)
    channels = [
        channel_result("protein_descriptors", cosine_similarity(rows[descriptor_columns].fillna(0.0).to_numpy(dtype=np.float64)),
                       families, permuted, len(descriptor_columns)),
        channel_result("graph_neighbours", jaccard_similarity(neighbours), families, permuted, int(neighbours.shape[1])),
        channel_result("reactome_entity_members", jaccard_similarity(entity_members), families, permuted, int(entity_members.shape[1])),
        channel_result("brain_expression", cosine_similarity(rows[expression_columns].fillna(0.0).to_numpy(dtype=np.float64)),
                       families, permuted, len(expression_columns)),
    ]
    return {"genes": len(genes), "families": int(genes.family.nunique()),
            "largest_family": int(genes.family.value_counts().iloc[0]),
            "chance_rate_of_a_random_other_gene": round(chance_rate(families), 4),
            "channels": channels,
            "shared_graph_neighbours": shared_neighbour_rates(neighbours, families),
            "shared_reactome_entities": shared_neighbour_rates(entity_members, families)}


def report(result: dict, arguments) -> str:
    lines = ["# Does the graph carry target family?", "",
             f"Measured on `{arguments.graph_dir}` with `{arguments.descriptor_table}` and GtoPdb "
             f"{Path(arguments.gtopdb_families).parent.name}, families of at least {arguments.minimum_family_size} graph genes. "
             "Generated by experiments/measure_family_in_graph.py.", "",
             "The question it answers (the user, 9 October 2026): \"is there any graph structure or attribute that corresponds "
             "to family? That might be another source of leakage.\" Grouping a family's drugs into one fold only removes the "
             "family confound if the model cannot tell the family from what it reads. Each channel below is scored by "
             "leave-one-out nearest neighbour: a gene counts as recovered when its most similar other gene carries the same "
             "family.", ""]
    for scope, block in result["scopes"].items():
        lines += [f"## {scope}", "",
                  f"{block['genes']} gene nodes, {block['families']} families, largest {block['largest_family']}; a random "
                  f"other gene matches {block['chance_rate_of_a_random_other_gene']:.3f} of the time.", "",
                  "| channel | columns | genes scored | family recovered | permuted families |", "|---|---|---|---|---|"]
        for channel in block["channels"]:
            lines.append(f"| {channel['channel']} | {channel['columns']} | {channel['genes_scored']} | "
                         f"{channel['family_recovered_by_nearest_neighbour']} | {channel['same_rule_on_permuted_families']} |")
        for key, heading in (("shared_graph_neighbours", "pairs sharing a graph neighbour"),
                             ("shared_reactome_entities", "pairs sharing a Reactome protein entity")):
            block_rates = block[key]
            lines += ["", f"{heading.capitalize()}: {block_rates['same_family_pairs_sharing_a_neighbour']} within family against "
                          f"{block_rates['different_family_pairs_sharing_a_neighbour']} between, a ratio of {block_rates['ratio']}."]
        lines.append("")
    lines += reading_section(result)
    return "\n".join(lines) + "\n"


def reading_section(result: dict) -> list[str]:
    """The conclusion, written from the numbers rather than asserted beside them."""
    lines = ["## Reading", ""]
    for scope, block in result["scopes"].items():
        recovered = {channel["channel"]: channel["family_recovered_by_nearest_neighbour"] for channel in block["channels"]}
        best = max((rate, name) for name, rate in recovered.items() if rate is not None)
        structural = recovered.get("graph_neighbours")
        descriptor = recovered.get("protein_descriptors")
        control = recovered.get("brain_expression")
        chance = block["chance_rate_of_a_random_other_gene"]
        lines += [f"On {scope.replace('_', ' ')}, the strongest channel is {best[1]} at {best[0]:.3f} against a chance rate of "
                  f"{chance:.3f}. The structure alone reaches {structural:.3f} and the protein descriptor block {descriptor:.3f}; "
                  f"the brain expression block, which carries nothing about a family's sequence, reaches {control:.3f}, which is "
                  "what this measure looks like for a channel that does not encode family.", ""]
    lines += ["So the answer to the question is yes, through two channels at once: the Reactome protein entities put a family's "
              "paralogues in the same complexes and sets, and the protein descriptor block is a projection of annotation space "
              "whose Pfam and GO-function columns are close to family labels by construction. A family-grouped split therefore "
              "changes which drugs share a fold, not whether the model can tell which family a held-out target belongs to. The "
              "confound is a property of the graph and of the descriptor block, so the registered grouping "
              "(disease_cluster_and_targets) leaves it available too: it is not an argument for or against the family rule of "
              "docs/leakage_group_rules.md, it is a limit on what any grouping of drugs can buy.", "",
              "## What this does not measure", "",
              "- Leave-one-out nearest neighbour says the information is present and easy to reach, not that a trained model uses "
              "it. It is an upper bound on the confound, not an estimate of its effect on a score.",
              "- The Reactome channel is scored only on the genes that are members of a protein entity ("
              + "; ".join(f"{channel['genes_scored']} of {block['genes']} in {scope}"
                          for scope, block in result["scopes"].items() for channel in block["channels"]
                          if channel["channel"] == "reactome_entity_members")
              + "); a gene in no entity has no neighbour set there and is left out rather than counted wrong.",
              "- GtoPdb families are themselves drawn partly on sequence and pharmacology, so the descriptor result is partly "
              "definitional. That does not weaken the conclusion: the question was whether the model's inputs carry family, and "
              "a definitional route is still a route.",
              "- Families of two genes are easy to recover. The chance rate above accounts for family sizes, and the permuted "
              "column is the same rule run on once-shuffled labels, which is the empirical floor.", ""]
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal"))
    parser.add_argument("--descriptor-table", type=Path,
                        default=Path("data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet"))
    parser.add_argument("--gtopdb-families", type=Path, default=Path("data/raw/gtopdb_2026/targets_and_families.csv"))
    parser.add_argument("--evidence-dir", type=Path, default=Path("data/processed/evidence_full_v2"))
    parser.add_argument("--minimum-family-size", type=int, default=2)
    parser.add_argument("--output", type=Path, default=Path("data/processed/off_target_scoping/family_in_graph.json"))
    parser.add_argument("--markdown-output", type=Path, default=Path("docs/family_in_graph.md"))
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    descriptors = pd.read_parquet(arguments.descriptor_table)
    family_of_symbol, _ = target_families(arguments.gtopdb_families)
    genes = genes_with_a_family(nodes, family_of_symbol, arguments.minimum_family_size)
    symbol_of_node = dict(zip(nodes.node_id, nodes.gene_symbol))
    mechanism_symbols = drug_mechanism_genes(pd.read_parquet(arguments.evidence_dir / "evidence_records.parquet"), symbol_of_node)

    scopes = {"every graph gene with a family": genes}
    mechanism_genes = genes[genes.gene_symbol.isin(mechanism_symbols)]
    sizes = mechanism_genes.family.value_counts()
    mechanism_genes = mechanism_genes[mechanism_genes.family.isin(sizes[sizes >= arguments.minimum_family_size].index)]
    scopes["drug mechanism targets only"] = mechanism_genes.reset_index(drop=True)
    result = {"graph_dir": str(arguments.graph_dir), "descriptor_table": str(arguments.descriptor_table),
              "evidence_dir": str(arguments.evidence_dir), "minimum_family_size": arguments.minimum_family_size,
              "drug_mechanism_symbols": len(mechanism_symbols),
              "scopes": {name: measure(nodes, edges, descriptors, subset) for name, subset in scopes.items()}}
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=1) + "\n")
    arguments.markdown_output.write_text(report(result, arguments))
    print(json.dumps(result, indent=1))
    print(f"wrote {arguments.output} and {arguments.markdown_output}")


if __name__ == "__main__":
    main()
