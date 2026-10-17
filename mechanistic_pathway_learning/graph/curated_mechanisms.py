"""Curated edges for perturbations whose established mechanism the imported sources leave without a route (decision
9.68 (a)).

Decision 9.68 measured 186 gene perturbations whose protein has no edge out, so the model can fit their labels only
through the gene's and the protein's own states. For a few of them the mechanism is established and sits inside layers
the graph already has: a co-chaperone that keeps the aromatic amino acid hydroxylases stable and active (DNAJC12), a
processing enzyme for a reaction Human-GEM leaves without a catalyst (MMACHC), and a cofactor the enzyme-regulation
source does not name (adenosylcobalamin for methylmalonyl-CoA mutase). Decision 9.71 found the same for a subunit:
GABRD belongs to no Reactome GABA-A receptor entity, so no delta-containing receptor conducted chloride, and its row
makes it a member of the GABA-bound heteropentamer entity that catalyses chloride entry. Each row of
docs/curated_dead_end_mechanisms.csv
names one edge with the papers behind it; the lysosomal release of cobalamin, a transport with no reaction to attach a
catalyst to, is a curated reaction in experiments/build_neuronal_graph_variant.py instead.

The rows are added to the neuronal variant, before the gene/protein split, so a gene end becomes the gene's protein node
as every imported edge's does (graph/gene_protein_split.py), and each edge's evidence_source starts with "curated" and
names its PubMed ids, which is how the edge evidence table reads a curated edge (graph/edge_evidence.py, decision 9.57).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DEFAULT_CURATED_MECHANISMS_PATH = Path("docs/curated_dead_end_mechanisms.csv")
REQUIRED_COLUMNS = ("source_id", "target_id", "relation_type", "sign", "decision", "pmids", "quotes", "reading")
# the sign each relation carries, and the node types its two ends must have before the split
SIGN_OF_RELATION = {"activates": 1.0, "inhibits": -1.0, "catalyzed_by": 1.0, "cofactor_of": 1.0, "member_of": 1.0}
NODE_TYPES_OF_RELATION = {"activates": ("gene", "gene"), "inhibits": ("gene", "gene"), "catalyzed_by": ("gene", "reaction"),
                          "cofactor_of": ("metabolite", "gene"), "member_of": ("gene", "protein_entity")}
EVIDENCE_PREFIX = "curated dead-end mechanism"


def read_curated_mechanisms(path: Path = DEFAULT_CURATED_MECHANISMS_PATH) -> pd.DataFrame:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    return table


def pubmed_ids_of(cell: str) -> list[str]:
    """The PubMed ids of a row, from a ';'-separated cell, each all digits."""
    identifiers = [part.strip() for part in str(cell).split(";") if part.strip()]
    for identifier in identifiers:
        if not identifier.isdigit():
            raise ValueError(f"not a PubMed id: {identifier!r}")
    return identifiers


def problems_of_row(row: pd.Series, node_type_of: dict[str, str], relation_types: list[str]) -> list[str]:
    """Why a row cannot be added as written (empty when it can)."""
    problems = []
    relation = row["relation_type"]
    if relation not in SIGN_OF_RELATION:
        return [f"relation {relation!r} is not one this table may add ({', '.join(sorted(SIGN_OF_RELATION))})"]
    if relation not in relation_types:
        problems.append(f"relation {relation!r} is not among the graph's relation types")
    try:
        sign = float(row["sign"])
    except ValueError:
        sign = float("nan")
    if sign != SIGN_OF_RELATION[relation]:
        problems.append(f"sign {row['sign']!r} does not match {relation} ({SIGN_OF_RELATION[relation]:+.0f})")
    for end, column, expected in (("source", "source_id", NODE_TYPES_OF_RELATION[relation][0]), ("target", "target_id", NODE_TYPES_OF_RELATION[relation][1])):
        node = row[column]
        if node not in node_type_of:
            problems.append(f"{end} {node} is not a node of the graph")
        elif node_type_of[node] != expected:
            problems.append(f"{end} {node} is a {node_type_of[node]} node, and {relation} needs a {expected} node there")
    try:
        if len(pubmed_ids_of(row["pmids"])) == 0:
            problems.append("no PubMed id")
    except ValueError as error:
        problems.append(str(error))
    if not str(row["quotes"]).strip():
        problems.append("no quote from the papers")
    return problems


def add_curated_mechanisms(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str],
                           table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with each row of the table added as one edge. A row that cannot be added as written raises, so the
    table cannot drift from the graph silently; a row whose edge the graph already holds (the same source, target and
    relation) is left out and reported, since the imported edge comes first."""
    node_type_of = dict(zip(nodes.node_id, nodes.node_type))
    problems = {index: problems_of_row(row, node_type_of, relation_types) for index, row in table.iterrows()}
    problems = {index: found for index, found in problems.items() if found}
    if problems:
        details = "; ".join(f"row {index + 2} ({table.loc[index, 'source_id']} -> {table.loc[index, 'target_id']}): {', '.join(found)}" for index, found in problems.items())
        raise ValueError(f"curated mechanisms that cannot be added: {details}")
    existing = set(zip(edges.source_id, edges.target_id, edges.relation_type))
    rows, already_present = [], []
    for row in table.itertuples():
        key = (row.source_id, row.target_id, row.relation_type)
        if key in existing:
            already_present.append(f"{row.source_id} -{row.relation_type}-> {row.target_id}")
            continue
        existing.add(key)
        pmids = "; ".join(f"PMID {identifier}" for identifier in pubmed_ids_of(row.pmids))
        rows.append({"source_id": row.source_id, "target_id": row.target_id, "relation_type": row.relation_type,
                     "sign": SIGN_OF_RELATION[row.relation_type], "evidence_source": f"{EVIDENCE_PREFIX}, decision {row.decision} ({pmids})"})
    new_edges = (pd.concat([edges, pd.DataFrame(rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
                 if rows else edges.reset_index(drop=True).copy())
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes = nodes.copy()
    if "degree" in new_nodes.columns:
        new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    summary = {"edges_added": len(rows), "by_relation": pd.Series([row["relation_type"] for row in rows], dtype=str).value_counts().to_dict(),
               "already_present": already_present}
    return new_nodes, new_edges, list(relation_types), summary
