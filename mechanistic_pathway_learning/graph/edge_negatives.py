"""IntAct's curated negative interactions as an anchor set for the graph's protein edges.

What it is. IntAct records, beside its interactions, the pairs curators found reported as *not* interacting
(intact_negative.txt, PSI-MITAB 2.7). They are the one curated statement of "these two proteins do not bind" the
project has, and the report of 15 October 2026 on evidence priors names such a set as what the reliability of an
edge source cannot be estimated without: every edge source asserts presences only, so its specificity is
unidentifiable until something asserts an absence. The set is small (908 human pairs on 14 January 2026), so what it
supports is an audit figure per source, an upper bound on how often a source asserts a pair curators recorded as
non-interacting, and not a per-edge weight.

Licence: IntAct is an EMBL-EBI resource; the EMBL-EBI terms of use (read 15 October 2026) impose no restriction
beyond the data owner's, and the IntAct data are released under CC BY 4.0 according to the resource's own
statements, which this session could not read from the script-rendered legal page; docs/data_sources.md names the
licence as read from the terms of use and asks the owner to confirm.
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

HUMAN_TAXID = "taxid:9606"
GENE_NAME = re.compile(r"uniprotkb:([^|()]+)\(gene name\)")
ACCESSION = re.compile(r"uniprotkb:([A-Z][0-9A-Z]{5,9})(?:-\d+)?")


@dataclass(frozen=True)
class NegativePair:
    accession_a: str
    accession_b: str
    gene_a: str | None
    gene_b: str | None
    interaction_type: str
    detection_method: str
    publications: str


def read_human_negative_pairs(path: Path) -> list[NegativePair]:
    """Every human-human row, with the isoform suffix of an accession dropped and the gene names read from the aliases."""
    pairs = []
    with open(path, encoding="utf-8", newline="") as negative_file:
        reader = csv.reader(negative_file, delimiter="\t")
        header = next(reader)
        column = {name.lstrip("#"): index for index, name in enumerate(header)}
        for row in reader:
            if HUMAN_TAXID not in row[column["Taxid interactor A"]] or HUMAN_TAXID not in row[column["Taxid interactor B"]]:
                continue
            accession_a = ACCESSION.search(row[column["ID(s) interactor A"]])
            accession_b = ACCESSION.search(row[column["ID(s) interactor B"]])
            if not accession_a or not accession_b:
                continue
            gene_a = GENE_NAME.search(row[column["Alias(es) interactor A"]])
            gene_b = GENE_NAME.search(row[column["Alias(es) interactor B"]])
            pairs.append(NegativePair(accession_a.group(1), accession_b.group(1), gene_a.group(1) if gene_a else None, gene_b.group(1) if gene_b else None,
                                      row[column["Interaction type(s)"]], row[column["Interaction detection method(s)"]], row[column["Publication Identifier(s)"]]))
    return pairs


def symbol_of_node(nodes: pd.DataFrame) -> dict[str, str]:
    """node id -> gene symbol for gene and protein nodes (a protein node's symbol is its gene's)."""
    symbols = {}
    for row in nodes.itertuples(index=False):
        symbol = getattr(row, "gene_symbol", None)
        if isinstance(symbol, str) and symbol and getattr(row, "node_type", "") in ("gene", "protein"):
            symbols[str(row.node_id)] = symbol
    return symbols


def contradicted_edges(pairs: list[NegativePair], edges: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """The graph edges whose two ends are a curated non-interacting pair, in either direction, with the pair's provenance."""
    negative_by_symbols: dict[frozenset[str], NegativePair] = {}
    for pair in pairs:
        if pair.gene_a and pair.gene_b:
            negative_by_symbols[frozenset((pair.gene_a, pair.gene_b))] = pair
    symbols = symbol_of_node(nodes)
    rows = []
    for edge in edges.itertuples(index=False):
        symbol_a, symbol_b = symbols.get(str(edge.source_id)), symbols.get(str(edge.target_id))
        if not symbol_a or not symbol_b or symbol_a == symbol_b:
            continue
        pair = negative_by_symbols.get(frozenset((symbol_a, symbol_b)))
        if pair is None:
            continue
        rows.append({"source_id": edge.source_id, "target_id": edge.target_id, "relation_type": getattr(edge, "relation_type", ""),
                     "sign": getattr(edge, "sign", None), "evidence_source": getattr(edge, "evidence_source", ""),
                     "intact_interaction_type": pair.interaction_type, "intact_detection_method": pair.detection_method, "intact_publications": pair.publications})
    return pd.DataFrame(rows, columns=["source_id", "target_id", "relation_type", "sign", "evidence_source", "intact_interaction_type", "intact_detection_method", "intact_publications"])


def audit(pairs: list[NegativePair], edges: pd.DataFrame, nodes: pd.DataFrame) -> dict:
    """Per edge source: how many of its protein-protein edges fall on a curated non-interacting pair, against how many
    of its edges join two proteins the negative set names at all (the exposure), so a source is judged on the pairs
    the anchor set can speak to."""
    named = {gene for pair in pairs for gene in (pair.gene_a, pair.gene_b) if gene}
    symbols = symbol_of_node(nodes)
    exposed: Counter = Counter()
    for edge in edges.itertuples(index=False):
        symbol_a, symbol_b = symbols.get(str(edge.source_id)), symbols.get(str(edge.target_id))
        if symbol_a and symbol_b and symbol_a != symbol_b and symbol_a in named and symbol_b in named:
            exposed[str(getattr(edge, "evidence_source", ""))] += 1
    contradicted = contradicted_edges(pairs, edges, nodes)
    by_source: dict[str, dict] = defaultdict(lambda: {"edges_on_named_proteins": 0, "edges_on_a_negative_pair": 0})
    for source, count in exposed.items():
        by_source[source]["edges_on_named_proteins"] = count
    for source, count in Counter(contradicted.evidence_source.astype(str)).items():
        by_source[source]["edges_on_a_negative_pair"] = int(count)
    return {"human_negative_pairs": len(pairs), "pairs_with_both_gene_names": sum(1 for pair in pairs if pair.gene_a and pair.gene_b),
            "proteins_named": len(named), "edges_on_a_negative_pair": int(len(contradicted)), "by_source": dict(sorted(by_source.items())),
            "contradicted": contradicted.to_dict("records")}


def write_audit(path: Path, result: dict) -> None:
    path.write_text(json.dumps(result, indent=1, default=str) + "\n", encoding="utf-8")
