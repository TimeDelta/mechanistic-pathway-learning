"""Build the per-edge evidence table of a graph for the edge priors (decision 9.57).

Reads the raw OmniPath and CollecTRI tables (their references and sources, which the graph build fetches and drops),
IntAct's curated non-interacting pairs and the graph's edges; writes <graph-dir>/edge_evidence.parquet aligned to
edges.parquet row for row (family, kind, records, contested) and edge_evidence_summary.json. The graph files themselves
are not touched, so no graph hash moves. The priors are computed from this table by the loader under --edge-prior-strength.

Usage:
  python experiments/build_edge_evidence.py --graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.edge_evidence import (
    DEFAULT_TIERS_PATH,
    edge_evidence_table,
    edge_priors,
    read_edge_tiers,
    signalling_references,
    summarise,
    transcription_references,
)
from mechanistic_pathway_learning.graph.edge_negatives import read_human_negative_pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--omnipath", type=Path, default=Path("data/raw/omnipath/omnipath_interactions.tsv"))
    parser.add_argument("--collectri", type=Path, default=Path("data/raw/omnipath/collectri_interactions.tsv"))
    parser.add_argument("--intact-negative", type=Path, default=Path("data/raw/intact/intact_negative.txt"))
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS_PATH)
    arguments = parser.parse_args()
    columns = ["source", "target", "source_genesymbol", "target_genesymbol", "sources", "references"]
    references = pd.concat([signalling_references(pd.read_csv(arguments.omnipath, sep="\t", usecols=columns)),
                            transcription_references(pd.read_csv(arguments.collectri, sep="\t", usecols=columns))], ignore_index=True)
    negative_pairs = {frozenset((pair.gene_a, pair.gene_b)) for pair in read_human_negative_pairs(arguments.intact_negative) if pair.gene_a and pair.gene_b}
    edges = pd.read_parquet(arguments.graph_dir / "edges.parquet")
    evidence = edge_evidence_table(edges, references, negative_pairs)
    evidence.to_parquet(arguments.graph_dir / "edge_evidence.parquet", index=False)
    gated, probability, strength = edge_priors(evidence, read_edge_tiers(arguments.tiers))
    summary = {"graph_dir": str(arguments.graph_dir), "tiers": str(arguments.tiers), **summarise(evidence, gated, probability, strength)}
    (arguments.graph_dir / "edge_evidence_summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {arguments.graph_dir / 'edge_evidence.parquet'}: {summary['gated']:,} of {summary['edges']:,} edges gated, {summary['contested']} contested")
    for key, value in summary["by_family_and_kind"].items():
        print(f"  {key:34s} {value['edges']:7,} edges, mean prior {value['mean_prior'] if value['mean_prior'] is None else round(value['mean_prior'], 3)}, mean strength {value['mean_strength']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
