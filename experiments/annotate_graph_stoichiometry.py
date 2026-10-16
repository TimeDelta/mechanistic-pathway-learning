"""Write the stoichiometric coefficient onto the reaction edges of a graph whose rows must not move.

Why. Item 20 (built 14 October 2026, decision 8.6) has the linear-response encoder divide a node's supply by how
many copies of it one turnover consumes, read from `stoichiometric_coefficient` on the substrate_of and product_of
edges. The build writes that column, but the rebuild copies the base graph from `data/releases/v0.4/graph_full`,
which predates it, to keep the release's row order (the degree-preserving rewiring of the H2 control swaps edges by
position, so another order is another rewired graph); every variant inherits the copy, so the physiology audit of
16 October 2026 found the column absent from the graph the model reads and every demand still 1. This script adds the
column in place, by joining Human-GEM's coefficients onto the existing rows: no row is added, dropped or moved, the
node table is untouched, and the file is replaced atomically.

What is written. abs(coefficient) from the SBML for every (metabolite, reaction) substrate_of or product_of edge, 1.0
on every other edge (a relation that carries no stoichiometry), as the build writes it. A graph that already has the
column is rewritten the same way (idempotent), and a reaction edge the model does not know keeps 1.0 and is counted.

Usage:
  python experiments/annotate_graph_stoichiometry.py --graph-dirs data/processed/graph_full data/processed/graph_full_neuronal ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mechanistic_pathway_learning.documentation import provenance_stamp  # noqa: E402

STOICHIOMETRIC_RELATIONS = ("substrate_of", "product_of")


def coefficient_map(model) -> dict[tuple[str, str, str], float]:
    """(source_id, target_id, relation_type) -> abs(coefficient) for every participant of every reaction of a cobra model,
    in the nominal direction the build uses: a negative coefficient is a substrate (metabolite -> reaction), a positive
    one a product (reaction -> metabolite)."""
    coefficients: dict[tuple[str, str, str], float] = {}
    for reaction in model.reactions:
        for metabolite, coefficient in reaction.metabolites.items():
            if coefficient < 0:
                coefficients[(metabolite.id, reaction.id, "substrate_of")] = abs(float(coefficient))
            else:
                coefficients[(reaction.id, metabolite.id, "product_of")] = abs(float(coefficient))
    return coefficients


def annotate(edges: pd.DataFrame, coefficients: dict[tuple[str, str, str], float]) -> tuple[pd.DataFrame, dict]:
    """The edge table with `stoichiometric_coefficient` written from the map; rows, order and every other column as
    they were. Returns the table and the counts."""
    annotated = edges.copy()
    is_reaction_edge = annotated.relation_type.isin(STOICHIOMETRIC_RELATIONS).to_numpy()
    values = []
    unknown = 0
    for source, target, relation, reaction_edge in zip(annotated.source_id, annotated.target_id, annotated.relation_type, is_reaction_edge):
        if not reaction_edge:
            values.append(1.0)
            continue
        coefficient = coefficients.get((str(source), str(target), str(relation)))
        if coefficient is None:
            unknown += 1
            values.append(1.0)
        else:
            values.append(coefficient)
    annotated["stoichiometric_coefficient"] = pd.Series(values, index=annotated.index, dtype="float64")
    stated = annotated.loc[is_reaction_edge, "stoichiometric_coefficient"]
    counts = {"edges": int(len(annotated)), "reaction_edges": int(is_reaction_edge.sum()), "reaction_edges_unknown_to_the_model": unknown,
              "reaction_edges_above_one": int((stated > 1.0).sum()), "reaction_edges_fractional": int((stated < 1.0).sum()),
              "largest_coefficient": float(stated.max()) if len(stated) else None}
    return annotated, counts


def write_atomically(table: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + ".tmp")
    table.to_parquet(temporary, index=False)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--graph-dirs", type=Path, nargs="+", required=True)
    parser.add_argument("--sbml", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--report", type=Path, default=Path("data/processed/rebuild_logs/graph_stoichiometry.json"))
    arguments = parser.parse_args()
    import cobra

    model = cobra.io.read_sbml_model(str(arguments.sbml))
    coefficients = coefficient_map(model)
    report = {"stamp": provenance_stamp("experiments/annotate_graph_stoichiometry.py", [str(arguments.sbml), *[str(path) for path in arguments.graph_dirs]]),
              "model_reactions": len(model.reactions), "coefficients": len(coefficients), "graphs": {}}
    for graph_directory in arguments.graph_dirs:
        edges_path = graph_directory / "edges.parquet"
        edges = pd.read_parquet(edges_path)
        before = edges.drop(columns=["stoichiometric_coefficient"], errors="ignore")
        annotated, counts = annotate(edges, coefficients)
        assert annotated.drop(columns=["stoichiometric_coefficient"]).equals(before), "no row, order or other column may change"
        write_atomically(annotated, edges_path)
        report["graphs"][str(graph_directory)] = counts
        print(f"{graph_directory}: {counts}")
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
