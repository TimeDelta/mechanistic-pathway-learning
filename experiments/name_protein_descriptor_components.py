"""Give each protein_rrr component a readable name: the function annotations its values correlate with most, on each
side (docs/node_descriptor_columns.md reads the result).

A component of the reduced-rank regression (mechanistic_pathway_learning/graph/protein_descriptors.py) is a direction
in the space of annotations best predicted from the ESM-2 embedding; it mixes many annotations, and its sign is the one
the fit gave. Over the proteins the fit used (those with any GO annotation), each component is correlated with every
annotation target column, as the targets entered the fit. The positive pole is the annotations with the largest
positive correlation, the negative pole those with the most negative. The name reads "high: <first term of the
positive pole>; low: <first term of the negative pole>": proteins with a high value look like the first, a low value
like the second. A term is skipped when its proteins nearly repeat those of a term already listed on that pole
(correlation of the two annotation columns above SAME_PROTEINS_CORRELATION), since GO terms are propagated to their
ancestors. The correlations show how well the name fits: a component whose top |r| is 0.2 is not "about" its top term.

Name sources (downloaded into their own directories under data/raw): GO term names from go-basic.obo; Pfam family
descriptions from Pfam-A.clans.tsv.gz (EBI Pfam current release); EC class names from ExPASy enzclass.txt.

    OMP_NUM_THREADS=2 python experiments/name_protein_descriptor_components.py
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.protein_descriptors import (
    annotation_target_matrix,
    go_annotations,
    read_go_ancestors,
    weighted_standardised_targets,
)

TERMS_PER_POLE = 4
SAME_PROTEINS_CORRELATION = 0.8
BLOCK_LABELS = {"ec_level_2": "EC", "go_function": "GO function", "go_component": "GO component", "location": "UniProt location", "pfam": "Pfam"}


def go_term_names(go_obo_path: Path) -> dict[str, str]:
    names, current = {}, None
    with open(go_obo_path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("id: GO:"):
                current = line[4:].strip()
            elif line.startswith("name: ") and current:
                names[current] = line[6:].strip()
                current = None
    return names


def pfam_family_names(clans_path: Path) -> dict[str, str]:
    """Pfam accession -> family description (column 5 of Pfam-A.clans.tsv)."""
    names = {}
    with gzip.open(clans_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                names[fields[0]] = " ".join(fields[4].split())
    return names


def ec_class_names(enzclass_path: Path) -> dict[str, str]:
    """'1.14' -> its subclass name, from lines such as '1.14. -.-  Acting on paired donors ...'."""
    names = {}
    for line in enzclass_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^\s*(\d+)\.\s*(\d+|-)\.\s*(\d+|-)\.-\s+(.*)$", line)
        if match and match.group(2) != "-" and match.group(3) == "-":
            names[f"{match.group(1)}.{match.group(2)}"] = match.group(4).strip().rstrip(".")
        elif match and match.group(2) == "-":
            names[match.group(1)] = match.group(4).strip().rstrip(".")
    return names


def readable(column: str, go_names: dict[str, str], pfam_names: dict[str, str], ec_names: dict[str, str]) -> str:
    block, term = column.split(":", 1)
    if block in ("go_function", "go_component"):
        return go_names.get(term, term)
    if block == "pfam":
        return f"{pfam_names.get(term, term)} ({term})"
    if block == "ec_level_2":
        top_class = ec_names.get(term.split(".")[0], "")
        return f"EC {term} {top_class.lower()}: {ec_names.get(term, '').lower()}".rstrip(": ")
    return term


def pole(correlation: pd.Series, targets: np.ndarray, column_position: dict[str, int], descending: bool) -> list[str]:
    chosen: list[str] = []
    for column in correlation.sort_values(ascending=not descending).index:
        if (correlation[column] > 0) != descending or len(chosen) == TERMS_PER_POLE:
            break
        candidate = targets[:, column_position[column]]
        if any(abs(np.corrcoef(candidate, targets[:, column_position[other]])[0, 1]) > SAME_PROTEINS_CORRELATION for other in chosen):
            continue
        chosen.append(column)
    return chosen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--gaf", type=Path, default=Path("data/raw/gene_ontology/goa_human.gaf.gz"))
    parser.add_argument("--go-obo", type=Path, default=Path("data/raw/gene_ontology/go-basic.obo"))
    parser.add_argument("--pfam-names", type=Path, default=Path("data/raw/pfam_names_2026/Pfam-A.clans.tsv.gz"))
    parser.add_argument("--ec-names", type=Path, default=Path("data/raw/expasy_enzclass_2026/enzclass.txt"))
    parser.add_argument("--per-entry-descriptors", type=Path, default=Path("data/processed/node_descriptors/protein_descriptors_per_entry.parquet"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/node_descriptors/protein_rrr_component_names.json"))
    arguments = parser.parse_args()

    descriptors = pd.read_parquet(arguments.per_entry_descriptors)
    uniprot = pd.read_csv(arguments.uniprot_table, sep="\t", usecols=["Entry", "Gene Names (primary)", "EC number", "Subcellular location [CC]", "Pfam"])
    uniprot = uniprot[uniprot.Entry.isin(descriptors.index)].reset_index(drop=True)
    targets, block_columns, observed = annotation_target_matrix(uniprot, go_annotations(arguments.gaf, read_go_ancestors(arguments.go_obo)))
    target_values = weighted_standardised_targets(targets, block_columns, observed)
    fit_rows = (observed["go_function"] | observed["go_component"]).to_numpy()
    fitted_targets = target_values[fit_rows]
    fitted_descriptors = descriptors.loc[uniprot.Entry[fit_rows]].to_numpy(dtype=np.float64)

    standardised_targets = (fitted_targets - fitted_targets.mean(axis=0)) / np.where(fitted_targets.std(axis=0) > 0, fitted_targets.std(axis=0), 1.0)
    standardised_descriptors = (fitted_descriptors - fitted_descriptors.mean(axis=0)) / fitted_descriptors.std(axis=0)
    correlations = standardised_descriptors.T @ standardised_targets / len(standardised_targets)  # [components, targets]

    go_names, pfam_names, ec_names = go_term_names(arguments.go_obo), pfam_family_names(arguments.pfam_names), ec_class_names(arguments.ec_names)
    column_position = {column: position for position, column in enumerate(targets.columns)}
    components = {}
    for index, component in enumerate(descriptors.columns):
        correlation = pd.Series(correlations[index], index=targets.columns)
        poles = {}
        for side, descending in (("positive", True), ("negative", False)):
            poles[side] = [{"annotation": column, "block": BLOCK_LABELS[column.split(":", 1)[0]],
                            "name": readable(column, go_names, pfam_names, ec_names), "correlation": round(float(correlation[column]), 3)}
                           for column in pole(correlation, fitted_targets, column_position, descending)]
        name = "; ".join(f"{label}: {entries[0]['name']}" for label, entries in (("high", poles["positive"]), ("low", poles["negative"])) if entries)
        components[f"protein_{component}"] = {"name": name, "largest_absolute_correlation": round(float(np.abs(correlation).max()), 3), **poles}
    result = {"per_entry_descriptors": str(arguments.per_entry_descriptors),
              "per_entry_descriptors_sha256": hashlib.sha256(arguments.per_entry_descriptors.read_bytes()).hexdigest(),
              "proteins_correlated": int(fit_rows.sum()), "annotation_columns": int(targets.shape[1]), "components": components}
    arguments.output.write_text(json.dumps(result, indent=1) + "\n")
    for column, entry in components.items():
        print(f"{column}: {entry['name']} (max |r| {entry['largest_absolute_correlation']})")
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
