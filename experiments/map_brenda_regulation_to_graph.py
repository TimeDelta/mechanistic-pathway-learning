"""Map BRENDA's human enzyme-regulation records onto Human-GEM metabolites and graph protein nodes.

Two routes beyond the exact-name match of the first pass (decision 9.26):

  the EC route, for the records BRENDA gives no accession for. UniProt's EC column names the human proteins of an
  EC class, but one EC class often names many isoforms - EC 1.14.14.1 names 25 - and a BRENDA statement about the
  class does not say which isoform, so fanning out would assert what the source does not. The route is therefore
  taken ONLY where the EC class names exactly one human accession, and the rest are reported as unresolved.

  the ChEBI route, for the ligand names that match no Human-GEM metabolite by name. A ligand name is resolved to a
  ChEBI term through the ontology's names and synonyms, and the term to a Human-GEM base through the project's own
  mapper, which walks conjugate acid and base pairs. Decision 9.17's guard applies to charge and not to identity,
  so a conjugate crossing is allowed here: an inhibitor named as an acid and held as its base is the same molecule
  acting on the same enzyme.
"""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import (
    ChemicalDefinition,
    human_gem_metabolites_by_chebi,
    map_definition_to_metabolites,
    read_chebi_relations,
)


def human_accessions_by_ec_number(uniprot_path: Path) -> dict[str, set[str]]:
    accessions_of_ec: dict[str, set[str]] = collections.defaultdict(set)
    with gzip.open(uniprot_path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for ec_number in [part.strip() for part in (row.get("EC number") or "").split(";") if part.strip()]:
                accessions_of_ec[ec_number].add(row["Entry"].strip())
    return accessions_of_ec


def chebi_id_by_lowercase_name(chebi_obo_path: Path) -> dict[str, str]:
    """Every ChEBI name and synonym, lowercased, to its term. First writer wins, so the primary name takes it."""
    name_to_term: dict[str, str] = {}
    opener = gzip.open if str(chebi_obo_path).endswith(".gz") else open
    current_term = None
    with opener(chebi_obo_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("id: CHEBI:"):
                current_term = line[4:].strip()
            elif current_term and line.startswith("name: "):
                name_to_term.setdefault(line[6:].strip().lower(), current_term)
            elif current_term and line.startswith("synonym: "):
                start, end = line.find('"'), line.rfind('"')
                if 0 <= start < end:
                    name_to_term.setdefault(line[start + 1:end].strip().lower(), current_term)
    return name_to_term


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--regulation-records", type=Path, required=True)
    parser.add_argument("--graph-dir", type=Path, default=Path("data/processed/graph_full_neuronal_split_binders_drugs_v4"))
    parser.add_argument("--uniprot", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_core.obo.gz"))
    parser.add_argument("--metabolites-table", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    arguments = parser.parse_args()

    regulation = pd.read_csv(arguments.regulation_records)
    nodes = pd.read_parquet(arguments.graph_dir / "nodes.parquet")
    accession_of_graph_protein = {str(value).upper() for value in nodes.uniprot_entry.dropna()}
    base_of_metabolite_display_name = {}
    for row in nodes[nodes.node_type == "metabolite"].itertuples():
        base_of_metabolite_display_name.setdefault(str(row.display_name).strip().lower(), row.base_metabolite_id)

    accessions_of_ec = human_accessions_by_ec_number(arguments.uniprot)
    term_of_name = chebi_id_by_lowercase_name(arguments.chebi_obo)
    chebi_names, chebi_neighbours, chebi_children = read_chebi_relations(arguments.chebi_obo)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(pd.read_csv(arguments.metabolites_table, sep="\t"))

    base_of_ligand_name: dict[str, str] = {}
    ligand_resolution_counts = collections.Counter()
    for ligand_name in sorted({str(value).strip() for value in regulation.ligand_name.dropna()}):
        lowered = ligand_name.lower()
        if lowered in base_of_metabolite_display_name:
            base_of_ligand_name[ligand_name] = base_of_metabolite_display_name[lowered]
            ligand_resolution_counts["matched a Human-GEM display name"] += 1
            continue
        term = term_of_name.get(lowered)
        if term is None:
            ligand_resolution_counts["no ChEBI term for the name"] += 1
            continue
        mapped, route = map_definition_to_metabolites(ChemicalDefinition("", (term,), 0, ""), chebi_names,
                                                      chebi_neighbours, chebi_children, metabolites_by_chebi)
        if mapped:
            base_of_ligand_name[ligand_name] = sorted(mapped)[0]
            ligand_resolution_counts[f"ChEBI route, {route}"] += 1
        else:
            ligand_resolution_counts[f"a ChEBI term but no Human-GEM metabolite ({route})"] += 1

    enzyme_resolution_counts = collections.Counter()
    output_rows = []
    for row in regulation.itertuples():
        base_metabolite = base_of_ligand_name.get(str(row.ligand_name).strip())
        if base_metabolite is None:
            continue
        accession = str(row.uniprot_accession).strip().upper()
        if accession and accession != "NAN":
            enzyme_route = "accession stated by BRENDA"
        else:
            candidates = accessions_of_ec.get(str(row.ec_number).strip(), set())
            if len(candidates) == 1:
                accession = next(iter(candidates))
                enzyme_route = "EC class names exactly one human accession"
            else:
                enzyme_resolution_counts["EC class names none or several human accessions"] += 1
                continue
        if accession not in accession_of_graph_protein:
            enzyme_resolution_counts["accession is not a graph protein node"] += 1
            continue
        enzyme_resolution_counts[enzyme_route] += 1
        output_rows.append({"base_metabolite_id": base_metabolite, "uniprot_accession": accession,
                            "relation_type": row.relation_type, "sign": row.sign, "ec_number": row.ec_number,
                            "brenda_field": row.field, "ligand_name": row.ligand_name,
                            "enzyme_route": enzyme_route,
                            "brenda_reference_count": row.brenda_reference_count})

    pairs = (pd.DataFrame(output_rows)
             .groupby(["base_metabolite_id", "uniprot_accession", "relation_type", "sign"], as_index=False)
             .agg(brenda_record_count=("ec_number", "size"),
                  ec_numbers=("ec_number", lambda values: ";".join(sorted(set(values))[:6])),
                  brenda_reference_count=("brenda_reference_count", "sum"),
                  enzyme_routes=("enzyme_route", lambda values: ";".join(sorted(set(values)))),
                  example_ligand_name=("ligand_name", "first")))
    pairs = pairs.sort_values(["relation_type", "base_metabolite_id", "uniprot_accession"])
    pairs.to_csv(arguments.output, index=False)

    summary = {"regulation_records_read": len(regulation),
               "distinct_ligand_names": int(regulation.ligand_name.nunique()),
               "ligand_names_resolved_to_a_metabolite": len(base_of_ligand_name),
               "ligand_resolution": dict(ligand_resolution_counts),
               "enzyme_resolution": dict(enzyme_resolution_counts),
               "pairs_written": len(pairs),
               "pairs_by_relation": pairs.relation_type.value_counts().to_dict(),
               "distinct_metabolites": int(pairs.base_metabolite_id.nunique()),
               "distinct_enzymes": int(pairs.uniprot_accession.nunique())}
    arguments.summary.write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
