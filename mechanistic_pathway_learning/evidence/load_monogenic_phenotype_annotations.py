"""Evidence class E1: monogenic gene-to-phenotype annotations from HPO (design section 4.2).

Input files (https://hpo.jax.org, annotation downloads):
  genes_to_phenotype.txt  tab-separated with a header line; columns include
                          ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, disease_id
  hp.obo                  the ontology, used to expand each target symptom to its descendants

Output: one EvidenceRecord per (gene, target symptom) with evidence_class "monogenic".
The grade A bar (two independent case series or an OMIM clinical synopsis) needs
disease-level provenance that this file does not carry; the loader records the
disease identifiers so assign_evidence_grades can be applied after OMIM lookup.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

from mechanistic_pathway_learning.evidence.assign_evidence_grades import EvidenceRecord


def load_hpo_is_a_parents_from_obo(obo_path: Path) -> dict[str, set[str]]:
    """Return term id -> set of direct parent ids from an OBO file (is_a lines only)."""
    parents_by_term: dict[str, set[str]] = defaultdict(set)
    current_term: str | None = None
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line == "[Term]":
                current_term = None
            elif line.startswith("id: HP:"):
                current_term = line[len("id: "):]
                parents_by_term.setdefault(current_term, set())
            elif line.startswith("is_a:") and current_term is not None:
                parent_id = line[len("is_a:"):].strip().split(" ")[0]
                parents_by_term[current_term].add(parent_id)
    return dict(parents_by_term)


def descendants_of(term_id: str, parents_by_term: Mapping[str, set[str]]) -> set[str]:
    """All terms whose is_a ancestry includes term_id, plus term_id itself."""
    children_by_term: dict[str, set[str]] = defaultdict(set)
    for child, parents in parents_by_term.items():
        for parent in parents:
            children_by_term[parent].add(child)
    collected = {term_id}
    frontier = [term_id]
    while frontier:
        current = frontier.pop()
        for child in children_by_term.get(current, ()):
            if child not in collected:
                collected.add(child)
                frontier.append(child)
    return collected


def parse_genes_to_phenotype(genes_to_phenotype_path: Path) -> list[dict[str, str]]:
    """Rows of genes_to_phenotype.txt as dictionaries keyed by the header names."""
    with open(genes_to_phenotype_path, encoding="utf-8") as annotation_file:
        header_line = annotation_file.readline().lstrip("#").strip()
        field_names = header_line.split("\t")
        reader = csv.DictReader(annotation_file, fieldnames=field_names, delimiter="\t")
        return [row for row in reader]


def monogenic_evidence_records(
    annotation_rows: Iterable[dict[str, str]],
    target_symptom_to_hpo_ids: Mapping[str, Iterable[str]],
    parents_by_term: Mapping[str, set[str]],
    genes_in_graph: set[str],
) -> list[EvidenceRecord]:
    """One record per (gene symbol, target symptom) for genes present in the physiology graph."""
    hpo_id_to_target_symptom: dict[str, str] = {}
    for target_symptom, hpo_ids in target_symptom_to_hpo_ids.items():
        for hpo_id in hpo_ids:
            for descendant_id in descendants_of(hpo_id, parents_by_term):
                hpo_id_to_target_symptom[descendant_id] = target_symptom
    disease_ids_by_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in annotation_rows:
        gene_symbol = row.get("gene_symbol", "")
        if gene_symbol not in genes_in_graph:
            continue
        target_symptom = hpo_id_to_target_symptom.get(row.get("hpo_id", ""))
        if target_symptom is None:
            continue
        disease_ids_by_pair[(gene_symbol, target_symptom)].add(row.get("disease_id", ""))
    return [
        EvidenceRecord(
            perturbation_identifier=gene_symbol,
            symptom_identifier=target_symptom,
            relation="induces",
            evidence_class="monogenic",
            source="HPO genes_to_phenotype; diseases=" + ";".join(sorted(disease_ids)),
            independent_case_series_count=len(disease_ids),  # provisional until OMIM lookup replaces it
        )
        for (gene_symbol, target_symptom), disease_ids in disease_ids_by_pair.items()
    ]
