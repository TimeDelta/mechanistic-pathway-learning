"""Gene-disease association types behind HPO annotations (design review v0.4, finding 4).

HPO's genes_to_phenotype.txt inherits every gene-disease association of genes_to_disease.txt:
OMIM associations typed MENDELIAN or POLYGENIC (susceptibility loci) and every Orphanet
association, which HPO types UNKNOWN. Grade A of design section 4.3 asks for a human loss
of function in a graph gene, so an annotation whose only gene-disease link is a polygenic
locus, a susceptibility factor, a candidate gene, a modifier or a biomarker is not grade A
evidence; it is a human association (grade C). Orphanet's own association types come from
Orphadata product 6 (en_product6.xml, CC BY 4.0): disease-causing germline or somatic
mutations count as causal, the rest do not.
"""
from __future__ import annotations

import csv
import xml.etree.ElementTree as ElementTree
from collections.abc import Mapping
from pathlib import Path

CAUSAL_ORPHANET_ASSOCIATION_TYPES = frozenset({
    "Disease-causing germline mutation(s) in",
    "Disease-causing germline mutation(s) (loss of function) in",
    "Disease-causing germline mutation(s) (gain of function) in",
    "Disease-causing somatic mutation(s) in",
})
NON_CAUSAL_ORPHANET_ASSOCIATION_TYPES = frozenset({
    "Major susceptibility factor in",
    "Candidate gene tested in",
    "Role in the phenotype of",
    "Part of a fusion gene in",
    "Modifying germline mutation in",
    "Biomarker tested in",
})
UNKNOWN_ASSOCIATION_TYPE = "unknown"
GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE = "Disease-causing germline mutation(s) (gain of function) in"
assert GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE in CAUSAL_ORPHANET_ASSOCIATION_TYPES


def load_omim_association_types(genes_to_disease_path: Path) -> dict[tuple[str, str], str]:
    """(gene symbol, disease id) -> association_type column of genes_to_disease.txt (MENDELIAN, POLYGENIC or UNKNOWN)."""
    types: dict[tuple[str, str], str] = {}
    with open(genes_to_disease_path, encoding="utf-8") as handle:
        header_line = handle.readline().lstrip("#").strip()
        reader = csv.DictReader(handle, fieldnames=header_line.split("\t"), delimiter="\t")
        for row in reader:
            types[(row["gene_symbol"], row["disease_id"])] = row["association_type"].strip()
    return types


def load_orphanet_association_types(product6_path: Path) -> dict[tuple[str, str], str]:
    """(gene symbol, ORPHA:<code>) -> Orphadata product 6 association type name."""
    types: dict[tuple[str, str], str] = {}
    root = ElementTree.parse(product6_path).getroot()
    for disorder in root.iter("Disorder"):
        orpha_code = disorder.findtext("OrphaCode")
        if not orpha_code:
            continue
        for association in disorder.iter("DisorderGeneAssociation"):
            gene = association.find("Gene")
            type_element = association.find("DisorderGeneAssociationType")
            symbol = gene.findtext("Symbol") if gene is not None else None
            type_name = type_element.findtext("Name") if type_element is not None else None
            if symbol and type_name:
                types[(symbol, f"ORPHA:{orpha_code}")] = type_name
    return types


def association_type_for(
    gene_symbol: str,
    disease_id: str,
    omim_types: Mapping[tuple[str, str], str] | None,
    orphanet_types: Mapping[tuple[str, str], str] | None,
) -> tuple[str, bool | None]:
    """The association type label behind one (gene, disease) annotation and whether it is causal.

    Returns (label, causal) with causal True for MENDELIAN or a disease-causing Orphanet type, False for
    POLYGENIC or a non-causal Orphanet type, and None when neither table types the pair.
    """
    if disease_id.startswith("ORPHA:") and orphanet_types is not None:
        orphanet_type = orphanet_types.get((gene_symbol, disease_id))
        if orphanet_type in CAUSAL_ORPHANET_ASSOCIATION_TYPES:
            return orphanet_type, True
        if orphanet_type in NON_CAUSAL_ORPHANET_ASSOCIATION_TYPES:
            return orphanet_type, False
        if orphanet_type:
            return orphanet_type, None
    if omim_types is not None:
        omim_type = omim_types.get((gene_symbol, disease_id))
        if omim_type == "MENDELIAN":
            return omim_type, True
        if omim_type == "POLYGENIC":
            return omim_type, False
        if omim_type and omim_type != "UNKNOWN":
            return omim_type, None
    return UNKNOWN_ASSOCIATION_TYPE, None


def genes_with_gain_of_function_associations_only(causal_association_types_by_gene: Mapping[str, set[str]]) -> set[str]:
    """The genes whose every causal association behind a report is Orphanet's gain-of-function type.

    causal_association_types_by_gene holds, per gene, the association types of its reports with a causal association
    (rubric_causal_association 1). A gene qualifies only when that set is the gain-of-function type alone: one causal
    report of another type (MENDELIAN, which OMIM gives without a direction, or an Orphanet loss-of-function or
    undirected type) leaves the gene as it was, because a gene is one perturbation with one sign and its reports would
    then disagree about it. Non-causal reports (a candidate gene, a polygenic locus) state no mechanism and do not count.
    """
    return {gene for gene, association_types in causal_association_types_by_gene.items()
            if set(association_types) == {GAIN_OF_FUNCTION_ORPHANET_ASSOCIATION_TYPE}}
