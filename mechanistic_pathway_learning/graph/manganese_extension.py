"""Manganese in the graph, connected the way Human-GEM connects every transported metabolite: compartment copies of
Mn2+, transport reactions between them catalysed by the transporter genes, and cofactor edges from the copy in each
dependent enzyme's compartment to that enzyme (graph variant for the ablation of design section 5.7;
experiments/build_manganese_graph_variant.py).

Human-GEM has no manganese metabolite and lists no metal cofactors, so the manganese dependence of MnSOD (SOD2),
arginase (ARG1, ARG2), prolidase (PEPD) or phosphoenolpyruvate carboxykinase (PCK1, PCK2) is absent, and so is the
route by which a transporter defect (SLC39A14, SLC30A10: hypermanganesemia; SLC39A8: deficiency) reaches them.

The variant adds:
  - metabolite nodes MN2e, MN2c, MN2m and MN2g (base metabolite MN2, Mn2+; extracellular, cytosol, mitochondrion,
    Golgi), none of them currency;
  - transport reactions in Human-GEM's form (substrate_of from the source copy, product_of to the destination copy,
    catalyzed_by from each transporter gene), with directions curated per transporter, since GO annotations
    (GO:0005384, GO:0071421) name the substrate but not the direction: uptake e -> c by SLC39A14 (ZIP14), SLC39A8
    (ZIP8), SLC11A2 (DMT1), SLC11A1 (NRAMP1, phagosome to cytosol) and TRPM2 (cation channel); export c -> e by
    SLC30A10; Golgi loading c -> g by the secretory pathway Ca2+/Mn2+ ATPases ATP2C1 and ATP2C2; mitochondrial entry
    c -> m without a catalyst in the graph (the uniporter MCU is not a gene node of the slice), so the mitochondrial
    pool follows the cytosolic one;
  - cofactor_of edges, sign +1 (less manganese, less activity), from the copy in the enzyme's UniProt location
    (mitochondrion, else Golgi, else cytosol) to each gene whose UniProt cofactor annotation requires Mn(2+): every
    Name= entry of its COFACTOR statements is Mn(2+). Proteins that accept Mn(2+) or another metal (Mg(2+), Co(2+);
    most glycosyltransferases and phosphatases) are left out, since a fall in manganese need not slow them.
A transporter missing from the graph is added as a gene node only when it has no other way in (SLC30A10, the exporter
whose loss causes hypermanganesemia with dystonia); other absent genes are reported. Degrees are recomputed as in plus
out edges.
"""
from __future__ import annotations

import re

import pandas as pd

MANGANESE_BASE_ID = "MN2"
COMPARTMENT_NAMES = {"e": "extracellular", "c": "cytosol", "m": "mitochondrion", "g": "Golgi apparatus"}
COFACTOR_RELATION = "cofactor_of"
# reaction id, display name, source compartment, destination compartment, catalysing genes
MANGANESE_TRANSPORT_REACTIONS = [
    ("MAR_MN2_UPTAKE", "transport of Mn2+ (extracellular to cytosol)", "e", "c", ("SLC39A14", "SLC39A8", "SLC11A2", "SLC11A1", "TRPM2")),
    ("MAR_MN2_EXPORT", "transport of Mn2+ (cytosol to extracellular)", "c", "e", ("SLC30A10",)),
    ("MAR_MN2_GOLGI", "transport of Mn2+ (cytosol to Golgi apparatus)", "c", "g", ("ATP2C1", "ATP2C2")),
    ("MAR_MN2_MITOCHONDRIA", "transport of Mn2+ (cytosol to mitochondrion)", "c", "m", ()),
]
GENES_ADDED_WHEN_ABSENT = ("SLC30A10",)
METAL_NAME_PATTERN = re.compile(r"Name=([^;]+);")


def requires_manganese(cofactor_text: str) -> bool:
    """True when every cofactor named in the UniProt COFACTOR statements is Mn(2+)."""
    if not isinstance(cofactor_text, str) or "Mn(2+)" not in cofactor_text:
        return False
    return {name.strip() for name in METAL_NAME_PATTERN.findall(cofactor_text)} == {"Mn(2+)"}


def enzyme_compartment(location_text: str) -> str:
    text = location_text if isinstance(location_text, str) else ""
    if "Mitochondrion" in text:
        return "m"
    if "Golgi" in text:
        return "g"
    return "c"


def manganese_enzymes(cofactor_table: pd.DataFrame) -> tuple[dict[str, str], set[str]]:
    """({gene requiring Mn(2+): compartment}, genes accepting Mn(2+) among other metals) from UniProt columns
    'Gene Names (primary)', 'Cofactor' and 'Subcellular location [CC]'."""
    required, alternative = {}, set()
    for gene, cofactor, location in zip(cofactor_table["Gene Names (primary)"], cofactor_table["Cofactor"], cofactor_table["Subcellular location [CC]"]):
        if not isinstance(gene, str):
            continue
        if requires_manganese(cofactor):
            required[gene] = enzyme_compartment(location)
        elif isinstance(cofactor, str) and "Mn(2+)" in cofactor:
            alternative.add(gene)
    return required, alternative


def add_manganese(nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], required_enzymes: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the Mn2+ compartment copies, the transport reactions and the cofactor edges; the summary lists the
    genes connected, added and absent."""
    gene_nodes = set(nodes.loc[nodes.node_type == "gene", "node_id"])
    blank = {column: None for column in nodes.columns}
    new_node_rows, edge_rows, added_genes, absent_transporters = [], [], [], set()
    for compartment, name in COMPARTMENT_NAMES.items():
        new_node_rows.append({**blank, "node_id": f"{MANGANESE_BASE_ID}{compartment}", "node_type": "metabolite", "display_name": "Mn2+",
                              "compartment": compartment, "base_metabolite_id": MANGANESE_BASE_ID, "is_currency": False, "in_metabolic_layer": True, "degree": 0})
    for gene in GENES_ADDED_WHEN_ABSENT:
        if f"GENE:{gene}" not in gene_nodes:
            new_node_rows.append({**blank, "node_id": f"GENE:{gene}", "node_type": "gene", "display_name": gene, "gene_symbol": gene, "in_metabolic_layer": True, "degree": 0})
            gene_nodes.add(f"GENE:{gene}")
            added_genes.append(gene)
    connected_transporters = []
    for reaction_id, display_name, source, destination, genes in MANGANESE_TRANSPORT_REACTIONS:
        catalysts = [gene for gene in genes if f"GENE:{gene}" in gene_nodes]
        absent_transporters |= set(genes) - set(catalysts)
        connected_transporters += catalysts
        new_node_rows.append({**blank, "node_id": reaction_id, "node_type": "reaction", "display_name": display_name, "compartment": f"{source};{destination}",
                              "is_transport": True, "reversible": False, "gene_reaction_rule": " or ".join(catalysts), "subsystem": "Transport reactions",
                              "is_currency": False, "degree": 0})
        edge_rows.append({"source_id": f"{MANGANESE_BASE_ID}{source}", "target_id": reaction_id, "relation_type": "substrate_of", "sign": 1.0, "evidence_source": "manganese extension"})
        edge_rows.append({"source_id": reaction_id, "target_id": f"{MANGANESE_BASE_ID}{destination}", "relation_type": "product_of", "sign": 1.0, "evidence_source": "manganese extension"})
        edge_rows += [{"source_id": f"GENE:{gene}", "target_id": reaction_id, "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "manganese extension (GO manganese transport)"}
                      for gene in catalysts]
    connected_enzymes = {gene: compartment for gene, compartment in required_enzymes.items() if f"GENE:{gene}" in gene_nodes}
    edge_rows += [{"source_id": f"{MANGANESE_BASE_ID}{compartment}", "target_id": f"GENE:{gene}", "relation_type": COFACTOR_RELATION, "sign": 1.0,
                   "evidence_source": "UniProt cofactor Mn(2+)"} for gene, compartment in sorted(connected_enzymes.items())]
    new_nodes = pd.concat([nodes, pd.DataFrame(new_node_rows, columns=nodes.columns)], ignore_index=True)
    new_edges = pd.concat([edges, pd.DataFrame(edge_rows, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    new_nodes["degree"] = counts.reindex(new_nodes.node_id).fillna(0).astype("int64").to_numpy()
    new_relations = list(relation_types) + ([COFACTOR_RELATION] if COFACTOR_RELATION not in relation_types else [])
    summary = {"transporters_connected": sorted(set(connected_transporters)), "transporters_absent_from_graph": sorted(absent_transporters),
               "gene_nodes_added": added_genes, "cofactor_enzymes_connected": {gene: compartment for gene, compartment in sorted(connected_enzymes.items())},
               "cofactor_enzymes_absent_from_graph": sorted(set(required_enzymes) - set(connected_enzymes))}
    return new_nodes, new_edges, new_relations, summary
