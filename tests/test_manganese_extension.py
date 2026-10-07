"""Tests for the manganese graph variant: cofactor parsing, compartments, transport reactions and cofactor edges."""
import pandas as pd

from mechanistic_pathway_learning.graph.manganese_extension import add_manganese, enzyme_compartment, requires_manganese

SOLE_MANGANESE = "COFACTOR: Name=Mn(2+); Xref=ChEBI:CHEBI:29035; Evidence={ECO:0000269}; Note=Binds 1 Mn(2+) ion per subunit."
EITHER_METAL = "COFACTOR: Name=Mg(2+); Xref=ChEBI:CHEBI:18420; Name=Mn(2+); Xref=ChEBI:CHEBI:29035; Note=Mn(2+) or Mg(2+)."


def toy_graph():
    nodes = pd.DataFrame({"node_id": ["GENE:SOD2", "GENE:ARG1", "GENE:SLC39A14", "MAR1"], "node_type": ["gene", "gene", "gene", "reaction"],
                          "display_name": ["SOD2", "ARG1", "SLC39A14", "x"], "compartment": [None, None, None, "c"], "base_metabolite_id": [None] * 4,
                          "degree": [1, 0, 0, 1], "gene_symbol": ["SOD2", "ARG1", "SLC39A14", None], "is_transport": [None, None, None, False],
                          "reversible": [None, None, None, False], "gene_reaction_rule": [None] * 4, "subsystem": [None] * 4, "is_currency": [None] * 4,
                          "in_metabolic_layer": [True, True, True, None]})
    edges = pd.DataFrame({"source_id": ["GENE:SOD2"], "target_id": ["MAR1"], "relation_type": ["catalyzed_by"], "sign": [1.0], "evidence_source": ["Human-GEM GPR"]})
    return nodes, edges


def test_cofactor_parsing_and_compartments() -> None:
    assert requires_manganese(SOLE_MANGANESE) and not requires_manganese(EITHER_METAL) and not requires_manganese(float("nan"))
    assert enzyme_compartment("SUBCELLULAR LOCATION: Mitochondrion matrix.") == "m"
    assert enzyme_compartment("SUBCELLULAR LOCATION: Golgi apparatus membrane.") == "g"
    assert enzyme_compartment("SUBCELLULAR LOCATION: Cytoplasm.") == "c"


def test_transporters_are_connected_through_transport_reactions_and_enzymes_to_their_compartment_copy() -> None:
    nodes, edges = toy_graph()
    new_nodes, new_edges, relations, summary = add_manganese(nodes, edges, ["catalyzed_by"], {"SOD2": "m", "ARG1": "c", "ABSENT": "c"})
    assert {"MN2e", "MN2c", "MN2m", "MN2g", "GENE:SLC30A10", "MAR_MN2_UPTAKE", "MAR_MN2_EXPORT"} <= set(new_nodes.node_id)
    def has(source, target, relation):
        return bool(((new_edges.source_id == source) & (new_edges.target_id == target) & (new_edges.relation_type == relation)).any())
    assert has("MN2e", "MAR_MN2_UPTAKE", "substrate_of") and has("MAR_MN2_UPTAKE", "MN2c", "product_of") and has("GENE:SLC39A14", "MAR_MN2_UPTAKE", "catalyzed_by")
    assert has("MN2c", "MAR_MN2_EXPORT", "substrate_of") and has("MAR_MN2_EXPORT", "MN2e", "product_of") and has("GENE:SLC30A10", "MAR_MN2_EXPORT", "catalyzed_by")
    assert has("MN2m", "GENE:SOD2", "cofactor_of") and has("MN2c", "GENE:ARG1", "cofactor_of")
    assert relations == ["catalyzed_by", "cofactor_of"]
    assert summary["gene_nodes_added"] == ["SLC30A10"] and summary["cofactor_enzymes_absent_from_graph"] == ["ABSENT"]
    counts = new_edges.source_id.value_counts().add(new_edges.target_id.value_counts(), fill_value=0)
    assert (new_nodes.set_index("node_id").degree == counts.reindex(new_nodes.node_id).fillna(0).to_numpy()).all()
