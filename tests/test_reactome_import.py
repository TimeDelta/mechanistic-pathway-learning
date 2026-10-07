"""The neuronal layer: compartment mapping, charge accounting across the plasma membrane, voltage gating, release
detection and the curated reactions (mechanistic_pathway_learning/graph/reactome_import.py)."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.reactome_import import (
    MEMBRANE_POTENTIAL_NODE,
    CuratedReaction,
    compartment_letter,
    import_reactome_layer,
    ions_moved_inward,
    moves_ions_across_a_membrane,
    neuronal_pool_id,
    parse_reactome_sbml,
    read_human_gem_participants,
    voltage_gating_sign,
)

SODIUM, POTASSIUM, CALCIUM, CHLORIDE, MAGNESIUM = "MAM02519", "MAM02200", "MAM01413", "MAM01442", "MAM02482"
ION_BASES = {"Na+": SODIUM, "K+": POTASSIUM, "Ca2+": CALCIUM, "Cl-": CHLORIDE, "Mg2+": MAGNESIUM}
ION_CHARGES = {SODIUM: 1, POTASSIUM: 1, CALCIUM: 2, CHLORIDE: -1, MAGNESIUM: 2}

REACTOME_SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
 <model id="test">
  <listOfCompartments>
   <compartment id="compartment_70101" name="cytosol"/>
   <compartment id="compartment_300" name="extracellular region"/>
   <compartment id="compartment_8" name="synaptic vesicle lumen"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="species_1000" name="Na+ [extracellular region]" compartment="compartment_300">
    <notes><p>Derived from a Reactome SimpleEntity</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:is xmlns:bqbiol="y">chebiId=CHEBI:29101</bqbiol:is></rdf:RDF></annotation>
   </species>
   <species id="species_1001" name="Na+ [cytosol]" compartment="compartment_70101">
    <notes><p>Derived from a Reactome SimpleEntity</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:is xmlns:bqbiol="y">chebiId=CHEBI:29101</bqbiol:is></rdf:RDF></annotation>
   </species>
   <species id="species_1002" name="dopamine [synaptic vesicle lumen]" compartment="compartment_8">
    <notes><p>Derived from a Reactome SimpleEntity</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:is xmlns:bqbiol="y">chebiId=CHEBI:18243</bqbiol:is></rdf:RDF></annotation>
   </species>
   <species id="species_1003" name="dopamine [extracellular region]" compartment="compartment_300">
    <notes><p>Derived from a Reactome SimpleEntity</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:is xmlns:bqbiol="y">chebiId=CHEBI:18243</bqbiol:is></rdf:RDF></annotation>
   </species>
   <species id="species_2000" name="Voltage gated Na+ channel [cytosol]" compartment="compartment_70101">
    <notes><p>Derived from a Reactome DefinedSet</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:hasPart xmlns:bqbiol="y">identifiers.org/uniprot:P35498</bqbiol:hasPart></rdf:RDF></annotation>
   </species>
   <species id="species_2001" name="Docked dopamine loaded synaptic vesicle [cytosol]" compartment="compartment_70101">
    <notes><p>Derived from a Reactome Complex</p></notes>
    <annotation><rdf:RDF xmlns:rdf="x"><bqbiol:hasPart xmlns:bqbiol="y">identifiers.org/uniprot:P21579</bqbiol:hasPart></rdf:RDF></annotation>
   </species>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="reaction_500" name="Na+ influx through voltage gated Na+ channels" reversible="false">
    <listOfReactants><speciesReference species="species_1000" stoichiometry="1" constant="true"/></listOfReactants>
    <listOfProducts><speciesReference species="species_1001" stoichiometry="1" constant="true"/></listOfProducts>
    <listOfModifiers><modifierSpeciesReference species="species_2000" sboTerm="SBO:0000013"/></listOfModifiers>
   </reaction>
   <reaction id="reaction_501" name="Release of docked dopamine loaded synaptic vesicle" reversible="false">
    <listOfReactants><speciesReference species="species_2001" stoichiometry="1" constant="true"/></listOfReactants>
    <listOfProducts><speciesReference species="species_1003" stoichiometry="1" constant="true"/></listOfProducts>
   </reaction>
  </listOfReactions>
 </model>
</sbml>
"""

HUMAN_GEM_SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml><model>
 <reaction id="R_MAR05429" name="sodium potassium pump" reversible="false">
  <listOfReactants>
   <speciesReference species="M_MAM02200e" stoichiometry="2" constant="true"/>
   <speciesReference species="M_MAM02519c" stoichiometry="3" constant="true"/>
  </listOfReactants>
  <listOfProducts>
   <speciesReference species="M_MAM02200c" stoichiometry="2" constant="true"/>
   <speciesReference species="M_MAM02519e" stoichiometry="3" constant="true"/>
  </listOfProducts>
 </reaction>
 <reaction id="R_MAR04956" name="sodium potassium chloride cotransport" reversible="false">
  <listOfReactants>
   <speciesReference species="M_MAM01442e" stoichiometry="2" constant="true"/>
   <speciesReference species="M_MAM02200e" stoichiometry="1" constant="true"/>
   <speciesReference species="M_MAM02519e" stoichiometry="1" constant="true"/>
  </listOfReactants>
  <listOfProducts>
   <speciesReference species="M_MAM01442c" stoichiometry="2" constant="true"/>
   <speciesReference species="M_MAM02200c" stoichiometry="1" constant="true"/>
   <speciesReference species="M_MAM02519c" stoichiometry="1" constant="true"/>
  </listOfProducts>
 </reaction>
 <reaction id="R_MAR00334" name="dopamine secretion via secretory vesicle" reversible="false">
  <listOfReactants><speciesReference species="M_MAM01736c" stoichiometry="1" constant="true"/></listOfReactants>
  <listOfProducts><speciesReference species="M_MAM01736e" stoichiometry="1" constant="true"/></listOfProducts>
 </reaction>
</model></sbml>
"""

NODE_COLUMNS = ["node_id", "node_type", "display_name", "compartment", "base_metabolite_id", "degree", "gene_symbol", "in_metabolic_layer", "is_transport",
                "reversible", "gene_reaction_rule", "subsystem", "is_currency"]
EDGE_COLUMNS = ["source_id", "target_id", "relation_type", "sign", "evidence_source"]


def base_graph() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """A graph with the ions in both compartments, cytosolic and extracellular dopamine, the pump, the cotransporter and
    Human-GEM's vesicular shortcut catalysed by VMAT2."""
    rows = []
    for base, name in ((SODIUM, "Na+"), (POTASSIUM, "K+"), (CALCIUM, "Ca2+"), (CHLORIDE, "Cl-"), ("MAM01736", "dopamine")):
        for letter in ("c", "e"):
            rows.append({"node_id": f"{base}{letter}", "node_type": "metabolite", "display_name": name, "compartment": letter, "base_metabolite_id": base,
                         "is_currency": base != "MAM01736", "in_metabolic_layer": True, "degree": 0})
    for reaction in ("MAR05429", "MAR04956", "MAR00334"):
        rows.append({"node_id": reaction, "node_type": "reaction", "display_name": reaction, "compartment": "c;e", "is_transport": True, "reversible": False,
                     "subsystem": "Transport reactions", "is_currency": False, "degree": 0})
    for gene in ("ATP1A3", "SLC12A2", "SLC18A2"):
        rows.append({"node_id": f"GENE:{gene}", "node_type": "gene", "display_name": gene, "gene_symbol": gene, "in_metabolic_layer": True, "degree": 0})
    nodes = pd.DataFrame(rows).reindex(columns=NODE_COLUMNS)
    edges = pd.DataFrame([{"source_id": "GENE:ATP1A3", "target_id": "MAR05429", "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM"},
                          {"source_id": "GENE:SLC12A2", "target_id": "MAR04956", "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM"},
                          {"source_id": "GENE:SLC18A2", "target_id": "MAR00334", "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM"}],
                         columns=EDGE_COLUMNS)
    return nodes, edges, ["catalyzed_by", "substrate_of", "product_of"]


def run_import(nodes=None, edges=None, **overrides):
    default_nodes, default_edges, relation_types = base_graph()
    nodes = default_nodes if nodes is None else nodes
    edges = default_edges if edges is None else edges
    arguments = dict(sbml_texts=[REACTOME_SBML], nodes=nodes, edges=edges, relation_types=relation_types,
                     gene_of_uniprot={"P35498": "SCN1A", "P21579": "SYT1"}, map_chebi_to_human_gem={"CHEBI:29101": SODIUM, "CHEBI:18243": "MAM01736"}.get,
                     ion_bases=ION_BASES, ion_charges=ION_CHARGES, human_gem_participants=read_human_gem_participants(HUMAN_GEM_SBML, set(ION_BASES.values())),
                     brain_expressed_genes={"ATP1A3", "SLC12A2", "SLC18A2"}, transmitter_bases={"MAM01736"}, vesicular_transporter_genes={"SLC18A2"})
    arguments.update(overrides)
    return import_reactome_layer(**arguments)


def test_compartment_letter_maps_vesicle_lumens_and_membranes():
    assert compartment_letter("synaptic vesicle lumen") == "v"
    assert compartment_letter("clathrin sculpted monoamine transport vesicle lumen") == "v"
    assert compartment_letter("plasma membrane") == "c"  # a membrane faces the cytosol
    assert compartment_letter("synaptic vesicle membrane") == "c"
    assert compartment_letter("extracellular region") == "e"
    assert compartment_letter("mitochondrial matrix") == "m"
    assert compartment_letter("mitochondrial intermembrane space") == "i"


def test_charge_counted_only_when_an_ion_crosses_the_plasma_membrane():
    sodium_entry = [("in", SODIUM, "e", 1.0), ("out", SODIUM, "c", 1.0)]
    assert ions_moved_inward(sodium_entry, set(ION_BASES.values())) == {SODIUM: 1.0}
    # the sodium-potassium pump: 3 Na+ out and 2 K+ in, one net positive charge leaving
    pump = [("in", POTASSIUM, "e", 2.0), ("in", SODIUM, "c", 3.0), ("out", POTASSIUM, "c", 2.0), ("out", SODIUM, "e", 3.0)]
    moved = ions_moved_inward(pump, set(ION_BASES.values()))
    assert moved == {SODIUM: -3.0, POTASSIUM: 2.0}
    assert sum(ION_CHARGES[base] * amount for base, amount in moved.items()) == -1
    # NKCC1: Na+, K+ and 2 Cl- inward, electroneutral
    cotransport = [("in", base, "e", amount) for base, amount in ((SODIUM, 1.0), (POTASSIUM, 1.0), (CHLORIDE, 2.0))]
    cotransport += [("out", base, "c", amount) for base, amount in ((SODIUM, 1.0), (POTASSIUM, 1.0), (CHLORIDE, 2.0))]
    moved = ions_moved_inward(cotransport, set(ION_BASES.values()))
    assert sum(ION_CHARGES[base] * amount for base, amount in moved.items()) == 0
    # receptor binding of a charged ligand outside the cell, and chemistry inside it, move nothing across
    assert ions_moved_inward([("in", CALCIUM, "e", 1.0), ("out", "RCTE_1", "c", 1.0)], set(ION_BASES.values())) == {}
    assert ions_moved_inward([("in", "MAM01371", "c", 1.0), ("out", SODIUM, "c", 1.0)], set(ION_BASES.values())) == {}


def test_voltage_gating_sign_by_channel_family():
    assert voltage_gating_sign({"SCN1A", "SCN2A"}) == 1.0
    assert voltage_gating_sign({"CACNA1A", "CACNB1"}) == 1.0
    assert voltage_gating_sign({"KCNA1", "KCNQ2", "KCNAB1"}) == 1.0
    assert voltage_gating_sign({"GRIN1", "GRIN2B"}) == 1.0  # the Mg2+ block is relieved by depolarisation
    assert voltage_gating_sign({"HCN1", "HCN4"}) == -1.0
    assert voltage_gating_sign({"KCNJ2", "KCNK3", "KCNN1"}) == 0.0  # inward rectifier, two-pore and SK: not voltage-gated
    assert voltage_gating_sign({"GABRA1", "SLC18A2"}) == 0.0
    # a reaction catalysed by a pore but not moving ions is not gated
    assert moves_ions_across_a_membrane([("in", SODIUM, "e", 1.0), ("out", SODIUM, "c", 1.0)], set(ION_BASES.values()))
    assert not moves_ions_across_a_membrane([("in", SODIUM, "c", 1.0), ("out", "MAM01371", "c", 1.0)], set(ION_BASES.values()))


def test_import_adds_pools_membrane_potential_and_gating():
    nodes, edges, relation_types, summary = run_import()
    assert summary["reactome_reactions"] == 2
    node_table = nodes.set_index("node_id")
    assert node_table.loc[MEMBRANE_POTENTIAL_NODE, "node_type"] == "membrane_potential"
    sodium_pool = neuronal_pool_id(SODIUM, "c")
    assert node_table.loc[sodium_pool, "is_currency"] is False or not node_table.loc[sodium_pool, "is_currency"]
    # the Reactome sodium entry: extracellular currency sodium in, the neuronal pool out, depolarising, Nav-gated
    assert ((edges.source_id == f"{SODIUM}e") & (edges.target_id == "RCT_500") & (edges.relation_type == "substrate_of")).any()
    assert ((edges.source_id == "RCT_500") & (edges.target_id == sodium_pool) & (edges.relation_type == "product_of")).any()
    depolarising = edges[(edges.source_id == "RCT_500") & (edges.target_id == MEMBRANE_POTENTIAL_NODE)]
    assert depolarising.relation_type.tolist() == ["changes_membrane_potential"] and depolarising.sign.tolist() == [1.0]
    gating = edges[(edges.source_id == MEMBRANE_POTENTIAL_NODE) & (edges.target_id == "RCT_500")]
    assert gating.relation_type.tolist() == ["voltage_gates"] and gating.sign.tolist() == [1.0]
    # the driving force of each pool on the potential
    driving = edges[edges.relation_type == "driving_force"].set_index("source_id").sign
    assert driving[neuronal_pool_id(SODIUM, "c")] == -1.0 and driving[neuronal_pool_id(POTASSIUM, "e")] == 1.0 and driving[neuronal_pool_id(CHLORIDE, "c")] == 1.0
    assert neuronal_pool_id(CALCIUM, "c") not in driving.index  # a signal, not a driving force
    assert "changes_membrane_potential" in relation_types and "driving_force" in relation_types


def test_human_gem_pump_hyperpolarises_and_cotransport_gives_no_edge():
    _, edges, _, summary = run_import()
    pump = edges[(edges.source_id == "MAR05429") & (edges.relation_type == "changes_membrane_potential")]
    assert pump.sign.tolist() == [-1.0]
    assert not ((edges.source_id == "MAR04956") & (edges.relation_type == "changes_membrane_potential")).any()
    # the pump moves sodium out of the cytosol and potassium into it
    pool_edges = edges[(edges.source_id == "MAR05429") & (edges.relation_type == "changes_ion_pool")].set_index("target_id").sign
    assert pool_edges[neuronal_pool_id(SODIUM, "c")] == -1.0
    assert pool_edges[neuronal_pool_id(POTASSIUM, "c")] == 1.0 and pool_edges[neuronal_pool_id(POTASSIUM, "e")] == -1.0
    assert summary["membrane_potential_edges"]["Human-GEM"] == 1


def test_a_reaction_whose_catalyst_is_not_brain_expressed_carries_no_current():
    _, edges, _, _ = run_import(brain_expressed_genes=set())
    assert not (edges.relation_type == "changes_membrane_potential").any() or not ((edges.source_id == "MAR05429")).any()


def test_release_takes_the_vesicle_lumen_transmitter_and_calcium():
    _, edges, _, summary = run_import()
    assert summary["release_reactions"] == {"RCT_501": ["MAM01736"]}
    # Reactome forms the loaded vesicle without consuming the lumen transmitter, so the lumen copy is added
    assert ((edges.source_id == "MAM01736v") & (edges.target_id == "RCT_501") & (edges.relation_type == "substrate_of")).any()
    calcium = edges[(edges.source_id == neuronal_pool_id(CALCIUM, "c")) & (edges.target_id == "RCT_501")]
    assert calcium.relation_type.tolist() == ["activates"] and calcium.sign.tolist() == [1.0]


def test_vesicular_shortcut_reaction_is_deleted_outright():
    nodes, edges, _, summary = run_import()
    # dropping the catalysis alone would leave an uncatalysed cytosol-to-extracellular route that still skips the vesicle
    assert "MAR00334" not in set(nodes.node_id)
    assert not ((edges.source_id == "MAR00334") | (edges.target_id == "MAR00334")).any()
    assert summary["nodes_deleted"] == 1 and len(summary["vesicular_shortcut_reactions_deleted"]) == 1
    assert "SLC18A2" in summary["vesicular_shortcut_reactions_deleted"][0]
    # the transmitter keeps no route out of the cytosol other than the vesicle cycle
    assert not ((edges.source_id == "MAM01736c") & edges.relation_type.eq("substrate_of")
                & edges.target_id.isin(edges[edges.target_id.eq("MAM01736e")].source_id)).any()
    # a non-vesicular transporter keeps its reaction
    assert ((edges.source_id == "GENE:ATP1A3") & (edges.target_id == "MAR05429")).any()


def test_shortcut_shared_with_another_transporter_keeps_the_reaction():
    """A reaction a non-vesicular transporter also carries loses only the vesicular catalysis."""
    nodes, edges, relation_types = base_graph()
    edges = pd.concat([edges, pd.DataFrame([{"source_id": "GENE:ATP1A3", "target_id": "MAR00334", "relation_type": "catalyzed_by", "sign": 1.0,
                                             "evidence_source": "Human-GEM"}], columns=EDGE_COLUMNS)], ignore_index=True)
    node_table, edge_table, _, summary = run_import(nodes=nodes, edges=edges)
    assert "MAR00334" in set(node_table.node_id) and summary["nodes_deleted"] == 0
    assert not ((edge_table.source_id == "GENE:SLC18A2") & (edge_table.target_id == "MAR00334")).any()
    assert ((edge_table.source_id == "GENE:ATP1A3") & (edge_table.target_id == "MAR00334")).any()
    assert len(summary["vesicular_shortcut_catalysis_removed"]) == 1


def test_curated_reaction_adds_its_missing_metabolite_copies():
    curated = [CuratedReaction("MAR_HISTAMINE_RELEASE", "release of histamine", ("MAM01736v",), ("MAM01736e",), ("SYT1",), "test", release="MAM01736")]
    nodes, edges, _, summary = run_import(curated_reactions=curated)
    node_table = nodes.set_index("node_id")
    assert node_table.loc["MAM01736v", "node_type"] == "metabolite" and node_table.loc["MAM01736v", "compartment"] == "v"
    assert node_table.loc["MAM01736v", "base_metabolite_id"] == "MAM01736"
    assert ((edges.source_id == "GENE:SYT1") & (edges.target_id == "MAR_HISTAMINE_RELEASE") & (edges.relation_type == "catalyzed_by")).any()
    assert "MAR_HISTAMINE_RELEASE" in summary["release_reactions"]


def test_every_edge_points_at_a_node_and_degrees_are_recomputed():
    nodes, edges, _, _ = run_import()
    assert not (set(edges.source_id) | set(edges.target_id)) - set(nodes.node_id)
    assert not nodes.node_id.duplicated().any()
    counts = edges.source_id.value_counts().add(edges.target_id.value_counts(), fill_value=0)
    assert nodes.set_index("node_id").degree[MEMBRANE_POTENTIAL_NODE] == counts[MEMBRANE_POTENTIAL_NODE]


def test_parse_reactome_sbml_reads_participants_and_modifier_roles():
    species, reactions = parse_reactome_sbml(REACTOME_SBML)
    assert species["species_1002"].compartment == "synaptic vesicle lumen" and species["species_1002"].chebi_ids == ["CHEBI:18243"]
    assert species["species_2000"].uniprot_ids == ["P35498"]
    entry = next(reaction for reaction in reactions if reaction.reaction_id == "500")
    assert entry.reactants == [("species_1000", 1.0)] and entry.catalysts == ["species_2000"] and not entry.inhibitors


def test_read_human_gem_participants_keeps_stoichiometry():
    participants = read_human_gem_participants(HUMAN_GEM_SBML, {SODIUM, POTASSIUM, CHLORIDE})
    assert ("in", SODIUM, "c", 3.0) in participants["MAR05429"]
    assert ("out", SODIUM, "e", 3.0) in participants["MAR05429"]
    assert "MAR00334" not in participants  # no ion among its participants


@pytest.mark.parametrize("variant_directory", ["data/processed/graph_neuronal"])
def test_built_variant_is_consistent_when_present(variant_directory):
    from pathlib import Path

    directory = Path(variant_directory)
    if not (directory / "nodes.parquet").exists():
        pytest.skip(f"{variant_directory} has not been built")
    nodes = pd.read_parquet(directory / "nodes.parquet")
    edges = pd.read_parquet(directory / "edges.parquet")
    summary = json.loads((directory / "neuronal_variant_summary.json").read_text())
    assert not (set(edges.source_id) | set(edges.target_id)) - set(nodes.node_id)
    assert (nodes.node_type == "membrane_potential").sum() == 1
    assert len(summary["release_reactions"]) >= 6
    assert summary["membrane_potential_edges"]["Reactome"] > 0 and summary["voltage_gating_edges"]["Reactome"] > 0


def test_curated_vesicle_chain_and_reverse_transport_in_the_built_variant():
    """In the built variant a transmitter reaches the cleft by exocytosis only through its vesicle lumen, and reverse
    transport through the plasma-membrane transporter is the one route that bypasses the vesicle."""
    from pathlib import Path

    directory = Path("data/processed/graph_neuronal")
    if not (directory / "nodes.parquet").exists():
        pytest.skip("the neuronal variant has not been built")
    nodes = pd.read_parquet(directory / "nodes.parquet")
    edges = pd.read_parquet(directory / "edges.parquet")
    summary = json.loads((directory / "neuronal_variant_summary.json").read_text())
    catalysts = edges[edges.relation_type == "catalyzed_by"]
    vesicular = {"SLC18A1", "SLC18A2", "SLC18A3", "SLC32A1", "SLC17A6", "SLC17A7", "SLC17A8"}
    compartment_of = dict(zip(nodes.node_id, nodes.compartment.astype(str)))
    for gene in sorted(vesicular):
        for reaction in catalysts[catalysts.source_id == f"GENE:{gene}"].target_id:
            letters = set(compartment_of.get(reaction, "").split(";"))
            assert not {"c", "e"} <= letters, f"{gene} still catalyses the cytosol-to-extracellular reaction {reaction}"
    # VMAT1 is neuroendocrine, so the curated brain loading reactions are VMAT2 alone (PMID 8643547)
    loading = [reaction for reaction in catalysts[catalysts.source_id == "GENE:SLC18A2"].target_id if reaction.startswith("MAR_")]
    assert loading and not [reaction for reaction in catalysts[catalysts.source_id == "GENE:SLC18A1"].target_id if reaction.startswith("MAR_")]
    # Release draws on the vesicle lumen copy wherever the loading reaction makes one (dopamine, serotonin, GABA,
    # histamine, glycine, noradrenaline); for glutamate and acetylcholine Reactome's loading consumes the cytosolic
    # transmitter straight into the loaded-vesicle complex, so the chain runs through that entity instead.
    node_ids = set(nodes.node_id)
    for reaction, released in summary["release_reactions"].items():
        for base in released:
            if f"{base}v" in node_ids:
                assert ((edges.source_id == f"{base}v") & (edges.target_id == reaction) & (edges.relation_type == "substrate_of")).any()
            else:
                loading = set(edges[(edges.source_id == f"{base}c") & (edges.relation_type == "substrate_of")].target_id)
                complexes = set(edges[edges.source_id.isin(loading) & (edges.relation_type == "product_of")].target_id)
                assert any(node.startswith("RCTE_") for node in complexes), f"{base} reaches release through no loaded vesicle"
    assert ((edges.source_id == "GENE:SLC6A3") & (edges.target_id == "MAR_DOPAMINE_EFFLUX")).any()


def test_a_lumped_human_gem_reaction_is_replaced_by_the_curated_families():
    """Human-GEM's MAR01527 lumps the Nav, HCN, NALCN, ENaC and TPCN families into one sodium reaction, so it can carry
    only one gating sign. The curated split gives each family its own, and the lumped reaction gets no edges at all."""
    curated = [CuratedReaction("MAR_NAV_SODIUM", "sodium entry through voltage-gated sodium channels", (f"{SODIUM}e",), (f"{SODIUM}c",), ("SCN1A",), "test"),
               CuratedReaction("MAR_HCN_SODIUM", "sodium entry through HCN channels", (f"{SODIUM}e",), (f"{SODIUM}c",), ("HCN1",), "test")]
    nodes, edges, relation_types = base_graph()
    nodes = pd.concat([nodes, pd.DataFrame([{"node_id": "MAR01527", "node_type": "reaction", "display_name": "Sodium Transport (Uniport)", "compartment": "c;e",
                                             "is_transport": True, "reversible": True, "subsystem": "Transport reactions", "is_currency": False, "degree": 0}],
                                           columns=NODE_COLUMNS)], ignore_index=True)
    edges = pd.concat([edges, pd.DataFrame([{"source_id": f"{SODIUM}e", "target_id": "MAR01527", "relation_type": "substrate_of", "sign": 1.0, "evidence_source": "Human-GEM"},
                                            {"source_id": "MAR01527", "target_id": f"{SODIUM}c", "relation_type": "product_of", "sign": 1.0, "evidence_source": "Human-GEM"},
                                            {"source_id": "GENE:ATP1A3", "target_id": "MAR01527", "relation_type": "catalyzed_by", "sign": 1.0, "evidence_source": "Human-GEM"}],
                                           columns=EDGE_COLUMNS)], ignore_index=True)
    node_table, edge_table, _, summary = run_import(nodes=nodes, edges=edges, curated_reactions=curated,
                                                    replaced_human_gem_reactions={"MAR01527": "lumps families whose gating differs"})
    assert "MAR01527" not in set(node_table.node_id)
    assert not ((edge_table.source_id == "MAR01527") | (edge_table.target_id == "MAR01527")).any()
    assert "MAR01527" in summary["human_gem_reactions_replaced_by_curated_ones"]
    gating = edge_table[edge_table.relation_type == "voltage_gates"].set_index("target_id").sign
    assert gating["MAR_NAV_SODIUM"] == 1.0 and gating["MAR_HCN_SODIUM"] == -1.0
    # both carry sodium inward, so both depolarise, however they are gated
    depolarising = edge_table[edge_table.relation_type == "changes_membrane_potential"].set_index("source_id").sign
    assert depolarising["MAR_NAV_SODIUM"] == 1.0 and depolarising["MAR_HCN_SODIUM"] == 1.0
    assert summary["voltage_gating_edges"].get("curated") == 2


def test_a_mixed_permeability_channel_keeps_its_net_current():
    """HCN passes sodium in and potassium out through one pore. Counted one for one the charges cancel and the
    current would vanish, so its net sign at rest is given explicitly; the pools still see both ions move."""
    sodium_pool, potassium_pool_in, potassium_pool_out = neuronal_pool_id(SODIUM, "c"), neuronal_pool_id(POTASSIUM, "c"), neuronal_pool_id(POTASSIUM, "e")
    mixed = CuratedReaction("MAR_HCN_CURRENT", "mixed current", (f"{SODIUM}e", potassium_pool_in), (sodium_pool, potassium_pool_out), ("HCN1",), "test",
                            membrane_potential_sign=1.0)
    counted = CuratedReaction("MAR_COUNTED", "the same ions without a sign", (f"{SODIUM}e", potassium_pool_in), (sodium_pool, potassium_pool_out), ("ATP1A3",), "test")
    _, edges, _, _ = run_import(curated_reactions=[mixed, counted])
    current = edges[(edges.relation_type == "changes_membrane_potential")].set_index("source_id").sign
    assert current["MAR_HCN_CURRENT"] == 1.0
    assert "MAR_COUNTED" not in current.index  # one sodium in, one potassium out: no net charge, no edge
    gating = edges[edges.relation_type == "voltage_gates"].set_index("target_id").sign
    assert gating["MAR_HCN_CURRENT"] == -1.0
    # the pools are reached through the reaction's own substrate and product edges, once each
    assert ((edges.source_id == "MAR_HCN_CURRENT") & (edges.target_id == sodium_pool) & (edges.relation_type == "product_of")).any()
    assert ((edges.source_id == potassium_pool_in) & (edges.target_id == "MAR_HCN_CURRENT") & (edges.relation_type == "substrate_of")).any()
    assert not ((edges.source_id == "MAR_HCN_CURRENT") & (edges.relation_type == "changes_ion_pool")).any()
