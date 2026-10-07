"""Release regulation by auto- and heteroreceptors (graph/neurotransmission_regulation.py) and the curated oxidant
sources and targets (graph/oxidative_regulation.py)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.neurotransmission_regulation import add_release_regulation
from mechanistic_pathway_learning.graph.oxidative_regulation import add_oxidative_regulation, oxidant_bases, oxidant_node

NODE_COLUMNS = ["node_id", "node_type", "display_name", "compartment", "base_metabolite_id", "degree", "gene_symbol", "in_metabolic_layer", "is_currency"]
EDGE_COLUMNS = ["source_id", "target_id", "relation_type", "sign", "evidence_source"]
DOPAMINE, SEROTONIN, HYDROGEN_PEROXIDE, SUPEROXIDE = "MAM01736", "MAM02897", "MAM02041", "MAM02631"


def base_graph() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    rows = [{"node_id": f"{base}{letter}", "node_type": "metabolite", "display_name": name, "compartment": letter, "base_metabolite_id": base,
             "is_currency": True, "in_metabolic_layer": True, "degree": 0}
            for base, name in ((DOPAMINE, "dopamine"), (SEROTONIN, "serotonin"), (HYDROGEN_PEROXIDE, "H2O2"), (SUPEROXIDE, "O2-")) for letter in ("c", "e", "m")]
    rows += [{"node_id": reaction, "node_type": "reaction", "display_name": reaction, "compartment": "c;e", "is_currency": False, "degree": 0}
             for reaction in ("RCT_380869", "RCT_380901")]
    rows += [{"node_id": f"GENE:{gene}", "node_type": "gene", "display_name": gene, "gene_symbol": gene, "in_metabolic_layer": False, "degree": 0}
             for gene in ("DRD2", "HTR1A", "KEAP1", "NFE2L2", "SOD2", "MAOA", "NDUFV1")]
    nodes = pd.DataFrame(rows).reindex(columns=NODE_COLUMNS)
    edges = pd.DataFrame([{"source_id": "RCT_380869", "target_id": f"{DOPAMINE}e", "relation_type": "product_of", "sign": 1.0, "evidence_source": "Reactome"}],
                         columns=EDGE_COLUMNS)
    return nodes, edges, ["product_of"]


def receptor_table() -> pd.DataFrame:
    return pd.DataFrame([
        {"receptor_gene": "DRD2", "endogenous_ligand": "dopamine", "population": "dopaminergic", "role": "autoreceptor", "site": "presynaptic_terminal", "effect_on_release": -1, "pmid": 29030431},
        {"receptor_gene": "HTR1A", "endogenous_ligand": "serotonin", "population": "serotonergic", "role": "autoreceptor", "site": "somatodendritic", "effect_on_release": -1, "pmid": 11821026},
        {"receptor_gene": "DRD2", "endogenous_ligand": "dopamine", "population": "serotonergic", "role": "heteroreceptor", "site": "presynaptic_terminal", "effect_on_release": -1, "pmid": 18064418},
        {"receptor_gene": "HTR1A", "endogenous_ligand": "serotonin", "population": "histaminergic", "role": "heteroreceptor", "site": "somatodendritic", "effect_on_release": -1, "pmid": 11111111},
    ])


def test_receptor_activation_nodes_and_signed_regulation():
    nodes, edges, relation_types = base_graph()
    new_nodes, new_edges, new_relations, summary = add_release_regulation(
        nodes, edges, relation_types, receptor_table(), {DOPAMINE: ["RCT_380869"], SEROTONIN: ["RCT_380901"]})
    node_table = new_nodes.set_index("node_id")
    assert node_table.loc["RCPT_DRD2", "node_type"] == "protein_entity"
    assert node_table.loc["RCPT_DRD2", "display_name"] == "DRD2 activated by dopamine"
    # the gene and its ligand outside the cell both feed the activation node
    assert ((new_edges.source_id == "GENE:DRD2") & (new_edges.target_id == "RCPT_DRD2") & (new_edges.relation_type == "member_of")).any()
    ligand = new_edges[(new_edges.target_id == "RCPT_DRD2") & (new_edges.relation_type == "activates")]
    assert ligand.source_id.tolist() == [f"{DOPAMINE}e"]
    # the autoreceptor loop and the heteroreceptor on another population
    regulation = new_edges[new_edges.relation_type == "regulates_release"].set_index(["source_id", "target_id"]).sign
    assert regulation[("RCPT_DRD2", "RCT_380869")] == -1.0 and regulation[("RCPT_DRD2", "RCT_380901")] == -1.0
    assert regulation[("RCPT_HTR1A", "RCT_380901")] == -1.0
    assert summary["receptors"] == 2 and summary["regulation_edges"] == 3 and summary["inhibitory_edges"] == 3
    # the histaminergic row has no release reaction here and is reported, not invented
    assert summary["rows_without_a_release_reaction"] == ["HTR1A on histaminergic (serotonin)"]
    # a transmitter that signals is no longer currency
    assert not node_table.loc[f"{DOPAMINE}e", "is_currency"]
    assert "regulates_release" in new_relations


def test_conflicting_receptor_signs_are_left_out():
    table = pd.concat([receptor_table(), pd.DataFrame([
        {"receptor_gene": "DRD2", "endogenous_ligand": "dopamine", "population": "dopaminergic", "role": "heteroreceptor", "site": "somatodendritic", "effect_on_release": 1, "pmid": 1}])])
    nodes, edges, relation_types = base_graph()
    _, new_edges, _, summary = add_release_regulation(nodes, edges, relation_types, table, {DOPAMINE: ["RCT_380869"], SEROTONIN: ["RCT_380901"]})
    assert summary["conflicting_signs_left_out"] == ["DRD2 on dopaminergic"]
    assert not ((new_edges.relation_type == "regulates_release") & (new_edges.target_id == "RCT_380869")).any()
    assert ((new_edges.relation_type == "regulates_release") & (new_edges.source_id == "RCPT_HTR1A")).any()  # the other receptor is unaffected


def test_oxidant_species_and_compartment_resolution():
    assert oxidant_bases("hydrogen peroxide") == [HYDROGEN_PEROXIDE]
    assert sorted(oxidant_bases("superoxide and hydrogen peroxide")) == sorted([HYDROGEN_PEROXIDE, SUPEROXIDE])
    assert oxidant_bases("reactive oxygen species (unspecified oxidant)") == []
    assert oxidant_bases("lipid hydroperoxides") == []
    node_ids = {f"{HYDROGEN_PEROXIDE}{letter}" for letter in ("c", "m")}
    assert oxidant_node(HYDROGEN_PEROXIDE, "mitochondrial inner membrane, matrix", node_ids) == f"{HYDROGEN_PEROXIDE}m"
    assert oxidant_node(HYDROGEN_PEROXIDE, "plasma membrane", node_ids) == f"{HYDROGEN_PEROXIDE}c"
    assert oxidant_node(HYDROGEN_PEROXIDE, "peroxisome", node_ids) == f"{HYDROGEN_PEROXIDE}c"  # no peroxisomal copy here


def test_oxidant_target_and_source_edges():
    nodes, edges, relation_types = base_graph()
    targets = pd.DataFrame([
        {"species": "hydrogen peroxide", "target_gene": "KEAP1", "effect_sign": -1, "compartment": "cytosol", "pmid": 14585973},
        {"species": "hydrogen peroxide", "target_gene": "NFE2L2", "effect_sign": 1, "compartment": "cytosol", "pmid": 14585973},
        {"species": "peroxynitrite", "target_gene": "SOD2", "effect_sign": -1, "compartment": "mitochondrion", "pmid": 12791589},
        {"species": "reactive oxygen species (unspecified oxidant)", "target_gene": "KEAP1", "effect_sign": -1, "compartment": "cytosol", "pmid": 1},
        {"species": "hydrogen peroxide", "target_gene": "PTPN1", "effect_sign": -1, "compartment": "cytosol", "pmid": 9624118},
    ])
    sources = pd.DataFrame([
        {"source_gene": "MAOA", "species_produced": "hydrogen peroxide", "compartment": "mitochondrial outer membrane", "pmid": 17158340},
        {"source_gene": "NDUFV1", "species_produced": "superoxide", "compartment": "mitochondrial inner membrane, matrix", "pmid": 21393237},
        {"source_gene": "ALOX15", "species_produced": "lipid hydroperoxides", "compartment": "cytosol", "pmid": 14607519},
    ])
    new_nodes, new_edges, new_relations, summary = add_oxidative_regulation(nodes, edges, relation_types, targets, sources)
    modifications = new_edges[new_edges.relation_type == "oxidatively_modifies"].set_index(["source_id", "target_id"]).sign
    assert modifications[(f"{HYDROGEN_PEROXIDE}c", "GENE:KEAP1")] == -1.0
    assert modifications[(f"{HYDROGEN_PEROXIDE}c", "GENE:NFE2L2")] == 1.0
    assert ("MAM02714m", "GENE:SOD2") not in modifications.index  # peroxynitrite is not in this small graph
    production = new_edges[new_edges.relation_type == "produces_oxidant"].set_index(["source_id", "target_id"]).sign
    assert production[("GENE:MAOA", f"{HYDROGEN_PEROXIDE}m")] == 1.0
    assert production[("GENE:NDUFV1", f"{SUPEROXIDE}m")] == 1.0
    assert summary["genes_absent_from_the_graph"] == ["PTPN1"]
    assert "ALOX15 -> lipid hydroperoxides" in summary["species_without_a_human_gem_metabolite"]
    assert summary["inactivating_edges"] == 1 and summary["activating_edges"] == 1
    # the oxidants stop being currency, so a change in them propagates
    assert not new_nodes.set_index("node_id").loc[f"{HYDROGEN_PEROXIDE}c", "is_currency"]
    assert "oxidatively_modifies" in new_relations and "produces_oxidant" in new_relations


def test_curated_tables_in_the_repository_are_well_formed():
    receptors = pd.read_csv(Path("docs/curated_presynaptic_receptors.csv"))
    assert set(receptors.effect_on_release) <= {-1, 1}
    assert receptors.role.isin({"autoreceptor", "heteroreceptor"}).all()
    assert receptors.pmid.notna().all() and receptors.doi.notna().all() and receptors.supporting_quote.str.len().min() > 20
    assert not receptors.duplicated(["receptor_gene", "population", "role", "site"]).any()
    targets = pd.read_csv(Path("docs/curated_oxidant_targets.csv"))
    assert set(targets.effect_sign) <= {-1, 1} and targets.pmid.notna().all() and targets.doi.notna().all()
    sources = pd.read_csv(Path("docs/curated_oxidant_sources.csv"))
    assert sources.pmid.notna().all() and sources.doi.notna().all() and sources.compartment.notna().all()
