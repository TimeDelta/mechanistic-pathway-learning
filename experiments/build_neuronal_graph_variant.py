"""Build the neuronal and oxidative graph variant on top of a graph directory, by default the manganese variant, into
<output-dir>: the Reactome layer with neuronal ion pools and a membrane potential node
(mechanistic_pathway_learning/graph/reactome_import.py), release regulation by auto- and heteroreceptors
(graph/neurotransmission_regulation.py) and curated oxidant sources and targets (graph/oxidative_regulation.py).

Inputs: the pinned Reactome exports (data/raw/reactome/v97, experiments/fetch_reactome_pathways.py), UniProt for
accession -> gene symbol, ChEBI for the ChEBI -> Human-GEM bridge of the laboratory labels, Human-GEM's SBML for ion
charges and transport stoichiometry, GTEx median TPM for the brain-expressed catalysts, and the curated tables
docs/curated_presynaptic_receptors.csv, docs/curated_oxidant_targets.csv and docs/curated_oxidant_sources.csv.

Usage:
  python experiments/build_neuronal_graph_variant.py --base-graph-dir data/processed/graph_manganese --output-dir data/processed/graph_neuronal
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.evidence.laboratory_abnormality_labels import (
    ChemicalDefinition,
    human_gem_metabolites_by_chebi,
    map_definition_to_metabolites,
    read_chebi_relations,
)
from mechanistic_pathway_learning.graph.build_physiology_graph import load_brain_expression
from mechanistic_pathway_learning.graph.neurotransmission_regulation import add_release_regulation
from mechanistic_pathway_learning.graph.node_descriptors import metabolite_charges_from_sbml
from mechanistic_pathway_learning.graph.brenda_enzyme_regulation import (
    DEFAULT_BRENDA_REGULATION_PATH,
    add_brenda_enzyme_regulation,
    load_brenda_regulation,
)
from mechanistic_pathway_learning.graph.oxidative_regulation import add_oxidative_regulation
from mechanistic_pathway_learning.graph.reactome_import import (
    read_chebi_charges,
    DEFAULT_MEMBRANE_POTENTIAL_BUILD,
    MEMBRANE_POTENTIAL_BUILDS,
    CuratedReaction,
    import_reactome_layer,
    read_human_gem_participants,
)
from mechanistic_pathway_learning.graph.redox_pools import add_redox_pools

# the inorganic ions whose movement across the plasma membrane carries current, with the Human-GEM base id of each
INORGANIC_IONS = {"Na+": "MAM02519", "K+": "MAM02200", "Ca2+": "MAM01413", "Cl-": "MAM01442", "Mg2+": "MAM02482"}
# Reactome's chloride is a second base id that no Human-GEM table holds, so it had no charge and was not counted as
# an ion moved (decision 9.4, the owner on 14 October 2026: "Definitely the second fix for chloride so it's separate
# from hcl"). MAM02144 is chloride: three of its four compartments are display-named "Cl-", and the cytosolic one
# reads "hydrochloride" only as a naming artefact - `RCT_975340` is "GABR heteropentamers:GABA transport **Cl-** from
# extracellular region to cytosol" and moves MAM02144e to MAM02144c. It is given chloride's charge and counted as an
# ion WITHOUT being merged into MAM01442 and without becoming a second chloride pool, which is the owner's choice of
# the two fixes: the two base ids stay separate entities and only the charge arithmetic is repaired.
# The charge is -1 per chloride and never -2 (the owner asked on 14 October 2026 whether it should be -1 rather than
# -2). Where -2 appears it is NKCC1's arithmetic, Na+ (+1) + K+ (+1) + 2 Cl- (2 x -1), which sums to 0 and so drops
# the current; the per-species charge below is -1. One limit now stated: Human-GEM annotates MAM02144 as
# `InChI=1S/ClH/h1H`, hydrogen chloride, which is formally neutral, so -1 is a judgement that the Reactome reactions
# use the entity as chloride and not a value read from a charge table. MAM01442, BiGG `cl`, is the entity that
# carries chloride's charge in Human-GEM, and the owner chose to keep the two separate rather than merge them.
IONS_WITHOUT_A_HUMAN_GEM_CHARGE = {"Cl- (Reactome)": ("MAM02144", -1)}
# The charge guard the owner approved on 14 October 2026 ("good call on prophylactically changing the mapper to
# ChEBl term Reactome names"). It turned out to be the substantive fix rather than a precaution, because NEITHER
# Human-GEM entity is chloride by its own annotation: MAM02144 is `InChI=1S/ClH/h1H`, hydrogen chloride, CHEBI:17883
# at charge 0, and MAM01442 is CHEBI:29311, "chlorine(.)", a radical, also at charge 0. Reactome names CHEBI:17996,
# chloride, at charge -1, and that is the charge the current needs. Every identifier and charge below was read from
# the pinned `data/raw/chebi/chebi_core.obo.gz` rather than recalled.
# Reading the charge from the term Reactome names, rather than from whichever Human-GEM entity the conjugate-acid
# step in `equivalent_forms` landed on, fixes the whole class: a mapping that crosses a conjugate acid or base
# changes the charge by one proton, and only chloride happens to be affected today because the other tracked ions
# are monatomic cations with no protonated form (decision 9.17).
ION_CHEBI_TERM_OF_BASE = {
    "MAM02519": "CHEBI:29101",  # sodium(1+)
    "MAM02200": "CHEBI:29103",  # potassium(1+)
    "MAM01413": "CHEBI:29108",  # calcium(2+)
    "MAM01442": "CHEBI:17996",  # chloride, where Human-GEM's own annotation is the chlorine radical
    "MAM02482": "CHEBI:18420",  # magnesium(2+)
    "MAM02144": "CHEBI:17996",  # Reactome's chloride, which Human-GEM holds only as hydrogen chloride
}
TRANSMITTER_BASES = {"MAM01736": "dopamine", "MAM02897": "serotonin", "MAM02617": "noradrenaline", "MAM01290": "adrenaline", "MAM01260": "acetylcholine",
                     "MAM02124": "histamine", "MAM01974": "glutamate", "MAM00970": "GABA", "MAM01986": "glycine"}
VESICULAR_TRANSPORTER_GENES = {"SLC18A1", "SLC18A2", "SLC18A3", "SLC32A1", "SLC17A6", "SLC17A7", "SLC17A8"}
# H2O2 everywhere (it diffuses and signals); glycine outside the cell and in vesicles (it is a transmitter)
SIGNALLING_METABOLITES = {"MAM02041": None, "MAM01986": {"e", "v"}}
MANGANESE_CHEBI = "CHEBI:29035"

# Vesicle loading, release and reverse transport Reactome's chosen pathways leave out, in Human-GEM's form.
# Loading is VMAT2 alone: central, peripheral and enteric neurons express only VMAT2, while VMAT1 (SLC18A1) is
# neuroendocrine, in chromaffin and enterochromaffin cells (Erickson et al. 1996, PMID 8643547). SLC18A1 keeps its
# gene node and its other reactions, since adrenal release reaches the brain through the circulation. Reactome loads dopamine,
# serotonin, GABA and acetylcholine into vesicles but not histamine, noradrenaline or glycine, and its loaded-vesicle
# complexes are formed without consuming the lumen transmitter, so release also takes the lumen copy as a substrate
# (reactome_import.import_reactome_layer). Histamine and noradrenaline are VMAT2 substrates (histamine with 30-fold
# higher affinity for VMAT2 than VMAT1; central, peripheral and enteric neurons express only VMAT2; Erickson et al.
# 1996, PMID 8643547). Glycine and GABA share the vesicular inhibitory amino acid transporter VGAT/SLC32A1, so
# glycinergic neurons have no separate vesicular glycine transporter (Wojcik et al. 2006, PMID 16701208).
# Noradrenaline is made inside the vesicle by dopamine beta-hydroxylase, a constituent of the catecholamine storage
# vesicles with its active site facing the lumen (Rush and Geffen 1980, PMID 6998654), so dopamine in the lumen, not
# noradrenaline in the cytosol, is its substrate. Release is Ca2+-triggered through synaptotagmin, with the Ca2+
# channels recruited to the active zone (Sudhof 2013, PMID 24183019).
CURATED_REACTIONS = [
    CuratedReaction("MAR_HISTAMINE_LOADING", "loading of histamine into secretory vesicles", ("MAM02124c",), ("MAM02124v",), ("SLC18A2",), "VMAT2 substrate (PMID 8643547)"),
    CuratedReaction("MAR_NORADRENALINE_LOADING", "loading of dopamine into noradrenergic vesicles", ("MAM01736c",), ("MAM01736v",), ("SLC18A2",), "VMAT2 substrate (PMID 8643547)"),
    # Human-GEM writes dopamine beta-monooxygenase (MAR06741) in the cytosol: ascorbate + dopamine + O2 ->
    # dehydroascorbate + H2O + noradrenaline. The same reaction with the lumen copies of the amines
    CuratedReaction("MAR_NORADRENALINE_VESICLE_SYNTHESIS", "dopamine beta-monooxygenase (vesicle lumen)", ("MAM01736v", "MAM02630c", "MAM01368c"),
                    ("MAM02617v", "MAM01655c", "MAM02040c"), ("DBH",), "dopamine beta-hydroxylase inside the storage vesicle (PMID 6998654)"),
    CuratedReaction("MAR_GLYCINE_LOADING", "loading of glycine into inhibitory vesicles", ("MAM01986c",), ("MAM01986v",), ("SLC32A1",), "shared vesicular inhibitory amino acid transporter (PMID 16701208)"),
    CuratedReaction("MAR_HISTAMINE_RELEASE", "release of histamine at the synapse", ("MAM02124v",), ("MAM02124e",), ("SYT1", "STX1A", "SNAP25", "VAMP2"),
                    "Ca2+-triggered exocytosis (PMID 24183019)", release="MAM02124"),
    CuratedReaction("MAR_GLYCINE_RELEASE", "release of glycine at the synapse", ("MAM01986v",), ("MAM01986e",), ("SYT1", "STX1A", "SNAP25", "VAMP2"),
                    "Ca2+-triggered exocytosis (PMID 24183019)", release="MAM01986"),
    # Reverse transport through the plasma-membrane transporter, the other half of the amphetamine mechanism: it raises
    # extracellular monoamine from the cytosolic pool without a vesicle and without depolarisation, while vesicular
    # redistribution empties the vesicle (Schmitz et al. 2001, PMID 11487614). Human-GEM writes these transporters as
    # uptake only, so without this reaction a drug acting on them can only lower extracellular transmitter
    CuratedReaction("MAR_DOPAMINE_EFFLUX", "reverse transport of dopamine (cytosol to extracellular)", ("MAM01736c",), ("MAM01736e",), ("SLC6A3",),
                    "reverse transport through the transporter (PMID 11487614)",
                    catalysis_relation="catalyzes_reverse_transport"),
    CuratedReaction("MAR_NORADRENALINE_EFFLUX", "reverse transport of noradrenaline (cytosol to extracellular)", ("MAM02617c",), ("MAM02617e",), ("SLC6A2",),
                    "reverse transport through the transporter (PMID 11487614)",
                    catalysis_relation="catalyzes_reverse_transport"),
    CuratedReaction("MAR_SEROTONIN_EFFLUX", "reverse transport of serotonin (cytosol to extracellular)", ("MAM02897c",), ("MAM02897e",), ("SLC6A4",),
                    "reverse transport through the transporter (PMID 11487614)",
                    catalysis_relation="catalyzes_reverse_transport"),
    # Human-GEM's MAR01527 "Sodium Transport (Uniport)" lumps every sodium-conducting channel into one reaction: the
    # voltage-gated Nav channels, the hyperpolarisation-activated HCN channels, the NALCN leak channel, the epithelial
    # ENaC subunits and the lysosomal two-pore channels. Their voltage dependence is opposite or absent, so one reaction
    # can carry only one gating sign, and the lumped reaction gave the HCN channels Nav's. Split by family: the Nav
    # channels carry the sodium current that generates the action potential (Catterall et al. 2020, PMID 33199904) and
    # open on depolarisation; the HCN channels carry I_h, open on hyperpolarisation and set the resting potential and
    # rhythmic firing (Biel et al. 2009, PMID 19584315), so their entry of sodium is depolarising while their gating is
    # not; NALCN and the atypical SCN7A are a leak; the ENaC subunits are epithelial and not voltage-gated. HCN conducts
    # potassium as well as sodium, which this reaction leaves out. TPCN1 and TPCN2 are lysosomal and keep Reactome's
    # lysosomal calcium reaction instead of a plasma-membrane sodium one
    # Products are the neuronal pool nodes, so these currents reach the pools as Reactome's ion species do.
    CuratedReaction("MAR_NAV_SODIUM", "sodium entry through voltage-gated sodium channels", ("MAM02519e",), ("NEURON_MAM02519c",),
                    ("SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A"), "voltage-gated sodium channel (PMID 33199904)"),
    # Calcium entry through the CaV1 (L-type) and CaV3 (T-type) channels, which no imported Reactome pathway and no
    # Human-GEM reaction carries: the graph's one voltage-gated calcium entry was Reactome's presynaptic RCT_210420,
    # catalysed by the CaV2 set (CACNA1A, CACNA1B, CACNA1E; decision 9.69). "The Ca(V)1 subfamily initiates contraction,
    # secretion, regulation of gene expression, integration of synaptic input in neurons, and synaptic transmission at
    # ribbon synapses in specialized sensory cells. ... The Ca(V)3 subfamily is important for repetitive firing of action
    # potentials in rhythmically firing cells such as cardiac myocytes and thalamic neurons" (Catterall 2011, PMID
    # 21746798). T-type channels "can be activated by small depolarizations of the plasma membrane", and "Ca2+ entry via
    # T-type channels can directly regulate intracellular Ca2+ concentrations" (Perez-Reyes 2003, PMID 12506128). Both
    # families open on depolarisation (DEPOLARISATION_ACTIVATED_PATTERN already names CACNA1A to CACNA1I and CACNA1S) and
    # carry Ca2+ inward, so they depolarise and raise the neuronal calcium pool, which drives release and the Ca2+-
    # activated channels. The auxiliary alpha2delta and beta subunits "enhance channel trafficking, increase the
    # expression of functional calcium channels at the plasma membrane and influence the channels' biophysical
    # properties" (Dolphin 2012, PMID 22805911), so they raise the CaV1 reaction's rate as activators rather than
    # catalysing it; the CaV3 alpha1 subunits are written without them. CACNA1F (retina) and CACNA1S (skeletal muscle)
    # are in the family as SCN4A and SCN5A are in the sodium one, and the cell-class weights give each class its share.
    CuratedReaction("MAR_CAV1_CALCIUM", "calcium entry through L-type voltage-gated calcium channels (CaV1)", ("MAM01413e",), ("NEURON_MAM01413c",),
                    ("CACNA1C", "CACNA1D", "CACNA1F", "CACNA1S"), "L-type voltage-gated calcium channel, CaV1 (PMID 21746798; auxiliary subunits PMID 22805911)",
                    activators=("CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNA2D4", "CACNB1", "CACNB2", "CACNB3", "CACNB4")),
    CuratedReaction("MAR_CAV3_CALCIUM", "calcium entry through T-type voltage-gated calcium channels (CaV3)", ("MAM01413e",), ("NEURON_MAM01413c",),
                    ("CACNA1G", "CACNA1H", "CACNA1I"), "T-type voltage-gated calcium channel, CaV3 (PMID 21746798; PMID 12506128)"),
    # Store-operated calcium entry (decision 9.71): "the identification of the STIM proteins as ER Ca(2+) sensors and the
    # Orai proteins as store-operated channels"; "Depletion of Ca(2+) from the ER causes STIM to accumulate at
    # ER-plasma membrane (PM) junctions where it traps and activates Orai channels" (Prakriya and Lewis 2015, PMID
    # 26400989). ORAI1 to ORAI3 are the pores and STIM1 and STIM2 their activators; ORAI is not voltage-gated, so the
    # reaction carries no gating edge, and its calcium entry depolarises and raises the neuronal pool.
    CuratedReaction("MAR_ORAI_CALCIUM", "store-operated calcium entry through ORAI channels", ("MAM01413e",), ("NEURON_MAM01413c",),
                    ("ORAI1", "ORAI2", "ORAI3"), "store-operated calcium channel, ORAI with STIM (PMID 26400989)", activators=("STIM1", "STIM2")),
    # P2X receptors (decision 9.71): "membrane ion channels that open in response to the binding of extracellular ATP",
    # and "All P2X receptors are permeable to small monovalent cations; some have significant calcium or anion
    # permeability" (North 2002, PMID 12270951). Written as the sodium entry every subtype carries, which depolarises and
    # raises the neuronal sodium pool, as the Nav and HCN families are; the calcium permeability of some subtypes is
    # left out, as HCN's potassium is. The ATP that opens them is a currency metabolite in the graph (it sends nothing),
    # so its gating cannot enter as an edge; a P2X gene's perturbation reaches the current through catalysis.
    CuratedReaction("MAR_P2X_CATION_ENTRY", "cation entry through ATP-gated P2X receptors", ("MAM02519e",), ("NEURON_MAM02519c",),
                    ("P2RX1", "P2RX2", "P2RX3", "P2RX4", "P2RX5", "P2RX6", "P2RX7"), "ATP-gated cation channel, P2X (PMID 12270951)"),
    # HCN: mixed sodium and potassium current; its net sign at rest is set explicitly (see CuratedReaction).
    CuratedReaction("MAR_HCN_CURRENT", "mixed sodium and potassium current through HCN channels", ("MAM02519e", "NEURON_MAM02200c"),
                    ("NEURON_MAM02519c", "NEURON_MAM02200e"), ("HCN1", "HCN2", "HCN3", "HCN4"),
                    "hyperpolarisation-activated cation channel, I_h (PMID 19584315)", membrane_potential_sign=1.0),
    CuratedReaction("MAR_SODIUM_LEAK", "sodium leak into the cytosol", ("MAM02519e",), ("NEURON_MAM02519c",), ("NALCN", "SCN7A"),
                    "sodium leak channel, not voltage-gated (PMID 33199904)"),
    CuratedReaction("MAR_ENAC_SODIUM", "sodium entry through epithelial sodium channels", ("MAM02519e",), ("NEURON_MAM02519c",),
                    ("SCNN1A", "SCNN1B", "SCNN1D", "SCNN1G"), "epithelial sodium channel, not voltage-gated (PMID 33199904)"),
]
# Human-GEM reactions the curated set replaces, deleted with every edge they carry
REPLACED_HUMAN_GEM_REACTIONS = {"MAR01527": "Sodium Transport (Uniport) lumps the Nav, HCN, NALCN, ENaC and TPCN families, whose voltage dependence differs; "
                                            "split into MAR_NAV_SODIUM, MAR_HCN_CURRENT, MAR_SODIUM_LEAK and MAR_ENAC_SODIUM"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base-graph-dir", type=Path, default=Path("data/processed/graph_manganese"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/graph_neuronal"))
    parser.add_argument("--reactome-dir", type=Path, default=Path("data/raw/reactome/v97"))
    parser.add_argument("--uniprot-table", type=Path, default=Path("data/raw/uniprot/uniprot_human_reviewed.tsv.gz"))
    parser.add_argument("--chebi-obo", type=Path, default=Path("data/raw/chebi/chebi_core.obo.gz"))
    parser.add_argument("--human-gem-metabolites", type=Path, default=Path("data/raw/Human-GEM/model/metabolites.tsv"))
    parser.add_argument("--human-gem-sbml", type=Path, default=Path("data/raw/Human-GEM/model/Human-GEM.xml"))
    parser.add_argument("--brain-expression", type=Path, default=Path("data/raw/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"))
    parser.add_argument("--minimum-brain-tpm", type=float, default=1.0)
    parser.add_argument("--receptor-table", type=Path, default=Path("docs/curated_presynaptic_receptors.csv"))
    parser.add_argument("--oxidant-target-table", type=Path, default=Path("docs/curated_oxidant_targets.csv"))
    parser.add_argument("--oxidant-source-table", type=Path, default=Path("docs/curated_oxidant_sources.csv"))
    parser.add_argument("--brenda-regulation-table", type=Path, default=DEFAULT_BRENDA_REGULATION_PATH,
                        help="BRENDA's human enzyme regulation as metabolite-to-enzyme pairs (decision 9.5): inhibits, activates and cofactor_of, "
                             "unweighted, CC BY 4.0. The graph held four metabolite-to-protein inhibits edges before this")
    parser.add_argument("--without-curated-layers", action="store_true",
                        help="leave out the curated receptor and oxidant edges, for the ablation that separates them from the Reactome layer")
    parser.add_argument("--membrane-potential", choices=MEMBRANE_POTENTIAL_BUILDS, default=DEFAULT_MEMBRANE_POTENTIAL_BUILD,
                        help="how VM_c gets its value. learned_from_reactions (the default, and the build registered up to 12 October 2026): one edge "
                             "from each of the 405 reactions that move net charge, which the encoder averages. computed_from_the_pools: the six ion "
                             "pools feed it with the weights of the linearised parallel-conductance equation and the 405 edges are not made "
                             "(mechanistic_pathway_learning/graph/membrane_potential.py, docs/membrane_potential_constants.md)")
    arguments = parser.parse_args()

    nodes = pd.read_parquet(arguments.base_graph_dir / "nodes.parquet")
    edges = pd.read_parquet(arguments.base_graph_dir / "edges.parquet")
    relation_types = json.loads((arguments.base_graph_dir / "relation_types.json").read_text())
    uniprot = pd.read_csv(arguments.uniprot_table, sep="\t")
    gene_of_uniprot = {accession: gene for accession, gene in zip(uniprot["Entry"], uniprot["Gene Names (primary)"]) if isinstance(gene, str)}
    names, neighbours, children = read_chebi_relations(arguments.chebi_obo)
    metabolites_by_chebi = human_gem_metabolites_by_chebi(pd.read_csv(arguments.human_gem_metabolites, sep="\t"))
    if "MN2" in set(nodes.base_metabolite_id.dropna()):
        metabolites_by_chebi[MANGANESE_CHEBI] = {"MN2"}  # Mn(2+) of the manganese variant

    def map_chebi_to_human_gem(chebi_id: str) -> str | None:
        mapped, _ = map_definition_to_metabolites(ChemicalDefinition("", (chebi_id,), 0, ""), names, neighbours, children, metabolites_by_chebi)
        return sorted(mapped)[0] if mapped else None

    human_gem_sbml = arguments.human_gem_sbml.read_text()
    ion_charges = {base: charge for base, charge in metabolite_charges_from_sbml(human_gem_sbml).items() if base in set(INORGANIC_IONS.values())}
    # the Reactome-only ions. They join ion_bases so `ions_moved_inward` counts them and ion_charges so the net charge
    # sees them; they are absent from NEURONAL_POOL_DRIVING_FORCE, so no pool node is built for them.
    ion_bases = dict(INORGANIC_IONS)
    for name, (base, charge) in IONS_WITHOUT_A_HUMAN_GEM_CHARGE.items():
        ion_bases[name] = base
        ion_charges[base] = charge
    # the charge guard: the term Reactome names wins over whatever entity the conjugate-acid step landed on
    chebi_charges = read_chebi_charges(arguments.chebi_obo)
    charges_taken_from_chebi = {}
    for base, chebi_term in ION_CHEBI_TERM_OF_BASE.items():
        if chebi_term in chebi_charges:
            charges_taken_from_chebi[base] = chebi_charges[chebi_term]
            ion_charges[base] = chebi_charges[chebi_term]
    print(f"ion charges read from the ChEBI term Reactome names: {charges_taken_from_chebi}")
    expression = load_brain_expression(arguments.brain_expression)
    brain_expressed_genes = set(expression.loc[expression.brain_median_tpm_max >= arguments.minimum_brain_tpm, "gene_symbol"])
    sbml_texts = [path.read_text() for path in sorted(arguments.reactome_dir.glob("*.sbml"))]
    new_nodes, new_edges, new_relations, summary = import_reactome_layer(
        sbml_texts, nodes, edges, relation_types, gene_of_uniprot, map_chebi_to_human_gem, ion_bases, ion_charges,
        read_human_gem_participants(human_gem_sbml, set(INORGANIC_IONS.values())), brain_expressed_genes, set(TRANSMITTER_BASES),
        VESICULAR_TRANSPORTER_GENES, CURATED_REACTIONS, REPLACED_HUMAN_GEM_REACTIONS, arguments.membrane_potential)

    release_reactions_of_transmitter: dict[str, list[str]] = {}
    for reaction, released in summary["release_reactions"].items():
        for base in released:
            release_reactions_of_transmitter.setdefault(base, []).append(reaction)
    receptor_summary = oxidative_summary = brenda_summary = "left out (--without-curated-layers)"
    if not arguments.without_curated_layers:
        new_nodes, new_edges, new_relations, receptor_summary = add_release_regulation(
            new_nodes, new_edges, new_relations, pd.read_csv(arguments.receptor_table), release_reactions_of_transmitter)
        new_nodes, new_edges, new_relations, oxidative_summary = add_oxidative_regulation(
            new_nodes, new_edges, new_relations, pd.read_csv(arguments.oxidant_target_table), pd.read_csv(arguments.oxidant_source_table))
        # before the gene/protein split the graph has gene nodes and no uniprot_entry column, so BRENDA's accessions
        # reach their gene nodes through the same UniProt accession -> gene name map the Reactome layer uses above
        new_nodes, new_edges, new_relations, brenda_summary = add_brenda_enzyme_regulation(
            new_nodes, new_edges, new_relations, load_brenda_regulation(arguments.brenda_regulation_table), gene_of_uniprot)
    new_nodes, new_edges, new_relations, redox_summary = add_redox_pools(new_nodes, new_edges, new_relations)
    for base, compartments in SIGNALLING_METABOLITES.items():
        rows = new_nodes.base_metabolite_id.eq(base) & (new_nodes.compartment.isin(compartments) if compartments else True)
        new_nodes.loc[rows, "is_currency"] = False

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    new_nodes.to_parquet(arguments.output_dir / "nodes.parquet")
    new_edges.to_parquet(arguments.output_dir / "edges.parquet")
    (arguments.output_dir / "relation_types.json").write_text(json.dumps(new_relations, indent=1) + "\n")
    summary.update(base_graph=str(arguments.base_graph_dir), reactome_pathways=json.loads((arguments.reactome_dir / "pathways.json").read_text()),
                   release_regulation=receptor_summary, oxidative_regulation=oxidative_summary, brenda_enzyme_regulation=brenda_summary,
                   redox_pools=redox_summary, minimum_brain_tpm=arguments.minimum_brain_tpm,
                   nodes_added_in_total=len(new_nodes) - len(nodes), edges_added_in_total=len(new_edges) - len(edges),
                   nodes_by_type=new_nodes.node_type.value_counts().to_dict(), edges_by_relation=new_edges.relation_type.value_counts().to_dict())
    (arguments.output_dir / "neuronal_variant_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in ("reactome_pathways", "small_molecules_without_human_gem_metabolite")}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
