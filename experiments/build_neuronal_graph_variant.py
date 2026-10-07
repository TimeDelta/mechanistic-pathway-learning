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
from mechanistic_pathway_learning.graph.oxidative_regulation import add_oxidative_regulation
from mechanistic_pathway_learning.graph.reactome_import import (
    CuratedReaction,
    import_reactome_layer,
    read_human_gem_participants,
)

# the inorganic ions whose movement across the plasma membrane carries current, with the Human-GEM base id of each
INORGANIC_IONS = {"Na+": "MAM02519", "K+": "MAM02200", "Ca2+": "MAM01413", "Cl-": "MAM01442", "Mg2+": "MAM02482"}
TRANSMITTER_BASES = {"MAM01736": "dopamine", "MAM02897": "serotonin", "MAM02617": "noradrenaline", "MAM01290": "adrenaline", "MAM01260": "acetylcholine",
                     "MAM02124": "histamine", "MAM01974": "glutamate", "MAM00970": "GABA", "MAM01986": "glycine"}
VESICULAR_TRANSPORTER_GENES = {"SLC18A1", "SLC18A2", "SLC18A3", "SLC32A1", "SLC17A6", "SLC17A7", "SLC17A8"}
# H2O2 everywhere (it diffuses and signals); glycine outside the cell and in vesicles (it is a transmitter)
SIGNALLING_METABOLITES = {"MAM02041": None, "MAM01986": {"e", "v"}}
MANGANESE_CHEBI = "CHEBI:29035"

# Vesicle loading and release Reactome's chosen pathways leave out, in Human-GEM's form. Reactome loads dopamine,
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
    CuratedReaction("MAR_HISTAMINE_LOADING", "loading of histamine into secretory vesicles", ("MAM02124c",), ("MAM02124v",), ("SLC18A2", "SLC18A1"), "VMAT2 substrate (PMID 8643547)"),
    CuratedReaction("MAR_NORADRENALINE_LOADING", "loading of dopamine into noradrenergic vesicles", ("MAM01736c",), ("MAM01736v",), ("SLC18A2", "SLC18A1"), "VMAT2 substrate (PMID 8643547)"),
    # Human-GEM writes dopamine beta-monooxygenase (MAR06741) in the cytosol: ascorbate + dopamine + O2 ->
    # dehydroascorbate + H2O + noradrenaline. The same reaction with the lumen copies of the amines
    CuratedReaction("MAR_NORADRENALINE_VESICLE_SYNTHESIS", "dopamine beta-monooxygenase (vesicle lumen)", ("MAM01736v", "MAM02630c", "MAM01368c"),
                    ("MAM02617v", "MAM01655c", "MAM02040c"), ("DBH",), "dopamine beta-hydroxylase inside the storage vesicle (PMID 6998654)"),
    CuratedReaction("MAR_GLYCINE_LOADING", "loading of glycine into inhibitory vesicles", ("MAM01986c",), ("MAM01986v",), ("SLC32A1",), "shared vesicular inhibitory amino acid transporter (PMID 16701208)"),
    CuratedReaction("MAR_HISTAMINE_RELEASE", "release of histamine at the synapse", ("MAM02124v",), ("MAM02124e",), ("SYT1", "STX1A", "SNAP25", "VAMP2"),
                    "Ca2+-triggered exocytosis (PMID 24183019)", release="MAM02124"),
    CuratedReaction("MAR_GLYCINE_RELEASE", "release of glycine at the synapse", ("MAM01986v",), ("MAM01986e",), ("SYT1", "STX1A", "SNAP25", "VAMP2"),
                    "Ca2+-triggered exocytosis (PMID 24183019)", release="MAM01986"),
]


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
    parser.add_argument("--without-curated-layers", action="store_true",
                        help="leave out the curated receptor and oxidant edges, for the ablation that separates them from the Reactome layer")
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
    expression = load_brain_expression(arguments.brain_expression)
    brain_expressed_genes = set(expression.loc[expression.brain_median_tpm_max >= arguments.minimum_brain_tpm, "gene_symbol"])
    sbml_texts = [path.read_text() for path in sorted(arguments.reactome_dir.glob("*.sbml"))]
    new_nodes, new_edges, new_relations, summary = import_reactome_layer(
        sbml_texts, nodes, edges, relation_types, gene_of_uniprot, map_chebi_to_human_gem, INORGANIC_IONS, ion_charges,
        read_human_gem_participants(human_gem_sbml, set(INORGANIC_IONS.values())), brain_expressed_genes, set(TRANSMITTER_BASES),
        VESICULAR_TRANSPORTER_GENES, CURATED_REACTIONS)

    release_reactions_of_transmitter: dict[str, list[str]] = {}
    for reaction, released in summary["release_reactions"].items():
        for base in released:
            release_reactions_of_transmitter.setdefault(base, []).append(reaction)
    receptor_summary = oxidative_summary = "left out (--without-curated-layers)"
    if not arguments.without_curated_layers:
        new_nodes, new_edges, new_relations, receptor_summary = add_release_regulation(
            new_nodes, new_edges, new_relations, pd.read_csv(arguments.receptor_table), release_reactions_of_transmitter)
        new_nodes, new_edges, new_relations, oxidative_summary = add_oxidative_regulation(
            new_nodes, new_edges, new_relations, pd.read_csv(arguments.oxidant_target_table), pd.read_csv(arguments.oxidant_source_table))
    for base, compartments in SIGNALLING_METABOLITES.items():
        rows = new_nodes.base_metabolite_id.eq(base) & (new_nodes.compartment.isin(compartments) if compartments else True)
        new_nodes.loc[rows, "is_currency"] = False

    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    new_nodes.to_parquet(arguments.output_dir / "nodes.parquet")
    new_edges.to_parquet(arguments.output_dir / "edges.parquet")
    (arguments.output_dir / "relation_types.json").write_text(json.dumps(new_relations, indent=1) + "\n")
    summary.update(base_graph=str(arguments.base_graph_dir), reactome_pathways=json.loads((arguments.reactome_dir / "pathways.json").read_text()),
                   release_regulation=receptor_summary, oxidative_regulation=oxidative_summary, minimum_brain_tpm=arguments.minimum_brain_tpm,
                   nodes_added_in_total=len(new_nodes) - len(nodes), edges_added_in_total=len(new_edges) - len(edges),
                   nodes_by_type=new_nodes.node_type.value_counts().to_dict(), edges_by_relation=new_edges.relation_type.value_counts().to_dict())
    (arguments.output_dir / "neuronal_variant_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in ("reactome_pathways", "small_molecules_without_human_gem_metabolite")}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
