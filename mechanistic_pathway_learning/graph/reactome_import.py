"""Neuronal and oxidative layer from Reactome, with neuronal ion pools and a membrane potential node (graph variant for
the ablation of design section 5.7; experiments/build_neuronal_graph_variant.py).

Human-GEM is metabolism only: it has no vesicles, exocytosis, ion channels, receptors, G proteins or redox signalling,
and its inorganic ions are currency metabolites, so nothing in the graph carries membrane potential. This module
imports the reactions of chosen Reactome pathways (experiments/fetch_reactome_pathways.py) in Human-GEM's form and adds
the electrical layer:

  - Entities. A Reactome small molecule maps to a Human-GEM metabolite through the ChEBI bridge of the laboratory
    labels (conjugate forms, manual entries, named stereo forms), in the Human-GEM compartment that corresponds to its
    Reactome compartment (vesicle lumens become a vesicle compartment v); one Human-GEM lacks becomes a metabolite node
    CHEBI<number><compartment>. A protein, complex or set becomes an entity node RCTE_<id> of its own, so the states of
    one protein stay apart (free NRF2, KEAP1:NRF2, nuclear NRF2), and every member gene feeds it (member_of, +1);
    genes missing from the graph are added.
  - Reactions. Each becomes a reaction node RCT_<id>: inputs substrate_of (consumed, so the linear-response encoder
    derives depletion), outputs product_of, catalysts catalyzed_by, positive regulators activates (+1), negative
    regulators inhibits (-1). Curated reactions (vesicle loading Reactome lacks; experiments/build_neuronal_graph_variant.py)
    are added in the same form.
  - Neuronal ion pools. Human-GEM's Na+, K+, Ca2+ and chloride are currency metabolites: they receive but do not pass
    a signal on, since through them every transporter would reach every other in two steps. The neuron's own pools are
    separate nodes NEURON_<base><compartment> (cytosolic Na+, K+, chloride and Ca2+; extracellular K+). Reactome ion
    species in those compartments map to the pools; a Human-GEM reaction that moves the ion across the plasma membrane
    reaches the pool through changes_ion_pool (+1 into the pool, -1 out of it).
  - Membrane potential. A node VM_c (the neuronal membrane potential, linearised around rest) receives from every
    reaction that moves net charge across the plasma membrane (changes_membrane_potential), with the sign of the charge
    it brings in: a cation entering or an anion leaving depolarises (+1), a cation leaving or an anion entering
    hyperpolarises (-1), electroneutral transport (NKCC1: Na+, K+, 2 Cl-) gives no edge. Only inorganic ions that cross
    count: an ion consumed outside and produced in the cytosol (or the reverse), with the reaction's stoichiometry
    (Reactome's, and Human-GEM's from its SBML, so the Na+/K+-ATPase moves 3 Na+ out and 2 K+ in). Receptor binding,
    exocytosis and chemistry on one side of the membrane move no charge across it. A Human-GEM reaction counts only when
    a catalysing gene is expressed in brain (GTEx), so kidney and gut transporters do not drive the neuron.
  - Driving force. Each pool feeds VM_c with the sign of its effect on the ion's Nernst potential: more cytosolic Na+
    or K+ lowers E_Na or E_K (-1), more extracellular K+ raises E_K (+1, hyperkalaemic depolarisation), more cytosolic
    chloride raises E_Cl (+1, the depolarising GABA-A response when NKCC1 dominates KCC2).
  - Voltage gating. VM_c -> each reaction that moves ions across a membrane and is catalysed by a voltage-gated pore:
    Nav, Cav, Kv, BK and Hv1 open on depolarisation (+1), the NMDA receptor's Mg2+ block is relieved by it (+1), HCN
    channels open on hyperpolarisation (-1).
  - Release. A reaction that releases a transmitter (the transmitter appears outside from a loaded vesicle or from the
    vesicle lumen) is Ca2+-triggered: cytosolic Ca2+ pool -> release, activates +1, so depolarisation reaches release
    through the voltage-gated Ca2+ channels (VM_c -> Cav2 -> Ca2+ pool -> release). Reactome forms its loaded vesicles
    of dopamine, serotonin, noradrenaline and GABA without consuming the lumen transmitter its loading reactions make,
    so the lumen transmitter is made a substrate of release, which connects VMAT2 and VGAT to the cleft.
Human-GEM's vesicular shortcuts (VMAT, VAChT, VGAT and VGLUT moving a transmitter straight from the cytosol to the
extracellular space in one ATP-driven step, including Recon3D's 5-hydroxytryptophan secretion) are deleted outright,
with every edge they carry, when no other transporter also catalyses them. Dropping their catalysis alone is not
enough: an uncatalysed reaction still propagates substrate -> reaction -> product, so the transmitter would still reach
the cleft without the vesicle. Deleting them leaves the vesicle cycle as the only route for a transmitter to leave the
cytosol by exocytosis, which is what makes a vesicular transporter defect, or a drug acting on one, reach the synapse.
"""
from __future__ import annotations

import re
import xml.etree.ElementTree as ElementTree
from collections import defaultdict
from dataclasses import dataclass, field

import pandas as pd

SBML_NAMESPACE = {"sbml": "http://www.sbml.org/sbml/level3/version1/core"}
SMALL_MOLECULE_KINDS = {"SimpleEntity", "ChemicalDrug"}
CATALYST_SBO, STIMULATOR_SBO, INHIBITOR_SBO = "SBO:0000013", "SBO:0000459", "SBO:0000020"
MEMBRANE_POTENTIAL_NODE = "VM_c"
UNIPROT_PATTERN = re.compile(r"identifiers\.org/uniprot:([A-Z0-9]+)")
CHEBI_PATTERN = re.compile(r"chebiId=CHEBI:(\d+)")
KIND_PATTERN = re.compile(r"Derived from a Reactome (\w+)")
LOADED_VESICLE_PATTERN = re.compile(r"loaded.*vesicle", re.IGNORECASE)
HUMAN_GEM_REACTION_PATTERN = re.compile(r'<reaction [^>]*?id="R_(\w+)"(.*?)</reaction>', re.DOTALL)
HUMAN_GEM_REFERENCE_PATTERN = re.compile(r'<speciesReference species="M_(MAM\d+|\w+?)([a-z])" stoichiometry="([\d.eE+-]+)"')

# Pore-forming subunits whose opening depends on membrane potential, with the sign of that dependence: Nav
# (SCN1A-SCN11A), Cav (CACNA1A-CACNA1S), Kv (KCNA, KCNB, KCNC, KCND, KCNF, KCNG, KCNH, KCNQ, KCNS, KCNV), BK
# (KCNMA1), Hv1 (HVCN1) and the NMDA receptor (GRIN1, GRIN2A-D, whose Mg2+ block depolarisation relieves) open on
# depolarisation; HCN1-4 open on hyperpolarisation. Two-pore (K2P), inward-rectifier (Kir), SK and ligand-gated
# channels are not voltage-gated in this sense, although Guide to Pharmacology files K2P and Kir under its
# voltage-gated-like superfamily.
DEPOLARISATION_ACTIVATED_PATTERN = re.compile(r"^(SCN(1|2|3|4|5|8|9|10|11)A|CACNA1[A-IS]|KCN[ABCDFGHQSV]\d+|KCNMA1|HVCN1|GRIN1|GRIN2[A-D])$")
HYPERPOLARISATION_ACTIVATED_PATTERN = re.compile(r"^HCN[1-4]$")

# Effect of a rise in each neuronal pool on the membrane potential through the ion's Nernst potential (cations:
# E = (RT/zF) ln([out]/[in]), so a cytosolic rise lowers it and an extracellular rise raises it; chloride the
# reverse). Cytosolic Ca2+ is a signal (release, Ca2+-activated channels) rather than a driving force here: E_Ca lies so
# far from rest that its relative change is small.
NEURONAL_POOL_DRIVING_FORCE = {("Na+", "c"): -1.0, ("K+", "c"): -1.0, ("K+", "e"): 1.0, ("Cl-", "c"): 1.0, ("Ca2+", "c"): 0.0}
COMPARTMENT_WORDS = {"c": "cytosol", "e": "extracellular"}


@dataclass
class ReactomeSpecies:
    species_id: str
    name: str
    compartment: str
    kind: str
    chebi_ids: list[str] = field(default_factory=list)  # identity (bqbiol:is)
    uniprot_ids: list[str] = field(default_factory=list)  # identity or members (bqbiol:is, bqbiol:hasPart)


@dataclass
class ReactomeReaction:
    reaction_id: str
    name: str
    reactants: list[tuple[str, float]]
    products: list[tuple[str, float]]
    catalysts: list[str]
    stimulators: list[str]
    inhibitors: list[str]


@dataclass
class CuratedReaction:
    """A reaction Reactome lacks, in graph node ids: substrates and products are metabolite nodes, catalysts gene
    symbols. release names the transmitter base a vesicle-to-cleft release reaction releases.
    membrane_potential_sign overrides the sign derived from the charge that crosses, for a current whose direction is
    known while its stoichiometry is not: a channel permeable to two ions at once carries a net current set by the
    permeability ratio and the driving forces, not by a whole-number stoichiometry, so deriving the sign from counted
    charge would either invent that ratio or cancel to zero and drop the current altogether.
    catalysis_relation names the relation of the gene-to-reaction edge, for a reaction whose catalysis the model should
    weigh apart from ordinary catalysis. Gains are learned per relation, so a transporter that this set gives a
    reaction in each direction has one shared gain over both unless the reverse direction carries a relation of its
    own, and the model then cannot weigh efflux differently from uptake. Only the catalysis edge is separated: the
    substrate and product edges of an efflux reaction carry the same mass action as any other transport."""
    reaction_id: str
    display_name: str
    substrates: tuple[str, ...]
    products: tuple[str, ...]
    catalysts: tuple[str, ...]
    evidence: str
    release: str | None = None
    membrane_potential_sign: float | None = None
    catalysis_relation: str = "catalyzed_by"


def compartment_letter(name: str) -> str:
    """Human-GEM compartment letter for a Reactome compartment name; vesicle lumens become v."""
    lowered = name.lower()
    if "vesicle" in lowered and "membrane" not in lowered:
        return "v"
    for words, letter in ((("extracellular",), "e"), (("intermembrane",), "i"), (("mitochondri",), "m"), (("endoplasmic reticulum", "sarcoplasmic"), "r"),
                          (("golgi",), "g"), (("lysosom", "endosom"), "l"), (("peroxisom",), "x"), (("nucle",), "n")):
        if any(word in lowered for word in words):
            return letter
    return "c"  # cytosol, cytoplasm, plasma membrane, cell junction and other membranes facing the cytosol


def parse_reactome_sbml(sbml_text: str) -> tuple[dict[str, ReactomeSpecies], list[ReactomeReaction]]:
    root = ElementTree.fromstring(sbml_text)
    model = root.find("sbml:model", SBML_NAMESPACE)
    compartment_names = {element.get("id"): element.get("name") for element in model.iter(f"{{{SBML_NAMESPACE['sbml']}}}compartment")}
    species = {}
    for element in model.iter(f"{{{SBML_NAMESPACE['sbml']}}}species"):
        raw = ElementTree.tostring(element, encoding="unicode")
        kind = KIND_PATTERN.search(raw)
        identity = raw.split("hasPart")[0]
        species[element.get("id")] = ReactomeSpecies(
            species_id=element.get("id"), name=re.sub(r"\s*\[[^\]]*\]$", "", element.get("name", "")),
            compartment=compartment_names.get(element.get("compartment"), "cytosol"), kind=kind.group(1) if kind else "Unknown",
            chebi_ids=[f"CHEBI:{number}" for number in CHEBI_PATTERN.findall(identity)], uniprot_ids=sorted(set(UNIPROT_PATTERN.findall(raw))))
    reactions = []
    for element in model.iter(f"{{{SBML_NAMESPACE['sbml']}}}reaction"):
        def references(list_name: str) -> list[tuple[str, float]]:
            container = element.find(f"sbml:{list_name}", SBML_NAMESPACE)
            if container is None:
                return []
            return [(reference.get("species"), float(reference.get("stoichiometry", 1))) for reference in container]
        modifiers = references("listOfModifiers")
        modifier_terms = {}
        container = element.find("sbml:listOfModifiers", SBML_NAMESPACE)
        if container is not None:
            modifier_terms = {reference.get("species"): reference.get("sboTerm") for reference in container}
        reactions.append(ReactomeReaction(
            reaction_id=element.get("id").replace("reaction_", ""), name=element.get("name", ""), reactants=references("listOfReactants"), products=references("listOfProducts"),
            catalysts=[species_id for species_id, _ in modifiers if modifier_terms.get(species_id) == CATALYST_SBO],
            stimulators=[species_id for species_id, _ in modifiers if modifier_terms.get(species_id) == STIMULATOR_SBO],
            inhibitors=[species_id for species_id, _ in modifiers if modifier_terms.get(species_id) == INHIBITOR_SBO]))
    return species, reactions


def read_chebi_charges(chebi_obo_path) -> dict[str, int]:
    import gzip

    opener = gzip.open if str(chebi_obo_path).endswith(".gz") else open
    charges, current = {}, None
    with opener(chebi_obo_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("id: CHEBI:"):
                current = line[4:].strip()
            elif current and line.startswith("property_value: chemrof:charge"):
                match = re.search(r'"(-?\d+)"', line)
                if match:
                    charges[current] = int(match.group(1))
    return charges


def read_human_gem_participants(sbml_text: str, bases: set[str]) -> dict[str, list[tuple[str, str, str, float]]]:
    """{reaction id: [(role 'in' or 'out', base metabolite id, compartment letter, stoichiometry)]} from Human-GEM's
    SBML for the reactions with a participant among bases (the graph keeps no stoichiometric coefficients)."""
    participants = {}
    for match in HUMAN_GEM_REACTION_PATTERN.finditer(sbml_text):
        body = match.group(2)
        reactant_part, _, product_part = body.partition("listOfProducts")
        rows = [("in", base, letter, float(amount)) for base, letter, amount in HUMAN_GEM_REFERENCE_PATTERN.findall(reactant_part)]
        rows += [("out", base, letter, float(amount)) for base, letter, amount in HUMAN_GEM_REFERENCE_PATTERN.findall(product_part)]
        if any(base in bases for _, base, _, _ in rows):
            participants[match.group(1)] = rows
    return participants


def ions_moved_inward(participants: list[tuple[str, str, str, float]], ion_bases) -> dict[str, float]:
    """Net amount of each ion the reaction moves from the extracellular space (e) into the cytosol (c), negative for
    the reverse. participants: (role 'in' for a reactant or 'out' for a product, base metabolite id, compartment
    letter, stoichiometry). An ion counts only when it is consumed on one side of the plasma membrane and produced on
    the other, so binding, release and chemistry on one side move nothing."""
    moved = {}
    # sorted, because ion_bases arrives as a set and the callers add one edge per entry in this dictionary's order,
    # so an unsorted loop gave the ion-pool edges a row order that changed with PYTHONHASHSEED
    for base in sorted(ion_bases):
        def amount(role: str, letter: str) -> float:
            return sum(stoichiometry for row_role, row_base, row_letter, stoichiometry in participants if row_role == role and row_base == base and row_letter == letter)
        net = min(amount("in", "e"), amount("out", "c")) - min(amount("in", "c"), amount("out", "e"))
        if net:
            moved[base] = net
    return moved


def moves_ions_across_a_membrane(participants: list[tuple[str, str, str, float]], ion_bases) -> bool:
    """True when an ion is consumed in one compartment and produced in another (a channel or transporter reaction)."""
    for base in ion_bases:
        consumed = {letter for role, row_base, letter, _ in participants if role == "in" and row_base == base}
        produced = {letter for role, row_base, letter, _ in participants if role == "out" and row_base == base}
        if consumed and produced and consumed != produced:
            return True
    return False


def voltage_gating_sign(catalyst_genes: set[str]) -> float:
    """+1 when a depolarisation-activated pore catalyses the reaction, -1 when only hyperpolarisation-activated ones do,
    0 otherwise."""
    if any(DEPOLARISATION_ACTIVATED_PATTERN.match(gene) for gene in catalyst_genes):
        return 1.0
    if any(HYPERPOLARISATION_ACTIVATED_PATTERN.match(gene) for gene in catalyst_genes):
        return -1.0
    return 0.0


def neuronal_pool_id(base: str, letter: str) -> str:
    return f"NEURON_{base}{letter}"


def import_reactome_layer(sbml_texts: list[str], nodes: pd.DataFrame, edges: pd.DataFrame, relation_types: list[str], gene_of_uniprot: dict[str, str],
                          map_chebi_to_human_gem, ion_bases: dict[str, str], ion_charges: dict[str, int], human_gem_participants: dict[str, list],
                          brain_expressed_genes: set[str], transmitter_bases: set[str], vesicular_transporter_genes: set[str] = frozenset(),
                          curated_reactions: list[CuratedReaction] = (),
                          replaced_human_gem_reactions: dict[str, str] = ()) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict]:
    """The graph with the Reactome entities and reactions, the curated reactions, the neuronal ion pools, the membrane
    potential node and their edges.

    map_chebi_to_human_gem(chebi_id) -> Human-GEM base metabolite id or None. ion_bases: {'Na+', 'K+', 'Ca2+', 'Cl-',
    'Mg2+': Human-GEM base id}; ion_charges: {base id: charge}. human_gem_participants: read_human_gem_participants for
    the ion bases. brain_expressed_genes: gene symbols whose Human-GEM reactions may drive the neuron. transmitter_bases:
    base ids of the transmitters whose release is Ca2+-triggered. vesicular_transporter_genes: genes whose Human-GEM
    cytosol-to-extracellular shortcuts lose their catalysis. replaced_human_gem_reactions: {reaction id: why}, the
    Human-GEM reactions the curated set replaces, deleted with every edge they carry (a reaction that lumps channel
    families whose voltage dependence differs)."""
    species: dict[str, ReactomeSpecies] = {}
    reactions: dict[str, ReactomeReaction] = {}
    for text in sbml_texts:
        file_species, file_reactions = parse_reactome_sbml(text)
        species.update(file_species)
        reactions.update({reaction.reaction_id: reaction for reaction in file_reactions})
    existing = set(nodes.node_id)
    blank = {column: None for column in nodes.columns}
    new_nodes: dict[str, dict] = {}
    new_edges: list[dict] = []
    removed_edges: set[tuple[str, str, str]] = set()
    deleted_nodes: set[str] = {reaction for reaction in dict(replaced_human_gem_reactions) if reaction in set(nodes.node_id)}
    node_of_species: dict[str, str] = {}
    base_of_node: dict[str, str] = dict(zip(nodes.node_id, nodes.base_metabolite_id))
    letter_of_node: dict[str, str] = dict(zip(nodes.node_id, nodes.compartment.astype(str)))
    name_of_node: dict[str, str] = dict(zip(nodes.node_id, nodes.display_name.astype(str)))
    unmapped_small_molecules = set()
    ion_names = {base: name for name, base in ion_bases.items()}
    ion_set = set(ion_bases.values())
    pools = {(ion_bases[name], letter): sign for (name, letter), sign in NEURONAL_POOL_DRIVING_FORCE.items() if name in ion_bases}

    def add_node(node_id: str, **values) -> None:
        if node_id not in existing and node_id not in new_nodes:
            new_nodes[node_id] = {**blank, "node_id": node_id, "degree": 0, **values}
            base_of_node[node_id] = values.get("base_metabolite_id")
            letter_of_node[node_id] = values.get("compartment")
            name_of_node[node_id] = values.get("display_name")

    def add_edge(source: str, target: str, relation: str, sign: float, evidence: str) -> bool:
        """True when the edge was added; a reaction the curated set replaces gets no edges of its own."""
        if source in deleted_nodes or target in deleted_nodes:
            return False
        new_edges.append({"source_id": source, "target_id": target, "relation_type": relation, "sign": sign, "evidence_source": evidence})
        return True

    add_node(MEMBRANE_POTENTIAL_NODE, node_type="membrane_potential", display_name="neuronal membrane potential", compartment="c", is_currency=False)
    for (base, letter), driving_force_sign in pools.items():
        pool = neuronal_pool_id(base, letter)
        add_node(pool, node_type="metabolite", display_name=f"{ion_names[base]} (neuronal pool, {COMPARTMENT_WORDS[letter]})", compartment=letter,
                 base_metabolite_id=base, is_currency=False, in_metabolic_layer=False)
        if driving_force_sign:
            add_edge(pool, MEMBRANE_POTENTIAL_NODE, "driving_force", driving_force_sign, "Nernst potential of the ion")

    for species_id, entity in species.items():
        letter = compartment_letter(entity.compartment)
        if entity.kind in SMALL_MOLECULE_KINDS or (entity.chebi_ids and not entity.uniprot_ids):
            chebi_id = entity.chebi_ids[0] if entity.chebi_ids else None
            base = map_chebi_to_human_gem(chebi_id) if chebi_id else None
            if base is None:
                unmapped_small_molecules.add(f"{chebi_id} {entity.name}" if chebi_id else entity.name)
                base = (chebi_id or f"RCTSM{species_id.split('_')[-1]}").replace(":", "")
            node_id = neuronal_pool_id(base, letter) if (base, letter) in pools else f"{base}{letter}"
            add_node(node_id, node_type="metabolite", display_name=entity.name, compartment=letter, base_metabolite_id=base, is_currency=False, in_metabolic_layer=True)
            node_of_species[species_id] = node_id
        else:
            node_id = f"RCTE_{species_id.split('_')[-1]}"
            add_node(node_id, node_type="protein_entity", display_name=entity.name, compartment=letter, is_currency=False)
            node_of_species[species_id] = node_id
            for gene in sorted({gene_of_uniprot[accession] for accession in entity.uniprot_ids if accession in gene_of_uniprot}):
                add_node(f"GENE:{gene}", node_type="gene", display_name=gene, gene_symbol=gene, in_metabolic_layer=False)
                add_edge(f"GENE:{gene}", node_id, "member_of", 1.0, "Reactome entity member")

    genes_of_node: dict[str, set[str]] = defaultdict(set)
    for edge in new_edges:
        if edge["relation_type"] == "member_of":
            genes_of_node[edge["target_id"]].add(edge["source_id"][5:])

    # reaction node -> (participants, catalyst genes, reactant nodes, product nodes, source)
    reaction_records: dict[str, tuple[list, set[str], list[str], list[str], str]] = {}
    for reaction_id, reaction in reactions.items():
        node_id = f"RCT_{reaction_id}"
        reactant_nodes = [node_of_species[s] for s, _ in reaction.reactants if s in node_of_species]
        product_nodes = [node_of_species[s] for s, _ in reaction.products if s in node_of_species]
        letters = {letter_of_node.get(node) for node in reactant_nodes + product_nodes}
        add_node(node_id, node_type="reaction", display_name=reaction.name, compartment=";".join(sorted(letter for letter in letters if letter)),
                 is_transport=len(letters) > 1, reversible=False, subsystem="Reactome", is_currency=False)
        for species_id, _ in reaction.reactants:
            if species_id in node_of_species:
                add_edge(node_of_species[species_id], node_id, "substrate_of", 1.0, "Reactome")
        for species_id, _ in reaction.products:
            if species_id in node_of_species:
                add_edge(node_id, node_of_species[species_id], "product_of", 1.0, "Reactome")
        for relation, sign, species_ids in (("catalyzed_by", 1.0, reaction.catalysts), ("activates", 1.0, reaction.stimulators), ("inhibits", -1.0, reaction.inhibitors)):
            for species_id in species_ids:
                if species_id in node_of_species:
                    add_edge(node_of_species[species_id], node_id, relation, sign, "Reactome regulator")
        participants = [("in", base_of_node.get(node_of_species[s]), letter_of_node.get(node_of_species[s]), n) for s, n in reaction.reactants if s in node_of_species]
        participants += [("out", base_of_node.get(node_of_species[s]), letter_of_node.get(node_of_species[s]), n) for s, n in reaction.products if s in node_of_species]
        catalyst_genes = set().union(*(genes_of_node.get(node_of_species.get(s), set()) for s in reaction.catalysts)) if reaction.catalysts else set()
        reaction_records[node_id] = (participants, catalyst_genes, reactant_nodes, product_nodes, "Reactome")

    curated_membrane_signs: dict[str, float] = {}
    for curated in curated_reactions:
        for node in curated.substrates + curated.products:  # a compartment copy of a metabolite the base graph lacks (a transmitter in the vesicle lumen)
            if node not in existing and node not in new_nodes:
                add_node(node, node_type="metabolite", display_name=name_of_node.get(f"{node[:-1]}c", node[:-1]), compartment=node[-1],
                         base_metabolite_id=node[:-1], is_currency=False, in_metabolic_layer=True)
        add_node(curated.reaction_id, node_type="reaction", display_name=curated.display_name,
                 compartment=";".join(sorted({letter_of_node.get(node) for node in curated.substrates + curated.products} - {None})),
                 is_transport=len({letter_of_node.get(node) for node in curated.substrates + curated.products}) > 1, reversible=False, subsystem="Curated neurotransmission", is_currency=False)
        for node in curated.substrates:
            add_edge(node, curated.reaction_id, "substrate_of", 1.0, curated.evidence)
        for node in curated.products:
            add_edge(curated.reaction_id, node, "product_of", 1.0, curated.evidence)
        for gene in curated.catalysts:
            add_node(f"GENE:{gene}", node_type="gene", display_name=gene, gene_symbol=gene, in_metabolic_layer=False)
            add_edge(f"GENE:{gene}", curated.reaction_id, curated.catalysis_relation, 1.0, curated.evidence)
        participants = [("in", base_of_node.get(node), letter_of_node.get(node), 1.0) for node in curated.substrates]
        participants += [("out", base_of_node.get(node), letter_of_node.get(node), 1.0) for node in curated.products]
        reaction_records[curated.reaction_id] = (participants, set(curated.catalysts), list(curated.substrates), list(curated.products), "curated")
        if curated.membrane_potential_sign is not None:
            curated_membrane_signs[curated.reaction_id] = curated.membrane_potential_sign

    # Human-GEM: the vesicular shortcuts go, brain-expressed ion transport joins the electrical layer
    catalysts_of: dict[str, set[str]] = defaultdict(set)
    for source, target, relation in zip(edges.source_id, edges.target_id, edges.relation_type):
        if relation == "catalyzed_by" and str(source).startswith("GENE:"):
            catalysts_of[target].add(source[5:])
    shortcuts_deleted, shortcut_catalysis_removed = [], []
    for reaction_node, genes in catalysts_of.items():
        if not (genes & set(vesicular_transporter_genes) and {"c", "e"} <= set(str(letter_of_node.get(reaction_node, "")).split(";"))):
            continue
        if genes <= set(vesicular_transporter_genes):
            # the reaction exists only as the shortcut, so it goes with every edge it has: dropping the catalysis alone
            # would leave an uncatalysed route from the cytosol to the extracellular space, which still skips the vesicle
            deleted_nodes.add(reaction_node)
            shortcuts_deleted.append(f"{reaction_node} ({name_of_node.get(reaction_node)}), catalysed by {', '.join(sorted(genes))}")
        else:  # a reaction a non-vesicular transporter also carries: only the vesicular catalysis goes
            for gene in sorted(genes & set(vesicular_transporter_genes)):
                removed_edges.add((f"GENE:{gene}", reaction_node, "catalyzed_by"))
                shortcut_catalysis_removed.append(f"{gene} -> {reaction_node} ({name_of_node.get(reaction_node)})")
            catalysts_of[reaction_node] = genes - set(vesicular_transporter_genes)
    for reaction_id, participants in human_gem_participants.items():
        genes = catalysts_of.get(reaction_id, set())
        if reaction_id in letter_of_node and genes & brain_expressed_genes:
            reaction_records[reaction_id] = (participants, genes, [], [], "Human-GEM")

    membrane_edges, pool_edges, gating_edges = defaultdict(int), 0, defaultdict(int)
    for reaction_node, (participants, catalyst_genes, _, _, source) in reaction_records.items():
        moved = ions_moved_inward(participants, ion_set)
        inward_charge = sum(ion_charges.get(base, 0) * amount for base, amount in moved.items())
        curated_sign = curated_membrane_signs.get(reaction_node)
        sign = curated_sign if curated_sign is not None else (1.0 if inward_charge > 0 else -1.0)
        evidence = ("net current of a channel permeable to more than one ion (curated)" if curated_sign is not None
                    else f"net charge across the plasma membrane ({source})")
        if (curated_sign is not None or inward_charge) and add_edge(reaction_node, MEMBRANE_POTENTIAL_NODE, "changes_membrane_potential", sign, evidence):
            membrane_edges[source] += 1
        # Reactome and curated reactions reach the pools through their own substrate and product edges, since their
        # ion species are the pool nodes; only Human-GEM, whose ions are currency nodes, needs these edges
        if source == "Human-GEM":
            for base, amount in moved.items():
                for letter, direction in (("c", 1.0), ("e", -1.0)):
                    if (base, letter) in pools and add_edge(reaction_node, neuronal_pool_id(base, letter), "changes_ion_pool",
                                                            direction * (1.0 if amount > 0 else -1.0), f"ion moved across the plasma membrane ({source})"):
                        pool_edges += 1
        gating = voltage_gating_sign(catalyst_genes) if moves_ions_across_a_membrane(participants, ion_set | {base for base in ion_charges}) else 0.0
        if gating and add_edge(MEMBRANE_POTENTIAL_NODE, reaction_node, "voltage_gates", gating, "voltage-gated pore"):
            gating_edges[source] += 1

    calcium_pool = neuronal_pool_id(ion_bases["Ca2+"], "c") if "Ca2+" in ion_bases else None
    release_reactions = {}
    for reaction_node, (participants, _, reactant_nodes, product_nodes, source) in reaction_records.items():
        if source == "Human-GEM":
            continue
        released = {base_of_node.get(node) for node in product_nodes if letter_of_node.get(node) == "e" and base_of_node.get(node) in transmitter_bases}
        from_vesicle = any(LOADED_VESICLE_PATTERN.search(str(name_of_node.get(node, ""))) for node in reactant_nodes) or any(
            letter_of_node.get(node) == "v" and base_of_node.get(node) in transmitter_bases for node in reactant_nodes)
        curated_release = next((curated.release for curated in curated_reactions if curated.reaction_id == reaction_node), None)
        if curated_release:
            released, from_vesicle = {curated_release}, True
        if not (released and from_vesicle):
            continue
        release_reactions[reaction_node] = sorted(released)
        for base in released:
            lumen = f"{base}v"
            if (lumen in existing or lumen in new_nodes) and lumen not in reactant_nodes:
                add_edge(lumen, reaction_node, "substrate_of", 1.0, "vesicle lumen content released by exocytosis")
        if calcium_pool and calcium_pool not in reactant_nodes:
            add_edge(calcium_pool, reaction_node, "activates", 1.0, "Ca2+-triggered exocytosis")

    node_table = pd.concat([nodes[~nodes.node_id.isin(deleted_nodes)], pd.DataFrame(list(new_nodes.values()), columns=nodes.columns)], ignore_index=True)
    kept = (~pd.Series(list(zip(edges.source_id, edges.target_id, edges.relation_type))).isin(removed_edges).to_numpy()
            & ~edges.source_id.isin(deleted_nodes).to_numpy() & ~edges.target_id.isin(deleted_nodes).to_numpy())
    edge_table = pd.concat([edges[kept], pd.DataFrame(new_edges, columns=edges.columns)], ignore_index=True).astype(edges.dtypes.to_dict())
    edge_table = edge_table.drop_duplicates(["source_id", "target_id", "relation_type"], ignore_index=True)
    counts = edge_table.source_id.value_counts().add(edge_table.target_id.value_counts(), fill_value=0)
    node_table["degree"] = counts.reindex(node_table.node_id).fillna(0).astype("int64").to_numpy()
    # a curated catalysis relation may be the default "catalyzed_by", which the literal list already names, so the
    # relations are deduplicated rather than concatenated; dict.fromkeys keeps the order the encoder indexes them by
    added_relations = ("substrate_of", "product_of", "catalyzed_by", "activates", "inhibits", "member_of",
                       "changes_membrane_potential", "voltage_gates", "changes_ion_pool", "driving_force"
                       ) + tuple(sorted({curated.catalysis_relation for curated in curated_reactions}))
    new_relations = list(relation_types) + [relation for relation in dict.fromkeys(added_relations)
                                            if relation not in relation_types]
    summary = {"reactome_reactions": len(reactions), "reactome_species": len(species), "curated_reactions": len(curated_reactions),
               "nodes_added": len(node_table) - len(nodes) + len(deleted_nodes), "nodes_deleted": len(deleted_nodes),
               "edges_added": len(edge_table) - len(edges) + int((~kept).sum()), "edges_removed": int((~kept).sum()),
               "membrane_potential_edges": dict(membrane_edges), "ion_pool_edges_from_human_gem": pool_edges, "voltage_gating_edges": dict(gating_edges),
               "release_reactions": release_reactions, "vesicular_shortcut_reactions_deleted": shortcuts_deleted,
               "human_gem_reactions_replaced_by_curated_ones": {reaction: why for reaction, why in dict(replaced_human_gem_reactions).items() if reaction in deleted_nodes},
               "vesicular_shortcut_catalysis_removed": shortcut_catalysis_removed,
               "small_molecules_without_human_gem_metabolite": sorted(unmapped_small_molecules)}
    return node_table, edge_table, new_relations, summary
