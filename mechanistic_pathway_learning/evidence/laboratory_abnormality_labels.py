"""Laboratory-abnormality labels: signed metabolite changes of monogenic disorders, from the HPO annotations of a gene
and the logical definitions of the HPO terms, mapped onto Human-GEM metabolites (design section 4, auxiliary labels).

An HPO term such as Hyperphenylalaninemia (HP:0004923) is defined in hp-base.owl as
  has_part some ('increased amount' and inheres_in some (phenylalanine and part_of some blood) and has_modifier abnormal),
so it names a direction (PATO:0000470 increased amount, PATO:0001997 decreased amount; PATO:0000070 amount gives no
direction), a ChEBI entity and a body fluid (UBERON: blood 0000178, urine 0001088, cerebrospinal fluid 0001359 and
others). hp.obo carries none of these definitions, hence the OWL file.

HPO names the neutral form of a compound (phenylalanine, CHEBI:28044) where Human-GEM often lists a charged or
stereo-specific one (L-phenylalanine, its zwitterion), so ChEBI identifiers are bridged in two tiers:
  - equivalent: the same compound in another protonation or tautomeric state (ChEBI relations is conjugate acid of,
    is conjugate base of, is tautomer of, followed transitively);
  - specific form: a direct is_a child of the named compound, or of one of its equivalents, whose name contains the
    named compound's name (L-phenylalanine under phenylalanine), with its own equivalents; kept apart in the report
    because a class such as 'amino acid' would otherwise reach every amino acid. Where a compound has both an L- and a
    D- child (ornithine, arginine), only the L- form is kept: clinical amino acid measurements are of the L- form.
Where HPO and Human-GEM name a measured element in different forms, a manual entry of MANUAL_HUMAN_GEM_MAPPINGS maps
it (route 'manual'): Human-GEM annotates K+ (MAM02200) with CHEBI:26216, the potassium atom, while the HPO potassium
terms (hypokalemia, hyperkalemia) name potassium(1+); conversely several HPO terms name the element or a generic
cation (calcium atom in hypercalciuria, iron atom and iron cation, sodium, zinc, copper and magnesium atoms) where
Human-GEM has the ion, and HPO's 'phosphate' is the clinical inorganic phosphate, Human-GEM's Pi. Serum iron counts
both oxidation states, so iron maps to Fe2+ and Fe3+; Human-GEM has no Cu+, so copper maps to Cu2+. Human-GEM has no
manganese metabolite; the manganese terms map to MN2, the Mn2+ node of the manganese graph variant
(graph/manganese_extension.py), and reach no node in the base graphs. The same table maps compounds that HPO names in a generic,
neutral or racemic form the ontology bridge cannot reach (glucose, lactic acid, cholesterol, adrenaline, galactose,
argininosuccinate, DOPAC and others below, each to the Human-GEM metabolite of that name), and a few class names to
the members a clinical test measures: ketone bodies (acetoacetate, 3-hydroxybutyrate, acetone), catecholamines
(adrenaline, noradrenaline, dopamine), LDL-, HDL- and VLDL-cholesterol to Human-GEM's lipoprotein particles, very
long-chain fatty acids to hexacosanoate (C26:0, the diagnostic one) and homocystine to homocystine and homocysteine
(total homocysteine). Classes without one clinical meaning (amino acids, organic acids, bile acids, porphyrins,
acylcarnitines, glycosaminoglycans) and peptide hormones, which Human-GEM does not contain, stay unmapped.
A term that still reaches more than MAXIMUM_METABOLITES_PER_TERM metabolites names a class (prostaglandins reach 32,
fatty acids 11) and is left out, as are generic pseudo-metabolites such as '[protein]'.
A label is a (gene, metabolite, direction, fluid) record; an unannotated pair is unlabelled, not normal.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

INCREASED_AMOUNT = "PATO_0000470"
DECREASED_AMOUNT = "PATO_0001997"
FLUID_BY_UBERON = {"UBERON_0000178": "blood", "UBERON_0001088": "urine", "UBERON_0001359": "cerebrospinal_fluid",
                   "UBERON_0001977": "serum", "UBERON_0002107": "liver"}
MAXIMUM_METABOLITES_PER_TERM = 3
MANUAL_HUMAN_GEM_MAPPINGS: dict[str, frozenset[str]] = {
    "CHEBI:29103": frozenset({"MAM02200"}),  # potassium(1+) -> K+ (Human-GEM lists the potassium atom, CHEBI:26216)
    "CHEBI:22984": frozenset({"MAM01413"}),  # calcium atom -> Ca2+
    "CHEBI:26020": frozenset({"MAM02751"}),  # phosphate (clinical inorganic phosphate) -> Pi
    "CHEBI:18248": frozenset({"MAM01821", "MAM01822"}),  # iron atom -> Fe2+, Fe3+
    "CHEBI:24875": frozenset({"MAM01821", "MAM01822"}),  # iron cation -> Fe2+, Fe3+
    "CHEBI:26708": frozenset({"MAM02519"}),  # sodium atom -> Na+
    "CHEBI:27363": frozenset({"MAM03157"}),  # zinc atom -> zinc (annotated zinc(2+))
    "CHEBI:28694": frozenset({"MAM01624"}),  # copper atom -> Cu2+
    "CHEBI:23378": frozenset({"MAM01624"}),  # copper cation -> Cu2+
    "CHEBI:25107": frozenset({"MAM02482"}),  # magnesium atom -> Mg2+
    "CHEBI:25155": frozenset({"MN2"}),  # manganese cation -> Mn2+ of the manganese graph variant (graph/manganese_extension.py); no node in the base graphs
    # the same compound in a generic, neutral or racemic form
    "CHEBI:17234": frozenset({"MAM01965"}),  # glucose -> glucose (Human-GEM: D-glucopyranose)
    "CHEBI:28358": frozenset({"MAM02403"}),  # rac-lactic acid -> L-lactate (the clinical lactate)
    "CHEBI:16113": frozenset({"MAM01450"}),  # cholesterol -> cholesterol
    "CHEBI:33568": frozenset({"MAM01290"}),  # adrenaline -> adrenaline ((R)-adrenaline)
    "CHEBI:28260": frozenset({"MAM01910"}),  # galactose -> galactose
    "CHEBI:17754": frozenset({"MAM01983"}),  # glycerol -> glycerol
    "CHEBI:20106": frozenset({"MAM03136"}),  # vanillylmandelic acid -> vanillylmandelate
    "CHEBI:18240": frozenset({"MAM03037"}),  # 4-hydroxy-L-proline -> trans-4-hydroxy-L-proline
    "CHEBI:24741": frozenset({"MAM03037"}),  # hydroxyproline -> trans-4-hydroxy-L-proline
    "CHEBI:28115": frozenset({"MAM03322"}),  # methylcobalamin -> methylcobalamin
    "CHEBI:61409": frozenset({"MAM00729"}),  # dihydroxyphenylacetic acid -> 3,4-dihydroxyphenylacetate (DOPAC)
    "CHEBI:17012": frozenset({"MAM02543"}),  # N-acetylneuraminic acid -> N-acetylneuraminate
    "CHEBI:89843": frozenset({"MAM03765"}),  # methylsuccinic acid -> methyl-succinate
    "CHEBI:184023": frozenset({"MAM01366"}),  # argininosuccinic acid -> argininosuccinate
    "CHEBI:46819": frozenset({"MAM03120"}),  # urate salt -> urate
    "CHEBI:26361": frozenset({"MAM02803"}),  # protoporphyrins -> protoporphyrin (IX)
    "CHEBI:12777": frozenset({"MAM02834", "MAM20001"}),  # vitamin A -> retinol (both Human-GEM entries)
    "CHEBI:33234": frozenset({"MAM01327"}),  # vitamin E -> alpha-tocopherol (the measured form)
    "CHEBI:23641": frozenset({"MAM01673"}),  # deoxyuridine phosphate, the entity of HP:0034277 'elevated circulating deoxyuridine' -> deoxyuridine
    "CHEBI:32797": frozenset({"MAM00653"}),  # (S)-2-hydroxyglutaric acid -> 2-hydroxyglutarate (the L2HGDH substrate; the R form is MAM20012)
    # a class name to the members a clinical test measures
    "CHEBI:17087": frozenset({"MAM01253", "MAM00157", "MAM01256"}),  # ketone -> ketone bodies: acetoacetate, (R)-3-hydroxybutanoate, acetone
    "CHEBI:33567": frozenset({"MAM01290", "MAM02617", "MAM01736"}),  # catecholamine -> adrenaline, noradrenaline, dopamine
    "CHEBI:47774": frozenset({"MAM03710"}),  # LDL cholesterol -> low density lipoprotein
    "CHEBI:47775": frozenset({"MAM03647"}),  # HDL cholesterol -> high density lipoprotein
    "CHEBI:47773": frozenset({"MAM04074"}),  # VLDL cholesterol -> very low density lipoprotein
    "CHEBI:27283": frozenset({"MAM01432"}),  # very long-chain fatty acid -> hexacosanoate (cerotic acid, C26:0)
    "CHEBI:17485": frozenset({"MAM03394", "MAM02133"}),  # homocystine -> L-homocystine and homocysteine (total homocysteine)
}
CHEBI_EQUIVALENCE_RELATIONS = {"RO:0018033", "RO:0018034", "RO:0018036"}  # is conjugate base of, is conjugate acid of, is tautomer of
HP_CLASS_PATTERN = re.compile(r'<owl:Class rdf:about="http://purl.obolibrary.org/obo/(HP_\d+)">')
EQUIVALENT_CLASS_PATTERN = re.compile(r"<owl:equivalentClass>(.*?)</owl:equivalentClass>", flags=re.S)


@dataclass(frozen=True)
class ChemicalDefinition:
    hpo_id: str  # HP:0004923
    chebi_ids: tuple[str, ...]  # CHEBI:28044
    direction: int  # +1 increased, -1 decreased, 0 abnormal amount without a direction
    fluid: str  # blood, urine, cerebrospinal_fluid, serum, liver, other or unspecified


def parse_hpo_chemical_definitions(owl_text: str) -> list[ChemicalDefinition]:
    """Every HPO class whose equivalence axiom names a ChEBI entity, with its direction and fluid."""
    starts = [(match.start(), match.group(1)) for match in HP_CLASS_PATTERN.finditer(owl_text)]
    definitions = []
    for index, (start, hpo_identifier) in enumerate(starts):
        end = starts[index + 1][0] if index + 1 < len(starts) else len(owl_text)
        equivalent = EQUIVALENT_CLASS_PATTERN.search(owl_text, start, end)
        if not equivalent:
            continue
        axiom = equivalent.group(1)
        chebi_ids = tuple(sorted({f"CHEBI:{number}" for number in re.findall(r"CHEBI_(\d+)", axiom)}))
        if not chebi_ids:
            continue
        direction = 1 if INCREASED_AMOUNT in axiom else -1 if DECREASED_AMOUNT in axiom else 0
        fluids = sorted({FLUID_BY_UBERON.get(uberon, "other") for uberon in re.findall(r"(UBERON_\d+)", axiom)})
        definitions.append(ChemicalDefinition(hpo_identifier.replace("_", ":"), chebi_ids, direction, ";".join(fluids) or "unspecified"))
    return definitions


def read_chebi_relations(chebi_obo_path: Path) -> tuple[dict[str, str], dict[str, set[str]], dict[str, set[str]]]:
    """(name by ChEBI id, equivalence neighbours by id, direct is_a children by id) from chebi_core.obo(.gz);
    alt_ids resolve to their primary id."""
    opener = gzip.open if str(chebi_obo_path).endswith(".gz") else open
    names: dict[str, str] = {}
    neighbours: dict[str, set[str]] = defaultdict(set)
    children: dict[str, set[str]] = defaultdict(set)
    alternative_of: dict[str, str] = {}
    current = None
    with opener(chebi_obo_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "[Term]" or line == "[Typedef]":
                current = None
            elif line.startswith("id: CHEBI:"):
                current = line[4:]
            elif current is None:
                continue
            elif line.startswith("name: "):
                names[current] = line[6:]
            elif line.startswith("alt_id: "):
                alternative_of[line[8:]] = current
            elif line.startswith("is_a: CHEBI:"):
                children[line[6:].split()[0]].add(current)
            elif line.startswith("relationship: "):
                relation, target = line.split()[1:3]
                if relation in CHEBI_EQUIVALENCE_RELATIONS:
                    neighbours[current].add(target)
                    neighbours[target].add(current)
    for alternative, primary in alternative_of.items():
        names.setdefault(alternative, names.get(primary, ""))
        neighbours[alternative].add(primary)
        neighbours[primary].add(alternative)
    return names, neighbours, children


def equivalent_forms(chebi_id: str, neighbours: dict[str, set[str]]) -> set[str]:
    """The compound in all its protonation and tautomeric states (transitive closure of the equivalence relations)."""
    seen, stack = {chebi_id}, [chebi_id]
    while stack:
        for neighbour in neighbours.get(stack.pop(), ()):
            if neighbour not in seen:
                seen.add(neighbour)
                stack.append(neighbour)
    return seen


def specific_forms(chebi_id: str, names: dict[str, str], neighbours: dict[str, set[str]], children: dict[str, set[str]]) -> set[str]:
    """Direct is_a children of the compound or its equivalents whose name contains the compound's name, with their
    own equivalent forms (L-phenylalanine and its zwitterion under phenylalanine)."""
    base_name = names.get(chebi_id, "").lower()
    if len(base_name) < 3:
        return set()
    named_children = {child for equivalent in equivalent_forms(chebi_id, neighbours) for child in children.get(equivalent, ())
                      if base_name in names.get(child, "").lower()}
    child_names = {names.get(child, "") for child in named_children}
    forms: set[str] = set()
    for child in named_children:
        child_name = names.get(child, "")
        if child_name.startswith("D-") and f"L-{child_name[2:]}" in child_names:
            continue  # both stereo forms are children: the measured one is the L- form
        forms |= equivalent_forms(child, neighbours)
    return forms


def human_gem_metabolites_by_chebi(metabolites_table: pd.DataFrame) -> dict[str, set[str]]:
    """ChEBI id -> Human-GEM base metabolite ids (metsNoComp), from the metChEBIID column (several ids separated by ';')."""
    mapping: dict[str, set[str]] = defaultdict(set)
    for base_id, chebi_field in zip(metabolites_table.metsNoComp, metabolites_table.metChEBIID):
        if pd.isna(chebi_field):
            continue
        for chebi_id in str(chebi_field).split(";"):
            chebi_id = chebi_id.strip()
            if chebi_id:
                mapping[chebi_id if chebi_id.startswith("CHEBI:") else f"CHEBI:{chebi_id}"].add(base_id)
    return mapping


def generic_metabolite_ids(nodes: pd.DataFrame) -> set[str]:
    """Base ids of Human-GEM pseudo-metabolites that stand for a class of molecules ('[protein]', '[acceptor]')."""
    metabolites = nodes[(nodes.node_type == "metabolite") & nodes.display_name.notna()]
    return set(metabolites.loc[metabolites.display_name.astype(str).str.startswith("["), "base_metabolite_id"].dropna())


def map_definition_to_metabolites(definition: ChemicalDefinition, names, neighbours, children, metabolites_by_chebi,
                                  excluded_metabolites: frozenset[str] = frozenset(),
                                  maximum_metabolites: int = MAXIMUM_METABOLITES_PER_TERM) -> tuple[set[str], str]:
    """(Human-GEM base metabolites, route): route is 'equivalent', 'manual' (MANUAL_HUMAN_GEM_MAPPINGS),
    'specific_form', 'class_term' (more than maximum_metabolites reached; no metabolites returned) or 'unmapped'."""
    def mapped_through(chebi_ids: set[str]) -> set[str]:
        return set().union(*(metabolites_by_chebi.get(chebi_id, set()) for chebi_id in chebi_ids)) - excluded_metabolites

    equivalent_ids = set().union(*(equivalent_forms(chebi_id, neighbours) for chebi_id in definition.chebi_ids))
    mapped = mapped_through(equivalent_ids)
    route = "equivalent"
    if not mapped:
        mapped = set().union(*(MANUAL_HUMAN_GEM_MAPPINGS.get(chebi_id, frozenset()) for chebi_id in equivalent_ids)) - excluded_metabolites
        route = "manual"
    if not mapped:
        mapped = mapped_through(set().union(*(specific_forms(chebi_id, names, neighbours, children) for chebi_id in definition.chebi_ids)))
        route = "specific_form"
    if not mapped:
        return set(), "unmapped"
    if len(mapped) > maximum_metabolites:
        return set(), "class_term"
    return mapped, route


def build_label_table(hpo_owl_text: str, gene_to_phenotype: pd.DataFrame, chebi_relations, metabolites_table: pd.DataFrame,
                      excluded_metabolites: frozenset[str] = frozenset()) -> pd.DataFrame:
    """One row per (gene, HPO term, Human-GEM base metabolite): gene_symbol, hpo_id, base_metabolite_id, direction
    (+1, -1 or 0) fluid and route; terms that do not map, or that name a class, give no rows."""
    names, neighbours, children = chebi_relations
    metabolites_by_chebi = human_gem_metabolites_by_chebi(metabolites_table)
    mapping = {}
    for definition in parse_hpo_chemical_definitions(hpo_owl_text):
        metabolites, route = map_definition_to_metabolites(definition, names, neighbours, children, metabolites_by_chebi, excluded_metabolites)
        if metabolites:
            mapping[definition.hpo_id] = (definition, metabolites, route)
    annotations = gene_to_phenotype.drop_duplicates(["gene_symbol", "hpo_id"])
    rows = []
    for gene_symbol, hpo_id in zip(annotations.gene_symbol, annotations.hpo_id):
        if hpo_id in mapping:
            definition, metabolites, route = mapping[hpo_id]
            rows += [(gene_symbol, hpo_id, metabolite, definition.direction, definition.fluid, route) for metabolite in sorted(metabolites)]
    return pd.DataFrame(rows, columns=["gene_symbol", "hpo_id", "base_metabolite_id", "direction", "fluid", "route"])
