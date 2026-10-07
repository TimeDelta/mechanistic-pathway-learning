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
    """(Human-GEM base metabolites, route): route is 'equivalent', 'specific_form', 'class_term' (more than
    maximum_metabolites reached; no metabolites returned) or 'unmapped'."""
    def mapped_through(chebi_ids: set[str]) -> set[str]:
        return set().union(*(metabolites_by_chebi.get(chebi_id, set()) for chebi_id in chebi_ids)) - excluded_metabolites

    mapped = mapped_through(set().union(*(equivalent_forms(chebi_id, neighbours) for chebi_id in definition.chebi_ids)))
    route = "equivalent"
    if not mapped:
        mapped = mapped_through(set().union(*(specific_forms(chebi_id, names, neighbours, children) for chebi_id in definition.chebi_ids)))
        route = "specific_form"
    if not mapped:
        return set(), "unmapped"
    if len(mapped) > maximum_metabolites:
        return set(), "class_term"
    return mapped, route
