"""The allele, not the gene, is the perturbation (decision 9.54 (d), the owner, 16 October 2026: "Huntington's proposal d
accepted").

Why. A gene was one perturbation with one sign, so a gene whose diseases arise from alleles that act in opposite
directions had to take one of them: HTT was seeded as a loss of function while Huntington disease is a gain of
function, and HTT's intellectual-disability label comes from another allele altogether (compound heterozygous
variants, Rodan et al. 2016, doi 10.1038/ejhg.2016.74). Here each gene's monogenic reports are split by the mechanism
of the disease they come from, and each mechanism class is a perturbation of its own:

  loss of function and untyped   the gene itself, seeded -1 as before (an untyped association is silence about the
                                 mechanism; it keeps the default and is marked untyped)
  gain of function (hypermorph)  "<gene>|gain_of_function", seeded +1 on the gene: more of the normal activity
  gain of function (neomorph)    a curated allele (configs/neomorphic_alleles.csv): a node of its own, present only where
                                 it is seeded, whose mechanism edges state what the mutant product does, each with a
                                 PubMed id and a quote

Every allele of a gene keeps the gene as its leakage group, so no split puts one allele in training and another in
the test.

Where a mechanism comes from. Orphanet types some gene-disease associations as loss or gain of function (Orphadata
product 6). The same disease is often recorded twice, once by Orphanet and once by OMIM, and OMIM's record carries no
mechanism: HTT's Huntington disease is ORPHA:399 (typed gain of function) and OMIM:143100 (untyped). Orphadata product 1
states how an Orphanet disorder maps to OMIM: E, the two concepts are equivalent; NTBT, the Orphanet code is narrower;
BTNT, it is broader. A typed mechanism is carried, for the same gene,

  - from an Orphanet disorder to every OMIM entry mapped E to it (one disease, two records), and
  - from such an OMIM entry to every Orphanet disorder mapped NTBT to it (a narrower form of the same disease, such as
    juvenile Huntington disease, ORPHA:248111 NTBT OMIM:143100),

and never through a BTNT mapping: a broader Orphanet code is a group (ORPHA:528084, non-specific syndromic
intellectual disability, maps BTNT to 108 OMIM entries in the release of 23 June 2026) and says nothing about any one
member's mechanism. A pair Orphanet types itself keeps its own type; when two inherited mechanisms disagree on a pair,
it stays untyped and is marked as a conflict.

The second source. Where Orphanet leaves a pair untyped, Gene2Phenotype's curated molecular mechanism is read (G2P,
all panels, CC BY 4.0; the owner, 15 October 2026: "Always try to find a second reputable source when one fails").
A G2P record speaks for an OMIM disease; it reaches an Orphanet disorder through the same E and NTBT mappings as
above. Only records whose gene-disease confidence is definitive, strong or moderate type a pair (the levels G2P's
curators name as used in clinical reporting; limited, disputed and refuted do not). "Loss of function" and "dominant
negative" give a loss of function, because a dominant-negative product lowers the activity of the wild-type one and
so the gene's net activity falls, which is what the seed -1 states; "gain of function" gives a gain; "undetermined"
and "undetermined non-loss-of-function" type nothing and are written into the basis. G2P records that disagree on one
pair leave it untyped as a conflict. The second source fills an absence and does not overrule the first: where both
type a pair and disagree, Orphanet's type stands and the disagreement is written into the basis for the owner to see.
"""
from __future__ import annotations

import csv
import xml.etree.ElementTree as ElementTree
from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

LOSS_OF_FUNCTION = "loss_of_function"
GAIN_OF_FUNCTION = "gain_of_function"
UNTYPED = "untyped"
NEOMORPH = "neomorph"
HYPERMORPH = "hypermorph"
MECHANISM_OF_ORPHANET_TYPE = {
    "Disease-causing germline mutation(s) (loss of function) in": LOSS_OF_FUNCTION,
    "Disease-causing germline mutation(s) (gain of function) in": GAIN_OF_FUNCTION,
}
EXACT, NARROWER, BROADER = "E", "NTBT", "BTNT"  # product 1's DisorderMappingRelation codes, Orphanet code relative to OMIM
GAIN_OF_FUNCTION_SUFFIX = "|" + GAIN_OF_FUNCTION
DEFAULT_NEOMORPHIC_ALLELES_PATH = Path("configs/neomorphic_alleles.csv")
SIGN_OF_EFFECT = {"increases": 1.0, "decreases": -1.0}
DEFAULT_GENE2PHENOTYPE_PATH = Path("data/raw/gene2phenotype/G2P_all_panels.csv")
GENE2PHENOTYPE_MECHANISM = {"loss of function": LOSS_OF_FUNCTION, "dominant negative": LOSS_OF_FUNCTION, "gain of function": GAIN_OF_FUNCTION}
GENE2PHENOTYPE_UNTYPED_MECHANISMS = ("undetermined", "undetermined non-loss-of-function")
GENE2PHENOTYPE_CONFIDENCES = ("definitive", "strong", "moderate")
GENE2PHENOTYPE_COLUMNS = ("g2p id", "gene symbol", "previous gene symbols", "disease mim", "confidence", "molecular mechanism", "allelic requirement")


def load_orphanet_omim_mappings(product1_path: Path) -> dict[str, list[tuple[str, str]]]:
    """ORPHA:<code> -> [(OMIM:<number>, relation code)] from Orphadata product 1 (en_product1.xml, CC BY 4.0)."""
    mappings: dict[str, list[tuple[str, str]]] = defaultdict(list)
    root = ElementTree.parse(product1_path).getroot()
    for disorder in root.iter("Disorder"):
        code = disorder.findtext("OrphaCode")
        if not code:
            continue
        for reference in disorder.iter("ExternalReference"):
            if reference.findtext("Source") != "OMIM":
                continue
            relation_element = reference.find("DisorderMappingRelation")
            relation = (relation_element.findtext("Name") or "") if relation_element is not None else ""
            mappings[f"ORPHA:{code}"].append((f"OMIM:{reference.findtext('Reference')}", relation.split(" ")[0]))
    return dict(mappings)


@dataclass(frozen=True)
class DiseaseMechanism:
    mechanism: str  # LOSS_OF_FUNCTION, GAIN_OF_FUNCTION or UNTYPED
    basis: str  # "orphanet typed", "same disease as ORPHA:x", "narrower than OMIM:y, same disease as ORPHA:x", "untyped" or "conflict: ..."


def disease_mechanisms(association_types: Mapping[tuple[str, str], str], mappings: Mapping[str, list[tuple[str, str]]],
                       pairs: Iterable[tuple[str, str]]) -> dict[tuple[str, str], DiseaseMechanism]:
    """The mechanism of each (gene, disease) pair in pairs, by the rules of the module docstring. association_types
    holds the Orphanet type of an (gene, ORPHA:x) association (product 6) and any other type for other pairs."""
    pairs = set(pairs)
    omim_to_narrower: dict[str, list[str]] = defaultdict(list)
    for orpha, targets in mappings.items():
        for omim, relation in targets:
            if relation == NARROWER:
                omim_to_narrower[omim].append(orpha)
    direct: dict[tuple[str, str], str] = {}
    reached: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)  # pair -> {inherited mechanism: basis}
    for (gene, disease), association_type in association_types.items():
        mechanism = MECHANISM_OF_ORPHANET_TYPE.get(association_type)
        if mechanism is None or not disease.startswith("ORPHA:"):
            continue
        direct[(gene, disease)] = mechanism
        for omim, relation in mappings.get(disease, []):
            if relation != EXACT:
                continue
            reached[(gene, omim)].setdefault(mechanism, f"same disease as {disease}")
            for narrower in omim_to_narrower.get(omim, []):
                if narrower != disease:
                    reached[(gene, narrower)].setdefault(mechanism, f"narrower than {omim}, same disease as {disease}")
    result: dict[tuple[str, str], DiseaseMechanism] = {}
    for pair in pairs:
        if pair in direct:
            result[pair] = DiseaseMechanism(direct[pair], "orphanet typed")
            continue
        found = reached.get(pair, {})
        if len(found) == 1:
            (mechanism, basis), = found.items()
            result[pair] = DiseaseMechanism(mechanism, basis)
        elif len(found) > 1:
            result[pair] = DiseaseMechanism(UNTYPED, "conflict: " + "; ".join(f"{mechanism} ({basis})" for mechanism, basis in sorted(found.items())))
        else:
            result[pair] = DiseaseMechanism(UNTYPED, "untyped")
    return result


@dataclass(frozen=True)
class Gene2PhenotypeRecord:
    g2p_id: str
    gene: str
    previous_symbols: tuple[str, ...]
    disease: str  # "OMIM:<number>"
    confidence: str
    mechanism: str  # G2P's own words, "loss of function", "dominant negative", ...
    allelic_requirement: str


def load_gene2phenotype(path: Path) -> tuple[list[Gene2PhenotypeRecord], dict]:
    """The records of G2P's all-panel download (www.ebi.ac.uk/gene2phenotype/api/panel/all/download) that name an
    OMIM disease, and a count of those that name none (a G2P disease with only a Mondo id cannot meet an HPO report)."""
    records: list[Gene2PhenotypeRecord] = []
    without_omim = 0
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [column for column in GENE2PHENOTYPE_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} lacks the G2P columns {missing}; the download format may have changed")
        for row in reader:
            disease_mim = row["disease mim"].strip()
            if not disease_mim:
                without_omim += 1
                continue
            previous = tuple(symbol.strip() for symbol in row["previous gene symbols"].split(";") if symbol.strip())
            records.append(Gene2PhenotypeRecord(row["g2p id"].strip(), row["gene symbol"].strip(), previous, f"OMIM:{disease_mim}", row["confidence"].strip().lower(),
                                                row["molecular mechanism"].strip().lower(), row["allelic requirement"].strip()))
    return records, {"records_with_an_omim_disease": len(records), "records_without_an_omim_disease": without_omim}


def gene2phenotype_mechanisms(records: Iterable[Gene2PhenotypeRecord], mappings: Mapping[str, list[tuple[str, str]]],
                              pairs: Iterable[tuple[str, str]]) -> dict[tuple[str, str], DiseaseMechanism]:
    """The mechanism G2P states for each (gene, disease) pair it covers, by the rules of the module docstring. A pair G2P
    does not cover is absent from the result. A gene is matched by its current symbol, and by a previous symbol that is
    no G2P gene's current symbol."""
    usable = [record for record in records if record.confidence in GENE2PHENOTYPE_CONFIDENCES]
    current_symbols = {record.gene for record in usable}
    by_gene_and_disease: dict[tuple[str, str], list[Gene2PhenotypeRecord]] = defaultdict(list)
    for record in usable:
        by_gene_and_disease[(record.gene, record.disease)].append(record)
        for symbol in record.previous_symbols:
            if symbol not in current_symbols:
                by_gene_and_disease[(symbol, record.disease)].append(record)
    result: dict[tuple[str, str], DiseaseMechanism] = {}
    for gene, disease in set(pairs):
        if disease.startswith("OMIM:"):
            routes = [(disease, "")]
        elif disease.startswith("ORPHA:"):
            routes = [(omim, f", for {omim}, the same disease" if relation == EXACT else f", for {omim}, of which {disease} is a narrower form")
                      for omim, relation in mappings.get(disease, []) if relation in (EXACT, NARROWER)]
        else:
            routes = []
        found = [(record, route) for omim, route in routes for record in by_gene_and_disease.get((gene, omim), [])]
        if not found:
            continue
        typed: dict[str, list[str]] = defaultdict(list)
        undetermined = []
        for record, route in found:
            description = f"gene2phenotype {record.g2p_id} ({record.mechanism}, {record.confidence}){route}"
            if record.mechanism in GENE2PHENOTYPE_MECHANISM:
                typed[GENE2PHENOTYPE_MECHANISM[record.mechanism]].append(description)
            else:
                undetermined.append(description)
        if len(typed) == 1:
            (mechanism, descriptions), = typed.items()
            result[(gene, disease)] = DiseaseMechanism(mechanism, "; ".join(sorted(set(descriptions))))
        elif len(typed) > 1:
            result[(gene, disease)] = DiseaseMechanism(UNTYPED, "conflict: " + "; ".join(sorted({description for descriptions in typed.values() for description in descriptions})))
        else:
            result[(gene, disease)] = DiseaseMechanism(UNTYPED, "untyped; " + "; ".join(sorted(set(undetermined))))
    return result


def first_source_unless_absent(first: Mapping[tuple[str, str], DiseaseMechanism], second: Mapping[tuple[str, str], DiseaseMechanism],
                               pairs: Iterable[tuple[str, str]]) -> dict[tuple[str, str], DiseaseMechanism]:
    """The first source's mechanism, with the second filling a pair the first leaves untyped. Where both type a pair
    and disagree the first stands and the basis says who disputes it; a conflict within the first is not settled by
    the second, whose view is added to the basis."""
    result: dict[tuple[str, str], DiseaseMechanism] = {}
    for pair in set(pairs):
        primary = first.get(pair, DiseaseMechanism(UNTYPED, "untyped"))
        secondary = second.get(pair)
        if secondary is None:
            result[pair] = primary
        elif primary.mechanism != UNTYPED:
            disputed = secondary.mechanism not in (UNTYPED, primary.mechanism)
            result[pair] = DiseaseMechanism(primary.mechanism, f"{primary.basis}; disputed by {secondary.basis}") if disputed else primary
        elif primary.basis.startswith("conflict"):
            result[pair] = DiseaseMechanism(UNTYPED, f"{primary.basis}; {secondary.basis}")
        else:
            result[pair] = secondary
    return result


@dataclass(frozen=True)
class AlleleMechanismEdge:
    target_symbol: str
    target_kind: str  # "protein" or "gene" (the gene node: a transcriptional effect)
    sign: float
    magnitude: float
    pmid: str
    quote: str


@dataclass(frozen=True)
class NeomorphicAllele:
    allele_id: str  # the perturbation id, "<gene>|<allele>"
    gene: str
    label: str
    edges: tuple[AlleleMechanismEdge, ...] = field(default_factory=tuple)

    @property
    def node_id(self) -> str:
        return "ALLELE:" + self.allele_id.replace("|", "_")


def read_neomorphic_alleles(path: Path = DEFAULT_NEOMORPHIC_ALLELES_PATH) -> dict[str, NeomorphicAllele]:
    """gene -> its curated neomorphic gain-of-function allele, from one row per mechanism edge. Every edge needs a PubMed
    id and a quote; an effect is "increases" or "decreases"; a magnitude is a stated quantity, 1 where the source states
    none, as a knockout's seed is 1."""
    edges_by_allele: dict[tuple[str, str, str], list[AlleleMechanismEdge]] = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["effect"] not in SIGN_OF_EFFECT:
                raise ValueError(f"{path}: effect must be one of {sorted(SIGN_OF_EFFECT)}, not {row['effect']!r}")
            if row["target_kind"] not in ("protein", "gene"):
                raise ValueError(f"{path}: target_kind must be protein or gene, not {row['target_kind']!r}")
            if not row["pmid"].strip().isdigit() or not row["quote"].strip():
                raise ValueError(f"{path}: every mechanism edge of {row['allele_id']} needs a PubMed id and a quote")
            magnitude = float(row["magnitude"])
            if not 0.0 < magnitude <= 1.0:
                raise ValueError(f"{path}: magnitude must be in (0, 1]")
            edges_by_allele[(row["allele_id"], row["gene"], row["label"])].append(
                AlleleMechanismEdge(row["target_symbol"], row["target_kind"], SIGN_OF_EFFECT[row["effect"]], magnitude, row["pmid"].strip(), row["quote"].strip()))
    alleles: dict[str, NeomorphicAllele] = {}
    for (allele_id, gene, label), edges in edges_by_allele.items():
        if not allele_id.startswith(gene + "|"):
            raise ValueError(f"{path}: allele id {allele_id!r} must start with its gene {gene!r} and a bar")
        if gene in alleles:
            raise ValueError(f"{path}: {gene} has two neomorphic alleles; key them by disease before adding a second")
        alleles[gene] = NeomorphicAllele(allele_id, gene, label, tuple(edges))
    return alleles


@dataclass(frozen=True)
class Allele:
    perturbation_id: str
    allele_class: str  # LOSS_OF_FUNCTION, UNTYPED, HYPERMORPH or NEOMORPH
    seed_sign: float | None  # -1 or +1 on the gene node; None for a neomorph, seeded through its own node
    label: str


def allele_of(gene: str, mechanism: DiseaseMechanism, neomorphic_alleles: Mapping[str, NeomorphicAllele]) -> Allele:
    """The perturbation a monogenic report of this gene and disease belongs to."""
    if mechanism.mechanism == GAIN_OF_FUNCTION:
        if gene in neomorphic_alleles:
            allele = neomorphic_alleles[gene]
            return Allele(allele.allele_id, NEOMORPH, None, allele.label)
        return Allele(gene + GAIN_OF_FUNCTION_SUFFIX, HYPERMORPH, 1.0, f"{gene} (gain of function)")
    return Allele(gene, mechanism.mechanism, -1.0, gene)


def gene_of_perturbation(perturbation_id: str) -> str:
    """The gene an allele perturbation belongs to: the part before the bar ("HTT|CAG_repeat_expansion" -> "HTT")."""
    return str(perturbation_id).split("|", 1)[0]
