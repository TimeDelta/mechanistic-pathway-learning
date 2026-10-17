"""Per-edge evidence for a prior on edge presence (decision 9.57, the owner, 16 October 2026).

The owner: "I've ALWAYS wanted differently weighted prior observations based on evidence strength (not just labels)".
The three rulings the record had generalised against this (the occupancy placeholder, the threshold grid, 9.20's
overexpression coupling) say that a confidence is not a transmission magnitude; this module never sets how much
signal an edge carries. It says how sure the sources are that the edge is there, as a prior probability and a
pseudo-count that a learnable edge gate is pulled toward (models/edge_gates.py).

What enters, per edge of a built graph, read from the edge's own evidence_source and node ids so that every graph
variant is served by one table:

  signalling (OmniPath)        the distinct PubMed ids among the raw row's references, counted apart for the
                               text-mined resource (ProtMapper assembles statements from the REACH and Sparser readers)
                               and for the rest; complexes are expanded to their members as the build expands them
  transcription (CollecTRI)    the distinct PubMed ids; an edge whose only supporting resources are the ExTRI text-
                               mining layers is text-mined only
  enzyme regulation (BRENDA)   the reference count the build wrote into evidence_source
  curated (this project)       one PubMed id, the one evidence_source names
  IntAct curated absence       a pair of proteins curators recorded as not interacting: an absence observation on
                               every edge between them, in either direction (graph/edge_negatives.py)

What does not enter: Human-GEM's metabolic edges (its reaction confidence score is 0 for 12,828 of 12,877 reactions,
so it discriminates nothing), the definitional and derived edges (encodes, the membrane-potential and ion-pool
layers, drug mechanism edges with their own magnitude, plasma binders, cofactors, Reactome's neuronal layer). A
family without a row in the tier table is ungated: its edges carry exactly what they carry today.

The tier table (configs/edge_prior_tiers.csv) gives each (family, kind) a prior probability that the edge is present
and a unit strength; the strength of an edge's prior is unit_strength * min(n, cap) ** tempering over its distinct
records, as on the label side (evidence/soft_priors.py: Cusick et al. 2009 support the step from one independent
publication to two, and little beyond). An IntAct absence is a Beta update: the prior's pseudo-counts are
alpha = p * n and beta = (1 - p) * n, the absence adds INTACT_ABSENCE_PSEUDO_COUNT to beta, and the edge's prior
becomes alpha / (alpha + beta) at strength alpha + beta.

Curated edges as edge labels (decision 9.63 (b), the owner, 16 October 2026). Each edge also carries the PubMed ids its
record count counts and the resources that state it, so that a paper's species can be read: experiments/
fetch_icite_species.py caches NIH iCite's human, animal and molecular/cellular scores per PubMed id, and an edge's
human records are its ids whose human score is above 0. An edge of a curated kind, not contested by an IntAct absence,
with at least two independent records (distinct PubMed ids) of which at least one is human evidence is an edge label:
taken as known, with no gate (curated_edge_labels). The owner: "Curated edges for other species should not be
considered labels unless there's corroborating evidence from humans". The species is the paper's, not the
interaction's: a paper that shows an interaction in mice and names a human cohort is human evidence. A transcription
edge's ids are CollecTRI's, which attributes most of them to itself rather than to the resource behind them, so they
can include ids of the ExTRI text-mining layer, as the record count of 9.57 always could.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_TIERS_PATH = Path("configs/edge_prior_tiers.csv")
TEXT_MINED_REFERENCE_RESOURCES = ("ProtMapper",)
TEXT_MINED_TRANSCRIPTION_SOURCES = ("ExTRI_CollecTRI", "ExTRI2_CollecTRI2")
AGGREGATOR_SOURCES = ("CollecTRI", "CollecTRI2", "DoRothEA")
INTACT_ABSENCE_PSEUDO_COUNT = 1.0
UNDIRECTED_RELATIONS = ("binds",)  # a directed edge never takes the evidence of its reverse, which is another claim
DEFAULT_TEMPERING = 0.5
DEFAULT_COUNT_CAP = 2
PMID = re.compile(r"^\d{1,9}$")
BRENDA_REFERENCES = re.compile(r"\((\d+) references?")
CURATED_PMID = re.compile(r"PMID[ :]?(\d{4,9})")
FAMILY_OF_EVIDENCE_SOURCE = (("OmniPath small_molecule_protein", "small_molecule"), ("OmniPath", "signalling"), ("CollecTRI", "transcription"),
                             ("BRENDA", "enzyme_regulation"), ("curated", "curated"))
LABEL_FAMILIES = ("signalling", "transcription", "curated")  # the families whose curated kind can become an edge label (decision 9.63)
DEFAULT_LABEL_MINIMUM_RECORDS = 2  # "two or more independent records" (decision 9.63 (b))
DEFAULT_LABEL_MINIMUM_HUMAN_RECORDS = 1  # "corroborating evidence from humans" (the owner, 16 October 2026)
# Decision 9.70 (3), the knowledge-graph practices review's M3: a PubMed id behind more than this many signalling or
# transcription edges is a bulk paper (114 on v5, among them a review of ubiquitin-conjugating enzymes behind 13,735
# edges and the NetworKIN predictor's paper behind 7,166), and the label rule can count an edge's bulk papers as one
# record ("one") or not at all ("none") instead of one record each ("count", the rule of 9.63).
BULK_PAPER_EDGES = 100
BULK_PAPER_FAMILIES = ("signalling", "transcription")
BULK_PAPER_MODES = ("count", "one", "none")
DEFAULT_HUMAN_SHARE_ABOVE = 0.0  # a paper is human evidence when iCite's human score is above this
SPECIES_COLUMNS = ("human", "animal", "molecular_cellular")


@dataclass(frozen=True)
class EdgeTier:
    family: str
    kind: str
    prior_probability: float
    unit_strength: float
    note: str = ""


def read_edge_tiers(path: Path = DEFAULT_TIERS_PATH) -> dict[tuple[str, str], EdgeTier]:
    """configs/edge_prior_tiers.csv as {(family, kind): EdgeTier}; a family without rows is ungated."""
    tiers = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            probability, strength = float(row["prior_probability"]), float(row["unit_strength"])
            if not 0.0 < probability < 1.0 or strength < 0.0:
                raise ValueError(f"{path}: {row['family']} {row['kind']} needs a probability strictly inside (0, 1) and a strength >= 0")
            tiers[(row["family"].strip(), row["kind"].strip())] = EdgeTier(row["family"].strip(), row["kind"].strip(), probability, strength, row.get("note", ""))
    return tiers


def family_of(evidence_source: str) -> str:
    text = str(evidence_source)
    for prefix, family in FAMILY_OF_EVIDENCE_SOURCE:
        if text.startswith(prefix):
            return family
    return "other"


UNIPROT_ACCESSION = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")


def symbol_of(node_id: str) -> str:
    """GENE:X and PROTEIN:X -> X, with the accession suffix of a gene's second protein product (PROTEIN:CALCA:P01258)
    dropped; any other node id is returned as it is."""
    text = str(node_id)
    for prefix in ("GENE:", "PROTEIN:"):
        if text.startswith(prefix):
            name = text[len(prefix):]
            head, _, tail = name.rpartition(":")
            return head if head and UNIPROT_ACCESSION.match(tail) else name  # GENE:CHEBI:15996 keeps its whole name
    return text


def accession_of(node_id: str) -> str | None:
    """The UniProt accession a protein node names: its suffix (PROTEIN:CALCA:P01258) or its whole name (PROTEIN:A0A087WWS6)."""
    text = str(node_id)
    if not text.startswith("PROTEIN:"):
        return None
    for part in reversed(text[len("PROTEIN:"):].split(":")):
        if UNIPROT_ACCESSION.match(part):
            return part.split("-")[0]
    return None


def pmids_of(references: str, excluded_resources: tuple[str, ...] = ()) -> tuple[set[str], set[str]]:
    """(PubMed ids from the kept resources, PubMed ids from the excluded ones) of a 'Resource:PMID;...' cell."""
    kept, excluded = set(), set()
    if not isinstance(references, str):
        return kept, excluded
    for reference in references.split(";"):
        if ":" not in reference:
            continue
        resource, identifier = reference.split(":", 1)
        identifier = identifier.strip()
        if not PMID.match(identifier):
            continue
        (excluded if resource.strip() in excluded_resources else kept).add(identifier)
    return kept, excluded


def resources_of(sources) -> set[str]:
    """The resources of a 'SIGNOR;HPRD;...' cell."""
    return {source.strip() for source in str(sources).split(";") if source.strip()} if isinstance(sources, str) else set()


def joined_pubmed_ids(pubmed_ids) -> str:
    """PubMed ids as one ';'-joined cell in numeric order ('' for none), the form edge_evidence.parquet stores."""
    return ";".join(sorted({str(identifier) for identifier in pubmed_ids}, key=int))


def split_pubmed_ids(cell) -> list[str]:
    """The PubMed ids of a ';'-joined cell; anything that is not an id is dropped."""
    return [identifier for identifier in cell.split(";") if PMID.match(identifier)] if isinstance(cell, str) and cell else []


def members(identifier: str) -> list[str]:
    """A complex 'A_B_C' (or 'COMPLEX:A_B_C') as its members; a single protein as itself; isoform suffixes dropped."""
    text = str(identifier)
    if text.startswith("COMPLEX:"):
        text = text[len("COMPLEX:"):]
    return [part.split("-")[0] if UNIPROT_ACCESSION.match(part) else part for part in text.split("_") if part]


def signalling_references(omnipath: pd.DataFrame) -> pd.DataFrame:
    """(source, target) -> curated and text-mined PubMed id sets, over OmniPath rows expanded to members, keyed by gene
    symbols and, where the table has them, by UniProt accessions too (a protein node named by accession matches there)."""
    records: dict[tuple[str, str], tuple[set[str], set[str], set[str]]] = {}
    has_accessions = "source" in omnipath.columns and "target" in omnipath.columns
    has_resources = "sources" in omnipath.columns
    for row in omnipath.itertuples(index=False):
        curated, text_mined = pmids_of(row.references, TEXT_MINED_REFERENCE_RESOURCES)
        resources = resources_of(row.sources) if has_resources else set()
        key_sets = [(members(row.source_genesymbol), members(row.target_genesymbol))]
        if has_accessions:
            key_sets.append((members(row.source), members(row.target)))
        for sources, targets in key_sets:
            for source in sources:
                for target in targets:
                    kept, mined, stated_by = records.setdefault((source, target), (set(), set(), set()))
                    kept |= curated
                    mined |= text_mined
                    stated_by |= resources
    return pd.DataFrame([{"family": "signalling", "source_symbol": source, "target_symbol": target, "curated_pmids": len(kept),
                          "text_mined_pmids": len(mined - kept), "text_mined_only": not kept and bool(mined),
                          "curated_pubmed_ids": joined_pubmed_ids(kept), "text_mined_pubmed_ids": joined_pubmed_ids(mined - kept),
                          "resources": ";".join(sorted(stated_by))}
                         for (source, target), (kept, mined, stated_by) in sorted(records.items())])


def transcription_references(collectri: pd.DataFrame) -> pd.DataFrame:
    """(TF, target) -> distinct PubMed ids, and whether every supporting resource is a text-mining layer."""
    records: dict[tuple[str, str], tuple[set[str], set[str]]] = {}
    for row in collectri.itertuples(index=False):
        pmids, _ = pmids_of(row.references)
        sources = resources_of(row.sources) - set(AGGREGATOR_SOURCES)
        for source in members(row.source_genesymbol):
            for target in members(row.target_genesymbol):
                kept, supporting = records.setdefault((source, target), (set(), set()))
                kept |= pmids
                supporting |= sources
    return pd.DataFrame([{"family": "transcription", "source_symbol": source, "target_symbol": target, "curated_pmids": len(pmids),
                          "text_mined_pmids": 0, "text_mined_only": bool(supporting) and supporting <= set(TEXT_MINED_TRANSCRIPTION_SOURCES),
                          "curated_pubmed_ids": joined_pubmed_ids(pmids), "text_mined_pubmed_ids": "", "resources": ";".join(sorted(supporting))}
                         for (source, target), (pmids, supporting) in sorted(records.items())])


def tempered(count, cap: int = DEFAULT_COUNT_CAP, exponent: float = DEFAULT_TEMPERING):
    return np.minimum(np.asarray(count, dtype=float), float(cap)) ** float(exponent)


def edge_evidence_table(edges: pd.DataFrame, references: pd.DataFrame, negative_pairs: set[frozenset[str]]) -> pd.DataFrame:
    """One row per graph edge, in edge order: family, kind, records (distinct PubMed ids or references), contested."""
    by_key = {(row.family, row.source_symbol, row.target_symbol): row for row in references.itertuples(index=False)} if len(references) else {}
    families, kinds, records, contested, pubmed, resources = [], [], [], [], [], []
    for edge in edges.itertuples(index=False):
        family = family_of(edge.evidence_source)
        source, target = symbol_of(edge.source_id), symbol_of(edge.target_id)
        kind, count, identifiers, stated_by = "none", 0, "", ""
        if family in ("signalling", "transcription"):
            undirected = str(edge.relation_type) in UNDIRECTED_RELATIONS  # a binding is one claim whichever way the build stored it
            row = None
            for first, second in ((source, target), (accession_of(edge.source_id) or source, accession_of(edge.target_id) or target)):
                row = by_key.get((family, first, second)) or (by_key.get((family, second, first)) if undirected else None)
                if row is not None:
                    break
            if row is None:
                kind = "unmatched"
            else:
                stated_by = str(getattr(row, "resources", "") or "")
                if row.text_mined_only:  # first: an ExTRI-only transcription edge carries PubMed ids that came from text mining
                    kind, count = "text_mined_only", int(row.text_mined_pmids or row.curated_pmids)
                    identifiers = str((getattr(row, "text_mined_pubmed_ids", "") if row.text_mined_pmids else getattr(row, "curated_pubmed_ids", "")) or "")
                elif row.curated_pmids > 0:
                    kind, count = "curated", int(row.curated_pmids)
                    identifiers = str(getattr(row, "curated_pubmed_ids", "") or "")
                else:
                    kind = "unreferenced"
        elif family == "enzyme_regulation":
            match = BRENDA_REFERENCES.search(str(edge.evidence_source))
            kind, count, stated_by = "in_vitro", int(match.group(1)) if match else 0, "BRENDA"
        elif family == "curated":
            found = set(CURATED_PMID.findall(str(edge.evidence_source)))
            kind, count, identifiers, stated_by = "curated", len(found), joined_pubmed_ids(found), str(edge.evidence_source)
        families.append(family)
        kinds.append(kind)
        records.append(count)
        contested.append(source != target and frozenset((source, target)) in negative_pairs)
        pubmed.append(identifiers)
        resources.append(stated_by)
    return pd.DataFrame({"family": families, "kind": kinds, "records": records, "contested": contested, "pubmed_ids": pubmed, "resources": resources},
                        index=edges.index)


def read_species(path: Path) -> pd.DataFrame:
    """The species scores experiments/fetch_icite_species.py cached, indexed by PubMed id (str): human, animal and
    molecular_cellular for each id iCite found."""
    table = pd.read_parquet(path)
    if "found" in table.columns:
        table = table[table["found"].astype(bool)]
    return table.assign(pmid=table["pmid"].astype(str)).drop_duplicates("pmid").set_index("pmid")[list(SPECIES_COLUMNS)]


def species_of_records(pubmed_ids: pd.Series, species: pd.DataFrame, human_share_above: float = DEFAULT_HUMAN_SHARE_ABOVE) -> pd.DataFrame:
    """Per edge, over its PubMed ids: how many iCite scored (species_known_records), how many are human evidence (human
    score above human_share_above: human_records), how many record animals and no human (animal_only_records), and each
    scored id with its human, animal and molecular/cellular scores ('11:0.20/0.00/0.80;...', pubmed_species), so that a
    reviewer of a closure sees the species paper by paper (decision 9.63 (c))."""
    known, human, animal_only, listed = [], [], [], []
    human_score, animal_score, cellular_score = species["human"].to_dict(), species["animal"].to_dict(), species["molecular_cellular"].to_dict()
    for cell in pubmed_ids:
        identifiers = [identifier for identifier in split_pubmed_ids(cell) if identifier in human_score]
        is_human = [float(human_score[identifier]) > human_share_above for identifier in identifiers]
        known.append(len(identifiers))
        human.append(int(sum(is_human)))
        animal_only.append(int(sum(1 for identifier, flag in zip(identifiers, is_human) if not flag and float(animal_score[identifier]) > 0.0)))
        listed.append(";".join(f"{identifier}:{float(human_score[identifier]):.2f}/{float(animal_score[identifier]):.2f}/{float(cellular_score[identifier]):.2f}"
                               for identifier in identifiers))
    return pd.DataFrame({"species_known_records": known, "human_records": human, "animal_only_records": animal_only, "pubmed_species": listed},
                        index=pubmed_ids.index)


def bulk_paper_ids(evidence: pd.DataFrame, threshold: int = BULK_PAPER_EDGES) -> set[str]:
    """The PubMed ids that stand behind more than threshold signalling or transcription edges of the table."""
    counts: dict[str, int] = {}
    for family, cell in zip(evidence["family"], evidence["pubmed_ids"]):
        if family in BULK_PAPER_FAMILIES:
            for identifier in split_pubmed_ids(cell):
                counts[identifier] = counts.get(identifier, 0) + 1
    return {identifier for identifier, count in counts.items() if count > threshold}


def records_with_bulk_papers(evidence: pd.DataFrame, mode: str, threshold: int = BULK_PAPER_EDGES) -> np.ndarray:
    """Each edge's record count with its bulk papers counted one each ("count", the table's own records), together as
    one record ("one") or not at all ("none"). Only the signalling and transcription families change, whose records
    are their distinct PubMed ids; every other family keeps its records."""
    if mode not in BULK_PAPER_MODES:
        raise ValueError(f"bulk paper mode must be one of {BULK_PAPER_MODES}")
    records = evidence["records"].to_numpy(dtype=float).copy()
    if mode == "count":
        return records
    bulk = bulk_paper_ids(evidence, threshold)
    for position, (family, cell) in enumerate(zip(evidence["family"], evidence["pubmed_ids"])):
        if family not in BULK_PAPER_FAMILIES:
            continue
        identifiers = split_pubmed_ids(cell)
        bulk_count = sum(identifier in bulk for identifier in identifiers)
        if bulk_count:
            records[position] = len(identifiers) - bulk_count + (1 if mode == "one" else 0)
    return records


def curated_edge_labels(evidence: pd.DataFrame, minimum_records: int = DEFAULT_LABEL_MINIMUM_RECORDS,
                        minimum_human_records: int = DEFAULT_LABEL_MINIMUM_HUMAN_RECORDS, minimum_sign_records: int | None = None,
                        exclude_disputed_signs: bool = False, bulk_papers: str = "count") -> np.ndarray:
    """Decision 9.63 (b) with the owner's human-evidence rule: True for an edge of a curated kind (LABEL_FAMILIES), not
    contested by an IntAct absence, with at least minimum_records independent records of which at least
    minimum_human_records are human evidence. Such an edge is an edge label: taken as known, with no gate.

    minimum_sign_records (decision 9.70, off by default): when given, a signed signalling or transcription edge is a
    label only if its sign itself rests on that many PubMed ids and, for a regulation, CollecTRI decided the sign from
    a paper (graph/sign_evidence.py), since a label's sign is one no gate can flip.

    exclude_disputed_signs (decision 9.72, off by default): a regulation whose two CollecTRI releases state opposite
    signs, or one of which states both, is not a label (graph/sign_evidence.py, sign_disputed).

    bulk_papers (decision 9.70 (3), "count" by default): how an edge's bulk papers count toward minimum_records,
    records_with_bulk_papers."""
    if "human_records" not in evidence.columns:
        raise ValueError("the edge evidence table has no human_records column; rebuild it with experiments/build_edge_evidence.py, which reads the species "
                         "experiments/fetch_icite_species.py caches")
    records = records_with_bulk_papers(evidence, bulk_papers)
    labels = ((evidence["kind"] == "curated") & evidence["family"].isin(LABEL_FAMILIES) & ~evidence["contested"].astype(bool)
              & (evidence["human_records"] >= minimum_human_records)).to_numpy() & (records >= minimum_records)
    if minimum_sign_records is not None:
        from mechanistic_pathway_learning.graph.sign_evidence import sign_is_evidenced

        labels &= sign_is_evidenced(evidence, minimum_sign_records)
    if exclude_disputed_signs:
        from mechanistic_pathway_learning.graph.sign_evidence import sign_is_undisputed

        labels &= sign_is_undisputed(evidence)
    return labels


def edge_priors(evidence: pd.DataFrame, tiers: dict[tuple[str, str], EdgeTier], global_strength: float = 1.0, tempering: float = DEFAULT_TEMPERING,
                count_cap: int = DEFAULT_COUNT_CAP, absence_pseudo_count: float = INTACT_ABSENCE_PSEUDO_COUNT) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(gated, prior probability, pseudo-count) per edge. An edge is gated when its (family, kind) has a tier; an
    ungated edge carries NaN and 0. A record count of 0 still gives the tier's unit strength once (min(n, cap) is read
    as max(n, 1) so an unreferenced tier is not a zero prior); an IntAct absence is a Beta update on the gated edge."""
    gated = np.zeros(len(evidence), dtype=bool)
    probability = np.full(len(evidence), np.nan)
    strength = np.zeros(len(evidence))
    for position, row in enumerate(evidence.itertuples(index=False)):
        tier = tiers.get((row.family, row.kind))
        if tier is None:
            continue
        count = global_strength * tier.unit_strength * float(tempered(max(int(row.records), 1), count_cap, tempering))
        alpha, beta = tier.prior_probability * count, (1.0 - tier.prior_probability) * count
        if row.contested:
            beta += global_strength * absence_pseudo_count
        gated[position] = True
        total = alpha + beta
        probability[position] = alpha / total if total > 0 else tier.prior_probability
        strength[position] = total
    return gated, probability, strength


def permute_edge_priors(gated: np.ndarray, probability: np.ndarray, strength: np.ndarray, relation: np.ndarray, random_seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """The control: the priors of the gated edges shuffled among the gated edges of the same relation, the wiring kept,
    so a relation keeps its distribution of priors and loses which edge holds which."""
    probability, strength = np.array(probability, copy=True), np.array(strength, copy=True)
    generator = np.random.default_rng(random_seed)
    relation = np.asarray(relation)
    for value in np.unique(relation[gated]):
        positions = np.flatnonzero(gated & (relation == value))
        order = generator.permutation(len(positions))
        probability[positions], strength[positions] = probability[positions][order], strength[positions][order]
    return probability, strength


def summarise(evidence: pd.DataFrame, gated: np.ndarray, probability: np.ndarray, strength: np.ndarray) -> dict:
    by = evidence.assign(gated=gated, probability=probability, strength=strength).groupby(["family", "kind"])
    return {"edges": int(len(evidence)), "gated": int(gated.sum()), "contested": int(evidence.contested.sum()),
            "by_family_and_kind": {f"{family}:{kind}": {"edges": int(len(group)), "gated": int(group.gated.sum()),
                                                        "mean_prior": None if not group.gated.any() else float(group.probability[group.gated].mean()),
                                                        "mean_strength": float(group.strength.mean())}
                                   for (family, kind), group in by}}
