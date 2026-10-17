"""The evidence behind an edge's sign, as distinct from the evidence behind the pair (decision 9.70, the knowledge-graph
practices review, `docs/knowledge_graph_practices_review.md` M1 and M2).

The edge evidence table (graph/edge_evidence.py) counts the distinct PubMed ids behind a pair of nodes. OmniPath keeps
the support for an interaction's positive sign, its negative sign, its direction and its mere occurrence apart (the
`evidences` field), and CollecTRI states for each regulation how its sign was decided: from a paper ("PMID"), by a
default it applied when no paper stated one ("default activation"; Mueller-Dott et al. 2023: "For the remaining 21730
TF-gene interaction an activating mode was assigned by default"), or from the factor's role across its regulon
("regulon"). A pair can therefore have two curated papers while its sign has none, and decision 9.63 takes such an edge
as known, with a sign no gate can flip. This module reads the two evidences downloads into, per (source, target):

  signalling      the PubMed ids behind the positive and behind the negative sign, over the curated resources (the
                  text-mining layers excluded, as edge_evidence.py excludes them), complexes expanded to members;
  transcription   the same from CollecTRI's own entries, with the set of its sign decisions.

and gives each graph edge the ids behind its stored sign (`sign_pubmed_ids`, `sign_records`) and, for a regulation,
how its sign was decided (`sign_basis`). Nothing here changes a label by itself: curated_edge_labels reads the columns
only when asked for a minimum (`minimum_sign_records`).

The CollecTRI download holds two releases of each regulation it covers (decision 9.72): the 2023 release's own entry
(resource CollecTRI, the one with a sign decision) and a later release's (resource CollecTRI2, dataset collectri2,
with no sign decision). They state opposite signs for 3,795 regulation edges of the v5 graph, and the sign the graph
stores is OmniPath's consensus between them. `sign_disputed` marks a regulation whose two releases state opposite
signs, or one of which states both; curated_edge_labels reads it only when asked (`exclude_disputed_signs`).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.edge_evidence import accession_of, family_of, joined_pubmed_ids, members, split_pubmed_ids, symbol_of

TEXT_MINING_RESOURCES = frozenset({"ProtMapper", "REACH", "Sparser", "ExTRI", "ExTRI2", "ExTRI_CollecTRI", "ExTRI2_CollecTRI2"})
SIGNALLING_DATASETS = frozenset({"omnipath", "pathwayextra", "ligrecextra"})
TRANSCRIPTION_RESOURCE = "CollecTRI"  # its own entry carries the sign decision; the resources it aggregates carry none
LATER_TRANSCRIPTION_RELEASE = "CollecTRI2"  # a later CollecTRI release in the same download, with no sign decision (9.72)
SIGN_DECISION_FROM_A_PAPER = "PMID"
SIGN_COLUMNS = ("sign_records", "sign_pubmed_ids", "sign_basis")
DISPUTE_COLUMN = "sign_disputed"
DEFAULT_SIGNALLING_EVIDENCES = Path("data/raw/omnipath/omnipath_interactions_evidences.json")
DEFAULT_TRANSCRIPTION_EVIDENCES = Path("data/raw/omnipath/collectri_interactions_evidences.json")


def read_evidences(path: Path) -> list[dict]:
    """An OmniPath interactions download with fields=evidences and format=json: a list of interactions."""
    rows = json.loads(Path(path).read_text())
    if not isinstance(rows, list):
        raise ValueError(f"{path} is not an OmniPath interactions list (format=json)")
    return rows


def references_of(entries, datasets: frozenset[str] | None = None, resource: str | None = None) -> set[str]:
    """The PubMed ids of the entries a curated resource made: text-mining layers left out, and, when given, only the
    named datasets or only the named resource's own entries (via none)."""
    found = set()
    for entry in entries or []:
        if datasets is not None and entry.get("dataset") not in datasets:
            continue
        if entry.get("resource") in TEXT_MINING_RESOURCES or entry.get("via") in TEXT_MINING_RESOURCES:
            continue
        if resource is not None and (entry.get("resource") != resource or entry.get("via") is not None):
            continue
        found.update(str(reference) for reference in entry.get("references") or [] if str(reference).isdigit())
    return found


def sign_decisions_of(entries) -> set[str]:
    """CollecTRI's sign decisions over its own entries; a decision can be a list."""
    decisions = set()
    for entry in entries or []:
        if entry.get("resource") != TRANSCRIPTION_RESOURCE or entry.get("via") is not None:
            continue
        decision = (entry.get("attrs") or {}).get("sign_decision")
        if isinstance(decision, list):
            decisions.update(str(value) for value in decision)
        elif decision:
            decisions.add(str(decision))
    return decisions


def signs_stated_by(evidences: dict, resource: str) -> set[int]:
    """The signs a resource's own entries (via none) state for one interaction: +1 for an entry among the positive
    evidence, -1 for one among the negative, both when the resource states both."""
    stated = set()
    for kind, sign in (("positive", 1), ("negative", -1)):
        if any(entry.get("resource") == resource and entry.get("via") is None for entry in evidences.get(kind) or []):
            stated.add(sign)
    return stated


def joined_signs(signs: set[int]) -> str:
    return ";".join(str(sign) for sign in sorted(signs))


def releases_dispute_the_sign(release_2023_signs: set[int], later_release_signs: set[int]) -> bool:
    """Decision 9.72: True when the two CollecTRI releases state opposite signs for a regulation, or either release
    states both. A release that is silent disputes nothing."""
    if len(release_2023_signs) > 1 or len(later_release_signs) > 1:
        return True
    return bool(release_2023_signs and later_release_signs and release_2023_signs != later_release_signs)


def sign_reference_table(rows: list[dict], family: str) -> pd.DataFrame:
    """(source symbol, target symbol) -> positive and negative PubMed ids (and, for transcription, the sign decisions
    and the signs each CollecTRI release states), complexes expanded to their members, keyed by gene symbols and
    UniProt accessions as edge_evidence.py keys them."""
    if family not in ("signalling", "transcription"):
        raise ValueError("family must be signalling or transcription")
    datasets = SIGNALLING_DATASETS if family == "signalling" else None
    resource = TRANSCRIPTION_RESOURCE if family == "transcription" else None
    table: dict[tuple[str, str], tuple[set[str], set[str], set[str], set[int], set[int]]] = {}
    for row in rows:
        evidences = row.get("evidences") or {}
        positive = references_of(evidences.get("positive"), datasets, resource)
        negative = references_of(evidences.get("negative"), datasets, resource)
        decisions, release_2023_signs, later_release_signs = set(), set(), set()
        if family == "transcription":
            for kind in ("positive", "negative", "directed", "undirected"):
                decisions |= sign_decisions_of(evidences.get(kind))
            release_2023_signs = signs_stated_by(evidences, TRANSCRIPTION_RESOURCE)
            later_release_signs = signs_stated_by(evidences, LATER_TRANSCRIPTION_RELEASE)
        key_sets = [(members(row.get("source_genesymbol", "")), members(row.get("target_genesymbol", "")))]
        if row.get("source") and row.get("target"):
            key_sets.append((members(row["source"]), members(row["target"])))
        for sources, targets in key_sets:
            for source in sources:
                for target in targets:
                    if source == target:
                        continue
                    kept = table.setdefault((source, target), (set(), set(), set(), set(), set()))
                    for kept_set, found in zip(kept, (positive, negative, decisions, release_2023_signs, later_release_signs)):
                        kept_set |= found
    return pd.DataFrame([{"family": family, "source_symbol": source, "target_symbol": target, "positive_pubmed_ids": joined_pubmed_ids(positive),
                          "negative_pubmed_ids": joined_pubmed_ids(negative), "sign_decisions": "+".join(sorted(decisions)),
                          "release_2023_signs": joined_signs(release_2023_signs), "later_release_signs": joined_signs(later_release_signs)}
                         for (source, target), (positive, negative, decisions, release_2023_signs, later_release_signs) in sorted(table.items())],
                        columns=["family", "source_symbol", "target_symbol", "positive_pubmed_ids", "negative_pubmed_ids", "sign_decisions",
                                 "release_2023_signs", "later_release_signs"])


def split_signs(cell) -> set[int]:
    return {int(value) for value in str(cell or "").split(";") if value.strip()}


BASIS_UNSIGNED = "unsigned"
BASIS_NOT_FOUND = "not in the evidences downloads"


def sign_evidence_columns(edges: pd.DataFrame, families: pd.Series, sign_tables: pd.DataFrame) -> pd.DataFrame:
    """Per graph edge, in edge order: the PubMed ids behind its stored sign (sign_pubmed_ids, sign_records) and how the
    sign was decided (sign_basis). The sign is the edge's own sign column (activates +1, inhibits -1, a regulation
    either). sign_basis is empty outside the signalling and transcription families, BASIS_UNSIGNED for an edge of sign
    0 (a binding), BASIS_NOT_FOUND for a signed edge the downloads do not hold (records NaN: unknown, not zero),
    "references" or "no reference" for a signalling edge found there, and CollecTRI's sign decisions joined by "+" for
    a regulation. sign_disputed is True only for a regulation whose two CollecTRI releases dispute its sign
    (releases_dispute_the_sign)."""
    by_key = {(row.family, row.source_symbol, row.target_symbol): row for row in sign_tables.itertuples(index=False)} if len(sign_tables) else {}
    records, identifiers, bases, disputed = [], [], [], []
    signs = edges["sign"].fillna(0.0).to_numpy(dtype=float) if "sign" in edges.columns else np.zeros(len(edges))
    for position, edge in enumerate(edges.itertuples(index=False)):
        family = families.iloc[position]
        if family not in ("signalling", "transcription"):
            records.append(np.nan)
            identifiers.append("")
            bases.append("")
            disputed.append(False)
            continue
        sign = signs[position]
        if sign == 0:
            records.append(np.nan)
            identifiers.append("")
            bases.append(BASIS_UNSIGNED)
            disputed.append(False)
            continue
        source, target = symbol_of(edge.source_id), symbol_of(edge.target_id)
        row = by_key.get((family, source, target)) or by_key.get((family, accession_of(edge.source_id) or source, accession_of(edge.target_id) or target))
        if row is None:
            records.append(np.nan)
            identifiers.append("")
            bases.append(BASIS_NOT_FOUND)
            disputed.append(False)
            continue
        behind_the_sign = split_pubmed_ids(row.positive_pubmed_ids if sign > 0 else row.negative_pubmed_ids)
        records.append(float(len(behind_the_sign)))
        identifiers.append(joined_pubmed_ids(behind_the_sign))
        if family == "transcription":
            bases.append(str(row.sign_decisions or "none stated"))
            disputed.append(releases_dispute_the_sign(split_signs(getattr(row, "release_2023_signs", "")),
                                                      split_signs(getattr(row, "later_release_signs", ""))))
        else:
            bases.append("references" if behind_the_sign else "no reference")
            disputed.append(False)
    return pd.DataFrame({"sign_records": records, "sign_pubmed_ids": identifiers, "sign_basis": bases, DISPUTE_COLUMN: disputed}, index=edges.index)


def sign_is_evidenced(evidence: pd.DataFrame, minimum_sign_records: int) -> np.ndarray:
    """True where a signed signalling or transcription edge's sign rests on at least minimum_sign_records PubMed ids and,
    for a regulation, CollecTRI decided the sign from a paper; a signed edge the downloads do not hold fails, since the
    rule asks for evidence that has been seen. True for every edge the rule does not concern: another family, or an
    unsigned edge such as a binding."""
    for column in SIGN_COLUMNS:
        if column not in evidence.columns:
            raise ValueError("the edge evidence table has no sign columns; rebuild it with experiments/build_edge_evidence.py and the OmniPath evidences "
                             "downloads (scripts/fetch_raw_sources.sh)")
    basis = evidence["sign_basis"].fillna("").astype(str)
    concerned = (evidence["family"].isin(["signalling", "transcription"]) & (basis != "") & (basis != BASIS_UNSIGNED)).to_numpy()
    enough = evidence["sign_records"].fillna(-1).to_numpy() >= minimum_sign_records
    is_regulation = evidence["family"].eq("transcription").to_numpy()
    from_a_paper = basis.str.split("+").map(lambda decisions: SIGN_DECISION_FROM_A_PAPER in decisions).to_numpy()
    return np.where(concerned, enough & (~is_regulation | from_a_paper), True)


def sign_is_undisputed(evidence: pd.DataFrame) -> np.ndarray:
    """True for every edge except a regulation whose two CollecTRI releases dispute its sign (decision 9.72)."""
    if DISPUTE_COLUMN not in evidence.columns:
        raise ValueError(f"the edge evidence table has no {DISPUTE_COLUMN} column; rebuild it with experiments/build_edge_evidence.py and the CollecTRI "
                         "evidences download (scripts/fetch_raw_sources.sh)")
    return ~evidence[DISPUTE_COLUMN].fillna(False).astype(bool).to_numpy()
