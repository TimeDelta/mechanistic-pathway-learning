"""Match SIDER drug names to MeSH chemical identifiers, the join key of PubTator3 and CTD chemicals (design section 4.2, E3).

First pass, by name only: a SIDER drug (drug_names.tsv, STITCH flat id and name) matches a MeSH chemical when the
MeSH descriptor name, or an entry term that the PubTator3 autocomplete reports as the matched synonym, equals the
SIDER name case-insensitively. The method behind every match is recorded (chemical_match_method) so a later pass
can replace name matching with a structure-based mapping (PubChem CID to MeSH through UniChem or the MeSH
registry numbers) and measure what changed. Salts, combination products and names that MeSH spells differently
("gamma-aminobutyric" in SIDER against "gamma-Aminobutyric Acid" in MeSH) do not match in this pass and are
counted as unmatched, never guessed.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

MATCH_METHODS: tuple[str, ...] = ("autocomplete_name", "autocomplete_synonym", "mesh_label_lookup", "ctd_chemical_name")
SYNONYM_MATCH_PATTERN = re.compile(r"<m>(.*?)</m>")


@dataclass(frozen=True)
class ChemicalMatch:
    stitch_flat_id: str
    drug_name: str
    chemical_mesh_id: str  # D- or C-number without the "MESH:" prefix
    chemical_mesh_name: str
    match_method: str


def normalize_name(name: str) -> str:
    return " ".join(name.strip().lower().split())


def load_sider_drug_names(sider_directory: Path) -> dict[str, str]:
    """STITCH flat id -> drug name from drug_names.tsv (one name per drug)."""
    names: dict[str, str] = {}
    with open(Path(sider_directory) / "drug_names.tsv", newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2 and row[0]:
                names[row[0]] = row[1]
    return names


def drug_names_by_normalized_name(drug_names: dict[str, str]) -> dict[str, list[tuple[str, str]]]:
    """Normalized name -> [(STITCH flat id, name)]; a name shared by two STITCH ids maps to both, which the caller reports."""
    index: dict[str, list[tuple[str, str]]] = {}
    for stitch_flat_id, drug_name in drug_names.items():
        index.setdefault(normalize_name(drug_name), []).append((stitch_flat_id, drug_name))
    return index


def strip_mesh_prefix(identifier: str) -> str:
    return identifier[5:] if identifier.upper().startswith("MESH:") else identifier


def match_autocomplete_candidates(stitch_flat_id: str, drug_name: str, candidates: list[dict]) -> ChemicalMatch | None:
    """The first chemical candidate of a PubTator3 autocomplete response whose name, or whose reported matched synonym,
    equals the drug name case-insensitively. Candidates of other biotypes (a gene named like a drug) are skipped."""
    wanted = normalize_name(drug_name)
    for candidate in candidates:
        if candidate.get("biotype") != "chemical" or candidate.get("db") != "ncbi_mesh":
            continue
        candidate_name = str(candidate.get("name", ""))
        if normalize_name(candidate_name) == wanted:
            return ChemicalMatch(stitch_flat_id, drug_name, strip_mesh_prefix(str(candidate["db_id"])), candidate_name, "autocomplete_name")
    for candidate in candidates:
        if candidate.get("biotype") != "chemical" or candidate.get("db") != "ncbi_mesh":
            continue
        for synonym in SYNONYM_MATCH_PATTERN.findall(str(candidate.get("match", ""))):
            if normalize_name(synonym) == wanted:
                return ChemicalMatch(stitch_flat_id, drug_name, strip_mesh_prefix(str(candidate["db_id"])), str(candidate.get("name", "")), "autocomplete_synonym")
    return None


def match_mesh_label(chemical_mesh_id: str, mesh_label: str, names_index: dict[str, list[tuple[str, str]]]) -> list[ChemicalMatch]:
    """Every SIDER drug whose name equals a MeSH chemical's preferred label case-insensitively (method mesh_label_lookup)."""
    hits = names_index.get(normalize_name(mesh_label), [])
    return [ChemicalMatch(stitch_flat_id, drug_name, strip_mesh_prefix(chemical_mesh_id), mesh_label, "mesh_label_lookup") for stitch_flat_id, drug_name in hits]


def match_ctd_chemical_name(chemical_mesh_id: str, chemical_name: str, names_index: dict[str, list[tuple[str, str]]]) -> list[ChemicalMatch]:
    """Every SIDER drug whose name equals a CTD ChemicalName case-insensitively (method ctd_chemical_name)."""
    hits = names_index.get(normalize_name(chemical_name), [])
    return [ChemicalMatch(stitch_flat_id, drug_name, strip_mesh_prefix(chemical_mesh_id), chemical_name, "ctd_chemical_name") for stitch_flat_id, drug_name in hits]


def chemical_matches_by_mesh_id(matches: list[ChemicalMatch]) -> dict[str, list[ChemicalMatch]]:
    index: dict[str, list[ChemicalMatch]] = {}
    for match in matches:
        index.setdefault(match.chemical_mesh_id, []).append(match)
    return index


MATCH_METHOD_PRIORITY: dict[str, int] = {method: rank for rank, method in enumerate(MATCH_METHODS)}


def pubchem_cid_of_stitch_flat_id(stitch_flat_id: str) -> int:
    """3386 from "CID100003386" (STITCH flat ids are "CID1" followed by the zero-padded PubChem CID)."""
    digits = stitch_flat_id[4:] if stitch_flat_id.upper().startswith("CID1") else stitch_flat_id[3:]
    return int(digits) if digits.isdigit() else 0


def collapse_matches_to_one_drug_per_mesh_id(matches: list[ChemicalMatch], preferred_stitch_ids: set[str] | None = None) -> tuple[list[ChemicalMatch], dict[str, list[str]]]:
    """One SIDER drug per MeSH chemical id, so a paper is never counted under a second STITCH id for the same chemical
    (SIDER lists brand names, salts and code names as separate drugs: olanzapine and Zyprexa, rizatriptan and MK-462).
    Order of preference: the match method (an exact descriptor name beats a synonym beats a label lookup beats a CTD
    name), then a STITCH id in preferred_stitch_ids (the ids the assembled table or the ChEMBL parent map already use),
    then the lowest PubChem CID. Returns the kept matches and, per MeSH id that lost entries, the dropped STITCH ids."""
    preferred = preferred_stitch_ids or set()
    kept: list[ChemicalMatch] = []
    collapsed: dict[str, list[str]] = {}
    for mesh_id, candidates in chemical_matches_by_mesh_id(matches).items():
        ranked = sorted(candidates, key=lambda match: (MATCH_METHOD_PRIORITY.get(match.match_method, len(MATCH_METHODS)), 0 if match.stitch_flat_id in preferred else 1, pubchem_cid_of_stitch_flat_id(match.stitch_flat_id)))
        kept.append(ranked[0])
        if len(ranked) > 1:
            collapsed[mesh_id] = [match.stitch_flat_id for match in ranked[1:]]
    return kept, collapsed
