"""Symptom constructs at two levels: the coarse crosswalk as registered, and a fine crosswalk derived from it.

Decision 9.58 (the owner, 16 October 2026): "Finely grained symptoms is the goal" and merging is to be an output of the
model, not an input. The coarse crosswalk (docs/symptom_crosswalk.csv) fuses several HPO roots into one construct in
seven rows. The fine level has one construct per HPO root, named by the root's own HPO name, so that no construct is
authored here: psychosis's four roots become four constructs, parkinsonism's six become six, and a construct with one
root keeps its coarse name. New constructs the coarse level never had (intellectual disability) come from an additions
table with the same columns as the crosswalk.

The MedDRA preferred terms and MeSH descriptors of a coarse row have to be shared out among its fine constructs. The
rule is deterministic and reviewable rather than clever: a term goes to the fine construct whose HPO name shares the
most word stems with it; a term that matches none goes to the coarse row's namesake (the root whose name slugs to the
coarse name, else the first root); an override table (term -> fine construct) wins over both. The generated page
lists every assignment so the owner can move any of them by adding an override.

The IMPC crosswalk is re-keyed the same way: a coarse row's proxy goes to the fine construct its MP term names match,
else to the namesake; every other fine construct of that coarse row takes a row with validity "none", because a proxy
stated for a coarse construct was not stated for its parts.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from mechanistic_pathway_learning.evidence.load_monogenic_phenotype_annotations import (
    descendants_of,
    load_hpo_is_a_parents_from_obo,
    load_hpo_term_names_from_obo,
)

COARSE_CROSSWALK_PATH = Path("docs/symptom_crosswalk.csv")
FINE_CROSSWALK_PATH = Path("docs/symptom_crosswalk_fine.csv")
ADDITIONS_PATH = Path("configs/symptom_construct_additions.csv")
TERM_OVERRIDES_PATH = Path("configs/fine_construct_term_overrides.csv")
IMPC_COARSE_PATH = Path("configs/impc_symptom_crosswalk.csv")
IMPC_FINE_PATH = Path("configs/impc_symptom_crosswalk_fine.csv")

CROSSWALK_COLUMNS = ("target_symptom", "hpo_ids", "hpo_names", "excluded_hpo_ids", "meddra_preferred_terms", "questionnaire_item", "rdoc_domain",
                     "identifier_status", "notes", "mesh_descriptors", "mesh_descriptor_level")
FINE_COLUMNS = CROSSWALK_COLUMNS + ("coarse_symptom",)
STOP_WORDS = {"of", "the", "and", "or", "in", "to", "a", "an", "disorder", "behaviour", "behavior", "behaviours", "behaviors", "abnormal", "increased", "decreased"}


def slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def stem(word: str) -> str:
    """A crude stem for matching a term to a construct name: lower case, no plural, no common suffixes."""
    word = word.lower()
    for suffix in ("ations", "ation", "ness", "ity", "ies", "ing", "ous", "ic", "al", "es", "s"):
        if word.endswith(suffix) and len(word) - len(suffix) >= 4:
            return word[: -len(suffix)]
    return word


BRITISH_SPELLINGS = {"behaviour": "behavior", "behaviours": "behaviors"}


def stems(text: str, keep_stop_words: bool = False) -> set[str]:
    words = [BRITISH_SPELLINGS.get(word, word) for word in re.findall(r"[a-z]+", text.lower())]
    return {stem(word) for word in words if (keep_stop_words or word not in STOP_WORDS) and len(word) > 2}


def phrase(text: str) -> frozenset[str]:
    """A phrase as its set of words, so that 'Hallucination, auditory' and 'Auditory hallucination' are one phrase."""
    return frozenset(stem(BRITISH_SPELLINGS.get(word, word)) for word in re.findall(r"[a-z]+", text.lower()))


def load_hpo_synonyms_from_obo(obo_path: Path) -> dict[str, list[str]]:
    """term id -> its synonym strings (every scope), from the synonym: lines of an OBO file."""
    synonyms: dict[str, list[str]] = {}
    current_term: str | None = None
    inside_term_stanza = False
    with open(obo_path, encoding="utf-8") as obo_file:
        for raw_line in obo_file:
            line = raw_line.strip()
            if line.startswith("["):
                inside_term_stanza = line == "[Term]"
                current_term = None
            elif inside_term_stanza and line.startswith("id: HP:"):
                current_term = line[len("id: "):]
                synonyms.setdefault(current_term, [])
            elif line.startswith("synonym:") and current_term is not None:
                match = re.match(r'synonym: "(.*?)"', line)
                if match:
                    synonyms[current_term].append(match.group(1))
    return synonyms


@dataclass(frozen=True)
class HpoVocabulary:
    """What the ontology calls each construct: the phrases (name and synonyms) of a root and its kept descendants, and
    the word stems of the root's own name and synonyms. Built once per obo file and consulted per term."""

    parents_by_term: dict[str, set[str]]
    names_by_term: dict[str, str]
    synonyms_by_term: dict[str, list[str]]

    @classmethod
    def from_obo(cls, obo_path: Path) -> "HpoVocabulary":
        return cls(load_hpo_is_a_parents_from_obo(obo_path), load_hpo_term_names_from_obo(obo_path), load_hpo_synonyms_from_obo(obo_path))

    def phrases_of(self, root_id: str, excluded_ids: list[str]) -> set[frozenset[str]]:
        kept = {root_id} | descendants_of(root_id, self.parents_by_term)
        for excluded in excluded_ids:
            kept -= {excluded} | descendants_of(excluded, self.parents_by_term)
        phrases: set[frozenset[str]] = set()
        for term in kept:
            for text in [self.names_by_term.get(term, "")] + self.synonyms_by_term.get(term, []):
                if text:
                    phrases.add(phrase(text))
        return phrases

    def ancestors_of(self, term_id: str) -> set[str]:
        """Every ancestor of a term through is_a, the term itself excluded."""
        ancestors, frontier = set(), list(self.parents_by_term.get(term_id, ()))
        while frontier:
            parent = frontier.pop()
            if parent not in ancestors:
                ancestors.add(parent)
                frontier.extend(self.parents_by_term.get(parent, ()))
        return ancestors

    def root_stems_of(self, root_id: str) -> set[str]:
        texts = [self.names_by_term.get(root_id, "")] + self.synonyms_by_term.get(root_id, [])
        return set().union(*[stems(text) for text in texts if text]) if texts else set()


@dataclass(frozen=True)
class FineConstruct:
    name: str
    coarse: str
    hpo_id: str
    hpo_name: str
    is_namesake: bool


def split_cell(cell: str | None) -> list[str]:
    return [part.strip() for part in (cell or "").split(";") if part.strip()]


def read_rows(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def fine_constructs_of(row: dict[str, str]) -> list[FineConstruct]:
    """One construct per HPO root of a coarse row; a single-root row keeps its coarse name. A row without an HPO root
    (psychomotor_retardation, abnormal_dreams, catatonia) is one construct under its coarse name."""
    coarse = row["target_symptom"]
    root_ids, root_names = split_cell(row.get("hpo_ids")), split_cell(row.get("hpo_names"))
    if len(root_ids) != len(root_names):
        raise ValueError(f"{coarse}: {len(root_ids)} hpo_ids but {len(root_names)} hpo_names")
    if len(root_ids) <= 1:
        return [FineConstruct(coarse, coarse, root_ids[0] if root_ids else "", root_names[0] if root_names else "", True)]
    names = [slug(name) for name in root_names]
    namesake_index = names.index(coarse) if coarse in names else 0
    constructs = []
    for index, (root_id, root_name, name) in enumerate(zip(root_ids, root_names, names)):
        constructs.append(FineConstruct(name, coarse, root_id, root_name, index == namesake_index))
    if len({construct.name for construct in constructs}) != len(constructs):
        raise ValueError(f"{coarse}: two HPO roots slug to one construct name: {names}")
    return constructs


def read_term_overrides(path: Path | None) -> dict[str, str]:
    """term (lower case) -> fine construct; the file has columns term, fine_construct and an optional note."""
    if path is None or not Path(path).exists():
        return {}
    return {row["term"].strip().lower(): row["fine_construct"].strip() for row in read_rows(Path(path)) if row.get("term", "").strip()}


def assign_term(term: str, constructs: list[FineConstruct], overrides: dict[str, str], vocabulary: HpoVocabulary | None = None,
                excluded_ids: list[str] | None = None) -> tuple[str, str]:
    """(fine construct, how). In order: an override; the construct whose HPO expansion names the term (name or synonym,
    as a set of word stems, so word order and plurals do not matter); the construct whose root name and synonyms share
    the most word stems with the term, after the stems every construct of the row shares are dropped (they carry no
    information about which part); else the row's namesake."""
    override = overrides.get(term.strip().lower())
    names = {construct.name for construct in constructs}
    if override is not None:
        if override not in names:
            raise ValueError(f"override sends {term!r} to {override!r}, which is not a construct of {constructs[0].coarse}: {sorted(names)}")
        return override, "override"
    namesake = next(construct for construct in constructs if construct.is_namesake)
    if len(constructs) == 1:
        return namesake.name, "single construct"
    if vocabulary is not None:
        term_phrase = phrase(term)
        holders = [construct for construct in constructs if construct.hpo_id and term_phrase in vocabulary.phrases_of(construct.hpo_id, excluded_ids or [])]
        # a root that is an ancestor of another holder's root holds the term only through that descendant; the most specific holder decides
        holders = [holder for holder in holders if not any(other is not holder and holder.hpo_id in vocabulary.ancestors_of(other.hpo_id) for other in holders)]
        if len(holders) == 1:
            return holders[0].name, "hpo vocabulary"
    for keep_stop_words in (False, True):  # a second pass keeps words such as behavior, which the first drops as generic
        root_stems = {construct.name: (vocabulary.root_stems_of(construct.hpo_id) if vocabulary is not None and construct.hpo_id else set())
                      | stems(construct.hpo_name or construct.name.replace("_", " "), keep_stop_words) for construct in constructs}
        shared = set.intersection(*root_stems.values()) if root_stems else set()
        term_stems = stems(term, keep_stop_words) - shared
        best, best_overlap, tied = None, 0, False
        for construct in constructs:
            overlap = len(term_stems & (root_stems[construct.name] - shared))
            if overlap > best_overlap:
                best, best_overlap, tied = construct, overlap, False
            elif overlap == best_overlap and overlap > 0:
                tied = True
        if best is not None and not tied:
            return best.name, "stem match"
    return namesake.name, "namesake"


def fine_rows_of(row: dict[str, str], overrides: dict[str, str], vocabulary: HpoVocabulary | None = None) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """The fine crosswalk rows of one coarse row and the term assignments made (term, construct, how)."""
    constructs = fine_constructs_of(row)
    assignments = []
    terms_by_construct: dict[str, list[str]] = {construct.name: [] for construct in constructs}
    excluded_ids = split_cell(row.get("excluded_hpo_ids"))
    for term in split_cell(row.get("meddra_preferred_terms")):
        name, how = assign_term(term, constructs, overrides, vocabulary, excluded_ids)
        terms_by_construct[name].append(term)
        assignments.append({"coarse_symptom": row["target_symptom"], "term": term, "kind": "meddra", "fine_construct": name, "how": how})
    namesake = next(construct for construct in constructs if construct.is_namesake)
    rows = []
    for construct in constructs:
        # a root that is an ancestor of a sibling construct's root would expand over the sibling's terms too; the sibling's root is excluded from it
        descending_siblings = [other.hpo_id for other in constructs if vocabulary is not None and other is not construct and other.hpo_id
                               and construct.hpo_id in vocabulary.ancestors_of(other.hpo_id)]
        derived = {column: "" for column in FINE_COLUMNS}
        derived.update({
            "target_symptom": construct.name, "coarse_symptom": construct.coarse, "hpo_ids": construct.hpo_id, "hpo_names": construct.hpo_name,
            "excluded_hpo_ids": ";".join(excluded_ids + descending_siblings), "meddra_preferred_terms": ";".join(terms_by_construct[construct.name]),
            "rdoc_domain": row.get("rdoc_domain", ""), "identifier_status": row.get("identifier_status", ""),
            "notes": (row.get("notes", "") if len(constructs) == 1 else f"derived from {construct.coarse} (decision 9.58); the coarse row's notes: {row.get('notes', '')}".rstrip(": ")),
        })
        if construct is namesake:  # the questionnaire item and the MeSH descriptors are stated at the coarse level and stay with its namesake
            derived["questionnaire_item"] = row.get("questionnaire_item", "")
            derived["mesh_descriptors"] = row.get("mesh_descriptors", "")
            derived["mesh_descriptor_level"] = row.get("mesh_descriptor_level", "")
        rows.append(derived)
    return rows, assignments


def build_fine_crosswalk(coarse_rows: list[dict[str, str]], additions: list[dict[str, str]], overrides: dict[str, str],
                         vocabulary: HpoVocabulary | None = None) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Every fine row of every coarse row, then the additions (each a construct of its own), and the term assignments."""
    rows, assignments = [], []
    for row in coarse_rows:
        derived, made = fine_rows_of(row, overrides, vocabulary)
        rows.extend(derived)
        assignments.extend(made)
    for addition in additions:
        if len(split_cell(addition.get("hpo_ids"))) > 1:
            raise ValueError(f"addition {addition['target_symptom']} names several HPO roots; add one construct per root")
        derived = {column: addition.get(column, "") for column in CROSSWALK_COLUMNS}
        derived["coarse_symptom"] = addition["target_symptom"]
        rows.append(derived)
    names = [row["target_symptom"] for row in rows]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"fine constructs named twice: {duplicates}")
    return rows, assignments


def coarse_of_fine(fine_rows: list[dict[str, str]]) -> dict[str, str]:
    return {row["target_symptom"]: row["coarse_symptom"] for row in fine_rows}


def fine_impc_rows(impc_rows: list[dict[str, str]], coarse_rows: list[dict[str, str]], additions: list[dict[str, str]]) -> list[dict[str, str]]:
    """The IMPC crosswalk re-keyed to the fine constructs: a coarse row's proxy goes to the fine construct its MP term
    names match, else to the namesake; the other fine constructs of that coarse row take a row saying no proxy was
    stated for them; an added construct takes a 'none' row until one is curated."""
    constructs_by_coarse = {row["target_symptom"]: fine_constructs_of(row) for row in coarse_rows}
    output = []
    for impc_row in impc_rows:
        constructs = constructs_by_coarse.get(impc_row["target_symptom"])
        if constructs is None:
            raise ValueError(f"IMPC crosswalk names {impc_row['target_symptom']!r}, which the coarse crosswalk does not hold")
        if len(constructs) == 1:
            output.append({**impc_row, "target_symptom": constructs[0].name})
            continue
        term_stems = stems(impc_row.get("mp_term_names", ""))
        best, best_overlap = None, 0
        for construct in constructs:
            overlap = len(term_stems & stems(construct.hpo_name))
            if overlap > best_overlap:
                best, best_overlap = construct, overlap
        chosen = best or next(construct for construct in constructs if construct.is_namesake)
        for construct in constructs:
            if construct is chosen:
                output.append({**impc_row, "target_symptom": construct.name})
            else:
                output.append({"target_symptom": construct.name, "relation": impc_row["relation"], "mp_term_ids": "", "mp_term_names": "", "negative_procedure_codes": "",
                               "validity": "none", "note": f"no mouse proxy stated at the fine level; the coarse row's proxy went to {chosen.name} (decision 9.58)"})
    seen = {(row["target_symptom"], row["relation"]) for row in output}
    for addition in additions:
        if (addition["target_symptom"], "induces") not in seen:
            output.append({"target_symptom": addition["target_symptom"], "relation": "induces", "mp_term_ids": "", "mp_term_names": "", "negative_procedure_codes": "",
                           "validity": "none", "note": "added construct; no mouse proxy curated yet"})
    return output


def write_rows(path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


IMPC_COLUMNS = ("target_symptom", "relation", "mp_term_ids", "mp_term_names", "negative_procedure_codes", "validity", "note")
