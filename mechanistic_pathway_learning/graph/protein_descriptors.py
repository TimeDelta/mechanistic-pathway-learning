"""Protein descriptors: ESM-2 embeddings reduced by reduced-rank regression onto function annotations that do not come
from phenotypes (design section 5.2, node descriptors).

Why not principal components: PCA keeps the directions along which the embeddings vary most, and variance is not
relevance. Here the kept directions are the ones that best predict what a protein does and where it sits, fitted over
the reviewed human proteome (about 20,000 proteins, none carrying a symptom label), so no symptom label enters and the
directions are the same for every split.

Targets, one multi-hot block each, every block standardised and scaled by 1 / sqrt(its number of columns) so that the
blocks weigh equally in the regression:
  - EC numbers to the second level (1.14, 2.7, ...), from UniProt;
  - GO molecular function and GO cellular component terms, propagated to their ancestors (is_a, part_of), from
    annotations whose evidence code is not phenotype-derived (IMP, IGI, HMP and HGI are dropped; biological process,
    which holds the behaviour terms, is not used) and that are not NOT-qualified;
  - UniProt subcellular locations and membrane topologies;
  - Pfam families.
A column is kept when between MINIMUM_PROTEINS_PER_TARGET proteins and MAXIMUM_TARGET_PREVALENCE of them carry it.
A protein without any GO function, GO component or location annotation is unannotated in that block, not negative for
every term: its targets in the block are set to the column mean (0 after standardisation). A missing EC number or Pfam
family is a negative (reviewed entries list them when they apply).

Reduced-rank regression (A. J. Izenman, J. Multivariate Analysis 5, 248-264, 1975): with standardised embeddings X and
targets Y, the ridge coefficients B = (X'X + alpha I)^-1 X'Y give fitted values X B; the top k right singular vectors
V_k of the fitted values span the k-dimensional annotation space best predicted from X, and the descriptors are
Z = X B V_k, standardised. alpha and k are chosen by cross-validated prediction of the targets (held-out R^2): alpha at
full rank, then the smallest k whose held-out R^2 reaches SUFFICIENT_FRACTION_OF_BEST of the best on the grid.

ICA rotation (optional, for reading): FastICA within the k-dimensional subspace; for the model's learned linear input
layer it is the same subspace, so it changes interpretability, not what the model can represent. Components are defined
up to sign and order, so they are matched across seeds and their stability is reported.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd

EXCLUDED_GO_EVIDENCE_CODES = frozenset({"IMP", "IGI", "HMP", "HGI"})
GO_ASPECTS = {"F": "go_function", "C": "go_component"}
MINIMUM_PROTEINS_PER_TARGET = 20
MAXIMUM_TARGET_PREVALENCE = 0.25
SUFFICIENT_FRACTION_OF_BEST = 0.95
RANK_GRID = (2, 4, 8, 12, 16, 24, 32, 48, 64)
RELATIVE_RIDGE_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)  # alpha = value * number of training proteins


def read_go_ancestors(go_obo_path) -> dict[str, set[str]]:
    """GO term -> its ancestors through is_a and part_of (the term itself included); obsolete terms map to themselves."""
    parents: dict[str, set[str]] = defaultdict(set)
    current = None
    with open(go_obo_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line in ("[Term]", "[Typedef]"):
                current = None
            elif line.startswith("id: GO:"):
                current = line[4:]
            elif current is None:
                continue
            elif line.startswith("is_a: GO:"):
                parents[current].add(line[6:16])
            elif line.startswith("relationship: part_of GO:"):
                parents[current].add(line[22:32])
    ancestors: dict[str, set[str]] = {}

    def ancestors_of(term: str) -> set[str]:
        if term not in ancestors:
            ancestors[term] = {term}  # placeholder against cycles
            result = {term}
            for parent in parents.get(term, ()):
                result |= ancestors_of(parent)
            ancestors[term] = result
        return ancestors[term]

    for term in list(parents):
        ancestors_of(term)
    return defaultdict(set, ancestors)


def go_annotations(gaf_path, go_ancestors: dict[str, set[str]]) -> dict[str, dict[str, set[str]]]:
    """aspect block name -> accession -> propagated GO terms, without phenotype-derived evidence or NOT qualifiers."""
    annotations: dict[str, dict[str, set[str]]] = {block: defaultdict(set) for block in GO_ASPECTS.values()}
    with gzip.open(gaf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != "UniProtKB" or fields[8] not in GO_ASPECTS:
                continue
            if fields[6] in EXCLUDED_GO_EVIDENCE_CODES or fields[3].startswith("NOT"):
                continue
            annotations[GO_ASPECTS[fields[8]]][fields[1]] |= go_ancestors.get(fields[4]) or {fields[4]}
    return annotations


def subcellular_locations(location_text: str) -> set[str]:
    """UniProt 'Subcellular location [CC]' text -> location and topology terms, without evidence, isoform labels and notes."""
    if not isinstance(location_text, str):
        return set()
    text = re.sub(r"\{[^}]*\}", "", location_text)
    text = re.sub(r"\[[^\]]*\]:", "", text)
    terms = set()
    for statement in text.split("SUBCELLULAR LOCATION:"):
        statement = statement.split("Note=")[0]
        for part in re.split(r"[.;]", statement):
            part = part.strip(" ,")
            if part:
                terms.add(part.lower())
    return terms


BLOCKS_WITH_MISSING_ANNOTATION = ("go_function", "go_component", "location")


def annotation_target_matrix(uniprot_table: pd.DataFrame, go_by_block: dict[str, dict[str, set[str]]]) -> tuple[pd.DataFrame, dict[str, list[str]], pd.DataFrame]:
    """Accession-indexed 0/1 target matrix, the columns of each block after the prevalence filter, and an
    accession-by-block table of whether the protein is annotated in the block at all."""
    accessions = uniprot_table["Entry"].tolist()
    sets_by_block: dict[str, dict[str, set[str]]] = {
        "ec_level_2": {accession: {".".join(ec.strip().split(".")[:2]) for ec in str(ec_field).split(";") if ec.strip() and ec.strip()[0].isdigit()}
                       for accession, ec_field in zip(accessions, uniprot_table["EC number"].fillna(""))},
        **{block: {accession: annotations.get(accession, set()) for accession in accessions} for block, annotations in go_by_block.items()},
        "location": {accession: subcellular_locations(text) for accession, text in zip(accessions, uniprot_table["Subcellular location [CC]"])},
        "pfam": {accession: {pfam for pfam in str(field).split(";") if pfam} for accession, field in zip(accessions, uniprot_table["Pfam"].fillna(""))},
    }
    observed = pd.DataFrame({block: [bool(sets[accession]) or block not in BLOCKS_WITH_MISSING_ANNOTATION for accession in accessions]
                             for block, sets in sets_by_block.items()}, index=pd.Index(accessions, name="accession"))
    columns, block_columns = {}, {}
    for block, sets in sets_by_block.items():
        counts: dict[str, int] = defaultdict(int)
        for terms in sets.values():
            for term in terms:
                counts[term] += 1
        kept = sorted(term for term, count in counts.items() if MINIMUM_PROTEINS_PER_TARGET <= count <= MAXIMUM_TARGET_PREVALENCE * len(accessions))
        block_columns[block] = [f"{block}:{term}" for term in kept]
        for term in kept:
            columns[f"{block}:{term}"] = np.fromiter((term in sets[accession] for accession in accessions), dtype=np.float32, count=len(accessions))
    return pd.DataFrame(columns, index=pd.Index(accessions, name="accession")), block_columns, observed


def weighted_standardised_targets(targets: pd.DataFrame, block_columns: dict[str, list[str]], observed: pd.DataFrame | None = None) -> np.ndarray:
    """Each column standardised over the proteins annotated in its block (unannotated ones get 0, the mean) and
    scaled by 1 / sqrt(number of columns in its block)."""
    values = targets.to_numpy(dtype=np.float64).copy()
    for block, columns in block_columns.items():
        if not columns:
            continue
        positions = [targets.columns.get_loc(column) for column in columns]
        rows = observed[block].to_numpy() if observed is not None else np.ones(len(targets), dtype=bool)
        block_values = values[np.ix_(rows, positions)]
        mean, deviation = block_values.mean(axis=0), block_values.std(axis=0)
        standardised = np.zeros((len(targets), len(positions)))
        standardised[rows] = (block_values - mean) / np.where(deviation > 0, deviation, 1.0)
        values[:, positions] = standardised / np.sqrt(len(positions))
    return values


@dataclass
class ReducedRankRegression:
    embedding_mean: np.ndarray
    embedding_scale: np.ndarray
    projection: np.ndarray  # [embedding dim, rank]: Z = standardised X @ projection
    descriptor_mean: np.ndarray
    descriptor_scale: np.ndarray

    def transform(self, embeddings: np.ndarray) -> np.ndarray:
        scores = ((embeddings - self.embedding_mean) / self.embedding_scale) @ self.projection
        return (scores - self.descriptor_mean) / self.descriptor_scale


def ridge_coefficients(standardised_embeddings: np.ndarray, targets: np.ndarray, alpha: float) -> np.ndarray:
    gram = standardised_embeddings.T @ standardised_embeddings
    return np.linalg.solve(gram + alpha * np.eye(gram.shape[0]), standardised_embeddings.T @ targets)


def top_target_directions(standardised_embeddings: np.ndarray, coefficients: np.ndarray, rank: int) -> np.ndarray:
    """V_k: the top right singular vectors of the fitted targets."""
    _, _, right_singular_vectors = np.linalg.svd(standardised_embeddings @ coefficients, full_matrices=False)
    return right_singular_vectors[:rank].T


def fit_reduced_rank_regression(embeddings: np.ndarray, targets: np.ndarray, rank: int, relative_alpha: float) -> ReducedRankRegression:
    mean, scale = embeddings.mean(axis=0), np.where(embeddings.std(axis=0) > 0, embeddings.std(axis=0), 1.0)
    standardised = (embeddings - mean) / scale
    coefficients = ridge_coefficients(standardised, targets, relative_alpha * len(embeddings))
    projection = coefficients @ top_target_directions(standardised, coefficients, rank)
    scores = standardised @ projection
    return ReducedRankRegression(mean, scale, projection, scores.mean(axis=0), np.where(scores.std(axis=0) > 0, scores.std(axis=0), 1.0))


def cross_validated_r_squared(embeddings: np.ndarray, targets: np.ndarray, ranks, relative_alpha: float, folds: int = 5, seed: int = 0) -> dict[int, float]:
    """Held-out R^2 of the targets (pooled over columns) for each rank, from folds over proteins."""
    order = np.random.default_rng(seed).permutation(len(embeddings))
    residual_by_rank, total = defaultdict(float), 0.0
    for fold in range(folds):
        test = order[fold::folds]
        train = np.setdiff1d(order, test)
        mean, scale = embeddings[train].mean(axis=0), np.where(embeddings[train].std(axis=0) > 0, embeddings[train].std(axis=0), 1.0)
        train_x, test_x = (embeddings[train] - mean) / scale, (embeddings[test] - mean) / scale
        target_mean = targets[train].mean(axis=0)
        train_y, test_y = targets[train] - target_mean, targets[test] - target_mean
        coefficients = ridge_coefficients(train_x, train_y, relative_alpha * len(train))
        directions = top_target_directions(train_x, coefficients, max(ranks))
        test_fitted_full = test_x @ coefficients
        total += float((test_y ** 2).sum())
        for rank in ranks:
            directions_k = directions[:, :rank]
            prediction = test_fitted_full @ directions_k @ directions_k.T
            residual_by_rank[rank] += float(((test_y - prediction) ** 2).sum())
    return {rank: 1.0 - residual_by_rank[rank] / total for rank in ranks}


def choose_rank(r_squared_by_rank: dict[int, float], sufficient_fraction: float = SUFFICIENT_FRACTION_OF_BEST) -> int:
    best = max(r_squared_by_rank.values())
    return min(rank for rank, value in r_squared_by_rank.items() if value >= sufficient_fraction * best)


def ica_rotation(descriptors: np.ndarray, seed: int):
    from sklearn.decomposition import FastICA

    return FastICA(n_components=descriptors.shape[1], whiten="unit-variance", random_state=seed, max_iter=2000, tol=1e-5).fit(descriptors)


def matched_component_stability(reference: np.ndarray, other: np.ndarray) -> float:
    """Mean absolute correlation of components matched one to one (Hungarian assignment on |correlation|)."""
    from scipy.optimize import linear_sum_assignment

    size = reference.shape[1]
    correlation = np.abs(np.corrcoef(reference.T, other.T)[:size, size:])
    rows, columns = linear_sum_assignment(-correlation)
    return float(correlation[rows, columns].mean())
