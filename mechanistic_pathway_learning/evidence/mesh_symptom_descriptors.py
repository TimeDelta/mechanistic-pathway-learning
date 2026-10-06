"""MeSH descriptors of the target symptoms, read from docs/symptom_crosswalk.csv (design sections 3.1 and 4.2).

The crosswalk carries two aligned semicolon lists per target symptom: mesh_descriptors (D-numbers) and
mesh_descriptor_level, each entry "symptom" or "diagnosis". Symptom-level descriptors (Anxiety D001007,
Depression D003863) name the sign or symptom itself; diagnosis-level descriptors (Anxiety Disorders D001008,
Depressive Disorder D003866) name a diagnostic category whose literature is about the disorder rather than the
symptom. Both are retrieved, because the literature on a diagnosis is a soft prior on its cardinal symptoms,
and the level travels with every report so the reliability model and the leakage audit can weight the two
apart (assumption A7 keeps diagnoses out of the targets, not out of the priors). A third level,
"mixed_polarity_diagnosis", marks a diagnosis whose cardinal features include the opposite of the target symptom
(Bipolar Disorder D001714 for elevated mood or mania): its literature is kept, but a cause or treat relation on it
is demoted to the undirected associated_with relation by the loaders, since "treats bipolar disorder" carries no
sign for mania.

Literature reports are grade C to E soft priors (design section 4.2): they never become evaluation positives,
and unobserved pairs stay unlabelled.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

DESCRIPTOR_LEVELS: tuple[str, ...] = ("symptom", "diagnosis", "mixed_polarity_diagnosis")  # the third: a diagnosis whose cardinal features include the opposite symptom (Bipolar Disorder for mania), whose directional relations are demoted
MESH_DESCRIPTOR_PATTERN = re.compile(r"^D\d{6,9}$")
LIST_SEPARATOR = ";"


@dataclass(frozen=True)
class SymptomMeshDescriptor:
    target_symptom: str
    mesh_descriptor: str  # D-number without the "MESH:" prefix
    level: str  # "symptom" or "diagnosis"


def split_list_cell(cell: str | None) -> list[str]:
    """Entries of a semicolon-separated crosswalk cell, blanks removed."""
    if cell is None:
        return []
    return [entry.strip() for entry in str(cell).split(LIST_SEPARATOR) if entry.strip()]


def parse_descriptor_columns(target_symptom: str, descriptors_cell: str | None, levels_cell: str | None) -> list[SymptomMeshDescriptor]:
    """The aligned (descriptor, level) pairs of one crosswalk row; raises when the lists disagree in length, an id
    is not a MeSH D-number or a level is not in DESCRIPTOR_LEVELS, so a malformed crosswalk fails loudly."""
    descriptors = split_list_cell(descriptors_cell)
    levels = split_list_cell(levels_cell)
    if len(descriptors) != len(levels):
        raise ValueError(f"{target_symptom}: {len(descriptors)} mesh descriptors but {len(levels)} levels")
    parsed: list[SymptomMeshDescriptor] = []
    for descriptor, level in zip(descriptors, levels):
        if not MESH_DESCRIPTOR_PATTERN.match(descriptor):
            raise ValueError(f"{target_symptom}: {descriptor!r} is not a MeSH descriptor id")
        if level not in DESCRIPTOR_LEVELS:
            raise ValueError(f"{target_symptom}: level {level!r} for {descriptor} is not one of {DESCRIPTOR_LEVELS}")
        parsed.append(SymptomMeshDescriptor(target_symptom, descriptor, level))
    return parsed


def load_symptom_mesh_descriptors(crosswalk_path: Path) -> list[SymptomMeshDescriptor]:
    """Every (symptom, descriptor, level) of the crosswalk, in file order. A descriptor may belong to one symptom only."""
    with open(crosswalk_path, newline="", encoding="utf-8") as crosswalk_file:
        rows = list(csv.DictReader(crosswalk_file))
    if rows and "mesh_descriptors" not in rows[0]:
        raise ValueError(f"{crosswalk_path} has no mesh_descriptors column")
    descriptors: list[SymptomMeshDescriptor] = []
    for row in rows:
        descriptors.extend(parse_descriptor_columns(row["target_symptom"], row.get("mesh_descriptors"), row.get("mesh_descriptor_level")))
    seen: dict[str, str] = {}
    for descriptor in descriptors:
        if descriptor.mesh_descriptor in seen and seen[descriptor.mesh_descriptor] != descriptor.target_symptom:
            raise ValueError(f"MeSH {descriptor.mesh_descriptor} is listed under both {seen[descriptor.mesh_descriptor]} and {descriptor.target_symptom}")
        seen[descriptor.mesh_descriptor] = descriptor.target_symptom
    return descriptors


def descriptor_lookup(descriptors: list[SymptomMeshDescriptor]) -> dict[str, SymptomMeshDescriptor]:
    """D-number -> its crosswalk entry, the filter used on the relation streams."""
    return {descriptor.mesh_descriptor: descriptor for descriptor in descriptors}
