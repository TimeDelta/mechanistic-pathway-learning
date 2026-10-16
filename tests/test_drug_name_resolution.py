"""Every drug label resolves to a curated name, and none is a fragment (the owner's ruling, 15 October 2026).

Why this exists. Four name-key failures landed in one session - `catalyzed_by` matched against `PROTEIN:` nodes, an
OnSIDES ATC lookup matched names against id columns, an exposure glob returned nothing, and `flupentixol` missed
`flupenthixol` on a British `h` - and each looked like a clean result until the key was checked. The owner's response
was to test the class rather than the instances: "Sounds like a good start".

What is checked. `data/processed/label_review/drug_names.csv` carries, per drug, the `perturbation_label` that
propagates into the evidence records and from there into the exposure table, beside `rxnorm_names`, a curated generic
list. A label that matches no curated name is either a real synonym (the project's name and RxNorm's differ) or a
defect (a truncated fragment, or a brand used where a generic belongs). The synonyms are enumerated below, so a new
mismatch fails loudly instead of joining a crowd of known ones.

What is deliberately NOT used. `pubchem_title` is not an authoritative generic name - it is often an IUPAC systematic
name ("13-Methyl-6,7,8,9,11,12,14,15,16,1...") or a brand ("Evoxac", "Celestone", "Didrate"). Testing the label
against it flags 77 of 295 rows, almost all of them correct labels. That detector was written first and discarded.

The defect pins follow the convention of tests/test_chloride_charge_sum.py: they assert the broken state so it cannot
be forgotten, and they are inverted when the rebuild fixes it.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DRUG_NAMES_PATH = REPOSITORY_ROOT / "data" / "processed" / "label_review" / "drug_names.csv"

# Real synonym pairs: the project's label and RxNorm's generic are both correct names for one drug. Each is an INN
# against a USAN, an older name against the current one, or an endogenous substance against its medicinal name.
ACCEPTED_SYNONYMS: dict[str, str] = {
    "CID100003640": "cortisol is hydrocortisone",
    "CID100004536": "norethisterone is the INN, norethindrone the USAN",
    "CID100005195": "deprenyl is the older name for selegiline",
    "CID100010413": "hydroxybutyrate is the acid, oxybate the salt dispensed",
    "CID100000861": "triiodothyronine is the hormone, liothyronine the medicine",
    "CID100121892": "retigabine is the INN, ezogabine the USAN",
}

# Labels that are wrong rather than synonymous, pinned until the rebuild fixes the name resolution. Since 16 October 2026
# the table is built by experiments/label_review/build_drug_names.py from the registered v4 drug set (222 drugs), so a
# pinned drug outside that set (dronabinol and etonogestrel entered through the approved-variant tables of 10 October)
# is checked only when the table holds it; the tests say which pins were not checked.
KNOWN_BAD_LABELS: dict[str, tuple[str, str]] = {
    "CID109966051": ("Lu", "vortioxetine"),      # truncated from the development code Lu AA21004
    "CID100002978": ("delta", "dronabinol"),     # truncated from a delta-9- prefixed name
    "CID100040976": ("Implanon", "etonogestrel"),  # a brand name standing where a generic belongs
}

SHORTEST_PLAUSIBLE_DRUG_LABEL = 5


def normalise(value) -> str:
    return str(value).strip().lower()


@pytest.fixture(scope="module")
def drug_names() -> pd.DataFrame:
    if not DRUG_NAMES_PATH.exists():
        pytest.skip(f"{DRUG_NAMES_PATH} is not built in this checkout")
    return pd.read_csv(DRUG_NAMES_PATH)


def labels_matching_no_curated_name(drug_names: pd.DataFrame) -> dict[str, tuple[str, str]]:
    """perturbation_id -> (label, curated names) for every label that matches no name in rxnorm_names."""
    unmatched: dict[str, tuple[str, str]] = {}
    for row in drug_names.itertuples(index=False):
        curated = normalise(row.rxnorm_names)
        if curated in ("nan", ""):
            continue
        label = normalise(row.perturbation_label)
        curated_names = [name.strip() for name in curated.split(";") if name.strip()]
        if not any(label == name or label in name or name in label for name in curated_names):
            unmatched[row.perturbation_id] = (str(row.perturbation_label), str(row.rxnorm_names))
    return unmatched


def test_no_new_label_fails_to_resolve(drug_names):
    """A label matching no curated name is a synonym or a defect; both are enumerated, so anything else is new."""
    unmatched = labels_matching_no_curated_name(drug_names)
    accounted = set(ACCEPTED_SYNONYMS) | set(KNOWN_BAD_LABELS)
    surprises = {identifier: names for identifier, names in unmatched.items() if identifier not in accounted}
    assert not surprises, (
        "a drug label resolves to no curated RxNorm name and is neither a listed synonym nor a known defect; "
        f"check whether it is a fragment or a spelling variant before adding it: {surprises}")


def test_the_known_bad_labels_are_still_the_only_bad_ones(drug_names):
    """Pins the defect. When the rebuild fixes the name resolution this fails, and the entries move out of
    KNOWN_BAD_LABELS rather than the assertion being loosened."""
    unmatched = labels_matching_no_curated_name(drug_names)
    present = {identifier for identifier in KNOWN_BAD_LABELS if identifier in set(drug_names.perturbation_id)}
    still_broken = {identifier for identifier in present if identifier in unmatched}
    assert still_broken == present, (
        "a pinned bad label now resolves, which is the fix landing; remove it from KNOWN_BAD_LABELS. "
        f"Still broken: {sorted(still_broken)}; pins outside this table's drug set, not checked: {sorted(set(KNOWN_BAD_LABELS) - present)}")


def test_no_label_is_a_bare_fragment(drug_names):
    """`Lu` is four characters and `delta` five, so length alone does not separate them - but a label shorter than a
    plausible drug name is always worth stopping on. Only the pinned defects may be this short."""
    short = {row.perturbation_id: str(row.perturbation_label)
             for row in drug_names.itertuples(index=False)
             if len(str(row.perturbation_label).strip()) < SHORTEST_PLAUSIBLE_DRUG_LABEL}
    unexpected = {identifier: label for identifier, label in short.items() if identifier not in KNOWN_BAD_LABELS}
    assert not unexpected, f"drug labels too short to be a drug name: {unexpected}"


def test_the_curated_name_is_available_for_every_known_defect(drug_names):
    """The fix needs no new data: the correct name is already in the same row. This asserts that, so the fix cannot be
    argued to need a fetch."""
    indexed = drug_names.set_index("perturbation_id")
    checked = 0
    for identifier, (_, expected_generic) in KNOWN_BAD_LABELS.items():
        if identifier not in indexed.index:  # a pin for a drug outside this table's drug set
            continue
        checked += 1
        row = indexed.loc[identifier]
        available = normalise(row.rxnorm_names) + " " + normalise(row.pubchem_title)
        assert expected_generic in available, (
            f"{identifier}: the correct generic {expected_generic!r} is not in the row's own rxnorm_names or "
            f"pubchem_title, so the fix would need a lookup after all")
    assert checked >= 1, "no pinned drug is in the table, so the defect class is unchecked"
