"""Decision 9.17, the owner on 14 October 2026: "good call on prophylactically changing the mapper to ChEBl term
Reactome names".

It is not a precaution. Neither Human-GEM entity that Reactome's chloride can land on is chloride by its own
annotation: MAM02144 is hydrogen chloride (CHEBI:17883, charge 0) and MAM01442 is the chlorine radical
(CHEBI:29311, charge 0). Reactome names CHEBI:17996, chloride, charge -1, and that is the charge the current needs.
"""
from pathlib import Path

import pytest

from experiments.build_neuronal_graph_variant import ION_CHEBI_TERM_OF_BASE
from mechanistic_pathway_learning.graph.reactome_import import read_chebi_charges

CHEBI_ONTOLOGY_PATH = Path("data/raw/chebi/chebi_core.obo.gz")
EXPECTED_CHARGE_OF_TRACKED_ION = {"MAM02519": 1, "MAM02200": 1, "MAM01413": 2,
                                  "MAM01442": -1, "MAM02482": 2, "MAM02144": -1}


@pytest.mark.skipif(not CHEBI_ONTOLOGY_PATH.exists(), reason="the pinned ChEBI ontology is not present")
def test_every_tracked_ion_takes_its_charge_from_the_chebi_term_reactome_names():
    chebi_charges = read_chebi_charges(CHEBI_ONTOLOGY_PATH)
    for base_metabolite_id, chebi_term in ION_CHEBI_TERM_OF_BASE.items():
        assert chebi_term in chebi_charges, f"{chebi_term} carries no charge in the pinned ontology"
        assert chebi_charges[chebi_term] == EXPECTED_CHARGE_OF_TRACKED_ION[base_metabolite_id], base_metabolite_id


@pytest.mark.skipif(not CHEBI_ONTOLOGY_PATH.exists(), reason="the pinned ChEBI ontology is not present")
def test_the_two_candidate_chloride_entities_are_both_neutral_by_their_own_annotation():
    """The reason the guard is needed rather than a choice between the two Human-GEM entities."""
    chebi_charges = read_chebi_charges(CHEBI_ONTOLOGY_PATH)
    assert chebi_charges["CHEBI:17883"] == 0, "hydrogen chloride, which MAM02144 is annotated as"
    assert chebi_charges["CHEBI:29311"] == 0, "the chlorine radical, which MAM01442 is annotated as"
    assert chebi_charges["CHEBI:17996"] == -1, "chloride, which Reactome names"


@pytest.mark.skipif(not CHEBI_ONTOLOGY_PATH.exists(), reason="the pinned ChEBI ontology is not present")
def test_the_guard_gives_the_nkcc1_and_gaba_a_arithmetic_the_owner_verified():
    """GABA_A moves one Cl- inward, so -1. NKCC1 is Na+ + K+ + 2Cl-, which sums to 0 and drops the current."""
    chebi_charges = read_chebi_charges(CHEBI_ONTOLOGY_PATH)
    charge_of = {base: chebi_charges[term] for base, term in ION_CHEBI_TERM_OF_BASE.items()}
    gaba_a_inward_charge = charge_of["MAM02144"] * 1
    nkcc1_inward_charge = charge_of["MAM02519"] * 1 + charge_of["MAM02200"] * 1 + charge_of["MAM01442"] * 2
    kcc2_outward_charge = charge_of["MAM02200"] * 1 + charge_of["MAM01442"] * 1
    assert gaba_a_inward_charge == -1
    assert nkcc1_inward_charge == 0
    assert kcc2_outward_charge == 0
