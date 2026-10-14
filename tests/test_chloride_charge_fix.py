"""Decision 9.4, the chloride charge fix of 14 October 2026 (the owner: "Definitely the second fix for chloride so
it's separate from hcl").

`MAM02144` is Reactome's chloride and no Human-GEM table holds it, so it had no charge and `ions_moved_inward` did
not count it. It now takes chloride's charge and is counted as an ion, **without** being merged into `MAM01442` and
without becoming a second chloride pool: the two base ids stay separate entities and only the charge arithmetic is
repaired.

These hold the arithmetic the fix produces. The graph itself changes only when it is rebuilt (decision 1.2), and
`tests/test_chloride_charge_sum.py` holds the defect until then.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "experiments" / "build_neuronal_graph_variant.py"
specification = importlib.util.spec_from_file_location("build_neuronal_graph_variant", SCRIPT)
build = importlib.util.module_from_spec(specification)
specification.loader.exec_module(build)

from mechanistic_pathway_learning.graph.reactome_import import ions_moved_inward  # noqa: E402

CHLORIDE = build.IONS_WITHOUT_A_HUMAN_GEM_CHARGE["Cl- (Reactome)"][0]
CHARGES = {"MAM02519": 1, "MAM02200": 1, CHLORIDE: -1}
ION_BASES = {"Na+": "MAM02519", "K+": "MAM02200", "Cl- (Reactome)": CHLORIDE}


def net_charge(participants) -> float:
    moved = ions_moved_inward(participants, set(ION_BASES.values()))
    return sum(CHARGES.get(base, 0) * amount for base, amount in moved.items())


def test_reactome_chloride_now_carries_chloride_s_charge():
    assert build.IONS_WITHOUT_A_HUMAN_GEM_CHARGE["Cl- (Reactome)"] == (CHLORIDE, -1)


def test_the_gaba_a_receptor_becomes_hyperpolarising():
    """One Cl- from extracellular to cytosol: net -1 inward, so the edge is -1. `RCT_975340` carried none at all."""
    gaba_a = [("in", CHLORIDE, "e", 1.0), ("out", CHLORIDE, "c", 1.0)]
    assert net_charge(gaba_a) == pytest.approx(-1.0)


def test_nkcc1_becomes_electroneutral_because_reactome_states_two_chlorides():
    """Na+ + K+ + 2Cl- inward is zero charge, and Reactome's SBML gives stoichiometry 2 for the chloride, so the
    sum is exactly zero and the reaction stops carrying a current. It read +1 before, from Na+ and K+ alone."""
    nkcc1 = [("in", "MAM02519", "e", 1.0), ("in", "MAM02200", "e", 1.0), ("in", CHLORIDE, "e", 2.0),
             ("out", "MAM02519", "c", 1.0), ("out", "MAM02200", "c", 1.0), ("out", CHLORIDE, "c", 2.0)]
    assert net_charge(nkcc1) == pytest.approx(0.0)


def test_kcc2_becomes_electroneutral():
    """K+ + Cl- outward is zero charge. It read -1 before, from K+ alone."""
    kcc2 = [("in", "MAM02200", "c", 1.0), ("in", CHLORIDE, "c", 1.0),
            ("out", "MAM02200", "e", 1.0), ("out", CHLORIDE, "e", 1.0)]
    assert net_charge(kcc2) == pytest.approx(0.0)


def test_the_two_chloride_pools_are_not_merged_and_no_second_pool_is_built():
    """The owner's choice of the two fixes: keep MAM02144 separate from MAM01442 and from hydrochloride."""
    assert CHLORIDE not in set(build.INORGANIC_IONS.values()), "the Human-GEM chloride base is untouched"
    assert build.INORGANIC_IONS["Cl-"] == "MAM01442"
    from mechanistic_pathway_learning.graph.reactome_import import NEURONAL_POOL_DRIVING_FORCE
    assert not any(name == "Cl- (Reactome)" for name, _ in NEURONAL_POOL_DRIVING_FORCE), \
        "no driving force is declared for it, so no pool node is built"
