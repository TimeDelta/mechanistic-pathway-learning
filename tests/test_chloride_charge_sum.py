"""The chloride charge-sum defect of 14 October 2026, held so it cannot drift while it waits for the rebuild.

Found by the physiology review of 13-14 October. The graph carries chloride under two identifiers: `MAM01442*`,
which is Human-GEM's chloride (CHEBI:17996, SBML charge -1), and `MAM02144*`, three of whose four compartments are
display-named "Cl-" and which is in **neither** the Human-GEM SBML nor its metabolites table. `INORGANIC_IONS` keys
chloride on `MAM01442` alone, `ion_charges` is built by filtering SBML charges to those bases, and
`reactome_import.py` computes `inward_charge = sum(ion_charges.get(base, 0) * amount ...)`. So a Reactome reaction
whose chloride is `MAM02144` contributes nothing to the charge sum and is not even counted as moving an ion.

These tests assert the defect as it stands rather than the behaviour that is wanted, because the fix requires the
graph to be rebuilt and the rebuild is decision 1.2, the owner's. Each one names what it should say after the fix,
so whoever rebuilds the graph can invert them and watch them pass.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

GRAPH = Path("data/processed/graph_full_neuronal")
ELECTRICAL = ("changes_membrane_potential", "changes_ion_pool")

pytestmark = pytest.mark.skipif(not (GRAPH / "edges.parquet").exists(),
                                reason="the confirmatory graph is not built in this checkout")


@pytest.fixture(scope="module")
def graph() -> tuple[pd.DataFrame, pd.DataFrame]:
    return (pd.read_parquet(GRAPH / "nodes.parquet"),
            pd.read_parquet(GRAPH / "edges.parquet", columns=["source_id", "target_id", "relation_type", "sign"]))


def test_two_chloride_identifiers_exist_and_only_one_is_in_the_charge_table(graph):
    """The root of it: `MAM02144` is named "Cl-" in three compartments and is in no Human-GEM table."""
    nodes, _ = graph
    human_gem_chloride = sorted(nodes.node_id[nodes.node_id.str.startswith("MAM01442")])
    reactome_chloride = sorted(nodes.node_id[nodes.node_id.str.startswith("MAM02144")])
    assert human_gem_chloride, "Human-GEM chloride must be present"
    assert reactome_chloride, "if this is empty the rebuild has merged the two pools and this file can be deleted"
    named_chloride = nodes[nodes.node_id.isin(reactome_chloride) & (nodes.display_name == "Cl-")]
    assert len(named_chloride) >= 3, "three of MAM02144's compartments are display-named Cl-"


def test_the_ligand_gated_chloride_channels_have_no_electrical_edge(graph):
    """GABA_A and the glycine receptor move only chloride, so their whole charge is dropped.

    AFTER THE FIX both should carry `changes_membrane_potential` at -1 (an inward anion current hyperpolarises)."""
    _, edges = graph
    for reaction, what in (("RCT_975340", "GABA_A heteropentamer"), ("RCT_975389", "glycine receptor")):
        electrical = edges[(edges.source_id == reaction) & edges.relation_type.isin(ELECTRICAL)]
        assert electrical.empty, f"{what} ({reaction}) now has an electrical edge; invert this test"


def test_two_electroneutral_cotransporters_are_read_as_opposing_currents(graph):
    """NKCC1 moves Na+ + K+ + 2Cl- inward and KCC2 moves K+ + Cl- outward. Both are electroneutral.

    With chloride dropped, NKCC1's Na+ and K+ leave +1 and KCC2's K+ leaves -1, which is what the graph carries.
    AFTER THE FIX neither should carry a `changes_membrane_potential` edge at all."""
    _, edges = graph
    signs = {}
    for reaction in ("RCT_426086", "RCT_426155"):
        row = edges[(edges.source_id == reaction) & (edges.relation_type == "changes_membrane_potential")]
        assert len(row) == 1, f"{reaction} should carry exactly one membrane-potential edge while the defect stands"
        signs[reaction] = float(row.sign.iloc[0])
    assert signs["RCT_426086"] == pytest.approx(1.0), "NKCC1 reads as depolarising, from Na+ and K+ alone"
    assert signs["RCT_426155"] == pytest.approx(-1.0), "KCC2 reads as hyperpolarising, from K+ alone"
    assert signs["RCT_426086"] == -signs["RCT_426155"], "the two electroneutral transporters oppose each other"


def test_no_reaction_moving_reactome_chloride_feeds_the_neuronal_chloride_pool(graph):
    """The pool node exists and the reactions that should fill it cannot reach it."""
    nodes, edges = graph
    assert "NEURON_MAM01442c" in set(nodes.node_id), "the neuronal chloride pool is built"
    reactome_chloride = set(nodes.node_id[nodes.node_id.str.startswith("MAM02144")])
    touching = edges[edges.source_id.isin(reactome_chloride) | edges.target_id.isin(reactome_chloride)]
    reactions = set(touching.source_id) | set(touching.target_id)
    into_the_pool = edges[(edges.target_id == "NEURON_MAM01442c") & edges.source_id.isin(reactions)]
    assert into_the_pool.empty, "a MAM02144 reaction now reaches the chloride pool; invert this test"
