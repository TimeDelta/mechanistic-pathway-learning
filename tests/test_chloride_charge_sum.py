"""The chloride charge-sum defect of 14 October 2026, held so it cannot drift while it waits for the rebuild.

Found by the physiology review of 13-14 October. The graph carries chloride under two identifiers: `MAM01442*`,
which is Human-GEM's chloride (CHEBI:17996, SBML charge -1), and `MAM02144*`, three of whose four compartments are
display-named "Cl-" and which is in **neither** the Human-GEM SBML nor its metabolites table. `INORGANIC_IONS` keys
chloride on `MAM01442` alone, `ion_charges` is built by filtering SBML charges to those bases, and
`reactome_import.py` computes `inward_charge = sum(ion_charges.get(base, 0) * amount ...)`. So a Reactome reaction
whose chloride is `MAM02144` contributes nothing to the charge sum and is not even counted as moving an ion.

These tests asserted the defect as it stood rather than the behaviour that is wanted, because the fix requires the
graph to be rebuilt and the rebuild is decision 1.2, the owner's. Each one named what it should say after the fix,
so whoever rebuilds the graph could invert them and watch them pass.

The rebuild of 15 October 2026 (a new container, scripts/rebuild_data_layer.sh, with the charge guard of commit
090f001 in the code) gave the post-fix graph: GABA_A and the glycine receptor each carry `changes_membrane_potential`
at -1, and NKCC1 and KCC2 carry no membrane-potential edge. The two tests below that pinned those defects are
inverted to that, as their docstrings asked. The two chloride identifiers still exist and no MAM02144 reaction
reaches the neuronal chloride pool, so the first and last tests still hold the defect they pin.
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


def test_the_ligand_gated_chloride_channels_carry_an_inward_anion_current(graph):
    """GABA_A and the glycine receptor move only chloride. While chloride's charge was dropped they had no electrical
    edge at all; with it read from the ChEBI term Reactome names, each carries `changes_membrane_potential` at -1
    (an inward anion current hyperpolarises). Inverted after the rebuild of 15 October 2026."""
    _, edges = graph
    for reaction, what in (("RCT_975340", "GABA_A heteropentamer"), ("RCT_975389", "glycine receptor")):
        electrical = edges[(edges.source_id == reaction) & edges.relation_type.isin(ELECTRICAL)]
        assert list(electrical.relation_type) == ["changes_membrane_potential"], f"{what} ({reaction}) should carry one membrane-potential edge"
        assert float(electrical.sign.iloc[0]) == pytest.approx(-1.0), f"{what} ({reaction}): an inward chloride current hyperpolarises"


def test_two_electroneutral_cotransporters_carry_no_current(graph):
    """NKCC1 moves Na+ + K+ + 2Cl- inward and KCC2 moves K+ + Cl- outward. Both are electroneutral.

    While chloride's charge was dropped, NKCC1's Na+ and K+ left +1 and KCC2's K+ left -1, and the graph carried the
    two as opposing currents. With chloride counted, neither carries a `changes_membrane_potential` edge at all.
    Inverted after the rebuild of 15 October 2026."""
    _, edges = graph
    for reaction, what in (("RCT_426086", "NKCC1"), ("RCT_426155", "KCC2")):
        row = edges[(edges.source_id == reaction) & (edges.relation_type == "changes_membrane_potential")]
        assert row.empty, f"{what} ({reaction}) is electroneutral and should carry no membrane-potential edge"


def test_no_reaction_moving_reactome_chloride_feeds_the_neuronal_chloride_pool(graph):
    """The pool node exists and the reactions that should fill it cannot reach it."""
    nodes, edges = graph
    assert "NEURON_MAM01442c" in set(nodes.node_id), "the neuronal chloride pool is built"
    reactome_chloride = set(nodes.node_id[nodes.node_id.str.startswith("MAM02144")])
    touching = edges[edges.source_id.isin(reactome_chloride) | edges.target_id.isin(reactome_chloride)]
    reactions = set(touching.source_id) | set(touching.target_id)
    into_the_pool = edges[(edges.target_id == "NEURON_MAM01442c") & edges.source_id.isin(reactions)]
    assert into_the_pool.empty, "a MAM02144 reaction now reaches the chloride pool; invert this test"
