"""Decision 9.10(a), the owner on 14 October 2026: "Go ahead with fixing the data for v4".

OmniPath records is_stimulation=True for the G-protein-coupled receptor kinases and the arrestins acting on a
receptor. Phosphorylating a receptor's intracellular face uncouples it from its G protein, so the propagation,
which reads a sign as the effect on the destination's activity, needs -1 there.

The second test is the one that matters over time: it fails if the source gains a stimulation edge from a
desensitising enzyme into a receptor that the curated table does not cover, so a refresh of OmniPath cannot
silently reintroduce the defect.
"""
import csv
import re
from pathlib import Path

import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.build_physiology_graph import (
    build_signaling_layer,
    interaction_sign,
    load_curated_edge_sign_overrides,
)

CURATED_OVERRIDES_PATH = Path("configs/curated_edge_sign_overrides.csv")
OMNIPATH_INTERACTIONS_PATH = Path("data/raw/omnipath/omnipath_interactions.tsv")
DESENSITISING_ENZYMES = {f"GRK{index}" for index in range(1, 8)} | {"ARRB1", "ARRB2", "ARR3", "SAG"}
RECEPTOR_FAMILY_PATTERN = re.compile(
    r"^(HTR|DRD|ADRA|ADRB|CHRM|HRH|OPR[MDKL]|GRM|GABBR|ADORA|CNR|TAAR|MTNR|GPR|P2RY|AVPR|OXTR|SSTR|EDNR|F2R|AGTR|PTGER|LPAR|S1PR)")


def test_every_curated_override_turns_a_stated_activation_into_an_inhibition():
    overrides = load_curated_edge_sign_overrides(CURATED_OVERRIDES_PATH)
    assert overrides, "the curated override table is empty"
    for (source_symbol, target_symbol), override in overrides.items():
        corrected_relation, corrected_sign, sourced_relation, sourced_sign = override
        assert (sourced_relation, sourced_sign) == ("activates", 1.0), f"{source_symbol}->{target_symbol}"
        assert (corrected_relation, corrected_sign) == ("inhibits", -1.0), f"{source_symbol}->{target_symbol}"
        assert source_symbol in DESENSITISING_ENZYMES, f"{source_symbol} is not a desensitising enzyme"


def test_the_override_is_applied_and_inverts_the_sign_in_the_built_layer(tmp_path):
    interactions_path = tmp_path / "interactions.tsv"
    with open(interactions_path, "w", newline="", encoding="utf-8") as interactions_file:
        writer = csv.DictWriter(interactions_file, delimiter="\t",
                                fieldnames=["source_genesymbol", "target_genesymbol", "is_directed",
                                            "is_stimulation", "is_inhibition", "consensus_stimulation",
                                            "consensus_inhibition"])
        writer.writeheader()
        writer.writerow({"source_genesymbol": "GRK5", "target_genesymbol": "ADRB2", "is_directed": "True",
                         "is_stimulation": "True", "is_inhibition": "False",
                         "consensus_stimulation": "True", "consensus_inhibition": "False"})
        writer.writerow({"source_genesymbol": "GRK5", "target_genesymbol": "MAPK1", "is_directed": "True",
                         "is_stimulation": "True", "is_inhibition": "False",
                         "consensus_stimulation": "True", "consensus_inhibition": "False"})
    nodes = pd.DataFrame([{"node_id": "GENE:GRK5", "node_type": "gene", "gene_symbol": "GRK5"}])
    overrides = {("GRK5", "ADRB2"): ("inhibits", -1.0, "activates", 1.0)}

    _, edges, counts = build_signaling_layer(nodes, interactions_path, curated_sign_overrides=overrides)

    corrected = edges[(edges.source_id == "GENE:GRK5") & (edges.target_id == "GENE:ADRB2")]
    untouched = edges[(edges.source_id == "GENE:GRK5") & (edges.target_id == "GENE:MAPK1")]
    assert corrected.sign.tolist() == [-1.0]
    assert corrected.relation_type.tolist() == ["inhibits"]
    assert untouched.sign.tolist() == [1.0], "an edge with no override keeps the sourced sign"
    assert counts["edges_corrected"] == 1
    assert counts["overrides_never_matched"] == 0


def test_an_override_whose_source_no_longer_agrees_is_skipped_rather_than_applied(tmp_path):
    interactions_path = tmp_path / "interactions.tsv"
    with open(interactions_path, "w", newline="", encoding="utf-8") as interactions_file:
        writer = csv.DictWriter(interactions_file, delimiter="\t",
                                fieldnames=["source_genesymbol", "target_genesymbol", "is_stimulation",
                                            "is_inhibition", "consensus_stimulation", "consensus_inhibition"])
        writer.writeheader()
        writer.writerow({"source_genesymbol": "GRK5", "target_genesymbol": "ADRB2", "is_stimulation": "False",
                         "is_inhibition": "True", "consensus_stimulation": "False", "consensus_inhibition": "True"})
    nodes = pd.DataFrame([{"node_id": "GENE:GRK5", "node_type": "gene", "gene_symbol": "GRK5"}])
    overrides = {("GRK5", "ADRB2"): ("inhibits", -1.0, "activates", 1.0)}

    _, edges, counts = build_signaling_layer(nodes, interactions_path, curated_sign_overrides=overrides)

    assert counts["edges_corrected"] == 0
    assert counts["overrides_skipped_because_the_source_no_longer_agrees"] == 1
    assert edges[(edges.source_id == "GENE:GRK5")].sign.tolist() == [-1.0], "the source already says inhibition"


@pytest.mark.skipif(not OMNIPATH_INTERACTIONS_PATH.exists(), reason="the pinned OmniPath download is not present")
def test_the_source_has_no_uncovered_stimulation_from_a_desensitising_enzyme_into_a_receptor():
    overrides = load_curated_edge_sign_overrides(CURATED_OVERRIDES_PATH)
    uncovered_pairs = set()
    with open(OMNIPATH_INTERACTIONS_PATH, encoding="utf-8") as interactions_file:
        for row in csv.DictReader(interactions_file, delimiter="\t"):
            relation, sign = interaction_sign(row)
            if (relation, sign) != ("activates", 1.0):
                continue
            for source_symbol in row["source_genesymbol"].split("_"):
                if source_symbol not in DESENSITISING_ENZYMES:
                    continue
                for target_symbol in row["target_genesymbol"].split("_"):
                    if RECEPTOR_FAMILY_PATTERN.match(target_symbol) and (source_symbol, target_symbol) not in overrides:
                        uncovered_pairs.add((source_symbol, target_symbol))
    assert not uncovered_pairs, (
        "the source states a stimulation from a desensitising enzyme into a receptor that the curated table does "
        f"not cover, so the sign defect of decision 9.10(a) would re-enter the graph: {sorted(uncovered_pairs)}")
