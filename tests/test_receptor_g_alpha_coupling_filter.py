"""Decision 9.23, the owner on 14 October 2026: "yes on the primary and secondary with dropping the 352 unlisted
since it exactly reproduces the true EVIDENCED physiology".

The three cases that matter are kept apart deliberately: an annotated receptor's unlisted coupling is dropped, an
unannotated receptor keeps everything because absence from a curated sheet is not evidence against, and a source
that is not a receptor is never touched, since this rule says nothing about what else regulates a G alpha
(decision 9.27, where the nine RGS and AGS edges were checked and found correctly signed).
"""
import csv
from pathlib import Path

import pandas as pd
import pytest

from mechanistic_pathway_learning.graph.build_physiology_graph import (
    build_signaling_layer,
    load_admitted_receptor_g_alpha_coupling,
    receptor_g_alpha_edge_is_dropped,
)

ADMITTED_COUPLING_PATH = Path("configs/admitted_receptor_g_alpha_coupling.csv")


def test_the_admitted_table_holds_primary_and_secondary_and_expands_families_to_genes():
    admitted_pairs, annotated_receptors = load_admitted_receptor_g_alpha_coupling(ADMITTED_COUPLING_PATH)
    assert admitted_pairs, "the admitted coupling table is empty"
    # HTR1A is primary Gi/o with no secondary, so every Gi/o member is admitted and no Gs or Gq member is
    assert ("HTR1A", "GNAI1") in admitted_pairs
    assert ("HTR1A", "GNAO1") in admitted_pairs
    assert ("HTR1A", "GNAS") not in admitted_pairs
    assert ("HTR1A", "GNAQ") not in admitted_pairs
    # ADRB2 is primary Gs with secondary Gi/o, so both families are admitted
    assert ("ADRB2", "GNAS") in admitted_pairs
    assert ("ADRB2", "GNAI1") in admitted_pairs
    assert "HTR1A" in annotated_receptors and "ADRB2" in annotated_receptors


def test_only_an_annotated_receptors_unlisted_coupling_is_dropped():
    admitted_pairs, annotated_receptors = load_admitted_receptor_g_alpha_coupling(ADMITTED_COUPLING_PATH)
    # an annotated receptor reaching a G alpha the annotation does not list
    assert receptor_g_alpha_edge_is_dropped("HTR1A", "GNAS", 1.0, admitted_pairs, annotated_receptors)
    # the same receptor reaching its own primary family
    assert not receptor_g_alpha_edge_is_dropped("HTR1A", "GNAI1", 1.0, admitted_pairs, annotated_receptors)
    # a secondary coupling survives
    assert not receptor_g_alpha_edge_is_dropped("ADRB2", "GNAI1", 1.0, admitted_pairs, annotated_receptors)
    # a source that is not a receptor is never filtered, however implausible it looks
    assert not receptor_g_alpha_edge_is_dropped("RGS4", "GNAS", -1.0, admitted_pairs, annotated_receptors)
    assert not receptor_g_alpha_edge_is_dropped("PRKACA", "GNAS", 1.0, admitted_pairs, annotated_receptors)
    assert not receptor_g_alpha_edge_is_dropped("MYC", "GNAQ", 1.0, admitted_pairs, annotated_receptors)
    # an unsigned edge is left alone, and a target that is not a G alpha is left alone
    assert not receptor_g_alpha_edge_is_dropped("HTR1A", "GNAS", 0.0, admitted_pairs, annotated_receptors)
    assert not receptor_g_alpha_edge_is_dropped("HTR1A", "MAPK1", 1.0, admitted_pairs, annotated_receptors)


def test_the_filter_drops_the_unlisted_edge_in_the_built_layer(tmp_path):
    interactions_path = tmp_path / "interactions.tsv"
    with open(interactions_path, "w", newline="", encoding="utf-8") as interactions_file:
        writer = csv.DictWriter(interactions_file, delimiter="\t",
                                fieldnames=["source_genesymbol", "target_genesymbol", "is_stimulation",
                                            "is_inhibition", "consensus_stimulation", "consensus_inhibition"])
        writer.writeheader()
        for source_symbol, target_symbol in (("HTR1A", "GNAI1"), ("HTR1A", "GNAS"), ("RGS4", "GNAS")):
            writer.writerow({"source_genesymbol": source_symbol, "target_genesymbol": target_symbol,
                             "is_stimulation": "True", "is_inhibition": "False",
                             "consensus_stimulation": "True", "consensus_inhibition": "False"})
    nodes = pd.DataFrame([{"node_id": "GENE:HTR1A", "node_type": "gene", "gene_symbol": "HTR1A"}])
    coupling = load_admitted_receptor_g_alpha_coupling(ADMITTED_COUPLING_PATH)

    _, edges, counts = build_signaling_layer(nodes, interactions_path, admitted_receptor_g_alpha_coupling=coupling)

    kept = {(row.source_id, row.target_id) for row in edges.itertuples()}
    assert ("GENE:HTR1A", "GENE:GNAI1") in kept, "the primary coupling must survive"
    assert ("GENE:HTR1A", "GENE:GNAS") not in kept, "the unlisted coupling must be dropped"
    assert ("GENE:RGS4", "GENE:GNAS") in kept, "a non-receptor source must not be filtered"
    assert counts["receptor_g_alpha_edges_dropped_as_unlisted"] == 1


def test_an_absent_table_keeps_every_coupling_edge(tmp_path):
    admitted_pairs, annotated_receptors = load_admitted_receptor_g_alpha_coupling(tmp_path / "missing.csv")
    assert (admitted_pairs, annotated_receptors) == (set(), set())
    assert not receptor_g_alpha_edge_is_dropped("HTR1A", "GNAS", 1.0, admitted_pairs, annotated_receptors)
