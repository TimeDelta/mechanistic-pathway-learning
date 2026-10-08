"""The degree used for strata and normalisation hops over encodes on a graph that splits genes from proteins."""
from pathlib import Path

import numpy as np
import pytest

from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data
from mechanistic_pathway_learning.evaluation.ranking_and_calibration_metrics import degree_strata

MERGED_GRAPH = Path("data/processed/graph")
SPLIT_GRAPH = Path("data/processed/graph_split")
EVIDENCE = Path("data/processed/evidence")


def load(graph_directory: Path):
    if not (graph_directory / "nodes.parquet").exists() or not (EVIDENCE / "evidence_records.parquet").exists():
        pytest.skip(f"{graph_directory} or {EVIDENCE} is not built here")
    return load_experiment_data(graph_directory, EVIDENCE, group_by="disease_cluster")


def test_the_seed_degree_collapses_on_the_split_graph_and_the_encodes_hop_does_not():
    split = load(SPLIT_GRAPH)
    assert len(np.unique(split.perturbation_degrees)) == 1, "the split graph is expected to give every perturbed gene degree 1"
    assert len(np.unique(split.perturbation_degrees_with_encoded_proteins)) > 1
    assert len(np.unique(degree_strata(split.perturbation_degrees_with_encoded_proteins))) > 1


def test_the_split_strata_are_the_merged_graph_strata():
    """The registered stratification is the merged graph's; the hop over encodes reproduces it."""
    merged, split = load(MERGED_GRAPH), load(SPLIT_GRAPH)
    assert list(merged.perturbation_ids) == list(split.perturbation_ids)
    assert np.array_equal(degree_strata(merged.perturbation_degrees), degree_strata(split.perturbation_degrees_for_strata))


def test_the_encodes_hop_recovers_the_merged_degree_itself_on_the_slice():
    """Not only the strata: on the slice the hop over encodes gives back the merged graph's degree value."""
    merged, split = load(MERGED_GRAPH), load(SPLIT_GRAPH)
    assert np.array_equal(split.perturbation_degrees_for_strata, merged.perturbation_degrees)


def test_a_merged_graph_keeps_the_seed_degree():
    merged = load(MERGED_GRAPH)
    assert np.array_equal(merged.perturbation_degrees_for_strata, merged.perturbation_degrees)
    assert np.array_equal(merged.perturbation_degrees_with_encoded_proteins, merged.perturbation_degrees)  # no encodes relation to hop over


def test_the_encodes_hop_is_at_least_the_seed_degree_everywhere():
    for graph_directory in (MERGED_GRAPH, SPLIT_GRAPH):
        data = load(graph_directory)
        assert (data.perturbation_degrees_with_encoded_proteins >= data.perturbation_degrees).all()
