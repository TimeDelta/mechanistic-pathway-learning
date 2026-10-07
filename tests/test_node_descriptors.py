"""Tests for the fixed node descriptors: reduced-rank regression, target handling, ICA matching and the per-type blocks."""
import numpy as np
import pandas as pd

from mechanistic_pathway_learning.graph.node_descriptors import assemble_node_descriptor_table, reaction_enzyme_classes
from mechanistic_pathway_learning.graph.protein_descriptors import (
    choose_rank,
    cross_validated_r_squared,
    fit_reduced_rank_regression,
    matched_component_stability,
    subcellular_locations,
    weighted_standardised_targets,
)


def planted_low_rank_problem(seed: int = 0, proteins: int = 600, embedding_dim: int = 30, true_rank: int = 3, targets: int = 40):
    """Targets depend on the embeddings only through 3 directions that carry little embedding variance, while the
    largest-variance directions carry none: principal components would keep the wrong ones."""
    rng = np.random.default_rng(seed)
    scales = np.ones(embedding_dim)
    scales[:5] = 10.0  # high-variance, irrelevant
    embeddings = rng.normal(size=(proteins, embedding_dim)) * scales
    relevant = embeddings[:, 10:10 + true_rank]
    targets_matrix = relevant @ rng.normal(size=(true_rank, targets)) + 0.3 * rng.normal(size=(proteins, targets))
    return embeddings, targets_matrix, relevant


def test_reduced_rank_regression_finds_low_variance_relevant_directions_and_the_rank() -> None:
    embeddings, targets, relevant = planted_low_rank_problem()
    r_squared = cross_validated_r_squared(embeddings, targets, [1, 2, 3, 4, 8], relative_alpha=1e-3)
    assert choose_rank(r_squared) == 3
    model = fit_reduced_rank_regression(embeddings, targets, rank=3, relative_alpha=1e-3)
    descriptors = model.transform(embeddings)
    # the descriptors span the relevant subspace: regressing the planted directions on them explains almost everything
    coefficients, *_ = np.linalg.lstsq(np.c_[descriptors, np.ones(len(descriptors))], relevant, rcond=None)
    residual = relevant - np.c_[descriptors, np.ones(len(descriptors))] @ coefficients
    assert 1 - residual.var(axis=0).sum() / relevant.var(axis=0).sum() > 0.98
    assert np.allclose(descriptors.std(axis=0), 1.0, atol=1e-6)
    # the first three principal components keep the high-variance directions and miss the relevant ones
    centred = embeddings - embeddings.mean(axis=0)
    principal_components = centred @ np.linalg.svd(centred, full_matrices=False)[2][:3].T
    coefficients, *_ = np.linalg.lstsq(np.c_[principal_components, np.ones(len(centred))], relevant, rcond=None)
    residual = relevant - np.c_[principal_components, np.ones(len(centred))] @ coefficients
    assert 1 - residual.var(axis=0).sum() / relevant.var(axis=0).sum() < 0.1


def test_unannotated_proteins_are_set_to_the_mean_not_counted_as_negatives() -> None:
    targets = pd.DataFrame({"location:a": [1.0, 0.0, 0.0, 1.0], "location:b": [0.0, 1.0, 0.0, 0.0], "pfam:x": [1.0, 0.0, 0.0, 0.0]})
    observed = pd.DataFrame({"location": [True, True, False, True], "pfam": [True, True, True, True]})
    values = weighted_standardised_targets(targets, {"location": ["location:a", "location:b"], "pfam": ["pfam:x"]}, observed)
    assert np.all(values[2, :2] == 0.0)  # the unannotated protein sits at the mean of the location block
    location_block = values[[0, 1, 3], :2] * np.sqrt(2)  # undo the block weight
    assert np.allclose(location_block.mean(axis=0), 0.0) and np.allclose(location_block.std(axis=0), 1.0)


def test_subcellular_location_text_is_parsed_into_terms() -> None:
    text = ("SUBCELLULAR LOCATION: [Isoform 1]: Membrane {ECO:0000255}; Single-pass type I membrane protein {ECO:0000255}.; "
            "SUBCELLULAR LOCATION: [Isoform 2]: Secreted {ECO:0000269|PubMed:26221034}. Note=Translocates to the nucleus.")
    assert subcellular_locations(text) == {"membrane", "single-pass type i membrane protein", "secreted"}
    assert subcellular_locations(float("nan")) == set()


def test_matched_component_stability_ignores_sign_and_order() -> None:
    rng = np.random.default_rng(0)
    components = rng.normal(size=(500, 4))
    assert matched_component_stability(components, -components[:, [2, 0, 3, 1]]) > 0.999
    assert matched_component_stability(components, rng.normal(size=(500, 4))) < 0.2


def test_blocks_are_zero_outside_their_node_type_and_ec_classes_are_read() -> None:
    nodes = pd.DataFrame({"node_id": ["MAM1c", "MAM1e", "MAR1", "GENE:A", "GENE:B"], "node_type": ["metabolite", "metabolite", "reaction", "gene", "gene"],
                          "base_metabolite_id": ["MAM1", "MAM1", None, None, None], "gene_symbol": [None, None, None, "A", "B"]})
    metabolite_table = pd.DataFrame({"log_molecular_weight": [5.0], "logp": [0.5], "net_charge": [0.0], "log_polar_surface_area": [4.0],
                                     "log_hydrogen_bond_donors": [1.0], "log_hydrogen_bond_acceptors": [1.0], "log_rotatable_bonds": [1.0],
                                     "log_rings": [0.7], "has_structure": [1.0], "partial_structure": [0.0]}, index=["MAM1"])
    sbml = '<reaction metaid="m" id="R_MAR1" name="x"> identifiers.org/ec-code/1.14.16.1 identifiers.org/ec-code/2.6.1.1 </reaction>'
    reaction_table = reaction_enzyme_classes(sbml)
    assert reaction_table.loc["MAR1", ["ec_class_1", "ec_class_2", "ec_class_3", "has_ec"]].tolist() == [1.0, 1.0, 0.0, 1.0]
    protein_table = pd.DataFrame({"rrr_1": [0.5], "rrr_2": [-1.0]}, index=pd.Index(["A"], name="gene_symbol"))
    table = assemble_node_descriptor_table(nodes, metabolite_table, reaction_table, protein_table)
    assert (table.loc[["MAR1", "GENE:A", "GENE:B"], table.columns.str.startswith("metabolite_")] == 0).all().all()
    assert table.loc["MAM1e", "metabolite_has_structure"] == 1.0  # every compartment copy gets the base metabolite's values
    assert table.loc["GENE:A", ["protein_rrr_1", "protein_rrr_2", "protein_has_protein_descriptors"]].tolist() == [0.5, -1.0, 1.0]
    assert table.loc["GENE:B", "protein_has_protein_descriptors"] == 0.0  # no reviewed entry: zeros and the flag at 0
    assert table.loc["MAR1", "reaction_ec_class_1"] == 1.0 and table.loc["MAM1c", "reaction_ec_class_1"] == 0.0
