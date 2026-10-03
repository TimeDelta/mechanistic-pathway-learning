"""GTEx brain expression on gene nodes, reaction inheritance and the regulon restriction (design 3.3 and 4.1)."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.graph.build_physiology_graph import (
    annotate_brain_expression,
    load_brain_expression,
    propagate_brain_expression_to_reactions,
    restrict_transcription_edges_to_brain_expressed,
)


def write_gct(path: Path) -> None:
    header = ["Name", "Description", "Liver", "Brain_Cortex", "Brain_Hypothalamus"]
    rows = [
        ["ENSG00000000419.12", "DPM1", "50", "0.2", "3.5"],   # expressed in one brain tissue
        ["ENSG00000001084.13", "GCLC", "80", "0.1", "0.4"],   # known, below threshold
        ["ENSG00000999999.1", "TFHIGH", "1", "12", "20"],    # matched by symbol only
        ["ENSG00000888888.1", "TFLOW", "9", "0.0", "0.3"],
    ]
    path.write_text("#1.2\n4\t3\n" + "\n".join("\t".join(row) for row in [header, *rows]) + "\n")


def test_annotation_matches_ensembl_then_symbol_and_restricts_regulons(tmp_path: Path) -> None:
    gct = tmp_path / "median_tpm.gct"
    write_gct(gct)
    expression = load_brain_expression(gct)
    assert expression.brain_tissue_count.iloc[0] == 2
    nodes = pd.DataFrame([
        {"node_id": "GENE:DPM1", "node_type": "gene", "gene_symbol": "DPM1", "ensembl_gene_id": "ENSG00000000419"},
        {"node_id": "GENE:GCLC", "node_type": "gene", "gene_symbol": "GCLC", "ensembl_gene_id": "ENSG00000001084"},
        {"node_id": "GENE:TFHIGH", "node_type": "gene", "gene_symbol": "TFHIGH", "ensembl_gene_id": None},
        {"node_id": "GENE:TFLOW", "node_type": "gene", "gene_symbol": "TFLOW", "ensembl_gene_id": None},
        {"node_id": "GENE:UNKNOWN", "node_type": "gene", "gene_symbol": "UNKNOWN", "ensembl_gene_id": None},
        {"node_id": "MAR00001", "node_type": "reaction", "gene_symbol": None, "ensembl_gene_id": None},
    ])
    annotated = annotate_brain_expression(nodes, expression, minimum_tpm=1.0)
    by_id = annotated.set_index("node_id")
    assert by_id.loc["GENE:DPM1", "brain_expressed"] and by_id.loc["GENE:DPM1", "brain_median_tpm_max"] == 3.5
    assert by_id.loc["GENE:GCLC", "brain_expression_known"] and not by_id.loc["GENE:GCLC", "brain_expressed"]
    assert by_id.loc["GENE:TFHIGH", "brain_expressed"]  # symbol match for a gene without an Ensembl id
    assert not by_id.loc["GENE:UNKNOWN", "brain_expression_known"] and not by_id.loc["GENE:UNKNOWN", "brain_expressed"]
    edges = pd.DataFrame([
        {"source_id": "MAR00001", "target_id": "GENE:GCLC", "relation_type": "catalyzed_by", "sign": 0.0},
        {"source_id": "MAR00001", "target_id": "GENE:DPM1", "relation_type": "catalyzed_by", "sign": 0.0},
        {"source_id": "GENE:TFHIGH", "target_id": "GENE:DPM1", "relation_type": "regulates_transcription_of", "sign": 1.0},
        {"source_id": "GENE:TFLOW", "target_id": "GENE:DPM1", "relation_type": "regulates_transcription_of", "sign": 1.0},
        {"source_id": "GENE:UNKNOWN", "target_id": "GENE:DPM1", "relation_type": "regulates_transcription_of", "sign": -1.0},
        {"source_id": "GENE:TFLOW", "target_id": "GENE:DPM1", "relation_type": "activates", "sign": 1.0},
    ])
    kept, counts = restrict_transcription_edges_to_brain_expressed(annotated, edges)
    assert counts["transcription_edges_before"] == 3 and counts["transcription_edges_dropped_factor_not_brain_expressed"] == 1
    assert counts["transcription_edges_kept_factor_expression_unknown"] == 1 and counts["transcription_factors_dropped"] == 1
    assert len(kept) == 5 and not ((kept.source_id == "GENE:TFLOW") & (kept.relation_type == "regulates_transcription_of")).any()
    assert ((kept.source_id == "GENE:TFLOW") & (kept.relation_type == "activates")).any()  # only regulon edges are restricted
    with_reactions = propagate_brain_expression_to_reactions(annotated, kept)
    assert with_reactions.set_index("node_id").loc["MAR00001", "brain_median_tpm_max"] == 3.5  # largest over its genes


def test_structural_features_carry_brain_expression_when_present() -> None:
    import numpy as np

    from mechanistic_pathway_learning.evaluation.experiment_data import ExperimentData

    def make(brain_expression, brain_expressed):
        return ExperimentData(
            node_ids=["a", "b"], node_index={"a": 0, "b": 1}, node_types=np.array(["gene", "reaction"]), is_currency=np.array([False, False]),
            node_degree=np.array([1.0, 3.0]), edge_source=np.array([0]), edge_target=np.array([1]), edge_relation=np.array([0]), relation_types=["catalyzed_by"],
            symptoms=["anxiety"], perturbation_ids=["g"], perturbation_labels=["g"], perturbation_types=["gene"], group_ids=["g"],
            perturbation_seeds=[np.array([0])], perturbation_signs=[np.array([-1.0])], perturbation_magnitudes=[np.array([1.0])],
            outcomes=np.ones((1, 1)), weights=np.ones((1, 1)), in_metabolic_layer=np.array([True]),
            node_brain_expression=brain_expression, node_brain_expressed=brain_expressed,
        )

    without = make(None, None).structural_node_features()
    with_expression = make(np.log1p(np.array([3.5, 0.0])), np.array([True, False])).structural_node_features()
    assert without.shape == with_expression.shape  # the feature width is fixed so models built on either graph have the same input size
    assert np.all(without[:, -2:] == 0.0)
    assert np.isclose(with_expression[0, -2], np.log1p(3.5)) and with_expression[0, -1] == 1.0 and with_expression[1, -1] == 0.0
