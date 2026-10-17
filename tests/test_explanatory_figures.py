"""The edge layers the explanatory figures report (experiments/make_explanatory_figures.py, figures 1 and 10): each
source names one layer, the first rule that matches wins, and an electrical edge whose source names Reactome at its
end is electrical, not Reactome."""
from __future__ import annotations

import pandas as pd

from experiments.make_explanatory_figures import OTHER_CURATED_LAYER, edge_layers


def test_each_edge_source_names_one_layer():
    sources = {
        "Human-GEM": "metabolic (Human-GEM)", "Human-GEM GPR": "metabolic (Human-GEM)", "OmniPath": "signalling (OmniPath)",
        "OmniPath small_molecule_protein": "signalling (OmniPath)", "CollecTRI": "transcription (CollecTRI)",
        "Reactome entity member": "Reactome complexes, sets and regulators", "Reactome regulator": "Reactome complexes, sets and regulators",
        "BRENDA CC BY 4.0 (1 references, accession stated by BRENDA)": "enzyme regulation (BRENDA)",
        "net charge across the plasma membrane (Reactome)": "electrical (ions and charge)",
        "ion moved across the plasma membrane (Human-GEM)": "electrical (ions and charge)",
        "drug mechanism (ChEMBL), moved from the perturbation's seeds": "drug mechanism (ChEMBL)",
        "plasma carrier of the drug (both; label)": "drug plasma sequestration",
        "curated plasma carriage": "plasma carriage of metabolites",
        "neomorphic allele mechanism (PMID 11264541)": "neomorphic allele mechanism",
        "glutathione consumed": OTHER_CURATED_LAYER,
    }
    edges = pd.DataFrame({"evidence_source": list(sources), "relation_type": "activates"})
    assert list(edge_layers(edges)) == list(sources.values())
    encodes = pd.DataFrame({"evidence_source": ["gene_protein_split", None], "relation_type": ["encodes", "encodes"]})
    assert list(edge_layers(encodes)) == ["gene to protein (encodes)", "gene to protein (encodes)"], "the relation decides, whatever the source"


def test_the_share_figure_is_computed_by_the_encoders_and_a_closed_gate_is_the_graph_without_its_edge():
    """Figure 14's numbers: linear response takes each input's share of the gated weights at the target, message
    passing divides by the gated in-degree; with b's gate at 0 both equal the graph without b; a sole input's value is 1
    for any gate above 0 and 0 at 0, and one of three equal inputs takes g / (g + 2)."""
    import numpy as np

    from experiments.make_explanatory_figures import share_reading_values

    values = share_reading_values()
    linear, message = values["linear_response"], values["message_passing"]
    assert np.allclose(linear["every_gate_open"], [0.25, 0.5, 0.25])
    assert np.allclose(linear["b_closed"], [0.5, 0.0, 0.5]) and np.allclose(linear["b_closed"], linear["graph_without_b"])
    assert np.allclose(message["every_gate_open"], [1 / 3, 2 / 3, 1 / 3])
    assert np.allclose(message["b_closed"], [0.5, 0.0, 0.5]) and np.allclose(message["b_closed"], message["graph_without_b"])
    grid = values["gate_grid"]
    assert values["sole_input"][0] == 0.0 and np.allclose(values["sole_input"][1:], 1.0), "a sole input is a cliff at 0"
    assert np.allclose(values["one_of_three_inputs"], grid / (grid + 2.0))


def test_the_wilson_interval_matches_known_values():
    from experiments.make_explanatory_figures import wilson_interval

    low, high = wilson_interval(0, 10)
    assert abs(low) < 1e-12 and abs(high - 0.2775) < 1e-3
    low, high = wilson_interval(5, 10)
    assert abs(low - 0.2366) < 1e-3 and abs(high - 0.7634) < 1e-3
    assert all(value != value for value in wilson_interval(0, 0)), "an empty total has no interval"


def simulated_run(closed: int, false: int, false_edges: int = 40, edges: int = 200, rows=None) -> dict:
    """A digest entry with the fields the simulation figures read."""
    precision = false / closed if closed else None
    return {"world": {"edges": edges, "false_edges": false_edges, "nodes": 30, "perturbations": 80, "symptoms": 4},
            "fit": {"closed": closed, "false_edges": false_edges, "precision": precision, "recall": false / false_edges},
            "closures_by_evidence": {f"at least {threshold} nats": {"closed": closed // 2, "false_edges": false_edges, "precision": 1.0 if closed // 2 else None}
                                     for threshold in ("0.0", "0.5", "1.0", "2.0", "3.0")},
            "closures_by_posterior_false": {f"at least {level}": {"closed": 1, "false": 1, "precision": 1.0} for level in ("0.5", "0.7", "0.8", "0.9", "0.95")},
            "stability": {"thresholds": {threshold: {"closed": 0, "false_edges": false_edges, "precision": None} for threshold in ("0.6", "0.7", "0.8", "0.9")}},
            "held_out_loss": {"fitted_graph": 0.5, "every_edge_open": 0.52, "true_graph": 0.48},
            "closed_edges": rows or [], "closed_edges_columns": ["tier", "prior_real", "real", "evidence_nats", "posterior_false"]}


def test_selector_points_pool_the_false_closures_of_every_world():
    from experiments.make_explanatory_figures import pooled_selector_points

    runs = {"one": simulated_run(10, 6), "two": simulated_run(30, 12)}
    points = pooled_selector_points(runs, ("one", "two"))
    read_out = points["read-out"][0]
    assert read_out["closed"] == 40 and abs(read_out["precision"] - 18 / 40) < 1e-12 and read_out["closures_per_world"] == 20
    assert read_out["low"] < read_out["precision"] < read_out["high"]
    assert points["stability"] == [], "a threshold with no closure in any world draws no point"


def test_posterior_calibration_bins_the_closed_edges_by_their_posterior():
    from experiments.make_explanatory_figures import posterior_calibration

    rows = [["text_mined_only", 0.525, 0, 0.5, 0.55], ["text_mined_only", 0.525, 1, 0.4, 0.57], ["curated", 0.9, 0, 5.0, 0.95], ["curated", 0.9, 0, 6.0, 1.0]]
    bins = posterior_calibration({"one": simulated_run(4, 3, rows=rows)}, ("one",))
    assert [(item["low"], item["edges"]) for item in bins] == [(0.5, 2), (0.9, 2)], "the last bin includes P(false) = 1"
    assert abs(bins[0]["share_false"] - 0.5) < 1e-12 and bins[1]["share_false"] == 1.0
    assert abs(bins[0]["mean_posterior_false"] - 0.56) < 1e-12


def test_a_condition_pools_its_runs_and_reports_the_share_of_the_attainable_held_out_gain():
    from experiments.make_explanatory_figures import condition_summary

    runs = {"one": simulated_run(10, 6), "two": simulated_run(30, 12)}
    summary = condition_summary(runs, ("one", "two", "absent"))
    assert summary["runs"] == 2 and summary["closed"] == 40 and abs(summary["precision"] - 0.45) < 1e-12
    assert abs(summary["recall"] - 18 / 80) < 1e-12 and abs(summary["base_rate"] - 0.2) < 1e-12
    assert abs(summary["held_out_gain_share"] - 0.5) < 1e-12, "(0.52 - 0.50) / (0.52 - 0.48)"
    assert summary["evidence_closed"] == 20 and summary["evidence_precision"] == 1.0
    assert condition_summary(runs, ("absent",)) is None


def test_the_digest_reads_a_closure_as_false_with_the_posterior_its_evidence_and_prior_give():
    """P(false) = sigmoid(E - logit pi): evidence equal to the closing cost leaves it at one half, and each further nat
    multiplies the odds of false by e."""
    import math

    from experiments.summarise_edge_gate_simulations import posterior_false, selected_scores

    prior = 0.9
    assert abs(posterior_false(math.log(prior / (1 - prior)), prior) - 0.5) < 1e-12
    assert abs(posterior_false(math.log(prior / (1 - prior)) + 1.0, prior) - 1 / (1 + math.exp(-1.0))) < 1e-12
    edges = [{"real": False, "posterior_false": 0.95}, {"real": True, "posterior_false": 0.92}, {"real": False, "posterior_false": 0.6}]
    assert selected_scores(edges, lambda edge: edge["posterior_false"] >= 0.9) == {"closed": 2, "false": 1, "precision": 0.5}
