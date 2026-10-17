"""The edge priors in the loader and the trainer (decision 9.57): off by default and then invisible, the evidence table
read through the loader's edge rows and refused when it is another graph's, the gates in an optimizer group of their own
without weight decay, the penalty and the gates recorded in results.json, the permutation control, and a gated run that
the rewired-graph evaluation rebuilds from its checkpoint."""
from __future__ import annotations

import json
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import torch

from experiments.evaluate_on_rewired_graph import evaluate_fold
from experiments.run_main_model import (
    ARGUMENTS_RECORDED_ONLY_WHEN_SET,
    EDGE_EVIDENCE_FILE,
    EDGE_PRIOR_SETTINGS,
    VALUES_BEFORE_THE_V4_DEFAULTS,
    build_argument_parser,
    build_edge_priors,
    configuration_fingerprint,
    optimizer_parameter_groups,
    write_closure_review,
)
from experiments.run_main_model_batch import CONFIGURATIONS
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, restrict_to_perturbations
from mechanistic_pathway_learning.graph.edge_evidence import curated_edge_labels, edge_priors, read_edge_tiers
from mechanistic_pathway_learning.models.edge_gates import EdgeGates
from tests.test_drug_entry_nodes import DRUG_SEEDS, evidence_row, variant, write_graph
from tests.test_soft_priors_in_the_trainer import write_toy_dataset

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TIERS = REPOSITORY_ROOT / "configs/edge_prior_tiers.csv"


def write_edge_evidence(graph_directory: Path) -> pd.DataFrame:
    """A table for the toy graph: product_of edges curated from one paper, catalyzed_by edges curated from two papers,
    text-mined or without a tier in turn, every seventh edge contested; every other two-paper curated edge has a human
    paper among its two, so it is an edge label under decision 9.63 (b) unless contested."""
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    rows = []
    for index, relation in enumerate(edges.relation_type):
        family, kind, records = ("curated", "curated", 1) if relation == "product_of" else [("signalling", "curated", 2), ("signalling", "text_mined_only", 1), ("other", "none", 0)][index % 3]
        pubmed_ids = ";".join(str(1000 + 10 * index + offset) for offset in range(records))
        human = int(records == 2 and (index // 3) % 2 == 0)
        rows.append({"family": family, "kind": kind, "records": records, "contested": index % 7 == 0, "pubmed_ids": pubmed_ids,
                     "resources": "SIGNOR" if family == "signalling" else "", "species_known_records": records, "human_records": human,
                     "animal_only_records": records - human, "pubmed_species": ""})
    table = pd.DataFrame(rows).assign(source_id=edges.source_id, target_id=edges.target_id, relation_type=edges.relation_type)
    table.to_parquet(graph_directory / EDGE_EVIDENCE_FILE, index=False)
    return table


def test_the_edge_prior_is_off_by_default_and_registered_as_absent_before_it_existed():
    arguments = build_argument_parser().parse_args(["--graph-dir", "g", "--evidence-dir", "e"])
    assert arguments.edge_prior_strength == 0.0 and arguments.permute_edge_priors is False and arguments.edge_gate_learning_rate is None
    assert arguments.edge_prior_tiers == Path("configs/edge_prior_tiers.csv") and arguments.edge_prior_tempering == 0.5 and arguments.edge_prior_count_cap == 2
    assert VALUES_BEFORE_THE_V4_DEFAULTS["edge_prior_strength"] == 0.0 and VALUES_BEFORE_THE_V4_DEFAULTS["permute_edge_priors"] is False
    assert EDGE_PRIOR_SETTINGS == {"edge_prior_tiers", "edge_prior_tempering", "edge_prior_count_cap", "permute_edge_priors", "edge_gate_learning_rate",
                                   "edge_gate_kind", "edge_gate_beta", "edge_gate_binary_epochs", "edge_gate_binary_learning_rate", "edge_gate_surrogate_width",
                                   "edge_gate_binary_full_batch", "edge_gate_closure_evidence", "edge_gate_flip_checks", "edge_gate_binary_learning_rate_decay",
                                   "curated_edge_labels", "curated_label_minimum_records", "curated_label_minimum_human_records",
                                   "curated_label_minimum_sign_records", "curated_label_exclude_disputed_signs", "curated_label_bulk_papers"}
    assert arguments.curated_edge_labels is True and arguments.curated_label_minimum_records == 2 and arguments.curated_label_minimum_human_records == 1
    assert arguments.curated_label_minimum_sign_records is None, "the sign-evidence rule of decision 9.70 is off unless the owner turns it on"
    assert arguments.curated_label_exclude_disputed_signs is False, "the disputed-sign rule of decision 9.72 is off unless the owner turns it on"
    assert arguments.curated_label_bulk_papers == "count", "a bulk paper counts one record, as in 9.63, unless the owner says otherwise"
    assert arguments.edge_gate_closure_evidence == 20, "the exact evidence for the 20 closures against the strongest priors; every closure gets the estimate"
    assert "edge_prior_strength" in ARGUMENTS_RECORDED_ONLY_WHEN_SET
    for name in ("v4_linear_response_edge_priors", "v4_linear_response_edge_priors_permuted", "v4_message_passing_edge_priors", "v4_message_passing_edge_priors_permuted"):
        flags = CONFIGURATIONS[name]
        assert flags[flags.index("--edge-prior-strength") + 1] == "1.0" and (("--permute-edge-priors" in flags) == name.endswith("_permuted"))


def test_the_fingerprint_names_the_edge_prior_only_while_it_is_on(tmp_path):
    data = Namespace(outcomes=np.zeros((2, 2)), weights=np.zeros((2, 2)), edge_source=np.array([0]), edge_target=np.array([1]), edge_relation=np.array([0]),
                     perturbation_seeds=[np.array([0]), np.array([1])], perturbation_signs=[np.array([-1.0]), np.array([-1.0])], perturbation_magnitudes=[np.array([1.0]), np.array([1.0])])
    (tmp_path / EDGE_EVIDENCE_FILE).write_bytes(b"table")
    base = dict(graph_dir=tmp_path, edge_prior_tiers=TIERS, edge_prior_tempering=0.5, edge_prior_count_cap=2, permute_edge_priors=True, edge_gate_learning_rate=0.1)
    off = configuration_fingerprint(Namespace(**base, edge_prior_strength=0.0), data, None)
    assert not [key for key in off if "edge_prior" in key or "edge_gate" in key or "edge_evidence" in key], "at strength 0 an old run's fingerprint stays as it was"
    on = configuration_fingerprint(Namespace(**base, edge_prior_strength=1.0), data, None)
    assert on["argument:edge_prior_strength"] == "1.0" and on["argument:permute_edge_priors"] == "True" and on["argument:edge_gate_learning_rate"] == "0.1"
    assert len(on["file:edge_prior_tiers"]) == 64 and len(on["file:edge_evidence"]) == 64


def test_the_loader_records_each_edges_row_in_the_graphs_edge_table(tmp_path):
    """drug_entry "targets" drops the rows the drug-entry build added; with those rows interleaved, the kept edges are a
    scattered subset and each must name its own row."""
    _, _, new_nodes, new_edges, relations, entry_table, _ = variant()
    shuffled = new_edges.sample(frac=1.0, random_state=3).reset_index(drop=True)
    write_graph(tmp_path / "drugs", new_nodes, shuffled, relations, entry_table)
    (tmp_path / "evidence").mkdir()
    rows = [evidence_row("knockout_B", "gene", "cluster:B", [["GENE:B", -1.0, 1.0]])]
    rows += [evidence_row(perturbation, "drug", perturbation, [list(triple) for triple in triples]) for perturbation, (_, triples) in DRUG_SEEDS.items()]
    pd.DataFrame(rows).to_parquet(tmp_path / "evidence/evidence_records.parquet", index=False)
    as_nodes = load_experiment_data(tmp_path / "drugs", tmp_path / "evidence", drug_entry="nodes")
    as_targets = load_experiment_data(tmp_path / "drugs", tmp_path / "evidence", drug_entry="targets")
    assert as_nodes.edge_row.tolist() == list(range(len(shuffled)))
    assert len(as_targets.edge_row) < len(shuffled) and as_targets.edge_row.tolist() != list(range(len(as_targets.edge_row))), "a scattered subset"
    for position, row in enumerate(as_targets.edge_row):
        assert shuffled.source_id[row] == as_targets.node_ids[as_targets.edge_source[position]] and shuffled.target_id[row] == as_targets.node_ids[as_targets.edge_target[position]]
    kept = np.zeros(len(as_targets.perturbation_ids), dtype=bool)
    kept[:2] = True
    assert np.array_equal(restrict_to_perturbations(as_targets, kept).edge_row, as_targets.edge_row), "the edges are the graph's, not the perturbations'"


def priors_arguments(paths: dict, **changes) -> Namespace:
    return Namespace(**{**dict(edge_prior_strength=2.0, encoder="message_passing", graph_dir=paths["graph"], edge_prior_tiers=TIERS, edge_prior_tempering=0.5,
                               edge_prior_count_cap=2, permute_edge_priors=False, seed=0), **changes})


def test_build_edge_priors_reads_the_table_through_the_loaders_rows_and_refuses_another_graphs(tmp_path):
    paths = write_toy_dataset(tmp_path, num_genes=12)
    table = write_edge_evidence(paths["graph"])
    data = load_experiment_data(paths["graph"], paths["evidence"])
    summary = build_edge_priors(data, priors_arguments(paths))
    gated, probability, strength = edge_priors(table, read_edge_tiers(TIERS), global_strength=2.0)
    np.testing.assert_array_equal(data.edge_prior_gated, gated)
    np.testing.assert_array_equal(data.edge_prior_probability, probability)
    np.testing.assert_array_equal(data.edge_prior_strength, strength)
    assert summary["gated"] == int(gated.sum()) and summary["strength"] == 2.0 and summary["permutation"] is None and len(summary["sha256"]) == 64
    assert set(summary["by_family_and_kind"]) == {"curated:curated", "signalling:curated", "signalling:text_mined_only", "other:none"}

    assert build_edge_priors(data, Namespace(edge_prior_strength=0.0)) is None and data.edge_prior_gated is None, "off clears what an earlier configuration set"
    permuted = build_edge_priors(data, priors_arguments(paths, permute_edge_priors=True))
    assert permuted["permutation"]["seed"] == 0 and len(permuted["permutation"]["priors_sha256"]) == 64
    for relation in np.unique(data.edge_relation):
        of_relation = gated & (data.edge_relation == relation)
        assert sorted(data.edge_prior_probability[of_relation]) == sorted(probability[of_relation]), "shuffled within the relation"
    assert not np.array_equal(data.edge_prior_probability[gated], probability[gated])

    with pytest.raises(ValueError, match="the none encoder reads none"):
        build_edge_priors(data, priors_arguments(paths, encoder="none"))
    table.iloc[::-1].to_parquet(paths["graph"] / EDGE_EVIDENCE_FILE, index=False)
    with pytest.raises(ValueError, match="row for row"):
        build_edge_priors(data, priors_arguments(paths))
    table.drop(columns=["source_id"]).to_parquet(paths["graph"] / EDGE_EVIDENCE_FILE, index=False)
    with pytest.raises(ValueError, match="before the endpoint columns"):
        build_edge_priors(data, priors_arguments(paths))
    (paths["graph"] / EDGE_EVIDENCE_FILE).unlink()
    with pytest.raises(FileNotFoundError, match="build_edge_evidence.py"):
        build_edge_priors(data, priors_arguments(paths))


def test_curated_edges_with_human_evidence_are_taken_as_known_and_the_permutation_leaves_them(tmp_path):
    paths = write_toy_dataset(tmp_path, num_genes=12)
    table = write_edge_evidence(paths["graph"])
    data = load_experiment_data(paths["graph"], paths["evidence"])
    summary = build_edge_priors(data, priors_arguments(paths, curated_edge_labels=True, curated_label_minimum_records=2, curated_label_minimum_human_records=1))
    gated, probability, strength = edge_priors(table, read_edge_tiers(TIERS), global_strength=2.0)
    labels = curated_edge_labels(table)
    assert labels.any() and summary["curated_edge_labels"]["edges"] == int(labels.sum()) == summary["curated_edge_labels"]["were_gated"]
    np.testing.assert_array_equal(data.edge_prior_gated, gated & ~labels)
    assert np.isnan(data.edge_prior_probability[labels]).all() and (data.edge_prior_strength[labels] == 0).all(), "a label carries no prior: it is known"
    np.testing.assert_array_equal(data.edge_prior_probability[gated & ~labels], probability[gated & ~labels])
    assert summary["gated"] == int((gated & ~labels).sum())
    assert len(data.edge_evidence) == len(data.edge_source) and data.edge_evidence.pubmed_ids.tolist() == table.pubmed_ids.tolist()
    stricter = build_edge_priors(data, priors_arguments(paths, curated_edge_labels=True, curated_label_minimum_records=3, curated_label_minimum_human_records=1))
    assert stricter["curated_edge_labels"]["edges"] == 0 and stricter["gated"] == int(gated.sum())
    permuted = build_edge_priors(data, priors_arguments(paths, curated_edge_labels=True, curated_label_minimum_records=2, curated_label_minimum_human_records=1,
                                                        permute_edge_priors=True))
    assert permuted["curated_edge_labels"]["edges"] == int(labels.sum()) and np.isnan(data.edge_prior_probability[labels]).all(), "the labels stay on their edges"
    table.drop(columns=["human_records"]).to_parquet(paths["graph"] / EDGE_EVIDENCE_FILE, index=False)
    with pytest.raises(ValueError, match="no-curated-edge-labels"):
        build_edge_priors(data, priors_arguments(paths, curated_edge_labels=True))
    assert build_edge_priors(data, priors_arguments(paths))["curated_edge_labels"] is None, "without the flag an old table still serves"


def test_a_label_whose_sign_the_collectri_releases_dispute_is_gated_only_when_asked(tmp_path):
    paths = write_toy_dataset(tmp_path, num_genes=12)
    table = write_edge_evidence(paths["graph"])
    labels = curated_edge_labels(table)
    disputed = np.zeros(len(table), dtype=bool)
    disputed[np.flatnonzero(labels)[::2]] = True
    table.assign(sign_disputed=disputed).to_parquet(paths["graph"] / EDGE_EVIDENCE_FILE, index=False)
    data = load_experiment_data(paths["graph"], paths["evidence"])
    label_arguments = dict(curated_edge_labels=True, curated_label_minimum_records=2, curated_label_minimum_human_records=1)
    plain = build_edge_priors(data, priors_arguments(paths, **label_arguments))
    assert plain["curated_edge_labels"]["edges"] == int(labels.sum()) and plain["curated_edge_labels"]["exclude_disputed_signs"] is False
    excluded = build_edge_priors(data, priors_arguments(paths, **label_arguments, curated_label_exclude_disputed_signs=True))
    assert excluded["curated_edge_labels"]["edges"] == int((labels & ~disputed).sum()) and excluded["curated_edge_labels"]["exclude_disputed_signs"] is True
    assert data.edge_prior_gated[labels & disputed].all(), "a disputed regulation goes back under its gate, at its tier's prior"
    table.to_parquet(paths["graph"] / EDGE_EVIDENCE_FILE, index=False)
    with pytest.raises(ValueError, match="sign_disputed"):
        build_edge_priors(data, priors_arguments(paths, **label_arguments, curated_label_exclude_disputed_signs=True))


def test_the_review_list_holds_every_closure_with_its_evidence(tmp_path):
    evidence = pd.DataFrame({"family": ["signalling", "transcription", "signalling"], "kind": ["curated", "text_mined_only", "unreferenced"],
                             "records": [1, 3, 0], "contested": [False, False, True], "species_known_records": [1, 3, 0], "human_records": [0, 1, 0],
                             "animal_only_records": [1, 0, 0], "pubmed_ids": ["11", "21;22;23", ""], "pubmed_species": ["11:0.00/0.50/0.50", "", ""],
                             "resources": ["SIGNOR", "ExTRI_CollecTRI", ""]})
    data = Namespace(node_ids=np.array(["A", "B", "C"], dtype=object), edge_source=np.array([0, 1, 2]), edge_target=np.array([1, 2, 0]),
                     edge_relation=np.array([0, 1, 0]), relation_types=["activates", "regulates_transcription_of"], edge_evidence=evidence)
    binary_pass = {"closed_edges": [{"edge": 2, "prior": 0.6}, {"edge": 0, "prior": 0.9}, {"edge": 1, "prior": 0.53}],
                   "closure_evidence": {"closed": 3, "evaluated": 1, "records": [{"gate": 0, "edge": 0, "prior": 0.9, "evidence_nats": 3.1, "posterior_false": 0.71}]}}
    summary = write_closure_review(binary_pass, data, tmp_path)
    table = pd.read_csv(tmp_path / summary["file"], keep_default_na=True)
    assert table.edge.tolist() == [0, 2, 1], "against the strongest priors first"
    assert summary["closed"] == 3 and summary["evaluated"] == 1 and summary["evidence_at_least_1_nat"] == 1
    assert table.evidence_nats.isna().tolist() == [False, True, True], "a closure left unmeasured is still listed"
    first = table.iloc[0]
    assert (first.source_id, first.target_id, first.relation, first.family, first.records, first.pubmed_ids, first.resources) == ("A", "B", "activates", "signalling", 1, "11", "SIGNOR")
    assert table.iloc[2].pubmed_ids == "21;22;23" and table.iloc[2].human_records == 1
    json.dumps(summary)
    assert summary["by_kind"] == {"signalling:curated": 1, "signalling:unreferenced": 1, "transcription:text_mined_only": 1}


class GatedEncoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.weight = torch.nn.Parameter(torch.zeros(2))
        self.edge_gates = EdgeGates(np.array([True, False]), np.array([0.8, np.nan]), np.array([1.0, 0.0]))


class HeadWithTimeScales(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.link = torch.nn.Parameter(torch.zeros(1))
        self.rest = torch.nn.Parameter(torch.zeros(1))

    def time_scale_parameter_groups(self):
        return {"links": [self.link], "leaks": [], "module_biases": [], "gates": []}


def test_the_edge_gates_are_an_optimizer_group_of_their_own_without_weight_decay():
    rates = dict(link_learning_rate=0.0, leak_learning_rate=0.0, module_bias_learning_rate=0.0, gate_learning_rate=0.0, descriptor_learning_rate=0.0, learning_rate=0.01)
    encoder = GatedEncoder()
    groups = optimizer_parameter_groups(encoder, torch.nn.Linear(2, 1), Namespace(**rates, edge_gate_learning_rate=None))
    assert len(groups) == 2 and groups[-1] == {"params": [encoder.edge_gates.logit], "lr": 0.01, "weight_decay": 0.0}
    assert all(parameter is not encoder.edge_gates.logit for parameter in groups[0]["params"]) and len(groups[0]["params"]) == 3
    groups = optimizer_parameter_groups(encoder, HeadWithTimeScales(), Namespace(**{**rates, "link_learning_rate": 0.02}, edge_gate_learning_rate=0.5))
    assert [group.get("lr") for group in groups] == [None, 0.02, 0.5] and groups[-1]["params"] == [encoder.edge_gates.logit]
    assert len(optimizer_parameter_groups(torch.nn.Linear(2, 2), torch.nn.Linear(2, 1), Namespace(**rates))) == 1, "without gates the groups are as they were"


def run_trainer(paths: dict[str, Path], run_directory: Path, *extra: str) -> dict:
    command = [sys.executable, "experiments/run_main_model.py", "--graph-dir", str(paths["graph"]), "--evidence-dir", str(paths["evidence"]), "--run-dir", str(run_directory),
               "--encoder", "message_passing", "--max-epochs", "2", "--num-bootstrap", "3", "--refit-on-validation", "--validation-fraction", "0.3",
               "--no-keep-large-groups-in-training", "--fold", "0", "--seed", "0", "--batch-size", "8", "--placebo-negative-weights-file", "none",
               "--group-by", "disease_cluster", "--label-selection", str(paths["selection"]), *extra]
    completed = subprocess.run(command, cwd=REPOSITORY_ROOT, capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin", "OMP_NUM_THREADS": "1"})
    assert completed.returncode == 0, completed.stderr[-3000:]
    split_directory = next(path for path in run_directory.iterdir() if path.is_dir())
    return {"directory": split_directory, "results": json.loads((split_directory / "results.json").read_text()),
            "predictions": np.load(split_directory / "test_predictions.npy")}


def test_end_to_end_the_edge_prior_off_is_the_plain_run_and_on_it_gates_trains_and_records(tmp_path):
    paths = write_toy_dataset(tmp_path)
    table = write_edge_evidence(paths["graph"])
    labels = curated_edge_labels(table)
    assert labels.any() and not labels.all()
    plain = run_trainer(paths, tmp_path / "plain")
    off = run_trainer(paths, tmp_path / "off", "--edge-prior-strength", "0", "--permute-edge-priors", "--edge-gate-learning-rate", "0.5")
    assert np.array_equal(plain["predictions"], off["predictions"]), "strength 0 with every other edge flag set is the registered run, bit for bit"
    assert plain["results"]["configuration_fingerprint"] == off["results"]["configuration_fingerprint"] and off["results"]["edge_priors"] is None
    assert not (off["directory"] / "edge_gate_values.npy").exists()

    on = run_trainer(paths, tmp_path / "on", "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05", "--edge-gate-kind", "soft")
    record = on["results"]["edge_priors"]
    gated, _, _ = edge_priors(table, read_edge_tiers(TIERS))
    gated &= ~labels  # the curated edges taken as known carry no gate (decision 9.63 (b), on by default with edge priors)
    assert record["curated_edge_labels"]["edges"] == int(labels.sum())
    assert record["gated"] == int(gated.sum()) and record["gates"]["gated_edges_the_encoder_reads"] == int(gated.sum()), "on the toy graph every gated edge is propagated along"
    assert set(record["penalty_denominator"]) == {"fit", "refit"} and record["penalty_denominator"]["refit"] > record["penalty_denominator"]["fit"] > 1.0
    assert record["gates"]["mean_absolute_move"] > 0 and record["gates"]["kl_penalty"] > 0
    assert all("edge_prior_penalty" in entry for entry in on["results"]["history"]) and all("edge_prior_penalty" in entry for entry in on["results"]["refit"]["history"])
    values = np.load(on["directory"] / "edge_gate_values.npy")
    assert np.isnan(values[~gated]).all() and ((values[gated] > 0) & (values[gated] < 1)).all()
    assert len(on["results"]["configuration_fingerprint"]["file:edge_evidence"]) == 64

    permuted = run_trainer(paths, tmp_path / "permuted", "--edge-prior-strength", "1.0", "--permute-edge-priors", "--edge-gate-kind", "soft")
    assert permuted["results"]["edge_priors"]["permutation"]["seed"] == 0
    # the toy message-passing model predicts the same for every perturbation (its field does not reach the head here),
    # so whether the gates move predictions is read on the linear-response model, whose predictions do depend on its field
    linear_plain = run_trainer(paths, tmp_path / "linear_plain", "--encoder", "linear_response")
    assert np.ptp(linear_plain["predictions"], axis=0).max() > 0
    linear = run_trainer(paths, tmp_path / "linear", "--encoder", "linear_response", "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05",
                         "--edge-gate-kind", "soft")
    assert linear["results"]["edge_priors"]["gates"]["gated_edges_the_encoder_reads"] == int(gated.sum())
    assert not np.array_equal(linear_plain["predictions"], linear["predictions"]), "gates that differ at a node change its shares and so the field"

    # the rewired-graph evaluation rebuilds the gated run from its checkpoint (gates attached before the weights load)
    evaluation = evaluate_fold(on["directory"], swaps_per_edge=2, rewiring_seed=0, data_cache={})
    assert evaluation["reproduced"], evaluation["reproduction_max_abs_difference"]


def test_end_to_end_hard_gates_train_run_the_binary_pass_and_save_the_scored_read_out(tmp_path):
    """Decision 9.62 (e): with the default gate kind the run samples hard concrete gates while it fits, runs the binary
    pass and the exact checks, scores the test with the binary read-out, saves it, and leaves the saved state at the
    scored model, so the rewired-graph evaluation reproduces the predictions."""
    paths = write_toy_dataset(tmp_path)
    table = write_edge_evidence(paths["graph"])
    gated, _, _ = edge_priors(table, read_edge_tiers(TIERS))
    labels = curated_edge_labels(table)
    gated &= ~labels
    for encoder in ("linear_response", "message_passing"):
        run = run_trainer(paths, tmp_path / encoder, "--encoder", encoder, "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05",
                          "--edge-gate-binary-epochs", "2", "--edge-gate-flip-checks", "3", "--edge-gate-closure-evidence", "50")
        record = run["results"]["edge_priors"]
        assert record["gates"]["kind"] == "hard_concrete" and record["gates"]["phase"] == "binary"
        assert len(record["binary_pass"]["history"]) == 2 and record["binary_pass"]["flip_checks"]["checked"] <= 3
        evidence = record["binary_pass"]["closure_evidence"]
        closures = pd.read_csv(run["directory"] / evidence["file"])
        assert evidence["evaluated"] == len(closures) == min(evidence["closed"], 50) == record["gates"]["closed_at_read_out"], "every closed edge measured"
        if len(closures):
            expected = 1 / (1 + np.exp(-(closures.evidence_nats - np.log(closures.prior / (1 - closures.prior)))))
            assert np.allclose(closures.posterior_false, expected), "P(false) = sigmoid(evidence - logit prior)"
            assert closures.prior.is_monotonic_decreasing, "the closures against the strongest priors are measured first"
            assert not closures.edge.isin(np.flatnonzero(labels)).any(), "an edge taken as known cannot close"
            assert {"family", "kind", "records", "pubmed_ids", "human_records", "resources"} <= set(closures.columns), "each closure with its evidence (9.63 (c))"
            assert closures.estimated_evidence_nats.notna().all() and closures.estimated_posterior_false.between(0, 1).all(), "every closure has its estimate"
        read_out = np.load(run["directory"] / "edge_gate_read_out.npy")
        assert set(np.unique(read_out)) <= {0, 1} and read_out[~gated].all(), "an ungated edge is always open"
        assert int((read_out[gated] == 0).sum()) == record["gates"]["closed_at_read_out"]
        probability = np.load(run["directory"] / "edge_gate_values.npy")
        assert np.isnan(probability[~gated]).all() and ((probability[gated] > 0) & (probability[gated] < 1)).all()
        evaluation = evaluate_fold(run["directory"], swaps_per_edge=2, rewiring_seed=0, data_cache={})
        assert evaluation["reproduced"], (encoder, evaluation["reproduction_max_abs_difference"])



def test_the_binary_pass_anneals_its_learning_rate_when_asked_and_keeps_it_by_default(tmp_path):
    """--edge-gate-binary-learning-rate-decay multiplies the binary pass's rate after each epoch, as the simulations'
    --binary-learning-rate-decay does (batch 4 of docs/edge_gate_simulation.md: 0.8 settled the pass in every world);
    the default 1 keeps the rate the trainer has always used."""
    paths = write_toy_dataset(tmp_path)
    write_edge_evidence(paths["graph"])
    common = ("--encoder", "linear_response", "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05", "--edge-gate-binary-epochs", "3",
              "--edge-gate-flip-checks", "0", "--edge-gate-closure-evidence", "0", "--edge-gate-binary-learning-rate", "0.1")
    annealed = run_trainer(paths, tmp_path / "annealed", *common, "--edge-gate-binary-learning-rate-decay", "0.5")["results"]["edge_priors"]["binary_pass"]
    assert annealed["learning_rate_decay"] == 0.5
    assert [round(epoch["learning_rate"], 6) for epoch in annealed["history"]] == [0.1, 0.05, 0.025]
    kept = run_trainer(paths, tmp_path / "kept", *common)["results"]["edge_priors"]["binary_pass"]
    assert kept["learning_rate_decay"] == 1.0 and [round(epoch["learning_rate"], 6) for epoch in kept["history"]] == [0.1, 0.1, 0.1]


def test_the_timing_mode_runs_with_hard_gates(tmp_path):
    """--timing-batches times training steps outside the training loop; with hard gates each step must draw its own
    sample, or the second backward runs through the first step's freed graph (found on the v5 graph, 17 October 2026)."""
    paths = write_toy_dataset(tmp_path)
    write_edge_evidence(paths["graph"])
    command = [sys.executable, "experiments/run_main_model.py", "--graph-dir", str(paths["graph"]), "--evidence-dir", str(paths["evidence"]),
               "--run-dir", str(tmp_path / "timing"), "--encoder", "linear_response", "--fold", "0", "--seed", "0", "--batch-size", "8",
               "--placebo-negative-weights-file", "none", "--group-by", "disease_cluster", "--label-selection", str(paths["selection"]),
               "--edge-prior-strength", "1.0", "--timing-batches", "3"]
    completed = subprocess.run(command, cwd=REPOSITORY_ROOT, capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin", "OMP_NUM_THREADS": "1"})
    assert completed.returncode == 0, completed.stderr[-3000:]
    timing = json.loads(next(line for line in completed.stdout.splitlines() if line.startswith("{") and "median_step_seconds" in line))
    assert timing["timing_batches"] == 3 and timing["median_step_seconds"] >= 0
