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
)
from experiments.run_main_model_batch import CONFIGURATIONS
from mechanistic_pathway_learning.evaluation.experiment_data import load_experiment_data, restrict_to_perturbations
from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers
from mechanistic_pathway_learning.models.edge_gates import EdgeGates
from tests.test_drug_entry_nodes import DRUG_SEEDS, evidence_row, variant, write_graph
from tests.test_soft_priors_in_the_trainer import write_toy_dataset

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TIERS = REPOSITORY_ROOT / "configs/edge_prior_tiers.csv"


def write_edge_evidence(graph_directory: Path) -> pd.DataFrame:
    """A table for the toy graph: product_of edges curated, catalyzed_by edges curated, text-mined or without a tier in
    turn, every seventh edge contested."""
    edges = pd.read_parquet(graph_directory / "edges.parquet")
    rows = []
    for index, relation in enumerate(edges.relation_type):
        family, kind, records = ("curated", "curated", 1) if relation == "product_of" else [("signalling", "curated", 2), ("signalling", "text_mined_only", 1), ("other", "none", 0)][index % 3]
        rows.append({"family": family, "kind": kind, "records": records, "contested": index % 7 == 0})
    table = pd.DataFrame(rows).assign(source_id=edges.source_id, target_id=edges.target_id, relation_type=edges.relation_type)
    table.to_parquet(graph_directory / EDGE_EVIDENCE_FILE, index=False)
    return table


def test_the_edge_prior_is_off_by_default_and_registered_as_absent_before_it_existed():
    arguments = build_argument_parser().parse_args(["--graph-dir", "g", "--evidence-dir", "e"])
    assert arguments.edge_prior_strength == 0.0 and arguments.permute_edge_priors is False and arguments.edge_gate_learning_rate is None
    assert arguments.edge_prior_tiers == Path("configs/edge_prior_tiers.csv") and arguments.edge_prior_tempering == 0.5 and arguments.edge_prior_count_cap == 2
    assert VALUES_BEFORE_THE_V4_DEFAULTS["edge_prior_strength"] == 0.0 and VALUES_BEFORE_THE_V4_DEFAULTS["permute_edge_priors"] is False
    assert EDGE_PRIOR_SETTINGS == {"edge_prior_tiers", "edge_prior_tempering", "edge_prior_count_cap", "permute_edge_priors", "edge_gate_learning_rate"}
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
    plain = run_trainer(paths, tmp_path / "plain")
    off = run_trainer(paths, tmp_path / "off", "--edge-prior-strength", "0", "--permute-edge-priors", "--edge-gate-learning-rate", "0.5")
    assert np.array_equal(plain["predictions"], off["predictions"]), "strength 0 with every other edge flag set is the registered run, bit for bit"
    assert plain["results"]["configuration_fingerprint"] == off["results"]["configuration_fingerprint"] and off["results"]["edge_priors"] is None
    assert not (off["directory"] / "edge_gate_values.npy").exists()

    on = run_trainer(paths, tmp_path / "on", "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05")
    record = on["results"]["edge_priors"]
    gated, _, _ = edge_priors(table, read_edge_tiers(TIERS))
    assert record["gated"] == int(gated.sum()) and record["gates"]["gated_edges_the_encoder_reads"] == int(gated.sum()), "on the toy graph every gated edge is propagated along"
    assert set(record["penalty_denominator"]) == {"fit", "refit"} and record["penalty_denominator"]["refit"] > record["penalty_denominator"]["fit"] > 1.0
    assert record["gates"]["mean_absolute_move"] > 0 and record["gates"]["kl_penalty"] > 0
    assert all("edge_prior_penalty" in entry for entry in on["results"]["history"]) and all("edge_prior_penalty" in entry for entry in on["results"]["refit"]["history"])
    values = np.load(on["directory"] / "edge_gate_values.npy")
    assert np.isnan(values[~gated]).all() and ((values[gated] > 0) & (values[gated] < 1)).all()
    assert len(on["results"]["configuration_fingerprint"]["file:edge_evidence"]) == 64

    permuted = run_trainer(paths, tmp_path / "permuted", "--edge-prior-strength", "1.0", "--permute-edge-priors")
    assert permuted["results"]["edge_priors"]["permutation"]["seed"] == 0
    # the toy message-passing model predicts the same for every perturbation (its field does not reach the head here),
    # so whether the gates move predictions is read on the linear-response model, whose predictions do depend on its field
    linear_plain = run_trainer(paths, tmp_path / "linear_plain", "--encoder", "linear_response")
    assert np.ptp(linear_plain["predictions"], axis=0).max() > 0
    linear = run_trainer(paths, tmp_path / "linear", "--encoder", "linear_response", "--edge-prior-strength", "1.0", "--edge-gate-learning-rate", "0.05")
    assert linear["results"]["edge_priors"]["gates"]["gated_edges_the_encoder_reads"] == int(gated.sum())
    assert not np.array_equal(linear_plain["predictions"], linear["predictions"]), "gates that differ at a node change its shares and so the field"

    # the rewired-graph evaluation rebuilds the gated run from its checkpoint (gates attached before the weights load)
    evaluation = evaluate_fold(on["directory"], swaps_per_edge=2, rewiring_seed=0, data_cache={})
    assert evaluation["reproduced"], evaluation["reproduction_max_abs_difference"]
