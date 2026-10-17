"""The digest of the edge-gate simulations (decision 9.62) that docs/edge_gate_simulation.md and figures 15 and 16 read.

Reads runs/simulations/<name>.json for the runs named in RUNS (scripts/run_edge_gate_simulations.sh writes them; the
directory is not tracked) and writes docs/edge_gate_simulation_results.json: per run, its settings and world, the
scores of the read-out, of the exact closure evidence and of stability selection, the held-out losses, and edge by edge
the closed edges with their tier, prior, truth and evidence, the closing frequencies over refits and the first-order
estimate beside the exact change. A run whose file is missing is left out and listed under "missing".

The v5 graph's labels and gated edges are added under "v5" when data/processed holds the tables, so that figure 16 can
set the simulated worlds beside the graph the gates will be fitted on; without them that entry is null.

The posterior that a closed edge is false. At a local optimum under single flips, a closed edge e is one whose
reopening alone would raise the data part of the objective, times W, by more than logit(pi_e) nats (its closure
evidence E_e, in nats of log-likelihood over the training labels). With pi_e the prior that e is real, the posterior log-odds that e is false,
given the other edges as fitted and the encoder's parameters, is E_e - logit(pi_e), so P(false) = sigmoid(E_e -
logit(pi_e)). The digest records it per closed edge; figure 15 reads its calibration.

Usage:
  PYTHONPATH=. python experiments/summarise_edge_gate_simulations.py
"""
from __future__ import annotations

import argparse
import json
import math
from datetime import date
from pathlib import Path

import numpy as np

RUN_DIRECTORY = Path("runs/simulations")
OUTPUT_PATH = Path("docs/edge_gate_simulation_results.json")
CONFIRMATORY_GRAPH_DIRECTORY = Path("data/processed/graph_full_neuronal_split_binders_drugs_v5")
EVIDENCE_DIRECTORY = Path("data/processed/evidence_full_v5")
EDGE_PRIOR_TIERS = Path("configs/edge_prior_tiers.csv")
LABEL_GRADES = ("A", "B")

# name -> (role, what the run changes from the calibrated oracle fit of 800 perturbations)
RUNS = {
    **{f"b3_s{seed}": ("calibrated", f"world {seed}: calibrated priors, the encoder's true parameters, 800 perturbations") for seed in range(5)},
    **{f"b3_s{seed}_nobinary": ("no_binary_pass", f"world {seed} without the binary pass and the flip checks") for seed in range(5)},
    **{f"b3_s{seed}_curated_fixed": ("curated_edges_as_labels", f"world {seed} with the curated tier taken as known: no gate, always open (decision 9.63)")
       for seed in range(5)},
    "b3_s0_p3000": ("more_data", "world 0 with 3,000 perturbations"),
    "b3_s0_optimistic_priors": ("optimistic_priors", "world 0 with fewer real edges than the priors say (curated 0.8, unreferenced 0.45, text-mined 0.3)"),
    "b3_s0_pessimistic_priors": ("pessimistic_priors", "world 0 with more real edges than the priors say (curated 0.97, unreferenced 0.85, text-mined 0.75)"),
    "b3_s0_pu": ("half_the_positives", "world 0 with each true positive observed with probability 1/2"),
    "b3_s0_joint_random_support": ("joint_random_support", "world 0, the encoder and readout learned from a fresh start, the readout on each symptom's true nodes"),
    "b2_s0_joint_random_dense": ("joint_random_dense", "world 0, the encoder and readout learned from a fresh start, the readout on every node"),
    "b2_s0_joint_oracle_dense": ("joint_oracle_dense", "world 0, the encoder and readout learned from the truth's values, the readout on every node"),
    "b2_realratio": ("v5_labels_per_edge", "a world of 2,000 nodes and 20 symptoms at a 5 percent positive rate: about v5's labelled positives per gated edge"),
    "b2_s0_nochecks": ("no_flip_checks", "world 0 with the binary pass but without the flip checks (batch 2)"),
    **{f"b2_s0_width{label}": ("surrogate_width", f"world 0 with the binary pass's surrogate width {width} (batch 2)")
       for label, width in (("025", "0.25"), ("05", "0.5"), ("1", "1"))},
}
KEPT_KEYS = ("world", "fit", "fit_by_tier", "held_out_loss", "closures_by_evidence", "flip_checks", "stability", "seconds")
# the settings kept per run: enough to rebuild its command line (tests/test_edge_gate_simulation.py checks the runner against them) and to
# compare the trainer's defaults with what was tested; a run made before a setting existed lacks it
SETTING_KEYS = ("world_seed", "nodes", "mean_in_degree", "perturbations", "symptoms", "nodes_per_symptom", "held_out_share", "label_bias", "truth_real_probability",
                "label_sensitivity", "ungated_tiers", "joint", "joint_start", "joint_readout", "beta", "batch_size", "gate_learning_rate", "learning_rate",
                "stochastic_epochs", "binary_epochs", "binary_learning_rate", "binary_full_batch", "surrogate_width", "flip_checks", "closure_evidence",
                "estimator_diagnostic", "stability_refits", "stability_scheme", "weight_decay", "patience", "joint_validation_share")


def logit(probability: float) -> float:
    return math.log(probability / (1.0 - probability))


def posterior_false(evidence_nats: float, prior_real: float) -> float:
    """P(the edge is false | the data) from its closure evidence and its prior, as the module docstring derives."""
    return 1.0 / (1.0 + math.exp(-(evidence_nats - logit(prior_real))))


def selected_scores(edges: list[dict], keep) -> dict:
    chosen = [edge for edge in edges if keep(edge)]
    false = sum(not edge["real"] for edge in chosen)
    return {"closed": len(chosen), "false": false, "precision": false / len(chosen) if chosen else None}


def digest_run(name: str, run: dict) -> dict:
    entry = {key: run.get(key) for key in KEPT_KEYS}
    entry["settings"] = {key: run["settings"].get(key) for key in SETTING_KEYS if key in run["settings"]}
    if entry["flip_checks"]:
        entry["flip_checks"] = {key: value for key, value in entry["flip_checks"].items() if key != "accepted_gates"}
    closed = run.get("closed_edges")
    if closed is not None:
        for edge in closed:
            edge["posterior_false"] = None if edge["evidence_nats"] is None else posterior_false(edge["evidence_nats"], edge["prior"])
        entry["closed_edges"] = [[edge["tier"], round(edge["prior"], 6), int(edge["real"]), None if edge["evidence_nats"] is None else round(edge["evidence_nats"], 5),
                                  None if edge["posterior_false"] is None else round(edge["posterior_false"], 5)] for edge in closed]
        entry["closed_edges_columns"] = ["tier", "prior_real", "real", "evidence_nats", "posterior_false"]
        entry["closures_by_posterior_false"] = {f"at least {level}": selected_scores(closed, lambda edge, level=level: edge["posterior_false"] is not None
                                                                                        and edge["posterior_false"] >= level)
                                                for level in (0.5, 0.7, 0.8, 0.9, 0.95)}
    frequencies = run.get("stability_frequencies")
    if frequencies is not None:
        entry["stability_frequencies"] = [[item["tier"], int(item["real"]), round(item["frequency"], 4)] for item in frequencies]
        entry["stability_frequencies_columns"] = ["tier", "real", "frequency"]
    diagnostic = run.get("estimator_diagnostic")
    if diagnostic:
        entry["estimator"] = {key: value for key, value in diagnostic.items() if key != "per_edge_data_part_nats"}
        entry["estimator"]["per_edge_data_part_nats"] = [[round(estimate, 5), round(exact, 5)] for _, estimate, exact in diagnostic.get("per_edge_data_part_nats", [])]
        entry["estimator"]["columns"] = ["estimate_nats", "exact_nats"]
    return entry


def v5_label_density() -> dict | None:
    """Labelled cells and positives against gated edges on the v5 tables, or None when they are not on this machine."""
    edge_evidence = CONFIRMATORY_GRAPH_DIRECTORY / "edge_evidence.parquet"
    records_path = EVIDENCE_DIRECTORY / "evidence_records.parquet"
    if not edge_evidence.exists() or not records_path.exists():
        return None
    import pandas as pd

    from mechanistic_pathway_learning.graph.edge_evidence import edge_priors, read_edge_tiers

    evidence = pd.read_parquet(edge_evidence)
    gated, _, _ = edge_priors(evidence, read_edge_tiers(EDGE_PRIOR_TIERS))
    records = pd.read_parquet(records_path)
    labels = records[records.grade.isin(LABEL_GRADES) & (records.relation == "induces")]
    positives = int(labels[["perturbation_id", "symptom"]].drop_duplicates().shape[0])
    perturbations, symptoms = int(records.perturbation_id.nunique()), int(records.symptom.nunique())
    gated_edges = int(gated.sum())
    return {"graph": CONFIRMATORY_GRAPH_DIRECTORY.name, "edges": int(len(evidence)), "gated_edges": gated_edges, "perturbations": perturbations, "symptoms": symptoms,
            "label_grade_positive_pairs": positives, "cells_per_gated_edge": perturbations * symptoms / gated_edges,
            "positives_per_gated_edge": positives / gated_edges,
            "note": ("every perturbation and every symptom of the label table, so an upper bound on what one training fold sees; the simulated worlds count "
                     "their training rows only")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs", type=Path, default=RUN_DIRECTORY)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    arguments = parser.parse_args()
    runs, missing = {}, []
    for name, (role, description) in RUNS.items():
        path = arguments.runs / f"{name}.json"
        if not path.exists():
            missing.append(name)
            continue
        runs[name] = {"role": role, "description": description, **digest_run(name, json.loads(path.read_text()))}
    digest = {"generated_by": "experiments/summarise_edge_gate_simulations.py", "generated_on": date.today().isoformat(),
              "written_by": "experiments/simulate_edge_gate_recovery.py through scripts/run_edge_gate_simulations.sh", "runs": runs, "missing": missing,
              "v5": v5_label_density()}
    arguments.output.write_text(json.dumps(digest, indent=1, default=float) + "\n", encoding="utf-8")
    print(f"wrote {arguments.output}: {len(runs)} runs, {len(missing)} missing")


if __name__ == "__main__":
    main()
