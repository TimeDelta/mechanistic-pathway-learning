"""Run experiments/run_main_model.py over a named configuration for several splits, sequentially and idempotently.

Each (configuration, split, seed) is one run directory under runs/<configuration>/; a split whose
DONE marker exists is skipped, so a batch can be resubmitted after an interruption. Configurations
are named sets of extra arguments (below) and can be extended on the command line with --extra.

Usage:
  python experiments/run_main_model_batch.py --configuration b6_default --folds 0 1 2 3 4
  python experiments/run_main_model_batch.py --configuration b6_default --holdout-subsystems "Tryptophan metabolism" ...
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

FULL_GRAPH_ARGUMENTS = ["--graph-dir", "data/processed/graph_full_neuronal", "--evidence-dir", "data/processed/evidence_full_v2",
                        "--label-selection", "data/processed/label_selection/better_v1_full_v2.parquet"]
FULL_GRAPH_NODE_PROPERTIES = "data/processed/node_descriptors/full_neuronal_descriptors_brain_expression.parquet"
FULL_GRAPH_CELL_CLASS_WEIGHTS = "data/processed/cell_class_weights/full_neuronal_cell_class_weights.parquet"
NOISY_OR_SETTINGS = ["--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02",
                     "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"]
MESSAGE_PASSING_ENCODER = ["--node-features", "typed"]
LINEAR_RESPONSE_ENCODER = ["--encoder", "linear_response", "--cofactor-relations", "--cell-class-weights", FULL_GRAPH_CELL_CLASS_WEIGHTS, "--extracellular-coupling"]

CONFIGURATIONS: dict[str, list[str]] = {
    "b6_default": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum"],
    "b6_absolute_mean": ["--head", "noisy_or", "--field", "absolute", "--pooling", "mean"],
    "b6_k1": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--num-modules", "1"],
    "b6_k16": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--num-modules", "16"],
    "b6_no_description_length": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--description-length-coefficient", "0"],
    "b6_leak_init": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--init-leak-from-base-rate"],
    "b6_two_entries_grade_policy": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--evidence-dir", "data/processed/evidence_two_entries"],
    "b6_mechanistic": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5"],
    "b6_mechanistic_expected_gate": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5"],
    "b6_mechanistic_time_scales": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5",
                                   "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002"],
    "b6_off_by_default": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5"],
    "b6_typed_nodes": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate"],
    "b3_linear_response": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response"],
    "b6_linear_response": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5"],
    "b6_linear_response_time_scales": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5",
                                       "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002"],
    "b3_linear_response_cofactors": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations"],
    # task 25: one global spectral scale and no per-node divisor, so stoichiometric counts survive; the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_spectral": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                              "--normalisation", "spectral"],
    # the counts the spectral arm keeps, without its global shrink: one divisor per node over every relation; the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_total_in_degree": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                     "--cofactor-relations", "--normalisation", "total_in_degree"],
    # the sign ablation: every edge +1 and every gain positive, same graph, split and head; the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_unsigned": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                              "--cofactor-relations", "--edge-signs", "all_positive"],
    # the sign permutation: the signs shuffled among edges, so the count of negative edges stays and their placement is lost;
    # the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_permuted_signs": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                    "--cofactor-relations", "--edge-signs", "permuted"],
    # the relation-typing ablation: one gain shared by every relation; the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_shared_gain": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                 "--cofactor-relations", "--relation-gains", "shared"],
    # task 24: a sparsemax-weighted mixture of odd statistics across relations in place of their mean; the only change from b3_linear_response_cofactors
    # the node descriptors (protein, metabolite and reaction properties) in the linear-response output gate, which until
    # now saw only the structural features; the only change from b3_linear_response_cofactors
    "b3_linear_response_cofactors_descriptors": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                 "--cofactor-relations", "--node-descriptors", "data/processed/graph/node_descriptors.parquet"],
    # brain region and brain cell-class expression (experiments/build_brain_expression_descriptors.py) appended to the
    # node descriptors, for genes and, through the gene rule, for reactions; the only change from
    # b3_linear_response_cofactors_descriptors
    "b3_linear_response_cofactors_descriptors_brain": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                       "--cofactor-relations", "--node-descriptors",
                                                       "data/processed/node_descriptors/slice_descriptors_brain_expression.parquet"],
    "b3_linear_response_cofactors_mixture": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                             "--cross-relation-aggregator", "softmax_mixture", "--mixture-weighting", "sparsemax"],
    "b6_linear_response_time_scales_cofactors": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                                 "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002"],
    "b3_linear_response_cofactors_log": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations", "--response-scale", "signed_log"],
    "b6_linear_response_time_scales_cofactors_log": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations", "--response-scale", "signed_log",
                                                     "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002"],
    # gate time scale: unsupported gates close within the 200 to 400 steps early stopping allows (docs/b6_module_diagnosis.md)
    "b6_mechanistic_gate_time_scales": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5",
                                        "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"],
    "b6_linear_response_gate_time_scales_cofactors": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                                      "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"],
    "b6_linear_response_gate_time_scales_cofactors_log": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations", "--response-scale", "signed_log",
                                                          "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"],
    # fixed node descriptors (experiments/build_node_descriptors.py) and the controls that credit propagation only with what it adds
    "b3_typed_nodes_descriptors": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--node-features", "typed",
                                   "--node-descriptors", "data/processed/graph/node_descriptors.parquet"],
    "b3_descriptors_only": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "local_descriptors",
                            "--node-descriptors", "data/processed/graph/node_descriptors.parquet"],
    "b3_local_structural": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "local_descriptors"],
    # auxiliary supervision of the field with measured metabolite directions (models/laboratory_readout.py)
    "b3_typed_nodes_laboratory": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--laboratory-label-weight", "1.0"],
    "b3_linear_response_cofactors_laboratory": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                                "--laboratory-label-weight", "1.0"],
    # manganese graph variant (graph/manganese_extension.py): ablation against the same configuration on the base graph
    "b3_linear_response_cofactors_laboratory_manganese": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                                          "--laboratory-label-weight", "1.0", "--graph-dir", "data/processed/graph_manganese"],
    # neuronal and oxidative graph variant (graph/reactome_import.py, neurotransmission_regulation.py, oxidative_regulation.py): the
    # ablation of design section 5.7 against the same configuration on the metabolic graph, with and without the curated layers
    "b3_linear_response_cofactors_laboratory_neuronal": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                                         "--laboratory-label-weight", "1.0", "--graph-dir", "data/processed/graph_neuronal"],
    "b6_linear_response_cofactors_neuronal": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                              "--graph-dir", "data/processed/graph_neuronal", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5",
                                              "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"],
    "b3_linear_response_cofactors_laboratory_neuronal_no_curation": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response",
                                                                     "--cofactor-relations", "--laboratory-label-weight", "1.0",
                                                                     "--graph-dir", "data/processed/graph_neuronal_no_curation"],
    # the full-graph runs (amendment of 7 October in docs/preregistration.md): graph_full_neuronal with the full evidence
    # and the better_v1 label selection, fixed before any of them was scored. full_linear_response_properties and
    # full_linear_response differ in one argument, the node properties (protein, metabolite and reaction descriptors and
    # brain region and cell-class expression); full_b3_typed_nodes_properties is the relational GNN comparator the primary
    # endpoint names, given the same node properties.
    "full_linear_response_properties": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                        *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES],
    "full_linear_response": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                             *FULL_GRAPH_ARGUMENTS],
    "full_b6_linear_response_properties": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations",
                                           *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES,
                                           "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02",
                                           "--leak-learning-rate", "0.0002", "--gate-learning-rate", "0.05"],
    "full_b3_typed_nodes_properties": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--node-features", "typed",
                                       *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES],
    # the confirmatory family (docs/preregistration.md, specification of 8 October 2026): two encoders by two heads, every one
    # on graph_full_neuronal with the node descriptors (protein, metabolite and reaction properties, brain region and
    # cell-class expression) and the better_v1 selection of the expanded labels; the linear-response encoder also
    # propagates once per cell class (the HPA classes, the dopaminergic class and all cells) with the extracellular
    # metabolites shared between classes. Run with --group-by disease_cluster_and_targets --lockbox configs/lockbox_v1.json.
    "confirmatory_message_passing_noisy_or": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", *MESSAGE_PASSING_ENCODER,
                                              *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES, *NOISY_OR_SETTINGS],
    "confirmatory_message_passing_sigmoid": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", *MESSAGE_PASSING_ENCODER,
                                             *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES],
    "confirmatory_linear_response_noisy_or": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", *LINEAR_RESPONSE_ENCODER,
                                              *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES, *NOISY_OR_SETTINGS],
    "confirmatory_linear_response_sigmoid": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", *LINEAR_RESPONSE_ENCODER,
                                             *FULL_GRAPH_ARGUMENTS, "--node-descriptors", FULL_GRAPH_NODE_PROPERTIES],
    "b3_typed_nodes_degree": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--degree-offset"],
    "b6_mechanistic_degree": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--node-features", "typed", "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--degree-offset"],
    "b3_degree_only": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "none", "--degree-offset"],
    "b3_linear_response_degree": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--degree-offset"],
    "b6_linear_response_time_scales_cofactors_degree": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--encoder", "linear_response", "--cofactor-relations", "--degree-offset",
                                                        "--init-leak-from-base-rate", "--module-bias-init", "-3", "--gate-init-noise", "0.5", "--link-learning-rate", "0.02", "--leak-learning-rate", "0.0002"],
    "b3_typed_nodes": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--node-features", "typed"],
    "b6_frequency_target": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--positive-target-from-frequency"],
    "b6_default_permuted": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--permute-labels"],
    "b3_sigmoid": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum"],
    "b3_sigmoid_absolute": ["--head", "sigmoid", "--field", "absolute", "--pooling", "mean"],
    "b3_sigmoid_permuted": ["--head", "sigmoid", "--field", "difference", "--pooling", "sum", "--permute-labels"],
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--configuration", choices=sorted(CONFIGURATIONS), required=True)
    parser.add_argument("--folds", type=int, nargs="*", default=[0, 1, 2, 3, 4])
    parser.add_argument("--holdout-modules", nargs="*", default=[])
    parser.add_argument("--holdout-subsystems", nargs="*", default=[])
    parser.add_argument("--seeds", type=int, nargs="*", default=[0])
    parser.add_argument("--group-by", choices=["gene", "disease_cluster", "disease_cluster_and_targets"], default="disease_cluster")
    parser.add_argument("--run-root", type=Path, default=Path("runs"))
    parser.add_argument("--lockbox", type=Path, default=None,
                        help="lockbox file (experiments/draw_lockbox.py): pilots run on the development set into <configuration>_<group-by>_development")
    parser.add_argument("--score-lockbox", action="store_true",
                        help="with --lockbox: one run per seed trained on the development set and scored on the lockbox, into <configuration>_<group-by>_confirmatory")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="further arguments passed through to run_main_model.py")
    arguments = parser.parse_args()
    if arguments.score_lockbox and arguments.lockbox is None:
        raise SystemExit("--score-lockbox needs --lockbox")
    suffix = "" if arguments.lockbox is None else ("_confirmatory" if arguments.score_lockbox else "_development")  # fold numbers differ once the lockbox is removed
    run_directory = arguments.run_root / f"{arguments.configuration}_{arguments.group_by}{suffix}"
    lockbox_arguments = [] if arguments.lockbox is None else ["--lockbox", str(arguments.lockbox)] + (["--score-lockbox"] if arguments.score_lockbox else [])
    jobs: list[list[str]] = []
    for seed in arguments.seeds:
        if arguments.score_lockbox:
            jobs.append(["--seed", str(seed)])
        elif arguments.holdout_modules or arguments.holdout_subsystems:
            jobs += [["--holdout-module", module, "--seed", str(seed)] for module in arguments.holdout_modules]
            jobs += [["--holdout-subsystem", subsystem, "--seed", str(seed)] for subsystem in arguments.holdout_subsystems]
        else:
            jobs += [["--fold", str(fold), "--seed", str(seed)] for fold in arguments.folds]
    for job in jobs:
        command = [sys.executable, "experiments/run_main_model.py", "--run-dir", str(run_directory), "--group-by", arguments.group_by, "--resume",
                   *CONFIGURATIONS[arguments.configuration], *lockbox_arguments, *job, *arguments.extra]
        print("running:", " ".join(command), flush=True)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            print(f"run failed with code {completed.returncode}; continuing with the next split", flush=True)


if __name__ == "__main__":
    main()
