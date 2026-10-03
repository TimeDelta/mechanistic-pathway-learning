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

CONFIGURATIONS: dict[str, list[str]] = {
    "b6_default": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum"],
    "b6_absolute_mean": ["--head", "noisy_or", "--field", "absolute", "--pooling", "mean"],
    "b6_k1": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--num-modules", "1"],
    "b6_k16": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--num-modules", "16"],
    "b6_no_description_length": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--description-length-coefficient", "0"],
    "b6_leak_init": ["--head", "noisy_or", "--field", "difference", "--pooling", "sum", "--init-leak-from-base-rate"],
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
    parser.add_argument("--group-by", choices=["gene", "disease_cluster"], default="disease_cluster")
    parser.add_argument("--run-root", type=Path, default=Path("runs"))
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="further arguments passed through to run_main_model.py")
    arguments = parser.parse_args()
    run_directory = arguments.run_root / f"{arguments.configuration}_{arguments.group_by}"
    jobs: list[list[str]] = []
    for seed in arguments.seeds:
        if arguments.holdout_modules or arguments.holdout_subsystems:
            jobs += [["--holdout-module", module, "--seed", str(seed)] for module in arguments.holdout_modules]
            jobs += [["--holdout-subsystem", subsystem, "--seed", str(seed)] for subsystem in arguments.holdout_subsystems]
        else:
            jobs += [["--fold", str(fold), "--seed", str(seed)] for fold in arguments.folds]
    for job in jobs:
        command = [sys.executable, "experiments/run_main_model.py", "--run-dir", str(run_directory), "--group-by", arguments.group_by, "--resume",
                   *CONFIGURATIONS[arguments.configuration], *job, *arguments.extra]
        print("running:", " ".join(command), flush=True)
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            print(f"run failed with code {completed.returncode}; continuing with the next split", flush=True)


if __name__ == "__main__":
    main()
