#!/bin/bash
#SBATCH --job-name=mpl-flux
#SBATCH --array=0-3000%40          # set the upper bound to (number of genes - 1) after gene_order.json exists
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00           # one gene per task fits inside a short wall-time window
#SBATCH --output=slurm-%A_%a.out
# Idempotent: a task whose output parquet exists returns immediately, so a partially
# finished array can be resubmitted as-is. Task 0 must finish first because it writes
# the baseline samples; submit it alone, then the rest with --dependency=afterok:<jobid>.
set -euo pipefail
source "${MPL_ENV_ACTIVATE:-$HOME/.venvs/mpl/bin/activate}"
cd "${MPL_REPO_DIR:-$HOME/mechanistic-pathway-learning}"
python -m mechanistic_pathway_learning.perturbation.flux_sampling_perturbation_features \
  --model-path data/raw/Human-GEM/model/Human-GEM.xml \
  --output-dir data/cache/flux_perturbation \
  --gene-index "${SLURM_ARRAY_TASK_ID}" \
  --num-samples 2000 --processes "${SLURM_CPUS_PER_TASK}"
