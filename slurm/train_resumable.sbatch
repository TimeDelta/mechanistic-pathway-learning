#!/bin/bash
#SBATCH --job-name=mpl-train
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --signal=B:USR1@600      # ten minutes before the wall-time limit, ask the trainer to checkpoint
#SBATCH --output=slurm-%j.out
# Resumable training for a runtime-limited cluster: the trainer writes a checkpoint on
# USR1 and on its own timer; this script requeues itself until the split's DONE marker appears.
#
# One split per job. Set MPL_RUN_NAME (run directory under runs/) and MPL_TRAIN_ARGS, for example
#   MPL_RUN_NAME=b6_default_disease_cluster MPL_TRAIN_ARGS="--head noisy_or --field difference --pooling sum --group-by disease_cluster --fold 0 --seed 0"
# or submit an array over folds with --array=0-4 and leave --fold out of MPL_TRAIN_ARGS: the task id is the fold.
# experiments/run_main_model_batch.py lists the named configurations and their arguments.
set -euo pipefail
source "${MPL_ENV_ACTIVATE:-$HOME/.venvs/mpl/bin/activate}"
cd "${MPL_REPO_DIR:-$HOME/mechanistic-pathway-learning}"
RUN_DIR="runs/${MPL_RUN_NAME:?set MPL_RUN_NAME}"
mkdir -p "${RUN_DIR}"
FOLD_ARGS=""
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  FOLD_ARGS="--fold ${SLURM_ARRAY_TASK_ID}"
fi
# shellcheck disable=SC2086
python experiments/run_main_model.py --run-dir "${RUN_DIR}" --resume ${MPL_TRAIN_ARGS:-} ${FOLD_ARGS} &
TRAIN_PID=$!
trap 'kill -USR1 "${TRAIN_PID}"; wait "${TRAIN_PID}"' USR1
wait "${TRAIN_PID}" || true
# The trainer names the split directory from its arguments; requeue while any split in the run lacks a DONE marker
# only when this job's own split is unfinished (the newest split directory without DONE is this job's).
if ! ls "${RUN_DIR}"/*/DONE >/dev/null 2>&1 || [[ -n "$(find "${RUN_DIR}" -mindepth 1 -maxdepth 1 -type d -newer "${RUN_DIR}" ! -exec test -e '{}/DONE' \; -print | head -1)" ]]; then
  sbatch --dependency=afterany:"${SLURM_JOB_ID}" "$0"
fi
