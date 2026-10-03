#!/bin/bash
#SBATCH --job-name=mpl-train
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --signal=B:USR1@600      # ten minutes before the wall-time limit, ask the trainer to checkpoint
#SBATCH --output=slurm-%j.out
# Resumable training for a runtime-limited cluster: the trainer writes a checkpoint on
# USR1 and on its own timer; this script requeues itself until a DONE marker appears.
set -euo pipefail
source "${MPL_ENV_ACTIVATE:-$HOME/.venvs/mpl/bin/activate}"
cd "${MPL_REPO_DIR:-$HOME/mechanistic-pathway-learning}"
RUN_DIR="runs/${MPL_RUN_NAME:?set MPL_RUN_NAME}"
mkdir -p "${RUN_DIR}"
trap 'kill -USR1 "${TRAIN_PID}"; wait "${TRAIN_PID}"' USR1
python -m experiments.run_main_model --config "${MPL_CONFIG:-configs/model_noisy_or_pathway_modules.yaml}" \
  --run-dir "${RUN_DIR}" --resume &
TRAIN_PID=$!
wait "${TRAIN_PID}" || true
if [[ ! -f "${RUN_DIR}/DONE" ]]; then
  sbatch --dependency=afterany:"${SLURM_JOB_ID}" "$0"
fi
