#!/bin/bash
# The edge-gate simulations of docs/edge_gate_simulation.md (decision 9.62), each written to runs/simulations/<name>.json,
# then the digest docs/edge_gate_simulation_results.json that the record and figures 15 and 16 read.
# Seconds to minutes per run on two cores; about 40 minutes in all.
#   bash scripts/run_edge_gate_simulations.sh            # every run
#   bash scripts/run_edge_gate_simulations.sh b3_s0      # the named runs only
#   bash scripts/run_edge_gate_simulations.sh --list     # each run's name and arguments, one per line, tab-separated
set -euo pipefail
cd "$(dirname "$0")/.."
export PYTHONPATH=.
mkdir -p runs/simulations
detail="--closure-evidence --per-edge-detail"
binary="--binary-full-batch --flip-checks 30"
declare -A RUN
for seed in 0 1 2 3 4; do
  RUN[b3_s${seed}]="--world-seed $seed $binary --estimator-diagnostic 300 --stability-refits 20 --stability-scheme complementary $detail"
  RUN[b3_s${seed}_nobinary]="--world-seed $seed --binary-epochs 0 --flip-checks 0 $detail"
done
RUN[b3_s0_p3000]="--world-seed 0 --perturbations 3000 $binary $detail"
RUN[b3_s0_pessimistic_priors]="--world-seed 0 $binary --truth-real-probability curated=0.97 unreferenced=0.85 text_mined_only=0.75 $detail"
RUN[b3_s0_optimistic_priors]="--world-seed 0 $binary --truth-real-probability curated=0.8 unreferenced=0.45 text_mined_only=0.3 $detail"
RUN[b3_s0_pu]="--world-seed 0 $binary --label-sensitivity 0.5 $detail"
RUN[b3_s0_joint_random_support]="--world-seed 0 $binary --joint --joint-start random --joint-readout support --stochastic-epochs 60 $detail"
RUN[b2_s0_joint_random_dense]="--world-seed 0 $binary --joint --joint-start random --stochastic-epochs 60 --closure-evidence"
RUN[b2_s0_joint_oracle_dense]="--world-seed 0 $binary --joint --stochastic-epochs 60 --closure-evidence"
RUN[b2_realratio]="--world-seed 0 $binary --nodes 2000 --symptoms 20 --label-bias -3 --stability-refits 10 --stability-scheme complementary --closure-evidence"
for seed in 0 1 2 3 4; do
  RUN[b3_s${seed}_curated_fixed]="--world-seed $seed $binary --ungated-tiers curated $detail"
done
# batch 4 (17 October 2026): does the binary pass settle with more epochs, or with its learning rate annealed?
for seed in 0 1 2 3 4; do
  RUN[b4_s${seed}_epochs32]="--world-seed $seed $binary --binary-epochs 32 $detail"
  RUN[b4_s${seed}_epochs32_decay08]="--world-seed $seed $binary --binary-epochs 32 --binary-learning-rate-decay 0.8 $detail"
done
RUN[b2_s0_nochecks]="--world-seed 0 --binary-full-batch --flip-checks 0 --closure-evidence"
for width in 0.25 0.5 1; do
  RUN[b2_s0_width${width/./}]="--world-seed 0 $binary --closure-evidence --surrogate-width $width"
done
if [ "${1:-}" = "--list" ]; then
  for name in "${!RUN[@]}"; do printf '%s\t%s\n' "$name" "${RUN[$name]}"; done | sort
  exit 0
fi
names=("$@")
if [ ${#names[@]} -eq 0 ]; then names=($(printf '%s\n' "${!RUN[@]}" | sort)); fi
for name in "${names[@]}"; do
  python experiments/simulate_edge_gate_recovery.py ${RUN[$name]} --output runs/simulations/$name.json > runs/simulations/$name.log 2>&1
  echo "done: $name"
done
python experiments/summarise_edge_gate_simulations.py
