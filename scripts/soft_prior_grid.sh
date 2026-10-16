#!/bin/bash
# The soft-prior grid on the development folds (docs/soft_priors.md, section 4): the global strength and the
# tempering exponent are chosen on development data before any registration, and the permuted control is run at
# the chosen point afterwards. This script only loops the batch runner over the grid; it does not choose.
#
# The data flags are the job's, as for every confirmatory run (docs/confirmatory_runbook.md: the job scripts under
# runs/full/ carry --graph-dir, --evidence-dir, --label-selection, --node-descriptors, --cell-class-weights and
# --combined-perturbations). Put them in MPL_DATA_ARGS; the report tables that become priors go in
# MPL_SOFT_PRIOR_REPORTS (defaults below to the rebuild's outputs and the release's literature table).
#
#   MPL_DATA_ARGS="--graph-dir data/processed/graph_full_neuronal_split_binders_drugs_v4 --evidence-dir data/processed/evidence_full_v4 ..." \
#   bash scripts/soft_prior_grid.sh linear_response            # or message_passing
#   bash scripts/soft_prior_grid.sh linear_response --permuted  # the control at one point: set MPL_GRID_STRENGTHS and MPL_GRID_EXPONENTS to it
#   MPL_DRY_RUN=1 bash scripts/soft_prior_grid.sh linear_response  # print the commands and run nothing
#
# One run directory per grid point: runs/soft_prior_grid/<encoder>/s<strength>_t<exponent>/v4_<encoder>_soft_priors[_permuted]_<group-by>_development.
set -euo pipefail
ENCODER="${1:?linear_response or message_passing}"
CONTROL=""
if [[ "${2:-}" == "--permuted" ]]; then CONTROL="_permuted"; fi
STRENGTHS="${MPL_GRID_STRENGTHS:-0.25 0.5 1.0}"
EXPONENTS="${MPL_GRID_EXPONENTS:-0 0.5 1.0}"
REPORTS="${MPL_SOFT_PRIOR_REPORTS:-data/processed/clinical_trials/clinical_trial_reports.parquet data/processed/faers/faers_reports.parquet data/processed/impc/impc_reports.parquet data/releases/v0.4/literature/literature_reports.parquet}"
LOCKBOX="${MPL_LOCKBOX:-configs/lockbox_v2.json}"
: "${MPL_DATA_ARGS:?set MPL_DATA_ARGS to the confirmatory data flags (docs/confirmatory_runbook.md)}"
RUN=""
if [[ -n "${MPL_DRY_RUN:-}" ]]; then RUN="echo"; fi
for strength in ${STRENGTHS}; do
  for exponent in ${EXPONENTS}; do
    run_root="runs/soft_prior_grid/${ENCODER}/s${strength}_t${exponent}"
    if [[ -z "${RUN}" ]]; then mkdir -p "${run_root}"; fi
    # shellcheck disable=SC2086
    ${RUN} python experiments/run_main_model_batch.py --configuration "v4_${ENCODER}_soft_priors${CONTROL}" --folds 0 1 2 3 4 \
      --group-by disease_cluster_and_targets --lockbox "${LOCKBOX}" --run-root "${run_root}" \
      --extra ${MPL_DATA_ARGS} --soft-prior-strength "${strength}" --soft-prior-tempering "${exponent}" --soft-prior-reports ${REPORTS}
  done
done
