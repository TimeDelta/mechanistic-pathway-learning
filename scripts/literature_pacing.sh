#!/bin/bash
# Gate for the paced literature appraisal. Prints "proceed" at most once per
# MINIMUM_SECONDS_BETWEEN_APPRAISALS and at most APPRAISAL_CAP times in total;
# otherwise prints the reason to stop and exits non-zero. Call it with "claim"
# before reading a source, which both checks and consumes the slot, so a firing
# cannot read two sources by calling twice.
set -euo pipefail

MINIMUM_SECONDS_BETWEEN_APPRAISALS=300
APPRAISAL_CAP=19
state_directory="runs/literature_appraisal"
state_file="$state_directory/pacing_state"
mkdir -p "$state_directory"
[ -f "$state_file" ] || printf 'last_appraisal_epoch=0\nappraisals_consumed=0\n' > "$state_file"
# shellcheck disable=SC1090
. "$state_file"

now_epoch="$(date -u +%s)"
seconds_since_last=$(( now_epoch - last_appraisal_epoch ))

if [ "$appraisals_consumed" -ge "$APPRAISAL_CAP" ]; then
  echo "stop: the cap of $APPRAISAL_CAP appraisals is spent; report and end the chain"
  exit 1
fi
if [ "$seconds_since_last" -lt "$MINIMUM_SECONDS_BETWEEN_APPRAISALS" ]; then
  echo "stop: only ${seconds_since_last}s since the last appraisal, minimum is ${MINIMUM_SECONDS_BETWEEN_APPRAISALS}s; read no source this firing"
  exit 1
fi
if [ "${1:-check}" = "claim" ]; then
  printf 'last_appraisal_epoch=%s\nappraisals_consumed=%s\n' "$now_epoch" "$(( appraisals_consumed + 1 ))" > "$state_file"
  echo "proceed: slot $(( appraisals_consumed + 1 )) of $APPRAISAL_CAP claimed, ${seconds_since_last}s since the last"
else
  echo "proceed: slot $(( appraisals_consumed + 1 )) of $APPRAISAL_CAP available, ${seconds_since_last}s since the last"
fi
