#!/bin/bash
# Relaunch unfinished, resumable background jobs after the cloud container was reclaimed (for example when a usage
# limit stopped the session). Every job is a file runs/jobs/<name>.job that sets:
#   JOB_COMMAND        the command, run from the repository root (for example "bash runs/literature_rerun.sh");
#                      it must be resumable: rerunning it continues from checkpoints and caches
#   JOB_DONE_MARKER    a path that exists once the job finished
#   JOB_FAILED_MARKER  optional; a path that exists when the job failed and needs attention (never relaunched)
#   JOB_DESCRIPTION    optional; one line for the status report
# A job is relaunched, detached with setsid nohup, when neither marker exists and no process runs its exact command.
#
# Usage:
#   bash scripts/resume_jobs.sh                 relaunch what is unfinished and report every job
#   bash scripts/resume_jobs.sh --status        report only, launch nothing
#   bash scripts/resume_jobs.sh --register NAME "COMMAND" DONE_MARKER [FAILED_MARKER] [DESCRIPTION]
#                                               write runs/jobs/NAME.job, then relaunch as above
# runs/ is gitignored: the registry and the markers live on the session's disk, which survives a container restart.
set -uo pipefail
repository_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repository_root"
job_directory="runs/jobs"
mkdir -p "$job_directory"

mode="resume"
if [ "${1:-}" = "--status" ]; then
  mode="status"
elif [ "${1:-}" = "--register" ]; then
  if [ "$#" -lt 4 ]; then
    echo "usage: $0 --register NAME \"COMMAND\" DONE_MARKER [FAILED_MARKER] [DESCRIPTION]" >&2
    exit 2
  fi
  {
    printf 'JOB_COMMAND=%q\n' "$3"
    printf 'JOB_DONE_MARKER=%q\n' "$4"
    printf 'JOB_FAILED_MARKER=%q\n' "${5:-}"
    printf 'JOB_DESCRIPTION=%q\n' "${6:-}"
  } > "$job_directory/$2.job"
  echo "registered job $2"
fi

job_files=("$job_directory"/*.job)
if [ ! -e "${job_files[0]}" ]; then
  echo "resume_jobs: no registered jobs"
  exit 0
fi

for job_file in "${job_files[@]}"; do
  job_name="$(basename "$job_file" .job)"
  JOB_COMMAND="" JOB_DONE_MARKER="" JOB_FAILED_MARKER="" JOB_DESCRIPTION=""
  # shellcheck disable=SC1090
  source "$job_file"
  label="job $job_name${JOB_DESCRIPTION:+ ($JOB_DESCRIPTION)}"
  if [ -z "$JOB_COMMAND" ] || [ -z "$JOB_DONE_MARKER" ]; then
    echo "$label: malformed job file, skipped"
    continue
  fi
  if [ -e "$JOB_DONE_MARKER" ]; then
    echo "$label: done"
  elif [ -n "$JOB_FAILED_MARKER" ] && [ -e "$JOB_FAILED_MARKER" ]; then
    echo "$label: failed ($(head -c 200 "$JOB_FAILED_MARKER" | tr '\n' ' ')), needs attention, not relaunched"
  elif pgrep -f -x -- "$JOB_COMMAND" > /dev/null; then
    echo "$label: running"
  elif [ "$mode" = "status" ]; then
    echo "$label: not running and unfinished"
  else
    setsid nohup bash -c "$JOB_COMMAND" > "$job_directory/$job_name.launch.log" 2>&1 < /dev/null &
    disown
    echo "$label: relaunched at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) relaunched" >> "$job_directory/$job_name.history"
  fi
done
