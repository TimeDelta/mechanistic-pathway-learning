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
#   bash scripts/resume_jobs.sh --stop NAME     stop job NAME and every process it started, launch nothing
# runs/ is gitignored: the registry and the markers live on the session's disk, which survives a container restart.
set -uo pipefail
repository_root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$repository_root"
job_directory="runs/jobs"
mkdir -p "$job_directory"

mode="resume"
stop_job_name=""
if [ "${1:-}" = "--status" ]; then
  mode="status"
elif [ "${1:-}" = "--stop" ]; then
  if [ "$#" -lt 2 ]; then
    echo "usage: $0 --stop NAME" >&2
    exit 2
  fi
  mode="stop"
  stop_job_name="$2"
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

# Stop a job and every process it started. setsid makes the job command a process group leader and its children
# inherit that group, so killing the group leaves nothing behind; killing only the pid that matches JOB_COMMAND leaves
# the trainer it spawned alive with no parent, and the next relaunch then has two processes writing one fold directory.
stop_job_process_group() {
  local job_name="$1" job_command="$2" leader_pid process_group surviving_pids
  leader_pid="$(pgrep -f -x -- "$job_command" | head -1)"
  if [ -z "$leader_pid" ]; then
    echo "job $job_name: not running, nothing to stop"
    return 0
  fi
  process_group="$(ps -o pgid= -p "$leader_pid" | tr -d ' ')"
  if [ "$process_group" != "$leader_pid" ]; then
    echo "job $job_name: pid $leader_pid is not its own process group leader (group $process_group), not stopped" >&2
    return 1
  fi
  echo "job $job_name: stopping process group $process_group (pids $(pgrep -g "$process_group" | tr '\n' ' '))"
  kill -TERM -- "-$process_group" 2>/dev/null
  for _ in 1 2 3 4 5 6 7 8 9 10; do
    surviving_pids="$(pgrep -g "$process_group" | tr '\n' ' ')"
    [ -z "$surviving_pids" ] && break
    sleep 1
  done
  if [ -n "$surviving_pids" ]; then
    echo "job $job_name: group $process_group ignored SIGTERM (pids $surviving_pids), sending SIGKILL"
    kill -KILL -- "-$process_group" 2>/dev/null
  fi
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) stopped process group $process_group" >> "$job_directory/$job_name.history"
}

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
  if [ "$mode" = "stop" ]; then
    [ "$job_name" = "$stop_job_name" ] && stop_job_process_group "$job_name" "$JOB_COMMAND"
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
