# The background job harness and the check-in routines

Long jobs run detached in a cloud container that is reclaimed a few minutes after the session goes idle, and again
whenever a usage limit stops the session. Two things keep them going: a registry with a resume script, and scheduled
check-ins that run the resume script.

## The registry

`scripts/resume_jobs.sh` relaunches every unfinished job registered as `runs/jobs/<name>.job`. A job file sets
`JOB_COMMAND` (resumable: rerunning it continues from checkpoints), `JOB_DONE_MARKER`, an optional
`JOB_FAILED_MARKER` and an optional `JOB_DESCRIPTION`. `runs/` is not in git, so the registry and the markers live on
the session's disk, which survives a container restart but not a new container; `.claude/hooks/session-start.sh` runs
the script when a container starts.

- `bash scripts/resume_jobs.sh` relaunches what is unfinished and reports every job.
- `bash scripts/resume_jobs.sh --status` reports only.
- `bash scripts/resume_jobs.sh --register NAME "COMMAND" DONE_MARKER [FAILED_MARKER] [DESCRIPTION]` registers, then
  relaunches.
- `bash scripts/resume_jobs.sh --stop NAME` stops job NAME and the whole process group it leads. Never kill the pid
  that matches `JOB_COMMAND` instead: that orphans the trainer it spawned, and the next resume puts two processes in
  one fold directory. To keep a job from relaunching, move its `.job` file to `runs/jobs/paused/`.

Every run but `--stop` ends with the line `resume_jobs: N active jobs` and writes the same count to
`runs/jobs/active_job_count`. A job is active when it is running, or unfinished with neither marker. A failed job is
not active: it is never relaunched and needs a person.

## The check-in routines, and when they are paused

Four scheduled routines fire into the working session, staggered so a reclaimed container is restarted within about
fifteen minutes:

| id | schedule | name |
|---|---|---|
| `trig_01E9fxBn6XU4RMdEjPQEyndm` | `3 * * * *` | Job keep-alive :03 |
| `trig_017eJmHnGUq6q2CzEz3kt1Y5` | `19 * * * *` | Job keep-alive :19 |
| `trig_01LJeguqEjTcRN7y2VWKj2qA` | `35 * * * *` | Job keep-alive :35 |
| `trig_01LoV65i2vYRtdnCBpV3dY2D` | `51 * * * *` | Hourly job check-in: resume reclaimed jobs |

**The pause rule** (the user's instruction of 9 October 2026: "Update the harness so that when there are no running
jobs, the check-ins get paused and only restarted when there are running jobs please"). The routines exist only to
relaunch registered jobs, so:

- **Disabled while no job is active.** A firing runs `scripts/resume_jobs.sh` first. If it reports `0 active jobs`,
  that firing disables all four routines (`update_trigger`, `enabled=false`), says so in one line and ends the turn
  without arming a wait. Each routine's prompt carries this as its first step, so the rule runs whether or not anyone
  remembers it.
- **Re-enabled when a job is registered.** A job is only ever registered from inside a session, so re-enabling
  belongs to whoever registers it: after `--register`, enable all four routines (`update_trigger`, `enabled=true`).
  The script prints a reminder on every `--register` for that reason.
- Disabling is a pause, not a deletion. The ids above stay valid, the prompts and schedules are kept, and
  `ended_reason` stays empty. Deleting a routine would lose its history and its id.
- The routines are not the only way back: a message from the user starts a turn whatever their state, and the session
  start hook relaunches jobs when a container starts.

First applied on 9 October 2026, with all 49 registered jobs done and nothing running.

## What a firing does when a job is active

The prompts hold the detail. In outline: pull, run the resume script, arm one background wait per unfinished job
against that job's own markers, and report a finished job against the twin the documents name. The working documents
win over any summary in a routine's prompt, which has gone stale before: docs/confirmatory_runbook.md,
docs/preregistration.md, docs/graph_content_null_results.md, docs/membrane_potential_reach.md and
docs/literature_appraisal_staging.md (owned by the separate literature session: read it, never edit it).
