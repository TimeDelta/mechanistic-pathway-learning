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

Four scheduled routines fire into the session that created them, staggered so a reclaimed container is restarted
within about fifteen minutes:

| id | schedule | name | fires into |
|---|---|---|---|
| `trig_01E9fxBn6XU4RMdEjPQEyndm` | `3 * * * *` | Job keep-alive :03 | `session_01M6JNfpSxQbjqBt6hD4RjMV` |
| `trig_017eJmHnGUq6q2CzEz3kt1Y5` | `19 * * * *` | Job keep-alive :19 | `session_01M6JNfpSxQbjqBt6hD4RjMV` |
| `trig_01LJeguqEjTcRN7y2VWKj2qA` | `35 * * * *` | Job keep-alive :35 | `session_01M6JNfpSxQbjqBt6hD4RjMV` |
| `trig_01LoV65i2vYRtdnCBpV3dY2D` | `51 * * * *` | Hourly job check-in: resume reclaimed jobs | `session_01M6JNfpSxQbjqBt6hD4RjMV` |

**A routine is bound to one session and cannot be repointed.** Each of the four carries
`persistent_session_id: session_01M6JNfpSxQbjqBt6hD4RjMV`, the session of 6 to 10 October 2026, and only
`create_trigger` sets that field: `update_trigger` takes a name, a schedule, an enabled state, a model and a prompt,
and no session. So a later session that enables one of these four sends the firing to the session in the last column,
not to itself, and is then left believing its jobs are covered when nothing is relaunching them. The four rows above
are a record of which routines exist, not instructions for a new session to enable.

**What a new session does when it registers a job.** It creates its own four routines with `create_trigger`, leaving
`persistent_session_id` out so each binds to the creating session, on the same four schedules, and adds a row per new
routine to the table above with its own session id. The four rows above stay disabled, and are never deleted: a pause
keeps the id, the prompt, the schedule and the run history, and a new routine's prompt is written by reading the
paused one it replaces (`get_trigger`, or `list_triggers`, reports the stored prompt). One substitution is needed
while copying a prompt: its pause-rule step names the four ids to disable, and those must become the new session's
own four, or the first firing with no active job disables another session's routines and leaves its own running. The
retired rows can be struck once the session they fire into is gone, but the ids stay in the table so a later reader
can still read the prompts.

**The pause rule** (the user's instruction of 9 October 2026: "Update the harness so that when there are no running
jobs, the check-ins get paused and only restarted when there are running jobs please"). The routines exist only to
relaunch registered jobs, so:

- **Disabled while no job is active.** A firing runs `scripts/resume_jobs.sh` first. If it reports `0 active jobs`,
  that firing disables the four routines bound to its own session (`update_trigger`, `enabled=false`), says so in one
  line and ends the turn without arming a wait. Each routine's prompt carries this as its first step, so the rule
  runs whether or not anyone remembers it.
- **Enabled while a job is registered.** A job is only ever registered from inside a session, so this belongs to
  whoever registers it: after `--register`, enable the four routines bound to that session, or create them as above
  when it has none. The script prints a reminder on every `--register` for that reason. Enabling a routine bound to
  another session is the failure the section above describes, so check the binding before enabling anything.
- Disabling is a pause, not a deletion. The ids above stay valid, the prompts and schedules are kept, and
  `ended_reason` stays empty. Deleting a routine would lose its history, its id and its prompt.
- The routines are not the only way back: a message from the user starts a turn whatever their state, and the session
  start hook relaunches jobs when a container starts, which is the one part of the harness a new session inherits
  with no setup.

The pause rule was first applied on 9 October 2026, with all 49 registered jobs done and nothing running. The binding
rule was written on 10 October 2026, at the user's question "The new session will automatically create its own check
in triggers though?": nothing creates them, and the instruction this section used to carry would have had a new
session enable four routines that wake a session it is not.

## What a firing does when a job is active

The prompts hold the detail. In outline: pull, run the resume script, arm one background wait per unfinished job
against that job's own markers, and report a finished job against the twin the documents name. The working documents
win over any summary in a routine's prompt, which has gone stale before: docs/confirmatory_runbook.md,
docs/preregistration.md, docs/graph_content_null_results.md, docs/membrane_potential_reach.md and
docs/literature_appraisal_staging.md (owned by the separate literature session: read it, never edit it).
