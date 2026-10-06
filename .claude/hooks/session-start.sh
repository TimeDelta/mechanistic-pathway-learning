#!/bin/bash
# SessionStart hook for Claude Code cloud sessions: make sure the package and the test tools import, then relaunch
# the unfinished background jobs registered under runs/jobs/ (scripts/resume_jobs.sh). A usage limit can stop the
# session, after which the container is reclaimed and every running process dies; the session's disk survives, so the
# jobs resume from their checkpoints and caches. The job report printed here becomes part of the session's context.
set -uo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"

if ! python -c "import mechanistic_pathway_learning, torch, pandas, pyarrow, pytest" > /dev/null 2>&1; then
  echo "session-start: installing the package with its test dependencies"
  pip install --quiet -e ".[dev]" pyarrow || echo "session-start: dependency install failed; jobs are still resumed"
fi

bash scripts/resume_jobs.sh
exit 0
