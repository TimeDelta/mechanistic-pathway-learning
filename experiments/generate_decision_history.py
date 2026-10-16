"""Generate docs/decision_history.md: when each decision in docs/open_decisions.md was recorded and when it closed,
read from the git history of the log rather than from the dates written in its rows (docs/documentation_process.md).

Why derived. The owner's rule is that "a decision is recorded here in the turn it is taken" (13 October 2026), and the
log itself admits that five entries were recorded after the decision behind them had been made and built. A row can
carry any date; the commit that first holds the row cannot. So this file states, for every id the log has ever
carried: the commit and date it first appeared, the commit and date its id was struck through (closed), whether it
is open now, and the defects a hand-kept log accumulates (an id used twice, an id that vanished, a literal "\\n"
that hides rows from a reader that splits on real newlines). The amendment history does the same for the
registration from `Amends:` trailers; this does it for the decision log from the log's own diffs, so no trailer is
needed and no session can forget to write one.

Usage:
  python experiments/generate_decision_history.py            # writes docs/decision_history.md
  python experiments/generate_decision_history.py --check    # fails if the committed file is out of date
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

from mechanistic_pathway_learning.documentation import (
    REPOSITORY_ROOT,
    closed_decision_ids,
    decision_ids_in_log,
    duplicated_decision_ids,
    git_output,
    literal_newline_escapes,
    provenance_stamp,
)

DEFAULT_OUTPUT_PATH = REPOSITORY_ROOT / "docs" / "decision_history.md"
LOG_RELATIVE_PATH = "docs/open_decisions.md"


@dataclass
class DecisionTimeline:
    identifier: str
    recorded_commit: str
    recorded_date: str
    closed_commit: str | None = None
    closed_date: str | None = None
    dropped_commit: str | None = None
    dropped_date: str | None = None
    recorded_closed: bool = False  # first seen already struck through: closed before it was ever recorded open

    @property
    def status(self) -> str:
        if self.dropped_commit:
            return "dropped"
        return "closed" if self.closed_commit else "open"


@dataclass
class LogHistory:
    timelines: dict[str, DecisionTimeline] = field(default_factory=dict)
    commits_read: int = 0
    duplicated_now: list[str] = field(default_factory=list)
    literal_newlines_now: list[int] = field(default_factory=list)


def commits_touching(path: str = LOG_RELATIVE_PATH, root: Path = REPOSITORY_ROOT) -> list[tuple[str, str]]:
    """(short sha, date) of every commit that touched the log, oldest first."""
    lines = git_output("log", "--reverse", "--format=%h %ad", "--date=short", "--", path, cwd=root).splitlines()
    return [tuple(line.split()) for line in lines if line.strip()]


def read_history(path: str = LOG_RELATIVE_PATH, root: Path = REPOSITORY_ROOT, include_working_tree: bool = True) -> LogHistory:
    history = LogHistory()
    seen_previously: set[str] = set()
    versions = [(sha, day, git_output("show", f"{sha}:{path}", cwd=root)) for sha, day in commits_touching(path, root)]
    if include_working_tree:
        working = (root / path).read_text(encoding="utf-8")
        if not versions or working != versions[-1][2]:
            versions.append(("working-tree", "uncommitted", working))
    for sha, day, text in versions:
        history.commits_read += 1
        present = set(decision_ids_in_log(text))
        closed = closed_decision_ids(text)
        for identifier in sorted(present, key=sort_key):
            timeline = history.timelines.get(identifier)
            if timeline is None:
                timeline = DecisionTimeline(identifier, sha, day, recorded_closed=identifier in closed)
                history.timelines[identifier] = timeline
            if timeline.dropped_commit and identifier in present:
                timeline.dropped_commit = timeline.dropped_date = None  # it came back
            if identifier in closed and timeline.closed_commit is None:
                timeline.closed_commit, timeline.closed_date = sha, day
        for identifier in seen_previously - present:
            timeline = history.timelines[identifier]
            if timeline.dropped_commit is None:
                timeline.dropped_commit, timeline.dropped_date = sha, day
        seen_previously = present
    if versions:
        history.duplicated_now = duplicated_decision_ids(versions[-1][2])
        history.literal_newlines_now = literal_newline_escapes(versions[-1][2])
    return history


def sort_key(identifier: str) -> tuple:
    match = re.match(r"(\d+)\.(\d+)([a-z]?)", identifier)
    return (int(match.group(1)), int(match.group(2)), match.group(3)) if match else (999, 0, identifier)


def render(history: LogHistory, stamp: str) -> str:
    timelines = sorted(history.timelines.values(), key=lambda timeline: sort_key(timeline.identifier))
    open_now = [timeline for timeline in timelines if timeline.status == "open"]
    closed = [timeline for timeline in timelines if timeline.status == "closed"]
    dropped = [timeline for timeline in timelines if timeline.status == "dropped"]
    same_commit = [timeline for timeline in closed if timeline.closed_commit == timeline.recorded_commit]
    lines = [
        "# Decision history (generated by experiments/generate_decision_history.py — do not edit by hand)",
        "",
        stamp,
        "",
        f"When each decision of `{LOG_RELATIVE_PATH}` was recorded and when it closed, read from the {history.commits_read} versions of the log in the git "
        "history (plus the working tree when it differs from HEAD). The commit that first holds a row is the earliest the decision can have been "
        "recorded; a date written inside the row is the author's claim and is not read here. A closed decision keeps its id, struck through, so "
        "\"dropped\" means an id the log once carried and no longer does, which is the defect this file exists to show.",
        "",
        f"**Now:** {len(open_now)} open, {len(closed)} closed, {len(dropped)} dropped. "
        f"{len(same_commit)} of the closed decisions were recorded and closed in the same commit, which is the recording lag the log's own header describes.",
        "",
    ]
    defects = []
    if history.duplicated_now:
        defects.append(f"ids used by more than one row: {', '.join(history.duplicated_now)}")
    if history.literal_newlines_now:
        defects.append(f"literal `\\n` on line{'s' if len(history.literal_newlines_now) > 1 else ''} {', '.join(map(str, history.literal_newlines_now))}, "
                       "which hides the rows after it from Markdown and from any reader that splits on real newlines")
    if dropped:
        defects.append("dropped ids: " + ", ".join(f"{timeline.identifier} (last seen before {timeline.dropped_commit}, {timeline.dropped_date})" for timeline in dropped))
    lines.append("**Defects in the log as it stands:** " + ("; ".join(defects) + "." if defects else "none."))
    lines.extend(["", "## Timeline", "", "| id | recorded | closed | status |", "| --- | --- | --- | --- |"])
    for timeline in timelines:
        recorded = f"{timeline.recorded_date} ({timeline.recorded_commit})" + (" already closed" if timeline.recorded_closed else "")
        closed_cell = f"{timeline.closed_date} ({timeline.closed_commit})" if timeline.closed_commit else "-"
        if timeline.status == "dropped":
            closed_cell = f"dropped at {timeline.dropped_commit} ({timeline.dropped_date})"
        lines.append(f"| {timeline.identifier} | {recorded} | {closed_cell} | {timeline.status} |")
    lines.append("")
    return "\n".join(lines)


def strip_stamp(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.startswith("Generated by `"))


PENDING_CELL = "uncommitted (working-tree)"


def only_pending_commits_resolved(committed_lines: list[str], rendered_lines: list[str]) -> bool:
    """True when the committed file differs from the rendering only where a row said a change was still in the
    working tree and the rendering now names the commit that holds it.

    The history is generated before the commit that closes a decision exists, so the committed copy says
    "uncommitted (working-tree)" for that closing and the next regeneration fills in the commit. That one-commit lag
    is the mechanism working as designed, not a stale document, so --check accepts it; any other difference fails.
    """
    if len(committed_lines) != len(rendered_lines):
        return False
    for committed, rendered in zip(committed_lines, rendered_lines):
        if committed == rendered:
            continue
        if PENDING_CELL not in committed or not committed.startswith("| ") or not rendered.startswith("| "):
            return False
        committed_cells = [cell.strip() for cell in committed.strip("|").split("|")]
        rendered_cells = [cell.strip() for cell in rendered.strip("|").split("|")]
        if len(committed_cells) != len(rendered_cells):
            return False
        for committed_cell, rendered_cell in zip(committed_cells, rendered_cells):
            if committed_cell != rendered_cell and not (PENDING_CELL in committed_cell and re.search(r"\(\w{7,40}\)", rendered_cell)):
                return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    history = read_history()
    rendered = render(history, provenance_stamp("experiments/generate_decision_history.py", [f"the git history of {LOG_RELATIVE_PATH}"]))
    if arguments.check:
        if not arguments.output.exists():
            print(f"{arguments.output} is missing; run experiments/generate_decision_history.py")
            return 1
        committed_lines = strip_stamp(arguments.output.read_text(encoding="utf-8")).splitlines()
        rendered_lines = strip_stamp(rendered).splitlines()
        if committed_lines != rendered_lines and not only_pending_commits_resolved(committed_lines, rendered_lines):
            print(f"{arguments.output} is out of date against the git history of {LOG_RELATIVE_PATH}; run experiments/generate_decision_history.py")
            return 1
        print(f"{arguments.output} is current")
        return 0
    arguments.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {arguments.output}: {len(history.timelines)} decisions over {history.commits_read} versions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
