"""Generate docs/status.md, the one page that says what is current, from things that cannot be stale by construction:
the release manifests, configs/registered_rules.json, the git history of the decision log, the document register and
the git log (docs/documentation_process.md). It also renders the release counts into the block README.md carries
between `<!-- generated:release_counts -->` markers, so the README restates no number by hand.

What goes in and what stays out. In: the data releases and their counts, the registered rules and their amendments,
the decisions open and closed, the state of the documents (which are generated and when, which living documents are
due for review, which refer to files that changed after they were last reviewed). Out, deliberately: anything about a
container, a session or a running job. That state lives on a machine's disk (`bash scripts/resume_jobs.sh --status`)
and the day it was written into a committed document it was already wrong somewhere else (decision 9.47).

Usage:
  python experiments/generate_status.py            # writes docs/status.md and the README block
  python experiments/generate_status.py --check    # fails if either is out of date against the sources
  python experiments/generate_status.py --review   # the review diagnostics, which depend on today and on HEAD and are never committed
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import date
from pathlib import Path

from experiments.generate_decision_history import read_history, sort_key
from mechanistic_pathway_learning.documentation import (
    REPOSITORY_ROOT,
    STAMP_POLICY_DATE,
    dangling_paths,
    git_output,
    head_commit,
    last_change_date,
    load_registry,
    parse_stamp,
    provenance_stamp,
    referenced_repository_paths,
    sha256_of,
    unregistered_and_missing,
)

DEFAULT_OUTPUT_PATH = REPOSITORY_ROOT / "docs" / "status.md"
README_PATH = REPOSITORY_ROOT / "README.md"
RELEASES_ROOT = REPOSITORY_ROOT / "data" / "releases"
REGISTERED_RULES_PATH = REPOSITORY_ROOT / "configs" / "registered_rules.json"
AMENDMENT_HISTORY_PATH = REPOSITORY_ROOT / "docs" / "amendment_history.md"
README_BLOCK_START = "<!-- generated:release_counts -->"
README_BLOCK_END = "<!-- /generated:release_counts -->"
REVIEW_INTERVAL_DAYS = 60


def thousands(value) -> str:
    return f"{value:,}" if isinstance(value, int) else str(value)


# ----------------------------------------------------------------------------------------------------------------------
# sources


def release_summaries() -> list[dict]:
    summaries = []
    for release_dir in sorted(RELEASES_ROOT.glob("v*")):
        manifest_path = release_dir / "MANIFEST.json"
        if not manifest_path.exists():
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        tables = {}
        for relative_path, record in manifest.get("files", {}).items():
            if relative_path.endswith(".parquet") and record.get("rows") is not None:
                tables[relative_path[: -len(".parquet")]] = record["rows"]
        graph = {}
        for graph_name in ("graph_full", "graph"):
            summary_path = release_dir / graph_name / "graph_summary.json"
            if summary_path.exists():
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                graph[graph_name] = {"nodes": sum(summary.get("nodes_by_type", {}).values()), "edges": sum(summary.get("edges_by_relation", {}).values())}
        evidence = {}
        for evidence_name in ("evidence_full", "evidence_slice_sweeps", "evidence_two_entries"):
            summary_path = release_dir / evidence_name / "evidence_summary.json"
            if summary_path.exists():
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                evidence[evidence_name] = {"observations": summary.get("observations"), "perturbations": summary.get("distinct_perturbations"),
                                           "symptoms": len(summary.get("by_symptom", {}))}
        summaries.append({"version": manifest.get("version", release_dir.name), "created": manifest.get("created"), "git_commit": manifest.get("git_commit", "")[:7],
                          "files": len(manifest.get("files", {})), "amendments": len(manifest.get("amendments", [])), "tables": tables, "graph": graph, "evidence": evidence})
    return summaries


def registered_rules() -> dict:
    rules = json.loads(REGISTERED_RULES_PATH.read_text(encoding="utf-8"))
    ids = [key for key in rules if not key.startswith("_")]
    amendments = 0
    if AMENDMENT_HISTORY_PATH.exists():
        text = AMENDMENT_HISTORY_PATH.read_text(encoding="utf-8")
        amendments = len(re.findall(r"^### ", text, flags=re.MULTILINE))
    return {"count": len(ids), "ids": ids, "amendments": amendments}


def decisions() -> dict:
    history = read_history()
    timelines = sorted(history.timelines.values(), key=lambda timeline: sort_key(timeline.identifier))
    open_ids = [timeline.identifier for timeline in timelines if timeline.status == "open"]
    closed_ids = [timeline.identifier for timeline in timelines if timeline.status == "closed"]
    latest = sorted(timelines, key=lambda timeline: (timeline.recorded_date, timeline.identifier))[-5:]
    return {"open": open_ids, "closed": closed_ids, "dropped": [timeline.identifier for timeline in timelines if timeline.status == "dropped"],
            "latest": [(timeline.identifier, timeline.recorded_date, timeline.status) for timeline in latest],
            "duplicated": history.duplicated_now, "literal_newlines": history.literal_newlines_now}


def commits_touching_since(paths: set[str], since_commit: str) -> list[str]:
    """Registered paths that have commits after `since_commit`."""
    changed = []
    for path in sorted(paths):
        if not (REPOSITORY_ROOT / path).exists():
            continue
        try:
            log = git_output("log", "--format=%h", f"{since_commit}..HEAD", "--", path)
        except subprocess.CalledProcessError:  # an as_of_commit the log no longer holds, after a rebase
            changed.append(f"{path} (as_of_commit {since_commit} is not in the log)")
            continue
        if log.strip():
            changed.append(path)
    return changed


def documents(today: date) -> dict:
    registry = load_registry()
    unregistered, missing = unregistered_and_missing(registry)
    by_class: dict[str, list[str]] = {}
    for path, entry in registry.items():
        by_class.setdefault(entry.document_class, []).append(path)
    unstamped, stamped = [], []
    for path in by_class.get("generated", []):
        file_path = REPOSITORY_ROOT / path
        if not file_path.exists():
            continue
        stamp = parse_stamp(file_path.read_text(encoding="utf-8", errors="replace"))
        changed = last_change_date(path)
        (stamped if stamp else unstamped).append((path, changed.isoformat() if changed else "?", stamp))
    unstamped.sort(key=lambda item: item[1])
    living_rows = []
    for path in by_class.get("living", []):
        entry = registry[path].fields
        reviewed = entry.get("reviewed")
        reviewed_date = date.fromisoformat(reviewed) if reviewed else None
        overdue = reviewed_date is not None and (today - reviewed_date).days > REVIEW_INTERVAL_DAYS
        text = (REPOSITORY_ROOT / path).read_text(encoding="utf-8", errors="replace") if (REPOSITORY_ROOT / path).exists() else ""
        changed_since = commits_touching_since(referenced_repository_paths(text), entry["as_of_commit"]) if entry.get("as_of_commit") else []
        living_rows.append({"path": path, "reviewed": reviewed or "?", "as_of_commit": entry.get("as_of_commit", "?"), "overdue": overdue,
                            "changed_since": changed_since, "dangling": dangling_paths(text)})
    frozen_rows = []
    for path in by_class.get("frozen", []):
        entry = registry[path].fields
        actual = sha256_of(REPOSITORY_ROOT / path) if (REPOSITORY_ROOT / path).exists() else ""
        frozen_rows.append({"path": path, "intact": actual == entry.get("sha256"), "sha256": actual[:16]})
    return {"counts": {name: len(paths) for name, paths in sorted(by_class.items())}, "unregistered": unregistered, "missing": missing,
            "unstamped": unstamped, "stamped": len(stamped), "living": living_rows, "frozen": frozen_rows,
            "drafts": by_class.get("draft", []), "records": len(by_class.get("record", []))}


# ----------------------------------------------------------------------------------------------------------------------
# rendering


def render_release_counts(summaries: list[dict]) -> list[str]:
    lines = []
    for summary in summaries:
        parts = []
        for graph_name, counts in summary["graph"].items():
            parts.append(f"`{graph_name}` {thousands(counts['nodes'])} nodes and {thousands(counts['edges'])} edges")
        for evidence_name, counts in summary["evidence"].items():
            parts.append(f"`{evidence_name}` {thousands(counts['observations'])} observations over {thousands(counts['perturbations'])} perturbations and {counts['symptoms']} symptoms")
        lines.append(f"- `data/releases/{summary['version']}` (created {summary['created']}, commit {summary['git_commit']}, {summary['files']} files, "
                     f"{summary['amendments']} manifest amendment{'s' if summary['amendments'] != 1 else ''}): " + "; ".join(parts) + ".")
    return lines


def render(stamp: str, summaries: list[dict], rules: dict, decided: dict, docs: dict, today: date) -> str:
    lines = ["# Status (generated by experiments/generate_status.py — do not edit by hand)", "", stamp, "",
             "What is current in this repository, read from the release manifests, `configs/registered_rules.json`, the git history of the decision log "
             "and the document register. Nothing here describes a container, a session or a running job: that state is on a machine's disk "
             "(`bash scripts/resume_jobs.sh --status`), and a committed statement about it is stale the moment another machine reads it. The commit "
             "this page was generated at is in the stamp above; `git log -1` says where the code is now.", "",
             "## Data releases", "", *render_release_counts(summaries), ""]
    for summary in summaries:
        lines.append(f"Tables of `{summary['version']}` by row count: " + ", ".join(f"{name} {thousands(rows)}" for name, rows in sorted(summary["tables"].items())) + ".")
    lines.extend(["", "## Registration", "",
                  f"{rules['count']} registered rules in `configs/registered_rules.json` ({', '.join(rules['ids'])}), asserted against the code by "
                  f"`tests/test_registration_matches_the_code.py`; {rules['amendments']} amendment{'s' if rules['amendments'] != 1 else ''} since the baseline in "
                  "`docs/amendment_history.md`.", "",
                  "## Decisions", "",
                  f"{len(decided['open'])} open and {len(decided['closed'])} closed in `docs/open_decisions.md` (timeline in `docs/decision_history.md`). "
                  f"Open: {', '.join(decided['open'])}.", ""])
    lines.append("Most recently recorded: " + "; ".join(f"{identifier} ({recorded}, {status})" for identifier, recorded, status in decided["latest"]) + ".")
    log_defects = []
    if decided["duplicated"]:
        log_defects.append(f"ids used twice: {', '.join(decided['duplicated'])}")
    if decided["dropped"]:
        log_defects.append(f"ids dropped from the log: {', '.join(decided['dropped'])}")
    if decided["literal_newlines"]:
        log_defects.append(f"literal `\\n` on line{'s' if len(decided['literal_newlines']) > 1 else ''} {', '.join(map(str, decided['literal_newlines']))}")
    lines.extend(["", "Log defects: " + ("; ".join(log_defects) + "." if log_defects else "none."), ""])
    lines.extend(["## Documents", "",
                  "By class (`configs/documentation_registry.json`): " + ", ".join(f"{name} {count}" for name, count in docs["counts"].items()) + ".", ""])
    if docs["unregistered"] or docs["missing"]:
        lines.append("Register defects: " + "; ".join(filter(None, [f"unregistered {', '.join(docs['unregistered'])}" if docs["unregistered"] else "",
                                                                     f"registered but untracked {', '.join(docs['missing'])}" if docs["missing"] else ""])) + ".")
        lines.append("")
    for row in docs["frozen"]:
        lines.append(f"Frozen: `{row['path']}` {'intact' if row['intact'] else 'CHANGED'} (sha256 {row['sha256']}...).")
    lines.extend(["", f"Generated documents: {docs['stamped']} carry a provenance stamp (script, date, commit); {len(docs['unstamped'])} were generated before the "
                  f"stamp existed ({STAMP_POLICY_DATE.isoformat()}) and say only which script wrote them. Their last change in git is the only date they have; "
                  "the oldest describe runs or tables that live on the disk of the container that made them (every container of this project is kept), "
                  "not in this repository:", ""])
    for path, changed, _ in docs["unstamped"][:12]:
        lines.append(f"- `{path}` last changed {changed}")
    if len(docs["unstamped"]) > 12:
        lines.append(f"- and {len(docs['unstamped']) - 12} more; `python experiments/register_documents.py --report` lists every one")
    lines.extend(["", "Living documents (reviewed = last read against the repository, at the commit named). Whether a review is due depends on "
                  f"today's date ({REVIEW_INTERVAL_DAYS} days) and on what changed after that commit, so it is not written here, where it would be stale "
                  "the next day; `python experiments/generate_status.py --review` prints it.", "",
                  "| document | reviewed | as of |", "| --- | --- | --- |"])
    for row in docs["living"]:
        lines.append(f"| `{row['path']}` | {row['reviewed']} | {row['as_of_commit']} |")
    lines.extend(["", f"Records: {docs['records']} dated documents that are corrected by a dated addition rather than an edit. Drafts awaiting the owner: "
                  + (", ".join(f"`{path}`" for path in docs["drafts"]) if docs["drafts"] else "none") + ".", ""])
    return "\n".join(lines)


def render_review(docs: dict, today: date) -> str:
    """The diagnostics that depend on today and on HEAD: which living documents are due for a review and why. Printed,
    never committed, because a committed copy is wrong by the next commit or the next day."""
    lines = [f"Review diagnostics as of {today.isoformat()} at commit {head_commit()} (a review is due after {REVIEW_INTERVAL_DAYS} days, "
             "or when a file the document names changed after its as_of_commit):", ""]
    for row in docs["living"]:
        due = "overdue" if row["overdue"] else ("due" if row["changed_since"] else "current")
        changed = ", ".join(row["changed_since"][:8]) + (f" and {len(row['changed_since']) - 8} more" if len(row["changed_since"]) > 8 else "")
        lines.append(f"  {due:8s} {row['path']}  (reviewed {row['reviewed']} at {row['as_of_commit']})" + (f": changed since: {changed}" if changed else ""))
        if row["dangling"]:
            lines.append(f"           names paths not in the tree: {', '.join(row['dangling'])}")
    return "\n".join(lines)


def strip_stamp(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.startswith("Generated by `"))


def readme_with_block(readme_text: str, block_lines: list[str]) -> str:
    if README_BLOCK_START not in readme_text or README_BLOCK_END not in readme_text:
        raise ValueError(f"README.md has no {README_BLOCK_START} ... {README_BLOCK_END} block")
    before, rest = readme_text.split(README_BLOCK_START, 1)
    _, after = rest.split(README_BLOCK_END, 1)
    return before + README_BLOCK_START + "\n" + "\n".join(block_lines) + "\n" + README_BLOCK_END + after


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--review", action="store_true", help="print the review diagnostics (due reviews, files changed since) and write nothing")
    arguments = parser.parse_args()
    today = date.today()
    if arguments.review:
        print(render_review(documents(today), today))
        return 0
    summaries = release_summaries()
    rendered = render(provenance_stamp("experiments/generate_status.py", ["data/releases/*/MANIFEST.json", "configs/registered_rules.json",
                                                                          "configs/documentation_registry.json", "the git log"]),
                      summaries, registered_rules(), decisions(), documents(today), today)
    readme_text = README_PATH.read_text(encoding="utf-8")
    readme_rendered = readme_with_block(readme_text, render_release_counts(summaries))
    if arguments.check:
        stale = []
        if not arguments.output.exists() or strip_stamp(arguments.output.read_text(encoding="utf-8")) != strip_stamp(rendered):
            stale.append(str(arguments.output))
        if readme_rendered != readme_text:
            stale.append("README.md release_counts block")
        if stale:
            print("out of date: " + ", ".join(stale) + "; run experiments/generate_status.py")
            return 1
        print("docs/status.md and the README block are current")
        return 0
    arguments.output.write_text(rendered, encoding="utf-8")
    if readme_rendered != readme_text:
        README_PATH.write_text(readme_rendered, encoding="utf-8")
        print("updated the README release_counts block")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
