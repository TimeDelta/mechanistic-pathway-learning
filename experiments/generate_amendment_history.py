"""Generate docs/amendment_history.md from Amends: trailers in the git log (the owner's proposal, 15 October 2026).

The design, in his words: "create a copy of the preregistration as it exists today that stays unmodified in the
repository as starting file for the current state. Then the commit trailer works fine without any additional code
unit stuff."

So the registration's history is a baseline plus a derived list, and neither is hand-maintained prose:

  docs/preregistration_baseline_2026-10-15.md   the document exactly as registered on that day, never edited. Its
                                                SHA-256 is recorded below and asserted by
                                                tests/test_amendment_history.py, so "unmodified" is checkable rather
                                                than promised. No header was added to it, because a frozen
                                                registration should be the document and not a copy with a note on it.
  docs/amendment_history.md                     every change since, generated from the log.
  docs/preregistration.md                       the living document, free to be overwritten, because the history no
                                                longer lives inside it.

Why generated rather than written. A hand-kept amendment list can claim a date the commits do not support; a derived
one can be checked against the repository. And the defect it replaces was real: on 15 October the floor rule appeared
on 44 lines of the registration, 16 of them normative, across four revisions, none marked superseded, so a session
reading the earliest took it for current.

Why the output is committed rather than produced on demand. The project pulls with --rebase, and rebasing rewrites
commits. If a registration record lived only in a commit message, a squash or a dropped commit would lose it in
silence. Committing the generated file turns that loss into a visible diff.

What this does not catch, stated because it is the hole: a rule changed in code with no trailer written leaves no
entry, and the generated list would look complete anyway. That is what
tests/test_registration_matches_the_code.py is for - it fails when the code and configs/registered_rules.json
disagree, whether or not anybody wrote a trailer. The two are a pair; neither closes the gap alone.

Trailer format, one or more per commit:

    Amends: floor_rule
    From: realised half-width plus 0.005
    To: none

Each Amends: line opens an entry and the From:/To: lines that follow belong to it, so one commit may amend several
rules. The rule id must appear in configs/registered_rules.json, which is how a typo becomes an error rather than a
silent omission.

Usage:
  python experiments/generate_amendment_history.py            # writes docs/amendment_history.md
  python experiments/generate_amendment_history.py --check    # fails if the committed file is out of date
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from mechanistic_pathway_learning.documentation import provenance_stamp

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = REPOSITORY_ROOT / "docs" / "preregistration_baseline_2026-10-15.md"
REGISTERED_RULES_PATH = REPOSITORY_ROOT / "configs" / "registered_rules.json"
DEFAULT_OUTPUT_PATH = REPOSITORY_ROOT / "docs" / "amendment_history.md"

# The commit the baseline was frozen at. Entries are generated from this commit forward; everything earlier is the
# in-document amendment list of the baseline itself, which stands as the record for that period rather than being
# reconstructed from commits that never carried trailers.
BASELINE_COMMIT_DESCRIPTION = "the freeze of 15 October 2026"
RECORD_SEPARATOR = "\x1e"
FIELD_SEPARATOR = "\x1f"


def read_commits(since_commit: str | None) -> list[dict]:
    """Every commit, newest first, as {sha, iso_date, subject, body}."""
    revision_range = [f"{since_commit}..HEAD"] if since_commit else []
    completed = subprocess.run(
        ["git", "log", "--no-merges", f"--format=%H{FIELD_SEPARATOR}%cI{FIELD_SEPARATOR}%s{FIELD_SEPARATOR}%b{RECORD_SEPARATOR}", *revision_range],
        cwd=REPOSITORY_ROOT, capture_output=True, text=True, check=True)
    commits = []
    for record in completed.stdout.split(RECORD_SEPARATOR):
        fields = record.strip("\n").split(FIELD_SEPARATOR)
        if len(fields) == 4 and fields[0]:
            commits.append({"sha": fields[0], "iso_date": fields[1], "subject": fields[2], "body": fields[3]})
    return commits


def amendments_in_commit(commit: dict, known_rule_ids: set[str]) -> tuple[list[dict], list[str]]:
    """The amendment entries a commit's trailers declare, and the complaints about trailers that do not parse."""
    entries: list[dict] = []
    complaints: list[str] = []
    current: dict | None = None
    for line in commit["body"].splitlines():
        stripped = line.strip()
        if stripped.lower().startswith("amends:"):
            rule_id = stripped.split(":", 1)[1].strip()
            if rule_id not in known_rule_ids:
                complaints.append(f"{commit['sha'][:8]} amends {rule_id!r}, which is not in configs/registered_rules.json")
            current = {"rule_id": rule_id, "from": "", "to": "", "sha": commit["sha"],
                       "iso_date": commit["iso_date"], "subject": commit["subject"]}
            entries.append(current)
        elif current is not None and stripped.lower().startswith("from:"):
            current["from"] = stripped.split(":", 1)[1].strip()
        elif current is not None and stripped.lower().startswith("to:"):
            current["to"] = stripped.split(":", 1)[1].strip()
    return entries, complaints


def render_history(entries: list[dict], complaints: list[str], baseline_digest: str, rules: dict) -> str:
    """The generated document. Oldest first, because an amendment history reads forward."""
    lines = [
        "# Amendment history (generated by experiments/generate_amendment_history.py — do not edit by hand)",
        "",
        provenance_stamp("experiments/generate_amendment_history.py", ["the git log", "configs/registered_rules.json"]),
        "",
        "Every change to a registered rule since the baseline, taken from `Amends:` trailers in the git log. Edits to "
        "this file are overwritten by the next run; change a commit message, or add a trailer to a new commit, instead.",
        "",
        f"**Baseline.** `docs/preregistration_baseline_2026-10-15.md` is the registration exactly as it stood at "
        f"{BASELINE_COMMIT_DESCRIPTION}, SHA-256 `{baseline_digest}`, and is never edited. Amendments before that date "
        "are the in-document list inside the baseline itself; they were made before this mechanism existed and are not "
        "reconstructed here, because the commits of that period carry no trailers and a reconstruction would assert "
        "dates the log cannot support.",
        "",
        f"**Current values** live in `configs/registered_rules.json` ({len(rules)} rules) and are asserted against the "
        "implementation by `tests/test_registration_matches_the_code.py`. This file says how they got there; that file "
        "says they are still true.",
        "",
    ]
    if complaints:
        lines += ["## Trailers that did not parse", "",
                  "These are errors rather than history: a trailer naming a rule the registry does not hold means "
                  "either a typo or a rule that was never registered.", ""]
        lines += [f"- {complaint}" for complaint in complaints] + [""]

    if not entries:
        lines += ["## Amendments", "",
                  "None since the baseline. An empty list here means no commit has carried an `Amends:` trailer, which "
                  "is not the same as no rule having changed — a change made without a trailer leaves no entry, and "
                  "`tests/test_registration_matches_the_code.py` is what catches that.", ""]
        return "\n".join(lines)

    lines += ["## Amendments", ""]
    for entry in sorted(entries, key=lambda item: item["iso_date"]):
        rule = rules.get(entry["rule_id"], {})
        lines.append(f"### {entry['iso_date'][:10]} — `{entry['rule_id']}`")
        lines.append("")
        # Always emitted, even when both are missing: a dropped line would leave a reader unable to tell an
        # amendment recorded without values from one the generator lost.
        lines.append(f"- **From** {entry['from'] or '(unstated)'} **to** {entry['to'] or '(unstated)'}")
        if rule.get("means"):
            lines.append(f"- Governs: {rule['means']}")
        lines.append(f"- Commit `{entry['sha'][:12]}`: {entry['subject']}")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--since-commit", default=None,
                        help="generate from this commit forward; the default reads the whole log, which is harmless "
                             "because only commits carrying trailers produce entries")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument("--check", action="store_true",
                        help="write nothing and exit non-zero if the committed file differs from what would be "
                             "generated, so a stale history fails rather than misleads")
    arguments = parser.parse_args()

    if not BASELINE_PATH.exists():
        raise SystemExit(f"{BASELINE_PATH} is missing: the baseline is what the history is relative to")
    baseline_digest = hashlib.sha256(BASELINE_PATH.read_bytes()).hexdigest()
    rules = json.loads(REGISTERED_RULES_PATH.read_text())["rules"]

    entries: list[dict] = []
    complaints: list[str] = []
    for commit in read_commits(arguments.since_commit):
        commit_entries, commit_complaints = amendments_in_commit(commit, set(rules))
        entries += commit_entries
        complaints += commit_complaints

    rendered = render_history(entries, complaints, baseline_digest, rules) + "\n"
    if arguments.check:
        existing = arguments.output.read_text() if arguments.output.exists() else ""
        # the stamp carries the date and commit of the run that wrote the file, which is not what is being compared
        def without_stamp(text: str) -> str:
            return "\n".join(line for line in text.splitlines() if not line.startswith("Generated by `"))

        if without_stamp(existing) != without_stamp(rendered):
            raise SystemExit(f"{arguments.output} is out of date; run experiments/generate_amendment_history.py")
        print(f"{arguments.output} is up to date ({len(entries)} amendments)")
        return
    arguments.output.write_text(rendered)
    print(f"wrote {arguments.output}: {len(entries)} amendments, {len(complaints)} unparsed trailers, "
          f"baseline SHA-256 {baseline_digest[:12]}")


if __name__ == "__main__":
    main()
