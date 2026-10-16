"""The documentation process, asserted (docs/documentation_process.md; the owner, 15 October 2026: sessions "keep
getting tripped up by stale things").

Two kinds of test. The first kind runs against the repository as it stands: every tracked document is registered and
confirmed, frozen documents are byte-exact, generated documents name their generator and carry a stamp once
regenerated, living documents describe the repository and not a machine, records are not rewritten, the decision log
keeps every id it ever had, the paths and decisions the documents name exist, and the generated status and history
pages match what their generators say now. The second kind exercises the parsers and rules on small fixtures, so a
regression in a check shows up as a failing fixture rather than as a stale document passing.
"""
from __future__ import annotations

import json
import subprocess
from datetime import date

import pytest

from mechanistic_pathway_learning import documentation
from mechanistic_pathway_learning.documentation import (
    DOCUMENT_CLASSES,
    REPOSITORY_ROOT,
    STAMP_POLICY_DATE,
    append_only_violation,
    closed_decision_ids,
    dangling_links,
    dangling_paths,
    decision_ids_dropped_or_renumbered,
    decision_ids_in_log,
    duplicated_decision_ids,
    item_ids_in_registration,
    last_change_date,
    literal_newline_escapes,
    load_registry,
    names_its_generator,
    parse_stamp,
    provenance_stamp,
    referenced_repository_paths,
    session_state_phrases,
    sha256_of,
    unregistered_and_missing,
    unresolved_decision_references,
    unresolved_item_references,
)

DECISION_LOG = REPOSITORY_ROOT / "docs" / "open_decisions.md"
REGISTRATION = REPOSITORY_ROOT / "docs" / "preregistration.md"


@pytest.fixture(scope="module")
def registry():
    return load_registry()


@pytest.fixture(scope="module")
def known_decisions():
    return set(decision_ids_in_log(DECISION_LOG.read_text(encoding="utf-8")))


@pytest.fixture(scope="module")
def known_items():
    return item_ids_in_registration(REGISTRATION.read_text(encoding="utf-8"))


def document_text(path: str) -> str:
    return (REPOSITORY_ROOT / path).read_text(encoding="utf-8", errors="replace")


# ----------------------------------------------------------------------------------------------------------------------
# the register


def test_every_tracked_document_is_registered_and_every_entry_is_tracked(registry):
    unregistered, missing = unregistered_and_missing(registry)
    assert not unregistered, f"register these in configs/documentation_registry.json (python experiments/register_documents.py --add-missing): {unregistered}"
    assert not missing, f"registered but not tracked by git (python experiments/register_documents.py --prune-missing): {missing}"


def test_every_entry_has_a_known_class_and_is_confirmed(registry):
    unconfirmed = [path for path, entry in registry.items() if entry.fields.get("confirmed") is False]
    assert not unconfirmed, f"a guessed class is not a classification; set the class and remove the flag: {unconfirmed}"
    assert all(entry.document_class in DOCUMENT_CLASSES for entry in registry.values())


def test_living_entries_carry_a_review_date_and_an_as_of_commit(registry):
    for path, entry in registry.items():
        if entry.document_class == "living":
            assert entry.fields.get("reviewed"), f"{path}: a living document says when it was last read against the repository"
            date.fromisoformat(entry.fields["reviewed"])
            assert entry.fields.get("as_of_commit"), f"{path}: a living document says which commit it was reviewed at"


def test_record_entries_carry_the_date_they_were_recorded(registry):
    for path, entry in registry.items():
        if entry.document_class == "record":
            assert entry.fields.get("recorded"), f"{path}: a record says when it was recorded"


def test_the_register_explains_itself():
    registry = json.loads((REPOSITORY_ROOT / "configs" / "documentation_registry.json").read_text(encoding="utf-8"))
    assert registry["_what_this_is"].startswith("The register of every tracked Markdown document")
    assert "tests/test_documentation_process.py" in registry["_how_it_is_enforced"]
    assert date.fromisoformat(registry["policy_date"]) == STAMP_POLICY_DATE


# ----------------------------------------------------------------------------------------------------------------------
# frozen, generated, living, record


def test_frozen_documents_are_byte_exact(registry):
    frozen = {path: entry for path, entry in registry.items() if entry.document_class == "frozen"}
    assert frozen, "the registration baseline is registered as frozen"
    for path, entry in frozen.items():
        assert sha256_of(REPOSITORY_ROOT / path) == entry.fields["sha256"], f"{path} changed; a frozen document is replaced by a new frozen copy, never edited"


def test_the_frozen_baseline_digest_agrees_with_the_amendment_generator(registry):
    """Two places pin the same bytes on purpose: the register and the generator. They must agree."""
    pinned = registry["docs/preregistration_baseline_2026-10-15.md"].fields["sha256"]
    history = (REPOSITORY_ROOT / "docs" / "amendment_history.md").read_text(encoding="utf-8")
    assert pinned in history


def test_generated_documents_name_their_generator(registry):
    """The ten mechanism cards of 3 October predate the convention and say so in the register (names_generator: false
    with a why); any other generated document names its script in its first lines."""
    nameless = [path for path, entry in registry.items() if entry.document_class == "generated"
                and entry.fields.get("names_generator") is not False and not names_its_generator(document_text(path), entry.generator)]
    assert not nameless, f"a generated document says in its first lines which script wrote it: {nameless}"
    excused = [path for path, entry in registry.items() if entry.fields.get("names_generator") is False]
    assert all(registry[path].fields.get("why") for path in excused), "an excused document says why"


def test_generated_documents_regenerated_after_the_policy_date_carry_a_stamp(registry):
    """A table nobody can date is the stale statement that bit hardest. From the policy date on, every generated
    document that changes carries the script, the date, the commit and the inputs."""
    unstamped = []
    for path, entry in registry.items():
        if entry.document_class != "generated":
            continue
        changed = last_change_date(path)
        if changed is not None and changed > STAMP_POLICY_DATE and parse_stamp(document_text(path)) is None:
            unstamped.append(f"{path} (changed {changed})")
    assert not unstamped, f"regenerated without a provenance stamp (mechanistic_pathway_learning.documentation.provenance_stamp): {unstamped}"


def test_generated_documents_registered_generators_exist(registry):
    missing = [f"{path} -> {entry.generator}" for path, entry in registry.items() if entry.document_class == "generated"
               and entry.generator not in (None, "unknown") and not (REPOSITORY_ROOT / entry.generator).exists()]
    assert not missing, f"the generator named in the register is gone: {missing}"


def test_living_documents_describe_the_repository_not_a_machine(registry):
    offenders = {path: session_state_phrases(document_text(path)) for path, entry in registry.items() if entry.document_class == "living"}
    offenders = {path: phrases for path, phrases in offenders.items() if phrases}
    assert not offenders, f"container or session state in a living document; it belongs in a dated record or on the machine: {offenders}"


def test_append_only_records_are_not_rewritten(registry):
    violations = {path: append_only_violation(path) for path, entry in registry.items()
                  if entry.document_class == "record" and entry.fields.get("append_only")}
    violations = {path: why for path, why in violations.items() if why}
    assert not violations, violations


# ----------------------------------------------------------------------------------------------------------------------
# the decision log


def test_the_decision_log_keeps_every_id_it_had_at_head():
    dropped = decision_ids_dropped_or_renumbered()
    assert not dropped, f"a closed decision keeps its id, struck through; these ids left the log: {dropped}"


def test_the_decision_log_carries_no_literal_newline_escapes_beyond_the_known_line():
    """Line 204 held three rows behind literal backslash-n on 15 October 2026 and is left to the owner to repair; a
    second such line is a regression."""
    lines = literal_newline_escapes(DECISION_LOG.read_text(encoding="utf-8"))
    assert lines in ([], [204]), f"literal \\n hides rows from Markdown and from readers that split on newlines: lines {lines}"


def test_the_decision_log_uses_no_id_twice_beyond_the_known_duplicate():
    duplicated = duplicated_decision_ids(DECISION_LOG.read_text(encoding="utf-8"))
    assert duplicated in ([], ["9.5"]), f"an id names one decision: {duplicated}"


def test_the_decision_history_is_current():
    """--check in test form: a committed history that no longer matches the log's git history is the stale artifact
    this mechanism exists to stop."""
    result = subprocess.run(["python3", "experiments/generate_decision_history.py", "--check"], cwd=REPOSITORY_ROOT,
                            capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin"})
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_status_page_and_the_readme_block_are_current():
    result = subprocess.run(["python3", "experiments/generate_status.py", "--check"], cwd=REPOSITORY_ROOT,
                            capture_output=True, text=True, env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin:/usr/local/bin"})
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_status_page_describes_no_machine():
    text = (REPOSITORY_ROOT / "docs" / "status.md").read_text(encoding="utf-8")
    body = text.split("Nothing here describes a container", 1)[1]
    assert not session_state_phrases(body)


# ----------------------------------------------------------------------------------------------------------------------
# references: paths, links, decisions, items


def test_paths_named_by_living_and_generated_documents_exist(registry):
    """A record may name a file that has since moved; a living or generated document may not."""
    dangling = {}
    for path, entry in registry.items():
        if path == "docs/status.md":
            continue  # it reports the absent paths the living documents name, so it names them too
        if entry.document_class in ("living", "generated"):
            allowed = set(entry.fields.get("names_absent_paths", {}))
            missing = [missing_path for missing_path in dangling_paths(document_text(path)) if missing_path not in allowed]
            if missing:
                dangling[path] = missing
    assert not dangling, f"paths named that are not in the tree (or say why in the register's names_absent_paths): {dangling}"


def test_absent_paths_allowed_in_the_register_are_still_absent_and_explained(registry):
    """An allowance that outlives its reason is itself stale: once the file exists, the entry goes."""
    for path, entry in registry.items():
        for absent_path, why in entry.fields.get("names_absent_paths", {}).items():
            assert why, f"{path}: {absent_path} needs a reason"
            assert not (REPOSITORY_ROOT / absent_path).exists(), f"{path}: {absent_path} exists now; drop it from names_absent_paths"
            assert absent_path in referenced_repository_paths(document_text(path)), f"{path} no longer names {absent_path}; drop the allowance"


def test_relative_links_in_every_document_resolve(registry):
    broken = {}
    for path in registry:
        missing = dangling_links(document_text(path), path)
        if missing:
            broken[path] = missing
    assert not broken, broken


def test_decision_references_resolve_to_ids_in_the_log(registry, known_decisions):
    unresolved = {}
    for path, entry in registry.items():
        if entry.document_class in ("living", "generated"):
            missing = unresolved_decision_references(document_text(path), known_decisions)
            if missing:
                unresolved[path] = missing
    assert not unresolved, f"decision ids named that the log does not hold (renumbered, or never recorded): {unresolved}"


def test_item_references_in_living_documents_resolve_to_registration_items(registry, known_items):
    unresolved = {}
    for path, entry in registry.items():
        if entry.document_class == "living" and path not in ("README.md", "CLAUDE.md"):
            missing = unresolved_item_references(document_text(path), known_items)
            if missing:
                unresolved[path] = missing
    assert not unresolved, f"item numbers named that docs/preregistration.md does not hold: {unresolved}"


def test_the_code_names_only_decisions_the_log_holds(known_decisions):
    """Docstrings and comments cite decisions by id; a renumbering that leaves the code pointing at nothing is a
    stale statement in the place sessions read first."""
    unresolved = {}
    for source in sorted((REPOSITORY_ROOT / "mechanistic_pathway_learning").rglob("*.py")):
        missing = unresolved_decision_references(source.read_text(encoding="utf-8", errors="replace"), known_decisions)
        if missing:
            unresolved[str(source.relative_to(REPOSITORY_ROOT))] = missing
    assert not unresolved, unresolved


# ----------------------------------------------------------------------------------------------------------------------
# the rules on fixtures


def test_the_stamp_carries_script_date_commit_and_inputs(tmp_path):
    stamp = provenance_stamp("experiments/x.py", ["a.json", "the git log"], today=date(2026, 10, 16))
    parsed = parse_stamp(stamp)
    assert parsed["generator"] == "experiments/x.py" and parsed["date"] == "2026-10-16"
    assert "from a.json, the git log" in stamp and stamp.endswith("do not edit by hand.")


def test_the_stamp_marks_a_dirty_tree_rather_than_claiming_a_clean_commit(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "--allow-empty", "-m", "x"], check=True)
    assert "-dirty" not in provenance_stamp("experiments/x.py", root=tmp_path)
    (tmp_path / "loose.txt").write_text("x")
    assert "-dirty" not in provenance_stamp("experiments/x.py", root=tmp_path), "an untracked file is not a change to tracked content"
    subprocess.run(["git", "-C", str(tmp_path), "add", "loose.txt"], check=True)
    assert "-dirty" in provenance_stamp("experiments/x.py", root=tmp_path)


def test_the_earlier_stamp_form_of_the_descriptor_column_document_is_read():
    parsed = parse_stamp("# Title\n\nVersion 4, generated 2026-10-09 by experiments/write_descriptor_column_doc.py from the code at commit 42f040b.")
    assert parsed == {"date": "2026-10-09", "generator": "experiments/write_descriptor_column_doc.py", "commit": "42f040b"}


def test_a_stamp_past_the_first_lines_does_not_count():
    assert parse_stamp("a\nb\nc\nd\ne\nf\nGenerated by `experiments/x.py` on 2026-10-16 at commit abcdef1; do not edit by hand.") is None


@pytest.mark.parametrize("head", ["# T\n\nGenerated by `experiments/x.py`", "# T (generated by experiments/x.py)", "# T\n\nMeasured by `experiments/a/b.py` on x",
                                  "# T\n\nBuilt by experiments/x.py from y", "# T\n\nRead from openFDA by `experiments/x.py`"])
def test_the_forms_in_which_generators_name_themselves_are_recognised(head):
    assert names_its_generator(head, None)


def test_prose_that_merely_mentions_a_script_is_not_a_generated_document():
    assert not names_its_generator("# T\n\nThe table experiments/x.py writes is discussed below.", None)


def test_repository_paths_are_read_from_backticks_and_bare_mentions_and_runtime_paths_are_left_out():
    text = "see `docs/a.md`, configs/b.json and `data/processed/graph_full/nodes.parquet`; also `runs/x` and `docs/<version>/c.md` and `docs/figures/`."
    assert referenced_repository_paths(text) == {"docs/a.md", "configs/b.json", "docs/figures"}


def test_dangling_paths_and_links_are_reported(tmp_path):
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "real.md").write_text("x")
    assert dangling_paths("`docs/real.md` and `docs/gone.md`", root=tmp_path) == ["docs/gone.md"]
    assert dangling_links("[a](real.md) [b](gone.md) [c](https://x.y/z) [d](#anchor) [e](real.md#part)", "docs/index.md", root=tmp_path) == ["gone.md"]


def test_decision_ids_are_read_struck_or_not_and_through_literal_newlines():
    log = "| 1.2 | open |\n| ~~1.3~~ | closed |\\n| 9.7a | open |\n| not | an id |\n"
    assert decision_ids_in_log(log) == ["1.2", "1.3", "9.7a"]
    assert closed_decision_ids(log) == {"1.3"}
    assert literal_newline_escapes(log) == [2]
    assert duplicated_decision_ids(log + "| 1.2 | again |\n") == ["1.2"]


def test_decision_and_item_references_are_matched_by_word():
    known = {"9.5", "4.2b", "9.10a", "1.3"}
    assert unresolved_decision_references("decision 9.5, decisions 4.2b and decision 9.99; Decision 9.10(a) too; decision 1.3(c) is a sub-point; version 9.5 is not one", known) == ["9.99"]
    assert unresolved_item_references("item 14 and item 5c, Item 99, item 14b (a sub-point of 14)", {"14", "5c"}) == ["99"]
    assert item_ids_in_registration("**14. A rule**\n- **23a. A sub-item in a bullet**\n**5c. Another**\n*not* 7. one") == {"14", "23a", "5c"}


def test_session_state_phrases_are_found_case_insensitively_and_dated_lines_are_scoped():
    assert session_state_phrases("This Container does not hold the rebuild") == ["this container"]
    assert session_state_phrases("the repository holds the rebuild") == []
    assert session_state_phrases("| 2026-10-10 | this container did not hold the rebuild |") == [], "a dated row is a record of that day"
    assert session_state_phrases("a claim this session made from memory on 14 October") == []


def test_append_only_accepts_additions_and_rejects_edits(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    record = tmp_path / "record.md"
    record.write_text("# Record\n\nfirst entry\n")
    subprocess.run(["git", "-C", str(tmp_path), "add", "record.md"], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x"], check=True)
    assert append_only_violation("record.md", root=tmp_path) is None
    record.write_text("# Record\n\nfirst entry\n\nsecond entry\n")
    assert append_only_violation("record.md", root=tmp_path) is None
    record.write_text("# Record\n\nfirst entry, corrected\n")
    assert "line 3" in append_only_violation("record.md", root=tmp_path)


def test_a_record_not_yet_committed_has_nothing_to_violate(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "new.md").write_text("anything")
    assert append_only_violation("new.md", root=tmp_path) is None


def test_dropped_decision_ids_are_measured_against_head(tmp_path, monkeypatch):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "docs").mkdir()
    log = tmp_path / "docs" / "open_decisions.md"
    log.write_text("| 1.1 | a |\n| 1.2 | b |\n")
    subprocess.run(["git", "-C", str(tmp_path), "add", "."], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x"], check=True)
    log.write_text("| 1.1 | a |\n| ~~1.2~~ | b, closed |\n| 1.3 | c |\n")
    assert decision_ids_dropped_or_renumbered(root=tmp_path) == []
    log.write_text("| 1.1 | a |\n| 1.4 | b renumbered |\n")
    assert decision_ids_dropped_or_renumbered(root=tmp_path) == ["1.2"]


def test_the_decision_history_reads_recording_and_closing_from_the_log_history(tmp_path):
    from experiments.generate_decision_history import read_history
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "docs").mkdir()
    log = tmp_path / "docs" / "open_decisions.md"

    def commit(text: str, message: str) -> None:
        log.write_text(text)
        subprocess.run(["git", "-C", str(tmp_path), "add", "."], check=True)
        subprocess.run(["git", "-C", str(tmp_path), "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", message], check=True)

    commit("| 1.1 | a |\n", "record 1.1")
    commit("| 1.1 | a |\n| ~~1.2~~ | recorded closed |\n", "record 1.2 closed")
    commit("| ~~1.1~~ | a, closed |\n", "close 1.1, drop 1.2")
    history = read_history(root=tmp_path, include_working_tree=False)
    assert history.commits_read == 3
    one = history.timelines["1.1"]
    assert one.status == "closed" and one.closed_commit != one.recorded_commit and not one.recorded_closed
    two = history.timelines["1.2"]
    assert two.recorded_closed and two.status == "dropped" and two.closed_commit == two.recorded_commit


def test_the_status_generator_reads_the_release_counts_from_the_manifest():
    from experiments.generate_status import release_summaries
    summaries = {summary["version"]: summary for summary in release_summaries()}
    manifest = json.loads((REPOSITORY_ROOT / "data" / "releases" / "v0.4" / "MANIFEST.json").read_text(encoding="utf-8"))
    assert summaries["v0.4"]["tables"]["graph_full/edges"] == manifest["files"]["graph_full/edges.parquet"]["rows"]
    assert summaries["v0.4"]["graph"]["graph_full"]["edges"] == manifest["files"]["graph_full/edges.parquet"]["rows"]
    assert summaries["v0.4"]["evidence"]["evidence_full"]["observations"] == manifest["files"]["evidence_full/evidence_records.parquet"]["rows"]


def test_the_readme_block_is_replaced_between_its_markers_only():
    from experiments.generate_status import readme_with_block
    readme = "before\n<!-- generated:release_counts -->\nold\n<!-- /generated:release_counts -->\nafter\n"
    assert readme_with_block(readme, ["new"]) == "before\n<!-- generated:release_counts -->\nnew\n<!-- /generated:release_counts -->\nafter\n"
    with pytest.raises(ValueError):
        readme_with_block("no markers", ["new"])


def test_the_readme_restates_no_release_count_by_hand():
    """The counts live in the generated block; the prose around it names the manifest instead."""
    from experiments.generate_status import README_BLOCK_END, README_BLOCK_START
    readme = (REPOSITORY_ROOT / "README.md").read_text(encoding="utf-8")
    before, rest = readme.split(README_BLOCK_START, 1)
    block, after = rest.split(README_BLOCK_END, 1)
    manifest = json.loads((REPOSITORY_ROOT / "data" / "releases" / "v0.4" / "MANIFEST.json").read_text(encoding="utf-8"))
    # a count stated on a line that names the release directory is scoped to that immutable release (the CTD notice
    # names the 5,191 rows of data/releases/v0.4), which cannot go stale; an unscoped count can
    unscoped_lines = [line for line in (before + after).splitlines() if "data/releases/v" not in line]
    for relative_path, record in manifest["files"].items():
        if relative_path.endswith(".parquet") and record["rows"] and record["rows"] > 1000:
            formatted = f"{record['rows']:,}"
            assert not any(formatted in line for line in unscoped_lines), f"{formatted} ({relative_path}) is typed into the README prose"


def test_the_process_document_and_the_claude_rules_point_at_each_other():
    process = (REPOSITORY_ROOT / "docs" / "documentation_process.md").read_text(encoding="utf-8")
    rules = (REPOSITORY_ROOT / "CLAUDE.md").read_text(encoding="utf-8")
    assert "docs/documentation_process.md" in rules and "docs/status.md" in rules
    for required in ("configs/documentation_registry.json", "experiments/generate_status.py", "experiments/generate_decision_history.py",
                     "experiments/register_documents.py", "tests/test_documentation_process.py"):
        assert required in process, f"the process document names {required}"


def test_the_documentation_module_states_the_classes_it_checks():
    assert set(DOCUMENT_CLASSES) == {"generated", "living", "record", "frozen", "draft"}
    for name in DOCUMENT_CLASSES:
        assert name in documentation.__doc__
