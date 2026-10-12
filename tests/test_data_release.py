"""Data releases: copies, CSV twins, MedDRA term labels blanked on label-derived rows, manifest checksums verified."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mechanistic_pathway_learning.data_release import build_release, strip_licensed_term_labels, verify_release


def write_inputs(root: Path) -> tuple[Path, Path]:
    graph = root / "graph"
    graph.mkdir()
    pd.DataFrame({"node_id": ["GENE:A", "MAR1"], "node_type": ["gene", "reaction"]}).to_parquet(graph / "nodes.parquet", index=False)
    pd.DataFrame({"source_id": ["MAR1"], "target_id": ["GENE:A"], "relation_type": ["catalyzed_by"]}).to_parquet(graph / "edges.parquet", index=False)
    (graph / "relation_types.json").write_text(json.dumps(["catalyzed_by"]))
    evidence = root / "evidence"
    evidence.mkdir()
    pd.DataFrame({"perturbation_id": ["A", "CID1"], "symptom": ["anxiety", "anxiety"], "weight": [1.0, 0.6]}).to_parquet(evidence / "evidence_records.parquet", index=False)
    pd.DataFrame({
        "source": ["HPO-OMIM", "SIDER-label", "SIDER-indication"],
        "source_term_id": ["HP:0000739", "C0003467", "C0003469"],
        "source_term_label": ["Anxiety", "Anxiety (MedDRA PT)", "Anxiety disorder (MedDRA PT)"],
        "symptom": ["anxiety", "anxiety", "anxiety"],
    }).to_parquet(evidence / "evidence_reports.parquet", index=False)
    (evidence / "evidence_summary.json").write_text(json.dumps({"observations": 2}))
    return graph, evidence


def test_strip_blanks_only_label_derived_rows() -> None:
    reports = pd.DataFrame({"source": ["HPO-OMIM", "SIDER-label", "OnSIDES-label"], "source_term_label": ["Anxiety", "x", "y"]})
    stripped, count = strip_licensed_term_labels(reports)
    assert count == 2 and stripped.source_term_label.tolist() == ["Anxiety", "", ""]
    assert reports.source_term_label.tolist() == ["Anxiety", "x", "y"]  # the input is not modified


def test_build_then_verify_and_detect_tampering(tmp_path: Path) -> None:
    graph, evidence = write_inputs(tmp_path)
    sources = tmp_path / "data_sources.md"
    sources.write_text("| Source | Use | Access | License note | Pinned version |\n|---|---|---|---|---|\n| HPO | E1 | public | open | release 2026-09-01 |\n")
    release = build_release("vtest", tmp_path / "releases", {"graph": graph}, {"evidence": evidence}, sources, None, notes="fixture")
    manifest = json.loads((release / "MANIFEST.json").read_text())
    assert manifest["licensed_term_labels_blanked"] == 2
    assert manifest["source_pins"] == [{"source": "HPO", "use": "E1", "pinned_version": "release 2026-09-01"}]
    assert (release / "evidence" / "evidence_reports.csv").exists() and (release / "graph" / "nodes.csv").exists()
    released_reports = pd.read_parquet(release / "evidence" / "evidence_reports.parquet")
    assert released_reports.source_term_label.tolist() == ["Anxiety", "", ""]
    assert manifest["files"]["evidence/evidence_records.parquet"]["rows"] == 2 and manifest["files"]["evidence/evidence_records.csv"]["rows"] == 2
    assert verify_release(release) == []
    (release / "graph" / "relation_types.json").write_text("[]")
    problems = verify_release(release)
    assert problems == ["graph/relation_types.json: checksum differs from the manifest"]
    (release / "stray.txt").write_text("x")
    assert any(problem.startswith("stray.txt") for problem in verify_release(release))


def test_extra_directories_drop_article_text_but_keep_metadata(tmp_path: Path) -> None:
    graph, evidence = write_inputs(tmp_path)
    literature = tmp_path / "literature"
    literature.mkdir()
    pd.DataFrame({"pmid": ["1"], "year": [2020], "title": ["A title"], "abstract": ["Long text"], "species_ids": ["9606"], "publication_types": ["Review"]}).to_parquet(literature / "documents.parquet", index=False)
    release = build_release("vdocs", tmp_path / "releases", {"graph": graph}, {"evidence": evidence}, None, None, extra_directories={"literature": literature})
    released = pd.read_parquet(release / "literature" / "documents.parquet")
    assert list(released.columns) == ["pmid", "year", "species_ids", "publication_types"]
    assert verify_release(release) == []


def test_ctd_terms_are_added_to_an_existing_release_without_touching_its_tables(tmp_path: Path) -> None:
    from mechanistic_pathway_learning.data_release import RELEASE_LICENSE, add_ctd_terms_to_release

    graph, evidence = write_inputs(tmp_path)
    literature = tmp_path / "literature"
    literature.mkdir()
    pd.DataFrame({"chemical_mesh_id": ["D000077443"], "mesh_descriptor": ["D001008"], "source": ["CTD-curated"]}).to_parquet(literature / "ctd_relations_filtered.parquet", index=False)
    (literature / "ctd_relations_summary.json").write_text(json.dumps({"ctd_pin": {"report_created": "Tue Sep 29 13:17:30 EDT 2026", "sha256": "abc", "checked_at": "2026-10-06T19:03:55Z"}}))
    release = build_release("vctd", tmp_path / "releases", {"graph": graph}, {"evidence": evidence}, None, None, extra_directories={"literature": literature})
    checksums_before = {name: entry["sha256"] for name, entry in json.loads((release / "MANIFEST.json").read_text())["files"].items()}

    notice = add_ctd_terms_to_release(release, today="2026-10-11")
    manifest = json.loads((release / "MANIFEST.json").read_text())
    text = notice.read_text()
    assert notice == release / "literature" / "CTD_TERMS.md" and "literature/CTD_TERMS.md" in manifest["files"]
    assert "You must notify CTD, and describe your use of our data" in text and "10.1093/nar/gkae883" in text
    assert "report created Tue Sep 29 13:17:30 EDT 2026" in text and "on 2026-10-06" in text and "It holds 1 rows" in text
    assert "http://ctdbase.org/detail.go?type=chem&acc=D000077443" in text and "http://ctdbase.org/detail.go?type=disease&db=MESH&acc=D001008" in text
    assert manifest["license"] == RELEASE_LICENSE and "CTD" in manifest["license"]
    assert manifest["amendments"][0]["date"] == "2026-10-11" and "previous_license" in manifest["amendments"][0]
    assert {name: entry["sha256"] for name, entry in manifest["files"].items() if name != "literature/CTD_TERMS.md"} == checksums_before
    assert verify_release(release) == []
    add_ctd_terms_to_release(release, today="2026-10-11")  # running it twice records one amendment
    assert len(json.loads((release / "MANIFEST.json").read_text())["amendments"]) == 1


def test_a_release_without_ctd_data_gets_no_ctd_notice(tmp_path: Path) -> None:
    from mechanistic_pathway_learning.data_release import add_ctd_terms_to_release

    graph, evidence = write_inputs(tmp_path)
    release = build_release("vplain", tmp_path / "releases", {"graph": graph}, {"evidence": evidence}, None, None)
    assert add_ctd_terms_to_release(release) is None and verify_release(release) == []


def test_ctd_page_addresses_drop_the_mesh_prefix() -> None:
    from mechanistic_pathway_learning.evidence.ctd_terms import ctd_chemical_page, ctd_disease_page

    assert ctd_chemical_page("MESH:D000077443") == "http://ctdbase.org/detail.go?type=chem&acc=D000077443"
    assert ctd_disease_page("D001008") == "http://ctdbase.org/detail.go?type=disease&db=MESH&acc=D001008"
