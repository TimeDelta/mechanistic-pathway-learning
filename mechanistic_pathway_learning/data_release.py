"""Versioned releases of the processed tables (graph and evidence), committed to git.

A release freezes the exact inputs of an experiment version so a checkout trains and
evaluates without any download, and so the pre-registration refers to checksummed
tables rather than to source pins alone. Each release directory holds copies of the
processed tables, CSV twins of the small tables for readable diffs, and MANIFEST.json
with the SHA-256, size and row count of every file, the git commit that produced it,
the source pins copied from docs/data_sources.md and a license statement.

Licensing. The report table carries MedDRA term names on rows derived from drug labels
(SIDER, OnSIDES); MedDRA terms may not be redistributed, so those names are blanked
before copying (the term identifier and the target symptom stay). Everything else in a
release is derivable from CC BY, CC BY-SA or public sources, and the SIDER-derived
rows make the release CC BY-NC-SA as a whole (design section 4.4).

Raw downloads are never part of a release; they are two gigabytes, several contain
licensed term tables, and docs/data_sources.md pins them for re-download.
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import date
from pathlib import Path

import pandas as pd

LABEL_DERIVED_SOURCES = ("SIDER-label", "SIDER-indication", "OnSIDES-label")
LICENSED_TERM_LABEL_COLUMN = "source_term_label"
GRAPH_FILES = ("nodes.parquet", "edges.parquet", "relation_types.json", "graph_summary.json")
EVIDENCE_FILES = ("evidence_records.parquet", "evidence_reports.parquet", "evidence_summary.json", "unmapped_records.parquet")
CSV_TWIN_MAX_ROWS = 10_000  # tables up to this size also get a CSV twin for readable diffs (evidence tables; not the node and edge tables)
MANIFEST_NAME = "MANIFEST.json"


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_count(path: Path) -> int | None:
    if path.suffix == ".parquet":
        import pyarrow.parquet as parquet_module

        return int(parquet_module.read_metadata(path).num_rows)
    if path.suffix == ".csv":
        with open(path, encoding="utf-8") as handle:
            return sum(1 for _ in handle) - 1
    return None


def strip_licensed_term_labels(reports: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Blank the MedDRA term names of label-derived report rows; returns the table and the number of rows blanked."""
    stripped = reports.copy()
    if LICENSED_TERM_LABEL_COLUMN not in stripped.columns or "source" not in stripped.columns:
        return stripped, 0
    mask = stripped["source"].isin(LABEL_DERIVED_SOURCES)
    stripped.loc[mask, LICENSED_TERM_LABEL_COLUMN] = ""
    return stripped, int(mask.sum())


def source_pins_from_table(data_sources_path: Path) -> list[dict[str, str]]:
    """Rows of the pinned-version table in docs/data_sources.md as dictionaries (source, use, pinned version)."""
    pins: list[dict[str, str]] = []
    if not data_sources_path.exists():
        return pins
    header: list[str] | None = None
    for line in data_sources_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if header is None:
            header = cells
            continue
        if all(set(cell) <= {"-", ":"} for cell in cells):
            continue
        if len(cells) == len(header):
            pins.append({"source": cells[0], "use": cells[1], "pinned_version": cells[-1]})
    return pins


def git_commit_hash(repository_root: Path) -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=repository_root, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def copy_table_with_twin(source: Path, destination: Path, transform=None) -> list[Path]:
    """Copy one table (optionally transformed) and write a CSV twin for small parquet tables; returns the written paths."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    if source.suffix == ".parquet":
        table = pd.read_parquet(source)
        if transform is not None:
            table = transform(table)
        table.to_parquet(destination, index=False)
        written.append(destination)
        if len(table) <= CSV_TWIN_MAX_ROWS:
            twin = destination.with_suffix(".csv")
            table.to_csv(twin, index=False)
            written.append(twin)
    else:
        shutil.copyfile(source, destination)
        written.append(destination)
    return written


def build_release(
    version: str,
    output_root: Path,
    graph_directories: dict[str, Path],
    evidence_directories: dict[str, Path],
    data_sources_path: Path | None = None,
    repository_root: Path | None = None,
    notes: str = "",
    extra_directories: dict[str, Path] | None = None,
) -> Path:
    """Write data/releases/<version>/ with copies of the named graph and evidence directories, any extra directories
    (every parquet, csv, json and md file at their top level, with licensed term labels blanked where the columns exist) and a manifest."""
    release_directory = output_root / version
    if release_directory.exists():
        shutil.rmtree(release_directory)
    release_directory.mkdir(parents=True)
    written: list[Path] = []
    stripped_rows_total = 0
    for name, directory in graph_directories.items():
        for file_name in GRAPH_FILES:
            source = directory / file_name
            if source.exists():
                written += copy_table_with_twin(source, release_directory / name / file_name)
    for name, directory in evidence_directories.items():
        for file_name in EVIDENCE_FILES:
            source = directory / file_name
            if not source.exists():
                continue
            if file_name == "evidence_reports.parquet":
                def strip_and_count(table: pd.DataFrame) -> pd.DataFrame:
                    nonlocal stripped_rows_total
                    stripped, count = strip_licensed_term_labels(table)
                    stripped_rows_total += count
                    return stripped
                written += copy_table_with_twin(source, release_directory / name / file_name, strip_and_count)
            else:
                written += copy_table_with_twin(source, release_directory / name / file_name)
    for name, directory in (extra_directories or {}).items():
        if not directory.exists():
            continue
        for source in sorted(directory.iterdir()):
            if source.suffix not in (".parquet", ".csv", ".json", ".md") or not source.is_file():
                continue
            if source.suffix == ".parquet":
                def strip_if_present(table: pd.DataFrame) -> pd.DataFrame:
                    nonlocal stripped_rows_total
                    stripped, count = strip_licensed_term_labels(table)
                    stripped_rows_total += count
                    return stripped
                written += copy_table_with_twin(source, release_directory / name / source.name, strip_if_present)
            else:
                written += copy_table_with_twin(source, release_directory / name / source.name)
    manifest = {
        "version": version,
        "created": date.today().isoformat(),
        "git_commit": git_commit_hash(repository_root) if repository_root is not None else None,
        "license": "Derived data; non-commercial and share-alike (CC BY-NC-SA) because SIDER-derived rows are included; MedDRA term names removed from label-derived report rows; see docs/experiment_design.md section 4.4 and docs/alternative_sources.md",
        "licensed_term_labels_blanked": stripped_rows_total,
        "notes": notes,
        "source_pins": source_pins_from_table(data_sources_path) if data_sources_path is not None else [],
        "files": {
            str(path.relative_to(release_directory)): {"sha256": sha256_of_file(path), "bytes": path.stat().st_size, "rows": row_count(path)}
            for path in sorted(written)
        },
    }
    (release_directory / MANIFEST_NAME).write_text(json.dumps(manifest, indent=1))
    return release_directory


def verify_release(release_directory: Path) -> list[str]:
    """Recompute every checksum in the manifest; returns a list of problems (empty when the release is intact)."""
    manifest_path = release_directory / MANIFEST_NAME
    if not manifest_path.exists():
        return [f"{manifest_path} is missing"]
    manifest = json.loads(manifest_path.read_text())
    problems: list[str] = []
    for relative_path, entry in manifest["files"].items():
        path = release_directory / relative_path
        if not path.exists():
            problems.append(f"{relative_path}: missing")
            continue
        if sha256_of_file(path) != entry["sha256"]:
            problems.append(f"{relative_path}: checksum differs from the manifest")
    for path in release_directory.rglob("*"):
        if path.is_file() and path.name != MANIFEST_NAME and str(path.relative_to(release_directory)) not in manifest["files"]:
            problems.append(f"{path.relative_to(release_directory)}: not in the manifest")
    return problems
