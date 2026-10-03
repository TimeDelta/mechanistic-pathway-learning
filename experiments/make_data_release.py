"""Build or verify a versioned release of the processed tables (mechanistic_pathway_learning/data_release.py).

Usage:
  python experiments/make_data_release.py --version v0.4                 # build data/releases/v0.4 from data/processed
  python experiments/make_data_release.py --version v0.4 --verify-only   # recompute checksums against MANIFEST.json

Training from a release needs no download:
  python experiments/run_main_model.py --graph-dir data/releases/v0.4/graph_full --evidence-dir data/releases/v0.4/evidence_full ...
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mechanistic_pathway_learning.data_release import build_release, verify_release


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="release name, for example v0.4")
    parser.add_argument("--output-root", type=Path, default=Path("data/releases"))
    parser.add_argument("--graph", nargs="*", default=["graph_full=data/processed/graph_full", "graph=data/processed/graph"], help="name=directory pairs")
    parser.add_argument("--evidence", nargs="*", default=["evidence_full=data/processed/evidence_full", "evidence_slice_sweeps=data/processed/evidence", "evidence_two_entries=data/processed/evidence_two_entries"], help="name=directory pairs")
    parser.add_argument("--data-sources", type=Path, default=Path("docs/data_sources.md"))
    parser.add_argument("--extra", nargs="*", default=["onsides=data/processed/onsides", "literature=data/processed/literature"], help="name=directory pairs whose top-level parquet, csv, json and md files are released")
    parser.add_argument("--notes", default="")
    parser.add_argument("--verify-only", action="store_true")
    arguments = parser.parse_args()
    release_directory = arguments.output_root / arguments.version
    if arguments.verify_only:
        problems = verify_release(release_directory)
        print("\n".join(problems) if problems else f"{release_directory}: every file matches the manifest")
        sys.exit(1 if problems else 0)
    pairs = lambda items: {item.split("=", 1)[0]: Path(item.split("=", 1)[1]) for item in items}  # noqa: E731
    built = build_release(arguments.version, arguments.output_root, pairs(arguments.graph), pairs(arguments.evidence), arguments.data_sources, Path("."), arguments.notes, pairs(arguments.extra))
    problems = verify_release(built)
    print(f"wrote {built}; files {len(list(built.rglob('*')))}; verification: {'ok' if not problems else problems}")


if __name__ == "__main__":
    main()
