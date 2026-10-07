"""Fetch and pin the sources of the fixed node descriptors and of the HPO laboratory-abnormality labels
(design section 5.2, node descriptors; section 4, auxiliary labels).

Files, each written to data/raw/<source>/ with a <file>.pin.json beside it (URL, bytes, SHA-256, Last-Modified or
release header, fetch time); a file already present with a pin is not fetched again unless --refresh is given:
  - uniprot/uniprot_human_reviewed.tsv.gz: the reviewed human proteome from the UniProt REST stream (accession,
    primary gene name, HGNC, EC numbers, subcellular location, Pfam, length, sequence); sequences feed the ESM-2
    embeddings, the annotations the reduced-rank regression that reduces them;
  - gene_ontology/goa_human.gaf.gz and gene_ontology/go-basic.obo: GO annotations with evidence codes (phenotype-derived
    codes are dropped downstream) and the ontology to propagate them to ancestor terms;
  - hpo/hp-base.owl: the HPO release whose logical definitions name the ChEBI entity and the direction of each
    laboratory abnormality (hp.obo carries none), pinned to the release of data/raw/hpo/hp.obo;
  - chebi/chebi_core.obo.gz: ChEBI with its conjugate acid and base relations, to bridge the neutral forms HPO names
    and the charged forms Human-GEM uses.
Raw files are never committed (data/raw is gitignored); the pins are recorded in docs/data_sources.md.

Usage:
  python experiments/fetch_node_descriptor_sources.py [--refresh]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path

from mechanistic_pathway_learning.evidence.cached_http import USER_AGENT

HPO_RELEASE = "v2026-09-01"  # the release of data/raw/hpo/hp.obo (data-version hp/releases/2026-09-01)
UNIPROT_FIELDS = "accession,gene_primary,xref_hgnc,ec,cc_subcellular_location,xref_pfam,length,sequence"
SOURCES = [
    ("uniprot", "uniprot_human_reviewed.tsv.gz",
     f"https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)%20AND%20(reviewed:true)&fields={UNIPROT_FIELDS}&format=tsv&compressed=true"),
    ("gene_ontology", "goa_human.gaf.gz", "https://current.geneontology.org/annotations/goa_human.gaf.gz"),
    ("gene_ontology", "go-basic.obo", "https://current.geneontology.org/ontology/go-basic.obo"),
    ("hpo", "hp-base.owl", f"https://github.com/obophenotype/human-phenotype-ontology/releases/download/{HPO_RELEASE}/hp-base.owl"),
    ("chebi", "chebi_core.obo.gz", "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_core.obo.gz"),
]
PINNED_HEADERS = ("last-modified", "content-length", "x-uniprot-release", "x-uniprot-release-date", "x-total-results", "etag")


def sha256_of_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_and_pin(raw_root: Path, source: str, file_name: str, url: str, refresh: bool) -> dict:
    directory = raw_root / source
    directory.mkdir(parents=True, exist_ok=True)
    target_path, pin_path = directory / file_name, directory / f"{file_name}.pin.json"
    if target_path.exists() and pin_path.exists() and not refresh:
        return json.loads(pin_path.read_text())
    partial_path = directory / f"{file_name}.part"
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=300) as response, open(partial_path, "wb") as handle:
        headers = {key.lower(): value for key, value in response.headers.items()}
        shutil.copyfileobj(response, handle, length=1 << 20)
    partial_path.replace(target_path)
    pin = {"url": url, "file": str(target_path), "bytes": target_path.stat().st_size, "sha256": sha256_of_file(target_path),
           "fetched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           **{header.replace("-", "_"): headers[header] for header in PINNED_HEADERS if header in headers}}
    pin_path.write_text(json.dumps(pin, indent=2) + "\n")
    return pin


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    parser.add_argument("--refresh", action="store_true", help="fetch again even when a pinned file is present")
    arguments = parser.parse_args()
    for source, file_name, url in SOURCES:
        pin = download_and_pin(arguments.raw_root, source, file_name, url, arguments.refresh)
        print(f"{source}/{file_name}: {pin['bytes']:,} bytes, sha256 {pin['sha256'][:12]}, "
              f"{pin.get('last_modified') or pin.get('x_uniprot_release') or 'no release header'}")


if __name__ == "__main__":
    main()
