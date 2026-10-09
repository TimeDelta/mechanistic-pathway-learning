"""Fetch and pin the UniProt annotations of the plasma substrate binders (docs/plasma_protein_binding.md, the user's
decision of 9 October 2026: "All of the plasma substrate binders should be modeled").

One file, data/raw/uniprot/uniprot_plasma_binders.tsv.gz, with a <file>.pin.json beside it: the reviewed human entries
of the binder genes with their binding-site features (ft_binding, whose ligand entries carry ChEBI identifiers that
mechanistic_pathway_learning/graph/plasma_binding.py maps to Human-GEM metabolites), their function text (cc_function),
their subcellular location and their primary gene name. The pinned proteome file
(uniprot_human_reviewed.tsv.gz) carries neither binding sites nor function text, which is why this is a separate fetch.

The query names the genes rather than a keyword, so the pinned file is exactly the binder set the build reads and a
later UniProt release cannot silently widen it. Raw files are never committed (data/raw is gitignored); the pin is
recorded in docs/data_sources.md.

Usage:
  python experiments/fetch_plasma_binder_annotations.py [--refresh]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.parse import quote

from experiments.fetch_node_descriptor_sources import download_and_pin

# The binders of docs/plasma_protein_binding.md. ORM1 and ORM2 are the two genes of orosomucoid; RBP4 carries retinol
# and circulates bound to transthyretin; APOA1 and APOB are the structural apolipoproteins of HDL and of VLDL/LDL.
BINDER_GENES = ("ALB", "ORM1", "ORM2", "SERPINA6", "SERPINA7", "TTR", "AFP", "SHBG", "RBP4", "APOA1", "APOB")
UNIPROT_FIELDS = "accession,gene_primary,gene_names,protein_name,ft_binding,cc_function,cc_subcellular_location"


def binder_query(genes=BINDER_GENES) -> str:
    """The UniProt query string: the reviewed human entries whose gene name is one of the binders."""
    return f"(organism_id:9606) AND (reviewed:true) AND ({' OR '.join(f'gene_exact:{gene}' for gene in genes)})"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    parser.add_argument("--refresh", action="store_true", help="fetch again even when a pinned file is present")
    arguments = parser.parse_args()
    url = (f"https://rest.uniprot.org/uniprotkb/stream?query={quote(binder_query())}"
           f"&fields={UNIPROT_FIELDS}&format=tsv&compressed=true")
    pin = download_and_pin(arguments.raw_root, "uniprot", "uniprot_plasma_binders.tsv.gz", url, arguments.refresh)
    print(json.dumps({key: value for key, value in pin.items() if key != "url"}, indent=2))
    print(f"genes asked for: {len(BINDER_GENES)} ({', '.join(BINDER_GENES)})")


if __name__ == "__main__":
    main()
