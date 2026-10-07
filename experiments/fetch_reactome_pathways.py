"""Fetch and pin the Reactome pathways the neuronal and oxidative graph variant imports
(mechanistic_pathway_learning/graph/reactome_import.py; design section 4.1).

Each pathway is exported as SBML by the Reactome content service (reactions with their reactants, products, catalysts
and positive or negative regulators, every species annotated with a ChEBI or UniProt identifier and a compartment)
and written to data/raw/reactome/v<release>/<stable id>.sbml with a pin beside it (URL, bytes, SHA-256, fetch time,
database release). Files already pinned are not fetched again unless --refresh is given.

Usage:
  python experiments/fetch_reactome_pathways.py
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fetch_node_descriptor_sources import download_and_pin  # noqa: E402

from mechanistic_pathway_learning.evidence.cached_http import USER_AGENT

PATHWAYS = {
    "R-HSA-112316": "Neuronal System (neurotransmitter release cycle, receptors, potassium channels, ionotropic and metabotropic transmission)",
    "R-HSA-983712": "Ion channel transport (stimuli-sensing channels, P-type ATPases, ligand-gated and other channels)",
    "R-HSA-426117": "Cation-coupled chloride cotransporters (NKCC1, KCC2 and relatives)",
    "R-HSA-8949215": "Mitochondrial calcium ion transport",
    "R-HSA-500792": "GPCR ligand binding",
    "R-HSA-418594": "G alpha (i) signalling events",
    "R-HSA-418555": "G alpha (s) signalling events",
    "R-HSA-416476": "G alpha (q) signalling events",
    "R-HSA-397795": "G-protein beta:gamma signalling",
    "R-HSA-3299685": "Detoxification of reactive oxygen species",
    "R-HSA-9755511": "KEAP1-NFE2L2 pathway",
    "R-HSA-9711123": "Cellular response to chemical stress",
    "R-HSA-202131": "Metabolism of nitric oxide: NOS3 activation and regulation",
    "R-HSA-1222556": "ROS and RNS production in phagocytes",
    "R-HSA-611105": "Respiratory electron transport",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    parser.add_argument("--refresh", action="store_true")
    arguments = parser.parse_args()
    request = urllib.request.Request("https://reactome.org/ContentService/data/database/version", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        release = response.read().decode().strip()
    source = f"reactome/v{release}"
    pins = {}
    for stable_id, description in PATHWAYS.items():
        pin = download_and_pin(arguments.raw_root, source, f"{stable_id}.sbml", f"https://reactome.org/ContentService/exporter/sbml/{stable_id}.xml", arguments.refresh)
        pins[stable_id] = {"description": description, "bytes": pin["bytes"], "sha256": pin["sha256"], "fetched_at": pin["fetched_at"]}
        print(f"{stable_id} ({description}): {pin['bytes']:,} bytes, sha256 {pin['sha256'][:12]}")
    (arguments.raw_root / source / "pathways.json").write_text(json.dumps({"reactome_release": release, "pathways": pins}, indent=1) + "\n")


if __name__ == "__main__":
    main()
