"""Physiology graph construction (design section 4.1).

Layers and sources:
  metabolic      Human-GEM (SBML/YAML from github.com/SysBioChalmers/Human-GEM) via cobra:
                 metabolite nodes (compartment-specific), reaction nodes, gene nodes,
                 edges substrate-of, product-of, catalyzed-by (from GPR rules),
                 transports-between (reactions spanning compartments)
  signaling      OmniPath (omnipath client): binds, activates, inhibits with signs
  transcription  CollecTRI (via decoupler or OmniPath): regulates-transcription-of, signed
  pharmacology   ChEMBL mechanisms and activities: targets edges with action type and pChEMBL
  attributes     brain expression weight per gene (HPA/GTEx), BBB penetration flag per drug,
                 currency tag per metabolite (tag_currency_metabolites)

Hard structural constraints are enforced here, not in the model: no other edge
types, no symptom-symptom edges, no disease nodes, compartment changes only via
transport reactions, pharmacological edges only with measured activity.

Outputs (data/processed/graph/):
  nodes.parquet   node_id, node_type, compartment, display_name, brain_expression_weight, is_currency
  edges.parquet   source_id, target_id, relation_type, sign, evidence_source, pchembl
  relation_types.json  ordered list used as the encoder's relation index

Not implemented in version 0.1.
"""
from __future__ import annotations

from pathlib import Path

RELATION_TYPES: list[str] = [
    "substrate_of",
    "product_of",
    "catalyzed_by",
    "transports_between",
    "binds",
    "activates",
    "inhibits",
    "regulates_transcription_of",
    "targets",
]


def build_physiology_graph(human_gem_path: Path, output_directory: Path, include_signaling: bool = True, include_transcription: bool = True):
    raise NotImplementedError("load Human-GEM with cobra; add OmniPath and CollecTRI layers; write nodes.parquet and edges.parquet")
