"""Interpretability deliverables (design section 6.6).

For each symptom: ranked consensus modules with support subgraph, compartment
annotations, supporting perturbations by grade and selection frequency across
seeds and folds; plus an overlap table against docs/curated_pathway_modules.csv.

Not implemented in version 0.1.
"""


def build_mechanism_card(symptom_identifier, consensus_modules, graph_nodes_table, evidence_records):
    raise NotImplementedError("render a markdown card per symptom under runs/<run>/mechanism_cards/")
