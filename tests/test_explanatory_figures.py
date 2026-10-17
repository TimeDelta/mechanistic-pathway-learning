"""The edge layers the explanatory figures report (experiments/make_explanatory_figures.py, figures 1 and 10): each
source names one layer, the first rule that matches wins, and an electrical edge whose source names Reactome at its
end is electrical, not Reactome."""
from __future__ import annotations

import pandas as pd

from experiments.make_explanatory_figures import OTHER_CURATED_LAYER, edge_layers


def test_each_edge_source_names_one_layer():
    sources = {
        "Human-GEM": "metabolic (Human-GEM)", "Human-GEM GPR": "metabolic (Human-GEM)", "OmniPath": "signalling (OmniPath)",
        "OmniPath small_molecule_protein": "signalling (OmniPath)", "CollecTRI": "transcription (CollecTRI)",
        "Reactome entity member": "Reactome complexes, sets and regulators", "Reactome regulator": "Reactome complexes, sets and regulators",
        "BRENDA CC BY 4.0 (1 references, accession stated by BRENDA)": "enzyme regulation (BRENDA)",
        "net charge across the plasma membrane (Reactome)": "electrical (ions and charge)",
        "ion moved across the plasma membrane (Human-GEM)": "electrical (ions and charge)",
        "drug mechanism (ChEMBL), moved from the perturbation's seeds": "drug mechanism (ChEMBL)",
        "plasma carrier of the drug (both; label)": "drug plasma sequestration",
        "curated plasma carriage": "plasma carriage of metabolites",
        "neomorphic allele mechanism (PMID 11264541)": "neomorphic allele mechanism",
        "glutathione consumed": OTHER_CURATED_LAYER,
    }
    edges = pd.DataFrame({"evidence_source": list(sources), "relation_type": "activates"})
    assert list(edge_layers(edges)) == list(sources.values())
    encodes = pd.DataFrame({"evidence_source": ["gene_protein_split", None], "relation_type": ["encodes", "encodes"]})
    assert list(edge_layers(encodes)) == ["gene to protein (encodes)", "gene to protein (encodes)"], "the relation decides, whatever the source"
