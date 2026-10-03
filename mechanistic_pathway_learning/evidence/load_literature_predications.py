"""Evidence class E3: literature predications as a soft prior only (design section 4.2).

The earlier in-house psychiatric literature knowledge graph is not used. Sources, in
order of preference:

1. PubTator3 (NCBI, public, no license): entity annotations (gene, chemical, disease
   including MeSH sign/symptom descriptors) and extracted relations such as
   associate, cause, treat, inhibit, stimulate, positive_correlate, negative_correlate.
   Query by target symptom MeSH descriptor and by graph gene symbols; keep PMID,
   relation type and year for the time split.
2. CTD (ctdbase.org, public, curated): chemical-gene and chemical-disease associations
   with direction (therapeutic versus marker/mechanism). Chemical-disease rows whose
   disease is a MeSH symptom descriptor supply induce/relieve priors for chemicals
   with known targets.
3. SemMedDB (NLM; requires a UMLS license, free for academic use): predications with
   predicate types CAUSES, AFFECTS, DISRUPTS, PREDISPOSES, ASSOCIATED_WITH between
   graph entities (by CUI) and symptom CUIs. Richer predicates than PubTator3; use
   when the license is in place.

Output: EvidenceRecord rows with evidence_class "literature", predication_type set and
evidence_available_date set to the publication date, so weights follow
assign_evidence_grades.DEFAULT_PREDICATION_WEIGHTS and the time split can exclude
late predications from training.

Not implemented in version 0.1.
"""
from __future__ import annotations


def fetch_pubtator3_relations(symptom_mesh_descriptors: list[str], gene_symbols: list[str]):
    raise NotImplementedError("use the PubTator3 API; cache raw responses under data/raw/pubtator3/")


def load_ctd_chemical_symptom_associations(ctd_directory):
    raise NotImplementedError("parse CTD_chemicals_diseases.tsv.gz; keep DirectEvidence and PubMedIDs")


def load_semmeddb_predications(semmeddb_directory, symptom_cuis: list[str], entity_cuis: list[str]):
    raise NotImplementedError("requires UMLS license; filter PREDICATION table by predicate and CUI lists")
