# Data sources, versions and licenses

Pin the exact release of every source in this table before Phase 1 counts are run; the pinned
versions are part of the pre-registration. Raw downloads live under data/raw/ (gitignored).

| Source | Use | Access | License note | Pinned version |
|---|---|---|---|---|
| Human-GEM (github.com/SysBioChalmers/Human-GEM) | base reconstruction | public | open | 2.0.1 (main, 2026-10-02) |
| Recon3D (vmh.life; BiGG) | identifier cross-reference | public | open for academic use | to pin |
| Human Phenotype Ontology (hpo.jax.org) | E1 monogenic gene-phenotype annotations; hp.obo | public | open | release 2026-09-01 |
| OMIM (omim.org) | clinical synopses for the grade A bar | registration | academic use | to pin |
| SIDER 4.1 (sideeffects.embl.de) | E2 label events, 2015 time slice | public | MedDRA terms inside are licensed; do not redistribute term tables | 4.1 (files dated 2015-10-21) |
| OnSIDES (github.com/tatonetti-lab/onsides) | E2 label events, current slice | public | same MedDRA note | to pin |
| ChEMBL (ebi.ac.uk/chembl) | drug mechanisms, targets, pChEMBL | public | CC BY-SA | ChEMBL_37 (2026-05-01) |
| UniChem (ebi.ac.uk/unichem) | PubChem to ChEMBL identifier mapping for SIDER drugs | public API | open | queried 2026-10-02 |
| OmniPath (omnipathdb.org) | signaling edges with signs | public | mixed by resource; check per-resource licenses before redistribution | to pin |
| CollecTRI (via decoupler or OmniPath) | transcription factor regulons | public | open | to pin |
| Human Protein Atlas; GTEx | brain expression weights | public | open | to pin |
| PubTator3 (ncbi.nlm.nih.gov/research/pubtator3) | E3 literature relations | public API | open | query date |
| CTD (ctdbase.org) | E3 chemical-symptom associations with direction | public | academic use; cite | to pin |
| SemMedDB (NLM) | E3 predications with predicate types | UMLS license required | UMLS terms of use | optional |
| GWAS summary statistics (GWAS Catalog) | E4 symptom-level genetics for validation | public | per study | to pin |
| MAGMA (cncr.nl/research/magma) | gene-set enrichment | public binary | academic use | to pin |

Not used: the earlier in-house psychiatric literature knowledge graph (quality judged insufficient);
individual-level cohort data (no new data collection in version 1).
