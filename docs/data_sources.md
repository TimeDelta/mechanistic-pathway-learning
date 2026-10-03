# Data sources, versions and licenses

Pin the exact release of every source in this table before Phase 1 counts are run; the pinned
versions are part of the pre-registration. Raw downloads live under data/raw/ (gitignored).

| Source | Use | Access | License note | Pinned version |
|---|---|---|---|---|
| Human-GEM (github.com/SysBioChalmers/Human-GEM) | base reconstruction (SBML for stoichiometry and gene rules; yml for reaction subsystems) | public | open | 2.0.1 (model date 2026-03-26; main fetched 2026-10-02 and 2026-10-03) |
| Recon3D (vmh.life; BiGG) | identifier cross-reference | public | open for academic use | to pin |
| Human Phenotype Ontology (hpo.jax.org) | E1 monogenic gene-phenotype annotations with frequency qualifiers and OMIM or Orphanet provenance; hp.obo for the is_a closure | public | open | release 2026-09-01 (genes_to_phenotype.txt, 333,983 rows) |
| OMIM (omim.org) | clinical synopses for the grade A bar | registration | academic use | to pin |
| SIDER 4.1 (sideeffects.embl.de, HTTPS) | E2 label events, 2015 time slice | public | CC BY-NC-SA: derived tables are non-commercial and share-alike (docs/alternative_sources.md, licensing posture); MedDRA terms inside are licensed; do not redistribute term tables | 4.1 (files dated 2015-10-21; fetched 2026-10-02 and 2026-10-03) |
| OnSIDES (github.com/tatonetti-lab/onsides) | E2 label events, current slice | public | same MedDRA note | to pin |
| ChEMBL (ebi.ac.uk/chembl) | drug mechanisms, targets, pChEMBL | public | CC BY-SA | ChEMBL_37 (2026-05-01) |
| UniChem (ebi.ac.uk/unichem) | PubChem to ChEMBL identifier mapping for SIDER drugs | public API | open | queried 2026-10-03; 1,177 of 1,430 SIDER drugs mapped; 7,561 ChEMBL mechanisms; 337 targets with gene symbols |
| OmniPath (omnipathdb.org) | signaling edges with signs (datasets omnipath, pathwayextra, ligrecextra) and small_molecule_protein interactions | public web service | mixed by resource; check per-resource licenses before redistribution | web service queried 2026-10-02 (first build) and 2026-10-03 (this build); 126,124 interaction rows, 6,553 small-molecule rows |
| GTEx v10 median gene TPM (storage.googleapis.com/adult-gtex, public bucket) | brain expression attribute on gene and reaction nodes; regulon restriction to brain-expressed transcription factors (design 3.3, 4.1) | public | GTEx open-access data; data-use terms of the GTEx portal apply | GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz, object date 2024-11-08, 59,033 genes, 13 brain tissues; threshold median TPM 1.0 |
| Orphadata product 6 (orphadata.com, en_product6.xml) | gene-disease association types for Orphanet entries: disease-causing versus susceptibility, candidate, modifier, biomarker or fusion (grade A rule, section 4.3) | public | CC BY 4.0 | JDBOR date 2026-06-23; 8,503 associations; downloaded 2026-10-03 |
| OnSIDES v3.1.1 (github.com/tatonetti-lab/onsides, release asset) | E2 later label slice: current-label adverse-event statements, one report per deduplicated statement | public | CC BY 4.0 (OnSIDES); MedDRA term names inside are licensed, not redistributed | data release 3.1.1-20260422; 6,928,666 mentions, 41,119 labels, 1,866 ingredients; bridge: RxNav and UniChem queried 2026-10-03; 1277 qualifying statements |
| RxNav (rxnav.nlm.nih.gov) | RxNorm ingredient to UNII and ATC classes | public API | open | queried 2026-10-03 for 1,866 ingredients |
| CollecTRI (via OmniPath, dataset collectri) | transcription factor regulons with signs | public web service | CC BY 4.0 (CollecTRI); OmniPath terms apply to the service | queried 2026-10-03; 64,515 regulon rows |
| Human Protein Atlas; GTEx | brain expression weights | public | open | to pin |
| PubTator3 (ncbi.nlm.nih.gov/research/pubtator3) | E3 literature relations | public API | open | query date |
| CTD (ctdbase.org) | E3 chemical-symptom associations with direction | public | academic use; cite | to pin |
| SemMedDB (NLM) | E3 predications with predicate types | UMLS license required | UMLS terms of use | optional |
| GWAS summary statistics (GWAS Catalog) | E4 symptom-level genetics for validation | public | per study | to pin |
| MAGMA (cncr.nl/research/magma) | gene-set enrichment | public binary | academic use | to pin |

Not used: the earlier in-house psychiatric literature knowledge graph (quality judged insufficient);
individual-level cohort data (no new data collection in version 1).
