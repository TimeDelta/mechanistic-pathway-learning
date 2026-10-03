# Literature survey, evidence class E3: retrieval half (PubTator3 and CTD)

Status: retrieval implemented and run on 3 October 2026; the join into the assembled evidence table and the
reliability fit is the next step. Design sections 3.1, 4.2 (E3 paragraph), 4.3 and 6.1; report schema in
docs/evidence_reports_spec.md section 1.

Literature reports are soft priors only (grades C to E). They never become evaluation positives, a pair that no
paper mentions is unlabelled rather than negative, and every row keeps the exact extraction model behind it
(PubTator3 relation type, CTD direct-evidence value, descriptor level, chemical match method) so that each report can be
weighted by the model that produced it and by its limitations, not by its count. This is a first pass: chemicals are
matched to SIDER by name, genes by NCBI Gene id, and the only study-design information is the relation type; the
document appraisal (llm_evidence_appraisal.py) is the step that reads the papers.

## 1. What was built

| File | What it does |
|---|---|
| docs/symptom_crosswalk.csv, columns mesh_descriptors and mesh_descriptor_level (appended) | MeSH descriptor ids per target symptom, each flagged "symptom" or "diagnosis" |
| mechanistic_pathway_learning/evidence/mesh_symptom_descriptors.py | reads and validates the two columns (aligned lists, D-numbers, one symptom per descriptor) |
| experiments/verify_mesh_symptom_descriptors.py | verifies every id against the NLM MeSH RDF service and the PubTator3 autocomplete endpoint (section 2) |
| mechanistic_pathway_learning/evidence/gene_identifier_map.py | symbol, NCBI Gene id and Ensembl id for every gene node of the full graph (section 3) |
| mechanistic_pathway_learning/evidence/sider_chemical_matching.py | SIDER drug name to MeSH chemical id by name, with the match method recorded (section 4) |
| mechanistic_pathway_learning/evidence/cached_http.py | rate-limited (3 requests per second), retried, disk-cached client for NCBI and NLM |
| mechanistic_pathway_learning/evidence/pubmed_publication_dates.py | esummary dates for PMIDs, 200 per request, cached, for the time split |
| mechanistic_pathway_learning/evidence/load_pubtator_relations.py | bulk-file row parsing, relation-type mapping, partner filters, report construction |
| experiments/fetch_pubtator_relations.py | the PubTator3 fetch (section 5) |
| mechanistic_pathway_learning/evidence/load_ctd_relations.py | CTD direct-evidence mapping and expansion to one report per cited paper |
| experiments/fetch_ctd_chemical_disease.py | the CTD fetch (section 6) |
| tests/test_literature_relations.py | fixtures for relation-type mapping, descriptor-level flags, gene id precedence, CTD evidence mapping, name matching, schema |

The stub load_literature_predications.py stays in place; the two loaders above supersede it. SemMedDB is still
gated on the UMLS license.

Outputs (gitignored): data/processed/literature/pubtator_relations_filtered.parquet,
ctd_relations_filtered.parquet, gene_identifier_map.parquet, sider_chemical_mesh_map.parquet,
pubmed_publication_dates.json, pubtator_relations_summary.json, ctd_relations_summary.json; raw downloads and
response caches under data/raw/pubtator3, data/raw/ctd and data/raw/hgnc with a pin file next to each download.

## 2. MeSH descriptors per target symptom and their verification

Choice. Each target symptom gets the MeSH descriptor that names the sign or symptom itself (level "symptom") and,
flagged apart, the diagnostic category whose literature is a soft prior on its cardinal symptom (level "diagnosis").
Both levels are retrieved; the level column travels with every report so the two can be weighted apart, and the
limitations text names a diagnosis-level report as such. Assumption A7 keeps diagnoses out of the targets, not out
of the priors. Diagnosis descriptors were limited to the category nearest the symptom (Depressive Disorder and
Major Depressive Disorder, Anxiety Disorders, Psychotic Disorders, Bipolar Disorder, Dementia, Disorders of Excessive
Somnolence); Schizophrenia, the chronic fatigue syndrome and the named dementias were not added, because their
literature is about the disease and not the symptom.

Verification, run on 3 October 2026 (experiments/verify_mesh_symptom_descriptors.py, responses cached under
data/raw/pubtator3/mesh_lookup and data/raw/pubtator3/autocomplete, result in
data/raw/pubtator3/mesh_descriptor_verification.json). "MeSH valid" means the NLM MeSH RDF service returns the
descriptor with its preferred label and tree numbers; "PubTator3 entity" is the entity the autocomplete endpoint
returns for that label with the same db_id; "tree opinion" is what the tree numbers alone say (C23.888 Signs and
Symptoms, F01 Behavior and C10.597 Neurologic Manifestations for symptom; F03 Mental Disorders, C10.228 and C10.886
for diagnosis).

| Target symptom | MeSH id | Level | MeSH preferred label | Tree numbers | Tree opinion | MeSH valid | PubTator3 entity |
|---|---|---|---|---|---|---|---|
| depressed_mood | D003863 | symptom | Depression | F01.145.126.350; F01.470.282 | symptom | yes | none (not in the PubTator3 disease vocabulary) |
| depressed_mood | D003866 | diagnosis | Depressive Disorder | F03.600.300 | diagnosis | yes | @DISEASE_Depressive_Disorder |
| depressed_mood | D003865 | diagnosis | Major Depressive Disorder | F03.600.300.375 | diagnosis | yes | @DISEASE_Depressive_Disorder_Major |
| anhedonia | D059445 | symptom | Anhedonia | C10.597.606.057; C23.888.592.604.039; F01.700.039 | symptom | yes | @DISEASE_Anhedonia |
| anxiety | D001007 | symptom | Anxiety | F01.470.132 | symptom | yes | @DISEASE_Anxiety |
| anxiety | D001008 | diagnosis | Anxiety Disorders | F03.080 | diagnosis | yes | @DISEASE_Anxiety_Disorders |
| irritability_or_aggression | D007508 | symptom | Irritable Mood | F01.470.047.110 | symptom | yes | none (not in the PubTator3 disease vocabulary) |
| irritability_or_aggression | D000374 | symptom | Aggression | F01.145.126.050.500; F01.145.126.125; F01.145.813.045 | symptom | yes | none (not in the PubTator3 disease vocabulary) |
| insomnia | D007319 | symptom | Sleep Initiation and Maintenance Disorders | C10.886.425.800.800; F03.870.400.800.800 | diagnosis | yes | @DISEASE_Sleep_Initiation_and_Maintenance_Disorders |
| somnolence_or_hypersomnia | D000077260 | symptom | Sleepiness | C23.888.900; G11.561.803.754.917 | symptom | yes | @DISEASE_Sleepiness |
| somnolence_or_hypersomnia | D006970 | diagnosis | Disorders of Excessive Somnolence | C10.886.425.800.200; F03.870.400.800.200 | diagnosis | yes | @DISEASE_Disorders_of_Excessive_Somnolence |
| fatigue | D005221 | symptom | Fatigue | C23.888.369 | symptom | yes | @DISEASE_Fatigue |
| fatigue | D005222 | symptom | Mental Fatigue | C23.888.369.500; F01.145.126.937 | symptom | yes | @DISEASE_Mental_Fatigue |
| psychomotor_agitation | D011595 | symptom | Psychomotor Agitation | C10.597.350.600; C10.597.606.881.700; C23.888.592.350.600; C23.888.592.604.882.700; F01.145.126.050.625; F01.700.875.700 | symptom | yes | @DISEASE_Psychomotor_Agitation |
| psychomotor_retardation | none | | | | | | no MeSH descriptor names psychomotor retardation; Psychomotor Disorders D011596 is a broader category and is not used (crosswalk notes) |
| psychosis | D006212 | symptom | Hallucinations | C10.597.606.762.300; C23.888.592.604.764.300; F01.700.750.300 | symptom | yes | @DISEASE_Hallucinations |
| psychosis | D003702 | symptom | Delusions | F01.145.126.200 | symptom | yes | none (not in the PubTator3 disease vocabulary) |
| psychosis | D011618 | diagnosis | Psychotic Disorders | F03.700.675 | diagnosis | yes | @DISEASE_Psychotic_Disorders |
| cognitive_impairment | D060825 | symptom | Cognitive Dysfunction | F03.615.250.700 | diagnosis | yes | @DISEASE_Cognitive_Dysfunction |
| cognitive_impairment | D008569 | symptom | Memory Disorders | C10.597.606.525; C23.888.592.604.529; F01.700.625 | symptom | yes | @DISEASE_Memory_Disorders |
| cognitive_impairment | D003221 | symptom | Confusion | C10.597.606.337; C23.888.592.604.339; F01.700.250 | symptom | yes | @DISEASE_Confusion |
| cognitive_impairment | D003704 | diagnosis | Dementia | C10.228.140.380; F03.615.400 | diagnosis | yes | @DISEASE_Dementia |
| elevated_mood_or_mania | D000087122 | symptom | Mania | C10.597.606.483; C23.888.592.604.487; F01.700.548 | symptom | yes | @DISEASE_Mania |
| elevated_mood_or_mania | D001714 | diagnosis | Bipolar Disorder | F03.600.150.500 | diagnosis | yes | @DISEASE_Bipolar_Disorder |

Findings of the verification.

- Mania has the current descriptor id D000087122 (the autocomplete returns @DISEASE_Mania with that id); Sleepiness
  exists as D000077260 and is used as the symptom-level descriptor for somnolence, with Disorders of Excessive
  Somnolence D006970 flagged as a diagnosis because its children are named disorders (narcolepsy, Kleine-Levin).
- Four symptom-level descriptors are valid MeSH ids that PubTator3 does not use: Depression D003863, Irritable Mood
  D007508, Aggression D000374 and Delusions D003702. They sit only in the Behavior branch F01, and PubTator3
  normalizes disease mentions to the MEDIC vocabulary (MeSH branches C and F03 plus a few F01 behaviors such as Anxiety),
  so the autocomplete returns no entity, the relations endpoint returns an empty list and the bulk file has zero
  rows for them (section 5). CTD uses the same MEDIC vocabulary. They stay in the crosswalk as the correct MeSH ids
  for a source that uses full MeSH (SemMedDB, MeSH-indexed PubMed queries), and the consequence is recorded: in
  PubTator3 and CTD, depressed mood rests on diagnosis-level descriptors only, irritability or aggression has no
  channel and psychosis rests on Hallucinations and Psychotic Disorders.
- Two level flags deviate from the tree numbers on purpose. Sleep Initiation and Maintenance Disorders D007319 is
  flagged "symptom": MeSH has no sign-or-symptom descriptor for insomnia, "Insomnia" and "Sleeplessness" are entry terms
  of D007319, and it is the id every insomnia mention normalizes to; it has no child disorders. Cognitive Dysfunction
  D060825 is flagged "symptom" because the descriptor names the impairment itself (entry terms Cognitive Impairment,
  Cognitive Decline, Mild Cognitive Impairment) and is the MeSH counterpart of HPO Cognitive impairment, even though
  MeSH trees it under Neurocognitive Disorders; the mild-cognitive-impairment diagnosis folded into the same
  descriptor is a limitation stated here rather than hidden.
- Psychomotor retardation has no MeSH descriptor; the crosswalk notes say so and the symptom has no E3 channel.

## 3. Gene identifier map

Human-GEM model/genes.tsv supplies the NCBI id (geneEntrezID) and Ensembl id for the metabolic genes; the HGNC
complete set (pinned in section 7) supplies the rest by current symbol, then by a previous symbol carried by
exactly one gene, then by an alias carried by exactly one gene. Human-GEM comes first because the metabolic layer
was built from it; an ambiguous previous or alias symbol is left unmapped rather than guessed.

<!-- GENE_MAP_COUNTS -->

## 4. SIDER drug to MeSH chemical matching (first pass, by name)

A SIDER drug matches a MeSH chemical when the MeSH descriptor name, or the entry term the PubTator3 autocomplete
reports as the matched synonym, equals the SIDER name case-insensitively (methods autocomplete_name and
autocomplete_synonym, one cached query per drug name); chemical ids on kept rows that no SIDER name matched are then
looked up by preferred label at id.nlm.nih.gov and matched again by name (method mesh_label_lookup). CTD rows match
by ChemicalName first (method ctd_chemical_name) and otherwise by a MeSH id the PubTator3 pass already matched. The
method is a column of every drug report. Salts, combination products and spelling differences ("gamma-aminobutyric"
against "gamma-Aminobutyric Acid") do not match in this pass and are counted, not guessed; a structure-based mapping
(PubChem CID to MeSH through UniChem or MeSH registry numbers) is the second pass.

<!-- CHEMICAL_MATCH_COUNTS -->

## 5. PubTator3 bulk relations

The bulk file relation2pubtator3.gz was 297,478,945 bytes (under the 400 MB download cap), so it was downloaded
to data/raw/pubtator3 and streamed once. Columns: PMID, relation type, entity 1, entity 2, entities written
"Type|Identifier" (Disease|MESH:D001007, Gene|8813, Chemical|MESH:D005473, Variant and Species forms); no date.
Rows are kept when one entity is a crosswalk descriptor and the other is a graph gene (NCBI id in the map of
section 3) or a matched SIDER drug. Relation types map as follows; nothing else is kept and every drop is counted.

| PubTator3 relation type | Repository relation | report_value | Note |
|---|---|---|---|
| cause, positive_correlate | induces | 1 | |
| treat, prevent, negative_correlate | relieves | 1 | |
| associate | associated_with | 1 | kept as its own undirected relation, not merged into induces |
| stimulate, inhibit | dropped | | chemical-gene relations; counted |
| cotreat, compare, interact, drug_interact | dropped | | not perturbation-symptom relations; counted |

Report columns follow docs/evidence_reports_spec.md section 1 (report_id
"PubTator3|PMID:{pmid}|MESH:{descriptor}|{symptom}|{ordinal}", source "PubTator3", evidence_class "literature",
evidence_code "pubtator3_{relation type}", references "PMID:{pmid}", rubric one-hots 0, perturbation_nodes the
gene node with the loss-of-function sign as a documented convention or an empty list for drugs, which the assembler
joins to ChEMBL targets by STITCH flat id as it does for SIDER rows) plus the literature columns pmid, relation_type,
entity_role (disease_first or disease_second), perturbation_ncbi_gene_id, chemical_mesh_id, chemical_match_method,
mesh_descriptor, mesh_descriptor_level, publication_year and direct_evidence (CTD only). Publication dates come from
esummary (pubdate parsed to an ISO date, first of the month or year when the finer part is missing).

<!-- PUBTATOR_COUNTS -->

## 6. CTD chemical-disease statements

CTD_chemicals_diseases.tsv.gz rows are kept when DiseaseID is a crosswalk descriptor and the chemical matches a SIDER
drug (section 4). DirectEvidence marker/mechanism maps to induces and therapeutic to relieves; a row carrying both
yields two reports; rows with empty DirectEvidence are inferred through a gene and are dropped with a count. Each
PubMed id of a row is its own report (report_id "CTD|MESH:{chemical}|MESH:{descriptor}|{symptom}|{ordinal}",
evidence_code "ctd_direct_evidence_{marker_mechanism|therapeutic}"); a direct-evidence row that cites no paper is one
report with no reference. Marker/mechanism covers correlation as well as causation, which the limitations text says.

<!-- CTD_COUNTS -->

## 7. Pins

<!-- PINS -->

## 8. What the next step needs

- Joining: the two parquet files have the evidence_reports.parquet columns in order followed by the literature
  columns, so the assembler can concatenate them after filling perturbation_nodes for drugs from the ChEMBL mapping
  and applying the same inclusion filters as for SIDER drugs (dominant target, ATC N). associated_with is a third
  relation; experiment_data.py filters on relation and will ignore it until it is told what to do with it.
- Weights: grade E with the predication weights of assign_evidence_grades (CAUSES-like 0.10 for cause and
  positive_correlate, ASSOCIATED_WITH 0.05 for associate; treat and negative_correlate need their own entries) until
  the reliability fit has a third source; the descriptor level and the chemical match method are the first rubric
  features to add.
- Time split: evidence_date is set for every dated PMID; the share of undated rows is in section 5.
- Second pass on chemicals (structure-based mapping) and on genes (species check through the Species entity of each
  PMID, since a human NCBI Gene id can be mentioned in an animal study).
