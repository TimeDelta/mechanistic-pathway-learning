# Literature survey, evidence class E3: retrieval and appraisal (PubTator3 and CTD)

Status: both halves implemented and first run on 3 October 2026. The adversarial review of the component
(docs/design_review_v04.md, third table) produced fixes that went into the code on 6 October 2026, and the chain was
rerun the same day on the PubTator3 bulk file that NCBI published on 6 October (section 7); the counts in this
document and in release v0.4 are from that build unless a paragraph names the 3 October build. The rule-based appraisal and the survey
(experiments/survey_literature.py, docs/literature_survey.md) join the two relation tables with the document
metadata into literature_reports.parquet, which enters the assembled table and the reliability fit through
assemble_evidence_table.py --extra-reports as grade E soft priors. Design sections 3.1, 4.2 (E3 paragraph), 4.3
and 6.1; report schema in docs/evidence_reports_spec.md section 1.

Literature reports are soft priors only (grades C to E). They never become evaluation positives, a pair that no
paper mentions is unlabelled rather than negative, and every row keeps the exact extraction model behind it
(PubTator3 relation type, CTD direct-evidence value, descriptor level, chemical match method) so that each report can be
weighted by the model that produced it and by its limitations, not by its count. Chemicals are matched to SIDER by
name (one SIDER drug per MeSH chemical, section 4), genes by NCBI Gene id; the study design is appraised from
PubMed publication types and PubTator3 species annotations by rule (rule_based_evidence_appraisal.py, section 6a),
not from the text of the papers; the document-grounded language-model appraisal (llm_evidence_appraisal.py)
remains unconfigured.

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
| experiments/fetch_pubtator_documents.py | document metadata per PMID: species and year from the PubTator3 export, publication types from esummary; capped, directional PMIDs first, resumable (section 6a) |
| mechanistic_pathway_learning/evidence/rule_based_evidence_appraisal.py | study design, species flags, rubric features and rubric weight from the document metadata (section 6a) |
| experiments/survey_literature.py | joins relations, documents and appraisal into literature_reports.parquet, writes literature_summary.json and docs/literature_survey.md with the overlap against the assembled table |
| tests/test_literature_relations.py | fixtures for relation-type mapping, Disease-Gene signs, mixed-polarity demotion, descriptor-level flags, gene id precedence, one drug per chemical, CTD evidence mapping, name matching, source names, schema |
| tests/test_rule_based_evidence_appraisal.py | design precedence (review, then species, then publication types), weights, the unfetched-metadata design |

The stub load_literature_predications.py stays in place; the two loaders above supersede it. SemMedDB is still
gated on the UMLS license.

Outputs (gitignored under data/processed/literature; released under data/releases/<version>/literature without
article text): pubtator_relations_filtered.parquet, ctd_relations_filtered.parquet, gene_identifier_map.parquet,
sider_chemical_mesh_map.parquet, pubmed_publication_dates.json, pubtator_relations_summary.json,
ctd_relations_summary.json, documents.parquet (titles and abstracts stay inside data/processed), documents_summary.json,
literature_reports.parquet and literature_summary.json; raw downloads and response caches under data/raw/pubtator3,
data/raw/ctd and data/raw/hgnc with a pin file next to each download.

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
| elevated_mood_or_mania | D001714 | mixed_polarity_diagnosis | Bipolar Disorder | F03.600.150.500 | diagnosis | yes | @DISEASE_Bipolar_Disorder |

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
- Bipolar Disorder D001714 carries the level mixed_polarity_diagnosis (added after the review): the diagnosis has
  the opposite mood pole as a cardinal feature, so a treat or cause relation on it says nothing about the direction
  of mania (lamotrigine, approved for bipolar depression and maintenance, had 775 "relieves mania" reports through
  this descriptor against 11 through Mania itself). The loaders keep the rows but demote induces and relieves on
  this descriptor to associated_with, with a count and a limitation clause; the survey table shows the
  symptom-level, diagnosis-level and demoted shares per symptom. Depressed mood rests on diagnosis-level
  descriptors only (Depression D003863 is outside the MEDIC vocabulary), which the same table shows.

## 3. Gene identifier map

Human-GEM model/genes.tsv supplies the NCBI id (geneEntrezID) and Ensembl id for the metabolic genes; the HGNC
complete set (pinned in section 7) supplies the rest by current symbol, then by a previous symbol carried by
exactly one gene, then by an alias carried by exactly one gene. Human-GEM comes first because the metabolic layer
was built from it; an ambiguous previous or alias symbol is left unmapped rather than guessed.

Build of 3 October 2026: 12,537 of 12,627 gene nodes carry an NCBI Gene id and 90 do not.

| Mapping source | Gene nodes |
|---|---|
| hgnc_alias_symbol | 1 |
| hgnc_previous_symbol | 20 |
| hgnc_symbol | 9,681 |
| human_gem_ensembl | 2,835 |
| unmapped | 90 |

Two NCBI ids are shared by two graph nodes each, a gene present under its current and its previous symbol: GBA and GBA1, SLC22A18 and SLC67A1. After the review a relation on such an id yields one report per node; in the 3 October build it reached the first node only.

## 4. SIDER drug to MeSH chemical matching (first pass, by name)

A SIDER drug matches a MeSH chemical when the MeSH descriptor name, or the entry term the PubTator3 autocomplete
reports as the matched synonym, equals the SIDER name case-insensitively (methods autocomplete_name and
autocomplete_synonym, one cached query per drug name); chemical ids on kept rows that no SIDER name matched are then
looked up by preferred label at id.nlm.nih.gov and matched again by name (method mesh_label_lookup). CTD rows match
by ChemicalName first (method ctd_chemical_name) and otherwise by a MeSH id the PubTator3 pass already matched. The
method is a column of every drug report. Salts, combination products and spelling differences ("gamma-aminobutyric"
against "gamma-Aminobutyric Acid") do not match in this pass and are counted, not guessed; a structure-based mapping
(PubChem CID to MeSH through UniChem or MeSH registry numbers) is the second pass.

One SIDER drug per MeSH chemical. SIDER lists brand names, salts and code names as separate drugs (olanzapine and
Zyprexa, rizatriptan and MK-462), so one MeSH id matched several STITCH ids and every paper on it was counted under
each (22.5 percent of the PubTator3 rows and 1,884 pairs in the 3 October build). After matching,
sider_chemical_matching.collapse_matches_to_one_drug_per_mesh_id keeps one STITCH id per MeSH id, in this order of
preference: the match method (exact descriptor name, then matched synonym, then label lookup, then CTD name), then a
STITCH id with a ChEMBL mapping in data/raw/chembl/pubchem_to_chembl_with_parents.json (the drugs the assembler can
join to targets), then the lowest PubChem CID. The dropped ids are listed in pubtator_relations_summary.json
(collapsed_sider_drugs_by_mesh_id); filter_ctd_rows applies the same rule, preferring the id the PubTator3 map kept.

Build of 3 October 2026, before the one-drug-per-chemical rule: of 1,430 SIDER drugs, 1,073 matched by descriptor name, 137 by matched synonym and 220 not at all. Of 9,880 chemical ids on rows with a target descriptor, 8,897 matched no SIDER name after the autocomplete pass; all 8,897 were looked up by label and none added a match, so 137,110 rows were dropped as chemicals without a SIDER drug.

## 5. PubTator3 bulk relations

The bulk file relation2pubtator3.gz was 297,478,945 bytes (under the 400 MB download cap), so it was downloaded
to data/raw/pubtator3 and streamed once. Columns: PMID, relation type, entity 1, entity 2, entities written
"Type|Identifier" (Disease|MESH:D001007, Gene|8813, Chemical|MESH:D005473, Variant and Species forms); no date.
Rows are kept when one entity is a crosswalk descriptor and the other is a graph gene (NCBI id in the map of
section 3) or a matched SIDER drug. Relation types map as follows; nothing else is kept and every drop is counted.

| PubTator3 relation type | Repository relation | report_value | Note |
|---|---|---|---|
| cause | induces | 1 | on target rows the partner is a chemical, an SNP or a mutation; gene partners carry no cause rows |
| treat, prevent | relieves | 1 | |
| stimulate, inhibit on a Disease-Gene row | induces | 1 | the bulk file renders a gene's positive or negative correlation with a disease as these types; the perturbation sign is +1.0 for stimulate (gene activity rises with the disease) and -1.0 for inhibit (gene activity falls with the disease, the loss-of-function direction of the monogenic class); the assembler weights them as ASSOCIATED_WITH, a correlation, not a cause |
| stimulate, inhibit on any other row | dropped | | chemical-gene relations; counted (none reach a target descriptor with a chemical partner in the pinned file) |
| positive_correlate, negative_correlate | induces, relieves | 1 | mapped for completeness; no row of the pinned file carries them against a Disease entity |
| associate | associated_with | 1 | kept as its own undirected relation, not merged into induces |
| cotreat, compare, interact, drug_interact | dropped | | not perturbation-symptom relations; counted |

A directional relation on the mixed-polarity descriptor D001714 is recorded as associated_with (section 2).

Report columns follow docs/evidence_reports_spec.md section 1 (report_id
"PubTator3-{relation type}|PMID:{pmid}|MESH:{descriptor}|{symptom}|{ordinal}", source "PubTator3-{relation type}",
the form assemble_evidence_table.predication_type_for_literature_sources parses (cause gives CAUSES 0.10; treat and
prevent AFFECTS 0.08; associate, stimulate and inhibit ASSOCIATED_WITH 0.05) and which makes each extraction type its
own sensor in the reliability fit; evidence_class "literature"; evidence_code "pubtator3_{relation type}", kept
through the appraisal; references "PMID:{pmid}"; rubric one-hots 0; perturbation_nodes the gene node with the sign
of the relation type for stimulate and inhibit rows and the loss-of-function sign as a documented convention
otherwise, or an empty list for drugs, which the assembler joins to ChEMBL targets by STITCH flat id as it does for
SIDER rows; the drug model_description carries no "ChEMBL targets" marker, so the assembler's target parser returns
nothing for literature rows) plus the literature columns pmid, relation_type, entity_role (disease_first or
disease_second), perturbation_ncbi_gene_id, chemical_mesh_id, chemical_match_method, mesh_descriptor,
mesh_descriptor_level, publication_year and direct_evidence (CTD only). A gene whose NCBI id resolves to two graph
nodes (GBA and GBA1, SLC22A18 and SLC67A1) yields one report per node. Publication dates come from esummary (pubdate
parsed to an ISO date, first of the month or year when the finer part is missing).

Build of 6 October 2026 from the bulk file pinned in section 7: 40,411,999 rows, 390,996 with a target descriptor on either side (3 October build on the 17 August file: 40,094,383 and 387,756).

| Relation type, whole file | Rows |
|---|---|
| associate | 18,816,263 |
| treat | 7,372,138 |
| negative_correlate | 4,106,519 |
| cause | 3,784,942 |
| positive_correlate | 3,487,922 |
| stimulate | 794,630 |
| inhibit | 627,652 |
| cotreat | 569,147 |
| compare | 542,110 |
| interact | 268,431 |
| prevent | 34,414 |
| drug_interact | 7,831 |

| Filter outcome on target-descriptor rows | Rows |
|---|---|
| rows kept | 209,184 |
| dropped chemical not matched to sider | 138,435 |
| dropped gene not in graph | 32,867 |
| dropped partner type:SNP | 5,550 |
| dropped partner type:ProteinMutation | 2,990 |
| dropped partner type:DNAMutation | 1,967 |
| dropped partner type:Mutation | 3 |

| Kept rows needing a rule of section 2 or 4 | Rows |
|---|---|
| Disease-Gene stimulate rows, kept as signed induces (before the gene-in-graph check) | 8,715 |
| Disease-Gene inhibit rows, kept as signed induces (before the gene-in-graph check) | 5,441 |
| directional rows on D001714 demoted to associated_with: treat | 13,946 |
| directional rows on D001714 demoted to associated_with: cause | 1,014 |
| directional rows on D001714 demoted to associated_with: stimulate | 461 |
| directional rows on D001714 demoted to associated_with: inhibit | 386 |

Kept: 209,184 reports over 133,193 papers, 5,011 genes and 985 drug ids; 209,068 dated through esummary (85 PMIDs not returned). By relation: associated_with 85,895, induces 47,107, relieves 76,182. Gene rows: associate 54,151, stimulate 6,277, inhibit 3,492; after the demotion on D001714, 8,922 gene reports are signed induces and 54,998 associated_with. SIDER drugs matched by MeSH chemical: 1,210, of which 82 are collapsed into another drug on the same MeSH id (64 MeSH ids had several), leaving 1,128.

The 3 October build (17 August file, rules before the review) kept 226,347 reports over 127,866 papers, 4,741 genes and 1,062 drug ids, with the 8,629 stimulate and 5,400 inhibit rows dropped as chemical-gene relations. Drug reports fell from 172,743 to 145,264, mostly through the one-drug-per-chemical rule (the earlier build counted a paper once per matched SIDER id), and gene reports rose from 53,604 to 63,920 with the restored stimulate and inhibit rows; the replacement bulk file adds 0.8 percent rows on its own.

## 6. CTD chemical-disease statements

CTD_chemicals_diseases.tsv.gz rows are kept when DiseaseID is a crosswalk descriptor and the chemical matches a SIDER
drug (section 4). DirectEvidence marker/mechanism maps to induces and therapeutic to relieves; a row carrying both
yields two reports; rows with empty DirectEvidence are inferred through a gene and are dropped with a count. Each
PubMed id of a row is its own report (report_id "CTD-curated|MESH:{chemical}|MESH:{descriptor}|{symptom}|{ordinal}",
source "CTD-curated", which the assembler reads as CAUSES for induces and AFFECTS for relieves, evidence_code
"ctd_direct_evidence_{marker_mechanism|therapeutic}"); a direct-evidence row that cites no paper is one report with
no reference. Marker/mechanism covers correlation as well as causation, which the limitations text says. A row that
carries both values is counted once under rows_with_both_direct_evidence_values. Statements on the mixed-polarity
descriptor D001714 are demoted to associated_with as in section 5, and one SIDER drug is kept per chemical (section 4).

Build of 6 October 2026 (same CTD file as 3 October): 9,903,450 CTD rows, 75,093 on a target descriptor, 14,739 of those with a SIDER chemical; 13,003 of these were inferred through a gene and dropped; 1,680 rows had a chemical matching several SIDER ids and were collapsed to one. Kept: 5,191 reports from 3,356 papers on 457 drug ids, marker/mechanism 2,851, therapeutic 2,340; by relation induces 2,590, relieves 1,827 and associated_with 774 (166 CTD rows on D001714, 98 marker/mechanism and 68 therapeutic, demoted); 5,190 dated; by chemical match method autocomplete_synonym 239, ctd_chemical_name 4,952. The 3 October build kept 5,607 reports on 490 drug ids, counting a paper once per matched SIDER id.

## 6a. Documents and the rule-based appraisal

experiments/fetch_pubtator_documents.py fetches, for the distinct PMIDs of both relation tables, the title, abstract,
year, journal and species annotations from the PubTator3 biocjson export (100 PMIDs per request) and the publication
types from E-utilities esummary (200 per request), cached per batch and resumable. The fetch is capped (--max-pmids,
30,000): PMIDs that carry a directional relation on any of their rows in any table (cause, treat, prevent, stimulate,
inhibit, the correlation types, CTD marker/mechanism and therapeutic) come first in first-seen order, associate-only
PMIDs fill the rest, documents already in documents.parquet are kept without a request, and the number of PMIDs
left without metadata is written to documents_summary.json. Titles and abstracts never leave data/processed.

rule_based_evidence_appraisal.py assigns each report a study design from the metadata, in this order: a review,
systematic review or meta-analysis publication type gives review_or_secondary; a paper whose only annotated species
is non-human gives animal_pharmacological (drug rows) or animal_genetic_perturbation (gene rows) whatever its other
publication types say; then Randomized Controlled Trial gives human_randomized_trial, the clinical trial, cohort,
case-control and observational types human_cohort_or_case_control (Comparative Study and Multicenter Study describe
scope, not design, and are not in that set), Genome-Wide Association Study or Genetic Association Studies
human_genetic_association and Case Reports human_case_report; anything else is not_reported. A report whose PMID is
beyond the fetch cap gets metadata_not_fetched, distinct from not_reported, so a paper nobody looked at is not
weighted like one that was read. Weights: randomised trial 1.0, cohort or case-control 0.9, genetic association 0.8,
case report 0.6, animal 0.5, review or secondary 0.4, not_reported 0.35, metadata_not_fetched 0.3, times 0.7 on a
diagnosis-level or mixed-polarity descriptor, floor 0.05. The rubric features (human species, animal only, the design
one-hots, symptom-level descriptor, publication year known, metadata available) and rubric_weight travel with the
report into the assembled table and the reliability fit; study_design is its own column and is prefixed to
model_description, while evidence_code keeps the extraction type.

Build of 6 October 2026: 134,747 papers across the two relation tables, 94,410 of them with a directional relation. The
cap selected 30,000 papers, all directional; 27,353 of them were already in documents.parquet from the 3 October
build and 2,647 were fetched new. The 2,647 documents of the 3 October selection that fall outside this selection stay
in documents.parquet and are used, so metadata exists for 32,647 papers (26,473 with species annotations, 32,638 with
publication types); a build from an empty cache would have 30,000. 102,100 papers have no metadata, among them 61,763
with a directional relation. 214,375 reports appraised; median rubric weight 0.245 (quantiles 0.1 to 0.9: 0.21, 0.21,
0.245, 0.30, 0.42).

| Study design | Reports, 6 October build | Reports, 3 October build and rules |
|---|---|---|
| metadata_not_fetched | 159,342 | (in not_reported) |
| not_reported | 19,078 | 190,724 |
| review_or_secondary | 10,097 | 10,752 |
| human_randomized_trial | 8,602 | 9,311 |
| human_cohort_or_case_control | 4,384 | 8,208 |
| human_case_report | 6,128 | 6,807 |
| animal_pharmacological | 6,686 | 6,142 |
| animal_genetic_perturbation | 58 | 10 |

Under the 3 October rules not_reported merged papers never fetched with papers read and uninformative, and Comparative
Study counted as a human cohort design; the review's rules split the first and put species before publication types,
which is most of the fall in human_cohort_or_case_control and the rise in the animal designs. Three quarters of the
reports (159,342) rest on papers beyond the cap, so the appraisal grades a minority of the literature; raising
--max-pmids is the lever, at about one biocjson and one esummary request per 100 and 200 papers.

## 7. Pins

| Download | Pinned file | Notes |
|---|---|---|
| PubTator3 relation2pubtator3.gz | 299,925,917 bytes; Last-Modified Tue, 06 Oct 2026 15:44:13 GMT; 40,411,999 rows; sha256 ffde6b53c1798363547c2fdf351c3385c21a249fcd3fe4be8dfe8dad3203cd76 | downloaded 6 October 2026 16:59 UTC; the file behind the 6 October build and release v0.4 |
| PubTator3 relation2pubtator3.gz, superseded | 297,478,945 bytes; Last-Modified Mon, 17 Aug 2026 13:41:35 GMT; 40,094,383 rows; sha256 6fca7a4c6a7b6fb9b727d564b1fc1aff1b04c0dfcadbc53aba188fb844e3ac42 | the file behind the 3 October build; NCBI replaced it on 6 October 2026 and no longer serves it, and the local copy was overwritten, so the 3 October build cannot be regenerated |
| CTD_chemicals_diseases.tsv.gz | 162,343,480 bytes; Last-Modified Tue, 29 Sep 2026 17:21:06 GMT; report created Tue Sep 29 13:17:30 EDT 2026; 9,903,450 rows; sha256 e11e6dc36a27d36ed88576e529a92211186e3e0b3cea4271cc78dcf3d1b821ee | downloaded 3 October 2026; unchanged at the 6 October check (same Last-Modified and sha256) |
| HGNC hgnc_complete_set.txt | 16,963,116 bytes; Last-Modified Fri, 02 Oct 2026 13:45:08 GMT; 45,187 rows; sha256 2b4224ea847df2fc6982f5b2a52804c5d92fbb8810b134afb636f6029452dc03 | downloaded 3 October 2026 |

## 8. Status of the join and what is still open

- Joined: literature_reports.parquet is appended to the assembled report table and the reliability fit by
  assemble_evidence_table.py --extra-reports. A literature-only pair is a grade E record whose weight is the
  predication weight of its extraction type (section 5); a pair that a label source also covers keeps that source's
  grade, the literature rows counting only in the reliability fit. The experiment loader keeps grades A and B as
  labels and refuses D and E (experiment_data.SOFT_PRIOR_ONLY_GRADES), so no literature row can become a label or an
  evaluation positive. associated_with stays a third relation that the induces experiments never read.
- Not yet built: the training-side use of the soft priors (when it is, literature rows must be excluded for every
  perturbation of the evaluation fold, since a report on a held-out gene's symptom would leak its label), the
  structure-based chemical mapping, the species check on gene rows through the Species annotations (an animal
  paper can mention a human NCBI Gene id; the appraisal now records animal-only papers, the loader does not drop them),
  the document-grounded language-model appraisal and the SemMedDB predications gated on the UMLS license.
- Open question 13 of the design concerns how gene association reports should count; after this build the gene
  rows also carry induces reports with a sign from the stimulate and inhibit types, so the question narrows to the
  associate rows.
