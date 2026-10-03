# Mechanistic pathway learning: experiment design, version 0.4

Status: 2 October 2026. Version 0.1 was written from the project notes before the repository was public; version 0.2 incorporated the decisions taken since (no hard symptom-level constraints, with strong evidence expressed as a 0.99 target; the earlier in-house literature knowledge graph is not used; compute is a runtime-limited GPU cluster). Version 0.3 replaced the Phase 1 estimates with measured counts from the HPO 2026-09-01 release, Human-GEM 2.0.1, SIDER 4.1 and ChEMBL 37 (docs/phase1_counts.md) and verified the HPO identifiers in the symptom crosswalk. Version 0.4 (3 October 2026) audits what the HPO term expansion pulls into each symptom and excludes diagnosis-level and off-construct terms (section 3.1, assumption A7), replaces the grade A proxy with annotation provenance plus frequency weights (section 4.3), adds a disease-cluster leakage group and a reconstruction-subsystem pathway split after measuring that the gene-wise split leaks and the curated-module split is too thin (section 6.1), changes how a module reads the encoder (the perturbation difference field with sum pooling, section 5.3), repairs a pooled-scoring artifact (section 6.2) and records the first proposed-model results on the monogenic slice (section 12). Numbered citations are in order of first appearance; PMIDs and DOIs are given in the reference list for import.

Implemented and tested: the noisy-OR pathway module head with hard-concrete supports and sum or mean pooling, the soft constraint losses, a relational message passing encoder with absolute and difference field readings, the sigmoid-head baseline B3, grouped (gene or disease cluster), curated-module, subsystem and time splits, the HPO loader with frequency qualifiers and provenance, the evidence grading policies, the learned reliability model, the document-grounded appraisal rubric, currency tagging, the equifinality indices, the metrics, the label-permutation and rewiring controls, the graph build (Human-GEM metabolic layer with subsystems; OmniPath, CollecTRI and small-molecule layers when their files are present), the evidence assembly, the popularity and random-walk baselines, a resumable training harness with validation early stopping, a batch driver, an aggregator and a module-analysis script. Still stubs: OnSIDES and literature loaders, the diagnosis proxy audit on the pharmacological side, the flux-sampling job (written, not run), the embedding and flux-feature and language-model baselines and the MAGMA wrapper.

## 1. Problem statement

### 1.1 Clinical problem

Psychiatric symptoms cross diagnostic boundaries and arise from heterogeneous biology, but mechanistic research is organized by diagnosis. A diagnosis-level mechanism averages over patients whose shared symptom may come from different causes, and a symptom-level mechanism (what produces psychosis, whichever disorder it is filed under) has no systematic map. The clearest evidence that single biochemical lesions produce psychiatric symptoms comes from inherited metabolic disorders: acute intermittent porphyria, Wilson disease, urea cycle defects, homocysteine remethylation disorders, Niemann-Pick type C and others present with psychosis, depression, anxiety or catatonia, sometimes for years before neurological or somatic signs appear [1]-[4]. One review describes hereditary metabolic disorders as a "rare but important cause of psychiatric disorders in adolescents and adults" [2]. These disorders are rare, but each is a natural perturbation experiment with a known enzyme, a known reaction and a documented symptom profile. Pharmacology supplies a second perturbation class: drugs with characterized targets that induce or relieve specific symptoms. Symptom-level genetics supplies a third: item-level GWAS of the nine PHQ-9 items found nine loci with "no overlap in loci across symptoms" [5], consistent with distinct biology per symptom.

### 1.2 Data science problem

Given a typed, compartmentalized physiology graph G (metabolites, reactions, enzymes and transporters, receptors, transcription factors) and a sparse, noisy, evidence-graded set of observations of the form "perturbation p of G induces (or relieves) symptom s", learn a function from a localized perturbation of G to a distribution over symptoms, together with an explanation in the form of one or more subgraphs (pathways) through which the perturbation acts. The function must (a) respect hard structural constraints from biochemistry and pharmacology, (b) accept soft constraints from literature evidence of graded quality, (c) use no diagnosis labels as inputs or targets and (d) allow several independent pathways to be individually sufficient for the same symptom without requiring one another.

### 1.3 Innovation

1. Supervision from perturbations with known molecular entry points (monogenic enzyme defects, single-target drugs) rather than from literature co-mention. This turns the task from knowledge-graph completion into mechanism inference.
2. A heterogeneous graph that joins a genome-scale metabolic reconstruction (with organelle compartments) to receptor signaling and transcriptional regulation, with an optional mechanistic flux layer that encodes what moves downstream of an enzyme lesion.
3. A noisy-OR output over learned sparse pathway modules. This encodes equifinality [6] directly: each module is a candidate sufficient cause and modules are not required to co-occur [7].
4. Orthogonal validation on symptom-level GWAS and metabolomic associations that never enter training, plus a time split, so the test measures mechanism recovery rather than literature reproduction.
5. An evaluation protocol with perturbation-wise and pathway-wise splits, degree-stratified reporting and negative controls, addressing the known failure modes of biomedical graph learning [8], [9].

## 2. Assumptions under challenge

Each item names an assumption in the project notes that is wrong or underdetermined as stated, followed by the design response.

A1. The graph should be full physiology. A generic human reconstruction has organelle compartments but no cell types and no brain regions. Recon3D covers "13,543 metabolic reactions involving 4,140 unique metabolites" [10] across the whole body; most of those reactions are irrelevant to any psychiatric symptom and the ones that matter are region- and cell-type-specific (dopamine in striatum versus prefrontal cortex; glutamine-glutamate cycling between astrocyte and neuron). Cell-type-resolved brain metabolic models exist [11] but add a modelling layer that is not needed to answer the first question. Design response: version 1 uses a brain-relevant subnetwork (section 3.2) with the full reconstruction retained as the background graph for propagation; cell-type resolution is a version 2 extension.

A2. Receptors are metabolized, so they belong in the metabolic graph. Protein turnover does not make a receptor a metabolic node. What a receptor does is signal transduction, which is not mass-balanced and does not fit stoichiometric edges; representing receptor binding as a reaction would corrupt flux-based features. Design response: one typed graph with separate edge semantics: catalysis and transport edges from the reconstruction; binding, activation and inhibition edges from a signaling resource; transcriptional regulation edges from a regulon resource. The loop described in the notes (neurotransmitter, receptor, effector, transcription factor, enzyme expression) is represented as a path across three edge types rather than as a reaction.

A3. Strong mechanistic evidence in psychiatric papers justifies hard constraints. Much of what counts as strong evidence in psychiatry comes from animal models and has not transferred to humans; even the dopamine hypothesis has been revised repeatedly [12]. A hard constraint is irreversible under training: a wrong one cannot be unlearned. Design response (agreed): hard constraints only for facts that are biochemical or pharmacological (stoichiometry, compartment, gene-reaction assignment, measured drug-target affinity). Every symptom-level link, however strong, enters as a weighted training observation or a soft prior. The policy discussed for strong evidence, a 0.99 "no-change" weighting, is kept and given two precise meanings in section 5.4: at the outcome level, strong positives are trained toward a target of 0.99 rather than 1.0; at the parameter level, a designated link can be anchored to its prior with a penalty equal to 99 pseudo-observations (strength 0.99 maps to 0.99 / 0.01). Both remain finite, so enough contradicting evidence can still move a link, which is what distinguishes them from a hard constraint. The anchored version is off by default and on in one ablation arm so the two policies are compared on held-out data.

A4. Symptom-level supervision is available. The literature is indexed by diagnosis. Symptom-level supervision exists in three places: phenotype ontologies attached to monogenic disorders, drug-label adverse events (MedDRA preferred terms such as anhedonia, insomnia or psychomotor agitation) and item-level questionnaire GWAS. Each has a confound. Drug-label events are not necessarily brain-mediated (fatigue on a beta blocker is cardiac output), are confounded by indication (antipsychotic labels list psychiatric events because of the treated population) and follow reporting conventions. Phenotype ontologies have their own confound: a symptom term's is_a closure can contain diagnoses and other constructs (section 3.1). Design response: restrict pharmacological evidence to CNS-penetrant drugs with a dominant target, model induction and relief as signed relations, run a peripheral-event negative control (section 6.3) and audit every ontology expansion term by term before it is used (docs/hpo_term_audit.md).

A5. Literature-derived training and literature-derived testing. This measures literature reproduction. Design response: orthogonal validation (section 6.4) and a time split (section 6.1) are mandatory parts of the primary analysis.

A6. Multiple independent pathways lead to the same symptom. Equifinality is plausible, but the competing hypothesis is a final common pathway: diverse upstream lesions converging on one downstream node set, as proposed for psychosis and striatal dopamine [12]. The two hypotheses produce different graph structures (disjoint supports versus converging supports). Design response: the noisy-OR model permits either outcome; section 6.5 defines an independence index and a convergence index and pre-registers both outcomes as findings.

A7. Removing diagnosis labels makes the learning diagnosis-free. Diagnosis leaks through proxies: monogenic phenotype annotations are attached to disease entries, antipsychotic labels track their indication and symptom GWAS cohorts are population samples with diagnosed subgroups. Two further leaks were measured in version 0.4. First, the ontology: the is_a closure of Cognitive impairment (HP:0100543) contained Attention deficit hyperactivity disorder (HP:0007018), a diagnosis, which supplied 76 of the 262 Human-GEM genes for that target and 41 of them through no other term; Restlessness contained Restless legs and Hallucinations contained sleep-related hallucinations. Second, the annotation unit: HPO annotations are made per disease and inherited by every gene of that disease, so two genes of one disease carry identical symptom profiles and sit next to each other in the metabolic graph. With genes grouped by shared disease entry (340 clusters for 451 genes, largest 27 genes; the largest cluster supplies 5 to 17 percent of a symptom's positive genes and the three largest 11 to 27 percent, data/processed/evidence/evidence_summary.json) the random-walk baseline's macro AUPRC fell from 0.347 to 0.293 per fold, which is within noise of its own label-permutation control (0.284); under gene-wise grouping the same baseline looked informative (section 6.1). Design response: no disease nodes in the graph, symptom terms only as targets, an excluded_hpo_ids column in the crosswalk with an audit table that lists every expansion term and what it contributes, the disease-cluster grouping as the default grouped split and a proxy audit in Phase 1 that tests whether drug class or disease identity can be predicted from the evidence features.

A8. Experiments are tractable with existing data. Computation is not the constraint; curation is. Measured (docs/phase1_counts.md): 451 Human-GEM genes carry at least one HPO annotation in a target symptom after the version 0.4 exclusions (473 before; cognitive impairment falls from 262 to 202 genes), and 65 SIDER drugs pass the version 1 pharmacological bar (a single ChEMBL mechanism target and an ATC nervous-system code). Ten of twelve symptoms clear the go/no-go threshold on these two classes alone; anhedonia (0 genes, 0 drugs) and psychomotor retardation (0 genes, 4 drugs) do not. The 65-drug count is small because most centrally acting drugs have several mechanism targets; relaxing the rule to two targets (section 4.2) is the first lever if the pharmacological class needs to be larger, and the symptom set shrinks rather than the evidence bar.

A9. Hub nodes are a nuisance to be handled later. Currency metabolites (ATP, NADH, water, protons) and hub symptoms dominate any propagation- or embedding-based model, and biomedical knowledge-graph embeddings are known to produce "densely connected entities being highly ranked no matter the context" [8]. Design response: currency metabolites are removed from path-based features (kept for flux computations), negatives are degree-matched and all metrics are reported by degree bin.

## 3. Scope for version 1

### 3.1 Symptom set

Twelve symptoms chosen for availability across all three evidence sources and for coverage of RDoC domains. Each symptom is a target with a crosswalk to the vocabularies used by the evidence sources (docs/symptom_crosswalk.csv). HPO identifiers were verified against the 2026-09-01 release: Anhedonia is HP:0012154; Agitation is HP:0000713 (with Restlessness HP:0000711); Hypersomnia (HP:0100786) is obsolete and replaced by Excessive daytime somnolence (HP:0001262); Psychomotor retardation (HP:0025356) was retired into Global developmental delay (HP:0001263), which is a different construct, so that symptom has no monogenic channel. Measured evidence per symptom is in docs/phase1_counts.md.

| Symptom (target) | HPO term | MedDRA preferred terms (examples) | Questionnaire item | RDoC domain |
|---|---|---|---|---|
| Depressed mood | Depression (HP:0000716) | Depressed mood; Depression | PHQ-9 item 2 | Negative valence |
| Anhedonia | Anhedonia (HP:0012154); no Human-GEM gene annotated | Anhedonia | PHQ-9 item 1 | Positive valence |
| Anxiety | Anxiety (HP:0000739) | Anxiety; Nervousness | GAD-7 items | Negative valence |
| Irritability or aggression | Irritability (HP:0000737); Aggressive behavior (HP:0000718) | Irritability; Aggression; Agitation | none | Negative valence; arousal |
| Insomnia | Insomnia (HP:0100785) | Insomnia; Initial insomnia | PHQ-9 item 3 | Arousal and regulatory |
| Somnolence or hypersomnia | Excessive daytime somnolence (HP:0001262) | Somnolence; Hypersomnia | PHQ-9 item 3 | Arousal and regulatory |
| Fatigue | Fatigue (HP:0012378) | Fatigue; Asthenia | PHQ-9 item 4 | Arousal and regulatory |
| Psychomotor agitation | Agitation (HP:0000713); Restlessness (HP:0000711); excluding Restless legs (HP:0012452) | Psychomotor hyperactivity; Agitation; Restlessness | PHQ-9 item 8 | Arousal |
| Psychomotor retardation | none (term retired in HPO) | Psychomotor retardation; Bradyphrenia | PHQ-9 item 8 | Arousal |
| Psychosis | Psychosis (HP:0000709); Hallucinations (HP:0000738); Delusion (HP:0000746); excluding Sleep related hallucination (HP:4000063) | Psychotic disorder; Hallucination; Delusion | none | Cognitive systems; perception |
| Cognitive impairment | Cognitive impairment (HP:0100543); excluding ADHD (HP:0007018), Progressive neurologic deterioration (HP:0002344), Motor deterioration (HP:0002333), Psychomotor deterioration (HP:0002361) | Memory impairment; Disturbance in attention; Confusional state | PHQ-9 item 7 | Cognitive systems |
| Elevated mood or mania | Mania (HP:0100754) | Mania; Hypomania; Euphoric mood | none | Positive valence; arousal |

Term expansion. Each HPO term is expanded to its is_a descendants, and the expansion is audited term by term (experiments/audit_hpo_term_expansion.py, docs/hpo_term_audit.md): for every descendant, the Human-GEM genes it contributes, the genes reachable through no other term, the provenance of the annotation rows and their frequency qualifiers. Two kinds of descendant are excluded through the crosswalk column excluded_hpo_ids, with the subtree below them: diagnosis-level terms, because a diagnosis is what assumption A7 keeps out of the targets (Attention deficit hyperactivity disorder under Cognitive impairment), and terms that name a different construct from the target symptom (Restless legs, a sensorimotor sleep disorder, under Restlessness; hypnagogic and hypnopompic hallucinations, sleep-wake transition phenomena, under Hallucinations; the motor and neurologic deterioration terms that HPO files under Mental deterioration). Dementia, Mental deterioration, Memory impairment, Confusion and Delirium stay in cognitive impairment; Akathisia stays in psychomotor agitation. The exclusions are a curation decision recorded as data, not a model constraint, and the audit table makes every one of them reversible.

Excluded in version 1: suicidal ideation (label events are dominated by class-wide boxed warnings rather than mechanism, so pharmacological evidence is uninformative) and appetite change (peripheral mechanisms dominate). Both are candidates for version 2 with genetics-only supervision. After the Phase 1 counts, anhedonia and psychomotor retardation stay in the crosswalk but are flagged as genetics-only targets: they have no grade A or B supervision, so they cannot be primary endpoints and enter only the GWAS enrichment analysis (section 6.4).

### 3.2 Pathway modules with perturbation anchors

The background graph is the full reconstruction. The brain-relevant subnetwork below is where learned modules are expected and where curation effort is concentrated. Anchors are the monogenic disorders and drug classes that give each module its grade A or B evidence (section 4.3). Module boundaries are curation scaffolding, not model constraints; learned modules are free to cross them, and the overlap between learned and curated modules is an interpretability output (section 6.6).

| Module | Enzymes, transporters, receptors | Monogenic anchors (gene) | Pharmacological anchors |
|---|---|---|---|
| Catecholamine synthesis and degradation | TH, DDC, DBH, PNMT, MAOA, MAOB, COMT, SLC6A3, SLC6A2, SLC18A2 | DBH deficiency; MAOA (Brunner syndrome); TH deficiency; SLC6A3 (dopamine transporter deficiency) | Amphetamines; levodopa; MAO inhibitors; D2 antagonists; VMAT2 inhibitors |
| Serotonin and melatonin | TPH1, TPH2, DDC, MAOA, AANAT, ASMT, SLC6A4 | none established | SSRIs; triptans; melatonin receptor agonists |
| Tetrahydrobiopterin cofactor | GCH1, PTS, SPR, QDPR, PCBD1 | BH4 deficiencies (GCH1, PTS, SPR, QDPR) | Sapropterin (relief) |
| Phenylalanine and tyrosine | PAH, TAT, HPD, HGD | Phenylketonuria (PAH); tyrosinemias | none |
| Tryptophan-kynurenine | TDO2, IDO1, KMO, KYNU, HAAO, KYAT1, KYAT3 | none established | Interferon alfa (indirect; soft) |
| GABA and glutamate | GAD1, GAD2, ABAT, ALDH5A1, GLUL, GLS, SLC1A2, GABA-A subunits, NMDA receptor subunits | SSADH deficiency (ALDH5A1); GABA transaminase deficiency (ABAT); GLUL deficiency | Benzodiazepines; vigabatrin; ketamine; memantine; phencyclidine |
| One-carbon, folate and homocysteine | MTHFR, MTR, MTRR, CBS, CTH, MAT1A, AHCY, MTHFD1, SLC19A1, FOLR1 | Remethylation disorders (MTHFR, MTR, MTRR); CBS deficiency; cerebral folate deficiency (FOLR1) | Methotrexate (soft); folinic acid (relief) |
| Urea cycle and ammonia | CPS1, OTC, ASS1, ASL, ARG1, NAGS, SLC25A15, SLC25A13 | OTC deficiency and other urea cycle defects | Valproate (hyperammonemia; soft) |
| Heme and porphyrin | ALAS1, ALAD, HMBS, UROS, UROD, CPOX, PPOX, FECH | Acute intermittent porphyria (HMBS); variegate porphyria (PPOX); hereditary coproporphyria (CPOX) | Porphyrinogenic drugs such as barbiturates (soft) |
| Copper and metal handling | ATP7B, ATP7A, CP | Wilson disease (ATP7B) | Chelators (relief) |
| Sterol, bile acid and lysosomal lipid | CYP27A1, NPC1, NPC2, DHCR7, GBA1, ARSA, HEXA, HEXB, ABCD1 | Cerebrotendinous xanthomatosis; Niemann-Pick C; Smith-Lemli-Opitz; Gaucher; metachromatic leukodystrophy; Tay-Sachs; X-linked adrenoleukodystrophy | Chenodeoxycholic acid (relief); miglustat |
| Creatine and energy metabolism | GAMT, GATM, SLC6A8, SLC2A1, PDHA1, POLG, respiratory chain subunits | Creatine deficiency syndromes; GLUT1 deficiency; PDH deficiency; POLG disorders | none |
| Purine and pyrimidine | HPRT1, ADA, ADSL, DPYD | Lesch-Nyhan (HPRT1); ADSL deficiency | Allopurinol (soft) |
| Acetylcholine, histamine and adenosine | CHAT, ACHE, SLC5A7, HDC, HNMT, ADORA1, ADORA2A | none established | Anticholinergics; cholinesterase inhibitors; H1 antagonists; caffeine |
| Steroid and thyroid hormone | CYP21A2, CYP11B1, HSD11B1, DIO2, DIO3, SLC16A2 | Congenital adrenal hyperplasia; Allan-Herndon-Dudley (SLC16A2) | Corticosteroids; levothyroxine (relief) |
| Vitamin cofactor transport | SLC19A3, BTD, TCN2, SLC52A2 | Biotin-thiamine responsive basal ganglia disease; biotinidase deficiency | Vitamin repletion (relief) |

### 3.3 Signaling and transcriptional layer

Ligand-receptor, receptor-effector and kinase-substrate edges with signs from OmniPath [13]; transcription factor to target regulons with signs from CollecTRI [14], restricted to transcription factors with brain expression evidence. This layer supplies the only route by which a receptor-level drug perturbation reaches enzyme expression.

## 4. Data layer

### 4.1 Physiology graph

Base reconstruction: Human-GEM [15], which is versioned on GitHub and actively curated, with a Recon3D [10] identifier map retained for cross-reference. Node types: metabolite (compartment-specific), reaction, gene, protein (enzyme, transporter, receptor, transcription factor) and drug. Edge types: substrate-of, product-of, catalyzed-by (from gene-protein-reaction rules), transports-between (compartments), binds, activates, inhibits (signed), regulates-transcription-of (signed) and targets (drug to protein, with action type and affinity from ChEMBL [16]). Node attributes: compartment; brain expression weight from Human Protein Atlas or GTEx brain tissues (cell-type expression as a version 2 attribute); blood-brain barrier penetration flag for drugs. Currency metabolites are tagged and removed from path-based features but kept for flux computation.

Hard structural constraints encoded in the graph: (1) no edges other than those listed; (2) no symptom-to-symptom edges; (3) no disease nodes; (4) compartment changes only through transport reactions; (5) drug-target edges only with measured activity (pChEMBL of at least 6, or a curated mechanism annotation); (6) sign consistency for pharmacological edges (agonist versus antagonist).

### 4.2 Perturbation-symptom evidence

Evidence class E1, monogenic (grade A). HPO gene-to-phenotype annotations [17] restricted to genes present in the graph and to phenotypes in the symptom set or their descendants after the exclusions of section 3.1. Each gene is a loss-of-function perturbation; its label is the set of annotated symptoms. The frequency column of genes_to_phenotype.txt (an HPO frequency qualifier from Obligate to Very rare, a patient fraction such as 3/5, or a percentage) is parsed to a midpoint and kept per pair as the largest reported frequency; rows qualified as Excluded assert absence and are dropped. In the 2026-09-01 release 94 percent of the 861 pairs on Human-GEM genes carry a frequency. The disease identifiers behind each pair are kept with their provenance (OMIM or Orphanet entry counts), which the grading and the disease-cluster split both use.

Evidence class E2, pharmacological (grade B). Drug-adverse event pairs from SIDER [18], which holds "1430 drugs, 5880 ADRs and 140 064 drug-ADR pairs" from labels current to 2015, for the earlier time slice, and from OnSIDES [19], which holds "3.6 million drug-ADE pairs for 3,233 unique drug ingredient combinations" from current labels, for the later slice. Drug indications from the same sources supply the relief relation; in SIDER only indications found by its NLP_indication method are kept, since text_mention rows include conditions that are merely named in the label. Drugs are mapped to ChEMBL through UniChem (PubChem identifier, stereo forms as fallback, salt forms resolved to parent molecules) and to their mechanism targets from the ChEMBL mechanism table; targets carry gene symbols and an action type that sets the sign of the perturbation (experiments/fetch_chembl_drug_targets.py). Inclusion filters, version 1: a single ChEMBL mechanism target (a protein-family target such as the GABA-A receptor counts as one), an ATC nervous-system code as a crude proxy for central action until a measured blood-brain barrier flag replaces it, and an event term that maps to the symptom set. The affinity-margin definition of dominance is deferred until ChEMBL activities are joined. Measured with these filters: 1,168 of 1,430 SIDER drugs map to ChEMBL, 580 have a mechanism annotation and 65 pass the single-target, ATC N filter with at least one target-symptom event.

Evidence class E3, literature (grades C to E; soft prior only). The earlier in-house psychiatric literature knowledge graph is not used; its quality was judged insufficient. E3 is derived from public resources instead, in order of preference: PubTator3 entity annotations and extracted relations (associate, cause, treat, inhibit, stimulate and correlation types) queried by target symptom MeSH descriptor and by graph gene symbol, with PMID and year kept for the time split; CTD curated chemical-disease associations with a therapeutic-versus-marker direction, for chemicals with known targets and MeSH symptom descriptors; and SemMedDB predications [20] (CAUSES, AFFECTS, DISRUPTS, PREDISPOSES, ASSOCIATED_WITH) once a UMLS license is in place, which is free for academic use. Weights follow predication type, count and recency. A targeted large-language-model extraction over the symptom-by-module space is a version 2 option only if it is validated against a hand-labelled gold set first. E3 never enters any test set.

Evidence class E4, symptom-level genetics (validation only in version 1). Gene-level statistics from item-level GWAS of the PHQ-9 [5] and of anxiety and depressive symptom factors [21], computed with MAGMA [22].

Evidence class E5, metabolomics (validation only; optional). Published metabolite-symptom association tables. No new cohort access is required.

### 4.3 Evidence grading and learned weighting

The grade table below is the fallback weighting and one ablation arm. The default weighting is learned: every source is treated as a noisy sensor of the true perturbation-symptom link with its own sensitivity and specificity, estimated by expectation-maximization from the pattern of agreement across sources (a Dawid-Skene model [27] with Beta priors and a feature-dependent form after Raykar et al. [28]; evidence_reliability_model.py). The posterior probability that a link is real is the observation's soft label and loss weight, and a link anchor's pseudo-observation count is read off that posterior (0.99 maps to 99) rather than set by hand. One global scale on the weights is the only hyperparameter left to validation-set tuning. A two-class latent model needs at least three conditionally independent sources to be identified, which the monogenic, label and literature classes supply. Weights learned as bare multipliers on the training loss would collapse to zero; weights learned as noise parameters inside the likelihood of the reports do not, which is why the model is specified this way.

Strength of findings enters through features rather than counts. A document-grounded language-model appraisal (llm_evidence_appraisal.py) extracts a structured rubric from each supporting document: study design, species, sample size, measurement level (symptom item versus diagnosis), effect direction, cohort identity and mechanistic specificity. The rubric is never a judgment of whether the link is true. Its fields become features of the reliability model, so sensitivity and specificity vary by study design within a source (the feature-dependent form of the model), and independence is counted by cohort rather than by mention. Before any appraisal enters training it is scored against a gold set rated by two humans, field by field with Cohen's kappa; the model name, prompt version and temperature are pinned, outputs are cached by content hash, and the time-split leakage audit in section 6.4 checks whether appraisal features help post-cutoff predictions more than within-period ones.

| Grade | Definition | Role | Initial loss weight |
|---|---|---|---|
| A | Human loss of function in a graph gene with the symptom in a curated clinical synopsis (OMIM or Orphanet provenance of the HPO annotation) or in at least two independent case series | Training label | 1.0, scaled by the HPO frequency qualifier when present |
| B | CNS-penetrant drug with a dominant target and the event on the label | Training label | 0.6, scaled by label frequency when available |

Grade A policy. Version 0.3 approximated "two independent case series" by the number of distinct disease identifiers behind a pair. That proxy counted an OMIM entry and its Orphanet counterpart for the same disease as two series, and it gave weight zero to genes with a single entry, HMBS among them: under it 408 of the 861 monogenic pairs on the slice (47 percent) are grade C with weight zero, excluded from training while remaining in the test sets (the alternative evidence table is kept for the ablation, data/processed/evidence_two_entries). Version 0.4 grades by provenance: every HPO gene-phenotype annotation descends from a curated disease annotation (OMIM clinical synopsis or Orphanet clinical-sign table), so provenance is the measurable form of the "clinical synopsis" clause, and the frequency qualifier carries the strength: weight = 1.0 times clip(4 times frequency, 0.25, 1), so a feature present in a quarter or more of patients keeps the full weight, Occasional (midpoint 0.17) gives 0.68 and Very rare the floor; unknown frequency keeps the base weight. The old policy remains selectable (configs/evidence_assembly.yaml, grade_a_policy) for the ablation that compares the two. A caveat stands: frequency is also an estimate of P(symptom | lesion), which is the quantity the model predicts; using it as a weight rather than as a soft target is a choice recorded as open question 8.
| C | Human association: gene-level GWAS statistic or metabolomic association with the symptom | Validation (version 1); soft prior (version 2) | none in version 1 |
| D | Animal perturbation with a symptom analogue | Soft prior | 0.2 |
| E | Literature predication or co-mention without study-type evidence | Soft prior | 0.05 to 0.1 by predication type |

### 4.4 Access and licensing

Public without registration: Human-GEM, Recon3D (VMH and BiGG), HPO, OnSIDES, SIDER, ChEMBL, UniChem, OmniPath, CollecTRI, PubTator3, CTD and GWAS summary statistics (GWAS Catalog). Registration: OMIM, for clinical synopses. License required: UMLS (free for academic use), for SemMedDB and MedDRA mappings; MedDRA term redistribution is restricted, so the public repository commits mappings to target symptoms, not MedDRA tables. Not required in version 1: individual-level cohort data. Pinned so far (docs/data_sources.md): HPO release 2026-09-01, Human-GEM 2.0.1, SIDER 4.1 (2015-10-21), ChEMBL 37 (2026-05-01); the pins become part of the pre-registration.

## 5. Model

### 5.1 Formal setup

Let G = (V, E, tau) be the typed graph, S the symptom set and P the set of perturbations. A perturbation p is a set of (node, sign, magnitude) triples. An observation is (p, s, r, y, w) with relation r in {induces, relieves}, outcome y in {0, 1} and weight w from the evidence grade. The model defines P(y = 1 | p, s, r).

### 5.2 Perturbation encoder

Two routes, compared as an ablation.

Route 1, propagation. A relational graph neural network over G using hard edges only. The perturbation is written onto its nodes as input; L rounds of typed message passing produce node states. The perturbation representation h_p is the field of node states over V rather than a pooled vector, so that modules can read local regions of it. Two readings of the field exist. The absolute field is the learned base state of every node plus the propagated perturbation; a readout of it can memorize which nodes exist. The difference field is the absolute field minus the field with no perturbation written on it (one extra batch-of-one pass), so it carries only what the perturbation changed and is zero everywhere the perturbation did not reach. The difference field is the default from version 0.4; the absolute field is an ablation arm, because the comparison between them is the measure of how much of the task is node identity rather than mechanism.

Route 2, flux. For enzyme and transporter perturbations, the gene's reactions are constrained in the reconstruction and the shift in the feasible flux distribution is computed by flux sampling with COBRApy [23]. Sampling is used rather than flux balance analysis because the brain has no defensible single objective function. The output is a per-reaction shift statistic (difference in sampled flux medians and in flux variability ranges). Receptor-level perturbations have no flux route and reach the metabolic layer only through route 1 via signaling and transcription edges.

Route 3 (default): concatenation of routes 1 and 2 on reaction nodes.

### 5.3 Pathway modules and noisy-OR output

K modules, each with a learned sparse support m_k in [0, 1]^|V| (hard-concrete relaxation of an L0 penalty [26]). Module activation a_k(p) = sigmoid(w_k . pool(m_k * h_p) + b_k), where pool is the sum over nodes (default) or the support-mass-normalized mean (ablation). The sum is the right pooling for the difference field: a perturbation that reaches a few support nodes is read at full strength and a support that nothing reaches reads as zero. With the mean over an absolute field, which version 0.3 used, the base states of thousands of gated nodes swamp a change at a few, and the activation becomes a property of the support rather than of the perturbation. Output:

P(y = 1 | p, s, induces) = 1 - (1 - lambda_s) * prod_k (1 - pi_{k,s} * a_k(p))

with pi_{k,s} in [0, 1] the module-to-symptom link and lambda_s a leak term for unmodelled causes. The leak can be initialized at the training base rate of each symptom (--init-leak-from-base-rate), so the untrained model predicts the base rate and the early epochs do not penalize every positive against a probability of a few percent; whether that changes the result is one of the configured comparisons. The relieves relation uses a separate link matrix and a signed activation. Equifinality is a property of this functional form: any single module with pi_{k,s} * a_k(p) near 1 suffices and no cross-module interaction term is needed [7].

### 5.4 Constraints

Hard: the graph itself (section 4.1). Nothing else is hard in the default configuration.

Soft, version 1 implementation, outcome level: every observation enters a weighted binary cross-entropy with its grade weight (section 4.3); observed positives are smoothed to a target of 0.99, so the model is asked to keep strong links near certain rather than certain. Unobserved pairs are sampled as degree-matched negatives at reduced weight because they are unlabelled, not negative.

Soft, version 1 implementation, parameter level: a designated module-to-symptom link can be anchored to a prior probability with a Bernoulli KL penalty whose coefficient is a pseudo-observation count; a prior strength of 0.99 is 99 pseudo-observations. With learned weighting the count comes from the reliability posterior (section 4.3). Disabled by default.

Soft, version 2 implementation: for each literature triple (entity e, symptom s, weight w), a hinge penalty max(0, tau * w - A(e, s)) on the attribution A from node e to output s. This nudges the explanation locally without fixing a pathway.

Optional hard lower bound: for designated (module, symptom) pairs, pi_{k,s} >= pi_min. Off by default; on in one ablation arm to test the project's original hard-constraint policy against held-out data.

### 5.5 Parsimony and identifiability

K and the support sizes are selected by a minimum description length criterion: total code length equals model cost (number of active support nodes and nonzero links, each at a fixed cost) plus the weighted negative log-likelihood of the observations. This is the parsimony principle used elsewhere in the research program and it replaces an arbitrary choice of K. The decomposition into modules is not identifiable (permutation, split and merge); stability selection across seeds and folds reports consensus supports with their selection frequency, and only consensus modules are interpreted.

### 5.6 Baselines

B0, popularity: predict each symptom at its base rate. Required because of hub bias.
B1, network proximity: random walk with restart from the perturbed nodes to symptom-anchored nodes [24].
B2, knowledge-graph embedding: TransE, ComplEx and RotatE on G with symptoms as nodes and evidence as edges, with hyperparameters chosen as recommended for drug-discovery graphs [25].
B3, relational GNN link prediction with a sigmoid head and no modules: the same encoder, a sum pool over all nodes, one hidden layer and an independent sigmoid per symptom (implemented; it shares the training loop with B6).
B4, flux-feature logistic regression (route 2 only).
B5, language-model zero-shot: symptom prediction from the gene or drug name alone, to measure how much of the task is text recall.
B6, the proposed model.

### 5.7 Ablations

No soft constraints; hard lower bounds on; compartments collapsed; signaling and transcription layers removed; flux route removed; propagation route removed; sigmoid head instead of noisy-OR; K = 1 and K = 16; absolute field with mean pooling instead of the difference field with sum pooling; description-length penalty off; currency metabolites retained in path features; fixed grade weights instead of learned reliability; two-distinct-entries grade A policy instead of provenance; learned reliability without rubric features.

## 6. Evaluation

### 6.1 Splits

Grouped perturbation-wise cross-validation, following the disjoint-group protocol shown to remove optimistic bias in drug knowledge graphs [9]. Drugs sharing a dominant target are assigned to the same fold. For genes two groupings are implemented: gene-wise (each gene its own group) and disease-cluster (connected components of the gene-disease graph over the target-symptom annotations, so genes annotated to one disease leave the training data together). The disease-cluster grouping is the default from version 0.4: HPO annotations are made per disease, genes of one disease inherit identical symptom profiles and often share reactions, and under gene-wise grouping the random-walk baseline reached a per-fold macro AUPRC of 0.347 against 0.234 for popularity, while under disease-cluster grouping it reached 0.293 against its own label-permutation control of 0.284 (docs/phase2_baselines_disease_cluster.md). The gene-wise number was mostly leakage.

Pathway-wise split: all perturbations writing onto a pathway are held out together, which tests whether the model generalizes to an unseen mechanism rather than interpolating within a known one. Two pathway definitions are implemented. The curated modules of section 3.2 are the pre-registered definition, but on the monogenic slice only six of sixteen reach ten positive pairs and together they hold 37 perturbations, so per-symptom metrics inside one module are not interpretable (held-out sets of 4 to 8 genes with homogeneous profiles). The Human-GEM subsystems (146 in release 2.0.1; every gene whose reactions mostly belong to subsystem S is held out with S) are the reconstruction's own pathway partition and cover many more annotated genes; they are proposed as the pathway-wise split for the primary endpoint, with the curated modules kept for the interpretability overlap (section 6.6). Held-out sets are scored pooled, per-set ranking metrics are always reported, per-set macro metrics only for sets of at least 20 perturbations, and every pathway-wise result is paired with the same split on permuted labels.

Time split: training on SIDER (labels to 2015) and HPO annotations released before a cutoff; testing on OnSIDES pairs for drugs first approved after 2016 and on HPO annotations added after the cutoff. Implemented for the monogenic class in version 0.4 from phenotype.hpoa: each (gene, symptom) pair is dated by the earliest availability date of the OMIM disease annotations behind it, where an annotation's date is the publication date of the PubMed reference it cites (all 204 references behind the target pairs were dated from PubMed on 3 October 2026, docs/hpo_reference_publication_dates.json, DOIs kept) and, for annotations that cite only the OMIM entry, the HPO biocuration date. Orphanet annotations all carry the release import date and leave their pairs undated. On the monogenic slice 408 of 861 pairs are dated, 144 of them after 31 December 2015; 97 pairs carry the 2009 bulk-curation date of the OMIM text import because they cite no publication. Training positives are the pairs dated on or before the cutoff; the scored pairs are every pair that is neither a training positive nor an undated positive, labelled by whether it was dated afterwards, with the same evaluation on permuted new-positive labels. The split is mostly a new-gene split: most genes with a post-cutoff pair had no dated pair before it.

### 6.2 Metrics

Per-symptom AUROC and AUPRC with bootstrap confidence intervals over perturbations; mean reciprocal rank and hits-at-3 for symptom ranking per perturbation; expected calibration error. All metrics stratified by evidence grade and by node degree bin. Out-of-fold predictions are pooled for the bootstrap and the ranking metrics, and the macro metrics are also reported as mean and standard deviation over folds: pooling has an artifact for any predictor that is constant within a fold, because a fold with a low training base rate has a high test base rate, which pulled the pooled AUROC of the popularity baseline to 0.44 in version 0.3 where 0.50 is the truth. Ranking metrics compare scores across symptoms within one perturbation, so they reward calibrated cross-symptom scores and favor the popularity baseline (MRR 0.65, hits-at-3 0.81 on the monogenic slice) over the random walk (0.48, 0.58), whose anchor-normalized scores are not comparable across symptoms; B0 is their reference, not their floor.

### 6.3 Negative controls

1. Symptom labels permuted within degree strata (run for every model and split; the random walk kept a per-fold macro AUPRC of 0.284 against 0.235 for popularity on permuted labels, which is the degree signal the stratified permutation preserves by design).
2. Degree-preserving rewiring of G (double-edge swaps within each relation type; run for the random walk under the grouped split). Performance should collapse toward B0 if the biology carries the signal.
3. Evidence grades shuffled.
4. Peripheral-event control: adverse events with no plausible central mechanism (rash, injection-site reaction). Predictions from brain modules should be at chance.

### 6.4 Orthogonal validation

Genetics. For each symptom, the top-N genes by attribution, excluding any gene in that symptom's training evidence, form a gene set. MAGMA competitive gene-set tests on the symptom-level GWAS summary statistics [5], [21], compared with size-matched random sets and with gene sets derived from B2. Power is limited (nine loci across all PHQ-9 items), so this analysis is pre-registered as secondary with a stated minimum detectable enrichment.

Metabolomics. Predicted symptom-linked metabolites compared with published metabolite-symptom associations by rank-based enrichment.

Prospective literature. Perturbation-symptom pairs predicted at high probability with no literature support at the training cutoff, checked against subsequent literature at analysis time.

### 6.5 Equifinality test (primary hypothesis)

For symptom s let M_s = {k : pi_{k,s} > theta} be the active modules.

Independence index I_s = 1 - max over module pairs in M_s of the Jaccard overlap of their supports.

Convergence index C_s = Jaccard overlap of the downstream reachable node sets of the modules in M_s within d hops.

Sufficiency test: let P_k be the held-out perturbations whose attribution mass lies in module k. Ablating a module k' different from k should leave AUPRC on P_k unchanged within the bootstrap interval; ablating k should reduce it. Paired bootstrap over perturbations. Implemented post hoc in experiments/analyze_pathway_modules.py from the saved test module activations, links and leaks: the symptom probability is recomputed with one module removed from the noisy-OR, P_k is the set of held-out perturbations whose largest contribution link_{k,s} a_k comes from k, and groups need at least eight perturbations and three positives to be scored.

Pre-registered expectations: psychosis and cognitive impairment will show at least two modules with I_s above 0.9 (porphyrin, sterol or lysosomal, urea cycle, catecholamine and glutamate anchors); depressed mood and anxiety will show at least two; psychomotor retardation will show one. Either outcome of the convergence index is reported as a finding; a high C_s for psychosis would favor the final common pathway account [12].

### 6.6 Interpretability deliverables

For each symptom: a ranked list of consensus modules, each with its support subgraph, compartment annotations, supporting perturbations by grade and the fraction of seeds in which it appeared; an overlap table against the curated modules in section 3.2.

## 7. Statistical analysis plan

Pre-registration (OSF) before Phase 3 of the symptom set, module scaffolding, splits, metrics, baselines and the primary endpoint. Primary endpoint: macro-averaged AUPRC over symptoms under the pathway-wise split for B6 versus the best of B1 to B3, five seeds by five folds, with a difference of at least 0.05 and a 95 percent bootstrap interval excluding zero. Proposed in version 0.4 and to be confirmed before pre-registration: the pathway-wise split is defined by Human-GEM subsystems rather than by the curated modules (section 6.1), and the comparison is on pooled held-out predictions paired with the permuted-label run of the same split. Secondary endpoints: disease-cluster and drug-wise grouped splits, time split, calibration, GWAS enrichment and the independence and convergence indices. Per-symptom tests are corrected by Benjamini-Hochberg. Failure criteria are stated in advance: if B6 does not beat B0 and B1 under the disease-cluster grouped split (the gene-wise split was shown to leak, section 6.1), the architecture is abandoned before the ablations are run.

## 8. Compute budget

Compute is a university cluster with GPUs whose binding limit is job wall time, not capacity. Graph size after joining layers: about 20,000 nodes and 100,000 typed edges (estimate). Route 1 training fits on one GPU in hours, so every training job checkpoints on a timer and on the Slurm pre-termination signal and requeues itself until a completion marker exists (slurm/train_resumable.sbatch); seeds, folds and the module-count candidates for the description-length selection are independent jobs submitted as arrays rather than loops inside one job. Route 2 requires one flux-sampling run per perturbed gene (about 3,000 genes in the reconstruction; minutes each, estimate); it is a Slurm array with one gene per task, idempotent so a partial array is resubmitted as-is, computed once and cached (slurm/flux_sampling_array.sbatch). The graph build and evidence assembly run on CPU partitions. No individual-level data and no secure enclave.

## 9. Phases, milestones and go/no-go criteria

Phase 0 (two weeks): freeze symptom set and crosswalk; data access inventory; pre-registration draft.

Phase 1 (four weeks): graph build with quality checks (connectivity, compartment counts, currency tagging); evidence assembly; counts per symptom and grade; proxy audit (A7). Go/no-go: at least 15 grade A or B perturbations for at least eight symptoms; otherwise drop symptoms. Status: E1 and E2 counts are measured (docs/phase1_counts.md) and ten symptoms pass; the graph build is done for the metabolic layer with subsystems and for the signaling, transcription and small-molecule layers when their files are present; the HPO expansion audit and the provenance-based grade A policy replace the OMIM lookup; the OnSIDES slice and the proxy audit remain.

Phase 2 (four weeks): baselines B0 to B5, evaluation harness, negative controls; metrics frozen. Status: B0, B1 and B3 run with per-fold and pooled metrics, the permutation control and three splits (grouped by gene or disease cluster, curated-module and subsystem hold-outs); B2, B4 and B5 remain.

Phase 3 (six weeks): B6 and ablations. Status: B6 and B3 run on the monogenic slice (section 12); the ablation arms are configured in experiments/run_main_model_batch.py.

Phase 4 (four weeks): orthogonal validation and equifinality tests.

Phase 5 (three weeks): mechanism cards, manuscript, release.

## 10. Repository layout

The layout below is committed. Files marked in the status block at the top are implemented and tested; the rest are stubs whose docstrings state their input files, output tables and the design section they serve.

```
mechanistic-pathway-learning/
  docs/
    experiment_design.md                      (this document)
    symptom_crosswalk.csv                     (section 3.1 as data)
    curated_pathway_modules.csv               (section 3.2 as data)
    evidence_grading.md
    preregistration.md
  configs/
    graph_build.yaml
    evidence_assembly.yaml
    model_noisy_or_pathway_modules.yaml
    baselines/*.yaml
    ablations/*.yaml
  data/
    raw/                                      (gitignored; downloaded by scripts)
    processed/                                (gitignored; graph and evidence tables)
  mechanistic_pathway_learning/
    graph/
      build_physiology_graph.py               (Human-GEM + OmniPath + CollecTRI + ChEMBL)
      tag_currency_metabolites.py
      brain_expression_weights.py
    evidence/
      load_monogenic_phenotype_annotations.py (HPO)
      load_drug_label_events.py               (SIDER, OnSIDES)
      load_literature_predications.py         (SemMedDB or the literature KG pipeline)
      assign_evidence_grades.py
      audit_diagnosis_proxies.py
    perturbation/
      flux_sampling_perturbation_features.py  (COBRApy)
      map_drug_targets_to_graph_nodes.py
    models/
      relational_message_passing_encoder.py
      noisy_or_pathway_module_model.py
      soft_constraint_losses.py
      minimum_description_length_selection.py
      baselines/
        popularity_baseline.py
        random_walk_with_restart_baseline.py
        knowledge_graph_embedding_baseline.py
        relational_gnn_sigmoid_baseline.py
        flux_feature_logistic_baseline.py
        language_model_zero_shot_baseline.py
    evaluation/
      perturbation_wise_and_pathway_wise_splits.py
      time_split.py
      ranking_and_calibration_metrics.py
      negative_controls.py
      gwas_gene_set_enrichment.py             (MAGMA wrapper)
      equifinality_independence_and_convergence.py
      mechanism_cards.py
  experiments/
    run_phase1_counts.py
    audit_hpo_term_expansion.py               (section 3.1 term audit)
    run_baselines.py                          (B0, B1; three splits; permutation control)
    run_main_model.py                         (B6 and B3; one split per invocation; resumable)
    run_main_model_batch.py                   (configurations and ablation arms; idempotent)
    aggregate_main_model_runs.py              (pooled and per-fold tables)
    analyze_pathway_modules.py                (supports, links, equifinality indices, curated overlap)
  tests/
```

## 11. Decisions taken and questions still open

Decided: no hard symptom-level constraints, with strong evidence as a 0.99 target or a 99-pseudo-observation anchor (A3); the earlier literature knowledge graph is not used and E3 comes from PubTator3, CTD and, with a license, SemMedDB; the reconstruction is Human-GEM with a Recon3D identifier map; compute is a runtime-limited GPU cluster and every long job is resumable. Decided in version 0.4 (reversible, recorded as data): diagnosis-level and off-construct HPO descendants are excluded through the crosswalk; grade A follows annotation provenance with frequency-scaled weights; the disease-cluster grouping is the default grouped split; modules read the perturbation difference field with sum pooling.

Open:
1. UMLS and OMIM access status (both free for academic use; OMIM registration is needed for the grade A bar).
2. Should the relieves relation be in version 1 or deferred? The scaffold supports both (separate link matrix per relation); including it doubles the evidence assembly work for drug indications.
3. Which language model and prompt fix baseline B5 before pre-registration.
4. Whether the cluster's wall-time limit is short enough that flux sampling should use fewer samples per gene (2,000 is the default) or a coarser flux variability summary alone.
5. Whether to relax the dominant-target rule from one to two ChEMBL targets, which would raise the pharmacological class above 65 drugs at the cost of less clean mechanism attribution.
6. Whether anhedonia and psychomotor retardation stay as genetics-only targets or are dropped from version 1.
7. The relief relation is thin at the symptom level in SIDER (one to ten drugs per symptom) because indications are recorded as diagnoses (major depressive disorder) rather than symptoms; decide whether diagnosis-level indications are admissible for the relieves relation only, which would reintroduce diagnosis labels on one side of the model.
8. HPO frequency as a loss weight (version 0.4) or as a soft target: the frequency estimates P(symptom | lesion), which is what the model outputs, so a soft target is the more literal use; the 0.99 positive target of A3 and a frequency target cannot both hold for the same pair.
9. Whether the pre-registered pathway-wise split uses Human-GEM subsystems or the curated modules, given that the curated modules hold 37 perturbations on the monogenic slice (section 6.1).
10. Whether "Progressive neurologic deterioration" and the other deterioration terms should stay excluded from cognitive impairment: they are filed under Mental deterioration in HPO, and the exclusion removes 13 genes reachable through no other term.

## 12. Reproduction status of this revision

Environment. The revision was produced in a cloud container whose network policy allowed GitHub and PyPI only. HPO (release 2026-09-01) and Human-GEM (2.0.1, SBML and yml) downloaded; the SIDER (sideeffects.embl.de), OmniPath and CollecTRI (omnipathdb.org), ChEMBL and UniChem (ebi.ac.uk), CTD, PubTator3 and Zenodo hosts were denied, so the pharmacological evidence class, the signaling and transcription layers and the small-molecule layer could not be rebuilt and every number in this section is for the monogenic slice on the metabolic layer: 451 genes, 861 gene-symptom pairs, 10 symptoms (elevated mood has one gene and is not scored), 24,185 nodes and 94,786 edges. The measured E2 counts of 2 October are kept in docs/phase1_e2_counts_measured.json so the Phase 1 table retains both classes.

Baselines (docs/phase2_baselines.md, docs/phase2_baselines_disease_cluster.md). Grouped by disease cluster, per-fold macro AUPRC: popularity 0.231, degree-scaled popularity 0.259, random walk with restart 0.293; on permuted labels 0.225, 0.273 and 0.284. Grouped by gene: 0.234, 0.258 and 0.347, with the random walk at 0.284 on permuted labels. The curated-module hold-out covers six modules and 37 perturbations. The subsystem hold-out covers 20 pathway subsystems and 182 perturbations (transport, exchange and isolated bins are not pathways and their genes are routed to their most frequent pathway subsystem); pooled macro AUPRC: popularity 0.185, degree-scaled popularity 0.290, random walk 0.296, against 0.156, 0.236 and 0.250 on permuted labels, so under a pathway hold-out the proximity baseline adds little beyond degree. Pooled AUROC of the constant predictor is 0.22 on this split, not 0.50, because held-out pathways have different base rates: pooled metrics under hold-outs are read against the permuted and popularity references, never against 0.5. Time split (cutoff 31 December 2015; 264 training pairs, 144 new positives among 3,793 scored pairs): macro AUPRC popularity 0.060, degree-scaled popularity 0.062, random walk 0.101 against 0.076 on permuted new-positive labels; macro AUROC of the random walk 0.547 against 0.514 permuted.

Proposed model and sigmoid baseline (docs/phase3_main_model.md). B6 (difference field, sum pooling, K = 8, description-length coefficient 1e-6, validation early stopping) and B3 (same encoder, sigmoid head) were trained on the disease-cluster split on CPU; the aggregated table is generated by experiments/aggregate_main_model_runs.py and the module analysis by experiments/analyze_pathway_modules.py. The numbers are in those files rather than here so that the design document does not drift from the runs.

## References

[1] F. Sedel et al., "Psychiatric manifestations revealing inborn errors of metabolism in adolescents and adults," J. Inherit. Metab. Dis., vol. 30, no. 5, pp. 631-641, 2007. PMID 17694356.
[2] C. Demily and F. Sedel, "Psychiatric manifestations of treatable hereditary metabolic disorders in adults," Ann. Gen. Psychiatry, vol. 13, p. 27, 2014. DOI 10.1186/s12991-014-0027-x.
[3] N. van de Burgt et al., "Psychiatric manifestations of inborn errors of metabolism: a systematic review," Neurosci. Biobehav. Rev., vol. 144, 104970, 2023. PMID 36436739.
[4] M. Walterfang et al., "The neuropsychiatry of inborn errors of metabolism," J. Inherit. Metab. Dis., vol. 36, no. 4, pp. 687-702, 2013. PMID 23700255.
[5] J. G. Thorp et al., "Genetic heterogeneity in self-reported depressive symptoms identified through genetic analyses of the PHQ-9," Psychol. Med., vol. 50, no. 14, pp. 2385-2396, 2020. PMID 31530331. DOI 10.1017/S0033291719002526.
[6] D. Cicchetti and F. A. Rogosch, "Equifinality and multifinality in developmental psychopathology," Dev. Psychopathol., vol. 8, no. 4, pp. 597-600, 1996. DOI 10.1017/S0954579400007318.
[7] J. Pearl, Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference. San Mateo, CA: Morgan Kaufmann, 1988.
[8] S. Bonner et al., "Implications of topological imbalance for representation learning on biomedical knowledge graphs," Brief. Bioinform., vol. 23, no. 5, bbac279, 2022. PMID 35880623.
[9] R. Celebi et al., "Evaluation of knowledge graph embedding approaches for drug-drug interaction prediction in realistic settings," BMC Bioinformatics, vol. 20, 726, 2019. DOI 10.1186/s12859-019-3284-5.
[10] E. Brunk et al., "Recon3D enables a three-dimensional view of gene variation in human metabolism," Nat. Biotechnol., vol. 36, no. 3, pp. 272-281, 2018. PMID 29457794. DOI 10.1038/nbt.4072.
[11] N. E. Lewis et al., "Large-scale in silico modeling of metabolic interactions between cell types in the human brain," Nat. Biotechnol., vol. 28, no. 12, pp. 1279-1285, 2010. PMID 21102456.
[12] O. D. Howes and S. Kapur, "The dopamine hypothesis of schizophrenia: version III, the final common pathway," Schizophr. Bull., vol. 35, no. 3, pp. 549-562, 2009. PMID 19325164.
[13] D. Turei et al., "Integrated intra- and intercellular signaling knowledge for multicellular omics analysis," Mol. Syst. Biol., vol. 17, no. 3, e9923, 2021. PMID 33749993.
[14] S. Muller-Dott et al., "Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities," Nucleic Acids Res., vol. 51, no. 20, 2023. PMID 37843125.
[15] J. L. Robinson et al., "An atlas of human metabolism," Sci. Signal., vol. 13, no. 624, eaaz1482, 2020. PMID 32209698.
[16] B. Zdrazil et al., "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods," Nucleic Acids Res., vol. 52, no. D1, 2024. PMID 37933841.
[17] M. A. Gargano et al., "The Human Phenotype Ontology in 2024: phenotypes around the world," Nucleic Acids Res., vol. 52, no. D1, 2024. PMID 37953324.
[18] M. Kuhn, I. Letunic, L. J. Jensen and P. Bork, "The SIDER database of drugs and side effects," Nucleic Acids Res., vol. 44, no. D1, pp. D1075-D1079, 2016. DOI 10.1093/nar/gkv1075.
[19] Y. Tanaka et al., "OnSIDES database: extracting adverse drug events from drug labels using natural language processing models," Med, 2025. DOI 10.1016/j.medj.2025.100642.
[20] H. Kilicoglu et al., "SemMedDB: a PubMed-scale repository of biomedical semantic predications," Bioinformatics, vol. 28, no. 23, pp. 3158-3160, 2012. PMID 23044550.
[21] J. G. Thorp et al., "Symptom-level modelling unravels the shared genetic architecture of anxiety and depression," Nat. Hum. Behav., vol. 5, 2021. PMID 33859377.
[22] C. A. de Leeuw et al., "MAGMA: generalized gene-set analysis of GWAS data," PLoS Comput. Biol., vol. 11, no. 4, e1004219, 2015. PMID 25885710.
[23] A. Ebrahim et al., "COBRApy: COnstraints-Based Reconstruction and Analysis for Python," BMC Syst. Biol., vol. 7, 74, 2013. PMID 23927696.
[24] L. Cowen et al., "Network propagation: a universal amplifier of genetic associations," Nat. Rev. Genet., vol. 18, no. 9, pp. 551-562, 2017. PMID 28607512.
[25] S. Bonner et al., "Understanding the performance of knowledge graph embeddings in drug discovery," Artif. Intell. Life Sci., vol. 2, 100036, 2022. DOI 10.1016/j.ailsci.2022.100036.
[26] C. Louizos, M. Welling and D. P. Kingma, "Learning sparse neural networks through L0 regularization," in Proc. ICLR, 2018. arXiv:1712.01312.
[27] A. P. Dawid and A. M. Skene, "Maximum likelihood estimation of observer error-rates using the EM algorithm," J. R. Stat. Soc. C (Appl. Stat.), vol. 28, no. 1, pp. 20-28, 1979. DOI 10.2307/2346806.
[28] V. C. Raykar et al., "Learning from crowds," J. Mach. Learn. Res., vol. 11, pp. 1297-1322, 2010.
