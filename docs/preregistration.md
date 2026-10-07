# Pre-registration skeleton (to be completed before Phase 3 and posted to OSF)

Frozen in Phase 0 and 1: symptom set and HPO term exclusions (docs/symptom_crosswalk.csv, audited in
docs/hpo_term_audit.md), curated module scaffolding (docs/curated_pathway_modules.csv), pinned data
versions (docs/data_sources.md), evidence grading policy and weights (configs/evidence_assembly.yaml:
grade A by curated-synopsis provenance with frequency scaling), splits and metrics (configs/evaluation.yaml).

Primary hypothesis: see experiment_design.md section 6.5 (equifinality independence and
convergence indices; sufficiency test by module ablation).

Primary endpoint: macro-averaged AUPRC under the pathway-wise split, proposed model versus the
best of network proximity, knowledge-graph embedding and relational GNN with sigmoid head;
minimum difference 0.05; 95 percent bootstrap interval excluding zero; five seeds by five folds.
Proposed in design version 0.4 and to be confirmed here: the pathway-wise split is the Human-GEM
subsystem hold-out (pathway subsystems only; transport, exchange and isolated bins are not
pathways), scored on pooled held-out predictions and paired with the permuted-label run of the
same split; the curated modules hold 37 perturbations on the monogenic slice and serve the
interpretability overlap instead.

Failure criterion: if the proposed model does not beat the popularity and network proximity
baselines under the disease-cluster grouped split, the architecture is abandoned before ablations
run. (Design version 0.3 named the gene-wise split; it leaks through genes sharing a disease,
design section 6.1.)

Negative controls run with every split: labels permuted within degree strata; degree-preserving
rewiring; grade shuffle; peripheral-event control.

Precision of the primary comparison on the monogenic slice (451 perturbations, disease-cluster
grouped split): the paired-bootstrap 95 percent half-width of a pooled macro AUPRC difference between
two models is about 0.015 to 0.02 (docs/phase3_main_model.md), so the minimum difference of 0.05 is
detectable there; on the curated-module hold-out (37 perturbations) it is not, which is why the
subsystem definition of the pathway split is proposed above. The full data (drug perturbations and the
signaling layers) will shift these widths and the statement is to be recomputed before posting.

Amendment of 7 October 2026, written after the monogenic-slice results were seen. Nothing here was posted
before those results, so the amendment is disclosed as made after them (Nosek et al. 2018,
doi:10.1073/pnas.1708274114, on preregistration when "the data are preexisting").

- The monogenic slice is a pilot. It has one layer of the graph (the metabolic layer), 451 genes, no brain
  expression and no drug perturbations, and on it even a readout trained on the random walk's own distributions
  does not beat the untrained walk (docs/trained_diffusion_readout.md), so the slice cannot separate a weak mechanism
  from too few labels. Its result, that no trained configuration beats the random walk, does not trigger the failure
  criterion above, and no component (the linear-response encoder, the pathway modules, the node descriptors) is
  dropped on its strength. This is the user's decision of 7 October: the mechanistic encoder and the node properties
  it reads stay until the full graph with all layers and better (not more) training examples has been tested.
- The confirmatory test is the full graph with every layer built at the time of the run (graph_full_neuronal on
  7 October: Human-GEM, OmniPath, CollecTRI and the Reactome neuronal and oxidative layer), the node descriptors and
  the brain region and cell-class expression, under the disease-cluster grouped split. The failure criterion applies
  there unchanged.
- What makes a training example "better" is fixed and dated in this file before the first full-graph model is
  scored, from evidence properties only (grade, frequency qualifier, provenance, layer contact), never from any
  model's predictions.
- The ablations (node descriptors, brain expression, edge signs, relation typing) run whatever the primary result,
  because they are what says which part carries the signal (confirmed by the user on 7 October). The with and without
  descriptors comparison is one of them; the descriptors are not dropped before the final experiment.
- The confirmatory configuration is the full graph with the node descriptors, the brain region and cell-class
  expression, the cell-class propagation channels (including a dopaminergic class) and the better training examples.
  Full-graph runs made before all of these exist are full-graph pilots and are reported as such.

Better training examples, fixed on 7 October 2026 before any full-graph model was scored (experiments/build_label_selection.py;
data/processed/label_selection/better_v1_full.parquet, selection version better_v1):

- A gene pair (grade A, an HPO annotation) is kept when its frequency is at least 0.30, the lower bound of the HPO
  frequency term Frequent (HP:0040282, present in 30 to 79 percent of cases). Pairs annotated Occasional or Very rare,
  pairs with a counted frequency below 30 percent and pairs with no frequency are set aside.
- A drug pair (grade B, SIDER or OnSIDES label) is kept when its label frequency is at least 1 percent, the lower bound of
  the CIOMS frequency category common, or when both SIDER and OnSIDES list it (two separate label extractions). The rest,
  most of them label statements with no frequency, are set aside.
- A pair set aside is masked, not relabelled: it is neither a positive nor a negative in the loss or in any metric
  (load_experiment_data, label_mask). Unobserved pairs stay unlabelled as before.
- On evidence_full this keeps 1,735 of 3,260 positive pairs: 1,395 gene pairs over 844 genes and 340 drug pairs over
  59 drugs. Elevated mood or mania keeps 11 positive pairs and psychomotor retardation keeps 1, so psychomotor
  retardation leaves the macro average (five positives are needed) and elevated mood or mania is scored on few positives.
- Layer contact was considered and not used, because it does not select on graph_full_neuronal: 2,286 of 2,425 gene pairs
  and all 876 drug pairs already touch a non-metabolic layer.
- Every full-graph run, baseline and comparison is trained and scored under this selection. A run trained on another
  selection is refused by the scoring scripts (aggregate_main_model_runs.py, compare_twin_runs.py), which compare the
  SHA-256 of the selection file (ffb06fb4… for better_v1_full). Before this rule, training already weighted each gene
  positive by its frequency (mean weight 0.64 over the 796 gene pairs below 0.30, 1.0 at 0.30 or above and 1.0 with no
  frequency), and every pair counted fully in scoring; the selection removes the set-aside pairs from both.

Symptom list expanded on 7 October 2026 at the user's request, before any full-graph model was scored (docs/symptom_crosswalk.csv,
audited in docs/hpo_term_audit.md; data/processed/evidence_full_v2; the better_v1 rule above applied to it:
data/processed/label_selection/better_v1_full_v2.parquet, SHA-256 c691a9aa…). This replaces evidence_full and better_v1_full
(ffb06fb4…) for every full-graph run; the slice keeps its own evidence table.

- Frame: a symptom is added when it is an item or domain of a standard instrument (PHQ-9, GAD-7, the Neuropsychiatric
  Inventory, HAM-D item 14, Y-BOCS) and has an HPO term, MedDRA preferred terms or both. Diagnosis-level HPO terms stay
  excluded (assumption A7).
- New symptoms: apathy, suicidality, self_injury, compulsive_behavior, disinhibition_or_impulsivity, hyperactivity,
  emotional_lability, increased_appetite, abnormal_dreams, decreased_libido and catatonia. Each excluded descendant is named in
  the crosswalk with its reason. Two were found by looking at which genes a term supplies: Eye poking (all its annotated genes
  are Leber congenital amaurosis genes; it is the oculo-digital sign of blind children) and Inappropriate laughter (mostly the
  Angelman syndrome sign: UBE3A, SLC9A6, MECP2, CDKL5).
- Existing symptoms gained the MedDRA preferred terms drug labels use for the same construct, several of which the HPO side
  already held: Akathisia (psychomotor agitation; HPO Akathisia sits under Restlessness), Delirium, Disorientation and the
  amnesias (cognitive impairment; HPO Delirium sits under Confusion), Paranoia and the hallucination subtypes (psychosis, with
  HPO Paranoia added), Hostility and Anger, the insomnia subtypes, Sedation, Dysphoria and Elevated mood. Tension was added to
  anxiety on the reading that MedDRA files it with the anxiety symptoms; that placement was not checked against the MedDRA
  hierarchy, which the project does not hold.
- Not added: decreased appetite (on drug labels it follows nausea and other gastrointestinal effects); autism and attention
  deficit hyperactivity disorder (diagnosis-level); intellectual disability and speech delay (developmental, not a symptom
  state); tics (a motor disorder); lethargy and drowsiness (in HPO, signs of reduced consciousness); sleep apnoea.
- Counts: 4,439 positive pairs, 2,211 kept (1,793 gene pairs over 925 genes; 418 drug pairs over 61 drugs), against 1,735
  before. Kept pairs by symptom: cognitive impairment 373, fatigue 331, irritability or aggression 273, anxiety 254,
  depressed mood 189, psychomotor agitation 89, psychosis 85, somnolence or hypersomnia 84, hyperactivity 83, self-injury 75,
  emotional lability 66, disinhibition or impulsivity 62, insomnia 57, apathy 54, compulsive behaviour 46, increased appetite 37,
  abnormal dreams 17, suicidality 13, elevated mood or mania 11, decreased libido 11, psychomotor retardation 1, catatonia 0.
- Scored symptoms: every symptom with at least five kept positive pairs in the whole data (so psychomotor retardation and
  catatonia are not scored, and anhedonia has no pair). The macro average takes the symptoms with at least five positives in
  the test set, as before; micro AUPRC (below) takes every scored symptom.
- The pipeline rebuilt from the previous crosswalk reproduces evidence_full and better_v1_full exactly (selection SHA-256
  ffb06fb4… both times), so the changes come from the crosswalk alone.
- What the expansion changes: most of the added kept pairs are gene pairs, and four new symptoms (self-injury, hyperactivity,
  compulsive behaviour, disinhibition) come mostly from neurodevelopmental syndromes with no or few drug positives. Drug kept
  pairs rose from 340 to 418.

Micro AUPRC, added on 7 October 2026 at the user's request (ranking_and_calibration_metrics.micro_auprc): average precision over
every labelled perturbation-symptom pair ranked in one list. It keeps the symptoms the macro average leaves out in a test set
with fewer than five positives. It weighs frequent symptoms more and credits ranking symptoms by base rate, so it is compared
against the same baselines (popularity among them). Every scoring script reports it beside the macro average, in the same
three readings (per fold, pooled, within degree strata).

Leakage grouping, decided 7 October 2026 (experiments/check_drug_target_leakage.py, docs/drug_target_leakage.md): under the
disease-cluster grouping a drug is grouped by its target set, so it can sit in another fold than the loss of function of its
own target gene. On the expanded labels, 69 of the 73 drugs that target a labelled gene had such a gene in another fold, and 66
of the 418 kept drug positives were also kept positives of such a gene. Where a labelled target gene is positive for a symptom
the drug is positive in 67 of 242 pairs (0.28), against 0.17 over all labelled drug pairs. Every full-graph run therefore uses
the grouping disease_cluster_and_targets, which holds a drug out with every labelled gene it targets and every drug sharing a
target node. That closes the overlap (0 drugs, 0 shared positives) and keeps the fold sizes at 307 to 308 perturbations; the
largest group grows from 140 to 232 perturbations.

To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; number of flux samples per gene; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
