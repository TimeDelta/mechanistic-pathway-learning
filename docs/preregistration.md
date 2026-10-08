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

Confirmatory specification, 8 October 2026 (experiments/score_confirmatory.py implements it). It replaces the primary
endpoint and failure criterion at the top of this file; the reasons are given with each item. Full-graph training before
it: two message-passing pilots interrupted after one epoch, on 3 October (noisy-OR head, the full graph and evidence of
that date, its command not logged) and 6 October (sigmoid head, graph_full with evidence_full), which printed validation
macro AUPRC 0.366 and 0.310 on perturbations some of which are now in the lockbox, and timing runs of five training
steps on a development fold (experiments/run_main_model.py --timing-batches), which score nothing. No full-graph model
has scored a test set. Decisions taken with the user on 7 and 8 October: both heads, both encoders, the rewiring
control, macro and micro AUPRC both required, a degree control, a lockbox, five seeds, and no flux sampling in the
confirmation.

- Split. The grouped split by disease_cluster_and_targets with a lockbox (configs/lockbox_v1.json, SHA-256 ebd851a8…,
  drawn by experiments/draw_lockbox.py with seed 20261007): 317 of 1,539 perturbations (253 genes, 64 drugs, 27 of them with
  a kept positive), 481 kept positive pairs, 18 symptoms with five or more kept positives inside it. The draw takes whole
  leakage groups, per stratum (groups holding a drug with a kept positive, and the rest), in a seeded order, until each
  stratum's lockbox reaches 20 percent of its perturbations; the one group above 5 percent of all perturbations
  (cluster:ABCA3, 232 perturbations) stays in development. The pathway-wise subsystem split of the skeleton becomes a
  secondary reading on the development set: Human-GEM subsystems hold metabolic genes only, so it cannot test drugs or the
  signalling layers.
- Models (experiments/run_main_model_batch.py): confirmatory_message_passing_noisy_or, confirmatory_message_passing_sigmoid,
  confirmatory_linear_response_noisy_or and confirmatory_linear_response_sigmoid. Every one reads graph_full_neuronal
  (nodes 48ebf531…, edges 78cf6c07…), evidence_full_v2 (228a4afe…), the better_v1 selection of it (c691a9aa…) and the node
  descriptors (14033cb3…: protein, metabolite and reaction properties, brain region and cell-class expression). The
  message-passing encoder reads typed node features with the descriptors. The linear-response encoder adds the cofactor
  relations and propagates once per cell class (cell_class_weights 234224d8…: the ten HPA single-nucleus classes, the
  dopaminergic class and an all-cells class), with the 1,684 extracellular metabolites one pool shared by the classes other
  than all cells. The noisy-OR head starts its leaks at the training base rate, its module biases at -3 and its gate noise at
  0.5, with link, leak and gate learning rates 0.02, 0.0002 and 0.05; the sigmoid head uses the trainer defaults. Asymmetry:
  the message-passing encoder sees cell-class expression only as node descriptors, not as channels.
- Training. Each confirmatory run trains on the 1,222 development perturbations, of which a grouped 15 percent is the
  validation set for early stopping (validation loss, patience 8, at most 60 epochs), and scores the lockbox once
  (--lockbox configs/lockbox_v1.json --score-lockbox). Seeds 0 to 4. Per seed and model, two control runs: labels permuted
  within degree quintiles and within the lockbox and the development set apart (rows move whole with their weights and
  label mask; --permute-labels), and the graph rewired within each relation with 50 attempted swaps per edge, no self-loops
  and no repeated edges (--rewire-swaps-per-edge 50; 32 seconds; 0.5 percent of edges stay in place). A model missing any of
  its fifteen runs is not confirmatory.
- Baselines. popularity, degree_popularity, random_walk_with_restart and knowledge_graph_embedding_transe, fitted on the
  development set and scored on the lockbox for each seed (run_baselines.py --lockbox ... --score-lockbox --seed k
  --with-kg-embedding --rewiring-method integer_draws), which also scores them on that seed's permuted labels. The best
  baseline of a reading is the one with the highest seed-mean lockbox score, chosen without reference to any model. The
  relational GNN with sigmoid head, a comparator in the skeleton, is now one of the four tested models; popularity and
  degree_popularity are added because they score highest on the full graph (docs/phase2_baselines_full_v2.md).
- Readings (seed mean of model minus best baseline on the lockbox rows): macro AUPRC over the 18 lockbox symptoms; micro
  AUPRC over every symptom with five or more kept positives in the whole data (20); both after ranking every score inside
  degree strata; the permutation difference in differences (the advantage on the real labels minus the advantage when model
  and baseline are trained and scored on permuted labels); the rewiring difference (the model minus the same model on its
  rewired graph).
- Hypotheses and decision. H1 (prediction) for a model: macro, micro, both within-strata and both permutation readings
  above zero, and the macro and micro differences at least 0.05. H2 (graph content): H1 and both rewiring readings above
  zero. Each reading's one-sided p-value comes from a paired bootstrap over lockbox perturbations (4,000 resamples, the same
  rows for every model, seed and baseline); a hypothesis's p-value is the largest of its readings' (an intersection-union
  test, Berger 1982). Holm's procedure (Holm 1979) runs over the eight hypotheses at one-sided 0.025. A model is confirmed
  for H1 or H2 when Holm rejects that hypothesis and its macro and micro differences reach 0.05.
- Failure criterion, replacing the one at the top: if no model is confirmed for H1, the architecture has not been shown
  to predict better than the baselines. The ablations still run (amendment of 7 October).
- Secondary readings, outside the Holm family and reported with paired intervals: noisy-OR minus sigmoid within each
  encoder, linear response minus message passing within each head, per-symptom AUPRC, the subsystem split and every
  ablation on the development set.
- Order of work. Pilots of any configuration run on the development set only (--lockbox without --score-lockbox; runs land
  in <configuration>_<grouping>_development). If a pilot shows a confirmatory configuration failing to train (a defect, not
  a low score), the configuration may be amended here, dated, before its first lockbox run. The lockbox is scored once by
  experiments/score_confirmatory.py, which writes a SCORED marker; a second scoring is listed in its output.
- Precision. On the full graph a 20 percent grouped hold-out (one of five folds, 308 perturbations, 382 to 548 kept
  positive pairs, 13 to 20 symptoms with five positives) gave paired-bootstrap 95 percent half-widths for differences
  between two baselines of 0.019 to 0.063 in macro AUPRC (median 0.041) and 0.018 to 0.041 in micro (median 0.028)
  (experiments/estimate_holdout_precision.py). At Holm's smallest level, 0.025 / 8, the half-width grows by about 1.4
  (2.73 / 1.96 standard errors), so a macro difference near 0.06 and a micro difference near 0.04 are the smallest the
  lockbox can confirm. Two caveats: a model and a baseline may co-vary less than two baselines do, which widens the
  interval, and the seeds are averaged, not resampled, so seed variance is not in it.
- Compute. On this 4-core machine, at 3 threads and batch 16, a training step on a development fold takes 1.8 s
  (linear response, sigmoid), 2.1 s (linear response, noisy-OR), 2.5 s (message passing, sigmoid) and 2.7 s (message
  passing, noisy-OR), after both encoders were changed on 8 October to multiply only the 54,982 stacked adjacency rows that
  hold an edge instead of all 700,435 (identical output, tested; before the change 15.8, 16.9, 6.9 and 7.0 s). With the
  slice's epochs to early stopping (about 15 for the sigmoid head, about 35 for noisy-OR), the 60 confirmatory runs take
  about 70 hours. Run order: one development pilot per configuration (fold 0, seed 0) to check that it trains; the
  development baselines; the 60 lockbox runs, two at a time with 2 threads each; then the lockbox baselines for seeds 0 to 4
  (their lockbox scores are not read before the models finish) and score_confirmatory.py.
- Flux route. Not in the confirmation: measured on Human-GEM, one LP takes 0.12 seconds with every boundary open, a flux
  variability or sampling warm-up per knockout about 0.9 core-hours; the perturbations a Human-GEM knockout can represent
  hold 14 percent of the kept pairs (235 genes, 15 drugs); no brain medium is defined; and a drug is not a knockout. It stays exploratory.

Amendment, 8 October 2026, written before any lockbox run and before the slice runs it refers to had finished: node
descriptors. On the slice, adding the node descriptors without brain expression lowered both encoders' macro and micro
AUPRC in all six readings of experiments/compare_twin_runs.py (per fold, pooled and within degree strata). Message
passing (b3_typed_nodes_descriptors against b3_typed_nodes): macro per fold -0.030 [-0.056, -0.004], within degree strata
-0.018 [-0.039, -0.002], micro pooled -0.018 [-0.035, +0.003]; with the descriptors it scored the same as the
descriptors-only control without a graph (+0.002 [-0.023, +0.027] per fold). Linear response
(b3_linear_response_cofactors_descriptors against b3_linear_response_cofactors): micro pooled -0.025 [-0.043, -0.008],
micro per fold -0.027 [-0.051, -0.002], macro pooled -0.003 [-0.016, +0.014]. The confirmatory configurations read the
descriptors with brain region and cell-class expression, which the slice had not tested, so two slice arms now test them
(b3_linear_response_cofactors_descriptors_brain and b3_typed_nodes_descriptors_brain). Rule: for each encoder, take the six
differences of its brain-expression arm minus the same encoder without descriptors (b3_linear_response_cofactors;
b3_typed_nodes). If their mean is below zero, that encoder's two confirmatory configurations drop --node-descriptors, their
development pilots are moved aside and rerun before the first lockbox run, and the linear-response encoder keeps its
cell-class channels; otherwise they stay as specified. The rule uses the point estimates, not the intervals, because the
slice cannot resolve differences of this size (half-widths of 0.015 to 0.04) and an input that does not help on
development data adds parameters to a model trained on about 1,000 perturbations. Limitation: the slice is the metabolic
graph, without the signalling, receptor and electrical layers of graph_full_neuronal. Held on 8 October: the user is
deciding how the descriptors enter, so the rule is reported when its runs finish and applied only on their confirmation.

Amendment, 8 October 2026, before any lockbox run (the user's decision): the two permutation readings leave H1 and become
secondary readings, reported with their intervals and p-values outside both hypotheses. H1 is now the macro, micro and
both within-strata readings above zero with the macro and micro differences at least 0.05; H2 is H1 and both rewiring
readings above zero; Holm still runs over the eight hypotheses. Reason: the difference in differences subtracts a second
noisy difference, so if the two are nearly uncorrelated its standard error is about 1.4 times that of the plain
difference and it would decide H1 (a projection, not a measurement). The within-strata readings and the degree_popularity
baseline remain the degree control. experiments/score_confirmatory.py implements it. The lockbox runs start only on the
user's go-ahead.

To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
