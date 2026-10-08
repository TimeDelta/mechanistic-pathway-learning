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
- Training. Each confirmatory run trains on the 1,222 development perturbations, of which one grouped seventh is the
  validation set for early stopping (fold 0 of round(1 / 0.15) = 7 grouped folds; validation loss, patience 8, at most 60
  epochs; as amended on 8 October, drawn from the groups no larger than half of it, 142 perturbations, and followed by a
  refit on all 1,222 for the chosen number of epochs, which is the model scored), and scores the lockbox once
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

Amendment, 8 October 2026, before any lockbox run and before the slice runs it names had started: how the node
descriptors enter. The user keeps the descriptors in both encoders ("I really don't want to drop them") and asked how to
stop them overwhelming the graph, so the drop branch of the descriptor rule above is withdrawn: neither encoder's
confirmatory configurations drop --node-descriptors. The rule's readings are still reported (docs/descriptor_rule.md).
On 8 October they were: linear response, mean of the six differences +0.002 (macro per fold -0.003, pooled +0.006,
within strata +0.007; micro -0.004, -0.001, +0.007); message passing -0.018 (macro -0.022, +0.002, -0.007; micro -0.040
[-0.057, -0.022], -0.033 [-0.052, -0.012], -0.006). The message-passing reference, b3_typed_nodes, was trained before the
structural node features gained the two GTEx brain expression columns (commit cb5cf70, 3 October; 16 structural columns
against 18), so that comparison, and the message-passing readings of the descriptor amendment above, also change the
structural features. b3_typed_nodes is retrained under the current code (runs/encoder) and the message-passing readings
are reported again.

Why a treatment: scored on a degree-preserving rewiring of its graph at test time (experiments/evaluate_on_rewired_graph.py,
50 attempted swaps per edge, seed 0, the trained weights unchanged), the linear-response encoder without descriptors
(b3_linear_response_cofactors) loses 0.026 ± 0.019 macro AUPRC over five folds and with the brain-expression descriptors
0.001 ± 0.010 (micro 0.004 ± 0.019 and -0.002 ± 0.005), at about the same AUPRC on the real graph; the descriptors-only and local-structural controls lose 0.000, as
a model without edges must. The descriptors replace what the edges gave instead of adding to it. In message passing with
descriptors the best validation epoch is 2 to 5 on the slice and 0 in the full-graph development pilot, and at the best
epoch the descriptor weights are at their initial scale (root mean square 0.049 against 0.046 for the initial uniform
draw), so the descriptors enter as a fixed random projection that tells the nodes apart.

Treatments (mechanistic_pathway_learning/models/descriptor_treatments.py, --descriptor-treatment): seed_masked, under
which each perturbation's perturbed nodes read their structural columns only, in the message-passing base state (in the
perturbed pass and in that perturbation's own unperturbed reference) and in the linear-response output gate, while every
other node keeps its descriptors; and zero_init_slow, under which the descriptor columns get their own linear map,
initialised at zero and trained at 0.0002, a tenth of the main learning rate. Adam moves a parameter by about one
learning rate per step, so on the slice (20 steps per epoch) a descriptor weight moves at most about 0.02 in five epochs;
a weight penalty was not used because under Adam the step does not scale with the gradient (Loshchilov and Hutter, ICLR
2019: L2 regularization and weight decay "are equivalent for standard stochastic gradient descent (when rescaled by the
learning rate), but as we demonstrate this is not the case for adaptive gradient algorithms, such as Adam") and a
penalty acts through the optimum, which early stopping does not reach. With plain (the default) the encoders are
unchanged: the same parameters, random draws and outputs, bit for bit, as before the treatments existed.

Slice arms (runs/encoder, five disease-cluster folds each): b3_typed_nodes_descriptors_brain_seed_masked,
b3_typed_nodes_descriptors_brain_zero_init_slow, b3_linear_response_cofactors_descriptors_brain_seed_masked and
b3_linear_response_cofactors_descriptors_brain_zero_init_slow, each one change from its *_descriptors_brain arm (the
plain arm). Readings per arm: graph reliance, the mean over folds of the macro and micro AUPRC lost on the rewired graph,
averaged over the two (docs/graph_reliance.md); and the six differences of experiments/compare_twin_runs.py against the
plain arm (docs/descriptor_treatments.md). Rule, per encoder: a treatment qualifies if its graph reliance is larger than
the plain arm's and the mean of its six differences against the plain arm is at least -0.01. The qualifying treatment
with the larger graph reliance replaces plain in that encoder's two confirmatory configurations; if none qualifies,
plain stays. Point estimates, as in the descriptor rule, because the slice does not resolve differences of this size.
The margin of -0.01 is below the slice's resolution on purpose: the treatments exist to move the work from the
descriptors to the graph, which the user asked for, and the rule rejects only a treatment whose loss is large enough to
show. As with the descriptor rule, the outcome is reported and applied only on the user's confirmation
(docs/confirmatory_runbook.md, step 2b); a change means new development pilots of the changed configurations before the
first lockbox run. Cost: seed_masked propagates a reference per perturbation, so a message-passing step takes about 2.4
times as long (3.8 against 1.55 seconds on the slice), and so would the 30 message-passing lockbox runs; in the
linear-response encoder it costs one more gate per perturbation. Limitations: test-time rewiring measures how much a
trained model depends on its edges, not whether the real edges beat degree-matched ones (the retrained rewiring control
of H2 asks that); one rewiring seed; the slice is the metabolic graph only; and under seed_masked a perturbation's
neighbours still carry descriptors that describe it indirectly (a gene's reactions carry its EC numbers).

Amendment, 8 October 2026, before the arms it names had started: descriptor blocks per node type. The user suggested
removing the descriptors of some node types, genes first, because a gene's descriptors are largely the graph around it
(which reactions it catalyses and which proteins it binds), and asked for the search to continue on the slice until it
is no longer worth it, with full-graph tests where they say more. --drop-descriptor-blocks leaves whole blocks out
(graph/node_descriptors.py, DESCRIPTOR_BLOCK_PREFIXES: metabolite, reaction_ec, reaction_brain, protein, gene_brain). Three
slice arms per encoder, each one change from its plain *_descriptors_brain arm: without_gene (protein and gene_brain out,
the user's suggestion), without_gene_derived (also reaction_brain, the genes' expression carried onto their reactions by
the gene rule, so a gene seed's own reactions no longer carry its expression pattern) and without_protein (only the
protein function block out; the genes keep brain expression, which the graph does not encode for the message-passing
encoder). On the slice every seed is a gene, so without_gene removes what seed_masked removes there and more, at no extra
compute. These arms join seed_masked and zero_init_slow under the treatment rule of the previous amendment (graph reliance
above the plain arm's, mean of the six differences against the plain arm at least -0.01, the larger graph reliance among
those that qualify). Arms the iteration adds later are dated here before they run, and every arm tried is reported in
docs/graph_reliance.md and docs/descriptor_treatments.md, chosen or not: selecting among slice arms makes the slice
estimate of the chosen arm optimistic, which does not reach the lockbox, since the lockbox is scored once by a model
fixed before it is opened.

Amendment, 8 October 2026, before any lockbox run (the user's decisions): the test, the micro floor, a no-descriptor
branch and graph paring.
- Test. Each of the four models is a candidate and is tested on its own at one-sided 0.025: H1 first, and H2 only if
  that model's H1 is confirmed, at the same level. Holm's procedure over the eight hypotheses is withdrawn; nothing is
  divided across models. Consequence, stated beside the results by experiments/score_confirmatory.py: with four models
  and no correction across them, the chance that at least one is confirmed by luck is above 0.025, at most
  1 - 0.975^4 = 0.096 if the four were independent and less since they share the data, the lockbox and the baselines.
- Floors. The macro difference must be at least 0.041 and the micro difference at least 0.028, the median projected
  95 percent half-widths of a macro and a micro difference on a 20 percent hold-out (the precision statement above). The
  first version of this amendment kept the macro floor at 0.05; the user moved it to 0.041 the same day, before any
  lockbox run, so that both floors are the same kind of bar: what the lockbox can detect, not a judged smallest effect
  worth having (a micro difference is the less noisy of the two, so one value of 0.05 would have been a different bar for
  each). At one-sided 0.025 a difference at its floor has a 95 percent interval that just reaches zero
  when the realised width equals the projection, so a floor binds only when the lockbox interval comes out narrower than
  projected and the floors add little to the significance test; the effect sizes and intervals reported next to every
  decision carry the question of size. With Holm withdrawn, the 0.06 and 0.04 of the precision statement (Holm's smallest
  level) no longer apply: at one-sided 0.025 the smallest differences the lockbox can confirm are near the half-widths.
- Ablation branch. The four configurations without --node-descriptors (confirmatory_*_without_descriptors in
  experiments/run_main_model_batch.py) are an ablation, run as development pilots (runs/full/ablation_pilots_without_descriptors.sh)
  and reported against the descriptor configurations; they are not among the tested models.
- Graph paring. Pared graph variants are tested on the development set only if the development pilots do not beat the
  development baselines. After the lockbox has been scored a pared graph can only be an exploratory result, because the
  lockbox is not scored twice.

Amendment, 8 October 2026, before any lockbox run: the early-stopping validation set (a defect found by diagnostics of
the development pilots; the user allowed a change of the stopping rule the same day, if needed and if differences in
fold prevalence could be shown not to interfere).
- Defect. The trainer takes as its early-stopping validation set fold 0 of a grouped split of the training pool into
  seven folds, and assign_grouped_folds places the largest group first, in fold 0. With the lockbox removed the
  development data hold one leakage group of 232 of 1,222 perturbations (disease_cluster_and_targets joins 169 genes and
  63 drugs through shared disease clusters and drug targets). It holds 27 percent of the development positives and 61 to
  100 percent of the positives of seven symptoms (abnormal dreams, decreased libido, elevated mood or mania, insomnia,
  psychomotor retardation, somnolence or hypersomnia, suicidality). A lockbox run's training pool is every development
  perturbation, so for every seed its validation set would have been that group alone and no confirmatory model would
  have trained on it, while the baselines do. Development folds 1 to 4 do the same. Development fold 0, where the pilots
  ran, has the group as its test set, so the pilots were not affected; their test fold is atypical for the same reason
  (rich in drugs and in those seven symptoms), and the development baselines are compared on that fold.
- Change. The eight configurations of the confirmatory family (the four tested models and the four ablations) add
  --keep-large-groups-in-training (experiments/run_main_model.py, early_stopping_validation): a group larger than half
  the expected validation set (0.15 x pool / 2, 92 perturbations in a lockbox run) stays in training and the seven folds
  are drawn over the other groups. In a lockbox run this gives 1,080 training and 142 validation perturbations, the
  largest validation group 58 and 14.5 to 16.1 percent of the positives in validation (seeds 0 to 4), against 990, 232, 232
  and 27 percent. On development fold 0 the split is identical (no group in that pool exceeds 73), so the four pilots
  stand for the amended configurations. Other runs keep the earlier split; on the slice the largest group (27
  perturbations) is about half of each validation set but holds 14 percent of the positives.
- Precision unchanged. The lockbox is unchanged (configs/lockbox_v1.json, 317 perturbations); what shrank is the
  validation set, which is drawn from the development perturbations and scored for no hypothesis. The precision
  statement's five 20 percent folds were drawn over all 1,539 perturbations, and its fold 0 holds the group of 232 (548
  kept positive pairs and 20 scorable symptoms, the most of any fold). The lockbox cannot hold that group, since the
  group lies wholly in the development data, so the four folds without it are the closer match: their medians are 0.040
  macro (range 0.021 to 0.063) and 0.027 micro (0.018 to 0.038) over twelve baseline-pair readings, against 0.041 and
  0.028 over all fifteen. The floors stay at 0.041 and 0.028.
- The stopping metric stays the validation loss. Validation macro AUPRC was considered and not adopted. (i) Not needed:
  the message-passing sigmoid pilot, which the loss stopped at epoch 0, scores macro 0.142 and micro 0.233 on its
  development test at its last state (epoch 8), against 0.131 and 0.257 at the selected state. (ii) Too noisy to choose
  an epoch: on its 140 validation perturbations the paired 95 percent bootstrap interval of the validation macro AUPRC
  difference between those two states is [-0.052, +0.047], wider than the range of validation macro AUPRC over all nine
  epochs (0.036). (iii) Prevalence does interfere. Within a fold every epoch is scored on the same validation
  perturbations, so a symptom's validation prevalence sets the same chance level for every epoch and cancels from the
  comparison, and a per-symptom AUPRC depends only on the ranking within the symptom, so learning a base rate does not
  move it. But prevalence decides which symptoms the validation set can score at all (five positives): 12 on fold 0
  against 18 on its test set, 10 in common, and 11 or 12 in a lockbox run under the amended split. An AUPRC rule would
  stop on a symptom set that leaves out about a third of what the test scores; the loss counts every labelled pair.

Amendment, 8 October 2026, before any lockbox run (the user's objection the same day: the baselines fit on every
development perturbation while the models leave out their early-stopping validation set, so a null result would be
less meaningful): refit after early stopping.
- The eight configurations of the confirmatory family add --refit-on-validation (experiments/run_main_model.py). The
  early-stopping run on the validation set (142 perturbations in a lockbox run) only chooses the number of epochs; a
  fresh model built from the run's seed is then trained on the training and validation perturbations together (all
  1,222 development perturbations in a lockbox run, 977 in a development pilot) for best epoch + 1 epochs, and its test
  predictions are the ones scored. The baselines fit on the same perturbations, so models and baselines now fit on the
  same rows; before, the models fitted on 1,080 of 1,222 (12 percent fewer), which biased a comparison against them and
  would have made a null result partly a cost of the held-out data. This is the second-pass strategy of Goodfellow,
  Bengio and Courville (Deep Learning, 2016, section 7.8, algorithm 7.2). The refit keeps the number of epochs, not the
  number of optimizer steps (each refit epoch has about 13 percent more steps), a choice that book leaves open; the
  overfitting early stopping guards against comes from repeated passes over the same perturbations.
- Permuted-label and rewired-graph runs are refitted the same way, on their own labels and graph. The early-stopped
  model's test predictions are kept (test_predictions_early_stopped.npy) and are not scored. experiments/score_confirmatory.py
  accepts a lockbox run only if it carries the refit and --keep-large-groups-in-training.
- The four development pilots are refitted from their stored early-stopping checkpoints (runs/full/refit_pilots.sh;
  a run that finished without the refit gets only the refit on resume) and reported against the development baselines
  on fold 0 beside their early-stopped scores. The ablation pilots carry the flag from the start.
- Checks on the slice (linear-response encoder, sigmoid head, fold 0, three epochs): a refit interrupted mid-epoch and
  resumed gives the same test predictions as one run straight through (largest difference 0), and so does a refit
  added on resume to a run that had finished without one.
- Cost: best epoch + 1 more epochs per run.

Amendment, 8 October 2026, before any lockbox run: fixes from a code review of the study (the user's request the same
day). Each is a defect against text already in this specification, not a change of the test.
- The TransE baseline is now seeded by --seed (experiments/run_baselines.py); before, every seed fitted the same TransE,
  although experiments/score_confirmatory.py describes the seed as setting its initialisation. Seed 0 is unchanged, so
  the development baselines stand.
- The models' subsystem and module hold-outs run on the development perturbations (--lockbox without --score-lockbox);
  the trainer refused that combination, so the secondary subsystem reading of the models could not run without the
  lockbox perturbations.
- The popularity and degree-scaled popularity baselines take their base rates over the labelled pairs only
  (mechanistic_pathway_learning/models/baselines/popularity_baseline.py). They counted a positive pair the selection set
  aside as a negative, while this specification, the trainer's loss and its leak initialisation treat it as neither.
  Within a symptom (and within a degree stratum) these baselines rank by degree alone, so macro AUPRC and both
  within-strata readings are unchanged; the pooled micro ranking changes, by -0.012 to +0.006 in per-fold micro AUPRC on
  the five development folds (a check run during the review, development rows only).
- read_lockbox also refuses data loaded with another --group-by or another evidence table than the lockbox was drawn on
  (it caught only groups straddling the lockbox, which a finer grouping never does).
- The pathway-wise hold-outs (the secondary subsystem reading, and the curated modules) leave out of each fit the
  perturbations that share a leakage group with a held-out one; they select perturbations by seed gene, so a held-out
  gene's disease-cluster or drug-target partner trained with its labels. On the development set five of the six
  qualifying subsystems have such partners (227 to 343 perturbations, mostly cluster:ABCA3); the lockbox holds whole
  groups and the grouped folds are drawn by group, so neither changes. The --min-holdout-positives filter counts kept
  positive pairs (it counted every positive, so hold-outs qualified on pairs that are never scored): on the development
  set 6 subsystems qualify instead of 21, and no curated module reaches 10 kept positives (the largest has 9).
- A resumed run now trains as an uninterrupted one does (experiments/run_main_model.py). Before: the epoch order was
  reseeded with seed + epoch on resume, where an uninterrupted run draws it from one generator; the torch generator
  behind the noisy-OR gate noise restarted; a checkpoint written in the middle of an epoch restarted that epoch from
  its first batch with a new order, after its first batches had already been applied (a checkpoint on the last batch
  of an epoch trained that epoch twice); and without the refit a crash after early stopping trained further epochs on
  resume. The checkpoint now holds both generators' states and, mid-epoch, the epoch's order, the next batch and the
  partial sums; a run that stopped early is marked finished. Checkpoints are written to a temporary file and renamed.
  On a toy graph (52 perturbations; both heads; resume at an epoch boundary, a mid-epoch interrupt, an interrupt on an
  epoch's last batch, a crash while scoring) every resumed run gives the test predictions of the uninterrupted run
  (largest difference 0), and an uninterrupted run gives the predictions of the code before the change. Runs that
  resumed under the old code keep their results; their test scores are what that run trained to, not what a run from
  scratch would give (the two finished development sigmoid pilots resumed at end-of-epoch checkpoints, their optimizer
  step counts show no batch trained twice).
- A checkpoint now records the training arguments and digests of the loaded inputs (labels, weights, label mask,
  edges, seeds, and the descriptor, cell-class and laboratory files). A resume under other arguments or inputs prints a
  warning and lists the changes in the run's code provenance. A resume does not refuse such a checkpoint, because the
  refit and a longer --max-epochs are legitimate changes on resume and are left out of the comparison. A resume also
  warns when the message-passing node features restored from the checkpoint differ from those built from the current
  files. Before, the split signature covered only the split and the controls. On the toy graph an uninterrupted run
  gives the predictions of the code before the change, and checkpoints written before it resume without a warning.
- With --positive-target-from-frequency (the slice configuration b6_frequency_target only), the early-stopping
  validation loss now uses the frequency targets of the training loss. Before, it scored every positive against 0.99
  under the frequency mode's weights. No confirmatory configuration uses the flag.
- experiments/score_confirmatory.py refuses a run whose predictions are not finite or not one row per lockbox
  perturbation and one column per symptom, whose symptom list differs from the data's, or (rewired variant) whose
  rewiring used another number of swaps per edge than 50; before, a NaN score ranked as the lowest. If a run is missing
  or refused it writes nothing, so a failed run can be resumed without spending the single scoring
  (--allow-incomplete scores with that model not confirmed, only on the user's decision). The scoring record holds the
  git status beside the commit. The scorer also refuses a run whose recorded trainer arguments differ from those its
  configuration in run_main_model_batch.py gives (the resume controls and the variant's own arguments aside), whose
  graph, evidence, descriptor or cell-class file hashes differ from the files being scored, or (permuted variant) whose
  permutation differs from the scorer's draw for that seed. The trainer now records the graph and evidence file hashes
  and a digest of its permutation in results.json. Before, only the baselines' permutation was checked, and nothing
  tied a run to its configuration or to the graph and evidence it read.
- Per-symptom AUPRC is undefined without a positive: per_symptom_auprc returns NaN there, where scikit-learn 1.9
  returns 0, so the per-symptom bootstrap intervals (secondary readings) leave such resamples out instead of scoring
  them 0. The confirmatory scorer computes its own AUPRC, which already left them out.
- Open, for the user's decision: the training loss takes the logarithms of P clamped to [1e-6, 1 - 1e-6]
  (mechanistic_pathway_learning/models/soft_constraint_losses.py), so a pair beyond the clamp has a constant loss and
  no gradient: a positive the sigmoid head scores below logit -13.8 cannot be pulled back, and the validation loss that
  stops training caps it. In the fold-0 development pilots 3.7 percent (message passing) and 1.2 percent (linear
  response) of the sigmoid heads' test probabilities are below 1e-6, among them 13 and 3 labelled positives; on the slice
  0.5 and 0.04 percent. --bce-in-log-space computes both losses from the heads' log P and log(1 - P) (log-sigmoid of the
  logits; for the noisy-OR head the log(1 - P) it already forms); it is off by default, which reproduces every earlier
  run bit for bit (checked on a toy graph for both heads). Whether the confirmatory configurations take it, which would
  need new development pilots of the configurations that change, is the user's decision.
- Decided by the user on 8 October 2026 (amendment below, symmetric rewiring adopted): the rewiring of H2 swaps every stored edge within its relation on its own. binds is
  undirected and stored in both directions (30,682 of its 31,118 edges have their reverse in graph_full_neuronal), and
  after 50 swaps per edge 1,664 do: the null graph keeps every in- and out-degree per relation but not the symmetry,
  so the real-minus-rewired difference also measures the loss of reciprocity. --keep-reciprocated-relations-symmetric
  (trainer and experiments/evaluate_on_rewired_graph.py; negative_controls.reciprocated_relations) rewires a relation
  stored in both directions as undirected edges, by undirected double-edge swaps written back in both directions, the
  treatment Maslov and Sneppen gave undirected networks; on the full graph it keeps all 30,682 reverse pairs, every
  per-relation degree, no duplicate or self-loop, and leaves 0.45 percent of edges unchanged (32 seconds). It is off by
  default, which reproduces every rewired run so far. Not covered by it: a reversible reaction's substrate_of and
  product_of edges between the same metabolite and reaction (15,867 pairs) are two relations and still move apart.
- Decided by the user on 8 October 2026 (amendment below, recomputed from the rewired wiring): a rewired run keeps node features computed from the real gene-reaction associations
  (the reaction_brain descriptor block and the reactions' cell-class weights come from the Human-GEM gene rules), so a
  model trained on the rewired graph still sees which genes catalyse each reaction through its features. This makes
  the rewiring difference of H2 smaller than the wiring's whole contribution: it can cost H2 power, not make it
  confirm wrongly. The alternatives are recomputing those blocks from the rewired catalyzed_by edges, or stating it as a
  limitation.
- The random walk's rewiring control is corrected. The walk merges the relations into one undirected graph with each
  pair once; rewiring the stored directed edges per relation moves apart two edges joining one pair (binds both ways, a
  reversible reaction's two edges), so the walk's rewired graph was denser than its real one (slice: 112,920 adjacency
  entries against about 135,000). --rewiring-method walk_graph (experiments/run_baselines.py and
  experiments/run_rewiring_null_distribution.py; negative_controls.rewire_walk_graph) rewires the walk's own graph so
  every node keeps its number of neighbours. The development reading that the slice wiring carries signal for the walk
  (docs/rewiring_null_distribution.md, 0 of 20 rewirings as good as the real graph, p = 0.048) used the old null. Rerun
  with the new one (docs/rewiring_null_distribution_walk_graph.md), the 20 rewirings give 0.264 ± 0.011 (best 0.285)
  against 0.293 on the real graph, again none as good (p = 0.048), so the reading stands. No confirmatory reading uses
  the walk's rewiring.
- The grade-shuffle and peripheral-event controls the skeleton lists under negative controls (design section 6.3)
  were never implemented and are not run; the confirmatory test has the two controls of its specification (permuted
  labels, a secondary reading, and the rewired graph of H2). The peripheral-event control needs adverse events without
  a plausible central mechanism, which the symptom set (central symptoms only) does not hold; a grade shuffle under the
  better_v1 selection would permute the weights of the kept positives only. configs/evaluation.yaml says so, and now
  states that the preregistration supersedes its endpoint and correction entries.
- Disclosed and left as they are. Changing any of them would change the evidence table or the graph, which confirmatory
  jobs read. Each can move the absolute scores of models and baselines alike. None of the counts below reads a lockbox
  outcome: lockbox-side counts use membership, seeds and disease identifiers only.
  - Disease entries annotated to 20 or more genes do not link those genes into one disease cluster
    (assemble_evidence_table.py; design open question 11). Seven entries are excluded this way (ORPHA:442835 with 49
    genes, ORPHA:154 with 45, and others). Such an entry gives every listed gene the same phenotype annotations, so two
    genes in different leakage groups can share labels through it. 29 of the 253 lockbox genes share such an entry with
    a development gene. In the five development folds, 107 test-fold gene kept positives are supported only by such
    entries, and the same entry gives a training gene the same kept positive. The direction of the effect on the
    model-minus-baseline difference is not known.
  - A drug action of type modulator, substrate or other has sign 0 (perturbation/map_drug_targets_to_graph_nodes.py).
    The linear-response encoder's input is sign times magnitude, so for that encoder gabapentin (CID100003446),
    pregabalin (CID100125889) and brivaracetam (RXCUI:1739745) give a field of exactly zero, the same as no
    perturbation. All three are development drugs: gabapentin has 17 positive pairs (12 kept), pregabalin 16 (none kept)
    and brivaracetam 3 (none kept). No lockbox perturbation has an all-zero input. The message-passing encoder still
    sees them, because its injection layer has a bias. The trainer now prints a warning that lists such perturbations.
  - evidence_full_v2 placed drug seeds on graph_full, not on graph_full_neuronal, which the models read. Every seed is a
    node of both. Thirty drugs, all in development, are seeded on fewer target genes than graph_full_neuronal holds, with
    their magnitudes renormalised over the targets present. For example, 21 GABA-A modulators miss GABRG3, the NMDA
    blockers miss GRIN3B, the gabapentinoids miss CACNA2D2 and CACNG4, and amifampridine misses several potassium
    channels. Nine genes that exist only in graph_full_neuronal (AFG3L2, CLCA4, CLCNKB, COX18, FOXRED1, KCNC2, KCNH5,
    NDUFAF4 and UQCC2) have annotations to crosswalk symptoms but are not perturbations. No gene perturbation sits on a
    node that exists only in graph_full_neuronal, so the leakage grouping is unaffected.
  - Every gene perturbation has sign -1 (loss of function), but grade A accepts Orphanet associations marked gain of
    function. In development, 10 of 1,144 gene perturbations (19 of 1,515 kept gene pairs) have only gain-of-function
    causal associations, so the signed encoder reads the wrong sign for them.
  - The symptom columns come from all perturbations. Catatonia is a column only because one lockbox perturbation has a
    set-aside catatonia pair; development holds no catatonia row. Catatonia is not scored (no kept positive). The
    development micro AUPRC of run_main_model.py, run_baselines.py, the aggregation and the twin comparisons ranks all 22
    columns, while the confirmatory micro AUPRC ranks the 20 scored symptoms.
  - Grade C pairs (human association) are not labels in this specification. They count as negatives in scoring and as
    weak negatives in training (159 such pairs among development perturbations). The positive pairs that the selection
    sets aside are masked instead. Masking grade C pairs in the same way would need an amendment; the user decided it on 8 October 2026 (amendment below, better_v2).
  - A drug pair's label frequency, which the selection compares with 1 percent, is the largest over the pair's preferred
    terms of the mean treatment-arm frequency midpoint of that term's SIDER reports
    (assemble_evidence_table.pharmacological_label_frequency). Post-marketing and rare reports (midpoint 0.0005) lower
    the mean, so 6 development drug pairs are set aside although one of their reports is at least 1 percent. The
    selection entry above does not define the frequency; this is the definition the selection file was built with.
- Decided by the user on 8 October 2026 (second amendment below, lockbox_v2): the lockbox draw over-samples drugs (experiments/draw_lockbox.py draw_groups). Each
  stratum's target is 20 percent of the sizes of all its groups, including cluster:ABCA3 (232 perturbations), which stays
  in development for its size. In the stratum of groups holding a drug with a kept positive, the drawable groups summed
  to less than that target, so all of them went into the lockbox. The lockbox therefore holds 64 of the 142 drugs (45
  percent) and 27 of the 61 drugs with a kept positive (44 percent), not 20 percent. All 34 development drugs with a
  kept positive (215 kept positives) are in cluster:ABCA3. Consequences:
  - The development folds put ABCA3 whole into one test fold, so no development reading trains on drug positives and
    then tests on drug positives.
  - In lockbox runs ABCA3 is in training, and the validation set holds no drug kept positive.
  - The lockbox's drug predictions extrapolate from the drug labels of one leakage group.
  The draw does not leak labels. The options are:
  - keep the lockbox and state this as a limitation;
  - return a seeded subset of the lockbox's drug-stratum groups to development, which only removes lockbox members, so
    the precision estimate is recomputed;
  - redraw the lockbox. This is not recommended, because development perturbations that pilots and the descriptor
    decisions have used would enter it.
  A future draw should leave excluded groups out of the target.
- Decided by the user on 8 October 2026 (amendment below, rotated and stratified): the early-stopping validation set is the same in every seed. assign_grouped_folds sorts
  groups largest first and breaks ties toward fold 0, and validation is fold 0, so the seed only reorders groups of
  equal size. For seeds 0 to 4 it is cluster:ARNT2 (58 perturbations) plus 84 single perturbations, with no drug kept
  positive. ARNT2 holds 10.5 percent of the development kept positives, among them 19 of 51 apathy and 15 of 50
  disinhibition positives. After the refit the scored model trains on ARNT2, but every seed's number of epochs is
  chosen on this one set, so the seed variance leaves out the variance from the choice of validation set. Rotating the
  validation fold with the seed (fold = seed mod 7) leaves seed 0, and so the development pilots, unchanged.
- Decided by the user on 8 October 2026 (amendment below, the group bootstrap decides): the scorer's paired bootstrap resamples lockbox perturbations one by one. Perturbations
  of one leakage group share disease annotations or targets, so their scores co-vary, and a bootstrap of single
  perturbations understates the variance. A check on development rows only drew 10 pseudo-lockboxes by the lockbox rule
  from the 1,222 development perturbations (257 to 299 perturbations each), fitted the baselines on the rest, and
  bootstrapped the three baseline-minus-baseline differences under the scorer's macro and micro statistics (1,000
  resamples). Median 95 percent half-widths:
  - resampling perturbations: 0.030 (macro) and 0.028 (micro);
  - resampling whole groups: 0.037 and 0.047 (median ratios 1.23 and 1.76).
  Seven of the ten pseudo-lockboxes hold a group of 50 to 58 perturbations. The lockbox's largest group has 13, and 210
  of its 237 groups are single perturbations. On the three pseudo-lockboxes whose largest group has at most 16
  perturbations, the half-widths are 0.024 and 0.025 resampling perturbations, against 0.027 and 0.038 resampling groups
  (ratios 1.18 and 1.31).
  The floors (0.041 and 0.028) were estimated with symptoms chosen again in each fold and micro AUPRC over all columns.
  Under the scorer's statistics, resampling single perturbations gives half-widths at or below both floors. Resampling
  groups gives a micro half-width above the micro floor (0.047, or 0.038 on the three lockbox-like pseudo-lockboxes),
  and the scorer's p-values are smaller than a group bootstrap would give. The alternative is a bootstrap over the
  lockbox's 237 leakage groups, either for the confirmatory p-values or as a sensitivity reading beside them.
  Resampling whole groups keeps the dependence inside a group; Field and Welsh (2007) show that "the cluster bootstrap
  gives consistent estimates under both the transformation and the random-effect model" of clustered data.
- Open, for the user's decision: no confirmatory baseline is fitted per perturbation type. Over labelled development
  pairs the base rate is 0.159 for drugs and 0.063 for genes, higher for drugs on 16 of 22 symptoms, and the lockbox
  holds 64 drugs of 317 against 78 of 1,222 in development. A predictor that knows only whether a perturbation is a drug
  therefore earns a share of every reading. experiments/measure_perturbation_type_offset.py (docs/perturbation_type_offset.md)
  draws 200 pseudo-lockboxes of 39 drugs and 154 genes from the development set, fits the predictors on the rest and
  scores the four H1 readings as the scorer does. Against the better of popularity and degree_popularity in each draw
  (mean, 2.5 and 97.5 percentiles):
  - type_popularity (per-symptom base rate of the perturbation's type): macro -0.034 [-0.064, -0.002], micro +0.059
    [+0.018, +0.091], macro within degree strata +0.057 [+0.017, +0.096], micro within degree strata +0.034 [+0.009, +0.064];
  - type_degree_popularity (the same times log(1 + degree)): +0.033 [+0.013, +0.058], +0.087 [+0.044, +0.129], +0.021
    [-0.004, +0.048] and +0.013 [+0.000, +0.027].
  Degree strata do not separate the types: the top stratum holds 58 development drugs and 173 genes. Neither graph-free
  predictor passes H1 on its own, but each reaches the micro floor (0.028) and one reaches the macro floor (0.041) inside
  strata, so a model that learns the type offset beside a weak graph signal can pass H1 against these two baselines. The
  random walk and TransE were not drawn (each needs the full graph per draw), so whether either already carries the
  offset is open. The fold-0 development pilots show the effect in reverse: that test fold holds 63 of the 78 development
  drugs, the 15 training drugs have lower base rates than the genes, and the message-passing sigmoid pilot's macro
  AUROC is 0.343 over all test perturbations against 0.409 and 0.457 inside genes and drugs (type_popularity scores 0.336
  on that fold); the linear-response sigmoid pilot reads 0.671 overall against 0.482 and 0.538 inside the types. The
  options are:
  - add type_popularity and type_degree_popularity to the baselines the best baseline of a reading is chosen from (graph
    free, fitted on the development set like popularity), so that the type offset no longer counts toward H1;
  - score the four readings inside perturbation type as well (degree strata crossed with type), as further H1 readings
    or as secondary readings;
  - keep the baselines and state the offset as a limitation.
  The first is cheap and leaves the models untouched; it changes the decision rule, so it is the user's decision.
  run_baselines.py --with-type-popularity fits the two predictors (off by default, so no run so far and none registered
  changes); the scorer uses them only when its --baselines names them.

Amendment, 8 October 2026, before any lockbox run (the user's decisions late the same day on the open items above). It
defines a second confirmatory family (confirmatory_v2_* in experiments/run_main_model_batch.py). The running development
pilots of the first family stay as they are and are reported as checks that each configuration trains.
- Grade C pairs are masked. experiments/build_label_selection.py --mask-grades C writes better_v2_full_v2: its positive
  rows equal better_v1's, and each pair whose only evidence is grade C (a human association) gets a row with keep False
  and masks_a_negative True, which load_experiment_data reads as neither positive nor negative. It masks 172 such pairs
  of the loaded data, 159 of them among the development perturbations of lockbox_v1; before, they were negatives in
  scoring and weak negatives in training.
- The early-stopping validation set rotates with the seed. run_main_model.py --validation-draw rotated_stratified
  splits the training pool once into round(1 / 0.15) = 7 grouped folds with a fixed partition seed, separately for
  groups holding a drug and the rest (the second stratum's folds offset by one, so the largest group of each stratum
  does not share a fold), and validates on fold seed mod 7, so seeds 0 to 4 stop on five disjoint validation sets. Groups
  larger than a quarter of the expected validation set stay in training (half before). On the lockbox_v1 development
  set as a lockbox run's pool, this gives validation sets of 125 to 128 perturbations whose largest group has 9 to 16
  members, holding 7.1 to 8.8 percent of the development kept positives (123 to 153) and 8 to 11 symptoms with five
  or more; the fold-0 draw gave 142 perturbations with cluster:ARNT2 (58) in every seed and 14.5 to 16.1 percent of the
  positives. The cost: the three largest groups (232, 58 and 50 perturbations, which hold a large share of the positives)
  never validate, so each validation set carries about half the positives of the fold-0 draw and its loss is noisier.
  The epoch count still comes from the validation loss, and the refit trains on every development perturbation.
- The confirmatory p-values and intervals come from a paired bootstrap over the lockbox's leakage groups
  (score_confirmatory.py --bootstrap-unit group: whole groups drawn with replacement). The bootstrap over single
  perturbations is reported beside it as a sensitivity reading. Field and Welsh (2007) give the reason quoted in the
  open item above.
- H2's rewiring keeps the relations stored in both directions symmetric (--keep-reciprocated-relations-symmetric), and
  a rewired run recomputes its reactions' expression from its wiring (--recompute-reaction-expression-on-rewired-graph,
  mechanistic_pathway_learning/graph/rewired_expression_features.py). The user asked why the rewired runs should keep
  descriptors derived from the real graph; they should not. Each reaction's gene rule is rewritten onto its rewired
  catalysts (a kept catalyst maps to itself, the others pair in node order; and/or structure kept), and the
  reaction_brain descriptor block and the reactions' cell-class weights are recomputed from the rewritten rules with
  the functions that built them, from the expression tables of experiments/write_expression_tables.py. With the real
  rules those tables reproduce the stored reaction rows exactly (largest difference 0). Intrinsic node properties
  (protein embeddings, metabolite and reaction chemistry, a gene's own expression) stay, since the rewiring moves edges,
  not nodes. They still carry some of the wiring (a protein's embedding predicts some of its partners), which can only
  make the rewired model better and the rewiring difference smaller, so it can cost H2 power but not confirm it
  wrongly. The scorer refuses a rewired run without either flag. runs/full/confirmatory_lockbox.sh passes both.
- Two tested models instead of four: one head per encoder, chosen on the development set before any lockbox run by
  experiments/choose_heads.py, written before any run it reads existed. For each encoder, the noisy-OR and the sigmoid
  configuration run on the five grouped development folds with seed 0, and so do the development baselines. Per fold,
  the four H1 readings are computed as the scorer computes them on the lockbox, each against the best baseline of that
  reading (the highest five-fold mean, chosen without reference to any model). A head's margin is the smallest of: the
  mean macro difference minus the macro floor, the mean micro difference minus the micro floor (0.032 and 0.043, the
  amendment below) and the two mean within-strata differences.
  The head with the larger margin is tested; margins within 0.001 choose noisy-OR. The other two configurations and
  the four without descriptors stay development readings. With two models tested at one-sided 0.025 each, the chance
  that at least one passes by luck is at most 2 x 0.025 = 0.05 whatever the dependence between the two tests (0.049 if
  they were independent), against at most 0.10 (0.096) with four. The choice reads development folds only, so it does
  not bias the lockbox test of the model it picks.
- Still open: the lockbox's drug share (the user asked for the precision estimate first; a draft that returns lockbox
  drug groups to development until the lockbox holds 20 percent of the drugs and moves unread development gene groups
  in is kept outside the repository and registers nothing), --bce-in-log-space, the per-type popularity baselines and
  the secondary contrasts. The second family's development runs wait for the lockbox decision, because the lockbox
  decides which perturbations are development.

Amendment, 8 October 2026, later the same day, before any lockbox run: the lockbox and the floors (the user's
decisions, after the precision estimate they asked for).
- Lockbox. The confirmatory test scores configs/lockbox_v2.json, derived from lockbox_v1 by experiments/rebalance_lockbox.py
  (its docstring states the rule, written before it was run; seed 20261009). Step 1: lockbox groups holding a drug, in a
  seeded order, return to development while the lockbox keeps at least round(0.2 x 142) = 28 drugs. Step 2: development
  groups with no drug and no perturbation in the fold-0 test or early-stopping validation set of the development pilots
  or in the slice evidence (396 eligible groups, 435 perturbations), in a seeded order, enter until the lockbox holds 317
  perturbations again. Nine groups returned (CHEMBL2023, CHEMBL2094253, CHEMBL2094268, CHEMBL2095181, CHEMBL264 and the
  clusters of ABAT, ALDH5A1, SLC18A2 and SLC6A3) and 38 gene groups entered. lockbox_v2: 317 perturbations (289 genes, 28
  drugs, 11 of them with a kept positive), 266 leakage groups (245 single perturbations, largest 11), 397 kept positive
  pairs, 16 symptoms in the macro average. Development: 1,222 perturbations, 114 drugs, 50 of them with a kept positive
  (before: 78 and 34, all 34 in cluster:ABCA3), 1,814 kept positive pairs. No leakage group straddles the two sides and
  no lockbox perturbation shares a seed node with a development perturbation (checked on the file). The script
  reproduces the file exactly apart from its time stamp.
- What was known when it was derived. Membership, perturbation types and group sizes; the label counts the file records
  (they fix the macro symptoms, as in draw_lockbox.py), which the user saw with the precision estimate before deciding.
  No model output on a perturbation that enters the lockbox entered any decision: the entering genes were outside the
  fold-0 test and validation sets of the development pilots and outside the slice. They were training perturbations of
  the fold-0 development pilots (fitted, not scored) and rows of baseline-only readings: the pre-lockbox precision
  folds (drawn over all 1,539 perturbations, as were lockbox_v1's members), the pseudo-lockboxes of the group bootstrap
  check and those of the perturbation-type check. 85 perturbations of lockbox_v1 that stay are slice perturbations, as
  before. The drug groups that return were never read by a model either. The first family's development pilots trained
  on the entering genes, so they are not readings of the second family; the second family's development runs start
  fresh on the new development set.
- Precision estimate given to the user before the decision. Scaling the precision folds' half-widths (perturbation
  bootstrap) by the square root of the kept positive pairs gives 0.041 macro and 0.028 micro for lockbox_v2 (397) and
  about 0.037 and 0.026 for lockbox_v1 (481). The precision fold whose count is closest to lockbox_v1's (479) gives
  0.028 and 0.022 (3 readings of one fold), the three with 382 to 402 (like lockbox_v2) 0.042 and 0.029 (medians); the
  spread between folds is wider than the counts explain. Under the group bootstrap lockbox_v2, with
  more and smaller groups (266 groups, 245 of them single, against 237 and 210), was expected to come out between 11
  percent narrower and 10 percent wider than lockbox_v1. A lockbox_v1-shaped hold-out cannot be drawn from the
  development data, so the comparison is a projection.
- Floors, revised a third time the same day by the user, before any second-family run: each floor is computed in the
  single scoring from the lockbox itself. For each tested model and for the macro and the micro reading, the floor is
  the half-width of that difference's 95 percent interval under the deciding (group) bootstrap plus 0.005
  (experiments/score_confirmatory.py --floor-rule half_width_plus_margin --floor-margin 0.005, the defaults;
  model_floors). A difference that reaches its floor lies 0.005 above what the lockbox resolves for that comparison.
  For a symmetric bootstrap distribution this is the condition that the interval's lower end exceeds 0.005, so H1's
  macro and micro conditions become a test against a margin of 0.005 instead of zero; with a skewed distribution the
  two differ slightly, and the rule is the half-width form. The floor no longer depends on a projection, and it binds
  by 0.005 above the significance condition whatever width the lockbox gives. Nothing about the floors is computed
  before the scoring, so the blinding is unchanged. The projected floors below (0.032 and 0.043) remain the floors of
  the head choice, where no lockbox interval exists, and the scorer reports them beside the realised floors. Power is
  as stated below when the realised width equals the projected median.
- Projected floors. H1 needed a macro difference of at least 0.032 and a micro difference of at least 0.043 (the user's values):
  the median 95 percent half-width of a difference under the group bootstrap, which now decides, plus 0.005. The
  medians are measured under the scorer's statistics (macro symptoms fixed per hold-out, micro over the scored
  symptoms) on the three development pseudo-lockboxes shaped like the lockbox (largest group at most 16; the group
  bootstrap check above): 0.027 macro (range 0.015 to 0.095) and 0.038 micro (0.025 to 0.053) over nine readings, three
  baseline pairs each (popularity, degree_popularity and the random walk). The first version of this bullet, the same
  day, set 0.053 and 0.042 from a projection, the precision folds' medians (0.041 and 0.028) times the ratio of group to
  perturbation half-widths on those pseudo-lockboxes (1.18 and 1.31). The user asked why the macro floor was not the
  measured median plus 0.005, and it should be: the precision folds' macro bootstrap chooses the scored symptoms again
  in every resample (ranking_and_calibration_metrics.paired_bootstrap_macro_difference skips symptoms with fewer than
  five positives in a resample), so symptoms near the threshold drop in and out and widen the interval, while the
  scorer fixes the macro symptoms in the lockbox file. That earlier projection also stated wrongly that the precision
  folds included TransE; their fifteen readings are the same three pairs over five folds. The precision statement's
  0.041 shares this inflation. With a floor 0.005 above the median width, the floor binds whenever the realised interval
  is no wider than the median, so a difference must also reach a fixed size, not only exclude zero. Three caveats.
  - Nine readings from three pseudo-lockboxes, with a wide range for macro: the median is itself uncertain. If the
    lockbox's realised macro width is larger, the significance test, not the floor, decides H1's macro condition.
    lockbox_v2 averages 16 macro symptoms against 12 or 13 on the pseudo-lockboxes, which should narrow its macro
    interval somewhat. The pairs are baseline against baseline; a model and a baseline may co-vary less, which widens
    the interval.
  - Power. At the median widths (standard errors 0.027 / 1.96 = 0.014 macro and 0.038 / 1.96 = 0.019 micro), a model
    whose true differences equal the floors reaches each floor half the time; reaching a floor with probability 0.8
    needs a true difference near 0.032 + 0.84 x 0.014 = 0.044 macro and 0.043 + 0.84 x 0.019 = 0.059 micro, against
    0.039 and 0.054 for the significance test alone (2.80 standard errors), each on its own; H1 also needs both
    within-strata readings above zero.
  - The scores of the first family's fold-0 development pilots were visible when the floors moved. The floors follow
    from the measured widths by a fixed rule (plus 0.005), not from those scores, and no second-family run had started.
  experiments/score_confirmatory.py now defaults to lockbox_v2, better_v2_full_v2, the two configurations in
  configs/head_choice.json (refusing a choice made on another lockbox's development set), the floors above and baselines
  under runs/full/lockbox_v2_baselines (refusing baselines that scored another lockbox). experiments/choose_heads.py
  takes the same floors. The scorer's statement of the chance that a model passes by luck is now the bound that holds
  whatever the dependence (2 x 0.025 = 0.05 for two models).
- Rewired reaction expression on the full graph (timing check, 50 symmetric swaps per edge, seed 0, on the loaded
  4-core machine): 40 s for the rewiring and 708 s for the recomputation. 7,748 of the 7,759 reactions with a rule
  change rule; 23,410 rule genes are replaced, 347 of them by a catalyst without expression of its own; 8 rule genes
  without a catalyzed_by edge stay. reaction_brain_gtex_brain_max correlates 0.24 between the real and the rewired
  graph and changes on 56 percent of reaction rows; every other descriptor column is unchanged.
- Still open: --bce-in-log-space, the per-type popularity baselines and the secondary contrasts. The second family's
  development runs (its four configurations on the five grouped development folds of lockbox_v2 with seed 0, and the
  development baselines; runs/full/v2_development.sh) start when the first family's pilots end; experiments/choose_heads.py
  reads them. On that development set every fold's early-stopping validation set (92 to 128 perturbations) holds drug
  kept positives except fold 4's. The first family's no-descriptor ablation pilots were retired before they trained
  anything; the second family's ablation runs the two chosen heads without descriptors on the same five folds after the
  head choice.

To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
