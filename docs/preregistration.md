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

Readings of the treatment rule (8 October 2026, descriptor_treatments_message_passing and _linear_response ended; graph
reliance from docs/graph_reliance.md as the mean of macro and micro reliance, the six differences against the plain arm
from docs/descriptor_treatments.md). Not applied: it waits for the user's confirmation (docs/confirmatory_runbook.md,
step 2b).

| encoder | arm | graph reliance | mean of six differences | qualifies |
|---|---|---|---|---|
| message passing | plain (b3_typed_nodes_descriptors_brain) | +0.0045 | | |
| message passing | seed_masked | +0.0149 | +0.0032 | yes |
| message passing | zero_init_slow | +0.0092 | +0.0087 | yes |
| message passing | without_gene | +0.0165 | -0.0003 | yes |
| message passing | without_gene_derived | +0.0150 | -0.0016 | yes |
| message passing | without_protein | +0.0168 | +0.0247 | yes; the rule picks it |
| linear response | plain (b3_linear_response_cofactors_descriptors_brain) | -0.0007 | | |
| linear response | seed_masked | +0.0019 | -0.0095 | yes |
| linear response | zero_init_slow | +0.0072 | +0.0084 | yes |
| linear response | without_gene | +0.0093 | -0.0063 | yes |
| linear response | without_gene_derived | +0.0157 | -0.0047 | yes; the rule picks it |
| linear response | without_protein | -0.0030 | +0.0047 | no (reliance below plain) |

What the readings support. Every arm moves graph reliance by 0.01 to 0.02, against standard deviations over folds of
about 0.02 (docs/graph_reliance.md, one rewiring seed), so the order among arms is not resolved. In the message-passing
encoder without_protein leads without_gene by 0.0003 in reliance, a tie; it is also the only arm with a clear AUPRC
gain over the plain arm (micro +0.043 pooled, 95 percent interval +0.021 to +0.063; macro +0.018, -0.002 to +0.036),
and against b3_typed_nodes without descriptors it gains +0.021 macro pooled (+0.005 to +0.039). In the linear-response
encoder without_gene_derived leads without_gene by 0.006 in reliance and costs -0.005 on the six-difference mean. The
slice holds gene seeds only, so a block's effect on drug perturbations is not measured. If applied, the
message-passing configurations would add --drop-descriptor-blocks protein and the linear-response configurations
--drop-descriptor-blocks protein gene_brain reaction_brain; the latter takes the reaction_brain block out, so for
linear response the rewired reaction expression of H2 would recompute only the reactions' cell-class weights.

The user's decision, 8 October 2026: no arm that drops the protein block ("Which means don't drop proteins", after "the
point of this project is [not] just to predict the symptoms from a pathways graph. It's to be able to use the graph to
determine what the pathways for symptoms are"). Every block arm drops the protein block, so neither rule pick is applied.
The protein block is the 65 columns derived from ESM-2 embeddings that the gene nodes carry (12,643 of 12,810 gene nodes
of graph_full_neuronal; the 1,681 protein_entity nodes carry no descriptor). No arm removes a node or an edge, so under
every arm the modules can name the same genes. Restricted to the arms that keep every block (plain, seed_masked,
zero_init_slow), the rule picks seed_masked for message passing (reliance +0.0149, mean of six differences +0.0032) and
zero_init_slow for linear response (+0.0072, +0.0084). seed_masked hides the descriptors only at each perturbation's
seed nodes, which are gene nodes for gene and drug perturbations alike (evidence_full_v2); every other gene keeps its
protein block. Against the unrestricted message-passing pick, seed_masked keeps most of the reliance (0.0149 against
0.0168, inside one fold standard deviation) and gives up the AUPRC gain (six-difference mean +0.0032 against +0.0247).
Not applied: the restricted picks wait for the user's confirmation; until then the confirmatory_v2 configurations keep
the plain treatment.

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
- Still open: --bce-in-log-space, the per-type popularity baselines, the secondary contrasts and --start-at-weighted-optimum
  (docs/best_epoch_zero.md: on the first family's pilots the noisy-OR leaks stayed at the raw base rate, 0.061, against a
  loss-optimal constant of 0.204, and the modules supplied the difference; the sigmoid head reached the constant within
  the first epoch; no pilot's best validation loss fell below that of a constant per-symptom prediction). The second family's
  development runs (its four configurations on the five grouped development folds of lockbox_v2 with seed 0, and the
  development baselines; runs/full/v2_development.sh) start when the first family's pilots end; experiments/choose_heads.py
  reads them. On that development set every fold's early-stopping validation set (92 to 128 perturbations) holds drug
  kept positives except fold 4's. The first family's no-descriptor ablation pilots were retired before they trained
  anything; the second family's ablation runs the two chosen heads without descriptors on the same five folds after the
  head choice.

Amendment, 8 October 2026, before any second-family run (the user's decisions): the sigmoid head is dropped and the
modules are fixed first.
- Head. "I guess we should drop the sigmoid head then". The tested models are the noisy-OR configurations of the two
  encoders, confirmatory_v2_message_passing_noisy_or and confirmatory_v2_linear_response_noisy_or. No head is chosen on
  the development runs: experiments/choose_heads.py stays as a record and is not run. configs/head_choice.json records
  the two models as the user's decision, and the scorer and the lockbox job read it as before. With two tested models,
  the chance that one passes by luck stays at most 0.05. The second family's sigmoid configurations are not run. Reason
  (the user's statement of purpose that day, "It's to be able to use the graph to determine what the pathways for
  symptoms are"): only the noisy-OR head has modules that can be read as pathways.
- Order. "Unacceptable to start the v2 experiment runs yet until that is done. Ya. Fixing the modules first makes
  sense." The second family's development job (v2_development) was stopped before it trained anything and is paused.
  It runs only after the user decides the module fix. The module fix is developed on the slice under the rule of
  docs/module_health.md, written before its runs.
- The descriptor-treatment readings above come from slice arms with the sigmoid head (b3_*), so they bear on the
  noisy-OR head only indirectly.
- Leak start (the user, the same day: "Ya. Start at the optimum."). The second family's noisy-OR configurations and
  their no-descriptor ablations carry --start-at-weighted-optimum, so each leak starts at its symptom's loss-optimal
  constant and not at the raw base rate (docs/best_epoch_zero.md). The slice runs of docs/module_health.md now say
  whether the leak start is enough or more of the fix is needed. The v2 runs still wait for that reading.
- Not decided: the descriptor treatments (the restricted picks above); where the descriptors
  sit (the user: "The genes shouldn't need descriptors but the proteins should"; in graph_full_neuronal a gene node
  stands for the gene and its product, and 11,723 of the 12,810 gene nodes take part in a protein-level relation, while
  the 1,681 protein_entity nodes are Reactome complexes and sets); the confirmatory pathway test (options and the
  decisions it needs: docs/pathway_test_options.md).
- The user, the same day, on the gene and protein split (docs/gene_protein_split.md): "The rest of that suggested split
  between the genes and proteins makes sense I think except the registered degree strata SHOULD change to the new split",
  and "the brain expression is still provided via pseudo-gene descriptors". Then: "It's just that the descriptors would
  be limited to only the brain expression columns. Also, the several genes [to one] protein change is fine but I would
  also require an equivalent several proteins from one gene. Seed [masking] would then be included." Built into new
  graph directories, not applied to any confirmatory configuration: gene nodes carry the brain expression columns
  only, protein nodes the protein descriptors, one protein node per reviewed entry (37 genes have several), seed
  masking on the slice message-passing arms. Measured: the mapping is not one to one (57 gene nodes share a reviewed
  entry with another gene, 37 have several entries, and isoforms are not resolved by any source), the sources copy
  most edges of a gene onto all its entries (1,082 of 1,395 edge ends), and on the split 140 of 142 drugs fall in the
  top degree stratum. Open: edge ends that name no entry (sent to every entry by default), seed masking over encodes
  partners (on by default on a split graph), the encodes relation in the H2 rewiring, the 110 nodes without a reviewed
  entry, and whether the confirmatory_v2 configurations move to the split graph.

Amendment, 9 October 2026, before any second-family run (the user's decision): the confirmatory family trains on the
graph that splits genes from proteins. This closes the last question of the amendment above.
- The decision, verbatim: "the confirmatory config should definitely be the split graph and include the protein
  descriptors and brain expressions, compartments, everything that physiologically makes sense."
- What changes. The four `confirmatory_v2_*` configurations and their four no-descriptor ablations read
  `data/processed/graph_full_neuronal_split`; the node descriptors come from
  `full_neuronal_split_descriptors_brain_expression.parquet`, which puts the 64 protein components on the 12,709
  protein nodes and the 27 brain region and cell-class expression columns on the gene nodes, with the metabolite and
  reaction blocks and the compartment features as they already were (docs/node_descriptor_columns.md, version 3); the
  linear-response encoder reads `full_neuronal_split_cell_class_weights.parquet`; message passing gets a third layer,
  because a perturbed gene node reaches its protein one hop later (decision 9 of docs/gene_protein_split.md). The
  development baselines (runs/full/v2_development.sh) and the lockbox baselines and scoring
  (runs/full/confirmatory_lockbox.sh, which now passes `--graph-dir`) read the same graph, so the models and the
  baselines see one graph. experiments/score_confirmatory.py itself is unchanged; it takes each run's descriptor and
  cell-class files from that run's own recorded arguments.
- What does not change. The evidence table (`evidence_full_v2`), the better_v2 label selection, configs/lockbox_v2.json
  (drawn on the evidence, not on the graph), configs/head_choice.json, the two tested models, the decision rule, the
  one-sided 0.025 H1-then-H2 sequence and the floors. No run directory existed under any `confirmatory_v2` name when
  this was written, so no finished run changes meaning and nothing is re-scored.
- Degree strata, measured rather than asserted. On a split graph a perturbed gene node has given its product-level
  edges to the protein, so its seed degree is no longer the degree the merged node had, and on the slice it is 1 for
  every perturbation, which collapses the registered strata into one. `perturbation_degrees_for_strata` therefore
  returns, on a graph with protein nodes, the degree of the seeds and of the proteins they encode, less the `encodes`
  edges themselves. Measured on `graph_full_neuronal_split` against the merged graph's seed degree over the 1,539
  perturbations of this evidence table: the same degree value for 1,527, Spearman 0.998, and six perturbations change
  degree quintile (CASP2, FGF8, GBA1, H3-3A, OTX2, WT1); on the slice the value is identical for all 451. The
  stratification is therefore preserved, not redefined, and is not re-registered. The first implementation of this
  hop read every node one edge away; that agreed on the slice but moved 43 percent of the full graph's strata, and it
  was corrected before any run (docs/gene_protein_split_results.md carries both readings).
- Still open, and blocking the lockbox run. Whether `encodes` is rewired in the H2 control (rewiring it would break
  the gene-to-protein identity the split exists to express, so holding it fixed is the proposal, not yet a decision);
  the descriptor treatments; and the development pilots of the two split configurations with the split baselines,
  which must be run and read first.
- The slice reading of docs/gene_protein_split_results.md is not evidence against this change: the user does not
  accept the slice as a verdict on the split until the two confirmatory models are tested on it.

Amendment, 9 October 2026 (second), recording the user's decisions after the first. Each is quoted and then stated as
a rule. No lockbox run has started, and none of these changes a finished run.

- **`encodes` is not rewired in the H2 control** ("Encodes should not be rewired in H2."). The H2 rewiring swaps edges
  within every other relation and leaves the gene-to-protein relation as it is
  (`fast_degree_preserving_rewiring(..., fixed_relations=[encodes])`, which is what `--rewire-encodes` off already
  does, so runs/full/confirmatory_lockbox.sh needs no change). Rewiring it would send gene X's knockout to gene Y's
  protein and to that protein's descriptors, so the H2 reading would measure scrambled identity beside scrambled
  wiring, and a model could lose the H2 comparison for a reason the hypothesis does not name. This closes the first of
  the three items the first amendment of 9 October named as blocking.
- **The protein descriptors keep rank 64** ("I understand the data support more dimensions for the protein descriptors
  but remember that we don't want them to completely outcompete the graph structure. I guess 64 is fine for now."). 64
  is the ceiling of `RANK_GRID` in graph/protein_descriptors.py rather than an interior optimum: held-out R-squared of
  the annotations still rises at the last point (0.113 at rank 48, 0.127 at 64), so the prespecified rule "the smallest
  rank reaching 95 percent of the best on the grid" selects the largest rank on offer. The width is not what keeps the
  graph in charge; the descriptor treatments do (the amendment of 8 October on how the node descriptors enter, and the
  graph-reliance reading of docs/graph_reliance.md), and the no-descriptor ablation of step 1d measures what the
  descriptors are worth on their own. If the treatments leave the descriptors dominant, the answer is a treatment or a
  dropped block, not a narrower table.
- **Repeated affinity measurements of one drug-target pair are aggregated by their maximum** ("Oh. In that case, I
  guess Max is fine."), in place of the median that experiments/scope_off_target_binding.py used. The premise was
  corrected first: a pair with more than one value does not hold one assay run twice but two curators' aggregates of
  different papers (DrugCentral keeps one value per pair 96.5 percent of the time, GtoPdb already medians its own
  papers), and 2,367 pairs are in both sources with a median absolute difference of 0.10 log units, a 90th percentile
  of 1.13 and 12.1 percent differing by more than 1 log unit. The maximum is defensible by Cheng-Prusoff, since an
  IC50-derived pK is biased downward against a Ki and 34.3 percent of DrugCentral pairs carry no Ki or Kd. My own
  proposal was the type-aware rule (median over Ki and Kd when present, maximum otherwise); the user's decision
  overrides it, and the type-aware rule is recorded here as the untried alternative. This changes the off-target
  scoping only: the measurement-coverage strata count genes with an affinity, not their values, so they do not move.
- **Every plasma substrate binder is modelled, and the measured-occupancy extension waits for them** ("All of the
  plasma substrate binders should be modeled. Once AAG and albumin are added, there will be state to attach that
  separate extension to."). The scope, the data audit and what is missing are in docs/plasma_protein_binding.md. The
  graph already holds the binders as gene nodes (ALB, ORM1, ORM2, SERPINA6, SERPINA7, TTR, AFP, SHBG and the
  apolipoproteins) with no edge to anything they carry, so the first stage is binding edges for the substrates, and
  the second is the fraction unbound per drug, which no pinned source in this repository carries. Nothing in the
  confirmatory specification changes on this yet: the extension needs the binder state, and the binder state needs
  data that must be fetched and pinned first.
- **The off-target summary enters the drug input normalisation, not the feature vector** ("You can drop the off target
  summary into the normalization but mark it down as a note in case a situation comes up where it is needed to try.").
  So when off-targets enter a drug's seeds, the summary of their occupancies scales that drug's input alongside
  `drug_input_normalised` (experiments/run_main_model.py), and no drug-level scalar is added beside the graph. The
  note: the standalone feature (mean occupancy over the kept off-targets, which is what I would have taken of the
  three the user offered) stays available and untried, as do log(number of off-targets) and log(sum of occupancies),
  both of which are count-driven and closer to a research-attention proxy than to pharmacology.
- **The leakage-group rules are measured, and the registered grouping stays until the user decides** ("I don't like the
  alphabetical thing. Try the other version of family tie breaking and also make sure to include a measurement of the
  remaining leakage confound from this decision."). experiments/measure_leakage_group_rules.py, run on two evidence
  tables because the answer differs between them: docs/leakage_group_rules.md on the off-target scoping table
  (evidence_full_v3_max3_targets, up to three mechanism targets per drug) and
  docs/leakage_group_rules_evidence_full_v2.md on the table the confirmatory family reads (evidence_full_v2, one
  mechanism target per drug).
  - The tie-break by best-bound mechanism target does what it was asked to do where it applies. On the scoping table it
    takes the largest group from 0.18 to 0.14 of the drugs and moves 15 drugs, among them the antipsychotics, from a
    serotonin-receptor family to the dopamine-receptor family their own affinities name. On `evidence_full_v2` it
    changes nothing at all: no drug there has mechanism targets in more than one family, so the two tie-breaks give the
    same groups and the question is moot for the confirmatory split as it stands.
  - The remaining confound is not the tie-break but the join the family rule drops. On the scoping table, where a drug
    can carry three mechanism targets in different families, 0.49 of test drugs share a mechanism target with a
    training drug under either tie-break, carrying 0.57 of the test fold's kept drug positives, against 0.000 under the
    registered rule, which puts two drugs sharing a target node in one group by construction. On `evidence_full_v2`
    that part is 0.000 under every rule, and what remains is the drug-to-gene leak of docs/drug_target_leakage.md:
    0.377 of test drugs target a gene whose own perturbation is in training, carrying 0.411 of the test fold's kept
    drug positives, against 0.000 under the registered rule.
  - So the family rule's gain is bought with that leak, not with a better partition. Its largest group holds 0.18 of the
    drugs against 0.44 for the registered rule on `evidence_full_v2`; putting the drug-to-gene join back closes the leak
    and takes the largest group to 0.45, which is the registered rule's own figure. `group_by` therefore stays
    `disease_cluster_and_targets` for the confirmatory reading, and the family grouping is available as a development
    sensitivity reading, where the measured leakage is reported beside it rather than assumed away.
- **The Reactome entities get their own annotation vector** ("I like the separate annotation vector much better for the
  complexes"), not the mean of their members' descriptors. graph/complex_descriptors.py and
  experiments/build_complex_descriptors.py: the entity's annotation vector is the union (maximum) of its members'
  Pfam families, EC numbers and GO molecular function terms and the intersection (minimum) of their GO cellular
  component and UniProt location terms, with a DefinedSet or CandidateSet taking the mean over its members because it
  is a disjunction rather than a machine; the vector is standardised with the proteome's statistics and projected onto
  the same annotation-space directions the protein descriptors use. docs/complex_descriptors.md reports the one
  measurement that decides whether the entity block can share the protein block's columns: the correlation, over the
  proteome, between the projection of a protein's recorded annotations and the projection of the annotations predicted
  from its sequence. Wiring the table into a graph's node_descriptors.parquet is a build step, and whether the
  confirmatory graph carries it is the user's call at the lockbox go-ahead.
- **Two symptom mapping conflicts stay open** ("I'll have to look at the antipsychotic bradykinesia and the other one
  later"): `psychomotor_retardation` matches `\bbradykine`, so the relabelling of the DRD2-inhibition claim to
  bradykinesia lands back on the symptom it was taken off; and `compulsive_behavior` matches `stereotyp`, which
  `ANIMAL_READOUT_PATTERN` also matches, so human punding cannot enter as a label without an exception to the
  rodent-readout filter (docs/rodent_readout_human_evidence.md). Neither is changed here.

## Amendment, 9 October 2026 (third)

Three decisions and one question from the user, each quoted.

- **The Reactome entity annotation vectors go into the preregistered descriptor table, and the pilots do not re-run
  yet** ("I'm ok with regenerating the pre-registered input but don't rerun the same pilots until we've completely
  nailed down the design options"). `data/processed/node_descriptors/full_neuronal_split_descriptors_brain_expression.parquet`
  is rebuilt with `--complex-descriptors`: 138 columns against 137, SHA-256 `df38d78c78ffe6d4` against
  `27ed5e89ccb369d9`. The 64 protein components now carry the entity annotation vectors on 1,527 of the 1,681 protein
  entity nodes, and the new column `protein_has_complex_descriptors` marks those rows, because an entity's vector is
  the projection of its recorded annotations while a protein's is the projection of the annotations predicted from its
  sequence, and the two are not the same quantity (docs/complex_descriptors.md). What did not change: the split graph's
  `nodes.parquet`, `edges.parquet` and `gene_to_protein.parquet` reproduce byte-identically, gene rows still carry no
  protein column, protein rows no other column, and the protein nodes' own values are unchanged. No confirmatory or
  development run has been launched on the new table, per the second half of the instruction; the pilots are the
  user's call once the design options are settled.
  - One code change was needed and is worth recording, because it would have failed silently on a later build:
    `gene_protein_split.place_descriptors_on_split` stripped the `protein_` prefix off every column of the block and
    looked the result up in the per-entry protein table, so the new flag would have been sought there as
    `has_complex_descriptors`. It now places only the columns the per-entry table holds and zeroes the rest on the
    protein nodes, which is correct: a protein node is not a Reactome entity.
- **The plasma substrate binders are built, orosomucoid among them** ("Implement the plasma binders please, download
  and pin whatever data you need. orosomucoid should also be implemented"). Stage 1 of
  docs/plasma_protein_binding.md: 27 `binds` edges, sign 0, from the extracellular copy of each cargo metabolite to
  its binder, in `data/processed/graph_full_neuronal_binders` and `data/processed/graph_full_neuronal_split_binders`
  (on the split graph the carriage sits on the protein nodes, since binding is the protein's property). Two sources,
  kept apart in `evidence_source`: the newly pinned UniProt binding-site features
  (`data/raw/uniprot/uniprot_plasma_binders.tsv.gz`, release 2026_03, sha256 `eab6550251dd8056`) and
  docs/curated_plasma_carriage.csv, 22 rows each with its source, PMID, DOI and a quoted sentence. Nothing in the
  confirmatory specification changes: the edges are in variant directories, and folding them into the confirmatory
  graph would need its own amendment and the pilots re-run, which the same message defers.
  - Orosomucoid gets one edge per gene, progesterone, and that is the honest limit rather than an oversight. Its
    UniProt function statement names no specific endogenous cargo ("Binds various ligands in the interior of its
    beta-barrel domain. Also binds synthetic drugs and influences their distribution and availability in the body"),
    and a graph with no drug nodes has nowhere to put the basic drugs that are its real cargo. Orosomucoid therefore
    enters this study through the fraction unbound of stage 2, which still has no pinned source.
- **Target family is recoverable from the graph and from the descriptors, which bounds what any drug grouping can buy**
  (the user's question: "is there any graph structure or attribute that corresponds to family? That might be another
  source of leakage."). Yes, and strongly: experiments/measure_family_in_graph.py, reported in
  docs/family_in_graph.md, scores each channel by leave-one-out nearest neighbour over the graph's gene nodes with a
  GtoPdb family. On the 104 genes that are drug mechanism targets (20 families, chance rate 0.092), a held-out gene's
  family is recovered 0.904 of the time from the 64 protein components, 0.912 from its graph neighbours and 0.973 from
  its Reactome protein-entity memberships; the brain expression block, which carries nothing about a family's
  sequence, reaches 0.173. Over all 2,651 family-bearing gene nodes (498 families, chance 0.004) the figures are 0.790,
  0.576 and 0.821 against 0.043. Same-family gene pairs share a Reactome protein entity 0.732 of the time against
  0.023 for different-family pairs.
  - What follows for the grouping question of docs/leakage_group_rules.md: a family-grouped split changes which drugs
    share a fold, not whether the model can tell which family a held-out target belongs to. The confound is a property
    of the graph and of the descriptor block, so the registered grouping `disease_cluster_and_targets` leaves it
    available too. This is not an argument for or against the family rule; it is a ceiling on what grouping drugs can
    remove, and it belongs in the limitations of any family-grouped sensitivity reading. The measurement is an upper
    bound on the confound, not an estimate of its effect on a score: nearest-neighbour recovery says the information
    is present and easy to reach, not that a trained model uses it. The user has deferred the grouping decision ("I'll
    have to think about the leakage family grouping later").


## Amendment, 9 October 2026 (fourth)

The user's decision, verbatim: "yes, include binder edges in confirmatory graph". The second confirmatory family
trains on the plasma-binder variant of the split graph. No lockbox run and no second-family development run has
started, so no finished run changes meaning and nothing is re-scored.

- **What changes.** The four `confirmatory_v2_*` configurations and their four no-descriptor ablations read
  `data/processed/graph_full_neuronal_split_binders` instead of `data/processed/graph_full_neuronal_split`
  (experiments/run_main_model_batch.py, `FULL_GRAPH_ARGUMENTS_V2`). The development baselines
  (runs/full/v2_development.sh), the lockbox baselines and experiments/score_confirmatory.py are passed the same
  directory (runs/full/confirmatory_lockbox.sh), so the models, the baselines and the scorer read one graph. Because
  git ignores `data/processed` and `runs/`, the variant has to be rebuilt by
  experiments/build_plasma_binder_variant.py in a fresh container and the two job-script edits reapplied; both are
  recorded in docs/confirmatory_runbook.md, step 1b.
- **What the variant is.** The split graph plus 27 `binds` edges, sign 0, from the extracellular copy of a cargo
  metabolite to its binder protein: 19 cargo metabolites and the 11 binder proteins AFP, ALB, APOA1, APOB, ORM1,
  ORM2, RBP4, SERPINA6, SERPINA7, SHBG and TTR (docs/plasma_binder_graph.md). It is a strict superset of the split
  graph: the same node ids, the same `relation_types.json` and the same 279,757 edge rows as its prefix, checked
  through the trainer's loader with nothing trained. The node descriptor table
  (`full_neuronal_split_descriptors_brain_expression.parquet`, 138 columns, SHA-256 `df38d78c78ffe6d4`) and the
  cell-class weights are the registered ones unchanged, because descriptors are keyed by node id and the two
  graph-derived model inputs, the `reaction_brain` block and the reactions' cell-class weights, both follow
  `catalyzed_by`, which carriage leaves alone. In the assembled feature matrix exactly one column moves, the
  log-degree, on exactly the 30 nodes the edges touch.
- **What does not change, measured.** experiments/measure_carriage_graph_effect.py, reported in
  docs/plasma_binder_confirmatory_effect.md: the perturbation ids, the symptom list, the leakage group ids, the
  outcome matrix, the label mask and the weights are identical with and without the carriage edges, so
  configs/lockbox_v2.json, the better_v2 selection, configs/head_choice.json, the two tested models, the decision
  rule, the one-sided 0.025 H1-then-H2 sequence and the floors all stand.
- **The degree strata, and the one point that needs stating.** No perturbation's own seed node is touched, because
  the carriage sits on protein nodes while a gene knockout seeds a gene node. But
  `perturbation_degrees_for_strata` hops from a perturbed gene over `encodes` to its protein, so three binder genes
  that are perturbations change strata degree: ALB 49 to 59, SERPINA6 3 to 5 and TTR 17 to 20. Two move degree
  quintile, ALB from 3 to 4 and SERPINA6 from 0 to 1, and both are development perturbations. No lockbox
  perturbation moves quintile, and the lockbox quintile sizes are identical (46, 69, 64, 74, 64). The two
  within-strata readings of H1 and the secondary within-stratum label permutation are therefore unchanged over the
  lockbox and changed over two development rows.
- **TTR is a lockbox perturbation whose graph neighbourhood this amendment changes after the lockbox was drawn**
  (drawn 2026-10-08; TTR holds 3 selection rows and 1 kept positive, cognitive_impairment). That is stated here
  rather than left implicit. It is defensible on one ground only: the carriage edges come from curated carriage
  literature (docs/curated_plasma_carriage.csv, 22 rows each with a PMID, a DOI and a quoted sentence) and from
  UniProt binding-site features (release 2026_03, sha256 `eab6550251dd8056`), neither of which was consulted for,
  derived from or filtered by any outcome label, and the edge set was fixed before this amendment and is not revised
  after any lockbox score is read. Nothing in the lockbox was opened to choose it. The honest residual is that the
  edge set was chosen by a person reading physiology, not drawn by a rule, and TTR's amyloid neuropathy phenotypes
  are the kind of thing that reading could have been influenced by; a reader who rejects that should read the
  lockbox result for TTR as the one row this amendment could in principle have touched.
- **What is still open.** The second stage of docs/plasma_protein_binding.md, the fraction unbound per drug, is not
  in this amendment and is not in the graph: it has no pinned source, and orosomucoid, whose cargo is basic drugs
  that the graph has no node for, enters through that stage rather than through these edges. The development pilots
  of the two split configurations are not re-run here (the user, the same message: "No new pilots yet"), so the
  carriage edges are registered but not yet read on any run.

Amendment, 12 October 2026, written before any lockbox run, before the held-out set has been redrawn and before any
model has been trained on the evidence table it registers. Each item names the decision it comes from; the decisions
are in docs/research_summary.md under their dates and the measurements behind them in docs/label_source_review.md.
The code for items 7 to 9 is not written at the time of this amendment, so those are registered and have been read on
no run.

**1. The evidence table.** `data/processed/evidence_full_v4` replaces `evidence_full_v2`, carrying four changes.
(a) **No limit on a drug's mechanism targets.** The dominant-target admission test is removed (user, 12 October: "I
really don't like having a target limit"; "I meant drop orphenadrine and the target limit"). (b) **A drug whose
target set joins two otherwise separate leakage blocks is set aside**, found by removing each drug in turn and
measuring the largest block. Applied to this table the procedure names orphenadrine and nothing else. Measured: with
the limit, 222 drugs and 5,732 kept positive pairs in blocks of 114 and 71; without it and without the exclusion, 230
drugs with the largest block taking 194 of them and half the kept positive pairs; without the limit and with
orphenadrine set aside, 229 drugs and 5,799 kept positive pairs in blocks of 116 and 77, the largest block holding 51
percent of drugs and 27 percent of kept positive pairs, which is what the three-target rule held to two figures
(docs/label_source_review.md 26.1). (c) **Sign +1 for the 17 genes whose every causal association is a gain of
function**, replacing the disclosure in the registered text that the signed encoder reads the wrong sign for them.
(d) **A kept drug pair whose frequency never exceeds its placebo arm becomes a labelled negative** rather than being
masked (user, 11 October: "Why would working no better than the placebo not be counted as evidence against a
label?"). The table also carries two earlier unregistered changes, drugs with a target of any type with the salt-form
fix, and parkinsonism as the 24th symptom; this amendment registers them together.

**2. The held-out set.** Its membership is redrawn on the table of item 1 by the registered leakage groups. An added
drug that acts on a node of a held-out perturbation is held out as its own stratum, outside the primary score, and
does not train, validate or enter H1 or H2 (user, 11 October: "I accept put the three mechanism target pairs being
mostly in the held out set"). The registered group rule chains; for added perturbations a one-step rule is used
instead, because chaining would pull 149 further perturbations into the held-out set.

**3. The loss.** `--equal-drug-shares` is on in every confirmatory configuration and its ablations: each drug's
positives and, apart, its negatives are rescaled to the mean counts over the training drugs (user, 12 October: "let's
do the control version and the equal loss per share version"). The negatives of item 1(d) take a weight of **0.5 by
default**, and where the drug's label states the trial arm sizes the weight is `n_effective / (n_effective + 214.39)`
with `n_effective = 1 / (1/n_drug + 1/n_placebo)` and 214.39 the median effective size over the **25 drugs** that
state both arm sizes, so a drug of typical size weighs exactly the default (user, 12 October: "I think -.5 is best
unless you disagree" and "look for placebo arm sizes for those specific studies and then calculate negative weight
based on that then use .5 as default for nothing found"). 41 of the 61 pairs carry a measured weight and 20 the
default, and the weights run from 0.057 to 0.772 (`configs/placebo_negative_weights.csv`). A loss weight below zero
is refused by the trainer, and so is a file that gives one.

**Corrected 13 October 2026 to the constant the file was built with.** This item read "214 the median effective size
over the pairs that have one". Both halves were wrong: the median over the 41 pairs is 196.56, and the file's weights
solve to 214.3877, the median over the 25 drugs. The drug median is the right statistic, because the effective size
is a property of the drug's trial and a drug with many symptoms repeats the same number in several pairs, so a pair
median lets those drugs set the constant. The text is corrected to the file rather than the file to the text; the
difference in any weight is below 0.0005 and nothing has been run on either. `tests/test_placebo_negative_weights.py`
now holds the file to the formula, so the two cannot drift apart again.

**4. Reported beside the score, changing no model.** The input-size baseline, macro AUPRC within strata of a drug's
number of targets, and the permuted-target negative control. All three were named before any run on this table.

**5. The occupancy anchor.** Where an off-target ever enters a model input, its occupancy is read from the measured
unbound concentration at each target's own site and not from the placeholder that fixes the best-bound target at 90
percent (user, 12 October). No registered input reads an off-target, so this binds nothing yet and is written here so
the anchor is not chosen by default later.

**6. The comparative checks.** `configs/comparative_checks.csv`, 214 statements named in advance, read after scoring,
testing no hypothesis of the study and changing no decision rule. Only the two-drug kind is a test of a model; the
carrier, displacement and enzyme kinds go through relations no label varies and are reported without a verdict.

**7. The propagation runs to a fixed point.** Both encoders iterate until the relative change falls below **1e-3**,
in place of the registered fixed counts of three message-passing layers and eight linear-response steps (user, 12
October: "Use the convergence instead of 8 steps. 8 steps is arbitrary, not physiological"; "Propagation window is
settled. It goes to a fixed point to simulate steady state"). The perturbation is held on at every round and the
update is bounded so that it is a contraction, without which the iteration has no fixed point. Measured on the
structural propagation at the contraction bound, 1e-3 is reached at 23 steps, and at gains of 0.7 and 0.5 at 14 and
11 (docs/label_source_review.md 17.2). **The claim this supports is the subthreshold steady state**: a fixed-point
model cannot represent excitability, and the electrical layer's seven voltage-gated feedback loops are therefore read
at their settled value and not as an action potential.

**8. The encoder parameterisation.** One set of relation transforms applied at every round, in place of one set per
layer. Each relation owns a diagonal gain of 32 parameters; the nine relations that 2,100 or more development
training pairs reach — `encodes`, `catalyzed_by`, `substrate_of`, `activates`, `product_of`, `binds`,
`regulates_transcription_of`, `inhibits`, `member_of` — additionally own a rank-2 factorisation of 128 parameters.
1,792 parameters in all against the registered 61,440, 0.60 per kept pair against 20.65. The rank follows the number
of training examples that reach a relation rather than its edge count or its flow magnitude, because the flow is
distorted by in-degree; the counts break cleanly at nine relations with 2,100 to 2,468 pairs and eleven with 451 or
fewer, and they are measured on development perturbations only. **Rank 4 is named now as the only step up that may be
taken** if the development runs show the nine underfitting, so that it cannot be chosen after a result is seen. The
factors carry an explicit norm bound, for the same reason as item 7 (user, 12 October: "set the reduced rank bases
with the example share flow through each edge type and 1 repeated layer regardless of depth", and "just do what you
think makes sense instead"; docs/label_source_review.md 26.2).

**9. Message passing reads the edge signs.** Its relation adjacencies take the value `sign / in_degree` in place of
`1 / in_degree`, as the linear-response encoder already does (user, 12 October: "The edge signs are also settled: the
MP has to adhere to the edge signs the same way the linear encoder does"). Twelve relations have a constant sign that
the relation index already carried; the change affects the eight whose sign varies within the relation, 42,903 edges,
of which 7,321 are negative and were read as positive, 7,087 of them transcriptional repressions among
`regulates_transcription_of`'s 41,638 edges. One consequence is registered with it: with signs in the values the
in-degree mean can cancel, so a node with one activator and one repressor of equal weight receives zero.

**10. The Open Targets statements are their own evidence class**, kept apart from gene loss (user, 11 October: "I
accept all the suggestions under section 8 of the label review", whose row 7 states it). The weight that class
carries is not in this amendment and is named in "what is still open" below.

**11. A reading at the peak of the unbound tissue concentration** is for the confirmation stage only, is not trained
on and is not selected on (user, 11 October: "Peak tissue unbound concentration can be relegated to confirmation
only"). Nothing is built for it.

**12. The Open Targets evidence class carries weight 0.45, with 0.0 as a declared secondary.** The number is grade
B's 0.6 times the measured agreement of the class with the study's own labels: of the 118 human-worded Open Targets
triples, 71 name a target that a study drug acts on in the same direction, and 54 of those 71 (0.7606) carry the symptom
on a kept label of that drug, 9 more on a label the selection masked, and 3 (0.042) are contradicted by a placebo
comparison. 0.6 x 0.7606 = 0.4563, **rounded down to the nearest 0.05**, which is 0.45. The direction of the rounding
is not free and is fixed here: the agreement is measured on the triples whose target already has a study drug while
the statements that train are the ones no study drug label makes, so where the rounding has to go one way it goes
towards the lower weight. (Corrected 13 October 2026: this item read "0.6 x 0.761 = 0.457, rounded to 0.45", and
0.457 rounds to 0.46, not 0.45. The registered weight is unchanged at 0.45; what is corrected is the statement of how
it was reached, and `tests/test_open_targets_class_weight.py` now holds the number to the derivation.) This replaces
the 0.3 of the 12 October draft, which was a judgment with no measurement behind it and is withdrawn. **The limitation is registered with the number**: the testable triples
are the ones whose target already has a study drug, and the statements that would train are by construction the ones no
study drug label makes (docs/label_source_review.md 9.3), so the two subsets are disjoint and whether 0.761 carries over
is not something this measurement can say. For that reason a second arm is registered now, before any run: the same
configuration with the class at weight 0.0 and its statements reported as a held-out stratum only. The class is 1.4
percent of training pairs, so the two arms are expected to agree; registering both makes that checkable instead of
assumed. The user expressed no preference between them on 12 October and the standing instruction of that day was to
lock in what I have no objection to.

**12a. The weight is built; the class itself is not, and what the recomputation shows has to be decided before it
is** (13 October 2026). `loss_weight_for_record` now gives an `open_targets` record the fixed 0.45 whatever its
grade, with 0.0 available as the registered secondary, so the arm exists the moment the class does. The class does
not exist: no loader writes those records, and
`experiments/label_review/open_targets_class_composition.py` recomputes on the live table what item 12 states from
the table of 10 October. Two of the three numbers move, and the second is not a number at all.

| | as registered | recomputed on `evidence_full_v4`, after the crosswalk ruling of 13 October |
| --- | --- | --- |
| the class: no kept drug label confirms it and the target is a graph gene | 49 statements, 33 of them training | **65 statements** |
| share of kept pairs | 1.4 percent | **2.18 percent** |
| of the class, inhibition | 16 of 33 | 30 of 65 |
| of the class, activation | 17 of 33 | **35 of 65** |
| agreement of the class with the study's own labels | 0.7606 (54 of 71) | **0.6418 (43 of 67)** |
| the weight that agreement gives | 0.45 | **0.35** |

- **The agreement does not reproduce, and the weight turns on it.** Item 12 derives 0.45 from an agreement of 0.7606
  (54 of 71). Recomputed on the capped v4 table the agreement is 0.6418 (43 of 67), which gives **0.35**, one step
  of the weight scale below the registered number. Four constructions were tried and none returns 54 of 71: matching
  the direction and reading kept labels gives 44 of 70, reading every label the selection lists gives 57 of 70, and
  ignoring the direction gives 72 of 94 (an agreement of 0.766, but on 94 statements rather than 71). The likeliest
  cause is item 1's three-target limit, which the table still carries: a drug's recorded targets stop at three, so a
  drug acting on a statement's target through a fourth is neither testable nor confirming here, and lifting the
  limit adds drugs to both counts. **The weight is therefore provisional at 0.45 and the number that settles it is
  the one the script computes on the rebuilt table, not the one quoted here.** If the rebuild leaves the agreement
  near 0.64 the registered weight becomes 0.35 by the same derivation, and that is a change to make before the
  lockbox is drawn rather than after.
- **The share is 2.18 percent, not 1.4.** The reasoning that the two arms are expected to agree is weakened rather
  than overturned, and the figure moves again at the rebuild, so what is registered is the rule for recomputing it
  and not a number from a table that no longer exists: the script above, on the table the lockbox is drawn from, run
  before the class is built.
- **More than half the class is activation, and an activation is not a label on a row the study has.** The study
  perturbs a gene by knocking it out, bar the 17 seeded as gain of function. "Activation of X causes Y" is a
  different perturbation of the same node, so 35 of the 65 statements would arrive as **new perturbation rows**
  rather than as labels on existing ones. That changes the perturbation set, so it changes the leakage groups, so it
  changes the held-out draw. Item 12 does not say this, and it is why the class cannot be added once the lockbox is
  drawn: it enters at the same rebuild as item 1 or not at all.
- **Two crosswalk mappings are ruled out, and both symptoms leave the class with them** (the owner, 13 October 2026:
  "weight gain -> increased appetite is unsupported as a 1:1 causal inference as is sexual_dysfunction -> decreased
  libido"). A weight change is one outcome of appetite among others and is not a report of appetite; sexual
  dysfunction is wider than libido, covering arousal and orgasmic disorders the symptom does not claim. The two
  wordings were the whole of what those patterns caught: `increased_appetite` had 13 rows of "increased body weight"
  and 2 of "weight gain", `decreased_libido` 7 of "sexual dysfunction", and no row of either says appetite or libido.
  So 22 of 452 curated rows go and **neither symptom is reachable from this curation**. The rejected wordings are
  listed in `REJECTED_WORDINGS` rather than deleted, so a later crosswalk cannot pick them up without reopening the
  ruling. The rest of the crosswalk stands on the same reading ("decreased memory" to cognitive impairment,
  "hallucinations" to psychosis), and the script prints every remaining statement beside the wording it was mapped
  from, which is the review that has not been done.

**12b. No Open Targets statement lands on a scored row** (the owner, 13 October 2026, of the inhibition
statements: "oh ... they're weak evidence", then "ok. drop them"). The weight of 0.45 prices this class in
training and does nothing where it would matter most: **AUPRC reads a label as a label**, so a statement in the
lockbox would score exactly like a grade-A monogenic observation, with no frequency and no recorded absence behind
it. And its effect would exceed its one point, because the macro average takes the symptoms with five or more kept
positives in the lockbox: one such pair could decide **which** symptoms are scored.
`docs/label_source_review.md` 9.3 recommended this when the class looked five times larger: "In the lockbox they
would be scored against the weakest labels of the study".

The two directions resolve differently, because the study already represents one of them.

- **Inhibition carries nothing.** 30 statements. The study's knockout of that gene is a perturbation of the same
  node in the same direction and it is scored, so putting the statement there puts weak evidence into the primary
  score; giving it a row of its own would duplicate a perturbation the table already holds. The ten pairs it would
  have added are not worth either.
- **Activation becomes its own perturbation row, never scored** (the owner: "what about just adding the receptor
  activations as a new perturbation type?", and "what about limiting that perturbation type to the training data
  since it's rows are low and it doesn't have a full mechanism pathway?"). 35 statements on **23 rows**, typed
  `receptor_activation` and seeded +1 at the gene node. Here nothing is duplicated: no other perturbation in the
  table drives these receptors upward, and the activation direction is the one thing this curation uniquely
  carries. The rows are 1.22 percent of perturbations and carry 1.52 labels each.

**12c. What "never scored" means, which is not "always trained".** `run_main_model.withhold_never_scored_rows`
takes a `receptor_activation` row out of every test set, and out of training whenever its leakage group is held
out. The second half is the point: training the activation row of a gene whose knockout sits in the lockbox shows
the model that node perturbed, which is the leak the leakage groups exist to stop. So such a row is **dropped**,
neither trained nor scored, and trained only when its group is not held out. **The held-out draw over every other
perturbation does not move**, because these rows leave the test set whatever group they fall in.

Two reasons this type is kept out of the score rather than merely down-weighted. It carries **no pharmacology at
all** - no exposure, occupancy, route or plasma binding - so the row is identifiable from its inputs alone, while
**43 percent of its labels are one symptom** (cognitive impairment, against 14.4 percent of all kept positives):
"identify the type, rank cognitive impairment" would be a strong rule with no mechanism in it, and macro-AUPRC is
scored per symptom. The base rate is not the exposure it might appear to be: at 0.063 it is exactly the gene rate,
so the type offset of decision 4.6 gains nothing from it.

**12d. The weight is provisional, and the construction behind it is weaker than the number suggests.** Asked how
0.45 was decided (the owner, 13 October 2026: "how did you decide the .45 though?"), the honest answer is that no
session decided it on evidence that reproduces. Its provenance, in order:

1. **0.3**, on 12 October: half of grade B "because halving looked proportionate". Recorded at the time as "a
   judgment with no measurement behind it", and withdrawn.
2. **0.45**, as grade B's 0.6 times the class's measured agreement with the study's own labels, 54 of 71 testable
   triples, 0.761.
3. **The arithmetic does not give 0.45.** 54/71 is 0.7606 and 0.6 x 0.7606 is 0.4563, which rounds to 0.46. The
   item now states the rounding as down to the nearest 0.05, with the extrapolation as the reason for the
   direction, rather than claiming an arithmetic that does not hold.
4. **The agreement does not reproduce.** On the capped table it is 0.6418, 43 of 67, which would give 0.35. Four
   constructions were tried and none returns 54 of 71 (item 12a).

**And the construction multiplies quantities that are not commensurable.** The grade weights are ordinal by
design - `assign_evidence_grades.py`: "Grades are deliberately coarse; the weights are the initial values from the
design" - so 0.6 is a rank expressed as a number and not a probability. The agreement, 0.76 or 0.64, **is** a
probability. Multiplying one by the other treats the agreement as a discount on grade B, which would mean
something only if 0.6 were itself a reliability. The product has no interpretation, and that is a defect in the
derivation and not in the measurement.

**The route that would settle it is already the registered default.** `configs/evidence_assembly.yaml` sets
`weighting: mode: learned_reliability`, and `evidence_reliability_model.py` exists for exactly this reason: "Fixed
grade weights and count-based weights share a flaw: they are set by hand." The Dawid-Skene fit estimates each
source's sensitivity and specificity from the pattern of agreement across sources, with a Beta prior and an
interval. `open_targets` is not in its source list. Adding it there would estimate from the same overlap the 0.761
was hand-computed on, but properly: with a prior, an interval, and no multiplication by an ordinal. The
extrapolation would remain - the sensitivity is estimated on the overlapping triples and applied to the
disjoint ones that train - but it would be explicit in the model rather than buried in a product.

**Decided: the class enters the reliability fit** (the owner, 13 October 2026: "add it to the Dawid-Skene fit
estimates"). `load_open_targets_statements.reliability_reports` emits a report for a statement **only where
another source already reports that item**, and `configs/evidence_assembly.yaml` lists `open_targets` among the
reliability sources. Three things follow from the "only where the item exists" condition, and all three are the
reason for it:

- **It creates no pair.** The assembler makes one observation per (perturbation, symptom, relation) with at least
  one explicit report, so a report on a new item would conjure a scored pair out of a review statement, which 12b
  rules out. Joining an existing item conjures nothing.
- **It reaches the fit and not the grade**, because the grade reads the pair's primary evidence class only. This
  is the route the literature class already takes, so no new machinery is involved.
- **It is estimated exactly where the hand-computed 0.761 was measured.** Those are the statements an existing
  item can corroborate or contradict. The fit redoes that measurement properly: a sensitivity and a specificity
  with a Beta prior and an interval, rather than an ordinal grade weight multiplied by a probability.

**The limit is unchanged and is now explicit in the model rather than hidden in the product.** The statements with
an overlapping item are by construction not the statements that train, so the sensitivity is estimated on one
subset and applied to a disjoint one. That extrapolation was always in the 0.45; what changes is that the fit
states it and carries an interval for it.

**0.45 remains as the fallback** for `--weighting grade`, the fixed-grades arm, where it is a convention and is
labelled one. Under the registered default weighting the class's weight comes from the fit.

**13. `VM_c` is computed, not fitted, and the claim is the subthreshold steady state.** The neuronal membrane potential
is the linearised parallel-conductance equation over the ion pools: each pool feeds `VM_c` one edge whose weight is that
ion's share of the resting leak conductance times RT/zF, six fixed weights adding no free parameter
(`mechanistic_pathway_learning/graph/membrane_potential.py`). Two neuronal extracellular pools are added, sodium and
chloride, so every ion in the equation has both sides; calcium stays a signalling node and does not enter it. The
constants are Düsterwald et al. 2018, eLife 7:e39575, Table 1, under CC BY 4.0, quoted with the place each was read in
`docs/membrane_potential_constants.md`. The law is the parallel-conductance equation and not Goldman-Hodgkin-Katz
because no measured permeability ratio for a mammalian central neuron could be verified and that source states ohmic
conductances, which substituted into GHK miss its own -72.6 mV by 39 mV. The build registered is
`computed_beside_the_reactions`: the 405 `changes_membrane_potential` edges are kept, because deleting them was measured
to lower the field at `VM_c` and not raise it. Registered with it, as 23.3 required: the seven reactions `VM_c` gates
all feed charge back into it, an excitable membrane is not contractive, and a model run to a fixed point under the
contraction bound represents where the membrane settles under a sustained perturbation and not an action potential. It
must not be read as showing excitability.

**14. A drug outside the nervous-system code is admitted by a three-part rule, read rather than by group name**
(user, 12 October 2026: "Can't exclude all beta blockers … just ones that don't act centrally because beta receptors
are all throughout the brain"; "The revised rule seems good for the beta blockers"; "What about combining the rules?";
"Sounds good. Do it."). All three must hold, and each is a property of the drug established before any outcome is
looked at, which is what keeps the drug set from being selected on the symptoms the study predicts:

  a. **the drug's own label states a central action** in a pharmacology section — mechanism of action, clinical
     pharmacology, pharmacodynamics or description — in a sentence about that drug. The adverse-reaction, warning and
     indication sections are not read, because admitting a drug for the psychiatric events on its label is the
     circularity this rule exists to avoid. `experiments/read_label_central_action.py` finds the candidate sentences
     from openFDA and quotes each with its set identifier; the sentence is then **read** and not counted, because
     class boilerplate ("the mechanism of the antihypertensive effects … has not been established; factors may
     include … vasomotor centers in the brain") appears on every beta blocker's label and says nothing about the drug
     in hand, and a combination-product label can put a drug's only match on a sentence about the other ingredient;
  b. **a mechanism target of the drug is a brain-expressed node** of the confirmatory graph, by the graph's own GTEx
     annotation;
  c. **the drug's measured Kp,uu,brain is at least 0.1.** The threshold is read off the measured distribution of the
     class this came up over rather than chosen: nothing measured lies between nadolol's 0.037 and oxprenolol's 0.20
     (docs/label_source_review.md 29.1).

**A label that states the drug does not reach the brain excludes it, whatever the other two parts say.** That is
positive evidence and it is quoted: nadolol "only minimal amounts of nadolol were detected in the brain", labetalol
"Only negligible amounts of the drug crossed the blood-brain barrier", sotalol "Sotalol crosses the blood brain
barrier poorly", montelukast "minimal distribution across the blood-brain barrier", etoposide "the concentrations are
lower than in extracerebral tumors".

On the evidence read (`configs/label_central_action.csv`, `configs/centrally_acting_drugs.csv`), the rule **admits
metoprolol and propranolol** and holds pindolol (penetrant at 0.50, no label sentence of its own), oxprenolol
(penetrant at 0.20, no US label in openFDA), timolol (a central-action sentence, no measured Kp,uu,brain) and
carvedilol (neither). Nothing is added to a table by this registration: the drug set changes only in the v4 rebuild,
which waits on the held-out draw.

**This rule replaces the ATC-group list for additions.** `configs/admitted_atc_prefixes_draft.csv` was never signed
off and is superseded: its content was a group name plus, for six drugs, a hand-read label sentence, and part (a) is
that same evidence standard applied to every candidate. The three groups whose ATC name states central action
(C02A, M03B, A08AA) and the hormone groups (H, G03) are not admitted by group any more; each of their drugs is subject
to the rule. The list stays in the repository as the record of what was considered.

**What the rule cannot do, registered with it.** It is not an automated gate. Over the 924-drug pool parts (b) and (c)
alone admit metformin, warfarin, omeprazole, etoposide and six NSAIDs while leaving carvedilol out, because almost
every drug has a target expressed somewhere in the brain and because having a measured Kp,uu,brain selects on having
been a reference compound in a brain-penetration study (docs/label_source_review.md 29.2). Part (a) is what turns
those out, and part (a) requires a reading. The screen's recall is also imperfect: metoclopramide, which is in the
superseded list with a hand-read sentence, is missed by the screen at eight labels.

**15. The unbound brain exposure reads more spellings of its own quantity, and its species limit is disclosed.**
`experiments/build_drug_exposure_table.py` matches eight further spellings of the unbound-brain ratio that the corpus
uses — `Kp,uu brain`, `Kp,uu, whole brain`, `Kp,uu,WB`, `Kp,uu,ISF`, `Kpuu,ECF`, `Kpuu,BBB`, the `Kp,uu from AUC…` and
`Kp,uu from VT…` forms, and the cerebral-microdialysis form — with a fluid-or-regional exclusion beside it, so no
cerebrospinal-fluid or regional value can enter (user, 12 October: "Yes keep additional kp,uu values"). 89 drugs become
94, none are lost and 16 values change, five beyond 2-fold; the diff is in the decision log. Registered with it, because
it bears on every occupancy the study computes: **the values are predominantly rodent**, and Fagerholm 2022 reports a
rat-to-monkey R² of 0.74 against a rat-to-human R² of 0.01 (n=13), so the species translation is an unquantified error
on this input and no claim rests on its absolute scale (docs/brain_exposure_data_feasibility.md). The cerebrospinal-fluid
exclusion is kept and now has its own evidence: lumbar CSF runs more than 3-fold high for P-gp substrates and
buprenorphine's own measured CSF range spans 46-fold.

**16. The bridging rule is applied once, and only orphenadrine is set aside** (user, 12 October 2026: "I don't know
about excluding tapentadol. It's got NE effects which are very important. I'm leaning towards leaving only orphenadrine
out."). The rule registered in item 1 had no stopping condition, which was an omission in my drafting: applied again to
the rebuilt table it names tapentadol, which seeds `OPRM1` and `SLC6A2` and so joins the 18-drug opioid block to the
97-drug monoamine block, and a third pass names agomelatine. **One pass is the rule.** The reason given for stopping
there is pharmacological and is recorded as such: tapentadol's noradrenaline reuptake inhibition is part of what the
study is trying to represent, so removing the drug to tidy the split would remove a mechanism rather than an artefact.
The cost is stated with it: the opioid and monoamine blocks stay joined, so the largest block holds 116 of 231 drugs,
50.2 percent of the drugs and 27.0 percent of the kept positive pairs. That is no worse than the table in use, whose
largest block holds 51.4 percent and 27.2 percent, so the rebuild does not degrade the split it inherits
(docs/rebuild_leakage_cost.md). The held-out draw therefore cannot separate opioids from monoamines, and no claim may
rest on its doing so.

**17. The enzyme, carrier and two-drug statements are training labels, in both encoders** (user, 12 October 2026:
"Implement"). 106 combinations and 317 positive pairs from `configs/comparative_checks.csv`, carried as the combined
perturbation type (`mechanistic_pathway_learning/evidence/combined_perturbations.py`). Both confirmatory models train
on them, so the two stay symmetric: the asymmetry I warned of on 12 October rested on the `clears` relation being the
only route these examples inform, and that was wrong. Lowering a clearing enzyme also moves the 61 to 103 reactions it
catalyses, through `catalyzed_by`, a relation with 24,266 edges that every reaction-touching label already trains, so a
combination constrains the net of two routes and the message-passing encoder is not being asked to fit 3,072 fresh
weights from 140 pairs (docs/label_source_review.md 27.1). The 21 negative controls of the check table are **not**
labels (user: "That's fine"): their statement is that nothing is expected, and a labelled negative would assert the
symptom does not occur, which no source states. Under `--equal-drug-shares` a combination shares its drug's loss slot,
so a drug with a statement does not take two drugs' share.

**18. Multiplicity: H1 first, then H2 and P1 by Holm at one-sided 0.025 per model** (user, 12 October 2026, on the
three ways set out in docs/pathway_test_options.md: "Accepted"). The chance of a false claim stays at most 0.025 per
model and at most 0.05 across the two. H2 therefore needs p below 0.0125 unless P1 also passes. P1 is read only behind
the two gates, G1 (stability) and G2 (faithfulness), which are preconditions and not tests at a level: the claim passes
only if G1, G2 and the Holm-adjusted test all pass, and a precondition can only lower the chance of passing, so the
conjunction keeps the level. The option this rejects is P1 on its own at 0.025, which would let the chance that at least
one of the four claims (two models x H1 and P1) passes by luck reach 0.10. Registered before any confirmatory module
exists, which is the condition docs/pathway_test_options.md set for fixing it at all.

**19. The probability clamp stays off, and is checked on the development runs** (user, 12 October 2026: "Leave the
clamp off for now but make sure you check for the issue in the dev runs"). The training loss keeps the logarithms of P
clamped to [1e-6, 1 - 1e-6], so a pair beyond the clamp has a constant loss and no gradient, and `--bce-in-log-space`
is not taken. The evidence that motivated it was measured on the sigmoid head, which is no longer in the confirmatory
family (user: "No sigmoid head will be used"), so for the noisy-OR head the clamp is **untested** rather than known to
be harmless, and that is why the check is registered rather than the option. Every development run reports the number
of its test probabilities below 1e-6 and above 1 - 1e-6 and how many labelled positives and negatives are among them
(`clamped_probabilities` in the run record). If labelled positives sit below the clamp in any development run, the
option is reconsidered before the confirmatory runs and not after.

Amendment, 13 October 2026, written before any lockbox run and before the v4 graph is rebuilt. Items 20 to 23 are
built and their tests pass; item 24 lists what the owner's decisions of this date leave open.

**20. A cell class attenuates a response by its limiting step, not by a product along the path** (the owner, 13
October 2026: "min along the path makes WAY more sense to me. do that so it aligns with the physiology"). In place of
multiplying the messages arriving at a node by that node's weight in the class on every round, which charges a route a
factor at every step, the settled response of class k is multiplied once by
`b_k(t) = max over routes to t of min over nodes n on the route of w_k(n)`,
the weakest weight on the widest route, computed by a max-min recurrence over every route (`class_bottleneck`). The
throughput of a chain at steady state is set by its weakest step, which is the same rule the gene rule already applies
inside a reaction, where an "and" takes the minimum. A per-edge stoichiometric demand divides the bottleneck, so a step
consuming two copies of its source is limited at half its availability (the owner: "ACh needing to bind two molecules
to the receptor in order to activate it, which would cut that expression in half"); the graph does not carry those
coefficients yet, so every demand is 1 until the rebuild writes them. At an extracellular pool node shared by the
classes the bottleneck is the largest over the pooled classes, because the pool is one compartment, and an entry node
the perturbation does not seed carries no route. Three things this changes and one it does not. It conditions the
magnitude of a class's response on that class's expression and no longer the route, because a factor applied on every
round compounds whichever factor it is, so the bottleneck has to be factored out of the propagation. The pool then
mixes states no class has attenuated. And the two rules err in opposite directions: the product understates, charging
a route of expressed nodes at every step, while the bottleneck overstates, scaling the sum over every route by the
widest route's ceiling, since masking the thin routes would need a cut-off and the weights have no natural one. What
it does not change is the fixed point, which is reached before the bottleneck is applied.

**No measurement decides whether this rule is used** (the owner, 13 October 2026, on being told the rule was
registered without one: "that doesn't need measured for confirmation to be in the final models. It should just be
in there to mirror the physiology"). He is right, and an earlier draft of this item said otherwise. The rule is
chosen because the throughput of a chain at steady state is set by its weakest step, which is a fact about
pathways and not a hypothesis this study tests; a comparison against the product rule would say how much the
choice moves the predictions, which is not the same question and is not a reason to hold the rule back.

**What does need checking is whether the implementation is the rule**, which is a correctness question and not a
model comparison, and `experiments/measure_limiting_step_rule.py` is repurposed to it. Two properties of the code
are not properties of a limiting step. The bottleneck `b` is computed on the **unsigned, unweighted** structure,
so a route through an inhibitory edge raises the ceiling exactly as an activating route does - "a bottleneck cares
only whether a route exists" is true of a bottleneck and false of a signed one. And `b` scales the settled state,
which already sums every route, so a contribution that arrived through a barely expressed route is credited at the
widest route's ceiling. Both are stated in the encoder and neither has been exercised. The first is a unit test;
the second is a property of the rule as registered and is disclosed rather than fixed, since masking thin routes
would need a cut-off and the weights have no natural one.

**21. A perturbed node does not read its own descriptors: `seed_masked` in both encoders** (the owner, 13 October
2026: "it HAS to use seed_masked to be physiological and to avoid signal leakage"). This replaces the reading of the
treatment rule of 8 October, which picked `seed_masked` for message passing and `zero_init_slow` for the linear
response; under `zero_init_slow` a perturbed node still reads its own descriptors, through a map that starts at zero,
so it does not meet the decision. The reason is a shortcut and not leakage of held-out information: no descriptor block
is derived from the labels, and what `plain` permits is a readout that learns which protein was perturbed and skips the
propagation, which is what the descriptors-only control measures. On the split graph the mask also covers each seed's
`encodes` partner and, for a drug seeded on its entry node, its mechanism targets and their partners. **The seed's
cell-class weights are not masked, and that is deliberate** (the owner, 13 October: "whether the cell class expresses
that knocked out gene at all determines whether the effect even matters in that cell type"): the weights are a separate
input from the descriptors, they set how far a class carries the response rather than naming the node, and under item
20 they are the class's ceiling. The seed's own brain-expression *descriptor* columns (`gene_brain`, `reaction_brain`)
are masked with the rest. The limit already registered holds: a seed's reactions still carry its EC numbers and the
expression its gene rule copied onto them.

**21a. A perturbed node still reads its own regional expression** (the owner, 13 October 2026: "that brain gene
expression masking not acceptable for perturbations ... whether the cell class expresses that knocked out gene at all
determines whether the effect even matters in that cell type"). Of the 138 descriptor columns, the 16 that carry the
perturbed gene's regional and whole-tissue expression stay readable at its seed under `seed_masked`
(`SEED_READS_THESE_DESCRIPTOR_PREFIXES`: the 13 `gene_brain_region_` columns, the two GTEx maxima and the
has-expression flag); the other 122 are masked as item 21 states. The 11 `gene_brain_class_` columns are **not** among
them, because they are the same eleven cell classes as the cell-class weight table, which is already a propagation
input and the better form of them, so unmasking them would add shortcut capacity and no mechanism. The honest residual:
16 columns of expression are still a narrow fingerprint of a gene, so this reopens a small part of the route item 21
closes. What closes it properly is item 23a, which makes the regional axis a propagation input and leaves these
columns redundant in the way the class columns already are. Until that is built the columns are readable, and masking
by column rather than by count is what makes the distinction expressible at all.

**22. Every channel carries the perturbation with the same sign, and the channels start equal** (the owner, 13
October 2026: "each class should be consistently signed within itself"). Each channel's input weight is `exp` of a free
parameter, so it is positive, and all start at one value; the classes stay free to learn different sizes. Before this
each channel held one free weight drawn from a normal distribution, so the sign and the size of the perturbation in a
channel were set by the draw. A readout absorbs a flipped sign everywhere but the extracellular pool, whose mean across
classes is fixed: it averaged copies of one extracellular quantity carrying opposite signs. Measured on a toy chain
with every class weight at 1, the pool lost a factor of 1.2 to 3.1 to that cancellation over six seeds, read negative
for a positive perturbation at one of them, and the channel sizes differed by up to 26 times.

**22a. Compartment pH is a throughput ceiling, read from the settled proton balance** (the owner, 13 October 2026:
"i don't think we've incorporated compartment pH yet, which is important", and of this design: "perfect solution.
implement it"). Protons are already nodes, one per compartment, with 6,271 reaction edges: 2,048 producing and 1,281
consuming reactions at the cytosolic proton alone. Two things stopped that being used, both measured. The nine
Human-GEM proton nodes carry `is_currency`, and a currency node's outgoing edges are dropped, so each one receives and
its deviation is computed while it sends nothing; and the same species has two pools, Human-GEM's `MAM02039*` and
Reactome's `CHEBI15378*`, with different base ids, so they never merged, the Human-GEM copy holding 50 to 100 times
the reactions in every compartment both describe while only the Reactome copy propagates. So whatever pH signal
existed travelled through one to two percent of the reactions that touch a proton. **Resolution:** the Human-GEM node
governs each compartment it describes and Reactome's governs the vacuole, which only it names; and after the
propagation settles, that node's deviation becomes a multiplier on every node in the compartment,
`w_pH = exp(-(sensitivity * deviation) ** 2)`. On the v4 graph that is 18,928 of 49,574 nodes over ten compartments
(cytosol 8,098, extracellular 3,713, mitochondrion 2,104, ER 1,909, peroxisome 984, lysosome 790, Golgi 670, nucleus
568, vacuole 49, intermembrane 43). A weight does not propagate, which is the point: the cytosolic proton touches
3,329 reactions, and a node that sent across all of them would make every pair of reactions neighbours, which is why
the currency mask exists and why unmasking it would be worse than leaving it.
**The form is the log one, which is what a pH is.** `ΔpH = -log10(1 + sensitivity * deviation)` and
`w_pH = exp(-(ΔpH / width) ** 2)`, as `compartment_ph_weight` implements it. Acidifying by a given fraction and
alkalinising by it move pH by different amounts, so the weight is **not** even in the deviation and a
gain-of-function does not take the same weight as a loss-of-function at a node whose compartment the perturbation
moves. That asymmetry is the chemistry and not a defect.

*(An earlier draft of this item registered the even form `exp(-(sensitivity * deviation) ** 2)` and listed among
its properties a **gain-versus-loss mirror**: that a gain of function produces the exact negative of the field a
loss of function produces. The owner never accepted that property and it should not have been in the registration.
It is withdrawn, 13 October 2026. It appears nowhere else in the registration, so nothing rests on it.)*

**What the weight does and does not do.** It is a multiplier on each node's settled value, applied once after the
fixed point and after the bottleneck. It therefore does **not** attenuate the messages a node passes during the
propagation: a collapse of lysosomal pH cannot change the field at any cytosolic node. The cell-class bottleneck
does carry its ceiling along routes, through the max-min recurrence; this one does not.
`docs/physiological_adherence_review.md` section 5 holds its reach (18,928 of 49,796 nodes take a weight; 5,127
transport reactions carry a multi-compartment label and take none, as do all 25,741 gene, protein and drug nodes)
and the conditioning of the proton node the ceiling is read from.

**23. The registered model is the trainer's default** (the owner, 13 October 2026: "all of the v4 stuff should be ON by
default"). Every registered choice that is built is now the default of `experiments/run_main_model.py`: the grouped
split by disease cluster and targets, typed node features, `seed_masked`, the fixed point at 1e-3, the limiting-step
rule, positive input weights, carrier relations and extracellular coupling for the linear response, the graph's edge
signs for message passing, equal drug shares, the placebo negative weight of 0.5, the weighted-optimum start, the
noisy-OR initialisation and learning rates, and the refit on validation. Three safeguards come with it, because a
default that carries a registered choice is how commit `9df62ee` changed a registered model silently
(`docs/confirmatory_architecture_review.md`, row 8). `experiments/run_main_model_batch.py` states the same model as
explicit flags, `tests/test_registered_defaults.py` holds the two equal so a default cannot drift from the
registration without a test failing, and `VALUES_BEFORE_THE_V4_DEFAULTS` gives every such argument the value it held
before this date, so a finished run is rebuilt as the model it was trained as and not as the current default. Each
choice can still be turned off by name. No data path is a default: the graph, evidence, label selection, descriptors,
cell-class weights and combined perturbations are named per run and belong to the rebuild.

**24. What the decisions of 13 October leave open, registered as open so it is not settled by default later.**
- **23a. The channel axis becomes HPA's cluster types, which is the region-by-class cross already in the source**
  (the owner, 13 October 2026: "i thought the region already had cell type expression levels and that the cross just
  has 0s in some spots, which would mean that particular sparse matrix can be mapped to a reduced one for faster
  processing in the bottleneck pass"). Measured: the 11 cell classes are a pooling of HPA's single-nuclei table, whose
  `Cluster type` column is already region-qualified, and its 34 types are that reduced matrix. About 20 are
  region-specific neuronal types (amygdala excitatory; hippocampal CA1-3, CA4 and dentate gyrus; thalamic excitatory;
  cerebellar inhibitory and Bergmann glia; midbrain-derived inhibitory; mammillary body; medium spiny and eccentric
  medium spiny; seven cortical types) and about 10 are pan-regional glial, vascular and immune types that carry no
  region because they have none. With the Siletti dopaminergic cluster, which HPA has no type for, that is 35 channels
  against the present 12: three times the bottleneck pass rather than the thirteen times a full 12-by-13 cross would
  cost, and the pass is cached per perturbation, so it is paid once per perturbation and not once per epoch. Registered
  before it is built; the cost is to be measured on the rebuilt graph, and the 13 regional descriptor columns of item
  21a become redundant when it lands.
- **Message passing is not at its registered form.** Item 8 of 12 October (one repeated layer, rank-factored, run to a
  fixed point) is not built, so the encoder keeps its three layers with weights per layer and does not reach a fixed
  point. Until it is built the two encoders differ in more than their response law, and item 7 is met by the
  linear-response encoder only.
- **Three registered items have no code.** The per-pair placebo weights of item 3 (`configs/placebo_negative_weights.csv`
  is written but not read by the trainer, which applies the 0.5 default to every such negative), the Open Targets
  weight of 0.45 of item 12, and the Holm procedure of item 18. Each is registered and unbuilt, and a run made before
  they are built does not implement them.
- **A gene knockout and an inhibitor of that gene's protein are nearly the same input.** Both reach the protein with
  sign -1 and magnitude 1, differing only in the entry hop, although a knockout removes the protein and everything it
  does while an inhibitor blocks activity and leaves it in place, at a partial occupancy, in the cells the drug
  reaches. The owner's direction of 13 October is that all three differences are to be made (occupancy from item 5's
  measured chain, exposure through the unbound brain ratio, and amount against activity), and that the third needs no
  node of its own for the inhibitor. None is built, and the first two change the graph's inputs, so they belong to the
  rebuild and to the physiology review it waits on.
- **Partial agonism is a halved magnitude in the agonist's direction, which is not what a partial agonist does.**
  `ACTION_TYPE_SIGN` maps a partial agonist, a positive allosteric modulator, a positive modulator and a stabiliser to
  +0.5 and a negative allosteric modulator to -0.5; on the v4 graph that is 357 entry rows at +0.5 and 19 at -0.5, over
  41 of the 222 drugs, of which 33 have no other kind of target. Two distinct errors sit in that bucket. A partial
  agonist's direction depends on the ambient tone, lowering signalling where the tone is above its intrinsic efficacy
  and raising it below, so a fixed +0.5 states one direction for both: aripiprazole at `DRD2` and `HTR1A`, varenicline
  at `CHRNA4` and `CHRNB2`, butorphanol at the opioid receptors, buspirone at `HTR1A`, with nalmefene, vilazodone and
  memantine alongside. And a positive allosteric modulator has no effect without the endogenous ligand, since it
  multiplies an existing signal rather than supplying one, which covers the benzodiazepines and the Z-drugs. Both would
  be answered by the same change, the drug's effect at a receptor scaling the endogenous ligand's state instead of
  adding to it, and that change is not linear, so it cannot sit inside the linear-response encoder without ending the
  linearity items 20 and 22 rest on. Registered as open.
- **Ligand binding competition is surveyed, computable and unbuilt.** `docs/label_source_review.md` 13.4 measured ten
  sources on 12 October: none labels competition for these drugs, and one route computes it from a shared binding site
  plus the two concentrations and affinities, which item 5's table holds for about half the drugs. The widest piece is
  already on disk, DrugCentral's displaced tracer per drug-gene record, which gives 1,021 of 1,307 same-gene drug pairs
  a shared tracer, with ChEMBL's curated site for 21 drugs where the GABA-A sites must be separated. Three limits are
  registered with it. Competition is saturable and so not linear, so it can enter only as a model input computed before
  the propagation and not as a relation inside it. The drug-against-endogenous half needs the endogenous ligand's
  concentration at the synapse, which no source constrains to within orders of magnitude, while its binding site is
  known for 85 of 147 targets. And there is almost nothing to test it on: of the 80 positive two-drug triples only 2
  have opposing actions at one target, and of the two-drug combination statements 19 distinct pairs share a target,
  7 of them with an affinity for both drugs at every shared target. It would therefore be pharmacology built into an
  input, not a tested claim, which is why it is registered as open rather than added (the owner, 12 October: "I don't
  want to add extra confirmatory model variations because that starts increasing probability of winning by luck").

**25. The held-out set is split on mechanism targets, not on drugs, and two readings are registered beside it**
(the owner, 13 October 2026: "NOT splitting on drug itself but splitting on each target individually (so each drug can
have some labels in the test/validation and some in the held-out but all of the labels for a specific mechanism are
consistently split across test and held out)"; "the memorization vs generalization seems like a defensible middle
ground though and i can also report a separate stratum for the unseen mechanisms as an estimate on generalization ...
do it and register all 3").

**Why the drug split cannot stand.** Over 200 full redraws at a held-out fraction of one seventh, **84.4 percent** of
held-out drugs have every target already present in a training drug, because the drugs genuinely share mechanisms: 58
of the 192 targets are held by six or more drugs. No grouping fixes it. Grouping hard enough to remove the sharing
gives a single block of 90 percent of the drugs, and the rare-target rule mitigates without eliminating, as the owner
said ("it doesn't truly stop leakage ... just mitigates it").

**a. The split: a partition of the targets, and a drug is held out whole when any of its targets is in it** (the
owner, 13 October 2026: "i'm actually thinking that's a better fit"). The drug is never trained on, and at test time
it is presented complete, so **the label is exact**: the drug's own observed symptoms are scored against the drug's
own targets.

**What is principled here is the rule, and what is a convention is the size. The registration says which is which**
(the owner, 13 October 2026: "we need a principled way of doing it instead of arbitrarily choosing a number of targets
for hold out, right?", and then, of an attempt to derive the size from the registered held-out fraction: "that's not
principled; it's just pushing the arbitrariness further down the chain"). He is right twice. A target fraction chosen
to reproduce the size the drug split happened to give is convenience; and deriving it instead from the registered
one-seventh only relocates the choice to one seventh, which is itself a convention.

- **The rule is determined by the hypothesis.** A drug is held out when any of its targets is, because the claim under
  test is about mechanism rather than drug identity, and this is what puts a mechanism the model never trained on into
  every held-out drug. Nothing about it is tunable.
- **The size is a convention**, and no construction inside a one-shot lockbox design can make it otherwise. The
  lockbox is drawn once and scored once; only cross-validation removes the choice, and that is a different design.
  What makes one seventh legitimate is not that it is derived but that it was fixed long before this split rule
  existed and is not tuned here.

**The rule governs the drug side; the gene side keeps the grouping it has.** Measured on the live graph, the rule
holds out 14.5 percent of the 222 drugs and so only 2.0 percent of the 1,636 perturbations, because 1,414 of them are
genes and a gene is drawn into it only when its own node is a held-out target. That is as it should be: the leakage
this rule exists to stop is between drugs sharing a mechanism, and the gene side does not have it. Genes are already
in 943 distinct groups with the largest holding 12 percent of them, against the drugs' 17 groups with the largest
holding 51 percent. So the held-out set is drawn on **each side at its own registered fraction**: the genes by the
registered disease-cluster grouping and the drugs by the rule below, one seventh of each rather than one seventh of
the whole.

The target count follows from that fraction: the targets are taken in the order the split's seed fixes and added until
the held-out share of drugs is as near one seventh as a whole number of targets allows. Measured over 200 draws
(`experiments/label_review/target_split_options.py`), that is **3.2 percent of targets**, 6.1 of 192, giving **14.4
percent of drugs** against 14.3, with a spread of 2.3 points.

**The design that has no convention in it, costed and not taken.** Leave-one-target-out: 192 folds, one per target,
nothing chosen. Median 2 drugs held out per fold, mean 6.4, largest 38 (17.1 percent), and 40.7 percent of a held-out
drug's targets seen in training, which is the same split quality the construction above gives (38.8 percent). Its cost
is 192 folds against the registered 60 runs: 192 x 4 configurations x 5 seeds is 3,840 runs. It is recorded here so
that the convention above is a decision taken against a measured alternative rather than for want of one.

What that buys: a held-out drug carries at least one mechanism training never saw, and **61.2 percent of its
mechanisms are new on average**, against the drug split where 84.4 percent of held-out drugs have **every** target
already in training. The residual, stated: the other 38.8 percent of a held-out drug's targets do appear in training
through other drugs, so the model knows part of its mechanism. That is the price of keeping the label exact.

**The registered alternative, and why it is not the primary.** The same sentence also reads as a split on the (drug,
target) **mechanism edge**, where a drug keeps labels on both sides and only its held-out target's contribution is
tested. That guarantees more - the mechanism predicted from is wholly unseen, 0 percent covered - but it costs the
label: the observed symptom was produced by the whole drug, so a partial mechanism is scored against a whole-drug
label, and the noise is one-sided. Since the study's task is (drug to symptoms) with whole-drug labels, the exact
label is worth more than the stronger novelty guarantee, and the edge split is registered here as the alternative
rather than taken. Both readings, and the two weaker ones, were measured before choosing:

| rule | held-out share | fully covered by training mechanisms | has a mechanism training never saw |
|---|---|---|---|
| drug, as drawn today | 14.3% | 84.4% | 15.6% |
| every target held | 3.1% | 81.1% | 18.9% |
| primary target held | 7.1% | 75.8% | 24.2% |
| any target held | 46.0% | 0.0% | 100% |
| **mechanism edge** | **13.8%** | **0.0%** | **100%** |

(that table is drawn at a held-out target fraction of one seventh, which is why `any target held` reads 46 percent
there; the fraction is a free knob and 3 percent is what the primary rule takes.) Neither of the two weaker readings
is taken: holding a drug out only when **every** one of its targets is held out reaches 3.1 percent of drugs and still
leaves 81.1 percent of them fully covered, and splitting on the **primary** target leaves 75.8 percent covered, since
a drug's other targets go on appearing in training.

**b. The permuted-target negative control carries the evidence on memorisation.** Already registered as item 4 and
read as the primary answer to whether the model uses mechanism or drug identity, because no split can settle it
structurally. Permuting which targets each drug has and retraining must lower the score; a score that survives the
permutation is drug identity and not mechanism, whatever the split.

**b2. A symptom is scored only when both sides of the split support it** (the owner, 13 October 2026: "the target
version doesn't guarantee all symptoms being trained either"). He is right, and it is measurable. Over 100 draws of
the rule in (a), a mean of **0.12 symptoms lose every training positive** (at most 1) and **1.05 fall below five**
(at most 2). The cause is visible in the table: one symptom carries a single positive in the whole of it and two more
carry 14 and 15, so a draw that takes that one positive leaves the symptom with nothing. A symptom in that state is
not merely unlearnable, it is **degenerate**: with its links and leak at their initial values the noisy-OR head
returns one fixed ranking of the perturbations for it, identical for every untrained symptom, driven by total module
activation, which is the degree-like quantity the input-size and degree controls already account for. Item 26c gives
the reason, which is an identity of the head rather than a measurement, and the measurement that confirms it
(`experiments/label_review/untrained_symptom_degeneracy.py`: correlation 1.000000, zero perturbations ranked
differently). Scoring it would be scoring the degree baseline under
a symptom's name.

The registered rule already requires **five kept positives in the lockbox** for a symptom to be scored, and that
stays. On the training side the condition is **not five but one**, and the owner's reason is the right one (13 October
2026: "i'm ok with the 5-kept positives for held out but not for training because including it could still help train
parts of pathways for other symptoms"). Two things follow, and they are separate.

- **Every symptom trains, whatever its count.** A symptom with two positives still sends gradient through the modules
  and the encoder, which are shared, so it helps the symptoms that do have support. Dropping it from training would
  throw that away for nothing.
- **A symptom is excluded from the score only when it has no training positive at all**, which is the degenerate case
  and not a thin one: at zero the links and leak never move from their initial values and the head returns the one
  fixed ranking. At two they move, badly, and a bad estimate is still an estimate, so it is scored and its count is
  reported beside it. Measured over 100 draws, that is a mean of 0.12 symptoms excluded rather than the 1.05 a
  five-positive condition would have excluded.

**The draw is not repeated until the condition holds**: rejecting draws on a property of the split would bias which
targets end up held out, which is exactly what the registration exists to prevent. Exclusion costs a symptom;
rejection would cost the split's meaning.

**c. A pre-registered secondary stratum for unseen mechanisms.** Scored apart on the held-out drugs **none** of whose
targets appears in training: the strictest reading available under the primary split, and **21.2 percent of the
held-out side** under the construction above, against the 15.6 percent the drug split would have given. It is the study's only estimate of extrapolation to a
mechanism never trained on, it is underpowered, and it is fixed here before the draw so it cannot be chosen after a
result is seen.

The held-out set is redrawn under this rule, which the owner's direction of 13 October permits ("doesn't have to
maintain current held out set because of the number of changes since that set was drawn").

**26. The nosology readout is a descriptive secondary analysis, reported apart from every prediction claim** (the
owner, 13 October 2026: "i accept the separated metrics for prediction (the 'any target held' version we mentioned
with the closest prefix before the previously registered 1/7 of targets version; not the full target xval because
that's way too expensive) vs nosology readout"). Item 25 governs the prediction side. This item fixes the other side
before any confirmatory run, so the grouping cannot be chosen after it is seen.

**a. What the readout reads.** For each symptom, its column of the module-to-symptom link probabilities, one column
per relation, fitted on all the evidence: `sigmoid(module_symptom_link_logit)[relation, :, symptom]`, 8 modules by 2
relations. That is the symptom's profile over mechanisms in the model's own basis, and clustering the symptoms by it
is the nosology question stated in the model's terms.

**The registration is on the quantity, not on the parameterisation.** If the head is later refactored so that the
links are a product of a shared module-to-latent matrix and a learned symptom vector (the owner's direction of 13
October: "we would remove module_symptom_link_logit and opt instead for a smaller output matrix that is fully
connected across modules and then learn symptom representations alongside that output matrix"), the readout object is
the same profile, recovered as that product, and nothing in this item reopens. The refactor itself is **not**
registered here; it is a head change, and it is listed below as the owner's.

**b. It is fitted in sample, which is legitimate here and is bounded.** The question is descriptive: which symptoms
this model produces through the same mechanisms. For a description, the estimator that uses all the evidence is the
right one. What it forfeits is stated rather than hedged: the readout carries **no** claim that the grouping predicts
anything held out, **no** claim about a symptom outside the set and **no** claim that it would recur in another
cohort. It enters neither H1 nor H2 nor P1, takes no part in the one-sided 0.025 sequence, and is reported in its own
section labelled descriptive. A reader who takes it for a validated nosology has been misled by the report, not by
the analysis, so the label is part of the registration.

**c. A per-symptom fold predicts nothing, and that is an identity. The symptom folds registered below do not
predict.** The identity first, because it fixes what a symptom fold cannot be asked for. With a symptom held out of
training and its own output parameters left at their initial values, its links stay at one common initial logit, so
its link probability is a single constant `c` across modules and

    P(symptom | perturbation) = 1 - (1 - leak) * prod_k (1 - c * activation_k(perturbation)),

where `leak` is the noisy-OR leak of that symptom, the background probability with no module firing, and not the
resting leak conductance of item 13. The product does not depend on the symptom. Two such symptoms are positive
affine maps of one another whatever their leaks, so they correlate at 1 and rank the perturbations identically; the
leak scales and does not reorder. The shared ranking falls as any module activation rises, so it ranks by total
module activation, the degree-like quantity the input-size and degree controls already account for. Measured at the
registered 8 modules and 24 symptoms over three seeds
(`experiments/label_review/untrained_symptom_degeneracy.py`, `docs/untrained_symptom_degeneracy.json`): smallest
correlation between untrained symptoms 1.000000, **zero** perturbations ranked differently by any two of them,
against 0.47 to 0.80 once one symptom trains.

**This item previously read that a per-symptom fold is unavailable "at any parameterisation", and that was wrong.**
The identity above holds only while the symptom's own output parameters sit at their initial values. The
construction the owner specified does not leave them there, and the correction is his (13 October 2026): "i thought
we talked about only doing the symptoms xval in the way that trains a lower rank of output matrix params on the
whole data THEN uses those output matrix params for the per-symptom-folds xval, which is only valid because it's not
being used to test the ability to generalize to new symptoms", and "it doesn't matter how much it leaks held-out
labels. the point is to check for stability".

**d. Stability is measured by resampling perturbations. No symptom fold is registered** (the owner, 13 October
2026: "i don't know that the joint reading refit is worth it, tbh"). The folds were pursued through this item's
drafting and are recorded here as not taken, with the reasons, because the construction is an obvious one to
propose again.

**d1. Why no symptom fold.** Four findings, in the order they arrived, and each narrowed the fold until nothing
was left worth the runs.

- **A fold cannot predict the held-out symptom**, by the identity of (c). That rules out the first thing a symptom
  fold is usually wanted for.
- **The construction that does work predicts nothing either.** Fit the readout on all the evidence, freeze it,
  refit the encoder and modules on the other 23 symptoms, and read the agreement of the supports (decision 4.9).
  The leak of the held-out symptom's labels into the frozen readout is admissible precisely because nothing is
  scored, which is the owner's point ("it doesn't matter how much it leaks held-out labels. the point is to check
  for stability").
- **It does not test circularity**, which is what the fold was being kept for. A module driven by one symptom gives
  that symptom its own module and makes it an outlier, so it produces a vacuous attribution rather than a false
  grouping. Circularity is (e)'s question and only (e)'s.
- **What it does measure is covered, or thin.** Whether a module rests on evidence one symptom uniquely supplies
  decomposes into a degenerate module (module health and gate G1 catch it, with harder perturbations than removing
  one of 24 label columns), a module whose links sit on one symptom (the concentration below, exact and free) and a
  grouping resting on co-occurrence ((e), and the raw label column correlations). What no combination of those
  strictly implies is the joint reading, since the optimiser's behaviour is not a function of them, and that is a
  thin residue against 24 folds x 2 models x 5 seeds, 240 runs where the whole confirmatory is 60. Leave-one-out
  also perturbs weakly, 1 of 24 symptoms, so a module that holds says little and only a module that moves says
  anything.

**What is reported instead, at no cost beyond runs already registered**: the link concentration
`max_s link_ks / sum_s link_ks` of every module in the full-data fit, beside the same statistic from the
permuted-label runs, which is what concentration looks like when no symptom leans on anything and is the
comparison gate G1 already makes for supports. No threshold is drawn on it: (b) fixes that this item makes no
decision, so the concentration is reported as the continuous quantity it is. Were a module's dependence ever
wanted end to end, the scale for "it moved" is the seed-to-seed support Jaccard the five registered seeds already
produce, on the principle the registered floor rule uses, and the route that avoids refitting altogether is an
influence function, recorded and not taken because (basu2021influence) finds such estimates "often erroneous" for
deeper networks and the quantity here runs through the encoder rather than the head alone.

**Decision 4.9 falls with the folds.** Freezing the readout existed to pin module identity across symptom folds.
With no symptom folds, nothing in item 26 freezes anything, and nothing in it asks for a change to the registered
head.

**d2. The perturbation resampling, which is where symptom profile stability is measured** (the owner, 13 October
2026: "resampling on perturbations instead is fine", and "i guess i don't need profile stability under symptom folds
so much as symptom profile stability in general"). Complementary-pair subsamples of the training perturbations
(shah2013variable), refit with nothing frozen, and compare. The profiles move here because the readout is refitted,
which is the other reason a symptom fold could not have carried this reading: under (d1)'s frozen readout every
fold's profiles are identical by construction.

**The reading is at the level of the partition, not of the profile vectors, and that removes the alignment problem
rather than solving it.** The items being clustered are the 24 symptoms, which are the same in every subsample, so
two subsamples' clusterings are compared by one-to-one Jaccard over symptom membership (kuhn1955hungarian) with the
stability index and its interval of (nogueira2018stability), and no coordinate has to mean the same thing in both.
This matters because modules are exchangeable under a refit: a profile vector from one subsample is not comparable
to another's until the modules are matched, and the matching is imperfect whenever modules split or merge. A
partition reading never asks the question.

**The profile vectors are read only behind the module-stability gate.** Where the vectors themselves are wanted -
to say which mechanisms a cluster is built on rather than which symptoms are in it - the modules are first matched
one to one across the subsamples by support Jaccard, and the matched Jaccard is reported beside the reading. A
profile that moves when the module basis has moved says nothing about the symptom, so a low matched Jaccard makes
the vector reading uninterpretable and it is reported as such rather than as instability of the profile. The
partition reading of the paragraph above does not depend on this gate.

The permuted-label runs are the null for both readings, as they are for the supports.

**Neither fold is evidence that the nosology is right.** Both say the structure is reproducible under a perturbation
of the data. A stable grouping of the wrong thing is stable, which is what (e) is for.

**The procedure is fixed now.** Distance: one minus the Pearson correlation between symptoms over the concatenated
induces and relieves profile. Algorithm: hierarchical clustering with average linkage. Number of clusters: taken by
the stability criterion of (benhur2002stability), which also returns the absence of structure ("it provides a means
of rationally defining an optimum number of clusters, and can also detect the lack of structure in data", abstract),
so the readout is able to report that there is no grouping to report. The owner may change any of these before the
first confirmatory run. None changes after it.

**e. The control that decides whether the readout says anything: label co-occurrence.** A symptom's link column is
fitted only to reproduce which perturbations carry that symptom, so two symptoms reported for the same drugs get
similar columns whether or not they share a mechanism. Without a control the readout would recover the label
co-occurrence structure and present it as mechanism. So the same procedure of (d) runs on the symptom columns of the
raw outcome matrix, same distance, same linkage, same stability selection, and the two partitions are compared
cluster by cluster by the same one-to-one Jaccard. **No cluster that the label co-occurrence baseline also returns is
reported as a mechanism-derived grouping.** No threshold is chosen for that comparison: the scale is the mechanism
clustering's own resampling spread from (d), so an agreement inside that spread is read as agreement. What the graph
adds, if anything, is the difference between the two partitions, and that difference is the readout's only positive
claim.

**f. Thin symptoms are reported with their counts and cannot form a cluster by being untrained.** By (c) a column
that has barely moved from its initial value is close to every other such column, so symptoms with little support
would group together and the group would read as a nosological finding while being an artefact of count. The
clustering is therefore over the symptoms meeting the registered condition of five kept positives; symptoms below it
are placed afterwards at their nearest cluster and reported with their counts and no claim. The inclusive clustering
is a declared secondary, reported beside it.

**g. What this item does not settle, and stays the owner's.**
- **The head refactor of (a) is not needed by this item, and nothing else in it is either.** The refactor was
  proposed to make symptom folds possible, then looked load-bearing again as the only way to make profile stability
  measurable under those folds; the owner wants profile stability in general, which (d2) gives without freezing
  anything, and then the folds went altogether. So item 26 asks for no change to the registered head: no latent
  layer, no frozen readout, no new parameter and no free rank. What is left arguing for a refactor is shrinkage
  across symptoms for the thin ones, against an imposed rank, a latent basis in place of named modules, a free
  rank `d`, and another confirmatory variation unless it replaces the current head outright ("I don't want to add
  extra confirmatory model variations because that starts increasing probability of winning by luck"). Decisions
  4.7 and 4.8 keep the two latent forms and what separates them, set aside rather than refused.
- **Whether the readout is published with the confirmatory result or apart from it.** Registering it here fixes the
  procedure either way.

**27. A perturbation can drive a reversible reaction backwards: `product_opposes`, product -> reaction at -1** (the
owner, 13 October 2026: "The reaction should only be able to follow the reversal path under certain conditions",
and "I don't want to accept a perturbation not being able to reverse the reaction").

**Linearised mass action gives reversal without a nonlinearity.** The net rate of a reversible reaction is
`v = k_f [S] - k_r [P]`, so around the steady state `dv = k_f dS - k_r dP`: the product enters the rate at a
negative sign. A perturbation that accumulates the product therefore makes the rate deviation negative and the
reaction runs backwards against its baseline, regenerating the substrate and pulling the product back down; a
perturbation that raises the substrate drives it forwards as before. That is the conditional reversal the owner
described - it happens only when the local state calls for it - and it is linear, so it needs one derived edge and
no gate on the state. Linearity was never the obstacle; the obstacle was that I had framed reversal as switching
edges on and off, when in a deviation model it is the sign of the rate.

*(An earlier draft of this item argued that a nonlinearity was available anyway, "the softmax-mixture aggregator
runs inside the loop". Two things wrong with that, corrected 13 October 2026 when the owner asked why the
aggregator used softmax rather than a sparse maximum. The mixture is **not in the registered model** -
`cross_relation_aggregator` defaults to `mean`, so the registered linear-response loop is linear throughout and
its only nonlinearities are the cell-class bottleneck and the pH weight, both applied after it. And the mixture's
weighting is **sparsemax** by default and has been, returning exact zeros so that a node type either uses a
statistic or does not; only the command-line spelling `softmax_mixture` is historical. Neither changes this item's
conclusion, which needs no nonlinearity at all.)*

**What is built.** The linear-response encoder derives, for every `product_of` edge of a reaction flagged
reversible, the edge product -> reaction at sign -1 under its own relation `product_opposes`, whose gain is learned
and stands for k_r as the other gains stand for elasticities. An irreversible reaction gets none, because a
thermodynamically irreversible step is insensitive to its product. Carrier products, which have no depletion edge,
get none either, for the same consistency. `--reversible-reactions`, on for the linear-response encoder.
`tests/test_reversible_reactions.py` runs the encoder's own adjacency to a fixed point on a three-node chain and
reads the signs: a product rise gives a negative rate and a rising substrate, an irreversible reaction ignores it,
and the forward direction is unchanged.

**What it replaces.** `build_physiology_graph` wrote the mirror pair for every reversible reaction, substrate_of
and product_of in both directions at +1, so the reaction received +1 from its product instead of -1 and every
participant received +1 and -1 from the same reaction. 15,851 such pairs over the 5,713 reversible Human-GEM
reactions, and 1,809 metabolites whose field was identically zero under every perturbation. The build now writes
the nominal direction only, from the stoichiometric sign it already read, and the encoder **refuses** a graph that
still carries the mirror, because on that graph the nominal direction is lost and the derived edge would sit on top
of the cancellation. It therefore reaches a run only at the rebuild; `--no-reversible-reactions` runs the old graph
as it was.

**What it does not do.** It is linear, so it does not saturate: a reaction cannot run backwards faster than its
reverse capacity in the physiology, and here nothing bounds it. That is the limit of the linear response
everywhere and not particular to reversal. The resting direction is Human-GEM's nominal one; a thermodynamic or
flux source could set it per reaction before the propagation (decision 4.10). The message-passing encoder does not
yet carry the relation. And 16 pairs on Reactome reactions not flagged reversible still cancel: species listed on
both sides because one step consumes and regenerates them - GRIP1 and GRIP2, PICK1, the AP2 complex, CoQ and CoQH2,
neuronal K+ in channel gating. Their net stoichiometry is zero, and the right encoding is a modifier edge that
enables the reaction without being consumed or produced. Not built.

**What this amendment does not change.** The two tested models and the noisy-OR head of each, the decision rule, the
one-sided 0.025 sequence of H1 then H2, the floors, the degree strata, the node descriptor table, and the training
stopping rule, which stays validation loss with patience 8 and at most 60 epochs. The head choice is not reopened: the
sigmoid head was dropped from the confirmatory family on 8 October because only the noisy-OR head has modules to read
as pathways.

**What is still open, and is deliberately not registered here.**
- **Off-targets in a model input.** The user asked for them on 12 October ("Occupancy input should be included as
  well as off targets") and the measurement argues against: admitting off-targets at any threshold from 0.1 to 0.5,
  under any concentration scheme, takes the largest leakage group from 241 perturbations and 71 drugs to between 203
  and 209 of the 222 drugs, which ends the drug split that item 2 depends on (docs/label_source_review.md 14.3). The
  objection is recorded and the decision is the user's.
- **Which held drugs of item 14 are admitted once their missing part is supplied.** The rule itself is registered;
  pindolol, oxprenolol, timolol and carvedilol each lack one part of it, and whether a label sentence is read for
  pindolol and timolol, or a Kp,uu,brain is sought for timolol and carvedilol, is the user's to direct. The superseded
  `configs/admitted_atc_prefixes_draft.csv` records what was considered under the old group rule.
- **The reach of the electrical layer, which is an aggregation question and not a wiring one.** Item 13 below
  registers the computed potential. What it does not fix is the magnitude of what arrives at `VM_c`, and the
  measurement of 12 October withdrew the diagnosis this item used to carry: the binding divisor is the in-degree of the
  ion pool (364 cytosolic sodium, 339 extracellular) and not the 405 at `VM_c`, so deleting the reaction edges lowers
  the field there by a factor of 26 rather than raising it (docs/label_source_review.md 27.5). Both alternatives the
  encoder already carries were measured worse in magnitude, so this stays open and is disclosed rather than fixed.
- **The held-out set itself**, which is drawn only on the user's go-ahead ("Don't finalize the lockbox until I've
  given the go ahead").


To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
