# Research summary: aims, results and open questions

Working notes on the project's aims and results. Every number below names the generated document it comes from and
the date it was read; when that document is regenerated, update the number here or mark it stale. Numbers are
per-fold macro AUPRC under the disease-cluster grouped split unless a row says otherwise.

## 1. The question in one paragraph

Inborn errors of metabolism and psychoactive drugs both perturb known molecular pathways and both produce psychiatric
symptoms. The project asks whether a model that routes a perturbation through a mechanistic graph (Human-GEM
metabolism joined to receptor signaling and transcriptional regulation) and reads it through a small set of learned
pathway modules with a noisy-OR output predicts which symptoms a held-out gene or drug produces better than models
that use network proximity or node degree alone, and whether the learned modules match known mechanisms. The noisy-OR
form encodes equifinality: any one sufficient module can produce a symptom. Design: docs/experiment_design.md;
pre-registration skeleton: docs/preregistration.md.

## 2. Aims

### Aim 1. A leakage-controlled benchmark for predicting psychiatric symptoms from gene and drug perturbations

Status: mostly done. The full-data rewiring control has to be repeated at the mixing strength the slice showed is
needed (findings 2 and 8).

- Labels: monogenic gene-symptom pairs from HPO with provenance-based grades and a reliability model; drug-symptom
  pairs from SIDER 4.1 and OnSIDES labels for drugs with a single mechanism target and an ATC nervous-system code.
  Unobserved pairs are unlabelled, not negative.
- Splits: disease-cluster grouping (genes of one disease leave training together), Human-GEM subsystem hold-out,
  time split by evidence date.
- Controls run with every split: labels permuted within degree strata, degree-preserving rewiring, popularity and
  degree-scaled popularity. On the slice the rewiring control is a null distribution of 20 draws at 50 swaps per edge.
- The disease-cluster grouping blocks leakage inside a disease but not across diseases that share a pathway or a
  complex. Measured on the slice (finding 12): the split divides almost every pathway between training and test, and a
  held-out gene's symptoms resemble those of its same-pathway training genes more than its baseline, while complex
  partners carry no such excess. Whether that is leakage depends on the claim: a new gene in a pathway with known
  members is the realistic use, but a claim about learning mechanism needs the subsystem hold-out (finding 3) or a
  pathway-grouped split. Reported separately (finding 13), the trained models trail the random walk only on genes
  whose pathway is split, and tie it on the rest.
- Contribution if Aim 2 fails: a measured account of how much apparent signal in this task is leakage and node degree.

### Aim 2. Test whether the pathway-module model beats the best baseline

Status: slice done (negative); full-data test not run. The full-data runs need more compute than the 4-core cloud
container; cluster access is not confirmed.

- Hypothesis: the proposed model (B6) beats popularity and the random walk with restart under the disease-cluster
  split.
- Pre-registered failure criterion (docs/preregistration.md): if B6 does not beat both, the architecture is
  abandoned before the ablations run.
- Primary endpoint: macro AUPRC difference of at least 0.05 with a 95 percent paired-bootstrap interval excluding
  zero.
- Slice result: no trained configuration beats the random walk, and the criterion is met (finding 5).
- Module diagnosis: done (finding 6). The modules were unused because training stopped while the gates were still
  shrinking. With their own learning rate they carry links, but no module passes the sufficiency test yet.
- Open question: the graph carries signal the random walk reads (finding 8), so the question is why a signed,
  relation-typed propagation extracts no more of it than unsigned diffusion. Finding 11 is a candidate cause.

### Aim 3 (conditional on Aim 2)

- If Aim 2 succeeds: validation on the pharmacological time split (SIDER 2015 pairs train, the 100 later OnSIDES
  candidates test; docs/onsides_label_slice.md) and interpretation of the modules against curated pathways.
- If Aim 2 fails: test whether flux-based perturbation features (design section 5.2, route 2) carry signal the graph
  lacks, first as the flux-feature logistic regression baseline (B4), then a revised architecture.

## 3. Results so far

| # | Finding | Numbers | Source (date read) |
|---|---|---|---|
| 1 | Grouping by gene leaks labels through genes that share a disease, much more on the slice than on the full data | slice random walk: 0.347 by gene, 0.293 by disease cluster, 0.284 by disease cluster on permuted labels; full data: 0.390 by gene, 0.374 by disease cluster (pooled 0.346 and 0.347) | docs/experiment_design.md section 12; docs/phase2_baselines_disease_cluster.md; docs/phase2_baselines_full_gene.md; docs/phase2_baselines_full_disease_cluster.md (6 October) |
| 2 | The full-data rewiring control was too weak to say whether the random walk reads topology or only degree | random walk per fold 0.374 on the real graph and 0.348 on a graph rewired at 2 swaps per edge (pooled 0.347 on both); degree-scaled popularity 0.339; 0.275 on permuted labels. On the slice 2 swaps per edge left the graph under-mixed (finding 8), so the earlier reading of this row, that the signal is degree and not topology, is withdrawn until the control is rerun at 50 swaps per edge | docs/phase2_baselines_full_disease_cluster.md, section "Negative control: degree-preserving rewiring of the graph (2 swaps per edge)"; docs/graph_content_null_results.md (7 October) |
| 3 | Inside a held-out pathway the baselines rank at chance | full-data subsystem hold-out, stratified AUROC: random walk 0.503 (0.560 on permuted labels), degree-scaled popularity 0.557; slice 0.506 | docs/phase2_baselines_full_disease_cluster.md; design section 12 (6 October) |
| 4 | No baseline predicts later monogenic findings | time split at 2015-12-31, 565 new positives: random walk macro AUPRC 0.083 against 0.079 permuted; publication-dated new positives only (129): random walk AUROC 0.432 | docs/phase2_baselines_full_disease_cluster.md (6 October) |
| 5 | On the slice no trained configuration beats the random walk, most fall below it, and the best beat popularity | 18 trained configurations score 0.241 to 0.293; B3 typed nodes 0.293 ± 0.031, B6 mechanistic with a separate gate learning rate 0.292 ± 0.024, random walk 0.293 ± 0.014, popularity 0.231. Paired bootstrap of that B6 configuration, pooled: minus random walk -0.020 [-0.043, +0.006], minus popularity +0.035 [+0.022, +0.064]; within degree strata, minus random walk +0.007 [-0.014, +0.032]. On the pooled paired difference, the pre-registered endpoint, 14 of the 18 are reliably below the random walk (interval entirely under zero); B3 typed nodes ties it per fold but is below it pooled, -0.031 [-0.055, -0.011]. The failure criterion is met on the slice | docs/phase3_main_model.md, generated 7 October 12:44 (7 October). The laboratory-loss and manganese runs finished later and appear only as paired differences in finding 9 |
| 6 | The modules were unused because training stopped while the gates were still shrinking; with their own learning rate they carry links, but none passes the sufficiency test | shared learning rate: in all 30 five-fold runs of six configurations the median gate had moved 92 to 100 percent of the distance Adam allows by the best epoch, expected support 5,000 to 13,500 nodes per module, no link above 0.5. Separate gate learning rate: expected support 187 to 334 nodes (B6 mechanistic) and 61 to 173 (B6 linear response); a link above 0.5 in 13 of 40 and 40 of 40 module-folds. Sufficiency test (ablating the dominant module): the own-ablation interval excludes zero in 0 of 62 symptom-module rows and 2 of 44, one in each direction. The sparse supports miss many genes: with the separate gate learning rate, B6 linear response gives 55 to 63 percent of each fold's held-out genes a prediction within 1e-6 of one shared row (13 to 43 percent exactly equal), so it ranks them by a tie or float noise, against 1 percent for the same run without the gate learning rate; B6 mechanistic 1 to 30 percent | docs/b6_module_diagnosis.md; docs/experiment_design.md section 5.5; docs/b6_mechanistic_gate_time_scales_disease_cluster_modules.md; docs/b6_linear_response_gate_time_scales_cofactors_disease_cluster_modules.md; docs/b6_module_diagnosis.md, columns held_out_at_modal_prediction and held_out_within_tolerance_of_modal (7 October) |
| 7 | The flux route reaches a minority of perturbations | a knockout blocks at least one reaction for 239 of 451 slice genes (440 of 861 positive pairs) and for 225 of 1,282 labelled full-data genes and 15 of 142 drugs (523 of 3,272 positive pairs); every boundary reaction of Human-GEM is open for uptake | measured 6 October 2026 from Human-GEM gene-reaction rules; not yet in a generated document |
| 8 | On the slice the random walk reads the graph's wiring, not only its degree sequence | real graph 0.293 against 0.268 ± 0.010 over 20 rewirings at 50 swaps per edge; none of the 20 as good, p = 0.048, which is the smallest value 20 draws allow; the real score is 2.4 null standard deviations above the null mean. Limits: the margin over the best rewiring is 0.006, below the real graph's own fold spread of ±0.014, and rewiring keeps degree per relation but not the metabolite-reaction bipartite structure | docs/rewiring_null_distribution.md; docs/graph_content_null_results.md (7 October) |
| 9 | Adding content to the graph has not moved the score; the protein descriptors are the only addition with a sign of an effect, and it does not survive correction | every one-argument twin with five folds, 20 pairs, three readings each (per fold with a t interval, pooled with a paired bootstrap over perturbations, within degree strata), uncorrected 95 percent intervals: manganese layer -0.002 [-0.008, +0.005] per fold; neuronal and oxidative variant +0.012 [-0.028, +0.053] per fold, +0.004 [-0.009, +0.015] pooled; its curated layers +0.007 [-0.029, +0.043]; signed logarithm +0.002 [-0.012, +0.015] (B3) and -0.001 [-0.027, +0.024] (B6); laboratory loss -0.002 [-0.008, +0.005]; protein descriptors -0.030 [-0.056, -0.004] per fold (p = 0.034) and -0.018 [-0.039, -0.002] within strata, but -0.001 [-0.017, +0.018] pooled. Three of 20 per-fold intervals exclude zero, about one is expected by chance, and none survives Holm's correction (smallest adjusted p 0.18) | docs/twin_comparisons.md; docs/graph_content_null_results.md; docs/neuronal_variant_comparison.md (7 October) |
| 10 | Untrained, the linear response predicts the direction of measured metabolite changes worse than always guessing "increased" | 1,067 signed gene-metabolite pairs from HPO laboratory abnormalities; sign agreement 0.615 at best (8 steps, summed over compartment copies) against 0.806 for the majority direction; AUROC 0.574 | docs/linear_response_sign_check.md (7 October) |
| 11 | Candidate cause of the null, not yet tested: the encoder's per-node averaging cancels production against consumption | 7,811 of 21,337 destinations (36.6 percent) receive both signs; at equal gains the aggregate input sums to zero at all 7,811 metabolites that receive both signs (92 percent of 8,460; first recorded as 88 percent, an undercount) and at no reaction. It predicts that keeping the stoichiometric counts moves the score. The spectral arm keeps them but also shrinks the field nearly everywhere (running, one fold of five); the total-in-degree arm keeps them at the in-degree magnitude and is queued. The cross-relation mixture cannot break an exact cancellation, since every statistic in it is odd | docs/membrane_potential_reach.md; docs/graph_content_null_results.md (7 October) |
| 12 | The disease-cluster split divides almost every pathway between training and test, and the split pathways carry label information across it; protein complexes do not | of 308 held-out genes whose primary subsystem is a pathway, 271 (88 percent) have a training gene with the same one and 37 have none; a held-out gene's mean Jaccard similarity of positive symptom sets to its same-pathway training genes exceeds its mean to all training genes by +0.050 [+0.028, +0.072] (bootstrap over 271 held-out genes), to genes sharing any pathway subsystem by +0.028 [+0.013, +0.044], to genes sharing a catalysed reaction by +0.020 [-0.014, +0.055] and to complex partners by -0.011 [-0.069, +0.050] (58 genes, 13 percent, have a complex partner in training) | docs/fold_mechanism_sharing.md (7 October) |
| 13 | The trained models trail the random walk only where a gene's pathway is split across folds; elsewhere they tie it | split-pathway genes (271): B3 typed nodes minus random walk -0.049 [-0.075, -0.023], B6 linear response with gate time scale -0.037 [-0.068, -0.005], B3 linear response with cofactors -0.073 [-0.099, -0.048], B6 mechanistic with gate time scale -0.022 [-0.051, +0.023]; unsplit genes (180): every trained model within an interval including zero (B3 typed nodes -0.000 [-0.037, +0.034]). The random walk's own lift over popularity is +0.072 on split-pathway genes, +0.076 on the 37 genes whose pathway is absent from training (5 symptoms scorable, intervals near ±0.1) and +0.041 on the 143 genes with no pathway subsystem, so its edge follows being a pathway enzyme rather than the split itself. The random walk reads the labels of training neighbours at prediction time and the trained models do not, which is the candidate reason; pooled predictions, strata by primary subsystem, no correction for the number of comparisons | docs/score_by_pathway_split.md (7 October) |
| 14 | Scored only on the genes it gives a distinct prediction, the noisy-OR model still does not beat the random walk on those genes | per fold, a held-out gene is undecided when its prediction is within 1e-6 of the largest tied group's, so the model makes no ranking among those genes. B6 linear response with the gate learning rate decides 176 of 451 genes (39 percent): 48 percent of split-pathway genes against 21 percent of genes with no pathway subsystem, median metabolic-graph degree 4 against 2. On the decided genes it scores 0.251 against the random walk's 0.278 on the same genes, paired -0.027 [-0.060, +0.014], and beats popularity by +0.040 [+0.026, +0.077] and B3 linear response with cofactors by +0.045 [+0.019, +0.084]. B6 mechanistic with the gate learning rate decides 376 (83 percent): -0.017 [-0.040, +0.019] against the random walk on its decided genes; on the 75 it leaves undecided the random walk scores 0.351, a lift of +0.168 over popularity. Uncorrected intervals; the decided set is chosen by the model's own outputs, without labels | docs/score_by_module_coverage.md (7 October) |
| 15 | The random-walk baseline is undirected while both encoders follow edge direction; following direction neither costs nor adds to the walk | the encoders put each edge at [target, source] (the linear response adds only the depletion edge in reverse); the walk and every rewiring in row 8 symmetrise the graph. Three unsigned walks with one shared readout, minus the undirected one: the encoder's edge directions +0.005 [-0.003, +0.013] per fold, +0.006 [+0.001, +0.011] pooled, -0.002 [-0.011, +0.009] within degree strata; downstream only +0.006 [-0.015, +0.027], +0.008 [+0.003, +0.014], +0.012 [-0.002, +0.028]. One reading of three excludes zero for each, uncorrected. The sign ablation of B3 linear response with cofactors (every edge +1) is running | docs/directed_random_walk.md; docs/graph_content_null_results.md (7 October) |
| 16 | A readout trained on the walk's own distributions does not beat the untrained walk, so the trained models' deficit is weak evidence against their mechanism | per-symptom L2 logistic regression on the walk distribution (square-rooted, 32 SVD components fitted per fold), minus the walk: undirected -0.012 [-0.031, +0.006] per fold, -0.026 [-0.044, -0.011] pooled, -0.004 [-0.026, +0.015] within degree strata; downstream-only walk +0.001 [-0.036, +0.037], -0.020 [-0.037, -0.002], +0.005 [-0.018, +0.023]. One readout design, fixed before scoring; uncorrected intervals | docs/trained_diffusion_readout.md; docs/graph_content_null_results.md (7 October) |

Data layer:

- Graph: 33,964 nodes and 255,367 typed edges (Human-GEM 2.0.1 with compartments, OmniPath, CollecTRI regulons
  restricted to brain-expressed transcription factors, GTEx v10 brain expression on gene nodes).
- Labels: 2,425 grade A gene pairs, 876 grade B drug pairs over 142 drugs, 609 grade C non-causal associations;
  the full-data baselines score 1,417 perturbations (1,275 genes, 142 drugs) in 901 disease clusters.
- Literature priors: 214,375 PubTator3 and CTD reports over 134,747 papers, grade E, never labels.
- Release: data/releases/v0.4 with a sha256 manifest.

Superseded numbers to avoid:

- The first full build of 3 October (1,585 perturbations, random walk 0.326, rewired-graph margin 0.011) predates the
  OnSIDES statements and the association types; rows 1 and 2 replace it.
- B6 mechanistic 0.272 with a paired difference of -0.029 [-0.051, -0.010] against the random walk was the best
  configuration on 6 October; row 5 replaces it as the headline, though that configuration's numbers still stand.
- "65 of 23,709 destinations (0.3 percent) see both signs" was computed from the stored edge table without the
  encoder's derived depletes_substrate relation; the correct figure is in row 11.
- 0.269 as the slice random walk's score on the real graph was a misread of the permuted-label table; the real graph
  scores 0.293.
- Twin intervals from a percentile bootstrap over the five fold differences (manganese [-0.005, +0.002], neuronal
  variant +0.012 [-0.012, +0.036]) were about a third too narrow; row 9 gives t intervals. "The one clear effect" for
  the protein descriptors is withdrawn for the same reason and for the 20 comparisons.
- "Rewiring leaves the slice score unchanged" came from the 2-swaps-per-edge control and is withdrawn (row 8).

## 4. Figures

| Figure | Content | Status |
|---|---|---|
| F1 | Data and model schematic: evidence sources, graph layers, perturbation encoder, modules, noisy-OR output | not made |
| F2 | Leakage cascade: random walk under gene, disease-cluster and permuted-label conditions (finding 1) | numbers ready |
| F3 | Wiring against degree: real graph against the 20-rewiring null (finding 8); the full-data panel waits for the rerun of finding 2 | slice ready |
| F4 | Main comparison with paired-bootstrap intervals, slice then full data (finding 5 and Aim 2) | slice ready |
| F5 | Per-symptom AUPRC against base rate | slice ready |
| F6 | Module diagnostics: gate travel against its bound, support size, link matrix, sufficiency test (finding 6) | slice ready |
| F7 | Time-split results (finding 4 and Aim 3) | monogenic ready |
| F8 | Sign agreement of the untrained response against the majority direction (finding 10) | ready |

## 5. Likely objections and current answers

| Question | Current answer or where it lives |
|---|---|
| What if the mechanistic model does not win? | Aim 1 stands on its own; the failure criterion was fixed in advance; Aim 3 has a branch for it |
| How do you know the evaluation does not leak? | disease-cluster grouping, drugs grouped by dominant target, permuted-label and rewired-graph controls (design 6.1, 6.3). The split does divide pathways, and split pathways carry label information (finding 12); the subsystem hold-out is the test that removes it (finding 3) |
| HPO annotations are incomplete; are missing pairs negatives? | no: unobserved pairs are unlabelled and sampled at reduced weight (design 5.4) |
| Why noisy-OR? | it encodes equifinality without interaction terms (design 5.3) |
| Why not flux balance analysis? | the brain has no defensible single objective; sampling is the planned route, with the coverage limits of finding 7 |
| Why these symptoms and not diagnoses? | assumption A7 and the symptom crosswalk (docs/symptom_crosswalk.csv) |
| Why does degree explain so much? | hub structure in the graph and in the labels (design assumption A9); on the slice the wiring adds signal beyond degree (finding 8), and the full-data control needs rerunning (finding 2) |
| Why were the modules unused, and are they used now? | the gates were still shrinking when early stopping ended training; a separate gate learning rate fixes that, and the modules now carry links but none passes the sufficiency test (finding 6) |
| The graph carries signal, so why does the mechanistic encoder not beat diffusion? | open; finding 11 is a candidate cause, tested by the spectral normalisation and cross-relation mixture arms. The literature search has not yet found a published case of signed propagation on a mechanistic network beating unsigned diffusion for a phenotype (docs/literature_appraisal_staging.md, question Q1) |
| Does the graph give the right direction of change? | untrained, no better than the majority direction (finding 10) |

## 6. Decisions and provenance

| Date | Decision | Reason |
|---|---|---|
| 2026-10-03 | Disease-cluster grouping is the default split (design v0.4) | gene-wise grouping leaked through shared diseases (finding 1) |
| 2026-10-06 | PubTator3 pinned to the 6 October bulk file | NCBI no longer serves the 17 August file |
| 2026-10-06 | Flux route deferred until a B4 test shows signal (proposed, not decided) | coverage and medium problems of finding 7 |
| 2026-10-07 | Gates get their own learning rate (--gate-learning-rate) | with the shared rate, supports were still shrinking at early stopping (finding 6) |
| 2026-10-07 | Rewiring control at 50 swaps per edge with a null of 20 draws on the slice | 2 swaps per edge left the graph under-mixed (finding 8) |
| 2026-10-07 | The 0.3 percent sign-mixing figure retracted | computed without the derived depletes_substrate relation (finding 11) |

AI assistance: code, data processing, reviews and drafts in this repository were produced with Claude Code; commits
carry a Co-Authored-By line.
