# Options for the confirmatory pathway test

Written on 8 October 2026, before any second-family run and before any module of the confirmatory configurations
exists. The user asked for a summary of the options. Nothing here is decided or registered; the preregistration lists
the pathway test under "Not decided". The review read no lockbox output and no second-family score. Sources are the keys
in docs/references.bib; each entry's note carries the quote that was read.

## What a pathway claim needs

The user's statement of purpose: "It's to be able to use the graph to determine what the pathways for symptoms are."
A module here is a node set (the hard support, gate above 0.5) with a link vector over symptoms. A claim that module k
is a pathway for symptom s needs three separate things:

1. **The model uses the support for s** (faithfulness). Agreement with outside data does not show this (deyoung2020eraser:
   "the degree to which provided rationales influenced the corresponding predictions").
2. **The support recurs across seeds and training sets**, and more often than for modules trained on permuted labels.
   Recurrence is not identity in general (locatello2019challenging: "well-disentangled models seemingly cannot be
   identified without supervision").
3. **Evidence for s that the model never trained on**, beyond what proximity to the training genes predicts.

Stability and faithfulness alone say nothing about biology: the slice modules were stable (cross-split Jaccard 0.25 to
0.74) while "the eight modules of each fold concentrate on the same hub nodes". Outside evidence alone is open to two
known artefacts: enrichment that also appears on modules from permuted data (levi2021domino: "GO terms enriched in
modules detected on the real data were often also enriched on modules found on randomly permuted data") and modules
that follow node degree (lazareva2021limits: "do not produce biologically more meaningful candidate disease modules on
widely used PPI networks than on random networks with the same node degrees").

## What the design already has

- Section 6.4: MAGMA competitive tests on the top-N genes by attribution, secondary. The wrapper is a stub, and the
  preregistration still lists the GWAS power statement under "To fill in".
- Section 6.5, the "sufficiency test": removes one module from the noisy-OR using saved activations, on perturbations
  grouped by that module's own largest contribution. It is a necessity test at the head, the grouping is circular, and
  it never touches the encoder or the graph.
- Stability selection (design line 189): experiments/analyze_pathway_modules.py takes a greedy best match of Jaccard
  across splits, with no one-to-one matching and no null.
- Section 6.6: overlap with the 16 curated modules and the Human-GEM subsystems.
- Null modules at no extra cost: the lockbox job's permuted-label and rewired-graph runs (2 models x 5 seeds each).

## Ranked options

| rank | option | question | data (never trained on) | null | pass rule (per model) | cost |
|---|---|---|---|---|---|---|
| 1 | Conditional MAGMA gene-property test | Do the genes the modules weight for s carry more GWAS association with the nearest trait than graph genes with the same covariates? | Pan-UKB item GWAS (Tier A); PGC (Tier B) | covariates: gene size, SNP density, log degree, LOEUF, brain specificity, B1 random-walk score from s's training positives; same statistic on permuted-label and rewired modules | Stouffer Z > 1.96 over Tier A symptoms; the controls < 1.645 | one gene analysis per GWAS, several large downloads |
| 2 | Deletion and keep-only of the support | Does the prediction of s on lockbox perturbations depend on the support as a route through the graph? | lockbox labels, inside the scoring step | 200 to 1,000 random node sets matched on size, type, degree and distance from the seeds | Holm p < 0.05 for at least one (k, s); unlinked symptoms show no drop | inference only, hours |
| 3 | Stability with one-to-one matching | Do the same supports recur across the 5 seeds and resampled training sets? | existing runs | permuted-label and rewired runs; random supports matched on size and degree | at least one responding module matched at Jaccard 0.5 or more in 4 of 5 seeds; real mean above permuted | free for seeds |
| 4 | Mouse knockout behaviour (IMPC; MGI secondary) | Do module genes more often carry the knockout call mapped to s than calls for unmatched terms? | IMPC release 24.0, CC BY 4.0 | matched versus unmatched (freund2018phenotype); covariates degree, phenotype count, viability | odds ratio > 1, 95 percent interval above 1 | low |
| 5 | Planted modules on the real graph | Does this training procedure recover known supports and links at this label density? | simulated labels from an independent field (random walk) | permuted simulated labels | median matched Jaccard 0.5 or more; link AUROC 0.8 or more | 6 to 12 full-graph runs |
| 6 | Drugs outside training (FAERS, OFFSIDES) | Do drugs whose targets lie in a support linked to s have more reports of s? | openFDA (CC0) | matched drugs; shrunk PRR or propensity-matched estimates | Holm one-sided p < 0.05 for a symptom with 20 or more drugs | moderate |
| 7 | Loss-of-function burden (Genebass; SCHEMA) | Is pLoF burden in module genes associated with the matched UKB item? | Genebass; coverage of the items not verified | random sets matched on LOEUF, length and degree | pooled burden Z, one-sided p < 0.025 | low to moderate |
| 8 | Curated pathway overlap | Which known pathways does each support overlap? | GO BP, KEGG, Reactome outside the 15 imported | random connected sets; permuted and rewired modules | above the null's 95th percentile | low; descriptive only |
| 9 | Literature time split | Do module (gene, symptom) pairs gain later papers more often than matched genes? | a future PubTator3 file | genes matched on prior publication count (stoeger2018large) | rate ratio > 1 | clean only prospectively |

Notes on the ranked options:
- **Option 1** replaces the top-N set of section 6.4 with a continuous weight per gene,
  w_s(g) = max over k of link(k, s) x E[gate_k(g)], averaged over the 5 lockbox seeds. The weight does not depend on module
  order, so it needs no matching across seeds. lamparter2016pascal argues against a cut: "eliminates the problem of
  arbitrary threshold selection inherent in any binary membership based pathway enrichment approach". The conditioning
  answers deleeuw2018conditional: "one property may be statistically associated even if biologically irrelevant to the
  phenotype, by being correlated with gene properties that are relevant". The universe is graph gene nodes minus every
  development training seed gene, because Mendelian genes sit near the GWAS loci of matched traits (freund2018phenotype:
  "a significant enrichment of genes linked to phenotypically matched Mendelian disorders in GWAS gene sets").
- **Option 2** compares with matched random deletions because deletion without retraining mixes two effects
  (hooker2019benchmark: "it is unclear whether degradation in performance is due to the introduction of artifacts outside
  of the original training distribution or because we actually removed information"; zheng2024robust gives the graph
  form). Sanity checks from adebayo2018sanity: re-initialised encoder weights must change the result, and permuted-label
  modules must fail. Seeds are excluded from the deleted set, and each module's seed share is reported.
- **Option 3** uses Hungarian matching (kuhn1955hungarian), selection frequency with the threshold of
  meinshausen2010stability, complementary-pair subsamples (shah2013variable) and the stability index with an interval of
  nogueira2018stability. Splits and merges of modules are penalised by one-to-one matching, so the order-free correlation
  of w_s across seeds is reported beside it.
- **Option 4.** Symptom to MP term, from the review's OLS4 search (to verify before registration): anxiety MP:0001362,
  hyperactivity MP:0001399, psychosis via abnormal prepulse inhibition MP:0003088, cognitive impairment MP:0002063,
  insomnia and somnolence via MP:0011396, irritability or aggression MP:0001354, increased appetite MP:0011939. Live IMPC
  counts: anxiety 371 genes, prepulse inhibition 468, hyperactivity 835, learning and memory 591, sleep 80. IMPC calls are
  mostly independent of the literature behind HPO (meehan2017disease: "90% of our phenotype annotations were novel"). The
  limit is construct validity (nestler2010animal: "extremely challenging given the subjective nature of many symptoms").
- **Option 8** is circular for the 15 imported Reactome pathways and the Human-GEM subsystems (a connected support
  overlaps them by construction) and for GO function and component (the ESM-2 descriptors are fitted onto them). It names
  modules; it confirms nothing (timmons2015multiple: "extreme caution should be applied when using or interpreting
  functional enrichment analysis").
- Other items a reviewer may expect: cross-class transfer (modules from gene perturbations predict the lockbox drugs);
  a held-out symptom fitted on frozen supports; a timestamped registry of module gene lists scored on later HPO and OnSIDES
  releases; a rule that a support made only of training seeds is a lookup, not a pathway.

## The review's recommended rule

Option 1 as the confirmatory pathway test, read only if options 2 and 3 pass as gates, with option 4 as the
preregistered secondary. Per model at one-sided 0.025:
- **G1 (stability):** at least one responding module (docs/module_health.md) is matched at Jaccard 0.5 or more in 4 of
  5 seeds, and the real-seed mean matched Jaccard is above the permuted-label one.
- **G2 (faithfulness):** deleting that module's non-seed support lowers the lockbox likelihood of its linked symptom
  (link 0.1 or more) beyond 95 percent of 1,000 matched random deletions (Holm).
- **P1 (outside evidence):** the Stouffer Z of the conditional MAGMA coefficient of w_s over Tier A symptoms is above
  1.96, and the same Z for the model's permuted-label and rewired-graph modules is below 1.645.
- **S1 (secondary):** the IMPC matched-versus-unmatched odds ratio with a 95 percent interval above 1.

GWAS coverage of the 16 symptoms in the macro average (review's tiers):
- Tier A, item-level GWAS (10): depressed mood (field 2050), anxiety (GAD symptoms 2026, or 1980 and 1970), insomnia
  (1200), somnolence or hypersomnia (1220), fatigue (2080), psychomotor agitation (2070 or 20516), emotional lability
  (1920), irritability or aggression (1940), disinhibition or impulsivity (2040, risk taking, a weak proxy), self-injury
  (20480, 5,880 cases).
- Tier B, diagnosis proxy (5): psychosis (scz2022), compulsive behaviour (ocd2025), hyperactivity (adhd2022), cognitive
  impairment (cognitive-function GWAS or 20508), increased appetite (binge eating 2026 or the bidirectional 20511).
- Tier C, none adequate (1): abnormal dreams (ollila2024nightmares: "the GWAS did not reveal individual risk variants").

## Points against the rule as written

1. **Preregister it now.** No confirmatory module exists yet, so the rule can be fixed before anyone sees one. Once the
   v2 modules exist, every threshold above is open to choice after the fact.
2. **The gates are preconditions, not tests at a level.** The claim passes only if G1, G2 and P1 all pass. If each is
   a test at level 0.025, the conjunction keeps the level 0.025 (an intersection-union test); a precondition can only
   lower the chance of passing. G1 as written compares two means of 5 seeds with no test, so it is a precondition only.
3. **Multiplicity against H1 and H2 is the user's decision.** Three ways, per model:
   - H1, then H2 and P1 by Holm at 0.025. The chance of a false claim per model stays at most 0.025 and across the two
     models at most 0.05. H2 then needs p below 0.0125 unless P1 also passes.
   - H1, then H2, then P1 in fixed sequence. H2 keeps 0.025; P1 is read only if the real wiring beats the rewired one.
   - P1 on its own at 0.025. The chance that at least one of the four claims (2 models x H1 and P1) passes by luck rises
     to at most 0.10.
   The first or second keeps the stated bound. Gating P1 on H1 fits the purpose: modules of a model that does not predict
   better than the baselines are not worth reading.
4. **Power comes before choosing P1 as primary.** Supports of 1 to 20 genes leave few genes with weight, and item GWAS
   can be weak (deary2018tiredness: "GWAS identified one genome-wide significant hit"). A power simulation should not
   open the P1 GWAS: draw null gene Z with the reference panel's gene correlations and plant a shift. If power is low at
   a plausible effect, a P1 failure says nothing, and IMPC (option 4), the closest analogue of the model's knockouts,
   is the stronger candidate for co-primary.
5. **No development reading of the P1 GWAS or the IMPC calls.** Running MAGMA on development modules against the Tier A
   list would spend the test. If a development check is wanted, split the GWAS list in advance (choobdar2019assessment
   split "180 GWASs" into "a leaderboard set" and "a separate holdout set").
6. **Drug-target circularity in depression.** pgc2025mdd: "The associations are enriched for antidepressant targets".
   Training drug targets leave the universe with the seeds; the depressed-mood test uses the item GWAS (2050), not the
   case-control meta-analysis.
7. **G2 reads lockbox labels, so it belongs to the single scoring.** experiments/score_confirmatory.py is frozen, so G2
   needs its own script, written and tested on development folds before any lockbox run, then run once after
   runs/confirmatory/SCORED. P1 and S1 read no lockbox outcome but read module files under the blinded directories, so
   they also run after SCORED.
8. **Option 5 is method validity and costs full-graph runs.** It tells whether the procedure can recover a planted
   module at a base rate of 0.064 at all. Full-graph training outside the registered jobs needs the user's approval; a
   slice version is cheap.
9. **None of these shows the order of a route.** A support is a node set. The claims are "the model uses these nodes for
   s" and "these nodes carry outside evidence for s", not "signal flows A to B to C".

## Decisions for the user

1. Which pathway test to register: the review's rule (P1 primary, IMPC secondary), IMPC as co-primary or another.
2. How the pathway test sits beside H1 and H2 (point 3).
3. Whether to run the power simulation (point 4) and the slice version of option 5 before registering.
4. Whether to restrict P1 to the 16 symptoms of the macro average (the 6 others have few positives for the links).

## Data access and terms

- Pan-UKB: reachable; manifest confirms the item fields above; licence page renders by script, not verified.
- PGC: "Use of these data is NOT unrestricted"; non-commercial without permission and no redistribution. Kept under
  data/raw, never committed.
- IMPC and MGI: CC BY 4.0. openFDA: CC0. gnomAD constraint, Genebass, nSIDES: licence not verified.
- MAGMA v1.10: "standard copyright … may not be distributed or modified"; v1.0 is GPL v3. The binary stays out of the
  repository.
- LDSC reference files: the old URLs return 404.
