# Soft priors: the evidence below the label bar as a target and a weight

Record, written 16 October 2026. The design of `mechanistic_pathway_learning/evidence/soft_priors.py` and its wiring
into `experiments/run_main_model.py`, the reasons for each choice, what is provisional, and the tests that would
have to pass before the arm means anything. The research pass behind it (16 October 2026, "Curated-evidence soft
priors for a positive-only physiology graph", in the session transcript) is quoted where a choice rests on it. Open
decision 9.50 is the log entry; this page does not say what is current.

## 1. What the arm does

A label (grade A or B) trains a pair toward the positive target at its evidence weight. Everything below that bar
said nothing to the loss until now: grade-C records were set aside by the label selection, grades D and E never left
the report table, and the three sources that record absences (placebo-controlled trial arms, FAERS disproportionality,
IMPC knockouts; `docs/open_decisions.md` 9.48) had no path into training at all. With `--soft-prior-strength` above 0
each such pair on a row the fit trains on becomes

    target = the tier's prior probability (positive evidence) or 0 (an absence)
    weight = strength * unit_strength(tier) * min(n, cap) ** tempering

with `n` the number of independent records behind the pair: distinct PMIDs where the reports cite one, else distinct
source records (a trial, a knockout line). The loss is the weighted binary cross-entropy already in use; with a
fractional target `q` the per-pair loss is the expectation of the hard loss over a Bernoulli(`q`) truth, which is the
proper way to say "probably positive" without saying "positive". The weight is a power prior: Ibrahim et al. (2015,
*Stat Med* 34:3724) describe the power prior as one that "is intended to be an informative prior constructed from
historical data", and its effective sample size is the exponent times the historical count, which is what
`strength * tempered count` is here.

The tiers are in `configs/soft_prior_tiers.csv` and are keyed by evidence class and direction:

| class | direction | prior probability | unit strength | the evidence |
| --- | --- | --- | --- | --- |
| monogenic | positive | 0.75 | 0.5 | an HPO annotation below the grade-A bar (grade C) |
| human_association | positive | 0.70 | 0.5 | a trial arm above placebo, or a FAERS signal (PRR >= 2, chi-square >= 4, >= 3 reports) |
| human_association | negative | 0 | 0.5 | a trial arm with no excess at >= 100 exposed, or a FAERS absence at >= 2,000 drug reports |
| animal | positive | 0.60 | 0.3 | an IMPC call on a crosswalked MP term |
| animal | negative | 0 | 0.3 | an IMPC line tested on the crosswalked procedure without a call |
| literature | positive | 0.55 | 0.1 | a PubTator3 cause, stimulate or inhibit relation, or a CTD direct-evidence row |

An absence trains toward 0 at the tier's strength, as the registered placebo rule does (amendment item 3 gives an
evidenced negative 0.5 where the registered run gives silence 0.2); `read_tiers` refuses a negative tier with a
positive target, so the two rules cannot drift apart by a typo.

## 2. The rules, and why

**Discrete tiers by evidence type, not continuous weights by count.** Every count in this project's inputs measures
how much a pair was studied at least as much as whether it is true. Schaefer et al. (2015, *Front Genet* 6:260) found
that "the degree of a protein is a function of the number of times it has been screened for interaction partners"
(r = 0.52), and Gillis and Pavlidis (2012, *PLoS Comput Biol* 8:e1002444) that "about one-half of performance could be
attributed entirely to node degree effects" in guilt-by-association prediction. A weight that rises with a count would
raise the prior on the same perturbations that already carry the most labels, which is the popularity confound the
design's degree strata exist to control. ClinGen's gene-disease validity scheme (Strande et al. 2017, *Am J Hum
Genet* 100:895) is the template: discrete tiers, capped accumulation so that "a preponderance of weak evidence" cannot
inflate a classification, and a replication requirement for the top tier.

**A capped, tempered count.** The one step a count supports is the step from one independent publication to two:
Cusick et al. (2009, *Nat Methods* 6:39) re-curated singly supported interactions and found that "curation can be
error-prone and possibly of lower quality than commonly assumed", a finding disputed in its size and not in its
direction. So the cap defaults to 2 and the exponent to 0.5: the second record multiplies the weight by 1.41, the
fiftieth adds nothing. Both are flags (`--soft-prior-count-cap`, `--soft-prior-tempering`); exponent 0 is an
indicator, which the research pass recommended ("0 / 1 / >= 2 deduplicated PMIDs").

**Record-level deduplication.** OmniPath's 2025 update reports that "only 22 % of the records are unique to a single
resource", and the same papers feed HPO, CTD and PubTator. `record_key` counts a pair's support by PMID where the
report cites one, else by source record, so two arms of one trial, or one paper in two sources, are one record.

**One class decides a pair.** The highest-ranked class present (`assemble_evidence_table.EVIDENCE_CLASS_RANK`) sets
the pair's target and strength; classes are not averaged, because a human source and a mouse or a sentence are not
exchangeable evidence. A pair with evidence in both directions across the tiered classes is contested and takes no
prior (`ResolvedPriors.contested` lists them and `results.json` counts them).

**Where a prior may land.** Never on a labelled positive, a labelled negative, or a pair the selection set aside for
a controlled comparison (the never-above-placebo pairs); on an unlabelled pair, or on a pair set aside only for
lacking a label grade (`ExperimentData.masked_for_soft_evidence`, written by the loader from the selection's
`masks_a_negative` rows). The counts of priors that fell on ineligible cells are reported, not silently dropped.

**Fit rows only, and the lockbox never.** The priors are built after the lockbox is removed from the data, so no
lockbox perturbation receives one; they are written into the targets of the rows the first fit trains on and, at the
refit, of the training and validation rows together; the validation loss, the early stopping and every score read the
plain matrices. Under `--permute-labels` the priors move with the label rows, so the permuted-labels control sees no
pair-specific information from them either.

**Literature: induction relations only, and CTD direct evidence only.** The co-mention relation ("associated_with")
is not a direction and is left out; CTD's inferred associations are transitive closures over gene-level edges the
graph already holds and would count the graph twice, so only its curated rows (which all carry a direct-evidence
code) enter. PubTator3's end-to-end relation extraction scores 37.6 F1 in an independent evaluation (Almeida et al.
2024, *Database* baae057), which is why the literature tier is the weakest.

**The Dawid-Skene posterior is not the strength.** With sources that only ever report positives, the fit's
"sensitivity" is the probability that a pair was studied and curated, not that it is true, and in the v0.4 table the
posterior ranks every pair within a band of 0.015 (`docs/evidence_reports_spec.md`, decision 18). It stays in the table
as a diagnostic (`reliability_posterior`) and sets no weight.

**No LLM appraisal.** Agreement of language models with human raters on structured risk-of-bias grading is weak
(kappa 0.06 to 0.22 in the 2025 comparative studies the research pass read), no study validates such grades as
priors, and a model trained on later literature knows which claims replicated, which is the labels' own signal. An
appraisal of design facts (species, dose, knockout against knockdown) would be extraction, not grading, and is not
built.

## 3. The controls the arm must pass

A prior that is only popularity will raise scores on the well-studied perturbations and nowhere else. Zietz et al.
(2024, *GigaScience*, the XSwap edge prior) showed a degree-only prior reconstructing one network at AUROC 0.797 while
it "was not able to predict the separate PPI network better than by random guessing". So the arm is read against:

1. **The registered model** (`--soft-prior-strength 0`, bit-identical to the run before the flag existed; tested).
2. **The degree-popularity baselines** already in the confirmatory set.
3. **The permuted prior** (`--permute-soft-priors`): the same priors moved among eligible cells so that every
   perturbation and every symptom keeps its number of them; a prior's target and strength travel with it. The real
   arm must beat it, and the gain must hold in the lowest popularity bins, which `results.json` now reports as
   `macro_auprc_by_literature_popularity_bin` (distinct PMIDs per perturbation in the literature table, five
   equal-count bins over the test rows) beside the degree bins.

The batch runner names the arms (`v4_linear_response_soft_priors`, `_permuted`, and the message-passing twins) and
`configs/ablations/all_ablations.yaml` lists them as `soft_priors_tiered` and `soft_priors_permuted`. Nothing here is
a primary reading; the registered defaults are unchanged (`tests/test_registered_defaults.py`).

## 4. What is provisional

- **The tier numbers.** 0.75 / 0.70 / 0.60 / 0.55 and 0.5 / 0.5 / 0.3 / 0.1 are a first setting by evidence type,
  ordered as the research pass orders the classes, and not fitted. They are to be chosen on development data before
  any registration: a grid over the global strength (0.25, 0.5, 1) and the exponent (0, 0.5, 1), scored on the
  development folds against the three controls, with the permuted arm run at the chosen point. A strength chosen on
  the lockbox would be a reading of the lockbox.
- **Perturbations with soft evidence only** are not in the perturbation set (the loader keys the set on the labels),
  so the arm fills in symptoms of labelled perturbations and adds no rows. Adding them changes the splits and the
  lockbox checks and is a separate decision.
- **Contested pairs** are skipped; a rule that lets a human absence override a mouse call (or the reverse) is a
  decision, not a sort order.
- **Combining classes** (a trial and a knockout agreeing) adds nothing today; ClinGen-style capped accumulation across
  classes is the natural extension.

## 5. Tests

`tests/test_soft_priors.py` (15): the committed tier table, refusal of a positive absence target, tempering and the
cap, the record key, report evidence by direction with co-mention dropped, record evidence by grade, resolution by
rank with contested pairs set aside, summation across tables before tempering, eligibility, matrix alignment,
application to fit rows only, strength 0 as a no-op, the permutation's margins and determinism, literature popularity,
the summary. `tests/test_soft_priors_in_the_trainer.py` (6): the defaults and the registration map, the fingerprint
naming the prior only while it is on, the fit targets, the loader's `masked_for_soft_evidence`, `build_soft_priors`
on a toy dataset, and an end-to-end run of the trainer three times (plain, strength 0 with every flag set, strength 1,
then permuted with permuted labels) checking bit-identical predictions and an unchanged fingerprint at strength 0 and
the recorded summaries at strength 1.
