# Draft of an amendment to the pre-registration: the evidence table v4 and the checks

**A draft for the owner. It is not applied.** `docs/preregistration.md` is unchanged, no registered configuration reads
anything named here, no held-out set has been drawn on the new table and no model has been trained on it. The draft
states what the decisions of 10 and 11 October 2026 would change when the owner applies them, and which decisions are
still missing. Written on 11 October 2026.

The registered state it would amend: the evidence table `evidence_full_v2`, the `better_v2` label selection,
`configs/lockbox_v2.json`, and the graph `data/processed/graph_full_neuronal_split_binders` (amendment of 9 October
2026, fourth).

## 1. What the owner has decided

Each line gives the owner's words and where they are recorded (`docs/research_summary.md`, decision log).

1. **Three label changes** (10 October: "Up to three targets", "Fix gain-of-function sign", "Mask never-above-placebo";
   11 October: "I accept all the suggestions under section 8 of the label review."; "Please do the recommended label
   fixes").
2. **Where the added drugs go** (11 October: "I accept put the three mechanism target pairs being mostly in the held
   out set"): an added drug that acts on a node of a held-out perturbation is held out, as its own stratum outside the
   primary score (`docs/label_source_review.md`, section 9.2, option a).
3. **Combined statements are checks** (10 October: "Build, report as checks"; 11 October: "Add the valproate, etc
   displacements as checks then since those are important and survived"; "Add the enzyme labels as the carrier twin
   instead then").
4. **A reading at the peak of the unbound tissue concentration** (11 October: "Peak tissue unbound concentration can
   be relegated to confirmation only").

## 2. What would change

- **The evidence table.** `data/processed/evidence_full_v4` in place of `evidence_full_v2`, built by
  `scripts/rebuild_data_layer.sh` with `--max-drug-targets 3` and `--seed-gain-of-function-genes-positive`; the
  selection `label_selection/better_v2_full_v4.parquet` with `--mask-grades C --mask-never-above-placebo`. SHA-256,
  first 16: records `e0206dbe181747e8`, reports `50d847913fa8ef9f`, selection `1d0cde3b1292217c`. 1,636 perturbations
  (1,414 genes and 222 drugs), 24 symptoms, 2,975 kept pairs of which 1,076 are drug pairs. Against `evidence_full_v2`
  the table also carries two earlier changes that are not registered either: drugs with a target of any type and the
  salt-form fix (the table called v3, `docs/drug_targets_any_type.md`), and parkinsonism as the 24th symptom (decision
  of 9 October). The amendment would register all of them together. 316 of the 317 perturbations of
  `configs/lockbox_v2.json` are in the table.
- **The sign of 17 genes.** A gene whose every causal association is a gain of function is seeded with sign +1. The
  registered text discloses the old sign and leaves it ("so the signed encoder reads the wrong sign for them"); the
  amendment would replace that sentence.
- **The held-out set.** The membership of `configs/lockbox_v2.json`, carried to v4 by its leakage groups (320
  perturbations, 31 drugs), stays the set the primary score is read on. The 29 added drugs that act on a held-out
  node (`configs/held_out_additions_v4_draft.csv`; 160 kept pairs) are held out and reported as their own stratum.
  They do not train, do not validate and do not enter H1 or H2. The 39 other added drugs train. Under the rule 349 of
  1,636 perturbations, 60 of 222 drugs and 307 of 1,076 kept drug pairs are held out. The registered group rule chains
  and is not used to place the added drugs, because it would pull 149 further perturbations into the held-out set;
  the amendment would say so and state the one-step rule in its place for added perturbations.
- **Checks.** `configs/comparative_checks.csv` (214 statements) is read by `experiments/run_comparative_checks.py`
  after scoring. The checks are named in advance, test no hypothesis of the study, change no decision rule and are
  reported apart from every score. A finished development run is asked the development checks only.
- **The peak reading.** Not built. The owner's sentence is read here as: a reading at the peak concentration would be
  reported at the confirmation stage, if at all, and is not trained on or selected on. The owner should correct this
  reading if it is wrong.

## 3. What would not change

The two tested models and their heads (`configs/head_choice.json`), the decision rule, the one-sided 0.025 sequence of
H1 then H2, the floors, the degree strata, the confirmatory graph and the node descriptor table. No drug-node graph and
no clearance edge is registered: `graph_full_neuronal_split_binders_drugs_v4` and its variant with clearance edges
exist for the ablation the owner asked for on 10 October and for the checks, and whether a confirmatory configuration
reads either is a separate decision.

## 4. What is still missing before this can be applied

1. **The limit on a drug's targets.** Three stays the rule until the owner decides (`docs/label_source_review.md`,
   section 12). Without a limit eight drugs enter and one of them joins the two large drug groups into one unless a
   written rule sets it aside.
2. **The placebo pairs.** Set aside, as decided on 10 October, or counted as negatives, as the owner's question of
   11 October suggests ("Why would working no better than the placebo not be counted as evidence against a label?").
   Both selections are built; the second is `better_v2_full_v4_placebo_negative.parquet`, SHA-256 first 16
   `3f4decacd783ab08`.
3. **The draft list of admitted drug groups** (`configs/admitted_atc_prefixes_draft.csv`): its rows need the owner's
   sign-off, and one drug of the list joins the two large drug groups in the same way.
4. **The input of a multi-target drug** and the other questions of `docs/confirmatory_architecture_review.md`,
   section 8.
5. **The Open Targets statements**: their evidence grade, and whether held-out statements enter the primary score
   (`docs/label_source_review.md`, section 9.3).
6. **Development runs on v4.** The registered head choice and the stopping rule were made on `evidence_full_v2`. A new
   table needs its development runs before a held-out run, and the owner's instructions "No new pilots yet" and "don't
   rerun the same pilots until we've completely nailed down the design options" stand until the points above are
   settled.
7. **The owner's go-ahead for the held-out set** ("Don't finalize the lockbox until I've given the go ahead").

## 5. What the amendment should state plainly

- The label changes were chosen from properties of the evidence (the number of curated targets, the type of a gene's
  association, the frequency in a placebo arm) and never from a model's predictions. No model was trained on v4
  before the changes were fixed.
- The membership of the held-out set was read to place the added drugs; no held-out prediction or score was.
- The 29 held-out added drugs share targets outside the held-out set with 43 development drugs. That is why they are
  their own stratum and do not enter the primary score.
- Catatonia's one kept pair is in the held-out set, so catatonia has no training positive.
