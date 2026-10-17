# Edge-gate simulations with planted truth (decision 9.62), 16 October 2026

Recorded 16 October 2026 (Eastern time; the runs carry UTC dates of 17 October). A record: corrected only by dated
additions. The numbers are those of `docs/edge_gate_simulation_results.json`, which
`experiments/summarise_edge_gate_simulations.py` wrote from the runs of `scripts/run_edge_gate_simulations.sh`;
figures 15 and 16 of `docs/figures.md` draw them.

## 1. Why a simulation

On the real graph no one knows which edges are false, so how well a fit of the hard edge gates of decision 9.62
finds them can be measured only on graphs whose truth is planted. A fit to the real labels is a pilot, which the
registration allows on the development set only ("Pilots of any configuration run on the development set only",
`docs/preregistration.md`, order of work), and it measures fit, not recovery. The settings the trainer uses are taken
from these runs.

## 2. The worlds and the fit

`experiments/simulate_edge_gate_recovery.py` builds each world: 300 nodes with lognormal in-degrees (median 6, so
2,086 to 2,465 edges by seed), two relations (an edge inhibits with probability 0.3 and activates otherwise), and the
v5 graph's evidence tiers in its proportions: curated (p 0.9, n 1, so π = 0.900; 68 percent of edges), unreferenced
(p 0.7, n 0.5, π = 0.604; 20 percent) and text-mined only (p 0.55, n 0.5, π = 0.525; 12 percent), with π =
sigmoid(n logit p). Each edge is real with probability π, so the priors are calibrated by construction, the most
favourable case for them. The labels come from the real edges only: a linear-response encoder on the true graph, 800
perturbations each seeding one node with sign +1 or -1, and 12 symptoms each reading 6 nodes through a logistic link,
the readout scaled so the logits have a standard deviation of 2 before a bias of -1 (a positive rate of 0.27 to 0.29),
drawn as Bernoulli. A quarter of the perturbations are held out.

The fit is the same encoder on every edge with the truth's parameters (the oracle case, which isolates what the gates
can learn), hard concrete gates (γ -0.1, ζ 1.1, β 2/3; Louizos, Welling and Kingma 2018) started at P(open) = π, 20
stochastic epochs (batches of 32, Adam at 0.05 on the gates' logits), the binary pass (8 full-batch epochs, Adam at
0.1, surrogate width 0.5), 30 exact flip checks, and for each closed edge its exact evidence E: the rise in the
training data loss, times W, when that edge alone is reopened, in nats, where W is the number of training labels.
Stability selection refits on 20 halves of the perturbations, 10 complementary pairs (Shah and Samworth 2013).

## 3. Results

### 3.1 The read-out

Over the five calibrated worlds the read-out closed 234 edges, 46.8 per world of about 2,300; 64.1 percent of them
were false (95 percent Wilson interval 0.58 to 0.70) against 21.1 percent of all edges, and they were 6.2 percent of
the false edges. By tier: 11 curated closures (8 false; the tier's prior puts 0.10 of its edges false), 58 unreferenced
(47 false, 0.81, against 0.40) and 165 text-mined (95 false, 0.58, against 0.47). A text-mined edge closes when
closing it lowers the training data loss, times W, by more than logit(0.525) = 0.10 nats, a bar that noise clears
often, so most of the read-out's closures, and most of its errors, are in that tier.

### 3.2 Which closures to believe

Pooled over the five worlds (figure 15, left):

| rule | closures per world | false (precision) | 95 percent interval |
| --- | --- | --- | --- |
| the read-out, every closed edge | 46.8 | 0.641 | 0.58 to 0.70 |
| exact evidence E ≥ 0.5 nats | 24.6 | 0.813 | 0.74 to 0.87 |
| exact evidence E ≥ 1 nat | 13.4 | 0.910 | 0.82 to 0.96 |
| exact evidence E ≥ 2 nats | 7.6 | 0.895 | 0.76 to 0.96 |
| P(false) = sigmoid(E - logit π) ≥ 0.7 | 11.6 | 0.931 | 0.84 to 0.97 |
| P(false) ≥ 0.8 | 7.0 | 0.943 | 0.81 to 0.98 |
| P(false) ≥ 0.9 | 4.2 | 1.000 | 0.85 to 1.00 |
| P(false) ≥ 0.95 | 3.2 | 1.000 | 0.81 to 1.00 |
| closing frequency over the 20 refits ≥ 0.6 | 19.6 | 0.796 | 0.71 to 0.86 |
| closing frequency ≥ 0.7 | 13.6 | 0.853 | 0.75 to 0.92 |
| closing frequency ≥ 0.8 | 9.0 | 0.889 | 0.77 to 0.95 |
| closing frequency ≥ 0.9 | 4.6 | 0.870 | 0.68 to 0.95 |

At each similar number of closures the posterior from the exact evidence was more precise than the closing
frequency, though the intervals overlap, and it needs one fit where stability selection needs twenty-one.

### 3.3 Is the posterior calibrated?

At a local optimum under single flips every closed edge has E ≥ logit π, so P(false) ≥ 1/2; one closure of the 234
fell below, since the checks test only some flips. Binning the other 233 by P(false) (figure 15, centre):

| P(false) | closures | mean P(false) | share false | 95 percent interval |
| --- | --- | --- | --- | --- |
| 0.5 to 0.6 | 134 | 0.536 | 0.485 | 0.40 to 0.57 |
| 0.6 to 0.7 | 41 | 0.644 | 0.756 | 0.61 to 0.86 |
| 0.7 to 0.8 | 23 | 0.743 | 0.913 | 0.73 to 0.98 |
| 0.8 to 0.9 | 14 | 0.859 | 0.857 | 0.60 to 0.96 |
| 0.9 to 1 | 21 | 0.967 | 1.000 | 0.85 to 1.00 |

Each bin's mean lies within its interval or below it: the posterior was calibrated or conservative. This is the most
favourable case: the priors are calibrated by construction and the encoder has the truth's parameters, and neither
will hold on real data. E also holds every other gate as fitted, so P(false) is conditional on the rest of the fit.

### 3.4 The first-order flip estimate

On 300 random edges per world, 1,500 in all, the binary pass's first-order estimate of the data part of a flip's
effect correlated 0.75 with the exact change and agreed with it in sign 96 percent of the time; among the 309 edges
whose exact change was at least 0.1 nats, 0.77 and 92 percent (figure 15, right). The exact checks accepted 36 of the
122 flips they tested. Among the tested flips the estimate's correlation with the exact change ran from -0.79 to 0.46
by world. The checks test the flips with the largest estimated falls, and those are where the estimate's errors are
largest, a selection effect; it is why the checks are exact.

### 3.5 What each pass contributes

Without the binary pass and the checks (the stochastic pass's read-out alone) the five worlds closed 295 edges at
precision 0.637 (188 false); with them, 234 at 0.641 (150 false). By tier, without the pass and with it: curated 4
(3 false) and 11 (8 false), unreferenced 77 (56 false) and 58 (47 false), text-mined 214 (129 false) and 165 (95
false). The held-out loss was lower with the binary pass in
three of the five worlds and higher in two (differences -0.0031 to +0.0007 nats per held-out label, mean -0.0007), so
whether the pass improves prediction is not settled by five worlds. The pass did not settle within its 8 epochs: in
the last epoch it still flipped 2 to 5 edges per world, so the read-out depends on when it stops.

In world 0, with the code of batch 2: the surrogate widths 0.25, 0.5 and 1 gave the same read-out (53 closures, 31
false), and without the 30 flip checks the read-out was 53 closures with 32 false and a held-out loss of 0.5495,
against 0.5507 with them. One world does not show whether the checks help.

### 3.6 When the favourable assumptions are dropped (world 0 unless stated)

The held-out gain share is (loss with every edge open - loss of the fitted graph) / (loss with every edge open - loss
of the true graph), on the held-out perturbations. It is noisy: the attainable gain was 0.0013 to 0.0060 nats per
held-out label, and over the five calibrated worlds the share ran from -1.20 (world 0) to +0.84 (mean +0.17). A
world-0 condition is therefore read against world 0, and precision carries more weight than the share.

| condition | closures | precision (false share of all edges) | closures with E ≥ 1 nat (precision) | held-out gain share |
| --- | --- | --- | --- | --- |
| calibrated, five worlds pooled | 234 | 0.641 (0.211) | 67 (0.910) | +0.17 (mean) |
| calibrated, world 0 | 53 | 0.585 (0.211) | 14 (0.929) | -1.20 |
| 3,000 perturbations | 92 | 0.663 (0.211) | 45 (0.800) | +0.87 |
| priors too optimistic (curated 0.8, unreferenced 0.45, text-mined 0.3 real) | 63 | 0.857 (0.339) | 20 (0.950) | -0.24 |
| priors too pessimistic (0.97, 0.85, 0.75 real) | 47 | 0.298 (0.084) | 12 (0.750) | -0.24 |
| half the true positives observed | 46 | 0.500 (0.211) | 15 (0.667) | none* |
| joint fit, fresh start, readout on each symptom's true nodes | 28 | 0.714 (0.211) | 1 (1.000) | -0.17 |
| joint fit, fresh start, readout on every node | 9 | 0.444 (0.211) | 0 | -0.22 |
| joint fit from the truth's values, readout on every node | 35 | 0.686 (0.211) | 5 (1.000) | -1.12 |
| v5's positives per gated edge (2,000 nodes, 20 symptoms at a 5 percent positive rate) | 28 of 15,012 | 0.536 (0.207) | 5 (0.800) | -0.42 |

\* With half the positives observed, the true graph predicted the observed held-out labels worse than every edge
open (0.4783 against 0.4753), so the share is undefined; the fitted graph predicted them slightly better than every
edge open (0.4748).

Against world 0, more labelled perturbations raised precision a little (0.585 to 0.663), doubled recall (6.5 to 12.8
percent of the false edges) and gave the one clearly positive held-out share. Priors that overstate how often edges
are false (the pessimistic case) made most closures wrong: precision 0.30 against 0.08 of edges false, so the
closures were still enriched for false edges, but 33 of 47 closed a real edge; the closures with at least 1 nat of
evidence held up better (0.75). Under-observed positives lowered precision (0.50) and the precision of the closures
with strong evidence (0.67). The joint fits, which learn the encoder's parameters and the readout with the gates
(AdamW at the trainer's weight decay of 1e-4, early stopping on an inner validation share at best epochs of 2 to 11),
closed 9 to 35 edges, and in all three the fitted graph predicted the held-out labels worse than every edge open. A
fresh start also left the encoder short of the truth: with the true graph, its held-out loss was 0.578 with the
readout on every node and 0.554 with it on the true nodes, against 0.548 for the truth's parameters. The first joint
run, before early stopping and weight decay were added, overfitted to a held-out loss of 1.11 and is not used.

### 3.7 Curated edges taken as known (decision 9.63)

With the curated tier ungated (`--ungated-tiers curated`), the five worlds closed 225 edges at precision 0.640, against
234 at 0.641 with it gated. The other tiers barely moved: the same 165 text-mined closures with the same 95 false,
and 60 unreferenced closures (49 false) against 58 (47 false). The held-out loss was lower with the curated tier gated
in three worlds and higher in two (mean difference -0.0003 nats per held-out label). Taking curated edges as known
removed their 11 closures, 3 of which had closed real edges, and changed nothing measurable about the rest. In these
worlds curated edges are real with probability 0.9 by construction; the question 9.63 answers on real data, whether
a curated edge with two independent records including a human study is that reliable, is not tested here.

### 3.8 Stability selection's bounds are not error rates here

At a closing frequency of 0.7 the bound of Meinshausen and Bühlmann (2010, Theorem 1), E(V) ≤ q² / ((2π_thr - 1) p),
was exceeded in four of the five worlds: the real edges closed were 4, 2, 3 and 1 against bounds of 1.54, 1.69, 2.06
and 0.36 (in world 2, none against 1.01). The theorem assumes that "the distribution of {1{k∈Ŝλ}, k ∈ N} is
exchangeable for all λ ∈ Λ" and that the procedure is "not worse than random guessing"; the authors write that "For
real data, we have no guarantee that the assumption is fulfilled but the numerical examples in Section 4 show that the
bound holds up very well for real data" (arXiv 0809.2932, read 17 October 2026). When the selection is a closure, the
noise variables are the real edges, and the cost of a closure, logit π, differs by tier, so their closing
probabilities differ by tier and they are not exchangeable; here the bound did not hold. The complementary-pairs
bounds of Shah and Samworth (2013) "require no (e.g. exchangeability) assumptions on the underlying model" (arXiv
1105.5578, abstract, read 17 October 2026), but they bound how many variables are selected from the set that "have low
selection probability" under the procedure fitted to half the data, not how many real edges are closed: a real edge
that the half-data fit closes often is outside that set. Neither bound limits the number of real edges closed; the
measured precision above is the evidence.

## 4. What this says for a fit on v5

The v5 label table has 8,199 label-grade positive (perturbation, symptom) pairs against 163,276 gated edges, 0.050 per
gated edge counting every perturbation, an upper bound for what one training fold sees. The 800-perturbation worlds
had 0.80 to 0.99 training positives per edge and the v5-ratio world 0.043. In that world the read-out closed 28 of
15,012 edges, 0.2 percent, at a precision of 0.54 against 0.21 of edges false, and the fitted graph predicted the
held-out labels worse than every edge open (-0.42); the 5 closures with at least 1 nat of evidence were 4 false edges
and 1 real one. Expected on v5: the read-out keeps nearly every edge, and the closures worth reading are the few with
strong exact evidence.

## 5. Settings kept and what is put to the owner

The trainer (`experiments/run_main_model.py`) takes the simulated settings: surrogate width 0.5 (no effect found
across 0.25 to 1), learning rate 0.1 for the binary pass, full-batch steps, and, from this record on, 8 binary-pass
epochs (it had 5). Flip checks stay off by default (`--edge-gate-flip-checks 0`): they cost one pass over the fit rows
each and showed no gain in world 0. The closure evidence of decision 9.63 (c) exists as
`--edge-gate-closure-evidence N`, which writes E and P(false) for up to N closed edges, strongest priors first, to
`edge_gate_closures.csv`; listing every closure by default, with its full evidence, is part of the 9.63 build. Put to
the owner, not decided here: reading closures by P(false) or E rather than by the read-out alone, and reporting
closing frequencies as descriptive rather than as an error rate.

## 6. Limits

One encoder (linear response); random graphs rather than the v5 topology; one seed node per perturbation; 12
symptoms; priors calibrated by construction except where stated; labels drawn from the true model (thinned in one
run); the truth's parameters except in the three joint fits; one world per condition except the five calibrated
worlds and their variants; a binary pass that had not settled; a held-out gain share too noisy to compare single
worlds; and the joint fits' early stopping also ended their gates' stochastic pass early.

## 7. Reproduction

`bash scripts/run_edge_gate_simulations.sh` writes every run to `runs/simulations/` and then the digest. The runs named
`b3_` started at commit 93b97fe; the simulation script changed while they ran (684b6c2 reads closures per edge, which
is the same for a run that gates every edge; cf57d4b moved the closure evidence into `models/edge_gate_passes.py`,
with the same output). The `b3_s*_curated_fixed` runs were made at commits 684b6c2, which added `--ungated-tiers`,
and cf57d4b. The runs named `b2_` are from commit f6d7b0d; `b2_s0_width05` gave the same read-out and held-out losses
as `b3_s0`.

## Addition, 17 October 2026: batch 4, whether the binary pass settles

Section 6 lists "a binary pass that had not settled" among the limits: at 8 epochs and a constant learning rate the
pass still flipped 2 to 5 gates in its last epoch (1 to 8 over its last four), so the closed set depended on the epoch
it stopped at. Batch 4
repeated the five calibrated worlds of section 3 with 32 binary epochs, once at the constant rate of 0.1 and once with
the rate multiplied by 0.8 after each epoch (`--binary-learning-rate-decay 0.8`), everything else as in `b3_s*`.

| world | flips in the last four epochs: 8, constant (`b3`) | 32, constant | 32, x 0.8 | precision after the checks: 8, constant | 32, constant | 32, x 0.8 | held-out loss of the fitted graph: 8, constant | 32, x 0.8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 3, 1, 8, 5 | 2, 5, 6, 3 | 0, 0, 0, 0 | 0.585 | 0.596 | 0.604 | 0.55071 | 0.55073 |
| 1 | 2, 2, 6, 3 | 4, 4, 3, 5 | 0, 0, 0, 0 | 0.694 | 0.667 | 0.686 | 0.56256 | 0.56316 |
| 2 | 3, 4, 3, 3 | 3, 5, 4, 8 | 0, 0, 0, 0 | 0.628 | 0.628 | 0.628 | 0.56170 | 0.56170 |
| 3 | 5, 4, 6, 3 | 7, 8, 7, 9 | 0, 0, 0, 0 | 0.689 | 0.683 | 0.698 | 0.56796 | 0.56788 |
| 4 | 4, 4, 2, 2 | 6, 3, 3, 4 | 1, 0, 0, 0 | 0.571 | 0.607 | 0.571 | 0.54311 | 0.54311 |

What it shows.

- **More epochs alone do not settle the pass; the annealed rate does.** At a constant rate the 32-epoch pass flipped 2
  to 9 gates in each of its last four epochs in every world; with the rate times 0.8 an epoch it flipped none in the
  last three epochs of every world. That is what a step-size effect would give: gates near the threshold cross it back
  and forth while the step stays the same size.
- **Settling did not change what is recovered.** Pooled over the five worlds the closures' precision was 0.641 at 8
  epochs (150 of 234 closed edges false), 0.641 at 32 constant (152 of 237) and 0.647 annealed (154 of 238); the
  closures with at least 1 nat of exact evidence 0.91, 0.90 and 0.91; and the held-out loss of the fitted graph moved
  by at most 0.0006 nats in any world (world 1, annealed against 8 epochs). Figure 16 has both conditions as rows.
- **A settled pass is not a local optimum under single flips.** After the annealed pass the exact flip checks still
  accepted 2 to 6 of 16 to 30 candidate flips per world (6 of 27 in world 0). The pass follows the first-order
  estimate, whose correlation with the exact change is weak among the largest promised changes (section 3), so where
  it stops is where the estimate stops promising improvement, not where the objective does.

For the trainer. `run_main_model.py` takes the option as `--edge-gate-binary-learning-rate-decay`, default 1 (the
rate it has always kept); its help text gives these results. Put to the owner (decision 9.62): 0.8 over 32 epochs
as the default, which makes the closed set independent of the epoch the pass stops at at no measured cost to
recovery, against the present 8 epochs at a constant rate. The pass's cost scales with its epochs: on v5 each binary
epoch costs about a training epoch with the limiting-step cache rebuilt (`docs/data_layer_rebuild.md`).

Reproduction. The `b4_` runs were made on 17 October 2026 between 04:05 and 04:30 UTC at commit a1af139, which added
the decay option, and 2e5d7de, which added a function to `models/edge_gate_passes.py` that the simulation does not
call. `bash scripts/run_edge_gate_simulations.sh --list` gives their commands.
