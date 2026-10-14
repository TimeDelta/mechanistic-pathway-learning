# Decision 3.2: the label sentences that decide item 14, and a proposed reading of each

Item 14(a) admits a drug when **its own label states a central action** in a pharmacology section, "in a sentence
about that drug", and the item says the sentence "is then **read** and not counted, because class boilerplate ...
appears on every beta blocker's label and says nothing about the drug in hand, and a combination-product label can
put a drug's only match on a sentence about the other ingredient". That reading is decision 3.2 and it has not been
done. This file is the material for it.

Every matching sentence is in `docs/label_central_action_sentences.csv` (238 sentences, 16 drugs), written by
`experiments/read_label_central_action.py --every-sentence`. **The readings below are proposals, not decisions.**

## Two corrections to what the decision log said

**It is twelve drugs, not "five or six".** Twelve drugs outside the nervous-system code pass parts (b) and (c) of
the rule and are admitted or refused by part (a) alone: ondansetron, hydrocodone, cyclobenzaprine, propranolol,
nimodipine, granisetron, pindolol, felodipine, indomethacin, metoprolol, omeprazole and metformin. Four more have a
candidate sentence and fail (b) or (c), so a reading changes nothing for them unless their missing part is supplied:
bisoprolol, nebivolol, palonosetron and timolol.

**The per-drug screen table could not support the reading at all.** `configs/label_central_action.csv` stores one
quote per drug, and that quote is the **first** matching sentence in document order, not the most telling one
(`best = (negative or positive or adverse or [{}])[0]`). Ondansetron has 32 matches and pindolol 11, so for those
the stored quote was chosen by position. This is why the reading looked smaller than it is: judged on the stored
quote, several drugs look clear-cut in the wrong direction.

## Proposed: a sentence of the drug's own, so part (a) holds

| drug | the sentence | what it says |
| --- | --- | --- |
| hydrocodone | "Effects on the Central Nervous System The principal therapeutic action of hydrocodone is analgesia." | a central action, stated of the drug |
| cyclobenzaprine | "Such studies show that cyclobenzaprine acts primarily within the central nervous system at brain stem as opposed to spinal cord levels" | a central action, stated of the drug |
| nimodipine | "nimodipine had a greater effect on cerebral arteries than on arteries elsewhere in the body perhaps because it is highly lipophilic, allowing it to cross the blood-brain barrier" | a cerebral effect and penetration |
| metoprolol | "Metoprolol is also known to cross the blood brain barrier following oral administration and CSF concentrations close to that observed in plasma have been reported." | penetration, measured |
| propranolol | "Propranolol crosses the blood-brain barrier and the placenta, and is distributed into breast milk." | penetration |
| felodipine | "Animal studies have demonstrated that felodipine crosses the blood-brain barrier and the placenta." | penetration, in animals |
| indomethacin | "Indomethacin has been found to cross the blood-brain barrier and the placenta, and appears in breast milk." | penetration |

**One thing to settle with this group, because the rule's wording and its use have drifted apart.** Item 14(a) says
"states a central **action**". Four of the seven above state **penetration** instead - the drug crosses the
blood-brain barrier - which is not the same claim, and for felodipine it is penetration in animals. The two drugs
item 14 already records as admitted, metoprolol and propranolol, both rest on penetration sentences, so the
registered precedent reads (a) as "a central action **or** stated CNS penetration". Either that broader reading is
made explicit in the item, or metoprolol and propranolol have to be re-examined with the rest. It should not stay
implicit, because part (c) already measures penetration and (a) would then be a weaker copy of it rather than an
independent part.

## Proposed: no sentence of the drug's own, so part (a) fails

| drug | its best candidate | why it is not the drug's sentence |
| --- | --- | --- |
| pindolol | "The mechanism of the antihypertensive effects of beta-blocking **agents** has not been established, but several mechanisms have been postulated: 1) an effect on the central nervous system" | the class boilerplate item 14 names, and a statement that the mechanism is *not* established. None of its 11 matches names pindolol |
| bisoprolol | "Factors which may be involved include: ... 3) Diminution of tonic sympathetic outflow from the vasomotor centers in the brain." | the same boilerplate. None of its 9 matches names bisoprolol |
| nebivolol | "Possible factors that may be involved include: ... (3) diminution of tonic sympathetic outflow to the periphery from cerebral vasomotor centers" | the same boilerplate |
| granisetron | "Serotonin receptors of the 5-HT 3 type are located peripherally on vagal nerve terminals and centrally in the chemoreceptor trigger zone of the area postrema." | a statement about where the *receptor* is, not about the drug. None of its 12 matches names granisetron |
| palonosetron | "5-HT 3 receptors are located ..." | the same receptor-location statement. None of its 32 matches names palonosetron |

## Proposed: the sentence is negative, hedged, or about a different drug

These three are why the reading matters most, and all three would be admitted on their stored quote.

| drug | its only or leading sentence | the problem |
| --- | --- | --- |
| **metformin** | "In mice, no **glipizide** or metabolites were detectable autoradiographically in the brain or spinal cord" | its **one and only** matching sentence is about glipizide, the other ingredient of a combination-product label. It is also negative. This is exactly the failure item 14 warns of |
| **omeprazole** | "Other Effects Systemic effects of omeprazole in the CNS, cardiovascular and respiratory systems have **not been found** to date." | a negative statement. The screen's negative detector does not match "have not been found to date", so all 6 of its sentences were classed as positive |
| **ondansetron** | "It is **not certain** whether ondansetron's antiemetic action is mediated centrally, peripherally, or in both sites." | the label declines to say. 32 matches, none of which asserts a central action of the drug |

**This is not a defect in the screen.** Item 14 makes the screen a candidate-finder and the reading the authority,
for these cases precisely. What it does show is that the screen's verdict column cannot be used as the rule: on
`configs/label_central_action.csv` alone, metformin, omeprazole and ondansetron all read "states a central action".

## Timolol, which is a genuine judgement and not a category

"Possible mechanisms of action include reduction in cardiac output, reduction in plasma renin activity, and a
central nervous system sympatholytic action." It is in timolol's own label and its subject is timolol's mechanism,
so the "does it name the drug" test used above gets it wrong; but it is hedged ("Possible mechanisms") and it is the
same three-mechanism list the other beta blockers carry. Timolol fails part (c) anyway for want of a measured
Kp,uu,brain (decision 3.1), so nothing turns on it until that is supplied.

## What the reading changes

If the proposals hold, part (a) removes **pindolol, granisetron, metformin, omeprazole and ondansetron** from the
drugs the rule would otherwise admit, and leaves **hydrocodone, cyclobenzaprine, nimodipine, metoprolol,
propranolol, felodipine and indomethacin**. Five of twelve, which is why this reading is not a formality. None of
it is decided here.


## Codeine, and a test of whether the screen's filter misses receptor mechanisms

The owner, 13 October 2026: "Codeine also hits opioid receptors in the brain ... is that recorded?"

**Part (b) records it fully.** Codeine's mechanism targets are `OPRM1`, `OPRD1` and `OPRK1`, and the graph's GTEx
annotation marks **all three brain-expressed**. Its Kp,uu,brain is 1.00, so part (c) holds as well.

**The label states the mechanism too, and the screen cannot see that sentence.** From the Mechanism of Action
section of codeine's own label: *"Codeine is an opioid agonist relatively selective for the mu-opioid receptor, but
with a much weaker affinity than morphine."* The screen keeps only sentences matching `CENTRAL_ACTION`, which
requires a word naming central anatomy, and that sentence names none. So codeine passes part (a) on a different
sentence: *"Codeine produces respiratory depression by direct action on brain stem respiratory centers"* (12.2
Pharmacodynamics).

**That looked like a systematic blind spot and it was measured rather than assumed.** Of the 24 drugs the screen
records as "no such sentence", **20 have at least one sentence naming the drug and a target with no central-anatomy
word**. But reading them settles it the other way: they are statements of *what the target is*, not of acting
centrally - "Verapamil hydrochloride is a calcium antagonist or slow-channel inhibitor", "Betaxolol is a
beta-1-selective (cardioselective) adrenergic receptor blocking agent", "Celecoxib is a potent inhibitor of
prostaglandin synthesis in vitro". A calcium channel or a cyclo-oxygenase is not a central claim. **Whether that
target is in the brain is exactly what part (b) tests**, from the graph's structured annotation rather than from
sentence matching, and it tests it better.

So the filter is right for part (a)'s job: requiring central anatomy is what separates a statement of central
action from a statement of target identity, and the second belongs to (b). For codeine nothing is lost - the
receptor mechanism is in (b) with all three receptors brain-expressed, and (a) contributes the thing (b) cannot, a
stated central **effect**.

**One check on that sentence, since respiratory depression is an unwanted effect.** Item 14 refuses to read the
adverse-reaction, warning and indication sections, because admitting a drug for the psychiatric events on its label
is the circularity the rule exists to avoid. This sentence is in a pharmacology section, it states a mechanism
rather than reporting an incidence, and respiratory depression is not among the study's 24 symptoms, so no
circularity arises through it.
