# Plasma substrate binders

The user's decision of 9 October 2026: "All of the plasma substrate binders should be modeled. Once AAG and albumin
are added, there will be state to attach that separate extension to." This document is the scope, the audit of what
this repository holds and the build in two stages. Stage 1, the carriage edges, is built (see Built below) in new graph
directories; stage 2, the fraction unbound, still has no source. Nothing in the confirmatory specification depends on
either.

It answers an earlier request as well: a measured occupancy that carries albumin, "as that is an important
differentiator between i.e. intermittent substance use and chronic substance use disorder".

## Why the binders are in the model at all

A drug in this study is not a node. It is a perturbation that seeds its target nodes with a magnitude, and the
occupancy rule of docs/off_target_scoping.md turns affinities into those magnitudes from the free concentration at the
target. Plasma binding sets that free concentration, so it enters in two places, and they are not the same work:

1. **The binders as proteins that carry things.** Albumin, orosomucoid and the rest are in the graph as gene nodes with
   no edge to anything they carry, so their physiological job is absent: albumin's carriage of fatty acids and
   bilirubin, transthyretin's and thyroxine-binding globulin's carriage of thyroid hormone, the steroid-binding
   globulins' carriage of cortisol and the sex steroids. These are `binds` edges, and `binds` already carries sign 0 in
   every one of the full graph's 31,118 such edges, which is the right reading for a carrier: it moves a substrate
   without acting on it.
2. **The free fraction of a drug.** One fraction unbound per drug is a constant that rescales that drug's seed
   magnitudes, which is the occupancy extension. The state the user wants to attach to it, raised orosomucoid or low
   albumin, needs a second condition per drug rather than a descriptor, because a constant cannot express a difference
   between two patients.

## Which binders

| binder | genes | what it carries |
|---|---|---|
| albumin | ALB | long-chain fatty acids, unconjugated bilirubin, bile acids, calcium, magnesium, tryptophan, thyroxine, acidic drugs |
| orosomucoid (alpha-1-acid glycoprotein, AAG) | ORM1, ORM2 | basic and neutral lipophilic drugs, steroid hormones |
| transthyretin | TTR | thyroxine, triiodothyronine, retinol through retinol-binding protein | 
| thyroxine-binding globulin | SERPINA7 | thyroxine, triiodothyronine |
| corticosteroid-binding globulin | SERPINA6 | cortisol, progesterone |
| sex hormone-binding globulin | SHBG | testosterone, dihydrotestosterone, estradiol |
| retinol-binding protein | RBP4 | retinol, carried as a complex with transthyretin |
| lipoproteins | APOA1, APOB and the other apolipoproteins | cholesterol, triacylglycerols, lipophilic drugs |
| alpha-fetoprotein | AFP | fetal carriage of fatty acids and bilirubin |

Retinol-binding protein is on the list although the original scope did not name it: it is the specific retinol carrier and
transthyretin's partner, so leaving it out would have given transthyretin a cargo with no carrier beside it.

What the graph holds now: `GENE:ALB`, `GENE:ORM1`, `GENE:ORM2`, `GENE:SERPINA6`, `GENE:SERPINA7`, `GENE:TTR`,
`GENE:AFP`, `GENE:SHBG` and the apolipoprotein genes exist, and `GENE:ALB` carries 49 edges: 32
`regulates_transcription_of`, 10 `binds` to other proteins (LRP2, FCGRT, VCAM1, B2M and others), 3 `activates`, 2
`member_of` and 2 `catalyzed_by` (its synthesis MAR05151 and its degradation MAR05258). Not one of its edges reaches a
small molecule it carries.

The substrates are present as metabolite nodes with an extracellular compartment, which is where a plasma binder meets
them: bilirubin (9 nodes over compartments c, e, r), thyroxine (3), cortisol and its precursors (9), testosterone (27),
estradiol (42), progesterone (8), retinol (18), palmitate and the other long-chain fatty acids, tryptophan and
cholesterol. So stage 1 needs no new nodes, only edges.

## Why orosomucoid, not albumin, for this study's drugs

This is the correction to make before albumin is used as the differentiator. According to PubMed, orosomucoid is the
carrier of the basic drugs, which is what most psychotropics are: Wallace and Verbeeck report that "Most drugs which
exhibit increased binding (decreased free fraction) in elderly subjects are basic and tend to have a greater affinity
for alpha 1-acid glycoprotein than for albumin"
([DOI](https://doi.org/10.2165/00003088-198712010-00004), PMID 3545616), and Fournier and colleagues that orosomucoid
"has the ability to bind and to carry numerous basic and neutral lipophilic drugs from endogenous (steroid hormones)
and exogenous origin; one to seven binding sites have been described. AGP can also bind acidic drugs such as
phenobarbital" ([DOI](https://doi.org/10.1016/s0167-4838(00)00153-9), PMID 11058758). Albumin is the carrier of the
acidic drugs, valproate and phenytoin among the ones this study holds.

The same review is the reason orosomucoid, not albumin, carries the state the user is after: "AGP is one of the major
acute phase proteins in humans, rats, mice and other species. As most acute phase proteins, its serum concentration
increases in response to systemic tissue injury, inflammation or infection, and these changes in serum protein
concentrations have been correlated with increases in hepatic synthesis"
([DOI](https://doi.org/10.1016/s0167-4838(00)00153-9), PMID 11058758). A raised orosomucoid lowers the free fraction of
a basic drug; a low albumin, from liver disease or malnutrition, raises the free fraction of an acidic one. Chronic use
moves both, in opposite directions and for different drugs, which is why the extension needs both binders and not one.

Two further cautions for the extension, not for stage 1. "One to seven binding sites" means the one-site Langmuir form
of the occupancy rule is an approximation on the binder even where it is reasonable at the target. And a fraction
unbound measured in normal plasma already carries the normal binder concentrations, so a model that scales by it and
then also varies the binder state would count the same binding twice unless the fraction is re-derived from the
binder's affinity and concentration.

## Stage 1: the carriage edges

Built: `binds` edges, sign 0, from the extracellular copy of each cargo metabolite to its binder's node, in new graph
directories rather than over `graph_full_neuronal` or its split, so no preregistered input changes. Direction and sign
follow the graph's own convention for a small molecule that binds a protein: its 436 `small_molecule_protein` edges
from OmniPath run metabolite -> gene, and all 31,118 `binds` edges carry sign 0.

Source. No pinned file in this repository carried ligand annotations, which is why one fetch was needed:
`data/raw/uniprot/uniprot_human_reviewed.tsv.gz` holds eight columns (Entry, primary gene name, HGNC, EC number,
subcellular location, Pfam, length, sequence) and no binding site, ligand or function text, and the ChEMBL tool's ADMET
endpoint returns calculated properties only. experiments/fetch_plasma_binder_annotations.py now pins
`data/raw/uniprot/uniprot_plasma_binders.tsv.gz` (UniProt release 2026_03, 2,989 bytes, sha256 `eab6550251dd8056`),
the reviewed entries of the eleven binder genes with `ft_binding`, `cc_function` and `cc_subcellular_location`. Ligand
entries carry ChEBI identifiers, which the same route the Reactome import uses turns into Human-GEM metabolites
(`map_definition_to_metabolites`).

What that fetch does not give, which is the reason for the second source. The binding-site features are structurally
evidenced and therefore sparse: they give albumin bilirubin IXalpha, Ca(2+), Cu and Zn(2+), corticosteroid-binding
globulin cortisol, thyroxine-binding globulin and transthyretin thyroxine, alpha-fetoprotein Cu(2+), and **nothing at
all** for orosomucoid, sex hormone-binding globulin, retinol-binding protein or the apolipoproteins. Albumin's fatty
acids, the protein's best-known cargo, are in its function text but not in its features. So the carriage those features
miss is in `docs/curated_plasma_carriage.csv`, one row per (binder, cargo) with its source, PMID, DOI and a quoted
sentence, in the form the repository already uses for the presynaptic receptor and oxidant tables.

## Known limitation: no drug reaches its plasma carrier in the graph

**None of the drugs orosomucoid carries has a carriage edge, and the one that matters most is fentanyl.**

ORM1 and ORM2 have exactly one carriage edge each, to progesterone. Every one of the 27 carriage edges starts at a
metabolite node; **0 of the study's 142 drug perturbations is a node in the graph at all**, because a drug enters as
a sustained input at its protein targets rather than as a node of its own. So the physiology the carriage edges were
added for - a binder holding a ligand and lowering its free concentration - is represented for progesterone,
thyroxine, retinol, cortisol and the fatty acids, and is absent for every basic drug, which is the cargo that made
orosomucoid worth modelling in the first place.

**Fentanyl is the drug this costs.** Of the 23 kept positive pairs the eleven orosomucoid-carried drugs hold,
fentanyl holds 12, over abnormal_dreams, anxiety, cognitive_impairment, depressed_mood, elevated_mood_or_mania,
fatigue, increased_appetite, insomnia, irritability_or_aggression, psychomotor_agitation, psychosis and
somnolence_or_hypersomnia. It is the one orosomucoid-carried drug whose labels the model is really scored on, its
free fraction is the one most moved by an acute phase response, and the graph gives it no path to ORM1 or ORM2.
Lidocaine (4 kept positives) and reboxetine (5) are the next two; the other eight hold 0, 1 or 2 between them.

**Why the edge is not simply added.** Three things stand between the drug and the edge, and only the first is a
hard stop:

1. There is no node to attach it to. None of the eleven is in Human-GEM, because they are xenobiotics;
   `configs/drugs_acting_as_graph_compounds.csv` holds the one drug that is an endogenous metabolite, tryptophan.
   A carriage edge for fentanyl needs a fentanyl node, which the graph would have to gain.
2. A drug node with only a carriage edge would carry nothing. The perturbation still enters at the drug's protein
   targets, so a node joined to ORM1 and to nothing else is a dead end that no field reaches. Making the edge
   informative means seeding the drug node instead, which changes how every drug perturbation enters.
3. A `binds` edge does not sequester. The graph gives all 31,118 of them sign 0, which reads as "moves a substrate
   without acting on it". Representing that bound drug is unavailable needs a signed depletion edge from binder to
   cargo, in the way `depletes_substrate` was added because "the stoichiometric reading needs edges the graph does
   not store" (`linear_response_encoder.py`). That applies to the 27 existing carriage edges too: as they stand,
   albumin binding bilirubin does not lower free bilirubin in this model either.

So the limitation is not that the edges are unplaceable. It is that placing them needs drug nodes, a change to
where a drug perturbation is seeded and a signed sequestration relation, which is three changes to the confirmatory
graph and the perturbation model, not one. Recorded here as a known limitation of the graph rather than attempted.
docs/drug_entry_nodes.md sets out those three changes, what the record holds about an earlier intention to make them
and the two objections to weigh, at the user's question of 9 October 2026 about drug nodes as new graph entry points.

## Stage 2: the free fraction and the binder state

The user's request of 9 October 2026, "Can you please find the data you need for orosomucoid?", is answered here.
The data is found and pinned, it reaches a quarter of this study's drugs, and one source it rests on refutes the
plan this document gave for using it.

### What is pinned

`data/raw/plasma_protein_binding/fraction_unbound_database.xlsx`, fetched by
experiments/fetch_fraction_unbound_database.py: the supplementary database of Al-Qassabi and colleagues (figshare
`10.48420/25243138.v1`, CC BY 4.0; the paper is `10.1016/j.xphs.2024.02.024`, PMID 38417790). The fetch goes through
the figshare API so the file's md5 can be checked against the one figshare publishes for it, which it matches.

It holds 556 measurement rows over 209 drugs, each row carrying

- the **major binding protein**: albumin (347 rows), orosomucoid (122) or both (87), which is the carrier identity no
  other source here gives;
- the **fraction unbound in a reference population**;
- the **fraction unbound in a special population**: renal impairment, hepatic impairment, elderly, paediatric,
  inflammatory disease (Crohn's disease, rheumatoid arthritis) and ethnicity;
- a PMID or DOI for the measurement.

That is the state the user asked for, measured rather than modelled: inflammatory disease is the raised-orosomucoid
condition and hepatic impairment the low-albumin one, and the database gives both ends of each.

### What it reaches, measured

experiments/measure_fraction_unbound_coverage.py, reported in docs/fraction_unbound_coverage.md: **27 of this
study's 142 drug perturbations**, and 2 of the 28 in the lockbox (caffeine and diphenhydramine). Eleven of the
twenty-seven name orosomucoid as a major binding protein: alfentanil, desipramine, diphenhydramine, fentanyl,
lidocaine, meperidine, methadone, reboxetine, ropivacaine, sufentanil and tiagabine. Two rows show the effect the
extension is for: methadone's fraction unbound falls from 0.1062 to 0.0912 in rheumatoid arthritis, and alfentanil's
from 0.11 to 0.052. One hundred and fifteen of the study's drugs have no row, so for them stage 2 has no free
fraction from this source.

### The plan this document gave for stage 2 is wrong, and the correction

The cautions above proposed re-deriving the fraction unbound from the binder's affinity and concentration, so that a
model that also varies the binder state would not count the same binding twice. The paper behind the pinned database
tested that approach and it fails: "concordance correlation coefficients for predicted fold-change in fu for the same
dataset were <0.38 for all populations and sub-groups", and "the predictions of fu solely based on changes in protein
concentrations in plasma cannot explain the observed values in some special populations" (abstract, PMID 38417790;
[DOI](https://doi.org/10.1016/j.xphs.2024.02.024)). Their recommendation is the opposite of the plan: "PBPK models of
special populations for highly bound drugs should preferably use measured fu data". So the second condition per drug
is a **pair of measured fractions**, the reference one and the special-population one, not a reference fraction
rescaled by an orosomucoid concentration. The double-counting the plan feared does not arise, because neither number
is derived from the other; what is lost is any drug the database does not measure in that population, which cannot be
filled by prediction without reintroducing exactly the error the paper quantifies.

### Why orosomucoid rather than albumin carries this state, with the sources

- Orosomucoid can be the carrier that matters even though it is the rarer protein: "hundreds of drugs with diverse
  structures bind to this glycoprotein. Although plasma concentration of AAG is much lower than that of albumin, AAG
  can become the major drug binding macromolecule in plasma with significant clinical implications" (Israili and
  Dayton, abstract, PMID 11495502, [DOI](https://doi.org/10.1081/dmr-100104402)).
- The free fraction of a basic drug tracks the orosomucoid concentration over the whole physiological range: "As the
  alpha 1-AGP increased from 0.05 to 2.0 gm/l, free fraction fell from 92.40% to 8.80%", measured for methadone, whose
  free fraction over 29 healthy subjects was "10.62 +/- 1.43" per cent (Romach and colleagues, abstract, PMID 7193106,
  [DOI](https://doi.org/10.1038/clpt.1981.34)). The same paper bounds where it does not matter: "Less than 20% of
  naloxone, codeine, morphine, heroin, pentazocine, and diphenoxylate bound to alpha 1-AGP."
- The per-drug binding constants are in Kremer, Wilting and Janssen, Pharmacol Rev 40:1-47, 1988 (PMID 3064105,
  [DOI](https://doi.org/10.1016/S0031-6997(25)00015-8)), which PubMed indexes without an abstract, so nothing from it
  is quoted or used here; it is recorded as the place to go if a drug the database misses has to be derived.

These are PubMed records, and PubMed asks that it be cited and the article DOIs given, which this section does.

### What the orosomucoid drugs are actually worth, and the treatment that follows

The user, 9 October 2026: "If orosomucoid carries 11 of 27 drugs by itself, we NEED a physiologically and
computationally viable way to deal with it ... If it just means dropping some portion of its pairs, check which
pairs and how important they are for psych symptoms. Knockouts for the orosomucoid binder would also be fine as
training examples." Both halves are measured in docs/orosomucoid_label_stake.md
(experiments/measure_orosomucoid_label_stake.py), and both answers cut against the worry.

**The pairs at stake are 23 of the study's 2,211 kept positives, 1.0 percent.** Seven of the eleven
orosomucoid-carried drugs hold no kept positive at all (alfentanil, desipramine, diphenhydramine, meperidine,
methadone, sufentanil, tiagabine), and diphenhydramine, the one of the eleven in the lockbox, is among them.
Fentanyl alone holds 12 of the 23. The symptoms they reach are anxiety (4), psychomotor_agitation, insomnia and
somnolence_or_hypersomnia (3 each), none of which depends on these drugs: every one is carried by many other
perturbations. So dropping the orosomucoid-carried pairs outright would cost the study about one part in a hundred
of its scored signal, concentrated in one opioid.

**The knockouts are not available.** A gene perturbation's labels come from HPO disease annotations, and ORM1 and
ORM2 have none: zero annotation rows and zero diseases, as do SHBG and SERPINA7. The binders that do carry HPO
annotations are TTR (100 rows), APOB (55), APOA1 (49), ALB (39), RBP4 (22), SERPINA6 (13) and AFP (4), and three of
those are already perturbations. Orosomucoid cannot be a training example until some phenotype source other than
HPO names it, which is a data gap and not a modelling choice. No human orosomucoid deficiency is described.

Sources searched for one, 9 October 2026, at the user's question ("did you look for additional sources to find
potential ORM1 / ORM2 knockouts?"):

| source | what it returned |
|---|---|
| HPO `genes_to_phenotype`, the study's own label source | 0 annotation rows and 0 diseases for ORM1 and for ORM2 |
| PubMed, for an ORM1 or ORM2 loss-of-function phenotype | nothing psychiatric. The nearest record is a mouse Orm2 study of adipose browning: "Orosomucoid 2 (Orm2) is identified as an IF-induced hepatokine that stimulates adipose browning", with "an obesity-associated Orm2 variant (D178E), which shows decreased GP130/IL23R binding and impaired browning capacity in mice" (abstract, PubMed 39248328, [DOI](https://doi.org/10.1002/advs.202407789), `zhu2024intermittent`) |
| IMPC, for a systematic mouse knockout phenotype | no phenotyping planned for Orm1 (MGI:97443), so no behavioural panel exists to read |
| commercial knockout lines | the described Orm2 phenotype is spontaneous obesity and hepatic steatosis, which is the same metabolic reading |

So the answer to whether the binder can be a training example is no on two counts rather than one: HPO annotates
neither gene, and the external literature describes the one existing knockout phenotype as metabolic. A psychiatric
label for ORM1 or ORM2 would have to be invented, and inventing one is what this study's grading rules exist to
prevent. These sources were retrieved from PubMed, which requires that it be cited and that the article DOI be given,
as above.

**The treatment.** Given those two numbers, a free-fraction input is the wrong instrument and a stratifier is the
right one:

1. The carriage edges stay as they are. They say which binder carries which cargo, which is structural and is
   already in the confirmatory graph. They do not claim a free fraction and do not need one.
2. The measured fraction unbound does not enter the model as a per-drug feature. It exists for 27 of 142 drugs, and
   the drugs it exists for are the better-studied ones; docs/off_target_scoping.md already measures that kept
   positives rise with measurement coverage (273, 374, 417 across the three terciles). A feature present for 19
   percent of drugs, and present exactly where the labels are denser, would carry how well studied a drug is as
   much as how much of it is free. That is the confounding the coverage strata exist to control, and a feature
   would reintroduce it inside the model.
3. The fraction unbound becomes a reported stratum instead, read the same way as the coverage terciles: the
   confirmatory numbers with and without the drugs whose free fraction moves most with the binder state. This
   changes no input, no graph and no label, so it needs no amendment and cannot disturb the lockbox.
4. If the binder state is ever to be modelled per drug, it needs the measured pair of fractions per population, not
   a concentration-rescaled one (the correction above), and it needs the 115 missing drugs, which no pinned source
   supplies. That is a later study, not this one.

### What is still missing

- The 115 drugs with no measured fraction unbound. ChEMBL does not fill the gap: its human plasma-binding rows
  (`standard_type` PPB and Fu, 2,260 and 1,759 rows) reach 19 of this study's 142 drugs, and its ADMET endpoint
  returns calculated properties only, which a call for one of these drugs confirms. Lombardo and colleagues' 1,352
  compound compilation (`10.1124/dmd.118.082966`) is the next candidate, and its protein-binding column is itself
  compiled "from other available sources", so it would need its own appraisal.
- Orosomucoid's own plasma concentration per population, which the database does not carry because it does not need
  it; the correction above is the reason this is no longer on the critical path.
- Whether the extension's second condition is modelled per drug or per population, which is a design question for the
  user, not a data gap.

## Stage 3, 9 October 2026: the carrier identity of the other drugs, and what the special-population column affords

Two questions from the user, both answered by measurement.

**"For the fraction unbound, the 'in special populations' is the part that i am leaning on there. maybe it doesn't
afford the leniancy i think it does?"** It affords less than it looks
(`experiments/measure_special_population_carrier_signal.py`, docs/special_population_carrier_signal.md). Scored
against the database's own carrier labels, the direction a drug's free fraction moves does separate the carriers, but
**the sign is not one sign**: a higher free fraction means orosomucoid in hepatic impairment (AUROC 0.742), in the
newborn (0.815) and across ethnic groups (0.925), and albumin in renal impairment (0.113), in the elderly (0.182) and
in inflammatory disease (0.248). That is the physiology, since orosomucoid is a hepatic acute-phase protein that
rises in inflammation and renal disease and falls in liver failure and at birth, but it means a free fraction with no
population attached carries nothing. Pooling the populations a drug has, leave-one-drug-out, gets 155 of 193 drugs
right, 80 percent, against a majority-class rate of 76 percent: it recovers 35 of the 47 orosomucoid carriers, which
the trivial rule finds none of, at the cost of calling 26 of 146 albumin-carried drugs orosomucoid. That error rate
is too high to place a carriage edge on and useful only for deciding which drugs to look up.

And the coverage it could add is zero, for a reason that has nothing to do with its accuracy: the classifier needs a
measured pair of fractions, which exists for the same 27 drugs whose carrier the database already labels. A ratio and
a carrier label come together.

**"I just want better coverage and i think it's possible to obtain in a valid way."** It is, from a source that
states the carrier rather than implying it: the drug label
(`experiments/fetch_fda_label_protein_binding.py`, `experiments/scope_label_plasma_binder_coverage.py`,
docs/label_plasma_binders.md). openFDA serves each label with a `set_id` and an `effective_time`, so the exact
version behind a sentence is recorded and the sentence can be quoted, which is the evidence standard the curated
carriage table already uses. Of the study's 154 drug perturbations, 130 have an openFDA label and **31 name a
carrier** (5 orosomucoid, 20 albumin, 6 both), **22 of them drugs the fraction-unbound database does not hold**. The
check that matters: 9 drugs are named by both sources and **all 9 agree**. Carrier identity therefore reaches 49 of
154 drugs rather than 27.

Those counts are the second reading. The first was 34 carriers over 11 agreeing drugs, and four of them were wrong in
ways a sentence-level rule has to refuse rather than report: oxcarbazepine was called an orosomucoid drug by the
sentence "Oxcarbazepine and MHD do not bind to alpha-1-acid glycoprotein", prilocaine by a sentence stating
lidocaine's binding, triazolam by one saying it does not displace bilirubin from albumin, and every phenytoin
sentence that mentions albumin at all mentions it as hypoalbuminemia, a dosing condition. The script now refuses a
negated binding or displacement claim, a binding clause whose subject is another drug in this study and the disease
name, which costs three drugs and one flip and is the reason to prefer a quoted sentence over a count.

Three cautions belong with the number that survives. A label states the carrier when it matters clinically, so the
drugs it names are not a random sample and the 31 cannot be read as an estimate of how many drugs orosomucoid
carries. Seven of the 31 rest on a sentence that does not name the drug itself (oxcarbazepine's names only its active
metabolite), and are marked in the table rather than counted silently. And a carrier identity is not yet an edge: a
xenobiotic has no node to attach one to, which is the subject of docs/drug_entry_nodes.md.

## Built

experiments/build_plasma_binder_variant.py (library: mechanistic_pathway_learning/graph/plasma_binding.py) wrote two
directories, reported in docs/plasma_binder_graph.md and in each directory's `plasma_binder_summary.json`:

- `data/processed/graph_full_neuronal_binders`: 27 edges added (266,966 -> 266,993), carriage on the gene nodes;
- `data/processed/graph_full_neuronal_split_binders`: the same 27 edges on the **protein** nodes, because on the
  gene/protein split graph binding is the protein's property and the gene node holds only expression.

30 (binder, cargo) rows were read, 8 from the UniProt features and 22 curated. All 30 mapped: no binder was missing
from the graph, no cargo failed the ChEBI route and no cargo lacked an extracellular copy. Three were already in the
graph from OmniPath (cortisol -> SERPINA6, thyroxine -> SERPINA7, thyroxine -> TTR), which is a check worth recording:
where the two sources overlap they agree, and the duplicate is dropped rather than added twice. Cargo connected:

| binder | cargo |
|---|---|
| ALB | palmitate, stearate, oleate, linoleate, arachidonate, thyroxine, bilirubin, Ca2+, Cu2+, zinc |
| ORM1, ORM2 | progesterone |
| SERPINA6 | cortisol, progesterone |
| SERPINA7 | thyroxine, triiodothyronine |
| TTR | thyroxine, triiodothyronine, retinol |
| RBP4 | retinol |
| SHBG | testosterone, estradiol-17beta, 5-alpha-dihydrotestosterone |
| AFP | Cu2+ |
| APOA1 | cholesterol |
| APOB | cholesterol, cholesterol-ester pool |

Node descriptors are unchanged, because the node set is: the source graph's `node_descriptors.parquet` is copied into
the variant. Node set 36,865 (merged) and 49,574 (split); only `degree` moves.

### What was deliberately left out, and why

- **Albumin's low-affinity buffering of every steroid.** Hammond reports that albumin "binds steroids with limited
  specificity and low affinity, but its high concentration in blood buffers major fluctuations in steroid
  concentrations and their free fractions" ([DOI](https://doi.org/10.1530/JOE-16-0070), PMID 27113851). True, and it
  would make albumin a hub over every steroid node for a binding no single site accounts for. The table keeps albumin
  to cargo it is a principal carrier of.
- **The lipoproteins as particles.** The carrier of cholesterol is the particle, not a site on one protein. The two
  apolipoprotein rows are an approximation stated in the table's own `mechanism` column, resting on there being one
  apolipoprotein B per VLDL or LDL particle (Fisher, [DOI](https://doi.org/10.1016/j.bbalip.2012.02.001),
  PMID 22342675) and on apolipoprotein A-I scaffolding HDL. Triacylglycerol is absent for a different reason: no ChEBI
  triglyceride term maps to a Human-GEM extracellular metabolite (the reconstruction has pools, not a triglyceride).
- **Alpha-fetoprotein's fatty acids and steroids.** The only source I could verify is hedged ("AFP *may* bind and
  transport fatty acids, steroids, heavy metal ions, drugs", Yin and Wang, PMID 12561448), and alpha-fetoprotein is
  fetal, so its carriage does not set the free fraction in the adults this study's labels come from. Only its UniProt
  Cu(2+) feature is used.
- **Orosomucoid's drugs.** Its UniProt function statement is explicit that the cargo is unnamed and largely
  xenobiotic: "Functions as a transport protein in the blood stream. Binds various ligands in the interior of its
  beta-barrel domain. Also binds synthetic drugs and influences their distribution and availability in the body." A
  graph with no drug nodes has nowhere to put that, so stage 1 gives orosomucoid only progesterone, the one endogenous
  cargo with a measured constant, and orosomucoid's real role in this study waits on the fraction unbound of stage 2.

### Checked against the trainer's own loader, 9 October 2026

The variant was read through `load_experiment_data` and `node_feature_matrix` on
`data/processed/graph_full_neuronal_split_binders`, with nothing trained. Node ids and `relation_types.json` are
identical to the source graph, the first 279,757 edge rows are the source's rows unchanged, and the 27 added rows are
all `binds`, all sign 0, all from an extracellular metabolite node to a `PROTEIN:` node. In the assembled feature
matrix one column moves, the log-degree column (index 16, after the six node types and ten compartments), on exactly
the 30 nodes the new edges touch; the perturbation and symptom counts are unchanged.

The descriptor table is copied rather than rebuilt, which holds only because carriage adds nothing a descriptor reads.
Two model inputs follow the wiring instead of an annotation: the `reaction_brain` block and the reactions' cell-class
weights both carry gene expression onto reactions through `catalyzed_by`
(`mechanistic_pathway_learning/graph/rewired_expression_features.py`). Carriage adds `binds` edges only, so neither
changes. A later variant that touched catalysis edges would have to recompute both, as the rewired runs do. Rebuilding
both variants a second time reproduced all seven parquet files byte for byte.

The preregistered descriptor input (`full_neuronal_split_descriptors_brain_expression.parquet`, 138 columns) is keyed
by node id and the node ids do not change, so it applies to the variant unchanged; the 84-column table inside each
graph directory is the plainer one the slice configurations read, not the confirmatory input.

## Open

- Whether the carriage edges go into the confirmatory graph or stay a variant. They change the graph, so they need a
  dated amendment and the development pilots re-run, which is the user's call; the user's instruction of 9 October 2026
  was "don't rerun the same pilots until we've completely nailed down the design options", so the variant waits.
- Whether a binder's own expression should gate the carriage edge, as the cell-class channels gate others.
- Whether the lipoproteins belong as one node or as the apolipoprotein genes they are now.
- Whether the binders' carriage should also reach the cargo's other compartments. Only the extracellular copy is joined
  now, which is right for plasma and wrong for the intracellular fatty-acid binding proteins, which are a different
  mechanism and not these genes.
