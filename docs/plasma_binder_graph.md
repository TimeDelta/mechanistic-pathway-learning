# Plasma substrate binders in the graph

Built by experiments/build_plasma_binder_variant.py from `data/processed/graph_full_neuronal` into `data/processed/graph_full_neuronal_binders`. 27 `binds` edges, sign 0, each from the extracellular copy of a cargo metabolite to its binder, carried on gene nodes, which is the direction and sign the graph's own small-molecule edges use.

The same run built `data/processed/graph_full_neuronal_split_binders` from the matching source graphs. On a gene/protein split graph the carriage targets the protein node, because binding is the protein's property and the gene node there holds only expression, so the counts above hold but the edge endpoints differ.

## Cargo connected, by binder

| binder | cargo | edges |
|---|---|---|
| AFP | Cu(2+) | 1 |
| ALB | (4Z,15Z)-bilirubin IXalpha, Ca(2+), Cu cation, Zn(2+), arachidonate, linoleate, oleate, palmitate, stearate, thyroxine | 10 |
| APOA1 | cholesterol | 1 |
| APOB | cholesterol, cholesteryl ester | 2 |
| ORM1 | progesterone | 1 |
| ORM2 | progesterone | 1 |
| RBP4 | retinol | 1 |
| SERPINA6 | cortisol, progesterone | 2 |
| SERPINA7 | thyroxine, triiodothyronine | 2 |
| SHBG | 5alpha-dihydrotestosterone, estradiol-17beta, testosterone | 3 |
| TTR | L-thyroxine, retinol, triiodothyronine | 3 |

By source: {'pubmed': 17, 'uniprot_binding_site': 8, 'uniprot_function': 2}. Rows read: 30 (22 curated, 8 from UniProt binding-site features).

## Binder genes absent from the graph

None.

## Cargo with no Human-GEM metabolite

None.

## Cargo with no extracellular copy in this graph

None.

## Edges the graph already held

- MAM01615e -> GENE:SERPINA6
- MAM02998e -> GENE:SERPINA7
- MAM02998e -> GENE:TTR

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
Those three changes, what the record holds about an earlier intention to make them and what they would cost are set
out in docs/drug_entry_nodes.md, written on 9 October 2026 at the user's question; the carrier identity they would
need now reaches 50 of the 154 drugs rather than 27 (docs/label_plasma_binders.md).

