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

## Stage 2: the free fraction and the binder state

Needs data that is not here. A measured fraction unbound per drug is in neither DrugCentral's flat files
(`classes.tsv`, `drug-to-class.tsv`, `indications.tsv`, `targets.tsv`) nor the GtoPdb release, and DrugBank's
protein-binding field is in the full database rather than the pinned slim tables. Candidates to fetch and pin, in
order of how much they can be checked: ChEMBL assay rows with a plasma-protein-binding or fraction-unbound endpoint;
the open pharmacokinetic compilations that report fraction unbound with a source paper per drug; a hand-pinned table
from the labels, which is the least checkable and would need its own appraisal.

Only then does the extension the user named become possible, and it is a second condition per drug (normal plasma
against raised orosomucoid or low albumin) with its own evidence, not a descriptor column.

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

## Open

- Whether the carriage edges go into the confirmatory graph or stay a variant. They change the graph, so they need a
  dated amendment and the development pilots re-run, which is the user's call; the user's instruction of 9 October 2026
  was "don't rerun the same pilots until we've completely nailed down the design options", so the variant waits.
- Whether a binder's own expression should gate the carriage edge, as the cell-class channels gate others.
- Whether the lipoproteins belong as one node or as the apolipoprotein genes they are now.
- Whether the binders' carriage should also reach the cargo's other compartments. Only the extracellular copy is joined
  now, which is right for plasma and wrong for the intracellular fatty-acid binding proteins, which are a different
  mechanism and not these genes.
