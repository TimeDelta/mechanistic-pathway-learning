# Plasma substrate binders

The user's decision of 9 October 2026: "All of the plasma substrate binders should be modeled. Once AAG and albumin
are added, there will be state to attach that separate extension to." This document is the scope, the audit of what
this repository holds and the build in two stages. Nothing here is built yet, and nothing in the confirmatory
specification depends on it.

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
| lipoproteins | APOA1, APOB and the other apolipoproteins | cholesterol, triacylglycerols, lipophilic drugs |
| alpha-fetoprotein | AFP | fetal carriage of fatty acids and bilirubin |

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

Build: `binds` edges, sign 0, between each binder's node and the extracellular metabolite nodes it carries, in a new
graph directory rather than over `graph_full_neuronal` or its split, so no preregistered input changes.

Source: no pinned file in this repository carries ligand annotations. The audit:
`data/raw/uniprot/uniprot_human_reviewed.tsv.gz` holds eight columns (Entry, primary gene name, HGNC, EC number,
subcellular location, Pfam, length, sequence) and no binding site, ligand or function text; the ChEMBL tool's ADMET
endpoint returns calculated properties only. So stage 1 needs one fetch: the UniProt REST query for these accessions
with the binding-site and function columns (`ft_binding`, `cc_function`), whose ligand entries carry ChEBI
identifiers, which the Reactome import already maps to Human-GEM metabolites
(`map_chebi_to_human_gem`). Where UniProt names a ligand class rather than a compound (long-chain fatty acid), the
edge set is the class's members already in the graph.

## Stage 2: the free fraction and the binder state

Needs data that is not here. A measured fraction unbound per drug is in neither DrugCentral's flat files
(`classes.tsv`, `drug-to-class.tsv`, `indications.tsv`, `targets.tsv`) nor the GtoPdb release, and DrugBank's
protein-binding field is in the full database rather than the pinned slim tables. Candidates to fetch and pin, in
order of how much they can be checked: ChEMBL assay rows with a plasma-protein-binding or fraction-unbound endpoint;
the open pharmacokinetic compilations that report fraction unbound with a source paper per drug; a hand-pinned table
from the labels, which is the least checkable and would need its own appraisal.

Only then does the extension the user named become possible, and it is a second condition per drug (normal plasma
against raised orosomucoid or low albumin) with its own evidence, not a descriptor column.

## Open

- Whether the carriage edges go into the confirmatory graph or a later one. They change the graph, so they need a dated
  amendment and the development pilots re-run, which is the user's call.
- Whether a binder's own expression should gate the carriage edge, as the cell-class channels gate others.
- Whether the lipoproteins belong as one node or as the apolipoprotein genes they are now.
