# Alleles as perturbations (decision 9.54 (d), built 16 October 2026)

The owner, 16 October 2026 (19:27 EDT): "Huntington's proposal d accepted". This record says what was built, from
which sources, what it changed in the v5 tables and what is still wrong. The code is
`mechanistic_pathway_learning/evidence/allele_mechanisms.py`; the assembler's `--allele-perturbations` uses it.

## 1. Why

A gene was one perturbation with one sign, seeded −1 unless every causal association behind it was Orphanet's gain
of function type. HTT was seeded as a loss of function, while Huntington disease is a toxic gain of function of the
CAG expansion, and HTT's intellectual-disability label comes from another allele altogether: OMIM:617435, an
autosomal recessive neurodevelopmental disorder from compound heterozygous variants (Rodan et al. 2016,
doi 10.1038/ejhg.2016.74). A gene whose diseases arise from alleles acting in opposite directions could not be
represented with one sign.

## 2. The rule

Each gene's monogenic reports are split by the mechanism of the disease they come from, and each mechanism class is
a perturbation of its own:

| class | perturbation id | seed |
| --- | --- | --- |
| loss of function, and untyped | the gene | −1 on the gene node, as before; an untyped association is silence about the mechanism and keeps the default, marked untyped |
| gain of function, hypermorph | `<gene>\|gain_of_function` | +1 on the gene node: more of the normal activity |
| gain of function, neomorph | the curated allele, `HTT\|CAG_repeat_expansion` | its own entry node, present only where it is seeded, whose edges state what the mutant product does |

Every allele of a gene keeps the gene as its leakage group and its disease cluster, so no split puts one allele in
training and another in the test. Literature statements about a gene with several alleles cannot be attributed to
one of them and are set aside in the trainer's soft priors; IMPC knockouts stay with the loss-of-function allele.
`--allele-perturbations` replaces `--seed-gain-of-function-genes-positive`, and the assembler refuses both together.

## 3. Where a disease's mechanism comes from

In order, each step consulted only where the ones before it left the pair untyped:

1. **Orphanet's own type** of the (gene, ORPHA:x) association, Orphadata product 6: "Disease-causing germline
   mutation(s) (loss of function) in" or "(gain of function) in".
2. **The same disease's other record**, through Orphadata product 1's alignments, whose names in the file are
   "E (Exact mapping: the two concepts are equivalent)", "NTBT (ORPHAcode is narrower than the targeted code used to
   represent it)" and "BTNT (ORPHAcode is broader than the targeted code used to represent it)". A type is carried
   from an Orphanet disorder to each OMIM entry mapped E to it, and from such an entry to each Orphanet disorder
   mapped NTBT to it (a narrower form of the same disease: juvenile Huntington disease, ORPHA:248111, is NTBT to
   OMIM:143100). Never through BTNT: a broader Orphanet code is a group (ORPHA:528084, non-specific syndromic
   intellectual disability, maps BTNT to 108 OMIM entries). Two inherited types that disagree leave the pair untyped
   as a conflict.
3. **Gene2Phenotype**, the second source (the owner, 15 October 2026: "Always try to find a second reputable source
   when one fails"). G2P records its curators' mechanism in "broad categories depending on the protein level
   consequence of variants reported. These include 'Loss of function', 'Dominant negative' and 'Gain of function'.
   There are also categories for 'Undetermined' and 'Undetermined non-loss-of-function'" (Yates et al. 2024,
   doi 10.1186/s13073-024-01398-1). A record speaks for an OMIM disease and reaches an Orphanet disorder through the
   same E and NTBT alignments. Only definitive, strong and moderate gene-disease confidence types a pair. A
   dominant-negative product is read as a loss of function, because it lowers the activity of the wild-type product
   ("Dominant-negative (DN) mutations interfere with the activity of the wild-type protein", Badonyi and Marsh 2025,
   doi 10.1038/s41467-025-63234-3) and so the gene's net activity falls, which is what −1 states. Undetermined types
   nothing and is written into the basis. G2P records that disagree with each other leave the pair untyped as a
   conflict (KCNA2, MAB21L2 and RARB, each with a monoallelic gain and a loss for one OMIM entry).
4. **Adjudication from the primary literature**, for the pairs the two sources dispute
   (`configs/allele_mechanism_adjudications.csv`, one row per pair with a PubMed id and a quote). Where both sources
   type a pair and disagree, the first stands unless a row settles it. The six disputed pairs of 16 October were read:

| gene, disease | Orphanet | G2P | the primary paper | adjudicated |
| --- | --- | --- | --- | --- |
| CDC73, hyperparathyroidism-jaw tumour syndrome (ORPHA:99880) | gain | loss | "we identified thirteen different heterozygous, germline, inactivating mutations in a single gene in fourteen families with HPT-JT" (Carpten et al. 2002, PMID 12434154) | loss; Orphanet wrong |
| PRPS1, mild PRPP synthetase superactivity (ORPHA:411536) | loss | gain | "PRS1 concentrations in cells from six affected males exceeded those in normal cells by 2-6-fold" (Becker et al. 1996, PMID 8702702) | gain, a hypermorph; Orphanet wrong |
| ARHGAP31, Adams-Oliver syndrome (ORPHA:974, OMIM:100300) | gain | loss | "Mutant transcripts are stable and increase ARHGAP31 activity in vitro through a gain-of-function mechanism." (Southgate et al. 2011, PMID 21565291) | gain; G2P wrong |
| KCNC1, progressive myoclonic epilepsy 7 (OMIM:616187) | loss | gain | "Functional analysis of the Arg320His mutant channel showed a dominant-negative loss-of-function effect." (Muona et al. 2015, PMID 25401298) | loss; G2P wrong |
| SLC40A1, hemochromatosis type 4 (OMIM:606069) | gain | loss (biallelic) | "Hemochromatosis type 4 is caused by dominant gain-of-function mutations of ferroportin preventing hepcidin-ferroportin binding and leading to hepcidin resistance." (Piperno et al. 2020, a review, PMID 32258529); HPO gives the entry autosomal dominant inheritance | gain; G2P's biallelic record does not fit the dominant entry |

The quotes were checked against the abstracts PubMed returned in the session that wrote them; the abstracts are not
stored in the repository, so no test re-reads them.

## 4. The curated neomorph

`configs/neomorphic_alleles.csv` holds HTT's CAG expansion as an entry node, `ALLELE:HTT_CAG_repeat_expansion`, with
four `targets` edges, each with a PubMed id and a quote: CREBBP decreased (PMID 11264541), TP53-driven transcription
decreased (PMID 10823891), PPARGC1A transcription decreased (PMID 17018277, an edge to the gene node) and REST's
nuclear action increased (PMID 12881722, the wild-type function the mutant loses). Magnitudes are 1 where the source
states none, as a knockout's seed is 1. On a graph with entry nodes the loader seeds the allele's own node; on a graph
without them it seeds the four targets. The v5 drug and allele graph,
`data/processed/graph_full_neuronal_split_binders_drugs_v5`, holds 49,797 nodes and 268,958 edges: the 222 drug nodes
of v4 with their 1,058 mechanism and 101 sequestration edges, and the allele's node with its 4 edges (evidence family
other, so ungated by the edge priors).

## 5. What it changed in the v5 tables (measured 16 October 2026)

Each build was run by `scripts/rebuild_data_layer.sh fine`; the two intermediate builds are kept under
`data/processed/superseded/` with their hashes.

| | without alleles | Orphanet only | Orphanet and G2P | and the adjudications (current) |
| --- | --- | --- | --- | --- |
| perturbations (genes) | 2,696 (2,474) | 2,745 (2,523) | 2,761 (2,539) | 2,760 (2,538) |
| leakage groups, disease clusters | 2,563, 1,698 | 2,563, 1,698 | 2,563, 1,698 | 2,563, 1,698 |
| positive pairs, kept | 8,059, 4,418 | 8,116, 4,442 | 8,124, 4,449 | 8,123, 4,448 |
| gene perturbations with a kept positive | 1,702 | 1,729 | 1,741 | 1,740 |
| alleles: loss, untyped, hypermorph, neomorph | | 311, 2,131, 80, 1 | 844, 1,579, 115, 1 | 845, 1,578, 114, 1 |
| genes with several alleles | | 49 | 65 | 64 |
| perturbations whose only kept label is intellectual disability | 757 | 774 | 782 | 782 |

The groups and clusters do not move, as the rule requires. Of the 4,537 (gene, disease) pairs behind the reports, G2P
alone types 1,086: 1,023 loss of function (45 of them dominant negative), which leaves their seed as it was and their
class known, and 63 gain of function, which change sign. The 63 are in ABCC9, ACTB, ACTG1, AKT1, BRAF, CACNA1G, CBL,
CCND2, CERT1, CHRNA4, CSNK2A1, EEF1A2, FGF12, FGFR1, FGFR2, FGFR3, FLNA, GLUD1, GMNN, GNAQ, GNB1, IRF6, KCNK4,
KCNMA1, KCNN3, KCNQ2, KRAS, MAP2K1, MAPK1, MAST3, MECOM, NFE2L2, NLRP3, NRAS, PAK1, PCGF2, PDE4D, PDGFRB, PIK3R2,
PRKAR1A, PRPS1, SETD2, SHOC2, SMAD4, SMO, SOS1, TBL1XR1, TGFBR1, TSHR, TTR and TUBB4A.

Examples on the kept labels: HTT keeps agitation, bradykinesia and intellectual disability, all three from its
recessive OMIM:617435 (Lopes-Maciel-Rodan syndrome), and its CAG expansion carries the thirteen labels of Huntington
disease. FMR1 keeps the fragile X labels, and its gain-of-function allele the labels of the fragile X-associated
tremor/ataxia syndrome (OMIM:300623, ORPHA:93256). LRRK2's Parkinson labels are all on its gain-of-function allele.

## 6. What is still wrong or open

- **The untyped default.** 2,763 pairs remain untyped and are seeded −1. By the inheritance HPO gives each disease,
  517 come only from recessive diseases, where −1 is the expected mechanism; 549 from diseases that are never
  recessive; 1,661 from diseases with no inheritance stated, mostly Orphanet records, which HPO annotates without it.
  For dominant phenotypes Badonyi and Marsh predict, from protein structure rather than curation: "Single-phenotype
  AD genes exhibited the largest fraction of phenotypes with a LOF mechanism, at 54.6%. The remaining fraction was
  attributed to GOF and DN mechanisms occurring at similar frequencies, at 23.8% and 21.6%, respectively. In
  multi-phenotype AD genes, the fraction of phenotypes with a LOF mechanism was lower, at 48.1%, followed by GOF at
  32.5% and DN at 19.3%." A dominant negative keeps −1 right, so on that estimate −1 has the wrong sign for roughly a
  quarter to a third of the never-recessive pairs. SNCA is the clearest case: its Parkinson disease and
  dementia with Lewy bodies records (OMIM:168601, 605543, 127750) are untyped in both sources, because G2P has no
  adult neurodegeneration panel, and so its main labels are seeded as a loss of function. Orphadata product 9 states
  each Orphanet disorder's inheritance and would fill most of the 1,661. Options for the owner, none built: keep −1;
  set the never-recessive untyped pairs aside; or admit the predicted mechanism as a third, lower-tier source.
- **Every gain of function is a hypermorph unless curated.** Several of the 63 above are toxic or new functions,
  not more of the normal one (TTR amyloidosis, TUBB4A, and FMR1's premutation, whose FXTAS is attributed to RNA
  toxicity), so +1 on the gene node misstates them. They are the candidates for `configs/neomorphic_alleles.csv`,
  with SNCA, PACS1 and ACTB and ACTG1.
- **G2P is a download made on the fly**, pinned by its hash and date in `docs/data_sources.md`; a later fetch
  reports DRIFT. Its static release directory is closed to automated readers by its robots file.
- **A correction to the decision log.** Decision 9.54's text of 16 October said G2P's licence was "not yet verified".
  It had been verified on 11 October (`docs/label_source_review.md`, sections 11.2 and 32).
