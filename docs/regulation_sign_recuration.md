# Recuration of disputed regulation signs (17 October 2026)

A dated record for decision 9.72 in `docs/open_decisions.md`. It is not updated; a later reading goes in a section of
its own below the last one.

## Why

The graph takes a regulation's sign from OmniPath's consensus, and OmniPath's CollecTRI download holds two releases
of each regulation it covers: the 2023 release (resource CollecTRI, which states how it decided each sign) and a
later one (resource CollecTRI2, which states no decision). On the v5 graph they state opposite signs for 3,795 of the
41,638 regulation edges, 1,930 of them labels under decision 9.63, whose sign no gate can flip. The question for this
record is which sign the cited papers support when the releases disagree.

## The sample

The frame is the 1,930 labels whose releases state opposite signs, in two strata by the release the stored sign
follows: 1,502 follow CollecTRI2 and 428 the 2023 release. Eight labels were drawn from each stratum with pandas'
`sample(8, random_state=20261017)` over the session's table of regulation edges (one row per v5 edge, in edge order).
For each pair the PubMed ids came from the pinned evidences download
(`data/raw/omnipath/collectri_interactions_evidences.json`, SHA-256 `60b79abc…`): those each release puts behind its
own sign, split into the ids only one release cites and the ids both cite.

## How the papers were read

A research pass read the PubMed record of every cited paper (all of them when a list held eight or fewer, otherwise
the first eight of each list and the shared ones), 84 paper records over 73 distinct ids, 68 of them cited by a
release and 5 found by search (below), and classed each
as "activates", "represses", "both", "indirect or not transcriptional" (an intermediate is named, or the effect is on
the target protein's activity, or it is a correlation), "not about this regulation" or "no abstract", copying the
sentence that supports the class verbatim. A pair's verdict is the direct transcriptional sign the cited papers
support: activation, repression, both reported, or undetermined. Where no cited paper was informative the pass ran one
PubMed search for the pair and read up to five results, marked as found by search. Every cited id returned a record
with an abstract. Of the 84 records, 18 were classed "activates", 13 "represses", 1 "both", 20 "indirect or not
transcriptional" and 32 "not about this regulation".

## Verdicts

| Regulation | Stored sign follows | 2023 release | CollecTRI2 | Verdict | Stored right | 2023 right | CollecTRI2 right |
| --- | --- | --- | --- | --- | --- | --- | --- |
| THAP11 to CDKN1A | CollecTRI2 | +1 (paper) | -1 | undetermined | - | - | - |
| SMAD2 to MUC4 | CollecTRI2 | +1 (regulon) | -1 | activation | no | yes | no |
| CBFA2T3 to TCF12 | CollecTRI2 | +1 (default) | -1 | repression | yes | no | yes |
| MYC to CYP3A4 | CollecTRI2 | +1 (default) | -1 | repression | yes | no | yes |
| CREM to INS | CollecTRI2 | +1 (default) | -1 | repression | yes | no | yes |
| RUNX1 to LYZ | CollecTRI2 | +1 (default) | -1 | undetermined | - | - | - |
| FLI1 to MDM2 | CollecTRI2 | +1 (default) | -1 | activation | no | yes | no |
| SOX5 to RUNX2 | CollecTRI2 | +1 (paper) | -1 | repression | yes | no | yes |
| SP3 to SLC9A2 | the 2023 release | +1 (default) | -1 | both reported | - | - | - |
| STAT5A to NR3C1 | the 2023 release | +1 (default) | -1 | undetermined | - | - | - |
| SP1 to PRKCA | the 2023 release | +1 (default) | -1 | activation | yes | yes | no |
| RELA to AKT1 | the 2023 release | -1 (paper) | +1 | activation | no | no | yes |
| NFE2L2 to BACH1 | the 2023 release | +1 (paper) | -1 | activation | yes | yes | no |
| USF1 to SHBG | the 2023 release | +1 (default) | -1 | both reported | - | - | - |
| NFKB1 to AKT1 | the 2023 release | -1 (paper) | +1 | activation | no | no | yes |
| SMAD1 to HAND1 | the 2023 release | +1 (paper) | -1 | activation | yes | yes | no |

"(paper)" is the 2023 release's sign decision "PMID", "(default)" its "default activation" and "(regulon)" its regulon
inference. Of the 11 pairs with one supported sign, the stored sign is right for 7, the 2023 release for 5 and
CollecTRI2 for 6. Of the five where the 2023 release took its sign from a paper and the literature decides, it is
right for 2 (NFE2L2 to BACH1, SMAD1 to HAND1) and wrong for 3 (SOX5 to RUNX2, RELA to AKT1, NFKB1 to AKT1). The other
five pairs are undetermined (3) or report both signs by cell type or species (2): Sp3 activates the human NHE2
promoter and represses Sp1-driven activation in mouse collecting-duct cells, and the USF factors repress the human
SHBG promoter in Sertoli cells and induce the rat one. The SP1 to PRKCA verdict rests on a paper found by search, not
on a cited one.

## Papers cited behind the opposite of what they report

According to PubMed, in 11 of the 16 pairs a release cites at least one paper whose abstract reports the opposite of
the sign the release puts it behind:

- THAP11 to CDKN1A, cited by CollecTRI2 behind -1: "In addition, THAP11 negatively regulated the expression of c-Myc,
  decreased cyclinD1 protein, and increased p27 and p21 protein levels" (Zhang et al. 2020,
  [doi:10.1155/2020/7838924](https://doi.org/10.1155/2020/7838924)); "Moreover, overexpression of THAP11 significantly
  decreased the colony numbers, and also inhibited the expression of c-myc target genes such as Cyclin D1, ODC and
  induced the expression of p21(Cip1)" (Nakamura et al. 2011, [doi:10.1002/ijc.26065](https://doi.org/10.1002/ijc.26065)).
  Both attribute the rise to THAP11 repressing c-Myc, which is why the verdict is undetermined.
- SMAD2 to MUC4, cited by CollecTRI2 behind -1: "Smad2 positively cooperates with Smad4 to activate the promoter"
  (Jonckheere et al. 2004, human pancreatic cancer cell lines, [doi:10.1038/sj.onc.1207769](https://doi.org/10.1038/sj.onc.1207769)).
- CBFA2T3 to TCF12, cited by the 2023 release behind its default +1: "Our results show that the E-box gene HEB is a
  direct target of CBFA2T3-ZNF652-mediated transcriptional repression" (Kumar et al. 2008,
  [doi:10.1074/jbc.M709136200](https://doi.org/10.1074/jbc.M709136200)).
- MYC to CYP3A4, cited by CollecTRI2 behind -1: "Furthermore, inhibition of c-myc may indirectly influence the
  expression of CYP3A" (Arora et al. 2000, rat liver, PMID 10688605, no DOI in PubMed); c-myc antisense lowered CYP3A
  activity after partial hepatectomy.
- FLI1 to MDM2, both releases: CollecTRI2 behind -1, "In this paper, we demonstrate that the promoter of MDM2 contains a
  consensus binding site for Fli-1 that is bound by this transcription factor in vitro and in vivo, resulting in MDM2
  transcriptional regulation", with "a positive correlation in the expression of Fli-1 and MDM2" in mouse leukemic
  tissue (Truong et al. 2005, [doi:10.1038/sj.onc.1208323](https://doi.org/10.1038/sj.onc.1208323)); the 2023 release
  behind +1, "The introduction of EWS-Fli1 resulted in significant reduction of promoter activities and mRNA levels of
  p21 and mdm2" (Li et al. 2010, [doi:10.1016/j.canlet.2010.01.022](https://doi.org/10.1016/j.canlet.2010.01.022)),
  which concerns the EWS-Fli1 fusion acting through p53.
- SP3 to SLC9A2, both releases: the 2023 release behind +1, "Furthermore, in mIMCD-3 cells, Sp1-mediated transcriptional
  activation was repressed by expression of Sp3 and Sp4" (Bai et al. 2001,
  [doi:10.1152/ajpcell.2001.280.5.C1168](https://doi.org/10.1152/ajpcell.2001.280.5.C1168)); CollecTRI2 behind -1, "The
  results in the present study indicate that Sp1 and Sp3 are required for constitutive NHE2 expression" (Pearse et al.
  2007, [doi:10.1042/BJ20070364](https://doi.org/10.1042/BJ20070364)).
- RELA to AKT1 and NFKB1 to AKT1, the 2023 release behind -1: "This was due to a RhoA-dependent retention of NFkappaB to
  the cytosol at 12 h which led to a transcriptional down-regulation of the anti-apoptotic protein BCL-2 as well as
  reduced AKT1 mRNA production at 24 h and thus diminished levels of PKB/AKT protein" (Aberg et al. 2007, human
  MDA-MB-231 cells, [doi:10.1016/j.thromres.2007.09.017](https://doi.org/10.1016/j.thromres.2007.09.017)). CollecTRI2,
  behind +1 for RELA to AKT1, cites "TNF-α activated p65 and MAPKs (p38, JNK1/2 and ERK1/2), iNOS and COX-2 which
  worsened the insulin signaling expressions of IR, IRS-1, IRS-2, PI3K-α, Akt/PKB, GLUT-2, AMPK, GCK, pyruvate kinase,
  and PPAR-γ" (Chen et al. 2019, HepG2 cells, [doi:10.1016/j.jsps.2019.01.014](https://doi.org/10.1016/j.jsps.2019.01.014)).
  The activation verdict rests on "Interestingly, in addition to stimulating the phosphorylation of Akt, overexpression
  of p65 causes an increase in the expression of Akt mRNA and protein" (Meng et al. 2002,
  [doi:10.1074/jbc.M112464200](https://doi.org/10.1074/jbc.M112464200)) and "Increased Akt expression by p65 is
  mediated at the transcriptional level" (Meng and D'Mello 2003, mouse NIH3T3 and primary endothelial cells,
  [doi:10.1016/j.bbaexp.2003.09.001](https://doi.org/10.1016/j.bbaexp.2003.09.001)), both cited by CollecTRI2.
- NFE2L2 to BACH1, cited by CollecTRI2 behind -1: "In the present study, we have demonstrated that BACH1 itself is
  regulated by Nrf2 as it is induced by Nrf2 overexpression and by Nrf2-activating agents in an Nrf2-dependent manner"
  (Jyrkkänen et al. 2011, [doi:10.1042/BJ20110526](https://doi.org/10.1042/BJ20110526)); "In GDM cells, the lipid
  peroxidation product 4-hydroxynonenal (HNE) failed to induce nuclear Nrf2 accumulation and mRNA and/or protein
  expression of Nrf2 and its target genes NAD(P)H:quinone oxidoreductase 1 (NQO1), Bach1, cystine/glutamate transporter,
  and glutamate cysteine ligase" (Cheng et al. 2013, [doi:10.2337/db13-0169](https://doi.org/10.2337/db13-0169)); "Nrf2
  is under feedback regulation through its target gene, transcriptional repressor Bach1" (Hushpulian et al. 2023,
  [doi:10.1002/bies.202300176](https://doi.org/10.1002/bies.202300176)).
- USF1 to SHBG, both releases: the 2023 release behind +1, the USF factors "repress human SHBG transcription in Sertoli
  cells through an interaction with FP4 within its proximal promoter" (Selva et al. 2004,
  [doi:10.1074/jbc.M409616200](https://doi.org/10.1074/jbc.M409616200)); CollecTRI2 behind -1, "Together, these data
  indicate that increased USF levels induce the expression of Nr5a1 and Shbg during the differentiation of Sertoli
  cells" (Wood et al. 2011, rat, [doi:10.1095/biolreprod.111.093013](https://doi.org/10.1095/biolreprod.111.093013)).
- SMAD1 to HAND1, cited by CollecTRI2 behind -1: "Chromatin immunoprecipitation (ChIP) revealed a subset of the BIG
  signature, including Satb2, Smad6, Hand1, Gadd45γ and Gata3, that was bound by Smad1/5 in the developing mandible,
  revealing direct Smad-mediated regulation" of Bmp4-induced genes (Bonilla-Claudio et al. 2012,
  [doi:10.1242/dev.073197](https://doi.org/10.1242/dev.073197)); "In contrast, inhibition of SMAD1 signalling acted to
  maintain SOX2 expression and prevent non-neural differentiation via HAND1" (Greber et al. 2011,
  [doi:10.1038/emboj.2011.407](https://doi.org/10.1038/emboj.2011.407)).

## Other faults in the citations

- SOX5 to RUNX2: a paper both releases cite, and the one the 2023 release's paper-based activation rests on, reports
  Sox6, not SOX5, at the Runx2 enhancer: "Dlx5, Mef2c, Tcf7, Ctnnb1, Sp7, Smad1, and Sox6, which localized on the
  enhancer region in primary osteoblasts, synergistically upregulated the enhancer activity" (Kawane et al. 2014,
  [doi:10.1002/jbmr.2240](https://doi.org/10.1002/jbmr.2240)).
- RUNX1 to LYZ: both cited papers test the AML1/ETO fusion oncoprotein, not wild-type RUNX1.
- STAT5A to NR3C1: the cited papers concern the STAT5 and glucocorticoid receptor proteins interacting, not NR3C1
  transcription.
- RELA and NFKB1 to AKT1: by the pass's reading, 8 of the 9 papers the 2023 release cites describe the reverse
  direction, Akt acting on NF-kappaB.

## How this record was checked

The session re-read 21 of the quoted records through PubMed (every quote in the two sections above except the
fusion and protein-interaction notes) and found each quote verbatim in its abstract. It did not repeat the
classification of the papers it did not quote, and it did not re-make the verdicts; those are the pass's. The pass
reported no text addressed to an AI assistant in any abstract, and no request that went unanswered.

## What this does not show

Sixteen pairs is a sample, not an estimate with useful precision: 7 of 11 correct stored signs is consistent with
anything from roughly a third to nine tenths. The verdicts are of the cited papers' abstracts, not of the whole
literature, and abstracts state the cell type and species unevenly. What the sample does show is that neither release
is a dependable tie-breaker: each is wrong for 5 or 6 of the 11 decidable pairs, and the 2023 release's sign decided
from a paper is wrong for 3 of 5.

## Sources

The release counts and the sample come from the pinned OmniPath downloads named in `docs/data_sources.md`; the papers
are in `docs/references.bib` under decision 9.72.
