# Gene and protein nodes: the split and the gene-to-protein mapping

Written on 8 October 2026. Nothing here is built or registered. No graph directory was regenerated; the measurements
read graph_full_neuronal, evidence_full_v2 with better_v2 and the pinned reviewed UniProt table (20,431 entries). The
degree rankers below read development outcomes only (the 317 lockbox perturbations are left out); lockbox counts are
membership counts.

## The user's statements

- "The genes shouldn't need descriptors but the proteins should."
- "I'm thinking for number two that the seed [masking] isn't necessary for the message passing after gene descriptors
  are removed."
- "But doesn't that assume a 1:1 mapping between gene and protein?"
- "If the gene descriptors are dropped, I would specify that the brain expression is still provided via pseudo-gene
  descriptors. The rest of that suggested split between the genes and proteins makes sense I think except the
  registered degree strata SHOULD change to the new split."

## Is the mapping one to one?

No. The merged gene node of the current graphs and the split as first proposed (one PROTEIN:X per GENE:X) both assume
it. Departures among the 12,810 gene nodes:

| case | gene nodes | examples | perturbation seeds (membership only) |
|---|---|---|---|
| one reviewed entry encoded by several genes | 57 | histone H2A, H2B and H3.3 genes (H3-3A, H3-3B), defensins, CKMT1A and CKMT1B | 1 development knockout |
| one gene with several reviewed entries (distinct products) | 37 | NRXN1, NRXN2 and NRXN3 (alpha and beta neurexins), GNAS, CDKN2A (p16INK4a and p14ARF), CALCA, PRNP, TSPO | 4 development knockouts, 1 development drug, 5 lockbox knockouts |
| isoforms inside an entry | not in the pinned table | most multi-exon genes | all |
| one entry cut into several peptides | not counted | POMC: edges to MC1R to MC5R (the melanocortins) and to OPRD1 and OPRK1 (beta-endorphin) on one node; PDYN, PENK | not counted |
| no reviewed protein | 113 | 19 ChEBI compounds typed as gene; 58 unreviewed UniProt accessions used as symbols; 36 outdated or malformed symbols (GBA for GBA1; nodes named "HBA1; HBA2", "CCL3L1; CCL3L3" and "CCL4L1; CCL4L2") | none |

The 113 nodes without a reviewed protein take part in 2,170 edges (activates 1,141, inhibits 713, binds 254,
catalyzed_by 40, regulates_transcription_of 14, member_of 8).

Isoforms are the largest departure. pan2008deep: "transcripts from approximately 95% of multiexon genes undergo
alternative splicing". yang2016widespread: "The majority of isoform pairs share less than 50% of their interactions" and
"alternative isoforms tend to behave like distinct proteins rather than minor variants of each other".

The sources do not resolve products below the gene:
- The OmniPath interaction file holds no isoform accession (0 of 126,125 rows on either side).
- It copies a gene's interactions onto each of the gene's reviewed entries. The partner lists are identical for CDKN2A's
  two entries (80 rows each), NRXN1's two (45 each) and CALCA's two (44 each). p16INK4a's partners CDK4 and CDK6 and
  p14ARF's partner MDM2 therefore sit on both entries.
- Interactions are collected by gene symbol, CollecTRI regulates genes, and ChEMBL single-protein targets are canonical
  entries.

A protein node finer than the gene would carry the same edges as its siblings. The one-to-one mapping is therefore a
limit of the data as much as of the design, and the split cannot remove it.

Defect found: experiments/build_protein_descriptors.py groups entries by UniProt's primary gene-name field. 68 entries
list several genes ("DEFB4A; DEFB4B"), so the 57 gene nodes they cover get no protein descriptors in any current graph.
It is not fixed under the running jobs; the split build fixes it.

## Handling in the split (proposed defaults, for the user)

1. **One protein node per gene,** except that genes encoding one shared entry get one shared protein node, with an
   encodes edge from each gene (H3-3A and H3-3B to H3.3). A knockout of one of them then leaves the other's route to the
   protein.
2. **Genes with several entries keep one protein node.** The sources give the entries the same edges, and its protein
   descriptor stays the mean over the entries, as now. Isoforms and peptides of one precursor are not represented. All
   three are stated as limitations.
3. **The 113 nodes without a reviewed protein.**
   - ChEBI compounds are retyped as metabolites.
   - Unreviewed accessions and outdated symbols are mapped to the current HGNC symbol where HGNC's previous symbols or
     UniProt's cross-references resolve them.
   - Any that stay unresolved remain protein nodes without descriptors.
4. **Leakage groups.** A protein node maps to the genes that encode it, and groups join through that map, so folds and
   lockbox membership stay as drawn (a test checks this).
5. **Brain expression.** The 27-column brain expression block (GTEx v10, HPA brain regions and HPA single-nucleus cell
   classes) is mRNA measured per gene. The default reading of "pseudo-gene descriptors": this block moves to the protein
   node as gene-derived descriptors, a proxy for protein abundance (liu2016dependency: "transcript levels by themselves
   are not sufficient to predict protein levels in many scenarios").
   - A shared protein takes the sum of its genes' values before the log.
   - A complex, if complexes get descriptors, takes the minimum over its members, as reactions do now under "and" rules
     (machado2014systematic).
   - The documents will not use the term "pseudogene", which names a non-functional gene copy (UniProt lists RPL9P7 to
     RPL9P9 among the genes of the RPL9 entry).
   - If the block were meant to stay on the gene node, knockout seeds would again carry 27 columns that tell genes apart,
     and seed masking would matter again for knockouts.
6. **Seed masking.**
   - With the block on the protein node, a knockout's seed carries no descriptors. Its encoded protein, one hop away,
     carries them.
   - A drug seeds a protein node, so its seed carries descriptors. The user's statement holds for knockouts at the seed
     only.
   - The slice arms can measure the remaining shortcut before the confirmatory configurations change.
7. **Degree** (the user's decision: degree follows the split). Every use of degree reads the split graph's own degree:
   - the strata;
   - degree_popularity;
   - the label permutation within strata;
   - the log-degree node feature.

   This is also what the scorer does by default, since it reads degree from the graph directory of the runs.
8. **Rewiring for H2 holds the encodes relation fixed** (for the user's decision, since it changes the H2 rewiring).
   Rewiring swaps edges within each relation. Swapping encodes edges keeps every degree but joins gene X to protein Y,
   so a knockout would read Y's descriptors, and H2 would then also measure scrambled identity.
9. **Message passing gets one more layer,** since a knockout reaches the protein layer one hop later.

## Degree strata on the split (measured)

On the split, a knockout seeds its gene node, whose degree is its transcription in-degree plus the encodes edge. A drug
seeds protein nodes, whose degree is the merged degree minus the transcription in-degree plus 1.

Five quantile strata over the 1,539 perturbations (1,397 knockouts, 142 drugs):

| stratum | merged: drugs | merged: knockouts | split: drugs | split: knockouts |
|---|---|---|---|---|
| lowest | 1 | 279 | (empty: no degree 0) | |
| 2 | 1 | 330 | 0 | 564 |
| 3 | 8 | 301 | 0 | 354 |
| 4 | 34 | 274 | 2 | 289 |
| top | 98 | 213 | 140 | 190 |

Quantile edges: merged 4, 11, 25 and 57; split 1, 2, 4 and 11. The split's ties leave four nonempty strata.

Degree as a ranker for every symptom on development perturbations (macro and micro AUPRC as the scorer computes them):

| perturbations | random ranking (macro) | merged degree | split seed-node degree | split protein-node degree |
|---|---|---|---|---|
| development, all (1,222; 20 symptoms) | 0.087 | 0.138 and 0.115 | 0.173 and 0.125 | 0.137 and 0.114 |
| development knockouts (1,108; 16 symptoms) | 0.101 | 0.116 and 0.105 | 0.115 and 0.103 | 0.116 and 0.105 |

Readings:
- 40 percent of knockouts have split degree 1 (no transcription in-edge), and 140 of 142 drugs fall in the top
  stratum. Inside knockouts, split and merged degree correlate at 0.57 (Spearman), and the protein node's degree
  correlates with the merged degree at 0.97.
- The split seed degree ranks better than the merged degree over all perturbations only through the type offset.
  Drugs have more edges and a higher base rate, and inside knockouts all three definitions score the same.

Consequences of the user's decision:
- **The pooled readings get a stronger baseline.** degree_popularity now carries part of the type offset, which is an
  open item of the preregistration. This makes H1 harder, not easier.
- **The type offset earns no credit inside the knockout-only strata.** Strata 2 to 4 hold knockouts only (two drugs
  aside). The top stratum mixes 140 drugs with 190 knockouts, so the offset still earns credit there.
- **The strata no longer order knockouts by their product's connectivity** (design section 6.2: no credit for ordering
  perturbations by degree; on annotation bias, haynes2018gene). A model can still read the protein node's degree one hop
  from the seed. On these data that degree ranks knockouts little better than chance (0.116 against 0.101), so little
  credit is open to it.

## Not changed

- graph_full_neuronal, the slice graphs and the descriptors are untouched. module_fix_slice, confirmatory_pilots and
  refit_pilots read them.
- The confirmatory_v2 configurations still name graph_full_neuronal. Moving them to the split graph is the user's
  decision, after the slice check of item 6.
