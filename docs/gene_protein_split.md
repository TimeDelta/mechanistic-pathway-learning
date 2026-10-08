# Gene and protein nodes: the split and the gene-to-protein mapping

Written on 8 October 2026 and revised the same day after the user's corrections (gene descriptors limited to the brain
expression columns; one protein node per entry where one gene has several; seed masking included). The split is built
into new directories (section "Built"); no existing graph directory was regenerated, and no confirmatory configuration
uses the split. The measurements read graph_full_neuronal, evidence_full_v2 with better_v2 and the pinned reviewed
UniProt table (20,431 entries). The degree rankers below read development outcomes only (the 317 lockbox perturbations
are left out); lockbox counts are membership counts.

## The user's statements

- "The genes shouldn't need descriptors but the proteins should."
- "I'm thinking for number two that the seed [masking] isn't necessary for the message passing after gene descriptors
  are removed."
- "But doesn't that assume a 1:1 mapping between gene and protein?"
- "If the gene descriptors are dropped, I would specify that the brain expression is still provided via pseudo-gene
  descriptors. The rest of that suggested split between the genes and proteins makes sense I think except the
  registered degree strata SHOULD change to the new split."
- "No. That's not what I meant by pseudo gene descriptors. I [guess] I just shouldn't say pseudo. It's just that the
  descriptors would be limited to only the brain expression columns. Also, the several genes [to one] protein change is
  fine but I would also require an equivalent several proteins from one gene. Seed [masking] would then be included."

## Is the mapping one to one?

No. The merged gene node of the current graphs and the split as first proposed (one PROTEIN:X per GENE:X) both assume
it. Departures among the 12,810 gene nodes:

| case | gene nodes | examples | perturbation seeds (membership only) |
|---|---|---|---|
| one reviewed entry encoded by several genes | 57 | histone H2A, H2B and H3.3 genes (H3-3A, H3-3B), defensins, CKMT1A and CKMT1B | 1 development knockout |
| one gene with several reviewed entries (distinct products) | 37 (84 entries) | NRXN1, NRXN2 and NRXN3 (alpha and beta neurexins), GNAS, CDKN2A (p16INK4a and p14ARF), CALCA (calcitonin and CGRP), PRNP, TSPO | 10 knockouts, 1 drug |
| isoforms inside an entry | not in the pinned table | most multi-exon genes | all |
| one entry cut into several peptides | not counted | POMC: edges to MC1R to MC5R (the melanocortins) and to OPRD1 and OPRK1 (beta-endorphin) on one node; PDYN, PENK | not counted |
| no reviewed protein | 113 | 19 ChEBI compounds typed as gene; 58 unreviewed UniProt accessions used as symbols; 36 outdated or malformed symbols (GBA for GBA1; nodes named "HBA1; HBA2", "CCL3L1; CCL3L3" and "CCL4L1; CCL4L2") | none |

The 113 nodes without a reviewed protein take part in 2,170 edges (activates 1,141, inhibits 713, binds 254,
catalyzed_by 40, regulates_transcription_of 14, member_of 8).

Isoforms are the largest departure. pan2008deep: "transcripts from approximately 95% of multiexon genes undergo
alternative splicing". yang2016widespread: "The majority of isoform pairs share less than 50% of their interactions" and
"alternative isoforms tend to behave like distinct proteins rather than minor variants of each other".

Several reviewed entries of one gene are distinct proteins, not isoforms. quelle1995alternative, on CDKN2A: "An
unrelated protein (p19ARF) arises in major part from an alternative reading frame of the mouse INK4a gene". Many of the
37 genes' second entries are short peptides from an upstream or alternative open reading frame (DDIT3: 34 residues;
PRKCH: 26; PRNP: 73), or polyproteins (ERVK-6: five entries).

What the sources resolve, for the 1,395 edge ends at these 37 genes in graph_full_neuronal (an end is the source or
target of an edge, the target end of transcription excluded since it stays at the gene):
- 1,082 ends come from rows that name every entry of the gene. OmniPath and CollecTRI copy a gene's interactions onto
  each of its entries: the partner lists are identical for CDKN2A's two entries (80 rows each), NRXN1's two (45 each) and
  CALCA's two (44 each), so p16INK4a's partners CDK4 and CDK6 and p14ARF's partner MDM2 sit on both entries.
- 69 ends name a subset: 40 Reactome entity memberships (CALCA, CDKN2A, GNAS, MT-ND4, MT-RNR2, PRKCH and SEM1) and
  29 OmniPath ends. CDKN2A's two Reactome memberships name p14ARF only; GNAS's 14 name Gs alpha and XLas, not NESP55
  or ALEX.
- 244 ends name none of the gene's reviewed entries (mostly OmniPath rows whose accession is not one of them) or come
  from a source that does not resolve entries (Human-GEM gene rules, 5 ends; OmniPath small molecules, 2).
- The OmniPath interaction file holds no isoform accession (0 of 126,125 rows on either side).
- ChEMBL targets name one entry for each of the five such genes they reach: CGRP (P06881) for CALCA, P24723 for PRKCH,
  P03905 for MT-ND4, Q9Y3U8 for RPL36 and P60896 for SEM1. No gene's targets name different entries.

Defects found:
- experiments/build_protein_descriptors.py groups entries by UniProt's primary gene-name field. 68 entries list several
  genes ("DEFB4A; DEFB4B"), so 57 gene nodes get no protein descriptors in any current graph. It is not fixed under the
  running jobs; the split build fixes it.
- The sign of each reduced-rank direction is not fixed by the fit: a refit with another thread count reproduced the 64
  descriptor columns to 3.7e-12 up to sign, with 29 columns flipped. A learned input layer absorbs a flip, but a rerun
  of build_protein_descriptors.py would not reproduce the stored table byte for byte. The per-entry table takes the
  stored signs.

## Handling in the split

1. **Protein nodes per reviewed entry** (the user's decisions). A gene listed by one entry gets one protein node. Genes
   listed by one shared entry share its protein node, with an encodes edge from each gene (H3-3A and H3-3B to H3.3); a
   knockout of one leaves the other's route to the protein. A gene listed by several entries gets one protein node per
   entry, each with its own encodes edge (CDKN2A to PROTEIN:CDKN2A:P42771 and PROTEIN:CDKN2A:Q8N726). Isoforms inside an
   entry and peptides cut from one precursor are not represented, since no source resolves them.
2. **Edges of a gene with several entries** (default, for the user). An edge end goes to the entries its source names
   (OmniPath and CollecTRI accession columns, Reactome entity members); where the source names none of the gene's
   entries or does not resolve entries, it goes to every entry, since which product is meant is not known. A drug acting
   on such a gene seeds the entries its ChEMBL targets name, or every entry when they name none.
   - Consequence: the sources copy most edges (1,082 of 1,395 ends), so the protein nodes of one gene mostly differ in
     their protein descriptors, not in their edges. A knockout of such a gene reaches every copy, which counts a copied
     edge once per entry (twice for CDKN2A, four times for GNAS). Both encoders average a node's incoming messages per
     relation (typed mean in message passing, in-degree normalisation in linear response), and the copies count in the
     denominator too, so a neighbour's input changes by a factor below the copy count. The 10 knockouts and 1 drug on these genes are 0.7 percent of the perturbations.
   - The alternative, sending unresolved ends to one main entry, needs a rule for which entry is the main one; UniProt
     and HGNC give none across entries.
3. **The 113 nodes without a reviewed protein.** The three named after a whole gene-names field ("HBA1; HBA2",
   "CCL3L1; CCL3L3", "CCL4L1; CCL4L2") stand for that entry and join its protein node. The other 110 get one protein
   node each, without protein descriptors. Retyping the 19 ChEBI compounds as metabolites and remapping the unreviewed
   accessions and outdated symbols to current HGNC symbols remain open for the user.
4. **Leakage groups.** Groups are built from the evidence's own gene nodes, so folds and lockbox membership stay as
   drawn (a test and the full-graph check below confirm it).
5. **Descriptors** (the user's decision: the gene's descriptors "would be limited to only the brain expression
   columns"). Gene nodes carry the 27-column brain expression block (GTEx v10, HPA brain regions and HPA single-nucleus
   cell classes, with the dopaminergic class) and the expression columns of nodes.parquet (log brain TPM, expressed
   flag). Protein nodes carry the protein block, each its own entry's row (experiments/build_protein_descriptors.py
   --per-entry-only), and no expression. Reactions keep their own blocks, including the expression carried onto them by
   the gene rule.
   - Expression is mRNA measured per gene, so it sits where it is measured. liu2016dependency: "transcript levels by
     themselves are not sufficient to predict protein levels in many scenarios"; a protein node reads it one encodes
     edge away.
   - Cell-class weights of the linear-response encoder (where a node's signal propagates) stay on the gene and are
     copied to each of its protein nodes (the largest over the genes of a shared node): a protein acts in the cells that
     make it.
6. **Seed masking** (the user: "Seed [masking] would then be included"). With seed_masked, each perturbation's seed
   nodes read their structural columns only, as on the merged graph, and on a split graph the mask also covers the nodes
   joined to a seed by an encodes edge (default, for the user). Without it the knocked-out gene's protein descriptors, or
   a drug target's gene expression, would sit one hop from the seed, where the merged graph hid them on the seed itself.
   It applies only on graphs with an encodes relation, so no merged-graph run changes. The slice arms carry seed masking
   for message passing; the linear-response arms stay without a treatment, since the descriptor rule's pick there
   (zero_init_slow) still waits for the user.
7. **Degree** (the user's decision: degree follows the split). Every use of degree reads the split graph's own degree:
   - the strata;
   - degree_popularity;
   - the label permutation within strata;
   - the log-degree node feature.

   This is also what the scorer does by default, since it reads degree from the graph directory of the runs.
8. **Rewiring for H2 holds the encodes relation fixed** (open for the user, since it changes the H2 rewiring). Rewiring
   swaps edges within each relation. Swapping encodes edges keeps every degree but joins gene X to protein Y, so a
   knockout would reach Y's descriptors, and H2 would then also measure scrambled identity. --rewire-encodes swaps them.
9. **Message passing gets one more layer,** since a knockout reaches the protein layer one hop later.

## Degree strata on the split (measured)

On the split, a knockout seeds its gene node, whose degree is its transcription in-degree plus its encodes edges (one
per entry). A drug seeds protein nodes, whose degree is the merged degree minus the transcription in-degree plus the
encodes edges.

Five quantile strata over the 1,539 perturbations (1,397 knockouts, 142 drugs):

| stratum | merged: drugs | merged: knockouts | split: drugs | split: knockouts |
|---|---|---|---|---|
| lowest | 1 | 279 | (empty: no degree 0) | |
| 2 | 1 | 330 | 0 | 561 |
| 3 | 8 | 301 | 0 | 354 |
| 4 | 34 | 274 | 2 | 291 |
| top | 98 | 213 | 140 | 191 |

Quantile edges: merged 4, 11, 25 and 57; split 1, 2, 4 and 11. The split's ties leave four nonempty strata.

Degree as a ranker for every symptom on development perturbations (macro and micro AUPRC as the scorer computes them),
measured on the first build (one protein node per gene); the per-entry nodes change the seed degree of 10 knockouts:

| perturbations | random ranking (macro) | merged degree | split seed-node degree | split protein-node degree |
|---|---|---|---|---|
| development, all (1,222; 20 symptoms) | 0.087 | 0.138 and 0.115 | 0.173 and 0.125 | 0.137 and 0.114 |
| development knockouts (1,108; 16 symptoms) | 0.101 | 0.116 and 0.105 | 0.115 and 0.103 | 0.116 and 0.105 |

Readings:
- 40 percent of knockouts have split degree 1 (no transcription in-edge, one entry), and 140 of 142 drugs fall in the
  top stratum. Inside knockouts, split and merged degree correlate at 0.57 (Spearman), and the protein node's degree
  correlates with the merged degree at 0.97.
- The split seed degree ranks better than the merged degree over all perturbations only through the type offset.
  Drugs have more edges and a higher base rate, and inside knockouts all three definitions score the same.

Consequences of the user's decision:
- **The pooled readings get a stronger baseline.** degree_popularity now carries part of the type offset, which is an
  open item of the preregistration. This makes H1 harder, not easier.
- **The type offset earns no credit inside the knockout-only strata.** Strata 2 to 4 hold knockouts only (two drugs
  aside). The top stratum mixes 140 drugs with 191 knockouts, so the offset still earns credit there.
- **The strata no longer order knockouts by their product's connectivity** (design section 6.2: no credit for ordering
  perturbations by degree; on annotation bias, haynes2018gene). A model can still read the protein node's degree one hop
  from the seed. On these data that degree ranks knockouts little better than chance (0.116 against 0.101), so little
  credit is open to it.

## Not changed

- graph_full_neuronal, the slice graphs and the descriptors are untouched. module_fix_slice, confirmatory_pilots and
  refit_pilots read them.
- The confirmatory_v2 configurations still name graph_full_neuronal. Moving them to the split graph is the user's
  decision, after the slice runs below.

## Built (8 October 2026, new directories only)

experiments/build_gene_protein_split.py (library: mechanistic_pathway_learning/graph/gene_protein_split.py) wrote:
- data/processed/graph_split, from the slice graph:
  - 2,848 gene nodes and 2,842 protein nodes (4 shared by 13 genes; MOCS2, MT-ND4 and POLR1D with two entries each;
    9 without a reviewed entry);
  - 97,632 edges; 2,845 gene nodes of degree 1, because the slice has no transcription edges, so on the slice almost all
    knockouts share one degree stratum.
- data/processed/graph_full_neuronal_split:
  - 12,810 gene nodes and 12,834 protein nodes: 18 shared by 41 genes (H3-3A/H3-3B, SMN1/SMN2, CKMT1A/CKMT1B, the
    USP17L cluster, HBA1 with the node named "HBA1; HBA2" and others), 84 for the 37 genes with several entries, 110
    without a reviewed entry;
  - 281,496 edges (266,966 before, 1,683 more from ends sent to several entries, 10 merged onto shared nodes, plus
    12,857 encodes edges); 6,578 gene nodes of degree 1.
- The descriptor tables and cell-class weight tables of both, under data/processed/node_descriptors/*_split_* and
  data/processed/cell_class_weights/*_split_*:
  - gene rows carry no protein column, protein rows no other column, and every other column is unchanged on the nodes
    that existed before;
  - on the full graph, 12,724 protein nodes (every one with a reviewed entry) carry protein descriptors, against 12,643
    gene nodes before; on one-to-one nodes the protein row equals the gene's earlier row, except 19 that had none (the
    multi-gene entry defect);
  - gene nodes keep their cell-class weights, and protein nodes take their gene's.

Code changes, all inactive on a merged graph:
- **Loader.** experiment_data.load_experiment_data seeds a drug's target proteins when the graph directory holds
  gene_to_protein.parquet, and builds leakage groups from the evidence's own gene nodes. On the slice and on
  graph_full_neuronal with evidence_full_v2 and better_v2, the edited loader returned the same seeds, signs, magnitudes,
  outcomes, masks and groups as before (checked on the first build).
- **Checks on graph_full_neuronal_split.** The perturbations, outcomes and groups are identical to the merged graph's;
  all 142 drugs seed protein nodes and all 1,397 knockouts gene nodes, with the same seed counts; the one drug on a gene
  with several entries (RXCUI:2282660, targets CALCA and CALCB) seeds PROTEIN:CALCA:P06881 (CGRP) and PROTEIN:CALCB.
- **Seed masking.** descriptor_treatments.seed_unit_mask widens the mask over encodes partners; run_main_model passes the
  pairs only when the graph has an encodes relation and the treatment is seed_masked. The partner index is not saved in
  checkpoints, so a checkpoint loads either way.
- **Rewiring.** encodes is held fixed unless --rewire-encodes. It comes last in index order, so the other relations'
  random stream is unchanged.
- **Rewired runs.** The rewired reaction rules read a protein catalyst's Ensembl id, which the protein node keeps as an
  identifier, not as a descriptor.
- **Tests and smoke runs.**
  - tests/test_gene_protein_split.py (6 tests); the full suite passes.
  - Two-epoch smoke runs of both split slice arms trained without error; the message-passing run reports the mask over
    2,851 encodes edges.

## Slice runs (job split_slice, after module_fix_slice)

Five disease-cluster folds, seed 0, each arm with --start-at-weighted-optimum and the slice brain expression descriptors
(experiments/run_main_model_batch.py):
- message passing with seed masking on the merged slice (2 layers) against the split slice (3 layers, so knockouts reach
  as far into the protein layer);
- linear response on the merged slice against the split slice, without a descriptor treatment.

The arm without seed masking on the split slice is dropped: it tested the user's earlier statement that masking may not
be needed once the genes carry no descriptors, and the genes carry the expression block again.

The readings, fixed before the runs, are paired five-fold differences in macro and micro AUPRC
(experiments/compare_twin_runs.py) and the module-health readings of docs/module_health.md
(experiments/module_health.py). The slice is metabolic: its knockouts reach reactions through one encodes edge and no
protein interactions, so it tests the descriptor placement and the extra hop, not the protein layer of the full graph.
They are development readings and go to the user, who decides whether the confirmatory configurations move to the split.

The merged and split arms differ in three arguments (graph, layer count and descriptor table), so their difference
reads the split as a whole, not one change. The within-degree-strata reading is taken twice, with strata from the split
slice (the user's decision for the registered strata) and from the merged slice, each with its own cache, because
compare_twin_runs.py keys its cache by run names only:

    pairs="b6_mechanistic_gate_time_scales_weighted_start_descriptors_split_seed_masked:b6_mechanistic_gate_time_scales_weighted_start_descriptors_seed_masked \
      b6_linear_response_gate_time_scales_cofactors_weighted_start_descriptors_split:b6_linear_response_gate_time_scales_cofactors_weighted_start_descriptors \
      b6_mechanistic_gate_time_scales_weighted_start_descriptors_seed_masked:b6_mechanistic_gate_time_scales_weighted_start \
      b6_linear_response_gate_time_scales_cofactors_weighted_start_descriptors:b6_linear_response_gate_time_scales_cofactors_weighted_start"
    OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/compare_twin_runs.py --graph-dir data/processed/graph_split \
      --run-root runs/module_fix --cache-dir runs/module_fix/twin_comparisons_split_strata --pairs $pairs \
      --markdown-output docs/split_slice_comparisons.md --json-output runs/module_fix/split_slice_comparisons.json
    OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/compare_twin_runs.py --graph-dir data/processed/graph \
      --run-root runs/module_fix --cache-dir runs/module_fix/twin_comparisons_merged_strata --pairs $pairs \
      --markdown-output docs/split_slice_comparisons_merged_strata.md --json-output runs/module_fix/split_slice_comparisons_merged_strata.json

The last two pairs read the brain expression descriptors (with seed masking for message passing) on the merged slice
under the module fix.

## Correction of 8 October 2026 (late): complex accessions

OmniPath sorts a complex's accessions but not its gene symbols: the CDKN2A–MDM2 complex reads `CDKN2A_MDM2` with
`COMPLEX:Q00987_Q8N726`, MDM2's accession first. The split paired them by position, so CDKN2A took MDM2's accession
from that row and only P42771 (p16INK4a) from its twin row `COMPLEX:P42771_Q00987`. That one record was the only thing
telling CDKN2A's entries apart. OmniPath copies every CDKN2A row, the CDK4 and CDK6 complex rows included, to both
P42771 and Q8N726. interaction_members now gives each symbol the row's accessions that are its own entries.

Rebuilt graph_full_neuronal_split: CDKN2A collapses to one protein node; CALCA (P01258 calcitonin, P06881 CGRP) is the
only gene with several protein nodes. 12,709 protein nodes, 279,757 edges, 36 protein nodes of several entries.
graph_split and its tables come out identical, so split_slice was not touched.
