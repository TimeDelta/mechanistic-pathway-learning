# Brain region and cell-class expression as node descriptors (7 October 2026)

The slice graphs carry no brain expression: only graph_full and graph_full_neuronal have the GTEx column, and no graph
had regional or cell-type expression. Following the user's direction of 7 October (brain region and cell type are
central to psychiatric symptoms; add them without adding noise), they now enter as node descriptors, which add columns
and no nodes or edges. Code: mechanistic_pathway_learning/graph/brain_expression_descriptors.py; table:
experiments/build_brain_expression_descriptors.py, written to data/processed/node_descriptors/ (gitignored), outside
every graph directory, so no running job's graph changes.

## What the block holds

Per gene, log1p of: the largest GTEx v10 median TPM over the 13 brain tissues and over every other tissue; Human Protein
Atlas consensus nTPM in 13 brain regions; and Human Protein Atlas single-nucleus nCPM, the largest over the cluster
types of each of 10 cell classes (excitatory neurons, cortical interneurons, medium spiny neurons, other inhibitory
neurons, other neurons, astrocytes, oligodendrocyte lineage, immune, vascular, ependymal and choroid). The 34 cluster
types come from the adult human brain atlas of Siletti et al. 2023 (doi:10.1126/science.add7046) and are grouped to
keep the column count low for 451 training perturbations. Reactions take each value from their Human-GEM gene rule,
minimum over enzyme complexes and maximum over isozymes, the mapping Machado and Herrgård 2014 report as general
(doi:10.1371/journal.pcbi.1003580). Each block is standardised over its node type; metabolites get zeros. With the
dopaminergic class of the next section, 27 columns per block, 54 in all, on top of the 83 existing descriptor columns
(26 and 52 before it was added).

Coverage on the slice graph: 2,843 of 2,848 gene nodes are in the
HPA region table; 7,772 of 12,877 reaction nodes have a known gene in their rule
(7,776 have a rule at all; the rest are transport, exchange and spontaneous reactions).

Before the dopaminergic class was added:

| gene | GTEx brain max TPM | GTEx other max TPM | top HPA region | top HPA cell class |
|---|---|---|---|---|
| PAH | 2.7 | 272.4 | medulla oblongata | other inhibitory neuron |
| OTC | 0.5 | 49.6 | white matter | astrocyte |
| CPS1 | 8.3 | 312.2 | basal ganglia | oligodendrocyte lineage |
| GCH1 | 8.7 | 58.0 | midbrain | immune |
| TH | 36.9 | 5.3 | pons | medium spiny neuron |
| DDC | 2.4 | 55.3 | midbrain | other neuron |
| DBH | 2.5 | 23.4 | pons | other neuron |
| MAOA | 19.2 | 168.8 | pons | vascular |
| COMT | 39.8 | 52.7 | white matter | ependymal choroid |
| SLC6A3 | 30.7 | 1.8 | midbrain | other neuron |
| SLC18A2 | 1.2 | 49.5 | midbrain | other inhibitory neuron |
| GAD1 | 60.2 | 8.1 | hypothalamus | cortical interneuron |
| ALDH5A1 | 44.1 | 30.2 | midbrain | other inhibitory neuron |
| ATP1A3 | 346.8 | 89.9 | cerebral cortex | other inhibitory neuron |
| SLC12A5 | 206.6 | 5.1 | cerebral cortex | other inhibitory neuron |

The peripheral anchors read as peripheral (PAH 2.7 against 272 TPM, OTC 0.5 against 50, CPS1 8.3 against 312), which
the slice had no way to see before. The top cell class is not interpretable for genes at a few nCPM: TH ranks medium
spiny neurons first. The atlas's 34 cluster types include no dopaminergic one and SLC6A3 peaks at 0.6 nCPM, so its top
class ("other neuron", the splatter and miscellaneous superclusters) rests on under 1 nCPM, so HPA's cluster types
cannot represent dopamine neurons as a class. The GTEx columns still carry SLC6A3 (30.7 TPM in brain, from substantia
nigra). The next section adds the class from the atlas itself.

## The dopaminergic class (7 October 2026)

The user required a dopaminergic class. HPA's finer brain table (rna_single_nuclei_brain_cluster, 260 clusters) is
brain region crossed with the same 34 cluster types, so it has none either (its highest SLC6A3 is 5.0 nCPM, in a
midbrain splatter cluster). The atlas itself has one: cluster 395, "DA VGLUT2" in its tables/cluster_annotation.xlsx
(github.com/linnarsson-lab/adult-human-brain), 998 nuclei from 3 donors, 92.6 percent from midbrain dissections, with
SLC6A3, TH, EN1, SLC18A2, LMX1B and LMX1A among its top genes. The atlas's per-dissection files on CZ CELLxGENE
(collection 283d65eb-dd53-496d-adb7-7570c7caa443, CC BY 4.0, raw UMI counts) hold the nuclei. Four dissections, SN-RN,
SN, PAG-DR and PAG (1.24 GB), hold 872 of the 998, from all 3 donors; experiments/build_dopaminergic_expression.py
reads only those rows, sums their counts per gene, keeps HPA's 19,580 genes and divides counts per million by a trimmed
mean of log ratios against HPA's splatter cluster type (the supercluster cluster 395 sits in; the method of Robinson
and Oshlack 2010, doi:10.1186/gb-2010-11-3-r25, without precision weights). The factor came out at 1.04, as expected
for a profile from the same atlas. The column is class_dopaminergic_neuron, an 11th cell class.

| gene | dopaminergic class nCPM | highest HPA cluster type nCPM (type) |
|---|---|---|
| TH | 388.0 | 6.4 (eccentric medium spiny neuron) |
| SLC6A3 | 302.7 | 0.6 (splatter) |
| SLC18A2 | 706.8 | 32.5 (cerebellar inhibitory) |
| DDC | 23.0 | 1.3 (splatter) |
| KCNJ6 | 2,079.7 | 677.3 (hippocampal dentate gyrus) |
| LMX1B | 90.1 | 14.6 (ependymal cell) |
| EN1 | 99.6 | 3.5 (oligodendrocyte precursor cell) |
| GAD1 | 2.6 | 554.5 (MGE interneuron) |
| AQP4 | 1.7 | 951.4 (astrocyte) |
| MBP | 111.7 | 5,605.9 (oligodendrocyte) |

The markers are where a dopamine neuron profile puts them, and the glial and GABAergic markers are low. MBP at 112 nCPM
is about 2 percent of its oligodendrocyte value, consistent with the ambient RNA that single-nucleus data carry. The
limits: 872 nuclei from 3 donors; Kamath et al. 2022 (doi:10.1038/s41593-022-01061-1) have about 15,700 dopamine
nuclei from 8 control donors and 10 subtypes, a larger source for subtypes if the single class proves too coarse, but
from a different study, so batch effects against the HPA columns would be expected. The DDC value (23 nCPM) is low for
the enzyme that makes dopamine, and is left as measured.

With the class in the tables (slice and full neuronal, both rebuilt on 7 October), TH, DDC, SLC6A3 and SLC18A2 rank it
first among the 11 classes; DBH, the noradrenergic enzyme, stays with other neurons. Coverage is unchanged, since the
column uses HPA's gene set.

## What it can and cannot test

Two slice arms, each one argument from its twin:

- b3_linear_response_cofactors_descriptors against b3_linear_response_cofactors: the existing protein, metabolite and
  reaction descriptors in the linear-response encoder's output gate (field = (h W) * sigmoid(node_features W_g + b)).
  This has never been run: the only descriptor runs used the message-passing encoder (b3_typed_nodes_descriptors) and
  the descriptors-only control. So the slice has not yet tested the node properties inside the mechanistic encoder.
- b3_linear_response_cofactors_descriptors_brain against b3_linear_response_cofactors_descriptors: the brain block
  added.

Expression in the output gate reweights where the response is read; it does not change how the response propagates.
A cell-class layer that changes propagation is separate work. Graph copies per cell class joined by transport, as in
Lewis et al. 2010 (doi:10.1038/nbt.1711), would multiply the graph by the number of classes; the user's direction of
7 October is instead one propagation channel per cell class, each with the same edges weighted by that class's
expression (docs/research_summary.md, decisions). Machado and Herrgård found that "for many conditions, the predictions obtained by simple flux balance analysis
using growth maximization and parsimony criteria are as good or better than those obtained using methods that
incorporate transcriptomic data", so a null here would not be surprising and would not by itself say the regional or
cell-type information is useless; the slice labels are mostly metabolic-disease phenotypes reached in few steps.

## The full neuronal graph (7 October 2026)

The same table for data/processed/graph_full_neuronal (data/processed/node_descriptors/
full_neuronal_descriptors_brain_expression.parquet, after experiments/build_node_descriptors.py wrote the 83 base
columns for that graph). Most of its gene nodes come from OmniPath and CollecTRI, which name genes by symbol only:
9,962 of 12,810 have no Ensembl id, and 9,783 of those are found through the Human Protein Atlas symbol table (symbols
that name more than one Ensembl id are not used). In all, 12,626 of 12,810 gene nodes and 7,752 of the 7,759 reactions
with a gene rule carry expression. The 1,681 Reactome protein entities (complexes and sets) and the membrane-potential
node get zeros; giving a complex the minimum over its members would follow the reaction rule, and is not done yet.
