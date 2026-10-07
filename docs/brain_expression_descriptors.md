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
(doi:10.1371/journal.pcbi.1003580). Each block is standardised over its node type; metabolites get zeros. 26 columns
per block, 52 in all, on top of the 83 existing descriptor columns.

Coverage on the slice graph: 2,843 of 2,848 gene nodes are in the
HPA region table; 7,772 of 12,877 reaction nodes have a known gene in their rule
(7,776 have a rule at all; the rest are transport, exchange and spontaneous reactions).

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
class ("other neuron", the splatter and miscellaneous superclusters) rests on under 1 nCPM, and this block cannot
represent dopamine neurons as a class. The GTEx columns still carry SLC6A3 (30.7 TPM in brain, from substantia nigra).

## What it can and cannot test

Two slice arms, each one argument from its twin:

- b3_linear_response_cofactors_descriptors against b3_linear_response_cofactors: the existing protein, metabolite and
  reaction descriptors in the linear-response encoder's output gate (field = (h W) * sigmoid(node_features W_g + b)).
  This has never been run: the only descriptor runs used the message-passing encoder (b3_typed_nodes_descriptors) and
  the descriptors-only control. So the slice has not yet tested the node properties inside the mechanistic encoder.
- b3_linear_response_cofactors_descriptors_brain against b3_linear_response_cofactors_descriptors: the brain block
  added.

Expression in the output gate reweights where the response is read; it does not change how the response propagates.
A cell-class layer that changes propagation (graph copies per cell class joined by transport, as in Lewis et al. 2010,
doi:10.1038/nbt.1711) is separate work, held for the user's decision on building a pared graph derived from the full
graph. Machado and Herrgård found that "for many conditions, the predictions obtained by simple flux balance analysis
using growth maximization and parsimony criteria are as good or better than those obtained using methods that
incorporate transcriptomic data", so a null here would not be surprising and would not by itself say the regional or
cell-type information is useless; the slice labels are mostly metabolic-disease phenotypes reached in few steps.
