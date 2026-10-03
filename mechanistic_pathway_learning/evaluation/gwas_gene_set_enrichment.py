"""Orthogonal validation with symptom-level GWAS (design section 6.4).

For each symptom: take the top-N genes by model attribution, drop any gene in
that symptom's training evidence, write a MAGMA gene-set file and run the
competitive gene-set test against the symptom's GWAS summary statistics
(PHQ-9 items; anxiety and depressive symptom factors). Compare with size-matched
random gene sets and with gene sets from the knowledge-graph embedding baseline.

MAGMA is an external binary (https://cncr.nl/research/magma/). This module only
prepares inputs and parses outputs. Not implemented in version 0.1.
"""


def write_magma_gene_set_file(gene_sets_by_symptom, output_path):
    raise NotImplementedError("one line per set: <set name> <gene ids...>")


def parse_magma_gene_set_results(results_path):
    raise NotImplementedError("read the .gsa.out table; return beta, standard error and p per set")
