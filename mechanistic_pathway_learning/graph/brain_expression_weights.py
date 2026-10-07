"""Brain expression weights per gene (design section 4.1, node attributes): an unused stub.

What is implemented, and where: mechanistic_pathway_learning/graph/build_physiology_graph.py reads the GTEx v10
median-TPM table (load_brain_expression) and gives every gene node the largest median TPM over the 13 GTEx brain
tissues and a brain_expressed flag (annotate_brain_expression); reactions take the largest value over their catalysts
(propagate_brain_expression_to_reactions), and regulon edges of factors known not to be brain expressed are dropped
(restrict_transcription_edges_to_brain_expressed). load_experiment_data turns the value into the log1p brain-expression
column of the structural node features and the flag into the next column. These are node attributes, not filters,
so peripheral genes stay in the graph.

Expression per brain region (13 Human Protein Atlas regions) and per brain cell class (10 classes grouped from 34
single-nucleus cluster types) is implemented as node descriptors, not as structural features, in
mechanistic_pathway_learning/graph/brain_expression_descriptors.py (built by
experiments/build_brain_expression_descriptors.py); unlike the GTEx column it is also available on the slice graphs.

What is not implemented: the rank-normalised weight in [0, 1] that compute_brain_expression_weights below was meant
to return (nothing calls it), and cell-class layer copies of the graph (nodes per cell class, as in Lewis et al. 2010).
"""


def compute_brain_expression_weights(expression_table_path):
    raise NotImplementedError("not implemented and not called; see the module docstring for the implemented brain expression")
