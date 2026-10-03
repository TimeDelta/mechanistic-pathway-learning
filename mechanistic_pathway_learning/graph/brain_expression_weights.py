"""Brain expression weights per gene (design section 4.1, node attributes).

Version 1: a scalar weight in [0, 1] from Human Protein Atlas brain tissue
consensus or GTEx brain tissues (median TPM, rank-normalized). Version 2:
cell-type-specific weights. Weights are node attributes used by the encoder,
not filters, so peripheral genes stay in the background graph.

Not implemented in version 0.1.
"""


def compute_brain_expression_weights(expression_table_path):
    raise NotImplementedError("rank-normalize brain median expression per gene to [0, 1]")
