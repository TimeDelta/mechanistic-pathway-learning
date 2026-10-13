# What the GtoPdb re-derivation costs in leakage

Measured by `experiments/label_review/gtopdb_leakage_cost.py` on `data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes` with `data/processed/evidence_full_v4`. Measurement only; nothing is written by it.

| arm | drugs | kept positive pairs | drug-holding groups | two largest blocks | largest block, share of drugs | share of kept pairs |
| --- | --- | --- | --- | --- | --- | --- |
| as built (ChEMBL) | 222 | 2975 | 17 | [114, 71] | 0.514 | 0.272 |
| re-derived, narrowing only (no GtoPdb target outside a ChEMBL row) | 222 | 2975 | 17 | [114, 71] | 0.514 | 0.272 |
| re-derived, narrowing and brain-expressed additions only | 222 | 2975 | 14 | [197, 5] | 0.887 | 0.482 |
| re-derived, narrowing and every addition | 222 | 2975 | 13 | [200, 5] | 0.901 | 0.484 |
| re-derived, narrowing on one subunit and every addition (the rejected rule) | 222 | 2975 | 13 | [200, 5] | 0.901 | 0.484 |
