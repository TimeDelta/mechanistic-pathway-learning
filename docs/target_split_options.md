# Splitting the held-out set on mechanism targets instead of drugs

Measured by `experiments/label_review/target_split_options.py` over 200 full redraws at a held-out fraction of 0.143, on 222 drugs and 192 distinct targets, with the GtoPdb narrowing and every addition kept.

`fully covered` is the share of the held-out side whose mechanisms all appear on the training side: the quantity that condemned the drug split, where it is 84.5 percent. `unseen mechanism` is its complement, and is the share that supports a generalisation claim.

| rule | held-out share | spread | fully covered by training mechanisms | has a mechanism training never saw |
| --- | --- | --- | --- | --- |
| drug (as drawn today) | 14.3% | 0.000 | 84.4% | 15.6% |
| all_targets_held | 3.1% | 0.021 | 81.1% | 18.9% |
| primary_target_held | 7.1% | 0.032 | 75.8% | 24.2% |
| any_target_held | 46.0% | 0.086 | 0.0% | 100.0% |
| mechanism_edge | 13.8% | 0.034 | 0.0% | 100.0% |
