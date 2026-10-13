# What a drug's own label says about acting centrally

Read from openFDA by `experiments/read_label_central_action.py`, pharmacology sections only (mechanism of action, clinical pharmacology, pharmacodynamics, description). The adverse-reaction, warning and indication sections are deliberately not read: a drug admitted for the psychiatric events on its label would select the drug set on the outcome the study predicts. A positive verdict is a shortlist to read, not a decision: class boilerplate makes every beta blocker positive, and a combination-product label can put a drug's only match on a sentence about the other ingredient.

2 drugs read.

| verdict | drugs |
| --- | --- |
| states a central action | 2 |

## states a central action

| drug | section | quote | openFDA set |
| --- | --- | --- | --- |
| metoprolol | clinical_pharmacology | Metoprolol is also known to cross the blood brain barrier following oral administration and CSF concentrations close to that observed in plasma have been reported. | 00940cc5-d2eb-4841-9138-de97d7b1c674 |
| propranolol | clinical_pharmacology | Factors that may contribute to the antihypertensive action include: (1) decreased cardiac output, (2) inhibition of renin release by the kidneys, and (3) diminution of tonic sympathetic nerve outflow from vasomotor centers in the brain. | 013a485a-4227-46ae-b19f-c550302fcb96 |
