# Re-deriving a drug's targets from the Guide to PHARMACOLOGY

Measured by `experiments/measure_gtopdb_target_rederivation.py`. Nothing is written by that script; this is the cost and reach of the change the owner agreed to on 13 October 2026, before any rebuild.

## Reach

- 222 drugs in the graph; 171 matched a GtoPdb ligand (39 by name, 132 by pubchem_cid).
- 114 have at least one human row naming a single gene with an action the vocabulary covers: 474 such rows of 505 human rows in all.
- Mechanism target entries over every drug: 1058 under ChEMBL, 1212 with GtoPdb's rows replacing ChEMBL's where it covers the drug.

## The drugs this narrows most

| drug | matched by | ChEMBL targets | GtoPdb targets | examples of what goes |
| --- | --- | --- | --- | --- |
| clonazepam | pubchem_cid | 13 | 3 | GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3 |
| diazepam | pubchem_cid | 13 | 3 | GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3 |
| flunitrazepam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |
| triazolam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |
| alprazolam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |
| pramipexole | name | 3 | 2 | DRD4 |
| CID100005095 | pubchem_cid | 3 | 2 | DRD4 |
| milnacipran | name | 2 | 2 | - |
| diphenhydramine | pubchem_cid | 1 | 1 | - |
| citalopram | pubchem_cid | 1 | 1 | - |
| sertraline | name | 1 | 1 | - |
| rivastigmine | name | 2 | 2 | - |
| retigabine | pubchem_cid | 5 | 5 | - |
| maprotiline | pubchem_cid | 1 | 1 | - |
| naltrexone | name | 3 | 3 | - |

## The drugs this widens

| drug | ChEMBL targets | GtoPdb targets | examples of what is added |
| --- | --- | --- | --- |
| frovatriptan | 2 | 3 | HTR1A |
| trazodone | 3 | 4 | HTR2B |
| paroxetine | 1 | 2 | P2RX4 |
| acetaminophen | 4 | 5 | TRPV4 |
| cocaine | 10 | 11 | HTR3A |
| doxepin | 2 | 3 | SLC6A4 |
| CID100002170 | 7 | 8 | HTR6 |
| diflunisal | 2 | 3 | DHFR |
| CID100005665 | 1 | 2 | SLC32A1 |
| imipramine | 2 | 3 | CHRM2 |
| carbamazepine | 10 | 11 | FZD8 |
| CID100006018 | 1 | 2 | SLC18A1 |
| CID100023897 | 2 | 3 | HRH1 |
| lorazepam | 13 | 14 | GPR68 |
| CID100010413 | 2 | 3 | HCAR1 |

## Rows that modulate an endogenous ligand rather than supplying one

These are what preregistration item 24's modulator relation reads: a partial agonist, a positive or negative allosteric modulation, or a potentiation.

- `CID100002118`: GABRA1, GABRA2, GABRA3, GABRA5
- `CID100002215`: ADRA2A, DRD2, DRD3, DRD4, DRD5, HTR1A
- `CID100002477`: HTR1A
- `CID100002487`: OPRK1, OPRM1
- `CID100002554`: FZD8
- `CID100002802`: GABRA1, GABRA2, GABRA3
- `CID100002803`: ADRA2A
- `CID100003016`: GABRA1, GABRA2, GABRA3
- `CID100003345`: OPRK1
- `CID100003380`: GABRA1
- `CID100003958`: GPR68
- `CID100004253`: OPRK1
- `CID100004419`: OPRD1, OPRM1
- `CID100004422`: OPRD1
- `CID100004440`: HTR1B
- `CID100005070`: KCNK4
- `CID100005078`: HTR1B
- `CID100005358`: HTR1B
- `CID100005556`: GABRA1, GABRA2, GABRA3, GABRA5
- `CID100005731`: HTR1B
- `CID100005732`: GABRA1, GABRA2, GABRA3
- `CID100060795`: DRD2, HTR1A, HTR2A, HTR2C
- `CID100060853`: HTR1A
- `CID106918313`: HTR1A
- `CID109966051`: HTR1A, HTR1B
