# Re-deriving a drug's targets from the Guide to PHARMACOLOGY

Measured by `experiments/measure_gtopdb_target_rederivation.py`. Nothing is written by that script; this is the cost and reach of the change the owner agreed to on 13 October 2026, before any rebuild.

## Reach

- 222 drugs in the graph; 171 matched a GtoPdb ligand (39 by name, 132 by pubchem_cid).
- 124 have at least one human row naming a single gene with an action the vocabulary covers: 474 such rows of 505 human rows in all.
- Mechanism target entries over every drug: 1058 under ChEMBL, 1067 with GtoPdb's rows replacing ChEMBL's where it covers the drug.

## The drugs this narrows most

| drug | matched by | ChEMBL targets | GtoPdb targets | examples of what goes |
| --- | --- | --- | --- | --- |
| gabapentin | pubchem_cid | 22 | 2 | CACNA1A, CACNA1B, CACNA1C, CACNA1D, CACNA1E, CACNA1F |
| amisulpride | pubchem_cid | 19 | 1 | DRD3, HTR1A, HTR1B, HTR1D, HTR1E, HTR1F |
| lorazepam | pubchem_cid | 13 | 1 | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5, GABRA6 |
| propofol | pubchem_cid | 13 | 3 | GABRA1, GABRA2, GABRA3, GABRA4, GABRA5, GABRA6 |
| diazepam | pubchem_cid | 13 | 3 | GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3 |
| clonazepam | pubchem_cid | 13 | 3 | GABRA4, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3 |
| flunitrazepam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |
| ropivacaine | name | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN2A, SCN3A, SCN4A |
| cocaine | name | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN2A, SCN3A, SCN4A |
| triazolam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |
| lidocaine | pubchem_cid | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN2A, SCN3A, SCN4A |
| carbamazepine | pubchem_cid | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN2A, SCN3A, SCN4A |
| lacosamide | pubchem_cid | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN2A, SCN3A, SCN4A |
| phenytoin | pubchem_cid | 10 | 1 | SCN10A, SCN11A, SCN1A, SCN3A, SCN4A, SCN5A |
| alprazolam | pubchem_cid | 13 | 4 | GABRA4, GABRA6, GABRB1, GABRB2, GABRB3, GABRD |

## The drugs this widens

| drug | ChEMBL targets | GtoPdb targets | examples of what is added |
| --- | --- | --- | --- |
| imipramine | 2 | 3 | CHRM2 |
| doxepin | 2 | 3 | SLC6A4 |
| CID100125017 | 2 | 3 | SLC6A3 |
| fluvoxamine | 1 | 2 | SIGMAR1 |
| CID100016362 | 6 | 7 | CACNA1H, CACNA1I, HRH1, HTR1A |
| moclobemide | 1 | 2 | MAOB |
| CID100071360 | 2 | 3 | HTR1A, HTR6, HTR7 |
| caffeine | 4 | 5 | ITPR1 |
| CID100213046 | 2 | 3 | ADRA2A, ADRA2C, HTR7 |
| CID100006018 | 1 | 2 | SLC18A1 |
| frovatriptan | 2 | 3 | HTR1A |
| CID100005665 | 1 | 2 | SLC32A1 |
| paroxetine | 1 | 2 | P2RX4 |
| CID100004991 | 1 | 2 | BCHE |
| eletriptan | 3 | 5 | HTR1A, HTR1E |

## Where the narrowing needs a read before it is taken

GtoPdb names one gene for these, where ChEMBL's complex named five or more. The narrowing is right in kind - the complex expansion is wrong - but GtoPdb's recall per drug varies, and a drug that acts on several subunits of one receptor may hold a row for only one of them. Listed so the loss is a decision and not a side effect.

| drug | ChEMBL targets | the one GtoPdb names |
| --- | --- | --- |
| amisulpride | 19 | DRD3 |
| lorazepam | 13 | GPR68 |
| ropivacaine | 10 | SCN10A |
| cocaine | 10 | HTR3A |
| lidocaine | 10 | SCN10A |
| carbamazepine | 10 | FZD8 |
| lacosamide | 10 | SCN10A |
| phenytoin | 10 | SCN10A |
| retigabine | 5 | KCNQ1 |

## Rows that modulate an endogenous ligand rather than supplying one

These are what preregistration item 24's modulator relation reads: a partial agonist, a positive or negative allosteric modulation, or a potentiation.

- `CID100002118`: GABRA1, GABRA2, GABRA3, GABRA5
- `CID100002215`: ADRA2A, DRD2, DRD3, DRD4, DRD5, HTR1A
- `CID100002476`: OPRM1
- `CID100002477`: HTR1A
- `CID100002487`: OPRK1, OPRM1
- `CID100002554`: FZD8
- `CID100002802`: GABRA1, GABRA2, GABRA3
- `CID100002803`: ADRA2A
- `CID100002818`: CHRM1
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
