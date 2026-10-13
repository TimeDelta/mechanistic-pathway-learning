# What the v4 rebuild costs in leakage

Measured by `experiments/label_review/rebuild_leakage_cost.py` before the rebuild is run, with the definitions of `docs/label_source_review.md` 26.1 so the numbers compare with the table there. Nothing is rebuilt, no held-out set is drawn.

| scheme | drugs | kept positive pairs | drug-holding groups | two largest blocks, by drugs | largest block's share of drugs / of kept pairs |
| --- | --- | --- | --- | --- | --- |
| three targets, the table in use | 222 | 2975 | 17 | 114, 71 | 51% / 27% |
| no limit, every drug | 230 | 3011 | 15 | 194, 6 | 84% / 48% |
| no limit, orphenadrine set aside | 229 | 3010 | 16 | 116, 77 | 51% / 27% |
| the rebuild: no limit, orphenadrine out, the two beta blockers in | 231 | 3021 | 17 | 116, 77 | 50% / 27% |

## Bridging drugs after the rebuild (the rebuild: no limit, orphenadrine out, the two beta blockers in)

Each row is the number of drugs the largest block loses when that drug alone is removed, which is the stated rule for finding a bridge.

| drug | drugs removed from the largest block |
| --- | --- |
| tapentadol | 19 |
| agomelatine | 4 |
| betahistine | 2 |
| fenfluramine | 2 |
| bupropion | 1 |
| melatonin | 1 |
| methamphetamine | 1 |
| amisulpride | 1 |

## The rebuild costs nothing in leakage

Against the table in use the rebuild gives **9 more drugs and 46 more kept positive pairs**, and the split is slightly
better rather than worse: the largest block holds 50.2 percent of the drugs against 51.4 percent now, the two-block
structure survives at 116 and 77 against 114 and 71, and the number of drug-holding groups is the same at 17. The two
beta blockers the central-action rule admits add nothing to the largest block, which is what the graph predicted: no
drug already in the table acts on `ADRB1` or `ADRB2`, so they form a group of their own.

Lifting the limit without setting orphenadrine aside is the arm that would be ruinous, and it is the reason the
exclusion is a rule: 194 drugs in one block, 84 percent of the drugs and 48 percent of the kept pairs.

## But the bridging rule names a second drug on the rebuilt table

Applying the registered procedure again — remove each drug in turn and measure the largest block — names **tapentadol**,
and it is a bridge in the same sense orphenadrine was. It seeds `GENE:OPRM1` and `GENE:SLC6A2`, so it joins the
18-drug opioid block to the 97-drug monoamine block:

| drugs set aside | blocks by drugs | largest block's share of drugs |
| --- | --- | --- |
| orphenadrine | **116**, 77, 6, 6 | 50.2 percent |
| orphenadrine, tapentadol | **97**, 77, 18, 6 | 42.4 percent |
| orphenadrine, tapentadol, agomelatine | 93, 77, 18, 6 | 40.6 percent |

Orphenadrine joined blocks of 116 and 77; tapentadol joins blocks of 97 and 18. Both satisfy the rule as registered:
"a drug whose target set joins two otherwise separate leakage blocks is set aside, found by removing each drug in turn."

**The registered rule has no stopping condition, and that is an omission in my drafting of the amendment.** Applied once
it names orphenadrine. Applied again to the result it names tapentadol. Applied a third time it names agomelatine, which
moves 4 drugs, and the sequence continues down to drugs that move one. Nothing in the registration says whether the rule
is applied once, until no drug removes more than some number, or until the largest block falls below some share of the
drugs. Three ways to close it, each defensible and each giving a different drug set:

1. **Once**, on the table the rule was stated for: orphenadrine alone, 231 drugs, largest block 50.2 percent. The
   registration reads most naturally this way and it is the smallest change.
2. **Until no remaining drug joins two blocks of at least ten drugs each**: orphenadrine and tapentadol, 230 drugs,
   largest block 42.4 percent. The threshold is what makes it terminate, and ten is a number that has to be stated
   rather than found.
3. **Until the largest block is below half the drugs**: orphenadrine and tapentadol, the same 230 drugs here, but the
   condition is on the outcome of the split rather than on the drug, so it would keep removing drugs on another table.

Whichever is chosen has to be written down before the held-out set is drawn, because every one of them changes which
drugs are in the pool the lockbox is drawn from.
