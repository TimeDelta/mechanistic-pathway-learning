# Rebuilding the data layer in a new container: what reproduced and what did not

Written 10 October 2026, by the session that followed the one of 6 to 10 October. `data/raw`, `data/processed` and
`runs/` are not in git, so a new container starts with none of them, and the confirmatory inputs had to be rebuilt
before anything else could be done. This file records each step, the check it was held to and the result. The commands
are in `scripts/fetch_raw_sources.sh` and `scripts/rebuild_data_layer.sh`.

The container had 2 cores and 7 GB of memory, against the 4 cores of the earlier one.

## The short version

- Every graph table reproduces its recorded node and edge counts, and every evidence table reproduces every recorded
  count.
- Three tables with a recorded SHA-256 match it byte for byte: the neuronal graph's `nodes.parquet`, the OnSIDES report
  table of v3 and the evidence report table of v3.
- Four recorded hashes do not reproduce: the neuronal graph's `edges.parquet`, for a cause that is found and fixed, and
  the evidence record tables and label selections of v3 and of v3 with parkinsonism, for a cause that is not found.
- Three sources had moved since they were first fetched. Two are now pinned to an archived copy that carries the
  recorded hash or content. The third is a live service whose change touches four identifiers and no table.
- Two defects in the builds made a table depend on Python's hash seed. One is fixed here. The other is worked around
  by copying the committed release, because fixing it would change the pinned row order.
- The node descriptors and the cell-class weights are rebuilt, with every recorded count, and the descriptor table is
  not the registered file: its SHA-256 is `0ba1f817029f32e0` against the registered `df38d78c78ffe6d4` (last section).
- On this container the confirmatory message-passing configuration trains at 9.4 s a step (11.4 minutes an epoch) and
  the linear-response one at 4.9 s (6.0 minutes); message passing with 6 or 9 layers is killed for memory at batch 16.

## Library versions

The tables of `data/releases/v0.4` were written by pandas 2.3.3 and pyarrow 25.0.1 (their parquet metadata says so).
A new container installs pandas 3 and a later pyarrow, and a parquet file written by another pyarrow has another hash
whatever its content. With pandas 2.3.3 and pyarrow 25.0.1 installed, the neuronal graph's `nodes.parquet` came out at
the recorded hash, `48ebf531`. Install those two before a rebuild that is to be compared by hash.

## Sources

| source | recorded | fetched on 10 October | result |
| --- | --- | --- | --- |
| Human-GEM | 2.0.1, model date 2026-03-26, no hash | the `main` branch held 2.1.1, dated 2026-10-07 | **moved**; pinned to the tag `v2.0.1` |
| HGNC complete set | 16,963,116 bytes, sha256 `2b4224ea` | the rolling URL served 16,973,159 bytes | **moved**; the monthly archive copy of 2026-10-02 has the recorded hash |
| RxNav (OnSIDES bridge) | bridge of 5 October in `data/releases/v0.4` | 4 of 1,866 ingredients had lost their UNII codes | **moved**; no other field of the bridge differs |
| HPO | release 2026-09-01, no hash | from the release tag `v2026-09-01` | no hash to check; pinned to the tag in place of `latest` |
| OmniPath | 126,124, 64,515 and 6,553 rows | 126,124, 64,515 and 6,553 rows | equal |
| GTEx v10 median TPM | 8,846,936 bytes | 8,846,936 bytes | equal |
| Human Protein Atlas, two files | sha256 `108fb699`, `a287c1ab` | the same | equal |
| OnSIDES v3.1.1 archive | 84,862,297 bytes; 6,928,666 mentions, 41,119 labels, 1,866 ingredients | the same | equal |
| ESM-2 weights and contact regression | sha256 `7f21e80e`, `16641e05` | the same | equal |
| UniProt reviewed proteome | release 2026_03, sha256 `35f4b89f` | the same | equal |
| GO annotations and ontology | sha256 `db472faf`, `b08d45b2` | the same | equal |
| HPO `hp-base.owl` | sha256 `42c7e775` | the same | equal |
| ChEBI | sha256 `72a9574f` | the same | equal |
| Reactome | release 97, 15 pathways, 12,086,097 bytes | release 97, 12,086,097 bytes | equal |
| UniProt plasma binders | sha256 `eab65502` | the same | equal |
| ChEMBL | ChEMBL_37, 7,561 mechanisms, 1,177 of 1,430 SIDER drugs mapped | the same | equal |
| Fraction-unbound database | sha256 `9fe90e80` | the same | equal |
| FDA label sentences (openFDA) | `docs/label_plasma_binders.md` | 154 drugs fetched | the document regenerates with no difference |

Human-GEM is the one that would have gone unnoticed. The README fetched the `main` branch, and a build from `main` on
10 October gave 12,748 reactions and 8,374 metabolites against the 12,877 and 8,460 of the release, with no error,
because every layer file is optional in the build. `scripts/fetch_raw_sources.sh` fetches the tag.

Seven inputs had a host in `docs/data_sources.md` and no URL anywhere in the repository: the two HPO annotation files
the evidence assembler reads when present (`phenotype.hpoa`, `genes_to_disease.txt`), the GTEx object, the Orphadata
product, the OnSIDES asset, the ESM-2 weights and the adult human brain atlas dissections. The fetch script now holds
them. Two inputs had no script at all and are written down for the first time:

- `data/raw/uniprot/uniprot_human_manganese_cofactor.tsv.gz`, which the manganese variant reads. A UniProt query for
  reviewed human entries whose cofactor annotation names ChEBI 29035 returns 281 entries on release 2026_03, and the
  variant reads them as 93 proteins requiring manganese alone and 188 accepting it among other metals, the two counts
  the earlier build recorded.
- `data/raw/onsides/v3.1.1/parquet/`, the OnSIDES tables as parquet. `scripts/convert_onsides_to_parquet.py` makes
  them. The conversion has to keep the literal string "NA" in `label_section`, which is OnSIDES' value for a non-US
  label: a reader that takes "NA" for a missing value drops every EU, UK and JP statement. That the conversion is the
  one the earlier build used is shown by the OnSIDES report table below, which matches its recorded hash.

## Tables

| table | recorded | rebuilt | result |
| --- | --- | --- | --- |
| `graph_full` | release v0.4 | same rows; 16,056 `catalyzed_by` rows in another order | content equal; the release copy is used |
| `graph_full_manganese` | 93 and 188 proteins | 93 and 188; 9 nodes and 99 edges added | equal |
| `graph_full_neuronal` nodes | 36,865; sha256 `48ebf531` | 36,865; `48ebf531` | **hash equal** |
| `graph_full_neuronal` edges | 266,966; sha256 `78cf6c07` | 266,966; `b078e173` | count equal, hash not (below) |
| `graph_full_neuronal_split` | 49,574 nodes, 279,757 edges, 12,709 protein nodes, 22 shared, CALCA the one gene with several | the same | equal; no hash was recorded |
| `graph_full_neuronal_split_binders` | 279,784 edges; 27 `binds` edges, sign 0, 11 binders, 19 cargo | the same | equal; the generated part of `docs/plasma_binder_graph.md` regenerates |
| OnSIDES reports, v3 | sha256 `177557ed28a85c12` | `177557ed28a85c12` | **hash equal** |
| evidence reports, v3 | sha256 `35992a4155bb4429` | `35992a4155bb4429` | **hash equal** |
| evidence records, v3 | 5,604 rows (4,090 gene, 1,514 drug), 154 drugs, 22 symptoms, 1,066 clusters; `d4e81f7ca1f1b087` | the same counts; `ae73de729b4c6403` | counts equal, hash not (below) |
| label selection `better_v2` on v3 | 5,521 rows, 4,744 positive, 2,515 kept, 777 masked negatives; `2aebe4f7e677adbc` | the same counts; `a90d427696abfc13` | counts equal, hash not |
| evidence records, v3 with parkinsonism | `17955be91f06be8c` | `ce7391cc4b76c654` | `docs/parkinsonism_symptom.md` regenerates with these two hash rows as its only difference |
| label selection on v3 with parkinsonism | `8d61694e300f02a3` | `cc15a045d75ddb3d` | as above |

New hashes of this container, for the tables that had none: `graph_full_neuronal_split` nodes `cd6ad3e7f0f3974d`, edges
`cb70bc7890112946`, `gene_to_protein.parquet` `9cf02ccc84953b60`; `graph_full_neuronal_split_binders` nodes
`623321fec4c0b6b2`, edges `64d479b56934a6a6`. The split build gave the same three hashes under three hash seeds.

## Two builds depended on Python's hash seed

**The neuronal variant, fixed.** Two builds of the same inputs under `PYTHONHASHSEED=0` and `PYTHONHASHSEED=1` gave
two `edges.parquet` files with the same rows and 381 of them in another position: 311 `changes_redox_pool`, 53
`changes_ion_pool` and 17 `redox_capacity`. Both loops iterated over a set of strings
(`graph/redox_pools.py`, `set(consumed) | set(produced)`; `graph/reactome_import.py`, the ion set of
`ions_moved_inward`). They are sorted now, and three builds under three seeds give `b078e173`. The recorded `78cf6c07`
was one seed's order and cannot be reproduced except by chance. The rows are the same set either way.

Row order matters downstream in two places: the degree-preserving rewiring of the H2 control swaps edges by position,
so another order is another rewired graph for the same seed, and a sum over a node's incoming edges is taken in row
order, which moves the last digits of a float.

**The base graph, not fixed.** `build_physiology_graph.py` writes a reaction's `catalyzed_by` edges in the order
cobra's `reaction.genes` yields them, and that is a frozenset. A rebuild has the release's rows with 16,056 of them
moved. Sorting would make the build reproducible and would also give a third order that matches neither the release nor
any earlier run, so the rebuild copies `data/releases/v0.4/graph_full` and the builder is left alone. If the base graph
is ever rebuilt on purpose, the sort belongs in that change.

## The evidence record tables: the cause is not found

The report tables match their recorded hashes, so the inputs of the assembly are the earlier build's inputs byte for
byte. The record tables and the label selections reproduce every recorded count and not the recorded hash. What was
tried:

- The assembly is deterministic here: the same hash under two hash seeds and with one or two threads.
- The fitted column is not the whole cause. `reliability_posterior` comes from a scikit-learn logistic regression and
  could differ in its last digits between two installs. But the label selection does not read it: a record table with
  every posterior value perturbed gave the same selection file, `a90d427696abfc13`. So the earlier selection, whose
  hash differs, was built from records that differ from these in something the selection does read (a frequency, a
  grade or the row order).
- `docs/parkinsonism_symptom.md` compares the two tables column by column, and it regenerates with the two hash rows
  as its only difference, so no reported number moves.

Settling it needs the earlier files, which are on the disk of the earlier session and nowhere else. Until then the
statement that holds is narrower than "reproduced": the same reports, the same counts, another file.

This matters for one thing. `configs/lockbox_v2.json` pins the evidence table and the label selection by file hash, and
`read_lockbox` refuses a table with another hash. A lockbox drawn on v3 will pin files, so those files have to travel
with the repository as a data release. Rebuilding them in the next container will not give the pinned bytes.

## What the earlier session's disk still holds

Everything under `runs/` (the pilots, the module analyses, the job scripts the runbook describes in prose) and the
earlier `data/processed`, with the stored protein descriptors whose column signs every later descriptor step copies.
If that session can still be opened, a data release cut from its disk would settle the open hash question above and
save the descriptor rebuild below.

## Node descriptors and cell-class weights

Rebuilt on 10 October after the embeddings finished. The commands are the `descriptors` stage of
`scripts/rebuild_data_layer.sh`.

| step | recorded | rebuilt | result |
| --- | --- | --- | --- |
| ESM-2 embeddings | 20,431 reviewed proteins | 103 batches, 20,431 proteins; 2 hours 55 minutes on two cores | count equal |
| protein descriptors | 1,880 targets, 19,628 proteins in the fit, rank 64, coverage 12,460 of 12,627 | the same | the report regenerates with no difference from `docs/protein_descriptor_report.md`, every held-out R squared included |
| complex descriptors | 1,681 entities, 1,527 with a vector; median correlation with a single member's protein descriptor 0.818 | 1,527 with a vector; 0.816 | counts equal; the one fitted figure the report prints differs in the third decimal |
| node descriptors, merged neuronal graph | 84 columns | 84 columns | equal |
| brain atlas dissections | four files, SHA-256 prefixes in `docs/data_sources.md` | the four prefixes match | equal |
| dopaminergic class | 872 nuclei of cluster 395 | 872 (581, 192, 62 and 37 in the four dissections) | equal |
| brain expression descriptors | 9,962 gene nodes without an Ensembl id, 9,783 resolved by symbol, 12,626 of 12,810 gene nodes and 7,752 of 7,759 reactions with expression | the same | equal |
| expression tables | the script refuses to write unless the recomputed reaction rows equal the stored ones | largest difference 0.0 in both tables | equal |
| split descriptor table | 138 columns; SHA-256 `df38d78c78ffe6d4` | 49,574 rows, 138 columns; `0ba1f817029f32e0` | **shape equal, hash not** |
| split cell-class weights | none recorded | 49,574 rows, 12 classes; `84976a5e15a9867a` | no hash to compare |
| drug entry variant | new | 49,728 rows in both tables, a zero row per drug node in the descriptors and a row of ones in the weights | new |

**The descriptor table is not the registered file.** The preregistration names
`full_neuronal_split_descriptors_brain_expression.parquet` by its SHA-256. The rebuilt table has the registered shape
and every recorded count, and other bytes. The 64 protein columns are a fit to embeddings computed on another machine,
so their last digits differ, as the complex report's 0.816 against 0.818 shows, and the sign of each direction is not
fixed by the fit. Whether the rebuilt directions agree in sign with the registered ones cannot be checked without the
registered table, which is on the earlier session's disk. A model does not care about a column's sign, since the
first layer can absorb it, so a run on the rebuilt table is a run of the registered configuration on equivalent
inputs and not on the registered file. A confirmatory run needs one of two things: the registered table itself, or
an amendment that pins the table it will read.

**What a new container costs now.** About three hours for the embeddings, 1.24 GB for the brain atlas, 45 minutes for
the OnSIDES identifier bridge and a few minutes for each other step. A data release that carries the processed tables
would remove all of it and would settle which file a pin refers to.

## What training costs on this container

`run_main_model.py --timing-batches 6` on the drug entry variant with the rebuilt tables, the noisy-OR head and the
confirmatory fitting arguments, on 2 cores with about 6 GB usable. It runs six training steps and exits: no split
directory, no validation, no score.

| configuration | seconds a step | minutes an epoch (69 steps) |
| --- | --- | --- |
| message passing, 3 layers | 9.4 | 11.4 |
| message passing, 3 layers, mechanism also before the first layer | 9.7 | 11.7 |
| message passing, 6 layers | killed for memory at 6.1 GB | |
| message passing, 9 layers | killed for memory at 6.1 GB | |
| linear response, 8 steps | 4.9 | 6.0 |

No lockbox was passed, so the step count is for 1,089 training perturbations of the 1,568; with the held-out
perturbations left out an epoch is shorter. A deeper message-passing variant needs gradients through its last rounds
only, or a smaller batch, to run here at all.
