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
- The node descriptors and the cell-class weights were not finished in this session (see the last section).

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

## Not finished: node descriptors and cell-class weights

The ESM-2 embeddings of the 20,431 reviewed proteins take about three hours on two cores and were still running when
the rest of this session's work was done. The steps after them, in order, with the checks to hold them to:

1. `experiments/build_protein_descriptors.py`, then the same with `--per-entry-only` (1,880 targets, 19,628 proteins in
   the fit, rank 64, coverage 12,460 of 12,627; the sign of each of the 64 directions is not fixed by the fit, so a new
   fit matches the earlier table up to sign and not by hash);
2. `experiments/build_complex_descriptors.py` (1,681 entities, 1,527 with a vector);
3. `experiments/build_node_descriptors.py --graph-dir data/processed/graph_full_neuronal --complex-descriptors ...`
   (84 columns);
4. `experiments/build_dopaminergic_expression.py` on the four brain atlas dissections (872 nuclei of cluster 395);
5. `experiments/build_brain_expression_descriptors.py` and `experiments/build_cell_class_weights.py` on the merged
   neuronal graph, then `experiments/write_expression_tables.py`;
6. `experiments/build_gene_protein_split.py` again with the `--descriptor-table` and `--cell-class-weights` pairs (138
   columns), then the binder variant and the drug-entry variant with theirs.

Training needs these tables. Nothing built or measured in this session does.
