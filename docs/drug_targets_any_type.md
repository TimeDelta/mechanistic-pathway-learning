# Drug targets of any type, and the salt-form mechanism fix

Written 8 October 2026 at the user's request ("fix that bug and incorporate drugs for any type of target"). The new
tables are data/processed/onsides_v3, data/processed/evidence_full_v3 and
data/processed/label_selection/better_v2_full_v3.parquet. Nothing that names evidence_full_v2 or better_v2_full_v2 was
changed: the confirmatory_v2 configurations, configs/lockbox_v2.json and the scorer. Moving to v3 is the user's decision
(last section).

## The bug

ChEMBL records many mechanisms on a salt or prodrug form, with the parent molecule in `parent_molecule_chembl_id`
(1,626 of the 7,561 mechanism rows in the cache). Amitriptyline's mechanisms, for example, are on its hydrochloride
CHEMBL1200964, whose parent is CHEMBL629. The SIDER route (UniChem PubChem CID to ChEMBL id, then the mechanism table)
indexed mechanisms by `molecule_chembl_id` alone. It therefore found no target for a drug whose UniChem id is the
parent, and dropped the drug. The OnSIDES bridge already indexed by both keys, so the two label sources applied
different rules to one drug.

A second gap: UniChem returns no ChEMBL id for 100 SIDER CIDs that the OnSIDES bridge unified with a ChEMBL parent by
RxNorm, UNII or name.

## Changes

- `mechanisms_by_parent_and_molecule` (mechanistic_pathway_learning/perturbation/map_drug_targets_to_graph_nodes.py)
  indexes each mechanism under both keys. The SIDER route (`load_chembl_caches`), the target fetch
  (experiments/fetch_chembl_drug_targets.py), the OnSIDES bridge and the OnSIDES report builder all use it.
  `mechanism_targets` is the one function that turns ChEMBL ids into targets for both label sources.
- The SIDER route (`sider_drug_targets` in assemble_evidence_table.py) uses the bridge's ChEMBL parent only when the
  UniChem route gives no target (`--onsides-bridge`).
- `GraphNodeLookup.nodes_of_target` places every ChEMBL target type on the graph where a node exists:

  | ChEMBL target type | graph nodes |
  |---|---|
  | single protein, complex, family, selectivity group | the gene nodes of its gene symbols (unchanged) |
  | nucleic acid naming a gene (TTR mRNA, SMN2 pre-mRNA) | that gene's node, by the Ensembl id ChEMBL gives as the accession |
  | metal, small molecule, lipid, oligosaccharide | the Human-GEM metabolites of configs/non_protein_drug_targets.csv, every compartment copy, magnitude 1/n |
  | none, for a drug that is itself a graph metabolite | that metabolite, sign +1 (EXOGENOUS SUPPLY; configs/drugs_acting_as_graph_compounds.csv) |
  | another organism, DNA or RNA in general, antibodies, unnamed classes | none |

  configs/non_protein_drug_targets.csv lists the 39 non-protein ChEMBL targets that act on human molecules. 24 map to
  Human-GEM metabolites; the other 15 carry a note saying why not (not in Human-GEM, no single compound named, a gut-lumen
  substrate or another drug).
- New signs: RNAI INHIBITOR, ANTISENSE INHIBITOR, HYDROLYTIC ENZYME, OXIDATIVE ENZYME and REDUCING AGENT −1; POSITIVE
  MODULATOR +0.5; EXOGENOUS SUPPLY +1; BINDING AGENT and CROSS-LINKING AGENT 0.
- A drug seeded on metabolite nodes counts as acting in the metabolic layer (`in_metabolic_layer`).
- The graph-compound table holds only drugs with no mechanism under any ChEMBL form, and the builders refuse an entry
  whose ChEMBL form carries one (`check_graph_compounds_have_no_mechanism`). Lithium is therefore not listed. SIDER has
  it under two CIDs: lithium carbonate (CID 11125), which carries IMPA1 and GSK3 inhibition in ChEMBL and fails the
  single-target rule, and the lithium ion (CID 28486), which carries no mechanism. Listing the ion would have let one
  drug both enter and be excluded.

## Rebuilds

```
OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/build_onsides_reports.py --output-dir data/processed/onsides_v3 --markdown-output <scratch>.md --evidence-dir data/processed/evidence_full_v3
OMP_NUM_THREADS=1 PYTHONPATH=. python -m mechanistic_pathway_learning.evidence.assemble_evidence_table --graph-dir data/processed/graph_full --output-dir data/processed/evidence_full_v3 --extra-reports data/processed/onsides_v3/onsides_reports.parquet --onsides-bridge data/raw/onsides/v3.1.1/ingredient_identifier_bridge.json
OMP_NUM_THREADS=1 PYTHONPATH=. python experiments/build_label_selection.py --evidence-dir data/processed/evidence_full_v3 --mask-grades C --output data/processed/label_selection/better_v2_full_v3.parquet
```

The v2 tables rebuild byte-identically from the code before this change, so every difference below comes from it.
SHA-256 (first 16 hex): onsides_reports 177557ed28a85c12, evidence_records 3b16a40740d0eb79, evidence_reports
1a96d6814e6c98c3, better_v2_full_v3 2aebe4f7e677adbc.

| | v2 | v3 |
|---|---|---|
| observation rows | 5,269 | 5,605 |
| gene rows | 4,090 | 4,090, identical (weights identical; the reliability posterior moves by at most 8.5e-7 because the fit sees more drug reports) |
| drug rows | 1,179 | 1,515 |
| drugs | 142 | 155 |
| qualified SIDER-label pairs | 738 | 1,279 |
| OnSIDES qualifying reports | 1,578 | 1,578 (125 non-qualifying rows of non-nervous-system enzyme drugs now carry metabolite seeds) |
| OnSIDES later-slice candidates | 136 | 136 |
| better_v2 positive pairs / kept pairs / masked negatives | 4,439 / 2,211 / 777 | 4,744 / 2,515 / 777 |
| perturbations loaded on graph_full_neuronal | 1,539 (142 drugs) | 1,551 (154 drugs) |

Most of the change comes from 42 drugs already in v2. In v2 they had OnSIDES reports only, because the SIDER route had
lost their salt-form mechanisms; in v3 they regain their 2015 SIDER label reports. Among them are citalopram,
fluoxetine, paroxetine, fluvoxamine, clomipramine, donepezil, galantamine, zolpidem, morphine, hydromorphone, methadone,
pramipexole, ropinirole, rotigotine, memantine, amantadine, tiagabine, varenicline, pregabalin and buspirone. Their
seeds are unchanged. The qualified SIDER pairs therefore rise from 738 to 1,279. The OnSIDES later-slice candidates
stay at 136, because the raw SIDER 4.1 check had already set the missing pairs aside.

Drugs that enter or leave (group = leakage group under disease_cluster_and_targets on graph_full_neuronal):

| drug | route | ChEMBL target and action | seed nodes | group |
|---|---|---|---|---|
| apomorphine | salt-form index | CHEMBL2331075 D2-like receptor, agonist | DRD2, DRD3, DRD4 | CHEMBL2096905 (a lockbox_v2 group) |
| pergolide | salt-form index | CHEMBL2096905 dopamine receptor, agonist | DRD1 to DRD5 | CHEMBL2096905 (a lockbox_v2 group) |
| pyridostigmine | salt-form index | CHEMBL220 acetylcholinesterase, inhibitor | ACHE | CHEMBL2095233 (a lockbox_v2 group) |
| phenelzine | salt-form index | CHEMBL2095205 monoamine oxidase, inhibitor | MAOA, MAOB | cluster:MAOA (a lockbox_v2 group) |
| maprotiline | salt-form index | CHEMBL222 noradrenaline transporter, inhibitor | SLC6A2 | cluster:SLC6A3 |
| tramadol | salt-form index | CHEMBL233 mu opioid receptor, agonist | OPRM1 | CHEMBL2095181 |
| propoxyphene | salt-form index | CHEMBL233 mu opioid receptor, agonist | OPRM1 | CHEMBL2095181 |
| procaine | salt-form index | CHEMBL2331043 sodium channel, blocker | ten SCN genes | cluster:ABCA3 |
| methohexital | salt-form index | CHEMBL2093872 GABA-A receptor, positive allosteric modulator | 13 GABA-A subunits | cluster:ABCA3 |
| pentazocine | salt-form index | CHEMBL287 sigma-1 receptor, modulator (sign 0) | SIGMAR1 | cluster:ADH1C |
| benztropine | bridge parent CHEMBL1201203 | CHEMBL216 muscarinic M1, antagonist | CHRM1 | CHEMBL216 |
| scopolamine | bridge parent CHEMBL569713 | CHEMBL216 muscarinic M1, antagonist | CHRM1 | CHEMBL216 |
| etomidate | bridge parent CHEMBL681 | CHEMBL2093872 GABA-A receptor, positive modulator | 13 GABA-A subunits | cluster:ABCA3 |
| gamma-aminobutyric acid | graph compound | none in ChEMBL; MAM00970 (GABA), +1 | five compartment copies | not loaded: its one pair is an indication, which the selection does not hold |
| doxepin (leaves) | salt-form index | gains CHEMBL222 (noradrenaline transporter) from its hydrochloride beside CHEMBL231 (H1) | | was a lockbox_v2 drug; now fails the single-target rule |

No drug enters through a nucleic-acid or a metabolite target. The drugs with such targets are absent from the qualifying
data for other reasons:
- the TTR mRNA drugs (patisiran, inotersen, vutrisiran) and the SMN2 drugs (nusinersen, risdiplam) have no OnSIDES
  label statement on a crosswalk symptom and postdate SIDER;
- the chelators, cysteamine and the enzyme replacements have no nervous-system ATC code.

The code places them when they qualify, as the tests show.

## What stays out

Of the 206 SIDER CIDs with a nervous-system ATC code and a crosswalk label event, 127 qualify. 62 have more than one
mechanism target. 17 have no human mechanism target in ChEMBL: amobarbital, chloroquine, levomepromazine, lisuride, lithium (the ion
CID; the carbonate CID has two targets), lormetazepam, mianserin, nefopam, nitrazepam, nitrous oxide, paraldehyde,
phenobarbital, pizotifen, propericiazine, salicylate, stiripentol and zuclopenthixol. The mechanism table holds no row
for their ChEMBL ids, under either key. Amobarbital has no ChEMBL id from either route, and chloroquine's only row is
for a Plasmodium falciparum target. None of the 17 is a listed graph compound (GABA and tryptophan are). Nitrous oxide
and salicylate are not in Human-GEM. Every target of a qualifying drug has a graph node.

## Caveats

- Pentazocine enters with an all-zero input. ChEMBL records only a sigma-1 "MODULATOR" mechanism for it (on salts of
  CHEMBL100116, ChEMBL's PENTAZOCINE) and no opioid-receptor mechanism. It joins gabapentin, pregabalin and brivaracetam,
  which the preregistration lists as zero-sign drugs; the trainer warns about them.
- Metabolite seeds are placed on graph_full, as the v2 gene seeds were (preregistration, "evidence_full_v2 placed drug
  seeds on graph_full"). graph_full_neuronal adds a vesicle copy of GABA (MAM00970v) that a GABA seed would miss.
- Scopolamine and benztropine are seeded on CHRM1 alone because ChEMBL lists only M1. Both act on other muscarinic
  subtypes too.

## Effect on the confirmatory experiment (membership only; no lockbox outcome or output was opened)

configs/lockbox_v2.json records the SHA-256 of evidence_full_v2's records and of better_v2_full_v2. Under v3:
- doxepin, a lockbox_v2 drug, leaves the data;
- apomorphine, pergolide, pyridostigmine and phenelzine join lockbox_v2 groups;
- 10 of the 28 lockbox_v2 drugs are among the 42 that regain SIDER label reports;
- no v3 leakage group holds lockbox_v2 and development perturbations together, so the lockbox groups still separate
  cleanly.

Adopting v3 would need a lockbox file derived from the same groups on v3, and the development runs would need to be
repeated on v3. Both are the user's decision; nothing here changes configs/lockbox_v2.json, the confirmatory_v2
configurations or the scorer.
