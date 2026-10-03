# mechanistic-pathway-learning

Constrained graph learning of symptom-to-mechanism pathways over a compartmentalized
physiology graph. The model maps a localized molecular perturbation (a monogenic enzyme
defect, a single-target drug) to a distribution over psychiatric symptoms, explained by
one or more sparse pathway modules combined through a noisy-OR, with no diagnosis labels
anywhere in the inputs or targets.

The full experiment design, including the assumptions it argues against, is in
docs/experiment_design.md. Read that first.

## Status

Phase 1 is measured and Phase 2 has started. Implemented and tested: the noisy-OR pathway
module head, soft constraint losses, a sparse relational message passing encoder, grouped
and pathway-wise and time splits, evidence grading, the learned evidence reliability model,
the document-grounded LLM appraisal rubric, currency metabolite tagging, equifinality
indices, ranking and calibration metrics, two negative controls, the HPO and SIDER loaders,
the SIDER-to-ChEMBL drug-target mapping, the physiology graph build (Human-GEM plus
OmniPath signaling, small-molecule ligand and CollecTRI transcription layers), the
evidence table assembly, the popularity and random-walk baselines and a resumable
training harness for the proposed model. Measured counts are in docs/phase1_counts.md
and the first baseline numbers in docs/phase2_baselines.md. Still stubs: OnSIDES and
literature loaders, the diagnosis proxy audit, the flux-sampling job (written, not yet
run), the embedding and GNN-sigmoid and flux-feature and language-model baselines, the
MAGMA wrapper and the mechanism cards.

Reproduce the data layer (downloads about 100 MB; raw files stay out of git):

```
mkdir -p data/raw/hpo data/raw/Human-GEM/model data/raw/sider_4.1 data/raw/omnipath
curl -L -o data/raw/hpo/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo
curl -L -o data/raw/hpo/genes_to_phenotype.txt https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt
for f in Human-GEM.xml genes.tsv metabolites.tsv reactions.tsv; do curl -L -o data/raw/Human-GEM/model/$f https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/$f; done
for f in meddra_all_se.tsv.gz meddra_freq.tsv.gz meddra_all_indications.tsv.gz drug_names.tsv drug_atc.tsv; do curl -L -o data/raw/sider_4.1/$f http://sideeffects.embl.de/media/download/$f; done
curl -L -o data/raw/omnipath/omnipath_interactions.tsv "https://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,ligrecextra&genesymbols=1&fields=sources,references,type,curation_effort,n_references&organisms=9606&format=tsv"
curl -L -o data/raw/omnipath/collectri_interactions.tsv "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=sources,references,n_references&organisms=9606&format=tsv"
curl -L -o data/raw/omnipath/small_molecule_protein.tsv "https://omnipathdb.org/interactions?types=small_molecule_protein&genesymbols=1&fields=sources,references,type&organisms=9606&format=tsv"
python experiments/fetch_chembl_drug_targets.py          # UniChem and ChEMBL, cached, resumable
python -m mechanistic_pathway_learning.graph.build_physiology_graph
python -m mechanistic_pathway_learning.evidence.assemble_evidence_table
python experiments/run_phase1_counts.py
python experiments/run_baselines.py
python experiments/run_main_model.py --fold 0           # one fold; see slurm/train_resumable.sbatch
```

## Setup

```
python -m venv ~/.venvs/mpl && source ~/.venvs/mpl/bin/activate
pip install -e ".[dev]"            # model and evaluation code
pip install -e ".[graph]"          # add cobra, omnipath, decoupler, chembl client for Phase 1
pytest
```

## Layout

```
docs/            experiment design, symptom crosswalk, curated pathway modules, data sources
configs/         graph build, evidence assembly, model, evaluation, baselines, ablations
slurm/           array job for flux sampling; resumable training job for a runtime-limited cluster
mechanistic_pathway_learning/
  graph/         physiology graph build (Human-GEM + OmniPath + CollecTRI + ChEMBL), currency tagging
  evidence/      HPO, SIDER/OnSIDES, PubTator3/CTD/SemMedDB loaders; evidence grading; proxy audit
  perturbation/  flux-sampling features per gene knockout; drug target mapping
  models/        encoder, noisy-OR head, soft constraint losses, MDL selection, baselines
  evaluation/    splits, metrics, negative controls, GWAS enrichment, equifinality tests, mechanism cards
experiments/     phase runners
tests/
```

## Conventions

Descriptive names over comments; no hard symptom-level constraints anywhere in the model;
unobserved perturbation-symptom pairs are unlabelled, not negative; every long job is
idempotent and resumable; MedDRA term tables are never committed.
