# mechanistic-pathway-learning

Constrained graph learning of symptom-to-mechanism pathways over a compartmentalized
physiology graph. The model maps a localized molecular perturbation (a monogenic enzyme
defect, a single-target drug) to a distribution over psychiatric symptoms, explained by
one or more sparse pathway modules combined through a noisy-OR, with no diagnosis labels
anywhere in the inputs or targets.

The full experiment design, including the assumptions it argues against, is in
docs/experiment_design.md. Read that first.

## Status

Scaffold only. Implemented and tested: the noisy-OR pathway module head, soft constraint
losses, a pure-PyTorch relational message passing encoder, grouped and pathway-wise and
time splits, evidence grading, currency metabolite tagging, equifinality indices, ranking
and calibration metrics and two negative controls. Everything that touches real data is a
stub with a written contract (see the docstring at the top of each file).

Next step: Phase 1 of the design (graph build, evidence assembly, go/no-go counts).

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
