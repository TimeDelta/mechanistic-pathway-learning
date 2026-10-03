# mechanistic-pathway-learning

Constrained graph learning of symptom-to-mechanism pathways over a compartmentalized
physiology graph. The model maps a localized molecular perturbation (a monogenic enzyme
defect, a single-target drug) to a distribution over psychiatric symptoms, explained by
one or more sparse pathway modules combined through a noisy-OR, with no diagnosis labels
anywhere in the inputs or targets.

The full experiment design, including the assumptions it argues against, is in
docs/experiment_design.md. Read that first.

## Status

Phase 1 is measured and Phases 2 and 3 have started on the monogenic slice. Implemented and tested:
the noisy-OR pathway module head (sum or mean pooling over hard-concrete supports), soft constraint
losses, a relational message passing encoder with absolute and perturbation-difference field
readings, the sigmoid-head baseline B3, grouped splits by gene or by disease cluster, curated-module
and reconstruction-subsystem pathway hold-outs, the time split, the HPO loader with frequency
qualifiers, provenance and crosswalk term exclusions, the evidence grading policies, the learned
evidence reliability model, the document-grounded appraisal rubric, currency tagging, equifinality
indices, ranking and calibration metrics, the label-permutation and rewiring controls, the SIDER
loader, the SIDER-to-ChEMBL mapping, the physiology graph build (Human-GEM with subsystems plus the
OmniPath, small-molecule and CollecTRI layers when their files are present), the evidence table with
disease clusters, popularity and random-walk baselines, a resumable training harness with
validation early stopping, a batch driver, an aggregator with paired bootstrap comparisons, a
module-analysis script with the sufficiency test, and mechanism cards.

Measured so far: docs/phase1_counts.md (counts), docs/hpo_term_audit.md (what each HPO term
contributes to each symptom), docs/phase2_baselines.md and docs/phase2_baselines_disease_cluster.md
(baselines under three splits with permutation controls), docs/phase3_main_model.md (B6 and B3).
Two results shape the design: the gene-wise split leaks through genes that share a disease, so the
disease-cluster grouping is the default; and the sixteen curated modules hold too few perturbations
for a pathway-wise endpoint, so Human-GEM subsystems are proposed instead (design section 6.1). A
monogenic time split dates each gene-symptom pair by the publication of the PubMed reference behind
its OMIM annotation (docs/hpo_reference_publication_dates.json). docs/references.bib holds the design
document's reference list for import into a reference manager.

Still stubs: OnSIDES and literature loaders, the diagnosis proxy audit, the flux-sampling job
(written, not yet run), the flux-feature and language-model baselines, the MAGMA wrapper. The pharmacological evidence class and the signaling layers need the SIDER, ChEMBL and
OmniPath hosts, which were not reachable from the environment that produced the current tables.

Reproduce the data layer (downloads about 110 MB; raw files stay out of git):

```
mkdir -p data/raw/hpo data/raw/Human-GEM/model data/raw/sider_4.1 data/raw/omnipath
curl -L -o data/raw/hpo/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo
curl -L -o data/raw/hpo/genes_to_phenotype.txt https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt
for f in Human-GEM.xml Human-GEM.yml genes.tsv metabolites.tsv reactions.tsv; do curl -L -o data/raw/Human-GEM/model/$f https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/$f; done
for f in meddra_all_se.tsv.gz meddra_freq.tsv.gz meddra_all_indications.tsv.gz drug_names.tsv drug_atc.tsv; do curl -L -o data/raw/sider_4.1/$f http://sideeffects.embl.de/media/download/$f; done
curl -L -o data/raw/omnipath/omnipath_interactions.tsv "https://omnipathdb.org/interactions?datasets=omnipath,pathwayextra,ligrecextra&genesymbols=1&fields=sources,references,type,curation_effort,n_references&organisms=9606&format=tsv"
curl -L -o data/raw/omnipath/collectri_interactions.tsv "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=1&fields=sources,references,n_references&organisms=9606&format=tsv"
curl -L -o data/raw/omnipath/small_molecule_protein.tsv "https://omnipathdb.org/interactions?types=small_molecule_protein&genesymbols=1&fields=sources,references,type&organisms=9606&format=tsv"
python experiments/fetch_chembl_drug_targets.py          # UniChem and ChEMBL, cached, resumable
python -m mechanistic_pathway_learning.graph.build_physiology_graph
python experiments/audit_hpo_term_expansion.py            # docs/hpo_term_audit.md
python -m mechanistic_pathway_learning.evidence.assemble_evidence_table
python experiments/run_phase1_counts.py
python experiments/run_baselines.py --group-by disease_cluster --with-kg-embedding --output-dir runs/baselines_disease_cluster --markdown-output docs/phase2_baselines_disease_cluster.md
python experiments/run_main_model_batch.py --configuration b6_default --group-by disease_cluster     # five folds, resumable
python experiments/run_main_model_batch.py --configuration b3_sigmoid --group-by disease_cluster
python experiments/aggregate_main_model_runs.py --run-dirs runs/b6_default_disease_cluster runs/b3_sigmoid_disease_cluster --baseline-results runs/baselines_disease_cluster/results.json
python experiments/analyze_pathway_modules.py --run-dir runs/b6_default_disease_cluster
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
docs/            experiment design, symptom crosswalk (with HPO exclusions), curated pathway modules, data sources, generated tables
configs/         graph build, evidence assembly, model, evaluation, baselines, ablations
slurm/           array job for flux sampling; resumable training job for a runtime-limited cluster
mechanistic_pathway_learning/
  graph/         physiology graph build (Human-GEM + OmniPath + CollecTRI + ChEMBL), currency tagging
  evidence/      HPO, SIDER/OnSIDES, PubTator3/CTD/SemMedDB loaders; evidence grading; proxy audit
  perturbation/  flux-sampling features per gene knockout; drug target mapping
  models/        encoder, noisy-OR head, soft constraint losses, MDL selection, baselines
  evaluation/    splits, metrics, negative controls, GWAS enrichment, equifinality tests, mechanism cards
experiments/     phase runners, HPO audit, batch driver, aggregator, module analysis
tests/
```

## Conventions

Descriptive names over comments; no hard symptom-level constraints anywhere in the model;
unobserved perturbation-symptom pairs are unlabelled, not negative; every long job is
idempotent and resumable; MedDRA term tables are never committed.
