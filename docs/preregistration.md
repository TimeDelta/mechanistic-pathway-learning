# Pre-registration skeleton (to be completed before Phase 3 and posted to OSF)

Frozen in Phase 0 and 1: symptom set (docs/symptom_crosswalk.csv), curated module scaffolding
(docs/curated_pathway_modules.csv), pinned data versions (docs/data_sources.md), evidence grading
rules and weights (configs/evidence_assembly.yaml), splits and metrics (configs/evaluation.yaml).

Primary hypothesis: see experiment_design.md section 6.5 (equifinality independence and
convergence indices; sufficiency test by module ablation).

Primary endpoint: macro-averaged AUPRC under the pathway-wise split, proposed model versus the
best of network proximity, knowledge-graph embedding and relational GNN with sigmoid head;
minimum difference 0.05; 95 percent bootstrap interval excluding zero; five seeds by five folds.

Failure criterion: if the proposed model does not beat the popularity and network proximity
baselines under the gene-wise split, the architecture is abandoned before ablations run.

To fill in: Phase 1 counts per symptom and grade; final symptom set after go/no-go; B5 language
model and prompt; number of flux samples per gene; power statement for the GWAS enrichment test.
