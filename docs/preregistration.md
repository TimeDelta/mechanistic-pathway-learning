# Pre-registration skeleton (to be completed before Phase 3 and posted to OSF)

Frozen in Phase 0 and 1: symptom set and HPO term exclusions (docs/symptom_crosswalk.csv, audited in
docs/hpo_term_audit.md), curated module scaffolding (docs/curated_pathway_modules.csv), pinned data
versions (docs/data_sources.md), evidence grading policy and weights (configs/evidence_assembly.yaml:
grade A by curated-synopsis provenance with frequency scaling), splits and metrics (configs/evaluation.yaml).

Primary hypothesis: see experiment_design.md section 6.5 (equifinality independence and
convergence indices; sufficiency test by module ablation).

Primary endpoint: macro-averaged AUPRC under the pathway-wise split, proposed model versus the
best of network proximity, knowledge-graph embedding and relational GNN with sigmoid head;
minimum difference 0.05; 95 percent bootstrap interval excluding zero; five seeds by five folds.
Proposed in design version 0.4 and to be confirmed here: the pathway-wise split is the Human-GEM
subsystem hold-out (pathway subsystems only; transport, exchange and isolated bins are not
pathways), scored on pooled held-out predictions and paired with the permuted-label run of the
same split; the curated modules hold 37 perturbations on the monogenic slice and serve the
interpretability overlap instead.

Failure criterion: if the proposed model does not beat the popularity and network proximity
baselines under the disease-cluster grouped split, the architecture is abandoned before ablations
run. (Design version 0.3 named the gene-wise split; it leaks through genes sharing a disease,
design section 6.1.)

Negative controls run with every split: labels permuted within degree strata; degree-preserving
rewiring; grade shuffle; peripheral-event control.

Precision of the primary comparison on the monogenic slice (451 perturbations, disease-cluster
grouped split): the paired-bootstrap 95 percent half-width of a pooled macro AUPRC difference between
two models is about 0.015 to 0.02 (docs/phase3_main_model.md), so the minimum difference of 0.05 is
detectable there; on the curated-module hold-out (37 perturbations) it is not, which is why the
subsystem definition of the pathway split is proposed above. The full data (drug perturbations and the
signaling layers) will shift these widths and the statement is to be recomputed before posting.

To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; number of flux samples per gene; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
