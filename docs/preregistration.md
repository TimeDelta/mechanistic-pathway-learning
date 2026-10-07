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

Amendment of 7 October 2026, written after the monogenic-slice results were seen. Nothing here was posted
before those results, so the amendment is disclosed as made after them (Nosek et al. 2018,
doi:10.1073/pnas.1708274114, on preregistration when "the data are preexisting").

- The monogenic slice is a pilot. It has one layer of the graph (the metabolic layer), 451 genes, no brain
  expression and no drug perturbations, and on it even a readout trained on the random walk's own distributions
  does not beat the untrained walk (docs/trained_diffusion_readout.md), so the slice cannot separate a weak mechanism
  from too few labels. Its result, that no trained configuration beats the random walk, does not trigger the failure
  criterion above, and no component (the linear-response encoder, the pathway modules, the node descriptors) is
  dropped on its strength. This is the user's decision of 7 October: the mechanistic encoder and the node properties
  it reads stay until the full graph with all layers and better (not more) training examples has been tested.
- The confirmatory test is the full graph with every layer built at the time of the run (graph_full_neuronal on
  7 October: Human-GEM, OmniPath, CollecTRI and the Reactome neuronal and oxidative layer), the node descriptors and
  the brain region and cell-class expression, under the disease-cluster grouped split. The failure criterion applies
  there unchanged.
- What makes a training example "better" is fixed and dated in this file before the first full-graph model is
  scored, from evidence properties only (grade, frequency qualifier, provenance, layer contact), never from any
  model's predictions.
- The ablations (node descriptors, brain expression, edge signs, relation typing) run whatever the primary result,
  because they are what says which part carries the signal (confirmed by the user on 7 October). The with and without
  descriptors comparison is one of them; the descriptors are not dropped before the final experiment.
- The confirmatory configuration is the full graph with the node descriptors, the brain region and cell-class
  expression, the cell-class propagation channels (including a dopaminergic class) and the better training examples.
  Full-graph runs made before all of these exist are full-graph pilots and are reported as such.

To fill in: Phase 1 counts per symptom and grade (docs/phase1_counts.md); final symptom set after
go/no-go; B5 language model and prompt; number of flux samples per gene; power statement for the
GWAS enrichment test; the open questions 8 to 10 of design section 11 (frequency as weight or target;
subsystem versus curated pathway split; the deterioration-term exclusion).
