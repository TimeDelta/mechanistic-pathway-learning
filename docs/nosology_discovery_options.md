# Options for post-experiment nosology discovery (draft for the owner), 17 October 2026

A draft: options put to the owner, nothing decided or built. The owner, 16 October 2026: "please consider my
auto-psych-nosology GitHub repo and how to incorporate it with modifications for the post-experiment nosology
discovery. I'm open to anything even if it means not using that project at all. I just want options for that part".
Decision 9.64 in `docs/open_decisions.md` is the log entry.

Sources. auto-psych-nosology was read at its single public commit (7706137, 9 October 2026) in a read-only clone; the
file and line references below are to that commit. A research report on nosology discovery, delivered to the owner on
16 October, is cited as "the research report". Literature was read through the PubMed API or arXiv on 17 October 2026.

## 1. What the experiment will hand over

The fitted model is the noisy-OR pathway module head over a physiology-graph encoder
(`models/noisy_or_pathway_module_model.py`): K modules, each with a sparse support over graph nodes (one hard concrete
gate per module and node), an activation in [0, 1] per perturbation, and

    P(symptom s | perturbation p) = 1 - (1 - leak_s) * prod_k (1 - link_{k,s} * activation_k(p)),

with separate link matrices for "induces" and "relieves". So three objects can carry a nosology: the symptom
constructs (about 38 on v5, with their links to modules), the perturbations (gene alleles and drugs, with their module
activations) and the modules themselves (each a set of physiology nodes). The design already says that "The
decomposition into modules is not identifiable (permutation, split and merge); stability selection across seeds and
folds reports consensus supports with their selection frequency, and only consensus modules are interpreted"
(`docs/experiment_design.md`, section 5.5), and K is chosen by description length (section 5.5;
`models/minimum_description_length_selection.py`).

Two limits to state first. There are no patients here: a nosology from this model classifies symptom constructs and
causal mechanisms by shared physiology, not people, so it is a mechanistic dimension structure that a patient-level
study would have to test. And no full-graph model has been scored yet; on the slice, runs before the gates got their own
learning rate ended where "no symptom had an active module" (section 5.5, `docs/b6_module_diagnosis.md`). Every option below assumes the confirmatory runs produce stable modules; if they do
not, the nosology question waits.

## 2. What auto-psych-nosology is

Its aim, in its README: "a proof-of-concept approach for automated, data-driven psychiatric nosology" that "tests
whether a dimensional framework can be derived independently of existing diagnostic systems" (README.md:63-64). It
builds a literature knowledge graph from iKraph (README.md:195; create_graph.py), scores psychiatric relevance by
keyword and ontology rules (psychiatry_scoring.py) and partitions graph nodes:

- **RGCN-SCAE** (self_compressing_auto_encoders.py, train_rgcn_scae.py): a relational GCN whose nodes are softly
  assigned to up to K = min(N, max(8, round(sqrt N), 4R)) prototypes (245 in the reported run) by Sinkhorn
  iterations, with a block decoder and L0 gates on clusters. Reading the code and the plots, problems that matter for
  any reuse: every decoder logit lies in [b_r, b_r + 1], because assignments sit on the simplex and the block weights
  in [0, 1]; the entropy floor equals log K, the maximum, so the hinge pushes every assignment toward uniform; the
  cluster gates never closed in the main run (the active-cluster plot is flat at 245); and same-type negative
  sampling does not exclude true edges.
- **"SBM"** (create_sbm_partitions.py): despite the name, "Vanilla SBM partitions via Louvain/greedy modularity
  heuristics." (line 1), NetworkX Louvain with description length reported as minus the modularity; there is no
  stochastic block model and no graph-tool. The research report's advice to reuse "its SBM wrapper" therefore does
  not hold: there is none to reuse.
- **Sparse reduced-rank logistic regression** (sparse_rr_logistic.py, run_sparse_rr_baseline.py,
  docs/sparse_rr_logistic.md): P(Y_ps = 1) = sigmoid(b_s + [ZB]_ps) on frozen perturbation features Z, with group-lasso,
  nuclear-norm and ridge penalties, fitted by proximal gradient, with masks, leakage groups and snapshot pins; the SVD
  of B gives symptom loadings. It was built for this project's kind of input and is tested; per its README, "fitting
  the final research graph remains pending" (README.md:543).
- **Alignment and stability tools** (align_partitions.py, partition_stability.py, cluster_bimodality_probe.py):
  NMI, AMI, ARI, enrichment against inferred HiTOP and RDoC labels (the Benjamini-Hochberg step lacks its cumulative
  minimum), bootstrap coherence, and partition agreement across checkpoints.
- **A physiology adapter** (adapt_physiology_graph.py, relation_schema.py) that already reads this project's
  nodes.parquet and edges.parquet; it rejects the v5 graph as written, whose 376 drug-target edges carry signs of +-0.5,
  and it reduces signs to relation categories ("Signs are encoded as categories, not enforced as signed propagation.",
  adapt_physiology_graph.py:227).

Its own results, in its words: they "fail to achieve meaningful global alignment or robust, non-degenerate structure
under perturbation" (README.md:849), and "Detectable signal exists, but it is fragmented and lacks the connectivity
required for robust global structure." (README.md:983).

The main mismatch: both of its partitioners group **graph nodes by adjacency**. A nosology from this experiment groups
**symptoms, perturbations or modules by what the fitted model says they do**. Run on the physiology graph, they would
return groups of genes, proteins and metabolites, and the symptoms are not nodes of that graph at all.

## 3. The options

**(A) The modules are the dimensions (no new model).** Read the consensus modules of the confirmatory runs as the
nosology's dimensions: each symptom's loading on a module is its link, each perturbation's position its activation, each
module's content its support of physiology nodes. Modules are matched across seeds and folds by their link profiles
and supports (the stability selection section 5.5 already requires), and only consensus modules are read. Symptoms are
then grouped by their link profiles and perturbations by their activation profiles, which gives the transdiagnostic
reading directly: two symptoms belong together when the same mechanisms raise them. Cost: matching and reporting code,
small. Risk: the dimensions are only as good as the modules' stability, and noisy-OR links are non-negative, so a
module that lowers a symptom is not represented in the induces matrix (it is in the relieves matrix, read apart).

**(B) A hierarchy by a nested stochastic block model on the fitted model.** Build the weighted tripartite graph
perturbation - module - symptom from (A)'s objects and fit a nested, degree-corrected, weighted stochastic block model
(graph-tool), which selects the number of groups at every level by description length. Peixoto (2014, Phys. Rev. X 4,
011047, doi 10.1103/PhysRevX.4.011047): "Here we construct a nested generative model that, through a complete
description of the entire network hierarchy at multiple scales, is capable of avoiding this limitation" (the resolution
limit); weights per Peixoto (2018, Phys. Rev. E 97, 012306): "We give a comprehensive treatment of different kinds of
edge weights". New code, a thin wrapper; graph-tool's licence and installation are to be checked. Risk: the hierarchy is
fitted to estimated weights, so their uncertainty (across seeds and folds) must enter as resampling, not be ignored.

**(C) auto-psych-nosology's sparse reduced-rank logistic regression as an independent cross-check.** Use
sparse_rr_logistic.py as it is: features Z frozen from the fitted encoder (module activations or field summaries)
computed on the training folds, the symptom labels of the same folds, the same leakage groups, never the lockbox. Its
symptom loadings (the SVD of B) are a convex, rank-selected second estimate of the symptom dimensions, against which
(A)'s modules can be compared. This is the part of that repository built for this input. Risk: it re-fits the labels,
so it is a check of (A), not an independent truth.

**(D) Consensus clustering with a null, as a guard.** Any clustering in (A) to (C) is resampled and compared with a
reference null, since Șenbabaoğlu et al. found that consensus clustering "is able to divide randomly generated unimodal
data into apparently stable clusters for a range of K" (Sci. Rep. 2014, 4:6207, doi 10.1038/srep06207), and M3C
"simulates null distributions of stability scores for a range of K values thus enabling a comparison with real data to
remove bias and statistically test for the presence of structure" (John et al., Sci. Rep. 2020, 10:1816, doi
10.1038/s41598-020-58766-1). partition_stability.py and cluster_bimodality_probe.py can carry the resampling; the null
is new.

**(E) Not recommended: RGCN-SCAE or the Louvain script on the physiology graph, and selection by agreement with HiTOP
or RDoC.** The first two partition the wrong object and carry the problems listed in section 2; hyperband_partition.py
selects by HiTOP and RDoC agreement, which would make the existing nosology the target, the circularity the owner named
on 16 October (decision 9.63: "I think it could introduce circularity from existing nosology").

## 4. Validation, common to every option

- Held-out perturbations: the structure is read on the development folds and tested once on the lockbox (does a held-out
  perturbation's predicted profile fall in the dimension its mechanism puts it in?).
- Stability: refits across seeds and leakage-group resamples, with the permuted-label and rewired-graph controls the
  registered runs already carry; a dimension that survives the rewired graph is a degree effect, not physiology.
- External structure as validation only, never as a selection criterion: the four genetic factors of Grotzinger et al.
  ("neurodevelopmental, compulsive, psychotic and internalizing", Nat. Genet. 2022, 54:548-559, doi
  10.1038/s41588-022-01057-4), the three groups of the Cross-Disorder Group ("Genetic correlation analyses revealed a
  meaningful structure within the eight disorders, identifying three groups of inter-related disorders", Cell 2019,
  179:1469-1482, doi 10.1016/j.cell.2019.11.020), HiTOP, and DSM or ICD. These are disorder-level and genetic; the
  model's units are symptom constructs and mechanisms, so the comparison needs a stated map from constructs to their
  disorders, fixed before the comparison is run.

## 5. Recommendation and what it would take

(A) as the primary reading, since the modules are what the experiment already fits and section 5.5 already requires
their consensus; (C) as the cross-check, reusing auto-psych-nosology's sparse reduced-rank regression as it is; (B) for
the hierarchy once (A)'s modules are stable; (D) around every clustering. From auto-psych-nosology: reuse
sparse_rr_logistic.py, the metric functions of align_partitions.py (with the Benjamini-Hochberg step corrected),
partition_stability.py and cluster_bimodality_probe.py; leave the RGCN-SCAE, the Louvain script and the iKraph builder.
The adapter is not needed for (A) to (D); if it is used, it must accept fractional signs and carry them as magnitudes.

## 6. Questions for the owner

1. Which unit should the nosology classify first: symptom constructs (by shared mechanisms), perturbations (by
   mechanism class), or both together through (B)?
2. Is (A) as the primary reading right, with (C) as its check, or should (C) or (B) lead?
3. Should the external structures of section 4 be compared at all before the lockbox is read, or only after?
