# Building a signed physiology graph: published practice read against this build (record, 17 October 2026)

A record of a research pass, not a plan: nothing here is decided. The owner, 16 October 2026: "Can you please also
research best practices for this type of knowledge graph creation and how they might be best incorporated?" A
research agent read the literature and the build and measured the v5 graph; the session that started it checked the
result (next section) and put the ranked changes of section 4 to the owner as decision 9.70 in `docs/open_decisions.md`.

What the pass read: the nine documents the task names; `docs/edge_gate_simulation.md`; decisions 9.21, 9.22, 9.42, 9.57, 9.62 and 9.63; and the docstrings of every module in `mechanistic_pathway_learning/graph/`. Counts are on `graph_full_neuronal_split_binders_drugs_v5` unless marked v4. They come from the pinned raw files and from two OmniPath downloads of 17 October with `fields=evidences`, read with `experiments/read_public_page.py`. Articles were found and read with the PubMed tools; their DOIs are linked in the Sources. Nothing in the repository was changed.

## How this record was checked

By the session that started the pass, on 17 October 2026, before it was written here.

- **M2, recounted** from the pass's download of OmniPath's CollecTRI dataset with `fields=evidences` (64,516
  interactions): CollecTRI's own `sign_decision` is "default activation" for 34,993, "PMID" for 27,632 and "regulon"
  for 1,892. One interaction carries a list of two decisions ("PMID" and "default activation"), which the pass's count
  dropped, so its figures below are one lower in the first two places. Nothing that follows depends on the difference.
- **M3, recounted exactly**: 114 PubMed ids each behind more than 100 signalling or transcription edges, 48,887 edges
  citing one of them, the largest 34199813 (13,735 edges) and 17570479 (7,166), and 4,639 of the 38,435 labels with
  fewer than two records once they are set aside.
- **Quotes**: 84 of the 86 quoted passages below were found verbatim in a source the pass read (its PubMed abstracts
  and full texts, the pages it saved and its downloads), with the pass's own drafts excluded from the search. Of the
  other two, "A live service with no release tag" is a comment of `scripts/fetch_raw_sources.sh` (lines 81 and 82),
  and the pass's quotation of the Biolink aspect qualifier was not found verbatim in any saved page; the bullet that
  carried it now quotes the Biolink page itself, read on 17 October 2026 with `experiments/read_public_page.py`
  (SHA-256 0c843c94b7a3ce2e...).
- **Not recounted**: M1, M4 and M5 and the label counts of section 4, item 1. They rest on the pass's scripts over the
  same downloads.
- **One statement corrected.** Section 1 said Reactome complexes are nodes "with a minimum over members". The minimum
  is the message-passing encoder's alone (`member_of` is its conjunction relation), and it applies to Reactome's sets as
  well as its complexes, although a set is a list of alternatives: the pinned SBML says of each, "Derived from a
  Reactome DefinedSet. This is a list of alternative entities, any of which can perform the given function". On v5,
  257 entity nodes with `member_of` edges are sets (210 DefinedSets and 47 CandidateSets, 1,552 edges) and 763 are
  complexes (4,708 edges). The linear-response encoder sums the members under its divisor. Both readings are recorded
  in `docs/complex_conjunction_semantics.md` and `docs/physiological_adherence_review.md`.
- The sources are in `docs/references.bib` (decision 9.70), for the owner's reference manager.

## 1. What the build already does

- Pins with sizes, SHA-256 and quoted licences (`docs/data_sources.md`); release manifests (`mechanistic_pathway_learning/data_release.py`).
- Per-edge family, evidence kind, distinct PubMed ids, resources, iCite species scores and IntAct absences (`graph/edge_evidence.py`; `docs/edge_priors.md` sections 2 and 11). Counting distinct PubMed ids is CollecTRI's rule: "PMIDs can only count as evidence for a TF–gene interaction once, even if they were curated by several databases" [1]. It is stricter than OmniPath's curation effort, which pypath counts over resource and reference pairs, "(self.a, self.b, res, ref)" [2].
- Text-mined layers kept apart. Priors act as per-edge costs, tested in planted worlds, with a permuted-prior control and degree-preserving rewiring within each relation (`experiments/run_baselines.py`). That is Hetionet's control: "Hetnet permutation preserves node degree while eliminating edge specificity" [3].
- Curated sign corrections with a drift check (`load_curated_edge_sign_overrides`, 9.10a), sourced receptor to G-alpha coupling (9.23), no family-level inference (9.22).
- Compartments linked only by transport reactions, and a mass and charge audit (`docs/physiology_audit.md`). Currency tags by name and degree and per-edge carrier rules (`tag_currency_metabolites.py`, `cofactor_edges.py`) follow currency-metabolite partitioning [4] and main-compound pathfinding [5].
- UniProt entries split only where data tell them apart (`docs/gene_protein_split.md`); Reactome complexes as nodes, read with a minimum over members by the message-passing encoder only (corrected; see the checks above).
- A cited validation set (`configs/physiology_validation_set.csv`), a degree-scaled popularity baseline (`models/baselines/popularity_baseline.py`) and sign ablations (`models/linear_response_encoder.py`).

The gaps are in what an evidence record is attached to, and in representation choices made before any prior is set.

## 2. Measured for this report

M1. A pair's evidence is not its sign's evidence. OmniPath records positive, negative, directed and undirected support separately for each interaction. Rebuilt from that field, the build's record count matches `edge_evidence.parquet` for 84,313 of 84,945 matched signed signalling edges. Decision 9.63 takes 14,680 signed signalling edges as known. Of these, 4,924 have two or more PubMed ids behind the stored sign, 5,906 have one and 3,727 have none (123 unmatched). Where there is none, the sign comes mostly from pathwayextra, "directed and signed interactions, without literature references" [6]. The resources are the one OmniPath names Wang (3,212 edges), KEGG (393), Cui2007 (289) and KEGG-MEDICUS (274). Of all 87,402 signed signalling edges, 37,885 have no PubMed id behind the sign.

M2. Half of CollecTRI's signs are a default. Its paper says: "For the remaining 21730 TF–gene interaction an activating mode was assigned by default" (of 43,175). The reason given is that "a default activating mode of regulation outperformed a default repressing mode" in TF-activity benchmarks [1]. In OmniPath's CollecTRI dataset, the CollecTRI entry's `sign_decision` is "default activation" for 34,992 of 64,516 interactions, "PMID" for 27,631 and "regulon" for 1,892. Of the default-signed interactions, 22,687 have no CollecTRI2 entry, and the two entries state opposite signs for 4,489 interactions. The build does not fetch the field (`scripts/fetch_raw_sources.sh`) and turns 1,280 tied rows into +1 (`build_transcription_layer`). Of the 41,638 transcription edges, 11,481 rest only on default-signed interactions without a CollecTRI2 entry. Another 9,851 are default-signed and have a CollecTRI2 entry whose basis was not resolved. Of the 15,076 transcription labels, 10,779 are PMID-signed, 417 default only, 3,259 default with a CollecTRI2 entry, 619 regulon-signed or mixed and 2 other. A gate never flips a sign (`docs/edge_priors.md` section 8), so a wrong sign on a label is permanent.

M3. Papers behind thousands of edges count as one curated record each. 114 PubMed ids each support more than 100 signalling or transcription edges, and 48,887 edges cite at least one of them. The largest:
- a 2021 review of ubiquitin-conjugating enzymes (13,735 edges) [7]
- the paper introducing the NetworKIN predictor (7,166 edges, 824 labels) [8]
- a cryo-EM structure (3,175)
- two database papers (1,824 and 1,789)
- the BioPlex 3.0 affinity-purification map (1,455) [9]
- an ATM and ATR substrate screen of "over 700 proteins" (1,385, including 1,311 of ATM's 1,677 binds edges) [10]

Not counting the 114, 4,639 of the 38,435 labels would have fewer than two records.

M4. Complexes become Cartesian products. 29,593 of 118,116 OmniPath edges and 2,093 CollecTRI edges exist only because a heteromeric complex was expanded into all its member pairs (1,014 and 582 labels). On v4, SIGNOR's ubiquitin–E2 complexes gave UBA52, UBC, UBB and RPS27A 423 to 476 activates edges each.

M5. Direction. All 30,714 OmniPath binds edges are reciprocated; of the 15,816 unsigned raw rows, only 1,142 have a reverse row. 18,084 raw rows have `consensus_direction` False; those with a sign (1,717 stimulation and 640 inhibition flags) enter as directed edges. SIGNOR marks 5,743 of 43,623 rows indirect; 455 signed graph edges are indirect-only there (v4). Neither field is read.

M6. Versions. The only release, v0.4 of 6 October, holds `graph_full` (33,964 nodes), not the v4 or v5 graph the models read (`docs/status.md`). OmniPath is "A live service with no release tag" (`scripts/fetch_raw_sources.sh`), and `docs/data_sources.md` records no SHA-256 for its files.

## 3. Practices by topic

Each topic ends with what the build does now, the change, its benefit for the physiology, the cost and the risk.

### 3.1 Standards and frameworks

Biolink [11] is a schema of classes, predicates and edge attributes. Several of its terms bear on this graph:
- Knowledge level separates a `knowledge_assertion`, "A statement of purported fact that is put forth by an agent as true, based on assessment of direct evidence", from a `prediction`, "A statement of a possible fact based on probabilistic forms of reasoning over more indirect forms of evidence" [12].
- Agent types include `text_mining_agent` [13].
- `primary_knowledge_source` is "The most upstream source of the knowledge expressed in an Association that an implementer can identify" [14].
- `object_aspect_qualifier`: "The purpose of the aspect slot is to indicate what aspect is being affected in an 'affects' association. This qualifier specifies a change in the object of an association", with the example values stability, abundance, expression and exposure [15].

ECO "allows tracking of annotation provenance, establishment of quality control measures, and evidence-based data mining" [16]. PSI-MI [17] gained "the new term 'causal interaction', and related child terms" in CausalTAB [18]. MI2CAST adds "a comprehensive set of other contextual details" for causal statements [19]. MIscore weights "MI:0407 | direct interaction" 1.00, "MI:0915 | physical association" 0.66 and "MI:0914 | association" 0.33 [20]. An affinity-purification map reports association, not binding.

BioCypher aims to "transparently build biomedical knowledge graphs while preserving provenances of the source data" [21]. KGX checks that "nodes are categorized using Biolink classes, edges are labeled using valid Biolink relationship types" [22]. SPOKE "is built weekly by python scripts which download each resource, check for integrity and completeness" [23]. Monarch [24] and Translator [25] build on Biolink [11]. These are mostly unsigned association graphs, and some of their cleaning steps are wrong for a causal graph: PrimeKG "cleaned the graph by dropping NaN and duplicate edges, adding reverse edges" and kept "its largest connected component" [26]. FAIR asks that "(Meta)data are assigned a globally unique and persistent identifier" and "are associated with detailed provenance" [27, 28].

- Now: the evidence kind stands in for agent type, for a whole edge. There is no knowledge level, aspect, ECO term or PSI-MI term. Node keys are symbols (`GENE:<symbol>`), and on v4, 56 gene nodes are keyed on ChEBI, UniProt or Ensembl ids.
- Change: record knowledge level, agent type, primary and aggregator source and aspect per component (presence, direction, sign), mapped to ECO and PSI-MI. Key nodes on HGNC, UniProt and ChEBI identifiers.
- Benefit: a CollecTRI default (M2) or a NetworKIN edge (M3) becomes a prediction rather than an assertion. SIGNOR's "up-regulates activity" and "up-regulates quantity by expression", now merged in `activates`, stay apart.
- Cost: a few columns, filled from the evidences JSON.
- Risk: mapping judgement, logged per resource.

Migrating to BioCypher or KGX would not change the physiology and is not recommended now.

### 3.2 Confidence, curation error, bias and negatives

How the resources score confidence:
- STRING: its "benchmarked confidence scores ... generally correspond to the probability of finding the linked proteins within the same KEGG pathway", while curated imports get "a uniform confidence score per dataset" [29]. Since 2023, confidence is "estimated based on the detection method used" [30].
- SIGNOR: the score "combines four features" and "can never be 0 since all SIGNOR relations are supported by at least one reference" [31].
- DoRothEA: levels run "from A (highest confidence) to E (lowest confidence)". A means "supported in at least two curated resources", E "uniquely supported by computational predictions", and there is "an inverse relationship between the coverage and the accuracy" [32].

Only STRING, DoRothEA and CollecTRI calibrate their scores against a reference, and STRING's reference is co-membership of a pathway, not a signed edge.

The evidence on curation error is old and concerns binding:
- Recuration found curation "error-prone and possibly of lower quality than commonly assumed". Of singly supported human interactions, "45% of interactions were not validated"; "The most common errors here were wrong species and wrong protein name"; and "40% of the human LC-multiple interactions were no longer supported by multiple publications after recuration" [33].
- Two databases curating the same paper "fully agree on 42% of the interactions" [34].
- Three pathway databases agreed on 0 to 61 percent of gene pairs [35].
- 63 percent of INDRA's statements rest on "a single sentence or database entry", and "even high‐mention Statements do not approach an accuracy of 100%" [36].

Two PubMed searches and a web search found no measured error rate for current curated causal signalling resources.

On bias:
- "multifunctionality, rather than association, is a primary driver of gene function prediction" [37]
- "functional information within gene networks is typically concentrated in only a very few interactions" [38]
- interaction data carry "selection bias in favor of biologically interesting bait proteins" [39]
- annotation bias "leads biomedical researchers to focus on richly annotated genes instead of those with the strongest molecular data" [40]
- in graphs, "certain entities can get topologically overrepresented" [41]
- a systematic map "appears strikingly more homogeneous" than the literature [42]

On negatives, Negatome 2.0 is "derived by manual curation of literature and by analyzing three-dimensional structures of protein complexes" [43].

- Now: the tiers are set by hand and, in planted worlds, "weigh as much as the data" (`docs/edge_priors.md` section 10); 9.63 (d) proposes estimating them. IntAct absences mark 54 edges as contested. There is a degree-scaled popularity baseline but no covariate for literature attention on the edge side.
- Change: calibrate the tiers against perturbations and recuration (3.6), not against a curated pathway set. Down-weight bulk papers (M3). Add Negatome as a second absence source for the audit only, since both absence sets describe binding, not regulation.
- Benefit: tier values that are measured rather than assumed.
- Cost: moderate.
- Risk: a recuration sample small enough to be cheap gives wide intervals.

The human-evidence rule of 9.63 targets the most common error Cusick and colleagues name. But the loss of multiple support on recuration, together with M3, makes two records weaker evidence than they look.

### 3.3 Signed and directed causal networks

Sign-consistency methods "test the consistency of network topologies and experimental data" and "apply repair operations to the network model suggesting missing or wrong interactions" [44]. When ranking regulators, causal reasoning over literature graphs is "extremely robust to random noise and to missing or spurious information" [45]. CausalR [46], CARNIVAL [47] and COSMOS [48] fit signed prior networks to omics data. The benchmarks:
- networks were scored "in terms of causal validity with respect to unseen interventional data", and "incorporating known biology was generally advantageous" [49]
- "many interactions in literature-derived networks do not appear to be functional in the liver cells" that were fitted [50]
- "the choice of network and algorithm has a profound impact" [51]
- a review names "the general lack of ground truth and the limitations of prior knowledge" [52]
- in kinase inference, "the evidence supporting the interactions has a strong effect on the predictions" [53]

- Now: signs meet data only through symptom labels, about 0.05 label-grade positives per gated edge (`docs/edge_gate_simulation.md` section 4). Signalling ties become binds and transcription ties become +1. In SIGNOR, 263 of 30,461 signed pairs carry both signs, 122 of them with a different recorded cell context or mechanism. SIGNOR draws "direct interactions ... as solid lines and indirect ones as dashed lines" [54], and "Indirect" is one of Cusick's error classes [33]. In a propagation model, an indirect edge beside its mediating path counts one influence twice.
- Change: audit sign consistency against TF knockdowns and against phosphoproteomics after kinase perturbation. Keep dual-sign edges with their context, and flag indirect edges where the mediating path exists. KnockTF, CollecTRI's benchmark source, answered HTTP 502 on 17 October; that is a failure, not a robots refusal. The second source is Perturb-seq "targeting all expressed genes with CRISPR interference (CRISPRi) across >2.5 million human cells", redistributed by scPerturb [55].
- Benefit: measured sign precision per tier.
- Cost: moderate to high.
- Risk: K562 cells are not neurons, and knockdowns include indirect effects.

### 3.4 Context specificity

Tissue networks "predict lineage-specific responses to perturbation" [56], and "network edges (transcription factor to target gene connections) have higher tissue specificity than network nodes (genes)" [57]. A brain atlas of "107 distinct cell types and their specific utilization of 544,735 candidate cis-regulatory DNA elements" [58] could gate transcription edges per cell class.

For metabolism, tINIT, iMAT and GIMME extract functional subnetworks [59, 60, 61]. A comparison of such methods found that "model extraction method choice had the largest impact on model accuracy" and "Model content varied substantially across different parameter sets" [62]. ftINIT "contextualized models from 202 single-cell clusters across 19 human organs" [63]. MEMOTE tests "four general areas: annotation, basic tests, biomass reaction and stoichiometry" [64]. Human-GEM reports "100% stoichiometric consistency, 99.4% mass-balanced reactions, and 98.2% charge-balanced reactions" [65]. Energy-generating cycles "occur in over 85% of metabolic models without extensive manual curation" [66].

- Now: cell-class weights act on nodes, and regulons are kept for brain-expressed factors; neither can express specificity at the level of an edge. The continuous weights avoid thresholds, which is defensible given the sensitivity reported above, but lack the network consistency that task-driven extraction enforces. Mass and charge balance is audited on the Human-GEM SBML only. The Reactome, manganese and carrier reactions the build adds are not checked, and that is where the chloride defect was (`docs/missing_physiology_review.md` section 1).
- Change: run MEMOTE-style consistency, balance and blocked-reaction tests over the whole metabolic layer, at low cost. Later: use ftINIT reaction sets as a check of the weights, and add edge-level context for transcription.
- Cost and risk of the later steps: high cost, and enhancer-to-gene links are themselves predictions.

### 3.5 Identifiers and versions

MetaNetX is "a reconciliation of metabolites and biochemical reactions providing cross-links between major public biochemistry and Genome-Scale Metabolic Network (GSMN) databases" [67]. Bioregistry starts from "The standardized identification of biomedical entities is a cornerstone of interoperability" [68]. Complex Portal holds "over 1100" curated human complexes [69]. Curated sources diverge on "organism or splice isoforms" and on "alternative representations of multi-protein complexes" [34].

- Now: Human-GEM's pinned `metabolites.tsv` already carries `metMetaNetXID` and `metChEBIID`. Even so, the 60 metabolite nodes minted from Reactome ChEBI ids (v4) include copies of Human-GEM species (H+, 12(S)-HETE, DAG). M4 and M6 stand.
- Change: merge metabolites through MetaNetX, map OmniPath complexes to Complex Portal or Reactome complexes, and release v5 with a manifest.
- Cost: low to moderate.
- Risk: the build keeps some pools apart on purpose (the neuronal and redox pools), so a merge needs an exclusion list.

### 3.6 Validating the built graph

Linked-data quality falls into "Accessibility (ii) Intrinsic (iii) Contextual and (iv) Representational" groups, and "Semantic accuracy is defined as the degree to which data values correctly represent the real world facts" [70]. Semantic accuracy is the owner's priority. Network benchmarks reward volume: "performance scales with network size", although "Correcting for size, we find that the DIP network provides the highest efficiency (value per interaction)" [71].

Held-out edge recovery suits link predictors, not gates that can only close edges. What fits here is Cusick's method of sampled recuration per tier: read the cited papers for species, direction, sign and directness, and report precision with intervals. Stratifying by each gene's literature attention [40] separates the graph's accuracy from publication volume. The same sample is the reference that 9.63 (d) needs.

## 4. Changes, ranked

1. Attach evidence to each component of an edge. Keep OmniPath's positive, negative and directed evidence and CollecTRI's `sign_decision`, and apply 9.63's rule (two records, one of them human) to the stored sign (M1, M2). 3,727 signed signalling labels carry a sign that no cited paper states. 4,295 transcription labels carry a sign that CollecTRI set by default or by the factor's regulon role. Under a sign-level rule the labels fall to at most 4,924 and 10,779; that is the owner's call.
2. Treat default, regulon-role and tied signs as unknown, using a sign prior or an unsigned relation. A gate cannot flip them, and CollecTRI chose its default for inference accuracy, not physiology.
3. Count a paper behind more than a set number of edges (a review, database paper, prediction or screen) as one record of its own kind. At 100 edges this removes 4,639 labels (M3).
4. Represent heteromeric complexes as complex nodes with the member_of minimum, as Reactome's already are, instead of as all member pairs (M4).
5. Stop symmetrising unsigned directed rows. Carry `consensus_direction` and SIGNOR's directness flag, and flag indirect edges whose mediating path exists (M5).
6. Measure sign precision per tier by perturbation checks and sampled recuration, and feed it to 9.63 (d). Without this every tier value stays assumed.
7. Record knowledge level, agent type, primary and aggregator source and aspect per component, with ECO and PSI-MI terms.
8. Release v5 with `edge_evidence.parquet` and a manifest, and record the SHA-256 of the OmniPath files and the evidences JSON (M6).
9. Run MEMOTE-style tests over the whole metabolic layer and merge the minted metabolites through MetaNetX.
10. Later: edge-level cell-class context for transcription, and ftINIT as a check of the expression weights.

## Notes on sources

- No permission request was withdrawn, in this run or in the earlier one.
- No page text addressed to an AI assistant was found. A fetch summary of an OmnipathR commit patch reported "a block of instructions" that repeated the request's own wording. It was read as the tool echoing the prompt, and nothing was acted on.
- Licences stated in this report, with where each was read:
  - ECO "is released into the public domain under a CC0 1.0 Universal license" (abstract of [16]).
  - Complex Portal: "Our data license is now CC0 to encourage data reuse" (abstract of [69]).
  - Bioregistry: "source code and data are available under the MIT and CC0 Licenses" (abstract of [68]).
  - scPerturb's Zenodo record 7041849, which holds `ReplogleWeissman2022_K562_gwps.h5ad`, gives "Creative Commons Attribution 4.0 International" on its record page (read with the script). The pertpy and scPerturb pages state no licence, and the original deposit was not read. Zenodo's API was refused by its robots file and was left unread.
- Files written outside the repository for reuse, in the session's scratchpad (not kept with the repository). Each download is named by its SHA-256 prefix, with a JSON record of its address, time and full digest beside it:
  - /tmp/claude-0/-home-claude/9f7f0165-a20e-58d4-b29a-868e7d449496/scratchpad/dl_collectri_evidences/60b79abc6a0779ca.txt
  - /tmp/claude-0/-home-claude/9f7f0165-a20e-58d4-b29a-868e7d449496/scratchpad/dl_omnipath_evidences/e121c5665fd8c824.txt
  - the counting scripts, in /tmp/claude-0/-home-claude/9f7f0165-a20e-58d4-b29a-868e7d449496/scratchpad/myscripts/
- The live service was read 14 days after the pinned files were fetched. Every pinned CollecTRI row matched, as did 84,945 of 87,402 signed signalling edges.

## Sources

PubMed-retrieved articles carry their DOI links. Web pages marked "script" were read with `experiments/read_public_page.py` on 17 October 2026.

1. Müller-Dott et al., CollecTRI regulons, Nucleic Acids Res 2023. https://doi.org/10.1093/nar/gkad841
2. pypath, `pypath/core/interaction.py` (script; SHA-256 6f565e82c5c3ecd2…). https://raw.githubusercontent.com/saezlab/pypath/master/pypath/core/interaction.py
3. Himmelstein et al., Hetionet, eLife 2017. https://doi.org/10.7554/eLife.26726
4. Huss and Holme, currency metabolites, IET Syst Biol 2007. https://doi.org/10.1049/iet-syb:20060077
5. Faust et al., RPAIR pathfinding, J Mol Biol 2009. https://doi.org/10.1016/j.jmb.2009.03.006
6. OmnipathR workshop vignette. https://bioconductor.statistik.tu-dortmund.de/packages/3.22/bioc/vignettes/OmnipathR/inst/doc/bioc_workshop.html
7. Bui et al., Cells 2021. https://doi.org/10.3390/cells10061383
8. Linding et al., NetworKIN, Cell 2007. https://doi.org/10.1016/j.cell.2007.05.052
9. Huttlin et al., BioPlex 3.0, Cell 2021. https://doi.org/10.1016/j.cell.2021.04.011
10. Matsuoka et al., Science 2007. https://doi.org/10.1126/science.1140321
11. Unni et al., Biolink Model, Clin Transl Sci 2022. https://doi.org/10.1111/cts.13302
12. Biolink KnowledgeLevelEnum (script). https://biolink.github.io/biolink-model/KnowledgeLevelEnum/
13. Biolink AgentTypeEnum (script). https://biolink.github.io/biolink-model/AgentTypeEnum/
14. Biolink primary_knowledge_source (script). https://biolink.github.io/biolink-model/primary_knowledge_source/
15. Biolink object_aspect_qualifier (script, SHA-256 0c843c94b7a3ce2e...). https://biolink.github.io/biolink-model/object_aspect_qualifier/
16. Nadendla et al., ECO, Nucleic Acids Res 2022. https://doi.org/10.1093/nar/gkab1025
17. Hermjakob et al., PSI-MI, Nat Biotechnol 2004. https://doi.org/10.1038/nbt926
18. Perfetto et al., CausalTAB, Bioinformatics 2019. https://doi.org/10.1093/bioinformatics/btz132
19. Touré et al., MI2CAST, Bioinformatics 2021. https://doi.org/10.1093/bioinformatics/btaa622
20. Villaveces et al., MIscore, Database 2015. https://doi.org/10.1093/database/bau131
21. Lobentanzer et al., BioCypher, Nat Biotechnol 2023, https://doi.org/10.1038/s41587-023-01848-y; quoted from its preprint, https://arxiv.org/abs/2212.13543v1
22. KGX. https://www.github.com/biolink/kgx
23. Morris et al., SPOKE, Bioinformatics 2023. https://doi.org/10.1093/bioinformatics/btad080
24. Putman et al., Monarch, Nucleic Acids Res 2024. https://doi.org/10.1093/nar/gkad1082
25. Fecho et al., Translator, Clin Transl Sci 2022. https://doi.org/10.1111/cts.13301
26. Chandak et al., PrimeKG, Sci Data 2023. https://doi.org/10.1038/s41597-023-01960-3
27. Wilkinson et al., FAIR, Sci Data 2016. https://doi.org/10.1038/sdata.2016.18
28. GO FAIR principles. https://www.go-fair.org/FAIR-principles/
29. von Mering et al., STRING, Nucleic Acids Res 2005. https://doi.org/10.1093/nar/gki005
30. Szklarczyk et al., STRING 2023, Nucleic Acids Res 2023. https://doi.org/10.1093/nar/gkac1000
31. Licata et al., SIGNOR 2.0, Nucleic Acids Res 2020. https://doi.org/10.1093/nar/gkz949
32. Garcia-Alonso et al., DoRothEA, Genome Res 2019. https://doi.org/10.1101/gr.240663.118
33. Cusick et al., Nat Methods 2009. https://doi.org/10.1038/nmeth.1284
34. Turinsky et al., Database 2010. https://doi.org/10.1093/database/baq026
35. Soh et al., BMC Bioinformatics 2010. https://doi.org/10.1186/1471-2105-11-449
36. Bachman et al., INDRA, Mol Syst Biol 2023. https://doi.org/10.15252/msb.202211325
37. Gillis and Pavlidis, PLoS One 2011. https://doi.org/10.1371/journal.pone.0017258
38. Gillis and Pavlidis, PLoS Comput Biol 2012. https://doi.org/10.1371/journal.pcbi.1002444
39. Gillis et al., J Proteomics 2014. https://doi.org/10.1016/j.jprot.2014.01.020
40. Haynes et al., Sci Rep 2018. https://doi.org/10.1038/s41598-018-19333-x
41. Bonner et al., Brief Bioinform 2022. https://doi.org/10.1093/bib/bbac279
42. Rolland et al., Cell 2014. https://doi.org/10.1016/j.cell.2014.10.050
43. Blohm et al., Negatome 2.0, Nucleic Acids Res 2014. https://doi.org/10.1093/nar/gkt1079
44. Thiele et al., BMC Bioinformatics 2015. https://doi.org/10.1186/s12859-015-0733-7
45. Chindelevitch et al., Bioinformatics 2012. https://doi.org/10.1093/bioinformatics/bts090
46. Bradley and Barrett, CausalR, Bioinformatics 2017. https://doi.org/10.1093/bioinformatics/btx425
47. Liu et al., CARNIVAL, NPJ Syst Biol Appl 2019. https://doi.org/10.1038/s41540-019-0118-z
48. Dugourd et al., COSMOS, Mol Syst Biol 2021. https://doi.org/10.15252/msb.20209730
49. Hill et al., HPN-DREAM, Nat Methods 2016. https://doi.org/10.1038/nmeth.3773
50. Saez-Rodriguez et al., Mol Syst Biol 2009. https://doi.org/10.1038/msb.2009.87
51. Hosseini-Gerami et al., BMC Bioinformatics 2023. https://doi.org/10.1186/s12859-023-05277-1
52. Garrido-Rodriguez et al., Mol Syst Biol 2022. https://doi.org/10.15252/msb.202211036
53. Hernandez-Armenta et al., Bioinformatics 2017. https://doi.org/10.1093/bioinformatics/btx082
54. SIGNOR tutorial. https://www.signor.uniroma2.it/tutorial/index.html
55. Replogle et al., Perturb-seq, Cell 2022, https://doi.org/10.1016/j.cell.2022.05.013 ; pertpy data page, https://pertpy.readthedocs.io/en/stable/usage/data/pertpy.data.replogle_2022_k562_gwps.html ; scPerturb (script), https://projects.sanderlab.org/scperturb/ ; Zenodo record 7041849 (script; SHA-256 761e3daabc3cc5e3…), https://zenodo.org/records/7041849
56. Greene et al., GIANT, Nat Genet 2015. https://doi.org/10.1038/ng.3259
57. Sonawane et al., Cell Rep 2017. https://doi.org/10.1016/j.celrep.2017.10.001
58. Li et al., Science 2023. https://doi.org/10.1126/science.adf7044
59. Agren et al., tINIT, Mol Syst Biol 2014. https://doi.org/10.1002/msb.145122
60. Shlomi et al., iMAT, Nat Biotechnol 2008. https://doi.org/10.1038/nbt.1487
61. Becker and Palsson, GIMME, PLoS Comput Biol 2008. https://doi.org/10.1371/journal.pcbi.1000082
62. Opdam et al., Cell Syst 2017. https://doi.org/10.1016/j.cels.2017.01.010
63. Gustafsson et al., ftINIT, PNAS 2023. https://doi.org/10.1073/pnas.2217868120
64. Lieven et al., MEMOTE, Nat Biotechnol 2020. https://doi.org/10.1038/s41587-020-0446-y
65. Robinson et al., Human-GEM, Sci Signal 2020. https://doi.org/10.1126/scisignal.aaz1482
66. Fritzemeier et al., PLoS Comput Biol 2017. https://doi.org/10.1371/journal.pcbi.1005494
67. Moretti et al., MetaNetX, Nucleic Acids Res 2021. https://doi.org/10.1093/nar/gkaa992
68. Hoyt et al., Bioregistry, Sci Data 2022. https://doi.org/10.1038/s41597-022-01807-3
69. Meldal et al., Complex Portal, Nucleic Acids Res 2022. https://doi.org/10.1093/nar/gkab991
70. Zaveri et al., linked-data quality survey, Semantic Web 2016: https://semantic-web-journal.net/node/773 and https://boa.unimib.it/bitstream/10281/96841/5/DQ_Survey.pdf
71. Huang et al., Cell Syst 2018. https://doi.org/10.1016/j.cels.2018.03.001