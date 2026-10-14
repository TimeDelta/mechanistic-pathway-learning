# Superfluous physiology in the graph: what cannot carry a response, what is duplicated, what it costs

14 October 2026, at the owner's question "Do you think there's superfluous physiology in the graph?" This is the
companion of `docs/missing_physiology_review.md`, which asks what is absent. Every count here was measured by
`experiments/measure_superfluous_physiology.py` on `data/processed/graph_full_neuronal_split_binders_drugs_v4`, the v4
drug graph, with `data/processed/evidence_full_v4` as the perturbation set. Three other graphs are named where they
differ: `graph_full_neuronal_split_binders_drugs_v4_enzymes` (identical to v4 but for 23 `clears` edges, measured here),
`graph_full_neuronal` (the confirmatory graph, no drug layer, run with `--graph-dir`) and
`graph_full_neuronal_split_binders` (where the reused `measure_relation_flow.py` was run). Measurements are separated
from inferences, and section 9 states the limit of each. Nothing here changes a graph, a label, a model or a model input.

Every number below is printed by that script, so each is reproducible by re-running it; the two runs behind this review
are `data/processed/label_review/superfluous_physiology.json` (v4) and `superfluous_physiology_neuronal.json`. The
numbers were re-derived from a clean run after a first draft, and that pass corrected four of them: the currency-drop
decomposition in section 2, the count of proteins the dead `encodes` edges dilute in sections 1 and 5 (4,036 that hold a
field, not the 1,949 that hold none), the `binds` count into a gene node on `graph_full_neuronal` in section 4, and a
misquotation of the decay figure in section 9.

## 1. The answer

Yes, and it is one named block rather than the metabolic layer the question implies. 8,729 of 49,796 nodes (17.5%) can
hold no field for any registered perturbation, and 6,013 of them are gene nodes with no transcriptional regulator and no
label; their 6,894 `catalyzed_by` edges sit in the `total_in_degree` divisor of 3,587 reactions that do hold a field and
divide the live message there by a median factor of 1.400, and their 5,998 `encodes` edges sit in the divisor of 4,036
proteins that hold a field and divide it by a median factor of 1.500. Both reach a 90th-percentile factor of 2.000. The
metabolic layer itself is not superfluous on any of the five measurements: 0 of 148 Human-GEM subsystems is
unreachable, 3.1% of metabolites and 3.2% of reactions are, and compartment copies are almost never interchangeable, with
one field-identical pair across 2,483 multi-copy species. The inference I draw is that paring should take the orphan gene
layer and leave the metabolism alone, and that `VM_c` is a separate problem: its divisor of 409 is 408 live edges crowding
each other and 1 dead one, so no amount of paring fixes it.

## 2. What "superfluous" was taken to mean, and why it is harmful rather than inert

Three readings, measured apart.

- Unreachable. A node no perturbation seed can reach holds a field of exactly zero in every example, because the encoder
  is linear in its input and the unperturbed response is zero (`linear_response_encoder.py`: "The system is linear in the
  input, so the unperturbed response is zero and the field is a difference field by construction"). An edge leaving such
  a node always carries a message of zero.
- Diluting. The registered normalisation is item 28's `total_in_degree`: "one per-node divisor over all relations
  together", `values = sign * weight / total_in_degree[target]` where `total_in_degree` is
  `torch.zeros(num_graph_nodes).index_add(0, target, weight)`
  (`mechanistic_pathway_learning/models/linear_response_encoder.py`, `signed_stacked_adjacency`). An edge that always
  carries zero still counts in that divisor, so unreachable content divides the real signal down at every node it feeds.
  This is the mechanism by which superfluous content is harmful and not merely dead weight.
- Duplicated. Two nodes receiving from exactly the same sources under exactly the same relations hold the same field in
  every example, so the second buys parameters and no information.

The adjacency measured is the one the encoder builds, not the stored edge table. The 280,943 stored edges gain a derived
`depletes_substrate` relation (36,562 edges, reaction to metabolite, one per `substrate_of` edge) and a derived
`product_opposes` relation (15,851 edges, product to reaction, one per `product_of` edge of a reversible reaction,
preregistration item 27), giving 333,356; 17,592 of those are then dropped because their source is a currency metabolite
("edges leaving a non-propagating node dropped"), which leaves 315,764 propagated edges. The drop falls on the stored and
the derived blocks alike, so `product_opposes` appears as 11,272 edges in the relation table of section 6 and
`substrate_of` as 23,575: the counts before and after the drop are different quantities and the script now prints both
(`how that count is reached`). The 1,058 `targets` edges are outside every divisor, because a
drug's mechanism edges "are kept out of the normalised adjacency and added on their own": they propagate and dilute
nothing.

The readout side has no node set. Both registered heads read the field at every graph node. The noisy-OR head's support
is a hard-concrete gate over all N nodes and the pool is a sum over them
(`noisy_or_pathway_module_model.py`, `module_activation`: `torch.einsum("kn,bnd->bkd", module_support, node_state_field)`),
and the sigmoid baseline "pools it over all graph nodes (sum of the field, which for the difference field is the total
propagated change)" (`relational_gnn_sigmoid_baseline.py`). So the backward cone from the registered readout is the whole
graph, and the only structural exclusion is the forward cone from the seeds. Two stricter readouts are reported beside it
as counterfactuals, each labelled as one.

Existing measurements reused rather than re-derived:

- `data/processed/label_review/relation_flow.json`, from `experiments/measure_relation_flow.py`, for the per-relation
  flow share and the labelled examples behind it in section 6. It was run on `graph_full_neuronal_split_binders`, so it
  carries no drug layer.
- `docs/membrane_potential_reach.md` for the in-degree dilution of the laboratory readout nodes: "135 of 501 (26.9%)"
  diluted at least tenfold and "36 of 501 (7.2%)" at least fiftyfold, against 4.6% and 1.0% of metabolite nodes in
  general. That measures the size of the divisor, not the dead share of it, so it is complementary to section 5 and not
  the same quantity.
- `docs/missing_physiology_review.md` for the verified facts the task names: 405 `changes_membrane_potential` edges into
  `VM_c` with a measured median field of 2.9e-15, 600 of 266,966 edges running metabolite to protein, 173 of 60,838 edges
  into reactions being regulatory. These were not re-derived.
- `docs/graph_content_null_results.md` and `docs/score_by_layer_contact_neuronal.md` were read and answer none of the five
  items; both report scores rather than graph content, and both are on the 451-gene metabolic slice. The one reading that
  bears on this review is the score side of the same question: the 98 slice genes in contact with the layer
  `graph_neuronal` adds scored -0.008 [-0.056, +0.035] against the base graph, and "Individual additions to the graph's
  content have not moved the score, while removing the wiring wholesale does."

## 3. The perturbation seeds

Under the registered filter (relation `induces`, grades A and B, a positive report) the set is 1,414 gene perturbations
and 222 drug perturbations. A gene knockout seeds its gene node; a drug seeds its own drug node, and its `targets` edges
carry the sign and magnitude on, which is the registered `--drug-entry nodes` mode. That gives 1,636 seed nodes. Under
`--drug-entry targets` the drugs seed 148 target proteins instead and the seed set is 1,562 nodes; the forward cone then
falls from 41,067 to 40,907 nodes, a difference of 160, so the choice of drug-entry mode is not what decides reachability.

## 4. Reachability

| | v4 drug graph | graph_full_neuronal |
| --- | --- | --- |
| nodes | 49,796 | 36,865 |
| stored edges | 280,943 | 266,966 |
| propagated edges | 315,764 | 301,788 |
| nodes a seed can reach | 41,067 (82.47%) | 34,114 (92.54%) |
| nodes no seed can reach | 8,729 (17.53%) | 2,751 (7.46%) |
| edges that can never carry a message | 20,141 (6.38%) | 14,175 (4.70%) |
| of those, in some node's divisor | 20,141 | 14,175 |

Every unreachable node by kind, on the v4 graph:

| node type | unreachable | of | share |
| --- | --- | --- | --- |
| gene | 6,013 | 12,810 | 46.9% |
| protein | 1,949 | 12,709 | 15.3% |
| reaction | 444 | 13,793 | 3.2% |
| metabolite | 263 | 8,580 | 3.1% |
| protein_entity | 60 | 1,681 | 3.6% |
| drug | 0 | 222 | 0% |
| membrane_potential | 0 | 1 | 0% |

Three facts fix what this is.

First, 5,966 of the 6,013 unreachable genes have no incoming propagated edge at all. They are source nodes that no
perturbation seeds and no transcription factor reaches. 69 of the 6,013 gene symbols do appear in the evidence table
under a grade or a relation the registered filter drops, so 5,944 of them are outside the label set entirely.

Second, the gene-protein split is what makes the gene layer unreachable, and the count is a consequence of the split
rather than of the biology. On the v4 split graph the relations entering a gene node are `regulates_transcription_of`
(41,638), `inhibits` (607), `binds` (109) and `activates` (68): essentially the transcriptional one alone. On
`graph_full_neuronal`, where a gene node is also its protein, the same node receives `activates` (69,569),
`regulates_transcription_of` (41,638), `binds` (31,096), `inhibits` (18,090), `cofactor_of` (83) and
`oxidatively_modifies` (41), and the unreachable gene count is 1,980 rather than 6,013. The 1,949 unreachable proteins all
have an incoming edge; they are unreachable because their gene is, and because no reachable protein acts on them.

Third, 6,105 of the 8,729 unreachable nodes have at least one edge into a reachable node, so they dilute rather than sit
idle. The remaining 2,624 are inert: they cost node parameters and nothing else.

The backward cone, with the registered readout and two counterfactuals:

| readout set | readout nodes | backward cone | in neither cone | nodes on a seed-to-readout route |
| --- | --- | --- | --- | --- |
| every node (the registered symptom head) | 49,796 | 49,796 (100%) | 0 nodes, 0 edges | 41,067 (82.5%) |
| every metabolite (the laboratory readout's node kind) | 8,580 | 44,689 (89.7%) | 1,238 nodes, 1,451 edges | 37,198 (74.7%) |
| the 16 curated pathway modules' genes and proteins | 205 | 39,854 (80.0%) | 1,903 nodes, 3,752 edges | 33,028 (66.3%) |

Under the registered head nothing lies in neither cone, because the readout reads every node. The share of the graph that
cannot carry a response from any seed to any readout is therefore the forward-cone figure, 17.53% of nodes and 6.38% of
propagated edges. The two counterfactual rows say what a readout restricted to metabolites or to the curated module genes
would additionally exclude: 1,238 and 1,903 nodes, which is small beside the 8,729.

## 5. Dilution of the total_in_degree divisor

Measured over the 40,334 nodes that both carry a field and have a divisor, with `VM_c` excluded.

| | value |
| --- | --- |
| total divisor mass | 310,969 edges |
| divisor mass that can never carry anything | 16,812 (5.406%) |
| dead share of the divisor, median | 0.0000 |
| dead share, mean | 0.0672 |
| dead share, 90th percentile | 0.3333 |
| dead share, 99th percentile | 0.5602 |
| nodes whose divisor is at least half dead | 2,447 |
| nodes whose divisor is entirely dead | 9 |
| divisor size, median | 3 |
| divisor size, 90th percentile | 14 |
| divisor size, largest | 3,329 (cytosolic H+) |
| gain if the dead edges left the divisor, median | 1.000 |
| that gain, 90th percentile | 1.500 |
| that gain, 99th percentile | 2.167 |
| that gain, largest | 9.000 |

The median node on a route loses nothing, so the headline reading is that the dilution from unreachable content is not
general. It is concentrated, and it concentrates in two relations.

A dead edge only does harm where it lands on a node that holds a field; a dead edge into a node that is itself
unreachable dilutes nothing, because there is no live message at that node to divide. Splitting the destinations that
way is what separates the harmful case from the inert one, and the split is not small:

| | `catalyzed_by` | `encodes` |
| --- | --- | --- |
| dead edges | 6,894 | 5,998 |
| distinct destinations | 3,705 | 5,985 |
| destinations that hold a field, so are diluted | 3,587 reactions | 4,036 proteins |
| destinations that hold no field, so dilute nothing | 118 reactions | 1,949 proteins |
| dead share of the diluted nodes' divisor, median | 0.286 | 0.333 |
| that share, 90th percentile | 0.500 | 0.500 |
| gain if the dead edges left, median | 1.400 | 1.500 |
| that gain, 90th percentile | 2.000 | 2.000 |
| that gain, largest | 8.000 | 9.000 |

The 1,949 proteins in the `encodes` column are exactly the 1,949 unreachable proteins of section 4: their gene is dead,
so the `encodes` edge into them is dead, and the node has no live message to lose. They are the clearest case of
superfluous content that is inert rather than harmful. The 4,036 proteins beside them are the harmful case: they do hold
a field, and a dead `encodes` edge divides it.

The encoder's own reading of the divisor is what makes this a defect rather than a convention: "The in-degree
normalisation reads isozymes as partial redundancy (losing one of k catalysts removes 1/k of the input)". That reading is
correct physiology. At a typical affected reaction, between a quarter and a half of the k is an enzyme this study can
never perturb.

Where the dead divisor mass comes from, by the node type of the edge's source: protein 11,086, gene 6,220, reaction 2,301,
metabolite 478, protein_entity 56. By relation: `catalyzed_by` 6,894, `encodes` 5,998, `activates` 2,291, `product_of`
1,424, `member_of` 873, `depletes_substrate` 870, `regulates_transcription_of` 664, `inhibits` 632, `substrate_of` 298,
`product_opposes` 179, `binds` 10, `changes_redox_pool` 6, `produces_oxidant` 1, `changes_membrane_potential` 1.

Worst by the number of dead edges in the divisor:

| node | name | type | dead | divisor | share |
| --- | --- | --- | --- | --- | --- |
| MAM02039c | H+ | metabolite | 73 | 3,329 | 2.2% |
| MAM02040c | H2O | metabolite | 53 | 2,246 | 2.4% |
| RCTE_749455 | Ligand:GPCR complexes that activate Gi:H... | protein_entity | 33 | 214 | 15.4% |
| RCTE_749445 | Ligand:GPCR complexes that activate Gi:H... | protein_entity | 33 | 214 | 15.4% |
| PROTEIN:RBX1 | RBX1 | protein | 32 | 184 | 17.4% |
| RCTE_380091 | Ligand:GPCR complexes that activate Gi | protein_entity | 28 | 193 | 14.5% |
| MAM01285c | ADP | metabolite | 28 | 780 | 3.6% |
| MAM01371c | ATP | metabolite | 27 | 952 | 2.8% |
| PROTEIN:SKP1 | SKP1 | protein | 27 | 123 | 22.0% |
| MAM02552c | NAD+ | metabolite | 25 | 320 | 7.8% |
| MAM01334c | AMP | metabolite | 25 | 223 | 11.2% |
| MAM01623c | CTP | metabolite | 25 | 69 | 36.2% |
| MAR09490 | Ubiquitin-C-terminal-thiolester hydrolase | reaction | 23 | 57 | 40.4% |
| MAM02751c | Pi | metabolite | 23 | 759 | 3.0% |
| RCTE_3299629 | TAS2Rs | protein_entity | 22 | 26 | 84.6% |

Worst by dead share among nodes with a divisor of at least 20:

| node | name | type | dead | divisor | share |
| --- | --- | --- | --- | --- | --- |
| RCTE_3299629 | TAS2Rs | protein_entity | 22 | 26 | 84.6% |
| RCTE_3299631 | TAS2Rs:Bitter-tasting compounds | protein_entity | 22 | 26 | 84.6% |
| MAR06545 | 1-phosphatidyl-myo-inositol-3,5-bisphosphate... | reaction | 12 | 21 | 57.1% |
| MAR06591 | phosphatidylinositol-3,5-bisphosphate hydrolase | reaction | 12 | 21 | 57.1% |
| MAR01091 | N-acetyl-LTE4 to leukotriene E4 conversion | reaction | 21 | 37 | 56.8% |
| MAR00765 | glucosylceramide pool to LacCer pool conversion | reaction | 20 | 37 | 54.1% |
| MAR01088 | N-acetyl-LTE4 to leukotriene E4 conversion | reaction | 19 | 36 | 52.8% |
| MAR06916 | ATP phosphohydrolase | reaction | 11 | 21 | 52.4% |
| MAM01424c | CDP | metabolite | 14 | 28 | 50.0% |
| MAM02016c | GMP | metabolite | 11 | 22 | 50.0% |

The two TAS2R entities are the clearest single case: 22 of 26 incoming edges are from bitter-taste receptors with no route
from any perturbation, so the complex reads at 4/26 of what its live inputs would give, and bitter taste is not in the
24-symptom set.

The nine nodes whose divisor is entirely dead are worth naming because they show two different causes. Six are seeds that
carry their own input and nothing else: GENE:SLC16A2, GENE:SYNJ1, GENE:CELF2, GENE:MAGI2, GENE:NACC1, GENE:TOR1A. Three
are drug-target proteins reached only through a `targets` edge, which is outside the divisor: PROTEIN:KCNB2,
PROTEIN:KCNG2, PROTEIN:SCN7A.

The 222 currency metabolites hold 27,165 of 314,706 divisor mass (8.6%). Their own outgoing edges are dropped, so they
receive and pass nothing on, and the pooled head still reads the field at them. They are not superfluous under the
registered readout; they are where the largest divisors are.

### VM_c, reported apart

`VM_c` has a divisor of 409: 405 `changes_membrane_potential` edges and 4 `driving_force` edges. One of the 409 is dead
(0.2%). So the pathology recorded in `docs/missing_physiology_review.md` as a "median field of 2.9e-15" is not caused by
superfluous content and is not addressable by paring. It is 408 live edges each divided by 409 at the one node that
represents the plasma membrane, which is the registered `total_in_degree` divisor doing what it is defined to do on a node
with one state variable and 405 contributors. Any fix is a change to the normalisation at that node or a split of the
single membrane-potential node, not a removal of anything.

## 6. Relations

No relation type has zero live edges on the v4 graph. The three with the lowest live share are the ones the orphan gene
layer feeds: `encodes` 6,738 live of 12,736 (52.9%), `catalyzed_by` 17,372 of 24,266 (71.6%), `member_of` 5,894 of 6,767
(87.1%). Every other relation is above 96% live.

One relation is dead in a way the live count does not show. `sequesters` has 101 edges, from 3 carriers (ALB, ORM1, ORM2)
into 62 drug nodes. An entry node that is not the seeded perturbation is held at zero at every step
("the response of an entry node that is not a seed is held at zero at every step"), so each of those 101 edges can carry
something in exactly 1 of the 1,636 perturbations and is inert in the other 1,635. The same holds for the 23 `clears`
edges that `graph_full_neuronal_split_binders_drugs_v4_enzymes` adds.

Relations with a handful of edges, and what the flow measurement already says routes through them. Edges from this
session on the v4 graph; examples and flow share from `data/processed/label_review/relation_flow.json`, measured on
`graph_full_neuronal_split_binders`, where a perturbation counts for a relation when that relation carries at least 0.1%
of its flow:

| relation | edges | live edges | labelled examples | flow share |
| --- | --- | --- | --- | --- |
| changes_membrane_potential | 405 | 404 | 0 | 0.0012% |
| changes_ion_pool | 395 | 395 | 11 | 0.0052% |
| changes_redox_pool | 324 | 318 | 8 | 0.0058% |
| regulates_release | 89 | 89 | 309 | 0.0374% |
| cofactor_of | 83 | 83 | 27 | 0.1573% |
| oxidatively_modifies | 41 | 41 | 15 | 0.0223% |
| produces_oxidant | 28 | 27 | 12 | 0.0317% |
| redox_capacity | 19 | 19 | 20 | 0.0119% |
| voltage_gates | 7 | 7 | 14 | 0.0078% |
| driving_force | 4 | 4 | 1 | 0.0011% |
| catalyzes_reverse_transport | 3 | 3 | 451 | 0.2300% |

`changes_membrane_potential` is the one relation with 405 edges and zero labelled examples above the threshold: every one
of its edges is structurally live and no development perturbation sends a thousandth of its flow along them, which is the
same fact as the divisor of 409 seen from the evidence side. The electrical layer is the clearest case in the graph of
content that is present, correct and unable to carry anything the model can learn from.

## 7. Duplication

Compartment copies. 2,483 of 4,213 metabolite species have more than one compartment copy, covering 6,850 nodes, of which
4,367 are copies beyond the first. The median species with copies has 2, the largest has 13 (MAM02555, NADPH). This is not
superfluous content by the field criterion: exactly 1 field-identical group consists of copies of one species, the pair
MAM00755m and MAM00755c (3alpha,7alpha-dihydroxy-5beta-cholest-24-enoyl-CoA). Compartment copies have distinct
neighbourhoods and hold distinct fields. 404 of the 6,850 have a total propagated degree of 2 or less, and 29 metabolite
nodes touch only Exchange/demand reactions; all 29 are extracellular and none is reachable.

Pool nodes. 14 exist: `NEURON_MAM02519c`, `NEURON_MAM02200c`, `NEURON_MAM02200e`, `NEURON_MAM01442c`, `NEURON_MAM01413c`,
`REDOX_GSH_c`, `REDOX_NADPH_c`, `REDOX_GSH_m`, `REDOX_NADPH_m`, `REDOX_GSH_r`, `REDOX_NADPH_r`, `REDOX_NADPH_n`,
`REDOX_GSH_x`, `REDOX_NADPH_x`. Every one shares its `base_metabolite_id` with plain compartment nodes of the same
species, so the same chemical species is carried twice: `NEURON_MAM02519c` beside seven MAM02519 nodes,
`REDOX_NADPH_c` beside eight MAM02555 nodes. That is the design (a separately regulated neuronal or redox pool), not an
accident, and the two copies do not hold the same field: the divisors run from 1 to 367 and the neighbourhoods differ.
One pool node is superfluous on the reachability criterion: `REDOX_GSH_x` has a divisor of 1 and is unreachable.

Identical neighbourhoods. Grouping every node by the multiset of its incoming edges, which is the criterion for holding
an identical field in every example:

| | groups | nodes in a group | redundant nodes | largest group |
| --- | --- | --- | --- | --- |
| identical incoming neighbourhood | 851 | 3,606 | 2,755 | 195 |
| of those, with any reachable member | 809 | | 2,672 | |
| identical incoming and outgoing neighbourhood | 184 | 465 | 281 | 15 |
| of those, with any reachable member | 145 | | 232 | |

The five largest field-identical groups are all gene nodes whose only input is one transcription factor's
`regulates_transcription_of` edge: 195 genes from PROTEIN:MYC, 151 from PROTEIN:TP53, 86 from PROTEIN:SP1, 71 from
PROTEIN:ESR1 and 67 from PROTEIN:AR. This has a
consequence for the readout that is worth stating because it is not about parameters: both heads pool the field by a sum
over all nodes, so a transcription factor with 195 single-input targets contributes 195 copies of one number to the pool.
2,755 of 49,796 nodes (5.5%) hold a field that is an exact copy of another node's. They are not removable, since each gene
encodes a different protein and their outgoing edges differ, but the pooled readout cannot tell them apart and the sum
weights them by group size.

The 281 nodes with an identical incoming and outgoing neighbourhood are interchangeable in full. The largest groups are 15
mitochondrially named gene nodes (GENE:P03186, GENE:P03200, GENE:P03206 and 12 more), 7 USP17L paralogue genes, and three
drug groups of 13, 12 and 11.

Drugs. 166 of 222 drugs fall into 36 groups sharing an identical target set with identical signs, the largest 23. The four
largest: 23 benzodiazepines and related GABA_A agents (alprazolam, chlordiazepoxide, clobazam, clonazepam, diazepam,
estazolam and 17 more); 22 sodium-channel agents, cocaine among them (phenytoin, carbamazepine, cocaine, ethotoin,
lidocaine, lamotrigine and 16 more); 12 opioids (dihydrocodeine, fentanyl, hydromorphone, meperidine, methadone, morphine
and 6 more); 9 tricyclics with venlafaxine (amitriptyline, dothiepin, imipramine, nortriptyline, protriptyline,
venlafaxine and 3 more). This reproduces the figure `docs/missing_physiology_review.md` records as "164 of 222 drugs in 36
identical-input groups"; the small difference is that this count groups on the signed target set under the registered
drug-entry mode. Three quarters of the drug perturbations are 36 distinct model inputs, so 130 of the 166 are duplicates
of a group representative and the model cannot order the members of a group.

## 8. Human-GEM subsystems

148 subsystems are named on reaction nodes. 0 of them have no reaction a seed can reach, so the answer to item 5 is a
negative result: there is no subsystem the perturbation set cannot touch. The lowest reachable shares are small
subsystems, not whole arms of metabolism:

| subsystem | reactions reachable | reactions | edges touching |
| --- | --- | --- | --- |
| Vitamin B12 metabolism | 2 | 4 | 46 |
| Fatty acid metabolism | 35 | 51 | 376 |
| Vitamin A metabolism | 7 | 10 | 89 |
| Riboflavin metabolism | 5 | 7 | 44 |
| Pool reactions | 40 | 48 | 4,495 |
| Phenylalanine metabolism | 5 | 6 | 71 |
| Pantothenate and CoA metabolism | 21 | 25 | 207 |
| Nucleotide metabolism | 107 | 125 | 1,463 |
| Isolated | 136 | 157 | 1,950 |
| Metabolism of other amino acids | 21 | 23 | 272 |

Phenylalanine metabolism is the one the design validates by hand, and 5 of its 6 reactions are reachable. Under the
strictest counterfactual readout (the curated pathway modules' genes and proteins) exactly 1 subsystem has no reaction on
a route: Alkaloids biosynthesis, 2 reactions, 16 edges, divisor mass 5.

The four subsystems that are modelling devices rather than physiology, which is where I expected to find the superfluous
content and did not:

| subsystem | reactions | reachable | on a route to the curated-module readout | edges touching | divisor mass |
| --- | --- | --- | --- | --- | --- |
| Transport reactions | 4,037 | 3,844 | 3,252 | 39,583 | 17,323 |
| Exchange/demand reactions | 1,658 | 1,542 | 1,001 | 6,566 | 3,250 |
| Isolated | 157 | 136 | 86 | 1,950 | 973 |
| Pool reactions | 48 | 40 | 40 | 4,495 | 1,978 |

Exchange/demand reactions are boundary pseudo-reactions rather than physiology and they are the largest candidate for a
cut: 1,658 reactions, 6,566 edges, 3,250 of the divisor mass, and 657 of them with no route to the curated-module
readout. The 48 Pool reactions are the opposite case: 48 nodes carrying 4,495 edges, so each is a hub of about 94 edges,
and 40 of the 48 are reachable and on a route. The transport layer is reachable at 95%, which bears on
`docs/transporter_connectivity_audit.md`'s figure of "2,266" transport reactions with no catalysing gene: an uncatalysed
transport reaction is still reachable through its substrate.

## 9. The limits of each measurement

- Reachability is structural. It says a node can hold a nonzero field for some perturbation, not that the field is large
  enough to read. The encoder's own figure for the decay is "the change at L-dopa is about 1e-5 of the response at the
  perturbed gene and at dopamine about 1e-8" (`docs/membrane_potential_reach.md`), so a node at the edge of the forward
  cone is reachable and numerically
  silent. A measurement of the settled field, as `measure_relation_flow.py` makes, is the stronger instrument and this
  review does not repeat it per node.
- The unreachable count is contingent on the registered label filter. 69 of the 6,013 unreachable gene symbols appear in
  the evidence table under a grade or a relation the filter drops, so adding grade C or the `relieves` relation would move
  the count. It is also contingent on the gene-protein split: 6,013 on the v4 split graph against 1,980 on
  `graph_full_neuronal`.
- The dilution figures count edges, because the v4 `edges.parquet` carries no `edge_weight` column; decision 4.11 makes
  the divisor a weight sum rather than a count, so on a graph that does carry weights the divisor and the dead share would
  both change. With every weight 1 the two forms are identical, which is the case here.
- "Dead" means the source can never hold a field. An edge whose source holds a field only under one perturbation out of
  1,636 counts as live, which is why the 101 `sequesters` edges had to be reported separately.
- Field-identical grouping ignores cell-class weights. Under `--cell-class-weights` each node's state is multiplied by its
  class weight at every step, so two nodes with the same incoming neighbourhood and different brain-expression weights do
  not in fact hold the same field. The 2,755 figure is therefore an upper bound on duplication under a cell-class run and
  exact under a run without one.
- The backward cones from the metabolite set and from the curated module genes are counterfactuals, not the registered
  readout. The registered heads read every node, so they are reported to show what a narrower readout would exclude, and
  no decision should rest on them.
- Subsystem labels are Human-GEM's and cover reaction nodes only. 13,793 reactions carry one; genes, proteins and
  metabolites do not, so section 8 says nothing about the signalling layer.
- Everything here is the v4 drug graph. `graph_full_neuronal`, the confirmatory graph, has no drug layer and different
  gene-layer reachability, and the two columns of section 4 are the only place this review measures both.

## 10. What I would put to the owner

Separated from the measurements above; these are inferences and none is a decision.

1. Cut the orphan gene layer, or give it an input. 5,944 gene nodes are outside the label set and have no transcriptional
   regulator, and the cut would remove 6,013 gene nodes and 1,949 protein nodes (16.0% of the graph) and 20,141 propagated
   edges (6.4%), raise the live catalysis message at 3,587 reactions by a median factor of 1.400 and the live message at
   4,036 proteins by a median factor of 1.500. Giving them an input
   instead is the same fix from the other side and is what `docs/missing_physiology_review.md` item 2 asks for on the
   metabolite side. Either way the number to watch is the divisor at those 3,587 reactions, not the node count.
2. Decide the `VM_c` divisor on its own. It is not a paring question: 408 of its 409 incoming edges are live. The options
   are the normalisation at that node or more than one membrane-potential node, and the review that raised it already says
   so ("That is a normalisation defect, not a scope limit").
3. Exchange/demand reactions are the one metabolic cut worth measuring. 1,658 reactions and 3,250 of the divisor mass sit
   on boundary pseudo-reactions, and 657 of them have no route to the curated module genes. The cost of the cut is that an
   exchange reaction is how a metabolite leaves the system, so removing them changes what a depletion can do.
4. Report the 36 drug groups in the paper. 166 of 222 drugs are 36 distinct model inputs, so a per-drug result is a
   per-group result for three quarters of the drug set, and the equifinality endpoint reads differently when the inputs
   are 36 rather than 222.
5. Do not pare the metabolism. 0 of 148 subsystems is unreachable, 3.1% of metabolites and 3.2% of reactions are, and 1
   field-identical pair exists across 2,483 multi-copy species. The measurement does not support a cut there, and
   `docs/graph_content_null_results.md` already shows that removing the wiring wholesale costs score.

## 11. Correction of 14 October: the diluting layer is the proteins, not the orphan genes

Sections 1, 5 and 10 attribute the dilution at reactions to the orphan **gene** nodes: "their 6,894 `catalyzed_by`
edges sit in the `total_in_degree` divisor of 3,587 reactions". **Gene nodes almost never emit `catalyzed_by` in this
graph.** Measured here on the v4 drug graph: of 3,157 distinct `catalyzed_by` sources, 2,839 are `protein` nodes, 316
are `protein_entity` complexes and **2 are gene nodes** (GENE:LRTOMT and GENE:USP41). The gene and protein layers
were split, and catalysis attaches to the protein.

The chain is one hop longer than the review said, and naming it correctly changes what fixes it:

    orphan gene (no incoming edge, never seeded, zero at every step)
      -> encodes -> protein (zero at every step, because encodes is its only input)
      -> catalyzed_by -> the reaction's total_in_degree divisor

Measured: **no protein node lacks an incoming edge at all**, because `encodes` from its gene is always there, but
**3,447 proteins have `encodes` as their only incoming relation**, and those emit **11,314 `catalyzed_by` edges into
4,886 reactions**. That is the set whose share of the divisor is unreachable, and it is larger than the 6,894 and
3,587 the review reported, because the review counted from the wrong layer.

**This is why the owner was right to hold the paring** (14 October: "Let's just wait on dropping the orphans though
because you could find more labels from BRENDA"). The mechanism is not labels - BRENDA carries enzyme and ligand
data, not symptom associations - and it is not the genes either. It is that a metabolite-to-protein edge gives an
`encodes`-only protein a second incoming relation and makes it live, and metabolite-to-protein is exactly what
BRENDA supplies and exactly the gap `docs/missing_physiology_review.md` section 2 measures at 600 of 266,966 edges.
Paring the layer first would delete the nodes the download is meant to connect.

