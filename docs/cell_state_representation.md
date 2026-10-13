# The cell state the model carries

![What the model holds as the state of a cell](figures/cell_state.png)

Drawn by `experiments/make_explanatory_figures.py --only cell_state` from
`data/processed/graph_full_neuronal_split_binders_drugs_v4_enzymes` (49,796 nodes, 280,966 edges). Written 13 October
2026 at the owner's request ("please include an explanatory figure for things like representation of the ph and other
pieces of the cell state").

## The distinction the model rests on

A node's state does one of two things, and which one decides almost everything else about it.

- **It travels.** The state at a node is multiplied by the relation's gain and added to its neighbours, every round,
  until the propagation settles. This is what a reaction, a protein, a metabolite of interest and a transcript do.
- **It caps.** The state is read once and multiplies what a node passes, without moving along any edge. This is what a
  cell class's expression does, what the limiting-step rule does with it, and what compartment pH now does.

Everything that would otherwise be a hub enters the second way, and that is not a convenience. The cytosolic proton
touches 3,329 reactions. A node that sent across all of them would make every pair of reactions neighbours and the
graph would stop meaning anything. The same is true of ATP, water and CoA. So the graph tags 222 metabolite nodes as
currency, and a currency node's **outgoing** edges are dropped while its incoming edges are kept: it receives, its
deviation is computed, and it passes nothing on.

That last sentence is the opening the recent work walks through. The deviation is already there. It was only ever
thrown away.

## What is computed and discarded

The right panel of the figure counts the reaction edges at each group of state-carrying species. The quantities, by
total degree over their compartmental copies:

| state | species | edges | what it reports |
|---|---|---|---|
| pH | H+ | 6,181 | the acid-base state of each compartment |
| redox | NAD+, NADH, NADP+, NADPH, FAD, FADH2, GSH, GSSG, FMN, FMNH2 | ~4,600 | the reducing power available, and oxidative stress |
| energy charge | ATP, ADP, AMP | ~2,700 | how much of the adenylate pool is charged |
| acyl carriers | CoA, acetyl-CoA, L-carnitine | ~2,900 | acetyl supply, which is also the acetylcholine precursor |
| acid-base | HCO3-, CO2 | ~1,375 | the bicarbonate buffer, so part of the pH story |
| metabolic ion pools | Na+, K+, Ca2+, Mg2+, Cl- | ~1,310 | the bulk ionic state |
| oxygen | O2 | 877 | hypoxia |
| methylation | SAM, SAH, THF | 368 | the methylation index |
| nitrogen and sulfur | NH3, H2S | 154 | ammonia load |

## How pH enters, and why nothing else needs a constant

`graph/compartment_ph.py` and preregistration item 22a. After the propagation settles, the deviation at the proton
node governing a compartment becomes a multiplier on every node in it,

    w_pH = exp(-(sensitivity * proton deviation) ** 2)

so a compartment pushed away from its resting proton balance passes less. 18,928 of 49,574 nodes take such a weight,
over ten compartments. The form is even in the deviation, so a gain-of-function stays the exact negative of a
loss-of-function; it is applied after the fixed point, so the contraction is untouched; and it lies in (0, 1], so it
cannot raise a norm. It breaks positive homogeneity, which is the registered cost and the first departure from
linearity in the model.

No compartment pH constant enters, because each reaction's optimum is taken as its own compartment's resting state and
the resting value then cancels out of the form. What does not cancel is the **scale**: the model's deviations are in
arbitrary units, the settled field is of order 1e-11, and an uncalibrated sensitivity would leave the weight at
exactly 1 with no gradient. So `calibrate_ph_sensitivity` sets it once from one batch.

**The better form, for the quantities still to be built.** Most of the rest are naturally ratios: energy charge
`(ATP + 0.5 ADP) / (ATP + ADP + AMP)`, NAD+/NADH, NADPH/NADP+, GSH/GSSG, SAM/SAH. A ratio is dimensionless, so it
needs no calibration and no learned unit scale at all — the one ugly part of the pH weight disappears. Every
denominator exists as a node, GSSG included.

## Two pools for one species, three times over

The graph describes some species twice, from two sources, with no edge between the copies.

- **Protons.** Human-GEM's `MAM02039*` against Reactome's `CHEBI15378*`. The Human-GEM copy holds 50 to 100 times the
  reactions in every compartment both describe, and only the Reactome copy propagates. Resolved by item 22a: the
  Human-GEM node governs the nine compartments it describes and Reactome's governs the vacuole, which only it names.
- **Ions.** The membrane potential reads `NEURON_`-prefixed pools — `NEURON_MAM02519c` (Na+, 369 edges), K+ cytosolic
  and extracellular, Cl-, Ca2+ — which are not currency and do propagate. The metabolic pools are different nodes:
  `MAM02519c` (Na+, 601 edges), currency. **Zero edges join them.** So a perturbation that changes sodium through
  metabolism cannot move the membrane potential. This is defensible as a design, since the neuronal pool is the
  electrically relevant sub-pool near the membrane and the 338 curated `changes_ion_pool` edges are what should move
  it, but it is a reach limit and it bears on the open question about the electrical layer's reach.

## Two defects in how a node is tagged currency

`tag_currency_metabolites` tags a metabolite when its name is on a curated list of 63 cofactor names **or** when its
degree reaches the 0.995 quantile of metabolite degrees. Both parts misfire.

1. **A charge suffix defeats the name list.** It holds `nad`, `nadp`, `hco3` and `glutathione`/`gthrd`/`gthox`, while
   the display names are `NAD+`, `NADP+`, `HCO3-` and `GSH`. So no copy of those is tagged by name, and only the
   high-degree copies are caught by the degree rule. The result is arbitrary: cytosolic NAD+ sends nothing while
   peroxisomal NAD+, with 135 edges, does.
2. **Nodes added after the tag is computed are never tagged.** `tag_currency` runs on the base graph, and
   `reactome_import.py` then adds nodes with `is_currency=False` written in. Thirteen species end up with copies that
   disagree: H2O in the vacuole propagates while its other nine copies do not, and the same for O2, Na+ in the
   intermembrane space and mitochondrion, Ca2+ in five compartments, Mg2+ and FMN.

And one that is a modelling error rather than an inconsistency: **glycine's cytosolic copy is tagged currency** by the
degree rule, with 171 edges, while its other six copies propagate. Glycine is an inhibitory neurotransmitter and the
NMDA co-agonist, so it is not a cofactor in any sense the tag intends. Every other transmitter in the graph
propagates: glutamate at 296 edges, glutamine at 258, aspartate at 107, and dopamine, serotonin, acetylcholine,
noradrenaline, histamine, adenosine, taurine and D-serine below those.

### Fixed in the code on 13 October 2026, and reaching the graph only at the rebuild

The owner directed the third of these ("glycine can be removed as currency as long as there's still a mechanism for
the amount of glycine reaching a receptor to vary ... it should also be consistent with other NTs"), and consistency
with the other transmitters is what the fix gives it: a transmitter is now on `NEVER_CURRENCY_NAMES` and no rule can
tag it, so the amount of glycine reaching a receptor varies the way glutamate's does.

- A name is matched with its trailing charge stripped as well as whole, so `NAD+`, `NADP+`, `HCO3-` reach the list,
  and `GSH` and `GSSG` are named outright since they are not charge variants of `glutathione`.
- `retag_currency_by_name` re-applies the name rule once a graph is complete. `build_physiology_graph.tag_currency`
  and `reactome_import.import_reactome_layer` both call it, so a node added after the degree rule ran no longer keeps
  whatever `is_currency` its adder wrote.
- **A fourth thing, which fixing the third exposed.** The degree rule would have tagged `NEURON_MAM02519c`, the
  cytosolic sodium pool the membrane potential reads, because it carries 369 edges and sits above the 0.995 quantile.
  Until now it was protected only by accident, by the very defect that left late-added nodes untagged. Tagging it
  would have stopped it reaching `VM_c` and silenced the electrical layer. Identifiers a build creates in order to
  propagate are now on `NEVER_CURRENCY_IDENTIFIER_PREFIXES` and no rule may tag them.

Run over the node table as it stands, the fix takes the currency set from 222 to 290: 69 copies of H2O2, H+, NH4+,
O2-, NADP+, Ca2+, GSSG, NAD+, GSH, Fe2+, Cl-, Na+, Fe3+, FMN, H2O, Mg2+, HCO3- and O2 gained, glycine released, and
all five neuronal pools left propagating. None of this reaches the graph until it is rebuilt.

## What a reader should take from the figure

The model's picture of a cell is a field over reactions, proteins and metabolites, shaped by three ceilings: which
cell class expresses a node, how limited the least-limited route to it is, and now the compartment's acid-base state.
The ceilings are where physiology enters without the graph collapsing into a single hub, and the bar panel is the
measure of how much physiology is still sitting there computed and unused.
