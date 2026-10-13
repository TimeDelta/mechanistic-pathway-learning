"""Linear-response encoder: a time-invariant, signed state space over the graph (design section 5.2, route 1).

The message passing encoder (relational_message_passing_encoder.py) emulates time badly for a chronic perturbation:
each of its L layers has its own weights, so the "dynamics" change from step to step; the perturbation is only the
initial state, so its signal is a transient rather than a sustained input; the difference field is exactly zero more
than L edges from the perturbed nodes (two on the default configuration, which on the metabolic layer reaches an
enzyme's reactions and their direct products and stops); and edges only run substrate -> reaction -> product, so a
blocked reaction cannot raise its substrate.

This encoder treats the field as the steady-state response of a linear system to a sustained input, the reading that
metabolic control analysis gives a small enzyme-level change (Kacser and Burns 1973; Heinrich and Rapoport 1974):

    h(t + 1) = (1 - damping) * h(t) + damping * (sum_r gain_r * (S_r h(t)) + u)

where u carries the perturbation at its nodes at every step, S_r is the signed, in-degree-normalised adjacency of
relation r (entry [destination, source] = sign of the edge / in-degree of the destination under r) and gain_r is a
learned per-channel gain. The same transition applies at every step, so after num_propagation_steps the field
approaches the fixed point h = sum_r gain_r * (S_r h) + u (damping sets the transient only). Gains of relations whose
edges carry a biological sign are positive, so the sign comes from the graph; relations without a sign (binds) get a
gain of either sign. Each node averages within each relation feeding it and then across those relations (entry =
sign / (in-degree under r * number of relations feeding the node)), so every row of sum_r |S_r| sums to at most one;
gains are bounded by contraction < 1, so the transition is a contraction in the max norm and the iteration converges.
(Bounding gains by contraction / R instead shrank the signal by about 0.2 per step, so changes six edges away were
of order 1e-9.)

Currency metabolites (water, protons, ATP, NAD(P)H and the like, flagged in the graph) receive the response but do
not pass it on: through them every reaction would reach every other in two steps, which is the hub problem of design
assumption A9 in another form. Carrier edges (cofactor_edges, mechanistic_pathway_learning/graph/cofactor_edges.py:
biopterins, quinones, folate carriers, transamination pairs and the like) can be given relations of their own,
cosubstrate_of (carrier -> reaction, supply) and coproduct_of (reaction -> carrier), whose gains are learned apart from
those of main substrates and products, and carriers get no depletion edge. A carrier pool is recycled (tetrahydrobiopterin
used by PAH is regenerated through PCBD1 and QDPR), so one consumer slowing down barely moves the pool, and the depletion
edge made a PAH loss raise tetrahydrobiopterin and through it tyrosine hydroxylase flux, L-dopa and dopamine, the
opposite of phenylketonuria. A synthesis defect does drain the pool, so the supply edges stay: dropping carrier edges
altogether left a GCH1 loss lowering tetrahydrobiopterin with nothing downstream able to see it, and 14 slice genes with
30 positive pairs make or recycle carriers (GCH1, PTS, SPR, QDPR, PCBD1, MTHFR, MTRR, SLC46A1, COQ2, COQ5, COQ7,
ALDH7A1, SLC19A3, HLCS). Untrained on the slice, with only the depletion edge dropped, a PAH loss raises phenylalanine
and lowers tyrosine, L-dopa and dopamine; a GCH1 loss lowers tetrahydrobiopterin, L-dopa and serotonin (dopamine still
rises, through a route not yet traced); an MTHFR loss raises homocysteine and lowers methionine. Carrier gains start low
(logit offset -3) and are learned.

The stoichiometric reading needs edges the graph does not store: a reaction's flux change depletes its substrates,
so every substrate_of edge (metabolite -> reaction) gets a reverse edge reaction -> metabolite with sign -1 under its
own relation, depletes_substrate. With catalyzed_by (gene -> reaction, +1), product_of (reaction -> metabolite, +1)
and substrate_of (metabolite -> reaction, +1, mass action), a loss of function at an enzyme lowers its reactions'
flux and its products and raises its substrates, and both changes travel on to the reactions that use them. The
in-degree normalisation reads isozymes as partial redundancy (losing one of k catalysts removes 1/k of the input) and
buffers metabolites with many producers. Transcriptional edges keep their own sign (CollecTRI activation +1,
repression -1), which the message passing encoder's single unsigned relation dropped.

The response propagates on a few channels (propagation_channels, default 4), each with its own per-relation gains,
so the channels settle at different time scales; a linear map then expands them to the node_state_dim channels the
heads read. Propagating all node_state_dim channels was memory-bound: each step wrote and reread a
[relations, nodes, batch, channels] tensor (about 200 MB on the metabolic slice at batch 16 and 32 channels), 19
seconds per batch forward and backward on one core.

The response decays steeply with distance under in-degree averaging (on the slice, the change at L-dopa is about 1e-5
of the change at the perturbed gene and at dopamine about 1e-8), so a linear readout sees little beyond the first few
edges. response_scale "signed_log" maps each channel through sign(h) * log(1 + |h| / s_c) with a learned scale s_c
(starting at 1e-4), which keeps the sign and the order of magnitude and makes distant changes readable; it is odd, so
a gain of function still mirrors a loss of function. response_scale "linear" leaves the response as it is.

Cell classes enter as channels (cell_class_weights, [N, K] in [0, 1], from
mechanistic_pathway_learning/graph/brain_expression_weights.py): with them the response propagates on K groups of
channels_per_cell_class channels over the same adjacency, and in class k each node's new state, input included, is
multiplied by its weight w_k at every step, h_k(t + 1) = (1 - damping) * h_k(t) + damping * w_k * (messages + u). A
node a class does not express then neither holds nor passes that class's response, and the per-class diagonal scaling
by weights at most one keeps the contraction. extracellular_pool_nodes makes chosen nodes (extracellular metabolites)
one pool shared by the classes in shared_pool_classes: after each step their state is the mean over those classes, so
a transmitter released in one class reaches every other. With an all-ones class among the K, the original propagation
is one of the channel groups. The readout expands the K * channels_per_cell_class channels to node_state_dim as before.

The system is linear in the input, so the unperturbed response is zero and the field is a difference field by
construction: no base state, and no node identity, enters it. Node kinds enter only through an output gate,
field = (h W_expand) * sigmoid(node_features W + b), which lets a pooled readout tell where the response landed
without adding anything where it did not. With node descriptors among the features, descriptor_treatment
(mechanistic_pathway_learning/models/descriptor_treatments.py) decides how they enter the gate: plain (one map of all the
columns), seed_masked (at each perturbation's perturbed nodes, where the response is largest, the gate reads the
structural columns only) or zero_init_slow (a separate descriptor map in the gate's logit that starts at zero and
trains at its own learning rate).

Entry nodes (entry_node_mask; a drug node, mechanistic_pathway_learning/graph/drug_entry_nodes.py) exist only in the
perturbation seeded on them: the response of an entry node that is not a seed is held at zero at every step. The edges
of entry_relation_names (a drug's mechanism edges) are kept out of the normalised adjacency and added on their own,

    message into target = sign * weight * gain_entry * h_source

with one positive gain per channel for each entry relation. Left in the adjacency they would be divided by the
target's in-degree under the relation, which is the number of study drugs sharing the target, and they would add one
to the number of relations feeding the target, which divides every other relation's message there for every
perturbation, gene knockouts included (docs/drug_entry_nodes_measured.md). An entry node with no incoming edge holds
its own input at every step, so its mechanism edges deliver a sustained input sign * weight * gain_entry * u to its
targets: the same input a seed on the target delivers, scaled per channel by the gain and one step later.
"""
from __future__ import annotations

import math

import torch
from torch import Tensor, nn

from mechanistic_pathway_learning.models.descriptor_treatments import (
    check_descriptor_treatment,
    features_without_descriptors,
    seed_unit_mask,
)

SUBSTRATE_RELATION = "substrate_of"
PRODUCT_RELATION = "product_of"
DEPLETES_SUBSTRATE_RELATION = "depletes_substrate"
COSUBSTRATE_RELATION = "cosubstrate_of"
COPRODUCT_RELATION = "coproduct_of"
UNSIGNED_RELATIONS = ("binds",)
# Where the edge signs come from. "graph" uses the stored signs and the -1 of the derived depletes_substrate relation;
# "all_positive" sets every edge to +1 and keeps every gain positive, the sign ablation of the same pipeline on the
# same graph and split (the check docs/literature_appraisal_staging.md applies to every appraised source). With no
# negative entries nothing cancels at a node, so the ablation removes the sign information and the cancellation
# together; the total_in_degree normalisation is the arm that keeps the signs and removes most of the cancellation.
# "permuted" shuffles the signs among all edges, the derived ones included, with a fixed seed, so the number of
# negative edges and the cancellation they can cause stay while their placement by mechanism is lost; set beside
# "all_positive", it separates the sign information from the cancellation (suggested in the literature session's
# closing notes). The seed is fixed so a resumed run rebuilds the same permutation.
EDGE_SIGNS = ("graph", "all_positive", "permuted")
SIGN_PERMUTATION_SEED = 0
# Whether each relation learns its own gain. "per_relation" is the default; "shared" ties every relation to one gain
# per channel (and drops the carriers' weak starting offset), the relation-typing ablation of the same pipeline. The
# adjacency is still normalised per relation under "in_degree", so "shared" removes the learned part of the typing
# and keeps the structural part; "total_in_degree" with "shared" removes both.
RELATION_GAINS = ("per_relation", "shared")

# How the messages of the relations feeding a node are combined into one. "mean" is the original
# behaviour: each relation's normalised message is multiplied by its gain and the results are summed,
# with the per-node divisor relations_feeding_node making that sum a mean. It annihilates a node that
# receives both signs, which is 36% of destinations (docs/membrane_potential_reach.md): a metabolite
# both produced and consumed gets product_of at mean +1 and depletes_substrate at mean -1, and their
# mean is zero. "softmax_mixture" drops that divisor and combines the feeding relations through a
# softmax over the statistics below, with separate weights per node type, so a destination that the
# mean blanks still carries the largest, smallest or most dispersed of its messages. The mixture
# changes which statistic is read, not the direction at an individual node, since its weights are per
# node type and so shared by every node of that type; direction is the normalisation's business.
CROSS_RELATION_AGGREGATORS = ("mean", "softmax_mixture")
# Every statistic offered here is odd, f(-x) = -f(x), because the response must mirror: a gain of
# function has to produce the negative of the field a loss of function produces (section 5.2). A plain
# maximum is not odd and is worse than even, it is sign-biased: over negative messages it returns the
# least negative, which at the first step is the zero of an untouched neighbour, so it suppresses the
# whole falling branch of the response and a metabolite that should drop stays at exactly zero. The
# minimum fails in mirror image, and a standard deviation is even, so it reports the same dispersion
# whichever way the field moves.
#
# The list is short because the graph is. No destination on data/processed/graph is fed by three or
# more relations, 73.1% are fed by exactly two and the rest by one; on graph_neuronal three relations
# reach 1.5% of destinations. At two values the median and the midrange both equal the mean, and a
# third moment does not exist, so a skewness term is unavailable on this graph whatever one thinks of
# it. What differs at two values is the mean, which message dominates, and the dispersion. The
# dispersion is the interesting one: at two messages of opposite sign, signed_standard_deviation is
# sign(m1 + m2) * (|m1| + |m2|) / 2, so it keeps the whole magnitude where the mean cancels and takes
# its direction from whichever of production and consumption wins. median is kept because the full
# graph may yet feed some destination three relations, and it is documented as collapsing to the mean
# below that.
#
# power_mean covers the span between the mean and the dominant message continuously, as
# sum_r m_r |m_r|^(p-1) / sum_r |m_r|^(p-1) with a learned exponent per node type: p = 1 is the mean
# and p -> infinity is signed_maximum_magnitude. It is odd, bounded by the largest |m_r|, and unlike a
# discrete mixture it has a gradient everywhere along that span.
MIXTURE_STATISTICS = ("mean", "signed_maximum_magnitude", "median", "signed_standard_deviation", "power_mean")
# How the mixture weights are formed from the logits. "softmax" spreads weight over every statistic, so
# a statistic a node type does not need keeps a small share for ever. "sparsemax" (Martins and
# Astudillo 2016) is the Euclidean projection onto the simplex and returns exact zeros, which makes the
# learned choice reportable: a node type either uses a statistic or it does not. Both are non-negative
# and sum to one, so the convex-combination bound that keeps the iteration contracting holds for either.
# The cost of exact zeros is a zero gradient outside the support, so a statistic dropped early can
# never return; entmax with the exponent between one and two interpolates, and MIXTURE_EXPONENT_SPARSEMAX
# names the endpoint rather than hard-coding it.
MIXTURE_WEIGHTINGS = ("softmax", "sparsemax")
# Running to a fixed point (amendment item 7 of 12 October 2026). The tolerance is on the relative change of the state
# between rounds; 1e-3 is the registered value. The cap exists so the loop terminates even where the update is not a
# contraction in practice, and 60 is chosen from the measurement rather than for roundness: at the encoder's contraction
# bound of 0.9 the structural propagation reaches 1e-3 at 23 rounds, and at gains of 0.7 and 0.5 at 14 and 11
# (docs/label_source_review.md 17.2), so 60 leaves room for a slower graph while still ending. A run that reaches the
# cap has not settled, and propagation_steps_taken says so in the run record.
# How a cell class's node weights attenuate a response along a path (the owner, 13 October 2026: "min along the path
# makes WAY more sense to me. do that so it aligns with the physiology").
#   "product", what every run before 13 October 2026 did: the weight multiplies the messages arriving at a node on
#     every round, so a path carries the product of the weights of the nodes it passes through. That is exactly what a
#     per-class copy of the graph with edge weights set by the destination's expression would give, which is how the
#     rule was justified, but it is the arithmetic of independent per-step failure, and it compounds. How much it
#     compounds is to be measured on the rebuilt v4 graph by experiments/measure_limiting_step_rule.py; no figure is
#     stated here because the only measurement so far was taken on a retired graph.
#   "limiting_step": the throughput of a chain at steady state is set by its weakest step, not by the product of its
#     steps - five reactions each at 2 percent of capacity pass about 2 percent of the flux, not 0.02^5. So the
#     attenuation at a node is the bottleneck of its least-limited route from the perturbed nodes,
#         b(t) = max over routes to t of min over nodes n on the route of w(n),
#     computed by the max-min recurrence of class_bottleneck below. The same limiting-step logic is already used
#     inside a reaction's gene rule, where evaluate_gene_rule takes the minimum over "and".
#     Over which routes. The propagation runs to a fixed point because the model's claim is the steady state, and at
#     a steady state every route contributes, so b is taken over every route and the recurrence runs until it stops
#     changing. Bounding it by the round at which the propagation stopped would make the stopping rule, which is a
#     numerical tolerance and is decided for a whole batch at once, set a physiological ceiling, and would make one
#     perturbation's ceiling depend on the perturbations it happened to share a batch with. Taken over every route, b
#     depends on nothing but the graph, the class weights and the perturbation's own seeds, so it is computed once per
#     perturbation and kept (cached_class_bottleneck). A run that reads the state at a chosen round (num_steps, or a
#     fixed step count with no tolerance) takes b over that many rounds, so a transient reading stays coherent.
#     Two places the recurrence follows the propagation rather than the classes' weights alone. The extracellular pool
#     (--extracellular-coupling) is one compartment shared by the pooled classes, so a metabolite one class releases
#     reaches a class that cannot make it; b at a pool node is therefore the largest over the pooled classes, which is
#     what a route that enters the pool in one class and leaves it in another gives. Without that the rule would set a
#     ceiling of zero at the pool in every class that cannot make the metabolite and so undo the coupling. And an entry
#     node (a drug) that a perturbation does not seed is absent, as in present_nodes, so no route passes through it.
# What the limiting-step rule is not: a per-path minimum. The state at a node sums contributions from routes of
# different lengths and no per-node scalar can hold a different minimum for each, so b is the bottleneck of the best
# route and scales everything the node holds.
# Where the two rules apply their weights, which is forced and not a choice. The product rule multiplies the messages
# arriving at a node on every round, so the attenuation is inside the propagation. A minimum cannot be, because a
# factor applied on every round compounds by construction, whichever factor it is: multiplying the messages by b every
# round would carry the product of b along the path, which is the failure the rule exists to remove. So the bottleneck
# is factored out of the loop and multiplies the settled state once, which is the only local way to make the
# attenuation at a node the bottleneck of its route rather than a product along it.
# What that costs, stated plainly because it is not a drop-in swap. Under the product rule a class's weights are part
# of its propagation, so class k answers a perturbation over its own weighted copy of the graph and the classes differ
# in which routes carry the response. Under the limiting-step rule the weights leave the loop, so inside the loop the
# classes differ only in their learned per-channel gains and the class's expression enters afterwards as a per-node
# ceiling on throughput. The rule therefore conditions the magnitude of a class's response on its expression and no
# longer conditions the route. It also moves the extracellular pool: share_extracellular_pool runs inside the loop, so
# it now mixes states that no class has attenuated yet, which is closer to a pool that is physically shared than
# mixing per-class attenuated states was.
# Two further overstatements, in opposite directions, which is why this is a registered choice and not a correction.
# The product rule understates: a route of expressed nodes is charged a factor at every step. The limiting-step rule
# overstates: the settled state at a node already sums contributions from every route, routes through barely expressed
# nodes among them, and b scales that whole sum by the best route's bottleneck, so a contribution that arrived through
# a thin route is credited at the widest route's ceiling. Masking the thin routes instead would need a cut-off on the
# weights, and the weights have no natural one (a smooth decline, median 0.33, no gap: docs/label_source_review.md).
# Three properties the rule does keep, which is what makes it safe to put in place of the product: zero input still
# gives exactly zero output, so the unperturbed field stays zero and the difference field is still the field; b does
# not depend on the input's magnitude, only on the graph, the class weights and which nodes are perturbed, so the
# response stays positively homogeneous of degree one and the input-size controls still mean what they meant; and b
# lies in [0, 1], so multiplying by it cannot raise a norm. The fixed point is reached by the propagation before b is
# applied, so the contraction the stopping rule rests on is the unweighted one and is untouched.
CLASS_WEIGHT_RULES: tuple[str, ...] = ("product", "limiting_step")
# How the perturbation is scaled into each channel (the owner, 13 October 2026: "each class should be consistently
# signed within itself"). "signed", every run before that date: one free weight per channel drawn from a normal
# distribution, so each channel started with a random sign and size. The readout absorbs a flipped sign everywhere
# except the extracellular pool, whose mean across classes is fixed and averaged copies of one extracellular quantity
# carrying opposite signs: on a toy graph the pool lost a factor of 1.2 to 3.1 to cancellation and read negative for a
# positive perturbation at one seed. "positive": the weight is exp of a free parameter, so every channel of every class
# carries a perturbation with its own sign, and all channels start at the same value, so no class starts out weighing
# more in the pool than another. The classes stay free to learn different sizes.
INPUT_WEIGHT_FORMS: tuple[str, ...] = ("positive", "signed")
# Compartment pH as a throughput ceiling (graph/compartment_ph.py has the reasoning and the measurements; the owner,
# 13 October 2026, of this design: "perfect solution. implement it"). After the propagation settles, the deviation at
# the proton node governing each compartment becomes a multiplier on every node in it,
#     w_pH = exp(-(sensitivity * proton deviation) ** 2),
# so an acidified or alkalinised compartment passes less. A weight does not propagate, which is why this does not
# reintroduce the hub the currency mask exists to stop: the cytosolic proton touches 3,329 reactions, and a node that
# sent freely across them would make every pair of reactions neighbours. The sensitivity is learned per compartment
# because the model's deviations are in arbitrary units, and PH_SENSITIVITY_TARGET_ARGUMENT is what
# calibrate_ph_sensitivity sets the argument's root mean square to on the batch it is given, since an uncalibrated
# sensitivity against deviations of order 1e-11 would leave the weight at exactly 1 with no gradient.
# The weight is a bell curve in pH units, not in [H+] (the owner, 13 October 2026: "why not use true pH on the log
# scale so it fits nicely into equations?"). The deviation is read as a fractional change in [H+], so
#     pH = resting pH - log10(1 + sensitivity * deviation),   w_pH = exp(-(delta pH / width) ** 2)
# with the optimum taken as the compartment's resting state, which makes w_pH a function of delta pH alone. The log is
# why this needs the resting values at all and why the ratios of the other cell-state quantities do not: a ratio is
# dimensionless but says nothing about how many protons there are, and a reaction consuming H+ is limited by the count.
# The fraction is clamped above -1 so the logarithm stays defined; a deviation that would empty the compartment of
# protons is held just above empty rather than being allowed to produce a NaN.
PH_SENSITIVITY_TARGET_ARGUMENT = 0.5  # root mean square of (sensitivity * deviation) after calibration
PH_WEIGHT_EPSILON = 1e-30
DEFAULT_PH_ACTIVITY_WIDTH = 1.0  # pH units; the half-width of a typical enzyme pH-activity curve, and a learned floor
SMALLEST_PROTON_FRACTION = -0.999999
BOTTLENECK_CHUNK = 32  # columns of the flattened batch-by-class gather taken at once; see class_bottleneck
DEFAULT_FIXED_POINT_TOLERANCE = 1e-3
DEFAULT_MAX_PROPAGATION_STEPS = 60
FIXED_POINT_NORM_EPSILON = 1e-12
# how the signed adjacency is scaled so the iteration contracts: "in_degree" divides each entry by the
# destination's in-degree under its relation times the number of relations feeding it, which averages and so
# decays the response steeply with path length; "spectral" divides every entry by one global constant, the
# spectral radius of the unsigned aggregate, which contracts just as surely without the per-node divisor.
# The per-node divisor is what discards the stoichiometric counts that decide direction at a metabolite: five
# producing reactions against one consuming reaction become one mean against another, and 92% of metabolites
# end up with an aggregate row sum of zero (docs/membrane_potential_reach.md). The spectral arm keeps those
# counts. Against it, that constant is 36.3 on the neuronal graph while the in-degree divisor is below it at
# 98.8 percent of destinations, so it divides harder nearly everywhere and the field reaching the membrane
# potential fell to 0.7 of its in-degree value; magnitude is not the same question as direction, so the arm is
# open rather than refuted. That also makes the spectral arm two changes at once, counts kept and the field shrunk
# nearly everywhere, so a loss on symptoms could not say which one cost it. "total_in_degree" separates them: it
# divides each entry by the destination's in-degree summed over every relation, so the unsigned row sums are one
# as under "in_degree" (the same contraction, no global shrink) while the counts survive: five producing reactions
# against one consuming reaction give (5 - 1) / 6 rather than 0, and the sign at a metabolite is set by its own
# stoichiometry rather than by the one gain ratio every metabolite shares.
# A third arm, "in_degree_power", raised the in-degree to an exponent before dividing and then rescaled the
# whole matrix once. It is removed as a measured negative result: the global rescale cancels the softening
# wherever in-degrees are uniform, and on the real graph the median field at VM_c fell to 2.9e-15 at an
# exponent of 0.5 and 1.0e-18 at 0.25, against 3.2e-10 for plain in-degree.
NORMALISATIONS = ("in_degree", "spectral", "total_in_degree")
CARRIER_INITIAL_GAIN_LOGIT_OFFSET = -3.0
INITIAL_LOG_RESPONSE_SCALE = -9.2  # log(1e-4)


def project_onto_simplex(logit: Tensor, dim: int) -> Tensor:
    """Sparsemax: the Euclidean projection of each slice onto the probability simplex.

    Martins and Astudillo 2016 (arXiv:1602.02068). Sort descending, take the largest support size whose
    sorted value still exceeds the running threshold, then shift by that threshold and clamp at zero.
    The result is non-negative and sums to one, like a softmax, but with exact zeros outside the
    support, so the statistics a node type does not use can be read off rather than inferred.
    """
    sorted_logit, _ = logit.sort(dim=dim, descending=True)
    cumulative = sorted_logit.cumsum(dim=dim)
    support_size = torch.arange(1, logit.shape[dim] + 1, device=logit.device, dtype=logit.dtype)
    support_size = support_size.view([-1 if axis == dim else 1 for axis in range(logit.dim())])
    is_in_support = sorted_logit * support_size > cumulative - 1.0
    chosen_size = is_in_support.to(logit.dtype).sum(dim=dim, keepdim=True)
    threshold = (cumulative.gather(dim, chosen_size.long() - 1) - 1.0) / chosen_size
    return (logit - threshold).clamp_min(0.0)


class CrossRelationAggregator(nn.Module):
    """Combines the gain-weighted messages of the relations feeding each node through a softmax over
    simple statistics, with separate weights per node type and per propagation channel.

    The mean annihilates a node fed by relations of opposite sign, so the statistics that survive sign
    mixing are offered beside it and the weights are learned. Only the relations that actually feed a
    node take part, since a relation with no edge into the node would otherwise contribute a zero that
    the median and the signed maximum magnitude would both read as a real message. Every statistic is
    odd, for the reason given at MIXTURE_STATISTICS, so the combination is odd and the response still
    mirrors under a change of perturbation sign.

    The contraction survives. With the per-node divisor dropped, each relation's rows sum to at most one
    in absolute value, so every per-relation message is bounded by maximum_gain * max |h|. The maximum,
    the minimum and the mean over a set of such values are bounded by the same quantity, and so is the
    population standard deviation, since the dispersion of values inside an interval cannot exceed its
    half-width times two. Softmax weights sum to one, so the combination is convex and keeps the bound.
    Zero input still gives exactly zero output, so the unperturbed field stays zero and the response
    stays positively homogeneous of degree one; it loses additivity, which linearity had supplied.
    """

    def __init__(self, num_node_types: int, propagation_channels: int,
                 statistics: tuple[str, ...] = MIXTURE_STATISTICS, weighting: str = "sparsemax",
                 initial_mean_weight: float = 0.9) -> None:
        super().__init__()
        unknown = [statistic for statistic in statistics if statistic not in MIXTURE_STATISTICS]
        if unknown:
            raise ValueError(f"unknown mixture statistics {unknown}, expected from {MIXTURE_STATISTICS}")
        if not statistics:
            raise ValueError("at least one mixture statistic is required")
        if weighting not in MIXTURE_WEIGHTINGS:
            raise ValueError(f"weighting must be one of {MIXTURE_WEIGHTINGS}, not {weighting!r}")
        self.statistics = tuple(statistics)
        self.weighting = weighting
        if "power_mean" in self.statistics:
            self.log_power_mean_exponent = nn.Parameter(torch.zeros(num_node_types, propagation_channels))
        # starts near the mean, so an untrained mixture reproduces the original behaviour up to the divisor
        initial_logit = torch.full((num_node_types, len(self.statistics), propagation_channels), 0.0)
        if "mean" in self.statistics:
            initial_logit[:, self.statistics.index("mean"), :] = self.warm_start_gap(initial_mean_weight, len(self.statistics), weighting)
        self.mixture_logit = nn.Parameter(initial_logit)

    @staticmethod
    def warm_start_gap(initial_mean_weight: float, num_statistics: int, weighting: str) -> float:
        """Logit gap that puts initial_mean_weight on the mean and spreads the rest evenly, under either weighting.

        The gap is solved for rather than chosen so the two weightings start from the same mixture, which a
        comparison between them needs. It also keeps sparsemax off its own boundary: with logits [a, b, ..., b],
        sparsemax keeps every statistic in its support only while a - b < 1, independent of the number of
        statistics, and outside the support a logit has exactly zero gradient, so a gap of 1 or more freezes
        the mixture at the mean for the whole of training. Requiring initial_mean_weight < 1 is the same bound:
        the sparsemax gap (num_statistics * w - 1) / (num_statistics - 1) is below 1 exactly when w is.
        """
        if num_statistics < 2:
            return 0.0
        if not 1.0 / num_statistics < initial_mean_weight < 1.0:
            raise ValueError(f"initial_mean_weight must be in (1/{num_statistics}, 1) so the mean is favoured and "
                             f"every statistic keeps a gradient, not {initial_mean_weight}")
        if weighting == "sparsemax":
            return (num_statistics * initial_mean_weight - 1.0) / (num_statistics - 1)
        return math.log(initial_mean_weight * (num_statistics - 1) / (1.0 - initial_mean_weight))

    def forward(self, relation_messages: Tensor, relation_feeds_node: Tensor, node_type_index: Tensor) -> Tensor:
        """relation_messages [R, N, B, C] already multiplied by their gains, relation_feeds_node [R, N]
        true where the relation has an edge into the node, node_type_index [N]; returns [N, B, C]."""
        feeding = relation_feeds_node[:, :, None, None].to(relation_messages.dtype)
        # a node with no incoming edge of any relation has no statistics to take, so its output is zero;
        # without this guard the median gathers the padding value and returns an infinity
        node_has_input = (feeding.sum(dim=0) > 0).to(relation_messages.dtype)
        num_feeding = feeding.sum(dim=0).clamp_min(1.0)
        masked = relation_messages * feeding

        relation_mean = masked.sum(dim=0) / num_feeding
        computed: list[Tensor] = []
        for statistic in self.statistics:
            if statistic == "mean":
                computed.append(relation_mean)
            elif statistic == "signed_maximum_magnitude":
                # the message furthest from zero, keeping its sign; odd, and it survives a cancellation
                largest = masked.abs().argmax(dim=0, keepdim=True)
                computed.append(masked.gather(0, largest).squeeze(0))
            elif statistic == "median":
                # the median of the feeding relations only: non-feeding entries are pushed to one end by
                # replacing them with the largest feeding magnitude plus one, then the median is read off
                # the feeding count rather than off the full relation axis
                padding = masked.abs().amax(dim=0, keepdim=True) + 1.0  # sorts above every real message
                ordered, _ = torch.where(feeding == 0, padding.expand_as(masked), masked).sort(dim=0)
                # num_feeding is [N, 1, 1], and gather does not broadcast its index, so the index is
                # expanded over the batch and channel axes; without it the median reads batch row 0 of
                # channel 0 and returns it for every row (and the stack below fails outright)
                index_shape = (1, *masked.shape[1:])  # [1, N, B, C]
                half_way = (num_feeding - 1) / 2
                lower_index = half_way.floor().long()[None].expand(index_shape)
                upper_index = half_way.ceil().long()[None].expand(index_shape)
                computed.append(0.5 * (ordered.gather(0, lower_index) + ordered.gather(0, upper_index)).squeeze(0))
            elif statistic == "power_mean":
                # sum_r m_r |m_r|^(p-1) / sum_r |m_r|^(p-1): the mean at p = 1, the dominant message as p grows
                exponent = torch.exp(self.log_power_mean_exponent)[node_type_index][None, :, None, :]  # [1, N, 1, C]
                magnitude = masked.abs().clamp_min(1e-30) ** (exponent - 1.0) * feeding
                computed.append((masked * magnitude).sum(dim=0) / magnitude.sum(dim=0).clamp_min(1e-30))
            else:  # signed_standard_deviation: population dispersion over the feeding relations, signed by the mean
                squared_deviation = ((relation_messages - relation_mean[None]) ** 2) * feeding
                dispersion = (squared_deviation.sum(dim=0) / num_feeding).clamp_min(0.0).sqrt()
                computed.append(torch.sign(relation_mean) * dispersion)

        stacked = torch.stack(computed, dim=0)  # [S, N, B, C]
        simplex_weight = (torch.softmax(self.mixture_logit, dim=1) if self.weighting == "softmax"
                          else project_onto_simplex(self.mixture_logit, dim=1))
        weight = simplex_weight[node_type_index]  # [N, S, C]
        return (stacked * weight.permute(1, 0, 2)[:, :, None, :]).sum(dim=0) * node_has_input


class LinearResponseEncoder(nn.Module):
    def __init__(
        self,
        num_graph_nodes: int,
        relation_types: list[str],
        edge_source: Tensor,
        edge_target: Tensor,
        edge_relation: Tensor,
        edge_sign: Tensor,
        node_features: Tensor,
        node_state_dim: int,
        non_propagating_nodes: Tensor | None = None,
        cofactor_edges: Tensor | None = None,
        num_propagation_steps: int = 8,
        fixed_point_tolerance: float | None = DEFAULT_FIXED_POINT_TOLERANCE,
        class_weight_rule: str = "limiting_step",
        input_weights: str = "positive",
        cache_class_bottlenecks: bool = True,
        ph_proton_node_index: Tensor | None = None,
        ph_compartment_of_node: Tensor | None = None,
        ph_compartment_names: list[str] | None = None,
        max_propagation_steps: int = DEFAULT_MAX_PROPAGATION_STEPS,
        propagation_channels: int = 4,
        response_scale: str = "linear",
        damping: float = 0.5,
        contraction: float = 0.9,
        normalisation: str = "in_degree",
        cross_relation_aggregator: str = "mean",
        mixture_weighting: str = "sparsemax",
        node_type_index: Tensor | None = None,
        edge_signs: str = "graph",
        relation_gains: str = "per_relation",
        cell_class_weights: Tensor | None = None,
        channels_per_cell_class: int = 1,
        extracellular_pool_nodes: Tensor | None = None,
        shared_pool_classes: Tensor | None = None,
        num_descriptor_columns: int = 0,
        descriptor_columns_kept: Tensor | None = None,
        descriptor_treatment: str = "plain",
        seed_mask_partner_index: Tensor | None = None,
        entry_node_mask: Tensor | None = None,
        entry_relation_names: tuple[str, ...] = (),
        edge_weight: Tensor | None = None,
    ) -> None:
        super().__init__()
        self.register_buffer("seed_mask_partner_index", None if seed_mask_partner_index is None else seed_mask_partner_index.long(), persistent=False)
        entry_edges = None
        if entry_node_mask is None:
            if entry_relation_names:
                raise ValueError("entry relations need entry_node_mask, the nodes that exist only where they are seeded")
        else:
            if entry_node_mask.shape != (num_graph_nodes,):
                raise ValueError("entry_node_mask must have one entry per graph node")
            present_entry_relations = [name for name in entry_relation_names if name in relation_types]
            is_entry_edge = torch.zeros(len(edge_relation), dtype=torch.bool)
            for name in present_entry_relations:
                is_entry_edge |= edge_relation.long() == relation_types.index(name)
            if bool(is_entry_edge.any()) and not bool(entry_node_mask.bool()[edge_source.long()[is_entry_edge]].all()):
                raise ValueError("every edge of an entry relation must leave an entry node")
            weights = torch.ones(len(edge_relation)) if edge_weight is None else edge_weight.to(torch.float32)
            entry_edges = (edge_source.long()[is_entry_edge], edge_target.long()[is_entry_edge],
                           edge_sign.float()[is_entry_edge] * weights[is_entry_edge],
                           torch.as_tensor([present_entry_relations.index(relation_types[int(relation)]) for relation in edge_relation[is_entry_edge]], dtype=torch.long),
                           len(present_entry_relations))
            # the entry edges leave the arrays the normalised adjacency is built from; the relation list is unchanged
            edge_source, edge_target, edge_relation, edge_sign = edge_source[~is_entry_edge], edge_target[~is_entry_edge], edge_relation[~is_entry_edge], edge_sign[~is_entry_edge]
            if cofactor_edges is not None:
                cofactor_edges = cofactor_edges[~is_entry_edge]
        self.register_buffer("entry_node_mask", None if entry_node_mask is None else entry_node_mask.bool(), persistent=False)
        self.register_buffer("entry_edge_source", None if entry_edges is None else entry_edges[0], persistent=False)
        self.register_buffer("entry_edge_target", None if entry_edges is None else entry_edges[1], persistent=False)
        self.register_buffer("entry_edge_signed_weight", None if entry_edges is None else entry_edges[2], persistent=False)
        self.register_buffer("entry_edge_relation", None if entry_edges is None else entry_edges[3], persistent=False)
        num_entry_relations = 0 if entry_edges is None else entry_edges[4]
        if node_features.shape[0] != num_graph_nodes:
            raise ValueError("node_features must have one row per graph node")
        check_descriptor_treatment(descriptor_treatment, num_descriptor_columns, node_features.shape[1])
        self.descriptor_treatment = descriptor_treatment
        self.num_descriptor_columns = num_descriptor_columns
        if not 0.0 < damping <= 1.0 or not 0.0 < contraction < 1.0:
            raise ValueError("damping must be in (0, 1] and contraction in (0, 1)")
        self.num_graph_nodes = num_graph_nodes
        self.node_state_dim = node_state_dim
        self.num_propagation_steps = num_propagation_steps
        if fixed_point_tolerance is not None and fixed_point_tolerance <= 0:
            raise ValueError("fixed_point_tolerance must be positive, or None to run a fixed number of rounds")
        if max_propagation_steps < 1:
            raise ValueError("max_propagation_steps must be at least 1")
        self.fixed_point_tolerance = fixed_point_tolerance
        self.max_propagation_steps = max_propagation_steps
        self.propagation_steps_taken = 0  # the rounds the last call took, for the run record
        if class_weight_rule not in CLASS_WEIGHT_RULES:
            raise ValueError(f"class_weight_rule must be one of {CLASS_WEIGHT_RULES}, not {class_weight_rule!r}")
        self.class_weight_rule = class_weight_rule
        if input_weights not in INPUT_WEIGHT_FORMS:
            raise ValueError(f"input_weights must be one of {INPUT_WEIGHT_FORMS}, not {input_weights!r}")
        self.input_weights = input_weights
        self.cache_class_bottlenecks = cache_class_bottlenecks
        # compartment pH (see PH_SENSITIVITY_TARGET_ARGUMENT). Both tensors or neither; a node whose compartment index
        # is -1 (none, or several, as a transport reaction has) takes no pH weight.
        self.ph_compartment_names = list(ph_compartment_names or [])
        if (ph_proton_node_index is None) != (ph_compartment_of_node is None):
            raise ValueError("ph_proton_node_index and ph_compartment_of_node are given together or not at all")
        if ph_proton_node_index is None:
            self.ph_proton_node_index = self.ph_compartment_of_node = self.log_ph_sensitivity = None
        else:
            if ph_compartment_of_node.shape[0] != num_graph_nodes:
                raise ValueError("ph_compartment_of_node needs one entry per graph node")
            if int(ph_compartment_of_node.max()) >= ph_proton_node_index.shape[0]:
                raise ValueError("ph_compartment_of_node indexes a compartment with no proton node")
            if self.ph_compartment_names and len(self.ph_compartment_names) != ph_proton_node_index.shape[0]:
                raise ValueError("ph_compartment_names must name one compartment per proton node")
            self.register_buffer("ph_proton_node_index", ph_proton_node_index.long(), persistent=False)
            self.register_buffer("ph_compartment_of_node", ph_compartment_of_node.long(), persistent=False)
            self.log_ph_sensitivity = nn.Parameter(torch.zeros(ph_proton_node_index.shape[0]))
            # the half-width of each compartment's pH-activity curve, in pH units; one per compartment and learned,
            # because no per-reaction pH optimum exists at the scale of this graph
            self.log_ph_activity_width = nn.Parameter(
                torch.full((ph_proton_node_index.shape[0],), math.log(DEFAULT_PH_ACTIVITY_WIDTH)))
            self.ph_sensitivity_calibrated = False
            self.compartment_ph_sensitivity_report: dict[str, float] | None = None
        self.class_bottleneck_cache: dict[tuple[int, ...], Tensor] = {}  # seeds -> [N, K]; see cached_class_bottleneck
        self.damping = damping
        if cross_relation_aggregator not in CROSS_RELATION_AGGREGATORS:
            raise ValueError(f"cross_relation_aggregator must be one of {CROSS_RELATION_AGGREGATORS}, not {cross_relation_aggregator!r}")
        self.cross_relation_aggregator = cross_relation_aggregator
        uses_mixture = cross_relation_aggregator == "softmax_mixture"
        stacked_adjacency, relation_names = self.signed_stacked_adjacency(num_graph_nodes, relation_types, edge_source, edge_target, edge_relation, edge_sign,
                                                                          non_propagating_nodes, cofactor_edges, normalisation,
                                                                          divide_by_relations_feeding=not uses_mixture, edge_signs=edge_signs)
        self.relation_names = relation_names
        self.register_buffer("stacked_adjacency", stacked_adjacency, persistent=False)
        # the stacked rows that hold an edge, as their own matrix [K, N] with the relation and target of each row: with the
        # mean across relations the response sums the rows' messages times their relation's gain into their targets, so
        # the product needs only these K rows (54,982 of 700,435 on graph_full_neuronal) instead of all R * N
        coalesced = stacked_adjacency.coalesce()
        nonempty_rows, compressed_row = torch.unique(coalesced.indices()[0], return_inverse=True)
        self.register_buffer("compressed_adjacency", torch.sparse_coo_tensor(torch.stack([compressed_row, coalesced.indices()[1]]), coalesced.values(),
                                                                             (len(nonempty_rows), num_graph_nodes)).coalesce(), persistent=False)
        self.register_buffer("row_relation_index", torch.div(nonempty_rows, num_graph_nodes, rounding_mode="floor"), persistent=False)
        self.register_buffer("row_target_index", nonempty_rows % num_graph_nodes, persistent=False)
        # the edge structure the bottleneck pass walks, taken from the same adjacency so the two agree about which
        # edges exist; a bottleneck cares only whether a route exists, so signs, weights and the in-degree divisor
        # play no part in it
        self.register_buffer("bottleneck_edge_source", coalesced.indices()[1], persistent=False)
        self.register_buffer("bottleneck_edge_target", coalesced.indices()[0] % num_graph_nodes, persistent=False)
        # how many copies of the source a step consumes, which is that participant's stoichiometry. The graph does not
        # carry it on its edges yet, so this stays None and every demand is 1; see class_bottleneck.
        self.register_buffer("bottleneck_edge_demand", None, persistent=False)
        self.register_buffer("node_features", node_features.to(torch.float32), persistent=False)
        self.register_buffer("relation_is_unsigned", torch.tensor([name in UNSIGNED_RELATIONS and edge_signs == "graph" for name in relation_names]),
                             persistent=False)
        self.maximum_gain = contraction
        if cell_class_weights is not None:  # propagation_channels is then K * channels_per_cell_class; see the module docstring
            cell_class_weights = cell_class_weights.to(torch.float32)
            if cell_class_weights.dim() != 2 or cell_class_weights.shape[0] != num_graph_nodes:
                raise ValueError("cell_class_weights must be [num_graph_nodes, num_cell_classes]")
            if bool((cell_class_weights < 0).any()) or bool((cell_class_weights > 1).any()):
                raise ValueError("cell_class_weights must lie in [0, 1], or the iteration may not contract")
            self.num_cell_classes = cell_class_weights.shape[1]
            self.channels_per_cell_class = channels_per_cell_class
            propagation_channels = self.num_cell_classes * channels_per_cell_class
            self.register_buffer("channel_node_weight", cell_class_weights.repeat_interleave(channels_per_cell_class, dim=1), persistent=False)  # [N, P], class-major
            # the same weights before the per-class channels are repeated out of them, [N, K]: the limiting-step rule reads these,
            # because every channel of one class carries that class's weight and a bottleneck over identical weights is the same
            # bottleneck, so running the pass per class instead of per channel costs channels_per_cell_class times less
            self.register_buffer("cell_class_node_weight", cell_class_weights, persistent=False)
        else:
            self.num_cell_classes, self.channels_per_cell_class = 0, 0
            self.channel_node_weight = self.cell_class_node_weight = None
        if extracellular_pool_nodes is not None:
            if cell_class_weights is None:
                raise ValueError("extracellular_pool_nodes needs cell_class_weights")
            pooled_classes = torch.ones(self.num_cell_classes, dtype=torch.bool) if shared_pool_classes is None else shared_pool_classes.bool()
            if pooled_classes.shape[0] != self.num_cell_classes or int(pooled_classes.sum()) < 2:
                raise ValueError("shared_pool_classes must name at least two of the cell classes")
            self.register_buffer("extracellular_pool_nodes", extracellular_pool_nodes.bool(), persistent=False)
            self.register_buffer("shared_pool_classes", pooled_classes, persistent=False)
        else:
            self.extracellular_pool_nodes = None
            self.shared_pool_classes = None
        self.propagation_channels = propagation_channels
        if response_scale not in ("linear", "signed_log"):
            raise ValueError("response_scale must be 'linear' or 'signed_log'")
        self.response_scale = response_scale
        self.log_response_scale = nn.Parameter(torch.full((propagation_channels,), INITIAL_LOG_RESPONSE_SCALE)) if response_scale == "signed_log" else None
        if relation_gains not in RELATION_GAINS:
            raise ValueError(f"relation_gains must be one of {RELATION_GAINS}, not {relation_gains!r}")
        self.relation_gains = relation_gains
        initial_gain_logit = torch.randn(1 if relation_gains == "shared" else len(relation_names), propagation_channels)  # spread so channels start at different time scales
        for carrier_relation in (COSUBSTRATE_RELATION, COPRODUCT_RELATION):
            if carrier_relation in relation_names and relation_gains == "per_relation":  # carrier pools are recycled and buffered: coupling through them starts weak (gain about 0.04)
                initial_gain_logit[relation_names.index(carrier_relation)] += CARRIER_INITIAL_GAIN_LOGIT_OFFSET
        self.gain_logit = nn.Parameter(initial_gain_logit)
        if input_weights == "positive":  # one shared starting value, so no class starts out weighing more in the pool
            self.log_input_weight = nn.Parameter(torch.full((propagation_channels,), -0.5 * math.log(propagation_channels)))
        else:
            self.input_weight = nn.Parameter(torch.randn(propagation_channels) / propagation_channels**0.5)
        self.channel_expansion = nn.Parameter(torch.randn(propagation_channels, node_state_dim) / propagation_channels**0.5)
        if descriptor_treatment == "zero_init_slow":
            self.num_structural_columns = node_features.shape[1] - num_descriptor_columns
            self.output_gate = nn.Linear(self.num_structural_columns, node_state_dim)
            self.descriptor_gate = nn.Linear(num_descriptor_columns, node_state_dim, bias=False)
            nn.init.zeros_(self.descriptor_gate.weight)
        else:
            self.output_gate = nn.Linear(node_features.shape[1], node_state_dim)
            self.descriptor_gate = None
        if descriptor_treatment == "seed_masked":
            self.register_buffer("node_features_without_descriptors", features_without_descriptors(self.node_features, num_descriptor_columns, descriptor_columns_kept), persistent=False)
        if uses_mixture:
            rows_with_an_edge = torch.zeros(len(relation_names) * num_graph_nodes, dtype=torch.bool)
            rows_with_an_edge[stacked_adjacency.coalesce().indices()[0]] = True
            self.register_buffer("relation_feeds_node", rows_with_an_edge.view(len(relation_names), num_graph_nodes), persistent=False)
            # with no node types given, every node shares one set of mixture weights
            given_types = torch.zeros(num_graph_nodes, dtype=torch.long) if node_type_index is None else node_type_index.long()
            if given_types.shape[0] != num_graph_nodes:
                raise ValueError("node_type_index must have one entry per graph node")
            self.register_buffer("node_type_index", given_types, persistent=False)
            self.cross_relation_mixture = CrossRelationAggregator(int(given_types.max()) + 1, propagation_channels,
                                                                  weighting=mixture_weighting)
        else:
            self.cross_relation_mixture = None
        # created last and only with entry edges, so an encoder without them draws the parameters it always drew
        self.entry_gain_logit = nn.Parameter(torch.randn(num_entry_relations, propagation_channels)) if num_entry_relations else None

    @staticmethod
    def spectral_radius_of_unsigned_aggregate(num_graph_nodes: int, source: Tensor, target: Tensor, iterations: int = 200) -> float:
        """Largest eigenvalue magnitude of the unsigned aggregate adjacency M, where M[target, source] counts the edges
        from source to target over every relation. M is non-negative, so its Perron root is the largest eigenvalue and
        a power iteration from a positive start converges to it. This is what the spectral normalisation divides by:
        with entries scaled to give rho(sum_r |S_r|) = 1, the gain-weighted operator obeys
        rho(sum_r gain_r S_r) <= max_r |gain_r|, so the existing contraction bound on the gains still makes the
        iteration a contraction without dividing any node's input by its own in-degree."""
        aggregate = torch.sparse_coo_tensor(torch.stack([target.long(), source.long()]), torch.ones(len(source)),
                                            (num_graph_nodes, num_graph_nodes)).coalesce()
        vector = torch.full((num_graph_nodes,), num_graph_nodes ** -0.5)
        radius = 0.0
        for _ in range(iterations):
            product = torch.sparse.mm(aggregate, vector[:, None]).squeeze(1)
            norm = float(product.norm())
            if norm == 0.0:
                return 0.0
            vector, previous_radius = product / norm, radius
            radius = norm
            if abs(radius - previous_radius) <= 1e-9 * max(radius, 1.0):
                break
        return radius

    @staticmethod
    def signed_stacked_adjacency(num_graph_nodes: int, relation_types: list[str], edge_source: Tensor, edge_target: Tensor,
                                 edge_relation: Tensor, edge_sign: Tensor, non_propagating_nodes: Tensor | None = None,
                                 cofactor_edges: Tensor | None = None, normalisation: str = "in_degree",
                                 divide_by_relations_feeding: bool = True, edge_signs: str = "graph") -> tuple[Tensor, list[str]]:
        """All relations' signed, normalised adjacencies stacked into one sparse [R * N, N] matrix, with the derived
        depletes_substrate relation appended, carrier edges moved to relations of their own when cofactor_edges is
        given, and edges leaving a non-propagating node dropped; returns the matrix and the relation names in stack order."""
        edge_relation = edge_relation.long().clone()
        relation_names = list(relation_types)
        if cofactor_edges is not None and SUBSTRATE_RELATION in relation_types and PRODUCT_RELATION in relation_types:
            cofactor_edges = cofactor_edges.bool()
            for main_relation, carrier_relation in ((SUBSTRATE_RELATION, COSUBSTRATE_RELATION), (PRODUCT_RELATION, COPRODUCT_RELATION)):
                relation_names.append(carrier_relation)
                edge_relation[cofactor_edges & (edge_relation == relation_types.index(main_relation))] = len(relation_names) - 1
        sources, targets, relations, signs = [edge_source.long()], [edge_target.long()], [edge_relation], [edge_sign.float()]
        for forward_relation, depletion_relation in ((SUBSTRATE_RELATION, DEPLETES_SUBSTRATE_RELATION),):  # carriers get none (module docstring)
            if forward_relation not in relation_names:
                continue
            forward_edges = edge_relation == relation_names.index(forward_relation)
            relation_names.append(depletion_relation)
            sources.append(edge_target[forward_edges].long())
            targets.append(edge_source[forward_edges].long())
            relations.append(torch.full((int(forward_edges.sum()),), len(relation_names) - 1, dtype=torch.long))
            signs.append(-torch.ones(int(forward_edges.sum())))
        source, target, relation, sign = torch.cat(sources), torch.cat(targets), torch.cat(relations), torch.cat(signs)
        if non_propagating_nodes is not None:
            keep = ~non_propagating_nodes.bool()[source]
            source, target, relation, sign = source[keep], target[keep], relation[keep], sign[keep]
        for index, name in enumerate(relation_names):
            if name in UNSIGNED_RELATIONS:
                sign = torch.where(relation == index, torch.ones_like(sign), sign)  # the learned gain carries the sign
        if edge_signs not in EDGE_SIGNS:
            raise ValueError(f"edge_signs must be one of {EDGE_SIGNS}, not {edge_signs!r}")
        if edge_signs == "all_positive":  # the sign ablation; see EDGE_SIGNS
            sign = torch.ones_like(sign)
        elif edge_signs == "permuted":  # the sign permutation; see EDGE_SIGNS
            sign = sign[torch.randperm(len(sign), generator=torch.Generator().manual_seed(SIGN_PERMUTATION_SEED))]
        if normalisation not in NORMALISATIONS:
            raise ValueError(f"normalisation must be one of {NORMALISATIONS}, not {normalisation!r}")
        row = relation * num_graph_nodes + target
        if normalisation == "spectral":  # one global scale and no per-node divisor; see NORMALISATIONS
            radius = LinearResponseEncoder.spectral_radius_of_unsigned_aggregate(num_graph_nodes, source, target)
            values = sign / max(radius, 1.0)
        elif normalisation == "total_in_degree":  # one per-node divisor over all relations together; see NORMALISATIONS
            total_in_degree = torch.zeros(num_graph_nodes).index_add(0, target, torch.ones(len(target)))
            values = sign / total_in_degree[target]
        else:
            in_degree = torch.zeros(len(relation_names) * num_graph_nodes).index_add(0, row, torch.ones(len(row)))
            relations_feeding_node = (in_degree.view(len(relation_names), num_graph_nodes) > 0).sum(dim=0).clamp_min(1)
            cross_relation_divisor = relations_feeding_node[target] if divide_by_relations_feeding else 1.0
            values = sign / (in_degree[row] * cross_relation_divisor)
        stacked = torch.sparse_coo_tensor(torch.stack([row, source]), values, (len(relation_names) * num_graph_nodes, num_graph_nodes)).coalesce()
        return stacked, relation_names

    def relation_gain(self) -> Tensor:
        """[R, D] gains: positive for signed relations, either sign for unsigned ones, all below the contraction bound."""
        gain_logit = self.gain_logit.expand(len(self.relation_names), -1)  # a no-op unless relation_gains is "shared"
        positive = torch.sigmoid(gain_logit)
        either_sign = torch.tanh(gain_logit)
        return self.maximum_gain * torch.where(self.relation_is_unsigned[:, None], either_sign, positive)

    def sustained_input(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor) -> Tensor:
        """u: [B, N, C], sign * magnitude * the channel's input weight at each perturbed node (index -1 is padding)."""
        batch_size = perturbation_node_index.shape[0]
        valid = (perturbation_node_index >= 0).to(perturbation_sign_and_magnitude.dtype)
        signed_magnitude = perturbation_sign_and_magnitude[:, :, 0] * perturbation_sign_and_magnitude[:, :, 1] * valid
        injected = signed_magnitude[:, :, None] * self.injection_weights()[None, None, :]
        field = torch.zeros(batch_size, self.num_graph_nodes, self.propagation_channels, device=injected.device, dtype=injected.dtype)
        index = perturbation_node_index.clamp_min(0)[:, :, None].expand(-1, -1, self.propagation_channels)
        return field.scatter_add(1, index, injected)

    def compartment_ph_weight(self, state: Tensor) -> Tensor:
        """[N, B, C] multiplier from the settled proton deviation of each node's compartment, 1 where a node has none.

        state is the settled response, node-major. The deviation is read per channel, so each cell class carries its own
        proton balance, and it is taken as a fractional change in [H+] so that the shift in pH is its negative base-ten
        logarithm (see PH_SENSITIVITY_TARGET_ARGUMENT above). Bounded in (0, 1], so it cannot raise a norm, and applied
        after the fixed point, so the contraction is untouched.

        One property is weaker on the log scale than it was on the linear one, and it is the registered cost. The
        linear form was even in the deviation, so a gain-of-function was the exact negative of a loss-of-function. A
        pH shift is not even: acidifying by a given fraction and alkalinising by it move pH by different amounts,
        which is what the logarithm is for and is the physiology. So a gain and a loss now take different weights at a
        node whose compartment the perturbation moves, and the two stop being exact negatives there. Positive
        homogeneity was already broken by the linear form.
        """
        deviation = state[self.ph_proton_node_index]  # [compartments, B, C]
        fraction = torch.exp(self.log_ph_sensitivity)[:, None, None] * deviation  # the fractional change in [H+]
        change_in_ph = -torch.log10(torch.clamp(1.0 + fraction, min=1.0 + SMALLEST_PROTON_FRACTION))
        per_compartment = torch.exp(-(change_in_ph / torch.exp(self.log_ph_activity_width)[:, None, None]) ** 2)
        ones = torch.ones_like(per_compartment[:1])
        with_none = torch.cat([per_compartment, ones], dim=0)  # the last row is the weight of a node in no compartment
        has_none = torch.full_like(self.ph_compartment_of_node, per_compartment.shape[0])
        return with_none[torch.where(self.ph_compartment_of_node >= 0, self.ph_compartment_of_node, has_none)]

    def calibrate_ph_sensitivity(self, state: Tensor) -> dict[str, float]:
        """Set each compartment's sensitivity so that the root mean square of (sensitivity * deviation) on this settled
        state is PH_SENSITIVITY_TARGET_ARGUMENT, then let training move it.

        The model's deviations are in arbitrary units set by the input magnitudes and the gains, so there is no scale to
        derive one from; without this the weight starts at exactly 1 and has no gradient to learn from. Called once, on
        one batch, before training, and the values it sets are reported in the run record. A compartment whose proton
        node the batch does not reach keeps its sensitivity, because there is nothing to calibrate against."""
        if self.ph_proton_node_index is None:
            return {}
        with torch.no_grad():
            deviation = state[self.ph_proton_node_index]  # [compartments, B, C]
            root_mean_square = deviation.pow(2).mean(dim=(1, 2)).sqrt()
            reached = root_mean_square > PH_WEIGHT_EPSILON
            self.log_ph_sensitivity[reached] = torch.log(PH_SENSITIVITY_TARGET_ARGUMENT / root_mean_square[reached])
            self.ph_sensitivity_calibrated = True
            names = self.ph_compartment_names or [str(index) for index in range(len(root_mean_square))]
            self.compartment_ph_sensitivity_report = {
                name: float(torch.exp(self.log_ph_sensitivity[index]))
                for index, name in enumerate(names) if bool(reached[index])}
            return self.compartment_ph_sensitivity_report

    def injection_weights(self) -> Tensor:
        """[C]: the scale of the perturbation in each channel, positive under input_weights "positive"."""
        return torch.exp(self.log_input_weight) if self.input_weights == "positive" else self.input_weight

    def class_bottleneck(self, perturbation_node_index: Tensor, rounds: int | None = None, chunk: int = BOTTLENECK_CHUNK) -> Tensor:
        """[N, B, K]: per cell class, the bottleneck of the least-limited route from the perturbed nodes to each node.

        The max-min recurrence, run until it stops changing (rounds None, the steady state; see CLASS_WEIGHT_RULES) or
        for a given number of rounds:

            b(t)   = w(t) at a perturbed node, 0 elsewhere
            b(t) <- max( b(t), min( w(t), max over edges s -> t of b(s) / demand(s -> t) ) )

        then, at each extracellular pool node, b of every pooled class is the largest over the pooled classes, and an
        entry node the perturbation does not seed has w = 0. b only grows and is bounded by 1, and a widest route is a
        simple path, so the recurrence settles within num_graph_nodes rounds. demand is how many copies of the source
        a step consumes, which is the stoichiometry of that participant: a receptor that needs two molecules of its
        transmitter is limited at half the transmitter's availability (the owner, 13 October 2026: "we also have to
        incorporate things like ACh needing to bind two molecules to the receptor in order to activate it, which would
        cut that expression in half for the final response amount"). The graph does not carry those coefficients on its
        edges yet - Reactome's and Human-GEM's SBML both state them and the build reads them for the membrane
        potential, then discards them - so bottleneck_edge_demand is None and every demand is 1 until a build writes
        them (set_bottleneck_edge_demand).

        No parameter enters this, only the fixed buffers and which nodes are perturbed, so it carries no gradient and
        runs under no_grad. The pass is per class rather than per channel, because every channel of a class shares that
        class's weights, and the edge gather is taken chunk columns at a time to bound its memory.
        """
        if self.cell_class_node_weight is None:
            raise ValueError("class_bottleneck needs cell_class_weights")
        with torch.no_grad():
            weights = self.cell_class_node_weight  # [N, K]
            batch_size = perturbation_node_index.shape[0]
            num_classes = weights.shape[1]
            node_limit = weights[:, None, :].expand(-1, batch_size, -1)  # [N, B, K]
            presence = self.present_nodes(perturbation_node_index)  # [N, B, 1], or None without entry nodes
            if presence is not None:
                node_limit = node_limit * presence
            seeded = torch.zeros((batch_size, self.num_graph_nodes + 1), dtype=torch.bool, device=weights.device)
            seeded.scatter_(1, torch.where(perturbation_node_index >= 0, perturbation_node_index,
                                           torch.full_like(perturbation_node_index, self.num_graph_nodes)), True)
            seeded = seeded[:, : self.num_graph_nodes].permute(1, 0)[:, :, None]  # [N, B, 1]
            bottleneck = self.share_bottleneck_at_the_pool(torch.where(seeded, node_limit, torch.zeros_like(node_limit)))
            source, target = self.bottleneck_edge_source, self.bottleneck_edge_target
            demand = None if self.bottleneck_edge_demand is None else self.bottleneck_edge_demand[:, None]
            flat_width = batch_size * num_classes
            flat_limit = node_limit.reshape(self.num_graph_nodes, flat_width)
            flat = bottleneck.reshape(self.num_graph_nodes, flat_width)
            for _ in range(self.num_graph_nodes if rounds is None else rounds):
                arriving = torch.zeros_like(flat)
                for start in range(0, flat_width, chunk):
                    per_edge = flat[:, start : start + chunk][source]
                    if demand is not None:
                        per_edge = per_edge / demand
                    arriving[:, start : start + chunk] = arriving[:, start : start + chunk].index_reduce(
                        0, target, per_edge, "amax", include_self=True)
                updated = torch.maximum(flat, torch.minimum(flat_limit, arriving))
                updated = self.share_bottleneck_at_the_pool(updated.reshape(self.num_graph_nodes, batch_size, num_classes)
                                                            ).reshape(self.num_graph_nodes, flat_width)
                if torch.equal(updated, flat):
                    break
                flat = updated
            return flat.reshape(self.num_graph_nodes, batch_size, num_classes)

    def share_bottleneck_at_the_pool(self, bottleneck: Tensor) -> Tensor:
        """[N, B, K] with each extracellular pool node's pooled classes set to their largest value: the pool is one
        compartment, so a route into it in one class continues out of it in another (see CLASS_WEIGHT_RULES)."""
        if self.extracellular_pool_nodes is None:
            return bottleneck
        pooled_largest = bottleneck[:, :, self.shared_pool_classes].amax(dim=2, keepdim=True).expand_as(bottleneck)
        replace = self.extracellular_pool_nodes[:, None, None] & self.shared_pool_classes[None, None, :]
        return torch.where(replace, pooled_largest, bottleneck)

    def cached_class_bottleneck(self, perturbation_node_index: Tensor) -> Tensor:
        """[N, B, K]: class_bottleneck at the steady state, computed once per set of seeds and kept.

        The steady-state bottleneck depends only on the graph, the class weights, the demands and the perturbation's own
        seeds, all fixed while a model trains, so every epoch after the first reads it from here. Memory: one [N, K]
        float tensor per distinct set of seeds, on the device of the weights. The cache is emptied when the module moves
        device or dtype and when the demands change."""
        keys = [tuple(sorted({int(node) for node in row.tolist() if node >= 0})) for row in perturbation_node_index]
        missing = [key for key in dict.fromkeys(keys) if key not in self.class_bottleneck_cache]
        if missing:
            widest = max(1, max(len(key) for key in missing))
            missing_index = torch.full((len(missing), widest), -1, dtype=torch.long, device=perturbation_node_index.device)
            for row, key in enumerate(missing):
                missing_index[row, : len(key)] = torch.as_tensor(key, dtype=torch.long)
            computed = self.class_bottleneck(missing_index)  # [N, M, K]
            for row, key in enumerate(missing):
                self.class_bottleneck_cache[key] = computed[:, row, :].clone()
        return torch.stack([self.class_bottleneck_cache[key] for key in keys], dim=1)

    def clear_class_bottleneck_cache(self) -> None:
        self.class_bottleneck_cache.clear()

    def set_bottleneck_edge_demand(self, demand: Tensor | None) -> None:
        """Per adjacency entry (bottleneck_edge_source, bottleneck_edge_target), how many copies of the source the step
        consumes; at least 1. Empties the cache, since every bottleneck in it was taken under the old demands."""
        if demand is not None:
            if demand.shape != self.bottleneck_edge_source.shape:
                raise ValueError("demand needs one value per adjacency entry of the bottleneck pass")
            if bool((demand < 1).any()):
                raise ValueError("a demand below 1 would raise availability past a step, and b would leave [0, 1]")
            demand = demand.to(self.bottleneck_edge_source.device, torch.float32)
        self.bottleneck_edge_demand = demand
        self.clear_class_bottleneck_cache()

    def _apply(self, fn, *args, **kwargs):
        self.clear_class_bottleneck_cache()  # cached tensors do not move with the buffers
        return super()._apply(fn, *args, **kwargs)

    def response(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, num_steps: int | None = None) -> Tensor:
        """The propagated response h, [B, N, C], before expansion and gating.

        With fixed_point_tolerance set (amendment item 7 of 12 October 2026; the owner: "Use the convergence instead of
        8 steps. 8 steps is arbitrary, not physiological") the iteration runs until the relative change of the state
        falls below that tolerance, to at most max_propagation_steps rounds, and the number of rounds taken is kept in
        propagation_steps_taken for the run record. Otherwise it runs num_propagation_steps rounds, which is what every
        run before that date did. num_steps overrides both, for the measurements that read the state at a chosen round.

        The relative change is ||h(t+1) - h(t)|| / max(||h(t)||, eps) in the Frobenius norm over the whole batch, so one
        number decides for the batch and a perturbation with a large field cannot be declared settled by a neighbour with
        a small one. The update is a contraction by construction (the gains are bounded and the adjacency is normalised),
        without which there is no fixed point to reach; the bound is what makes the stopping rule a statement about the
        steady state rather than about the arithmetic running out."""
        sustained = self.sustained_input(perturbation_node_index, perturbation_sign_and_magnitude)
        batch_size = sustained.shape[0]
        gain = self.relation_gain()  # [R, C]
        num_relations = gain.shape[0]
        node_major_input = sustained.permute(1, 0, 2).contiguous()  # [N, B, C]
        # the product rule weights the input and the messages inside the loop; the limiting-step rule leaves the loop
        # unweighted and multiplies the settled state by the bottleneck once, below. See CLASS_WEIGHT_RULES.
        node_weight = None if self.channel_node_weight is None or self.class_weight_rule != "product" else self.channel_node_weight[:, None, :]  # [N, 1, C]
        if node_weight is not None:
            node_major_input = node_major_input * node_weight  # a class that does not express the perturbed node does not feel it there
        state = node_major_input  # h(0) = u, so after k steps the response reaches k edges from the perturbed nodes
        previous_state = state
        presence = self.present_nodes(perturbation_node_index)  # [N, B, 1], or None without entry nodes
        entry_edge_gain = None
        if self.entry_gain_logit is not None and self.entry_edge_source.numel():
            # sign * weight * gain per entry edge and channel; the gain is positive, so the sign is the mechanism's
            entry_edge_gain = (self.entry_edge_signed_weight[:, None] * torch.sigmoid(self.entry_gain_logit)[self.entry_edge_relation])[:, None, :]  # [E, 1, C]
        row_gain = gain[self.row_relation_index][:, None, :] if self.cross_relation_mixture is None else None  # [K, 1, C]
        rounds = (self.num_propagation_steps if num_steps is None else num_steps) if (num_steps is not None or self.fixed_point_tolerance is None) \
            else self.max_propagation_steps
        steps_taken = 0
        for _ in range(rounds):
            if row_gain is not None:  # the sum over relations of gain times the relation's messages, over the rows that hold an edge
                row_messages = torch.sparse.mm(self.compressed_adjacency, state.reshape(self.num_graph_nodes, batch_size * self.propagation_channels))
                relation_messages = torch.zeros_like(state).index_add(0, self.row_target_index, row_messages.view(-1, batch_size, self.propagation_channels) * row_gain)
            else:
                aggregated = torch.sparse.mm(self.stacked_adjacency, state.reshape(self.num_graph_nodes, batch_size * self.propagation_channels))
                per_relation = aggregated.view(num_relations, self.num_graph_nodes, batch_size, self.propagation_channels) * gain[:, None, None, :]
                relation_messages = self.cross_relation_mixture(per_relation, self.relation_feeds_node, self.node_type_index)
            if entry_edge_gain is not None:  # added on their own: no in-degree divisor, and no part in the mean or the mixture across relations
                relation_messages = relation_messages.index_add(0, self.entry_edge_target, state[self.entry_edge_source] * entry_edge_gain)
            if node_weight is not None:  # every mixture statistic is positively homogeneous, so scaling after it equals scaling the messages
                relation_messages = relation_messages * node_weight
            state = (1.0 - self.damping) * state + self.damping * (relation_messages + node_major_input)
            if self.extracellular_pool_nodes is not None:
                state = self.share_extracellular_pool(state)
            if presence is not None:  # an entry node that is not a seed is absent: it holds nothing and so sends nothing
                state = state * presence
            steps_taken += 1
            if self.fixed_point_tolerance is not None and num_steps is None:
                change = torch.linalg.vector_norm(state - previous_state) / torch.linalg.vector_norm(previous_state).clamp_min(FIXED_POINT_NORM_EPSILON)
                if float(change) < self.fixed_point_tolerance:
                    break
            previous_state = state
        self.propagation_steps_taken = steps_taken
        if self.class_weight_rule == "limiting_step" and self.cell_class_node_weight is not None:
            # at the steady state the bottleneck over every route, kept per perturbation; a transient reading (a chosen
            # round, or a fixed step count) takes it over the rounds the propagation took. repeat_interleave puts the
            # per-class value on that class's channels. See CLASS_WEIGHT_RULES.
            if num_steps is None and self.fixed_point_tolerance is not None:
                bottleneck = self.cached_class_bottleneck(perturbation_node_index) if self.cache_class_bottlenecks \
                    else self.class_bottleneck(perturbation_node_index)
            else:
                bottleneck = self.class_bottleneck(perturbation_node_index, rounds=steps_taken)
            state = state * bottleneck.repeat_interleave(self.channels_per_cell_class, dim=2)
        if self.ph_proton_node_index is not None:  # read from the settled state, so after the loop and after no bottleneck
            state = state * self.compartment_ph_weight(state)
        return state.permute(1, 0, 2)

    def present_nodes(self, perturbation_node_index: Tensor) -> Tensor | None:
        """[num_graph_nodes, batch_size, 1], zero for an entry node the perturbation is not seeded on; None without
        entry nodes. Padding (index -1) marks no node."""
        if self.entry_node_mask is None:
            return None
        batch_size = perturbation_node_index.shape[0]
        seeded = torch.zeros((batch_size, self.num_graph_nodes + 1), dtype=torch.bool, device=perturbation_node_index.device)
        seeded.scatter_(1, torch.where(perturbation_node_index >= 0, perturbation_node_index, torch.full_like(perturbation_node_index, self.num_graph_nodes)), True)
        present = ~self.entry_node_mask[None, :] | seeded[:, : self.num_graph_nodes]
        return present.permute(1, 0)[:, :, None].to(torch.float32)

    def share_extracellular_pool(self, state: Tensor) -> Tensor:
        """state [N, B, C] with the pool nodes' channels of the pooled classes replaced by their mean over those classes."""
        num_nodes, batch_size, _ = state.shape
        by_class = state.view(num_nodes, batch_size, self.num_cell_classes, self.channels_per_cell_class)
        pooled_mean = by_class[:, :, self.shared_pool_classes, :].mean(dim=2, keepdim=True).expand_as(by_class)
        replace = self.extracellular_pool_nodes[:, None, None, None] & self.shared_pool_classes[None, None, :, None]
        return torch.where(replace, pooled_mean, by_class).reshape(num_nodes, batch_size, -1)

    def forward(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:
        """The gated response field [B, N, D]; relation_adjacencies is accepted for the message passing signature and
        ignored, since the signed adjacency is built once from the edges given at construction."""
        gate = self.node_gate(perturbation_node_index)  # [1, N, D], or [B, N, D] under seed_masked
        response = self.response(perturbation_node_index, perturbation_sign_and_magnitude)
        if self.response_scale == "signed_log":
            response = torch.sign(response) * torch.log1p(response.abs() / torch.exp(self.log_response_scale))
        return (response @ self.channel_expansion) * gate

    def descriptor_parameters(self) -> list[nn.Parameter]:
        """The parameters that read the descriptor columns on their own (zero_init_slow), for their own learning rate."""
        return [] if self.descriptor_gate is None else list(self.descriptor_gate.parameters())

    def node_gate(self, perturbation_node_index: Tensor) -> Tensor:
        """sigmoid of the gate logit per node, [1, N, D]; under seed_masked [B, N, D], with each perturbation's perturbed
        nodes gated on their structural columns only."""
        if self.descriptor_gate is not None:
            logit = (self.output_gate(self.node_features[:, : self.num_structural_columns])
                     + self.descriptor_gate(self.node_features[:, self.num_structural_columns :]))
            return torch.sigmoid(logit)[None, :, :]
        gate = torch.sigmoid(self.output_gate(self.node_features))[None, :, :]
        if self.descriptor_treatment != "seed_masked":
            return gate
        gate_without_descriptors = torch.sigmoid(self.output_gate(self.node_features_without_descriptors))[None, :, :]
        seed_rows = seed_unit_mask(perturbation_node_index, self.num_graph_nodes, self.seed_mask_partner_index)
        return torch.where(seed_rows[:, :, None], gate_without_descriptors, gate)

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:
        """The response is linear in the input, so the unperturbed field is zero and the difference field is the field."""
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude)
