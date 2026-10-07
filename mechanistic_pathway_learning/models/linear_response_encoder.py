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

The system is linear in the input, so the unperturbed response is zero and the field is a difference field by
construction: no base state, and no node identity, enters it. Node kinds enter only through an output gate,
field = (h W_expand) * sigmoid(node_features W + b), which lets a pooled readout tell where the response landed
without adding anything where it did not.
"""
from __future__ import annotations

import math

import torch
from torch import Tensor, nn

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
EDGE_SIGNS = ("graph", "all_positive")

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
        propagation_channels: int = 4,
        response_scale: str = "linear",
        damping: float = 0.5,
        contraction: float = 0.9,
        normalisation: str = "in_degree",
        cross_relation_aggregator: str = "mean",
        mixture_weighting: str = "sparsemax",
        node_type_index: Tensor | None = None,
        edge_signs: str = "graph",
    ) -> None:
        super().__init__()
        if node_features.shape[0] != num_graph_nodes:
            raise ValueError("node_features must have one row per graph node")
        if not 0.0 < damping <= 1.0 or not 0.0 < contraction < 1.0:
            raise ValueError("damping must be in (0, 1] and contraction in (0, 1)")
        self.num_graph_nodes = num_graph_nodes
        self.node_state_dim = node_state_dim
        self.num_propagation_steps = num_propagation_steps
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
        self.register_buffer("node_features", node_features.to(torch.float32), persistent=False)
        self.register_buffer("relation_is_unsigned", torch.tensor([name in UNSIGNED_RELATIONS and edge_signs == "graph" for name in relation_names]),
                             persistent=False)
        self.maximum_gain = contraction
        self.propagation_channels = propagation_channels
        if response_scale not in ("linear", "signed_log"):
            raise ValueError("response_scale must be 'linear' or 'signed_log'")
        self.response_scale = response_scale
        self.log_response_scale = nn.Parameter(torch.full((propagation_channels,), INITIAL_LOG_RESPONSE_SCALE)) if response_scale == "signed_log" else None
        initial_gain_logit = torch.randn(len(relation_names), propagation_channels)  # spread so channels start at different time scales
        for carrier_relation in (COSUBSTRATE_RELATION, COPRODUCT_RELATION):
            if carrier_relation in relation_names:  # carrier pools are recycled and buffered: coupling through them starts weak (gain about 0.04)
                initial_gain_logit[relation_names.index(carrier_relation)] += CARRIER_INITIAL_GAIN_LOGIT_OFFSET
        self.gain_logit = nn.Parameter(initial_gain_logit)
        self.input_weight = nn.Parameter(torch.randn(propagation_channels) / propagation_channels**0.5)
        self.channel_expansion = nn.Parameter(torch.randn(propagation_channels, node_state_dim) / propagation_channels**0.5)
        self.output_gate = nn.Linear(node_features.shape[1], node_state_dim)
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
        positive = torch.sigmoid(self.gain_logit)
        either_sign = torch.tanh(self.gain_logit)
        return self.maximum_gain * torch.where(self.relation_is_unsigned[:, None], either_sign, positive)

    def sustained_input(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor) -> Tensor:
        """u: [B, N, C], sign * magnitude * input_weight at each perturbed node (index -1 is padding)."""
        batch_size = perturbation_node_index.shape[0]
        valid = (perturbation_node_index >= 0).to(perturbation_sign_and_magnitude.dtype)
        signed_magnitude = perturbation_sign_and_magnitude[:, :, 0] * perturbation_sign_and_magnitude[:, :, 1] * valid
        injected = signed_magnitude[:, :, None] * self.input_weight[None, None, :]
        field = torch.zeros(batch_size, self.num_graph_nodes, self.propagation_channels, device=injected.device, dtype=injected.dtype)
        index = perturbation_node_index.clamp_min(0)[:, :, None].expand(-1, -1, self.propagation_channels)
        return field.scatter_add(1, index, injected)

    def response(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, num_steps: int | None = None) -> Tensor:
        """The propagated response h after num_steps (default num_propagation_steps), [B, N, C], before expansion and gating."""
        sustained = self.sustained_input(perturbation_node_index, perturbation_sign_and_magnitude)
        batch_size = sustained.shape[0]
        gain = self.relation_gain()  # [R, C]
        num_relations = gain.shape[0]
        node_major_input = sustained.permute(1, 0, 2).contiguous()  # [N, B, C]
        state = node_major_input  # h(0) = u, so after k steps the response reaches k edges from the perturbed nodes
        for _ in range(self.num_propagation_steps if num_steps is None else num_steps):
            aggregated = torch.sparse.mm(self.stacked_adjacency, state.reshape(self.num_graph_nodes, batch_size * self.propagation_channels))
            per_relation = aggregated.view(num_relations, self.num_graph_nodes, batch_size, self.propagation_channels) * gain[:, None, None, :]
            relation_messages = (per_relation.sum(dim=0) if self.cross_relation_mixture is None
                                 else self.cross_relation_mixture(per_relation, self.relation_feeds_node, self.node_type_index))
            state = (1.0 - self.damping) * state + self.damping * (relation_messages + node_major_input)
        return state.permute(1, 0, 2)

    def forward(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        """The gated response field [B, N, D]; relation_adjacencies is accepted for the message passing signature and
        ignored, since the signed adjacency is built once from the edges given at construction."""
        gate = torch.sigmoid(self.output_gate(self.node_features))  # [N, D]
        response = self.response(perturbation_node_index, perturbation_sign_and_magnitude)
        if self.response_scale == "signed_log":
            response = torch.sign(response) * torch.log1p(response.abs() / torch.exp(self.log_response_scale))
        return (response @ self.channel_expansion) * gate[None, :, :]

    def perturbation_difference_field(self, perturbation_node_index: Tensor, perturbation_sign_and_magnitude: Tensor, relation_adjacencies=None) -> Tensor:  # noqa: ARG002
        """The response is linear in the input, so the unperturbed field is zero and the difference field is the field."""
        return self.forward(perturbation_node_index, perturbation_sign_and_magnitude)
