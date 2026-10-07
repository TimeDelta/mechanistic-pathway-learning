"""Cofactor (carrier) edges of the metabolic layer, for the linear-response encoder (design section 5.2).

In the stoichiometric graph a carrier such as tetrahydrobiopterin is an ordinary substrate of every reaction that uses
it, so a linear response couples reactions through it: a PAH loss of function leaves tetrahydrobiopterin unconsumed,
raises it and through it raises tyrosine hydroxylase flux, L-dopa and dopamine, the opposite of what phenylketonuria
shows. Currency metabolites (tag_currency_metabolites.py) already stop passing the response on; carriers are not
currency (they are specific and few reactions use them), so their edges get relations of their own instead
(cosubstrate_of and coproduct_of, with no depletion edge since a recycled pool barely moves when one consumer slows), whose gains the model learns apart from those of
the main substrates and products.

An edge between a metabolite and a reaction is a carrier edge when either
  - the metabolite's name is on the carrier list below (conserved moieties that recur in few reactions, so a count
    cannot find them: biopterins, quinones, cytochromes, thioredoxins, ferredoxins, ascorbate, lipoamide, B6 and
    one-carbon folate carriers), or
  - within that reaction the metabolite forms a substrate-product pair that occurs in at least
    minimum_pair_reactions reactions of the reconstruction (on the slice and the full graph: PAPS -> PAP in 92,
    2-oxoglutarate <-> glutamate in 40 and 26 transaminations, aminoacyl-tRNA -> tRNA in 24 each); a recurring pair
    is a carrier pair by construction, while a main pair such as phenylalanine -> tyrosine occurs in 4.
The rule is per edge, not per metabolite: glutamate is a carrier in a transamination and a main substrate elsewhere.

carrier_rule selects which of the two clauses applies, so the two can be told apart by measurement rather than by
argument (design section 5.7):
  "all"       both clauses, the behaviour the layer was built with;
  "all_but_transmitters"
              both clauses, except that the pair clause never demotes a metabolite the graph releases as a transmitter
              or reads out as a laboratory label (TRANSMITTER_METABOLITE_NAMES). On the neuronal graph the pair clause
              marks 1758 edges of which 936 are aminoacyl-tRNA and 184 PAPS <-> PAP, both carriers by any reading, and
              178 are 2-oxoglutarate <-> glutamate; this arm keeps the first two and releases the third, so it
              separates the heuristic's correct calls from the one that touches a readout;
  "names"     the curated list only, with the recurring-pair clause off. The pair clause is the heuristic half, and it
              catches 2-oxoglutarate <-> glutamate in 40 reactions, which makes glutamate a carrier in every
              transamination; glutamate is also a transmitter and a laboratory label, so this arm asks whether the
              heuristic suppresses a response the experiment needs;
  "neuronal"  the carriers that gate transmitter synthesis and the oxidative response only (biopterins for the
              hydroxylases, B6 for the decarboxylases, ascorbate for dopamine beta-hydroxylase, the one-carbon folate
              carriers for methyl supply, thioredoxins for the oxidative arm), pair clause off. Narrowing the list is
              not free in one direction: a metabolite off the list goes back to being a main substrate, so this arm
              reintroduces exactly the spurious coupling the layer exists to stop, in the respiratory chain
              (ubiquinone, the cytochromes) and in steroidogenesis (the ferredoxins). It is the arm that tests whether
              that coupling costs less than treating the transmitter carriers and the respiratory carriers alike.
"""
from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np

DEFAULT_CARRIER_METABOLITE_NAMES: frozenset[str] = frozenset(
    name.lower() for name in (
        "tetrahydrobiopterin", "dihydrobiopterin", "quinonoid dihydrobiopterin", "4alpha-hydroxytetrahydrobiopterin",
        "O2-4a-cyclic-tetrahydrobiopterin",
        "ubiquinone", "ubiquinol",
        "ferricytochrome C", "ferrocytochrome C", "ferricytochrome B5", "ferrocytochrome B5",
        "thioredoxin", "oxidized thioredoxin", "mitothioredoxin", "mitooxidized thioredoxin",
        "oxidized ferredoxin", "reduced ferredoxin", "oxidized adrenal ferredoxin", "reduced adrenal ferredoxin",
        "ascorbate", "monodehydroascorbate", "dehydroascorbic acid",
        "lipoamide", "dihydrolipoamide",
        "[ACP]", "acyl carrier protein",
        "pyridoxal-phosphate", "pyridoxamine-phosphate",
        "5,10-methylene-THF", "5-methyl-THF", "10-formyl-THF", "5,10-methenyl-THF", "5-formyl-THF",
        "5-formiminotetrahydrofolate", "dihydrofolate",
    )
)
DEFAULT_MINIMUM_PAIR_REACTIONS = 10
# the carriers that gate transmitter synthesis and the oxidative response; see carrier_rule in the module docstring
NEURONAL_CARRIER_METABOLITE_NAMES: frozenset[str] = frozenset(
    name.lower() for name in (
        # tetrahydrobiopterin, the cofactor of phenylalanine, tyrosine and tryptophan hydroxylase
        "tetrahydrobiopterin", "dihydrobiopterin", "quinonoid dihydrobiopterin", "4alpha-hydroxytetrahydrobiopterin",
        "O2-4a-cyclic-tetrahydrobiopterin",
        # pyridoxal phosphate, the cofactor of aromatic L-amino acid decarboxylase and of glutamate decarboxylase
        "pyridoxal-phosphate", "pyridoxamine-phosphate",
        # ascorbate, the cofactor of dopamine beta-hydroxylase
        "ascorbate", "monodehydroascorbate", "dehydroascorbic acid",
        # one-carbon folate carriers, the methyl supply of S-adenosylmethionine and so of catechol O-methyltransferase
        "5,10-methylene-THF", "5-methyl-THF", "10-formyl-THF", "5,10-methenyl-THF", "5-formyl-THF",
        "5-formiminotetrahydrofolate", "dihydrofolate",
        # thioredoxins, the arm of the oxidative response the glutathione and NADPH pools do not cover
        "thioredoxin", "oxidized thioredoxin", "mitothioredoxin", "mitooxidized thioredoxin",
    )
)
# metabolites the experiment reads out, as a transmitter the graph releases or as a laboratory label; the
# recurring-pair clause makes glutamate a carrier in every transamination, which removes the depletion edge and the
# main-substrate coupling of a node whose response the experiment needs. See the "all_but_transmitters" rule.
TRANSMITTER_METABOLITE_NAMES: frozenset[str] = frozenset(
    name.lower() for name in ("glutamate", "GABA", "glycine", "aspartate", "L-aspartate", "serine", "L-serine")
)
CARRIER_RULES = ("all", "all_but_transmitters", "names", "neuronal")


def cofactor_edge_mask(
    edge_source: np.ndarray,
    edge_target: np.ndarray,
    edge_relation: np.ndarray,
    relation_types: list[str],
    node_base_metabolite_ids: np.ndarray,
    node_display_names: np.ndarray,
    is_currency: np.ndarray,
    minimum_pair_reactions: int = DEFAULT_MINIMUM_PAIR_REACTIONS,
    carrier_names: frozenset[str] = DEFAULT_CARRIER_METABOLITE_NAMES,
    carrier_rule: str = "all",
) -> np.ndarray:
    """Boolean mask over edges: True for substrate_of and product_of edges that are carrier edges (module docstring).

    carrier_rule is one of CARRIER_RULES and overrides carrier_names for "neuronal"; it is the only argument that
    changes which clause of the rule applies."""
    if carrier_rule not in CARRIER_RULES:
        raise ValueError(f"carrier_rule must be one of {CARRIER_RULES}, not {carrier_rule!r}")
    if carrier_rule == "neuronal":
        carrier_names = NEURONAL_CARRIER_METABOLITE_NAMES
    use_recurring_pair_clause = carrier_rule in ("all", "all_but_transmitters")
    pair_clause_exempt_names = TRANSMITTER_METABOLITE_NAMES if carrier_rule == "all_but_transmitters" else frozenset()
    mask = np.zeros(len(edge_source), dtype=bool)
    if "substrate_of" not in relation_types or "product_of" not in relation_types:
        return mask
    substrate_relation, product_relation = relation_types.index("substrate_of"), relation_types.index("product_of")
    names = np.array([str(name).lower() for name in node_display_names])
    on_carrier_list = np.isin(names, list(carrier_names))
    exempt_from_pair_clause = np.isin(names, list(pair_clause_exempt_names)) if pair_clause_exempt_names else np.zeros(len(names), dtype=bool)
    substrate_edges = np.where(edge_relation == substrate_relation)[0]  # metabolite -> reaction
    product_edges = np.where(edge_relation == product_relation)[0]  # reaction -> metabolite
    substrates_of_reaction: dict[int, list[int]] = defaultdict(list)
    products_of_reaction: dict[int, list[int]] = defaultdict(list)
    for edge in substrate_edges:
        if not is_currency[edge_source[edge]]:
            substrates_of_reaction[int(edge_target[edge])].append(int(edge_source[edge]))
    for edge in product_edges:
        if not is_currency[edge_target[edge]]:
            products_of_reaction[int(edge_source[edge])].append(int(edge_target[edge]))
    base = np.asarray(node_base_metabolite_ids, dtype=object)
    pair_reaction_counts: Counter = Counter()
    for reaction, substrates in substrates_of_reaction.items():
        for substrate in substrates:
            for product in products_of_reaction.get(reaction, []):
                if base[substrate] != base[product]:  # a transport reaction moves one species between compartments
                    pair_reaction_counts[(base[substrate], base[product])] += 1

    def in_recurring_pair(metabolite: int, reaction: int, metabolite_is_substrate: bool) -> bool:
        if not use_recurring_pair_clause or exempt_from_pair_clause[metabolite]:
            return False
        partners = products_of_reaction.get(reaction, []) if metabolite_is_substrate else substrates_of_reaction.get(reaction, [])
        for partner in partners:
            pair = (base[metabolite], base[partner]) if metabolite_is_substrate else (base[partner], base[metabolite])
            if pair_reaction_counts.get(pair, 0) >= minimum_pair_reactions:
                return True
        return False

    for edge in substrate_edges:
        metabolite, reaction = int(edge_source[edge]), int(edge_target[edge])
        mask[edge] = bool(on_carrier_list[metabolite]) or in_recurring_pair(metabolite, reaction, metabolite_is_substrate=True)
    for edge in product_edges:
        reaction, metabolite = int(edge_source[edge]), int(edge_target[edge])
        mask[edge] = bool(on_carrier_list[metabolite]) or in_recurring_pair(metabolite, reaction, metabolite_is_substrate=False)
    return mask
