"""Cofactor (carrier) edges of the metabolic layer, for the linear-response encoder (design section 5.2).

In the stoichiometric graph a carrier such as tetrahydrobiopterin is an ordinary substrate of every reaction that uses
it, so a linear response couples reactions through it: a PAH loss of function leaves tetrahydrobiopterin unconsumed,
raises it and through it raises tyrosine hydroxylase flux, L-dopa and dopamine, the opposite of what phenylketonuria
shows. Currency metabolites (tag_currency_metabolites.py) already stop passing the response on; carriers are not
currency (they are specific and few reactions use them), so their edges get relations of their own instead
(cosubstrate_of, coproduct_of and the derived depletes_cosubstrate), whose gains the model learns apart from those of
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
        "ascorbate", "monodehydroascorbate",
        "lipoamide", "dihydrolipoamide",
        "pyridoxal-phosphate", "pyridoxamine-phosphate",
        "5,10-methylene-THF", "5-methyl-THF", "10-formyl-THF", "5,10-methenyl-THF", "5-formyl-THF",
        "5-formiminotetrahydrofolate", "dihydrofolate",
    )
)
DEFAULT_MINIMUM_PAIR_REACTIONS = 10


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
) -> np.ndarray:
    """Boolean mask over edges: True for substrate_of and product_of edges that are carrier edges (module docstring)."""
    mask = np.zeros(len(edge_source), dtype=bool)
    if "substrate_of" not in relation_types or "product_of" not in relation_types:
        return mask
    substrate_relation, product_relation = relation_types.index("substrate_of"), relation_types.index("product_of")
    names = np.array([str(name).lower() for name in node_display_names])
    on_carrier_list = np.isin(names, list(carrier_names))
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
