"""experiments/audit_transporter_connectivity.py: a transporter counts as connected through its protein node in a split
graph, its own node in a merged one, or a Reactome entity it is a member of."""
import pandas as pd

from experiments.audit_transporter_connectivity import transport_connectivity


def test_a_transporter_is_connected_through_its_protein_its_entity_or_its_own_gene_node():
    nodes = pd.DataFrame({"node_id": ["GENE:A", "PROTEIN:A", "GENE:B", "PROTEIN:B", "RCTE_1", "GENE:C", "PROTEIN:C", "PROTEIN:X", "GENE:D", "T1", "T2", "T3", "T4",
                                      "glucose_c", "glucose_e", "glutamate_c", "glutamate_v"],
                          "is_transport": [False, False, False, False, False, False, False, False, False, True, True, True, True, False, False, False, False],
                          "display_name": ["A", "A", "B", "B", "entity", "C", "C", "X", "D", "T1", "T2", "T3", "T4", "glucose", "glucose", "glutamate", "glutamate"]})
    rows = [("GENE:A", "PROTEIN:A", "encodes"), ("PROTEIN:A", "T1", "catalyzed_by"), ("glucose_e", "T1", "substrate_of"), ("T1", "glucose_c", "product_of"),
            ("GENE:B", "PROTEIN:B", "encodes"), ("PROTEIN:B", "RCTE_1", "member_of"), ("RCTE_1", "T2", "catalyzed_by"), ("glutamate_c", "T2", "substrate_of"),
            ("T2", "glutamate_v", "product_of"),
            ("GENE:C", "PROTEIN:C", "encodes"), ("PROTEIN:C", "PROTEIN:X", "binds"),
            ("GENE:D", "T3", "catalyzed_by")]
    edges = pd.DataFrame(rows, columns=["source_id", "target_id", "relation_type"])
    connectivity, transport_reactions, orphans = transport_connectivity(nodes, edges)
    assert connectivity["A"]["transport_reactions"] == {"T1"} and connectivity["A"]["metabolites"] == {"glucose"}, "through its protein node"
    assert connectivity["B"]["transport_reactions"] == {"T2"}, "through the Reactome entity its protein is a member of"
    assert not connectivity["C"]["transport_reactions"] and connectivity["C"]["other_edges"] == 1, "a binding edge only"
    assert connectivity["D"]["transport_reactions"] == {"T3"}, "a merged graph's gene node catalyses for its product"
    assert transport_reactions == 4 and orphans == 1, "T4 has no catalysing gene"
