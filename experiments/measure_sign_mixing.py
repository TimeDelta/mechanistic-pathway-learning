"""Signs as the encoder actually builds them, not as edges.parquet stores them."""
import sys
import pandas as pd

graph_directory = sys.argv[1]
edges = pd.read_parquet(f"{graph_directory}/edges.parquet")

# The encoder appends depletes_substrate: the reverse of every substrate_of edge, sign -1.
substrate_edges = edges[edges["relation_type"] == "substrate_of"]
derived_depletion = pd.DataFrame({
    "source_id": substrate_edges["target_id"].values,
    "target_id": substrate_edges["source_id"].values,
    "relation_type": "depletes_substrate",
    "sign": -1.0,
})
signed_edges = pd.concat([edges[["source_id", "target_id", "relation_type", "sign"]], derived_depletion],
                         ignore_index=True)

print(f"graph {graph_directory}")
print(f"  stored edges {len(edges)}, with the derived depletion relation {len(signed_edges)}")
print("  sign values per relation:")
for relation, group in signed_edges.groupby("relation_type"):
    print(f"    {relation:24s} {len(group):7d} edges, distinct signs {sorted(group['sign'].unique())}")

destination_signs = signed_edges.groupby("target_id")["sign"].agg(["min", "max"])
mixed_sign_destinations = int((destination_signs["min"] < 0).sum() & 1 if False else
                              ((destination_signs["min"] < 0) & (destination_signs["max"] > 0)).sum())
all_nodes = pd.unique(pd.concat([signed_edges["source_id"], signed_edges["target_id"]]))
print(f"  destinations receiving any message: {len(destination_signs)}")
print(f"  destinations receiving BOTH signs:  {mixed_sign_destinations} "
      f"({100 * mixed_sign_destinations / max(len(destination_signs), 1):.1f}% of destinations, "
      f"{100 * mixed_sign_destinations / len(all_nodes):.1f}% of {len(all_nodes)} nodes)")

within_relation_sign_variation = sum(
    int(group.groupby("target_id")["sign"].nunique().gt(1).sum())
    for _, group in signed_edges.groupby("relation_type"))
print(f"  destinations seeing both signs WITHIN one relation: {within_relation_sign_variation}")
