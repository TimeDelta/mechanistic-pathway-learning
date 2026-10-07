"""How much the encoder's two-stage mean cancels at each destination, at equal gains.

The encoder averages within each relation feeding a node, then across those relations. With all
gains equal, the net message at a destination is the mean over relations of the mean sign within
the relation. The cancellation ratio |net| / (mean over relations of mean |sign|) is 1 when no
cancellation occurs and 0 when it is total. Equal gains is the honest caveat: learned gains can
reweight relations, but they are shared across all nodes, so a gain cannot rescue one
destination without moving every other destination fed by the same relation.
"""
import sys
import numpy as np
import pandas as pd

graph_directory = sys.argv[1]
edges = pd.read_parquet(f"{graph_directory}/edges.parquet")
substrate_edges = edges[edges["relation_type"] == "substrate_of"]
signed_edges = pd.concat([
    edges[["source_id", "target_id", "relation_type", "sign"]],
    pd.DataFrame({"source_id": substrate_edges["target_id"].values,
                  "target_id": substrate_edges["source_id"].values,
                  "relation_type": "depletes_substrate", "sign": -1.0}),
], ignore_index=True)

within_relation = signed_edges.groupby(["target_id", "relation_type"])["sign"].mean().rename("signed_mean")
within_relation_magnitude = signed_edges.groupby(["target_id", "relation_type"])["sign"].apply(
    lambda signs: signs.abs().mean()).rename("magnitude_mean")
per_destination = pd.concat([within_relation, within_relation_magnitude], axis=1).groupby("target_id").mean()

cancellation_ratio = per_destination["signed_mean"].abs() / per_destination["magnitude_mean"]
print(f"graph {graph_directory}")
print(f"  destinations {len(cancellation_ratio)}")
for quantile in (0.1, 0.25, 0.5, 0.75, 0.9):
    print(f"  cancellation ratio, quantile {quantile:4.2f}: {cancellation_ratio.quantile(quantile):.3f}")
for threshold in (0.99, 0.5, 0.1):
    below = int((cancellation_ratio < threshold).sum())
    print(f"  destinations with ratio below {threshold:4.2f}: {below} "
          f"({100 * below / len(cancellation_ratio):.1f}%)")
print(f"  destinations with total cancellation (ratio below 0.01): "
      f"{int((cancellation_ratio < 0.01).sum())}")
