import pandas as pd

drivers = pd.read_csv("results/cluster_driver_genes.csv")

# Positive Difference means higher in Cluster 0
cluster0 = drivers[drivers["Difference"] > 0].copy()

# Negative Difference means higher in Cluster 1
cluster1 = drivers[drivers["Difference"] < 0].copy()

cluster0 = cluster0.sort_values("Difference", ascending=False)
cluster1 = cluster1.sort_values("Difference", ascending=True)

cluster0_top = cluster0.head(100)
cluster1_top = cluster1.head(100)

cluster0_top.to_csv("results/cluster0_signature_genes.csv", index=False)
cluster1_top.to_csv("results/cluster1_signature_genes.csv", index=False)

print("Cluster 0 signature genes:")
print(cluster0_top[["Gene", "Difference"]].head(20))

print("\nCluster 1 signature genes:")
print(cluster1_top[["Gene", "Difference"]].head(20))

print("\nSaved:")
print("results/cluster0_signature_genes.csv")
print("results/cluster1_signature_genes.csv")