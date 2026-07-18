import pandas as pd

INFO_FILE = "data/raw/string/9606.protein.info.v12.0.txt"
LINKS_FILE = "data/raw/string/9606.protein.links.v12.0.txt"

CLUSTER0_FILE = "results/cluster0_signature_genes.csv"
CLUSTER1_FILE = "results/cluster1_signature_genes.csv"

# Load STRING protein info
info = pd.read_csv(INFO_FILE, sep="\t")

print("STRING info columns:")
print(info.columns.tolist())
print(info.head())

# Load subtype signatures
cluster0 = pd.read_csv(CLUSTER0_FILE)
cluster1 = pd.read_csv(CLUSTER1_FILE)

cluster0_genes = set(cluster0["Gene"].dropna().astype(str))
cluster1_genes = set(cluster1["Gene"].dropna().astype(str))

# Map gene symbol -> STRING protein ID
# STRING usually uses preferred_name for gene symbols.
gene_to_string = dict(zip(info["preferred_name"], info["#string_protein_id"]))

cluster0_ids = {
    gene_to_string[g] for g in cluster0_genes
    if g in gene_to_string
}

cluster1_ids = {
    gene_to_string[g] for g in cluster1_genes
    if g in gene_to_string
}

print("\nCluster 0 genes:", len(cluster0_genes))
print("Cluster 0 genes mapped to STRING:", len(cluster0_ids))

print("\nCluster 1 genes:", len(cluster1_genes))
print("Cluster 1 genes mapped to STRING:", len(cluster1_ids))

# Load STRING links
# This file can be large, but your PC should handle it.
links = pd.read_csv(LINKS_FILE, sep=" ")

print("\nSTRING links columns:")
print(links.columns.tolist())
print("Total STRING interactions:", len(links))

# Keep only high-confidence interactions
links = links[links["combined_score"] >= 700]

print("High-confidence interactions:", len(links))

def build_subtype_network(ids, name):
    network = links[
        links["protein1"].isin(ids) &
        links["protein2"].isin(ids)
    ].copy()

    network.to_csv(
        f"results/{name}_string_network.csv",
        index=False
    )

    print(f"\n{name}:")
    print("Nodes in signature:", len(ids))
    print("Edges found:", len(network))
    print(f"Saved results/{name}_string_network.csv")

build_subtype_network(cluster0_ids, "cluster0")
build_subtype_network(cluster1_ids, "cluster1")