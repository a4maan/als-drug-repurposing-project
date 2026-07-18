import pandas as pd

# -------------------------
# Load Cluster 1 signature
# -------------------------

cluster1 = pd.read_csv(
    "results/cluster1_signature_genes.csv"
)

cluster1_genes = set(
    cluster1["Gene"].dropna().astype(str)
)

print("Cluster 1 genes:", len(cluster1_genes))

# -------------------------
# Load DGIdb
# -------------------------

dgidb = pd.read_csv(
    "data/raw/dgidb/interactions.tsv",
    sep="\t",
    low_memory=False
)

print("DGIdb interactions:", len(dgidb))

# -------------------------
# Keep interactions involving
# Cluster 1 genes
# -------------------------

hits = dgidb[
    dgidb["gene_name"].isin(cluster1_genes)
].copy()

print("Matching interactions:", len(hits))

# -------------------------
# Count how many subtype genes
# each drug hits
# -------------------------

drug_scores = (
    hits.groupby("drug_name")["gene_name"]
    .nunique()
    .reset_index()
)

drug_scores.columns = [
    "drug_name",
    "genes_targeted"
]

drug_scores = drug_scores.sort_values(
    "genes_targeted",
    ascending=False
)

drug_scores.to_csv(
    "results/cluster1_drug_rankings.csv",
    index=False
)

print("\nTop candidate drugs:")
print(drug_scores.head(30))