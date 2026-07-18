import pandas as pd

# Load cluster assignments
clusters = pd.read_csv(
    "results/phase2_patient_clusters.csv"
)

# Load expression matrix
expr = pd.read_csv(
    "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",
    sep="\t"
)

# Load annotations
annot = pd.read_csv(
    "data/raw/Human.GRCh38.p13.annot.tsv",
    sep="\t",
    low_memory=False
)

expr = expr.merge(
    annot[["GeneID", "Symbol", "GeneType"]],
    on="GeneID",
    how="left"
)

expr = expr[
    expr["GeneType"] == "protein-coding"
]


# Build sample lists
cluster0 = clusters[clusters["cluster"] == 0]["patient"].tolist()
cluster1 = clusters[clusters["cluster"] == 1]["patient"].tolist()

# Only keep samples that actually exist in expression matrix
cluster0 = [s for s in cluster0 if s in expr.columns]
cluster1 = [s for s in cluster1 if s in expr.columns]

print("Cluster 0 samples found:", len(cluster0))
print("Cluster 1 samples found:", len(cluster1))

missing0 = [s for s in clusters[clusters["cluster"] == 0]["patient"].tolist() if s not in expr.columns]
missing1 = [s for s in clusters[clusters["cluster"] == 1]["patient"].tolist() if s not in expr.columns]

print("Missing from cluster 0:", missing0)
print("Missing from cluster 1:", missing1)

results = []

for _, row in expr.iterrows():

    gene = row["Symbol"]

    if pd.isna(gene):
        continue

    mean0 = row[cluster0].mean()
    mean1 = row[cluster1].mean()

    fold_change = mean0 - mean1

    results.append(
        [gene, fold_change, mean0, mean1]
    )

results = pd.DataFrame(
    results,
    columns=[
        "Gene",
        "Difference",
        "Cluster0Mean",
        "Cluster1Mean"
    ]
)

results["AbsDiff"] = results[
    "Difference"
].abs()

results = results.sort_values(
    "AbsDiff",
    ascending=False
)

results.to_csv(
    "results/cluster_driver_genes.csv",
    index=False
)

print(
    results.head(30)
)