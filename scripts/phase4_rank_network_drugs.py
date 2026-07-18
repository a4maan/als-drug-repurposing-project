import pandas as pd

# -----------------------------
# Load cleaned inflammatory ALS module
# -----------------------------

module = pd.read_csv("results/cluster1_signature_clean.csv")

# Use top 50 cleaned genes for now
module = module.head(50).copy()

# Convert Difference to importance score
# More negative = higher in Cluster 1
module["importance"] = module["Difference"].abs()

module_scores = dict(
    zip(module["Gene"], module["importance"])
)

module_genes = set(module_scores.keys())

print("Module genes:", len(module_genes))

# -----------------------------
# Load DGIdb
# -----------------------------

dgidb = pd.read_csv(
    "data/raw/dgidb/interactions.tsv",
    sep="\t",
    low_memory=False
)

# Keep only interactions with module genes
hits = dgidb[
    dgidb["gene_name"].isin(module_genes)
].copy()

print("Drug-gene hits:", len(hits))

# -----------------------------
# Score drugs
# -----------------------------

drug_results = []

for drug, group in hits.groupby("drug_name"):

    targeted_genes = sorted(set(group["gene_name"]))

    score = sum(
        module_scores.get(gene, 0)
        for gene in targeted_genes
    )

    drug_results.append({
        "drug_name": drug,
        "network_score": score,
        "genes_targeted": len(targeted_genes),
        "targeted_genes": ",".join(targeted_genes)
    })

ranked = pd.DataFrame(drug_results)

ranked = ranked.sort_values(
    ["network_score", "genes_targeted"],
    ascending=False
)

ranked.to_csv(
    "results/cluster1_network_drug_rankings.csv",
    index=False
)

print("\nTop 30 network-aware drug candidates:")
print(ranked.head(30))