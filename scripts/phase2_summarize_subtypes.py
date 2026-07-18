import pandas as pd

drivers = pd.read_csv("results/cluster_driver_genes.csv")

def classify_gene(gene):
    gene = str(gene)

    if gene.startswith("HB") or gene in ["HBA1", "HBA2", "HBB", "HBD", "HBQ1"]:
        return "Hemoglobin / blood"

    if gene.startswith("HLA") or gene in [
        "S100A8", "S100A9", "LYZ", "B2M", "IFITM2", "IFITM3",
        "ISG15", "IFI6", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"
    ]:
        return "Immune / inflammatory"

    if gene in [
        "COX1", "COX2", "COX3", "ATP6", "ATP8", "ND1", "ND2",
        "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB"
    ]:
        return "Mitochondrial"

    if gene.startswith("RPL") or gene.startswith("RPS"):
        return "Ribosomal"

    return "Other"

drivers["Category"] = drivers["Gene"].apply(classify_gene)

top100 = drivers.head(100)

summary = top100["Category"].value_counts()

print("\nTop 100 driver gene categories:")
print(summary)

print("\nTop genes by category:")
for category in summary.index:
    subset = top100[top100["Category"] == category]
    print("\n" + category)
    print(subset[["Gene", "Difference", "Cluster0Mean", "Cluster1Mean"]].head(10))

summary.to_csv("results/subtype_driver_category_counts.csv")
top100.to_csv("results/top100_cluster_driver_genes_categorized.csv", index=False)

print("\nSaved:")
print("results/subtype_driver_category_counts.csv")
print("results/top100_cluster_driver_genes_categorized.csv")