import pandas as pd

importance = pd.read_csv(
    "results/subtype_classifier_gene_importance.csv"
)

annot = pd.read_csv(
    "data/raw/Human.GRCh38.p13.annot.tsv",
    sep="\t",
    low_memory=False
)

importance["GeneID"] = importance["GeneID"].astype(str)
annot["GeneID"] = annot["GeneID"].astype(str)

decoded = importance.merge(
    annot[["GeneID", "Symbol", "Description"]],
    on="GeneID",
    how="left"
)

print(decoded.head(30)[
    ["GeneID", "Symbol", "importance"]
])

decoded.to_csv(
    "results/subtype_classifier_gene_importance_named.csv",
    index=False
)