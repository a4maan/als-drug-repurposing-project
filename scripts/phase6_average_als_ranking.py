import pandas as pd
import numpy as np

# -------------------------
# Load expression matrix
# -------------------------

expr = pd.read_csv(
    "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",
    sep="\t"
)

# -------------------------
# Metadata
# -------------------------

samples = expr.columns[1:]

als_samples = samples[:96]
control_samples = samples[96:]

print("ALS:", len(als_samples))
print("Controls:", len(control_samples))

# -------------------------
# Average ALS signature
# -------------------------

expr["als_mean"] = expr[als_samples].mean(axis=1)
expr["control_mean"] = expr[control_samples].mean(axis=1)

expr["score"] = np.log2(
    (expr["als_mean"] + 1)
    /
    (expr["control_mean"] + 1)
)

# -------------------------
# Gene ID -> Symbol
# -------------------------

annot = pd.read_csv(
    "data/raw/Human.GRCh38.p13.annot.tsv",
    sep="\t",
    low_memory=False
)

expr = expr.merge(
    annot[["GeneID", "Symbol"]],
    on="GeneID",
    how="left"
)

expr = expr.dropna(subset=["Symbol"])

# -------------------------
# Top ALS genes
# -------------------------

expr["importance"] = expr["score"].abs()

top_genes = (
    expr.sort_values(
        "importance",
        ascending=False
    )
    .head(100)
)

top_genes.to_csv(
    "results/average_als_signature.csv",
    index=False
)

print(top_genes[
    ["Symbol", "score"]
].head(20))