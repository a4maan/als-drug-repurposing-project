import pandas as pd
import numpy as np

EXPRESSION_FILE = "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
ANNOTATION_FILE = "data/raw/Human.GRCh38.p13.annot.tsv"

# Load files
expr = pd.read_csv(EXPRESSION_FILE, sep="\t")
annot = pd.read_csv(ANNOTATION_FILE, sep="\t")

# Add gene symbols
expr = expr.merge(
    annot[["GeneID", "Symbol", "Description", "GeneType"]],
    on="GeneID",
    how="left"
)

# Sample columns are GSM columns
sample_cols = [c for c in expr.columns if c.startswith("GSM")]

als_samples = sample_cols[:96]
control_samples = sample_cols[96:]

print("Total samples:", len(sample_cols))
print("ALS samples:", len(als_samples))
print("Control samples:", len(control_samples))

# Build control baseline
control_mean = expr[control_samples].mean(axis=1)

# Pick one ALS patient
patient = als_samples[0]
patient_value = expr[patient]

# Abnormality score
score = np.log2((patient_value + 1) / (control_mean + 1))

profile = pd.DataFrame({
    "GeneID": expr["GeneID"],
    "Symbol": expr["Symbol"],
    "Description": expr["Description"],
    "GeneType": expr["GeneType"],
    "score": score,
    "abs_score": score.abs()
})

# Remove genes with missing symbols
profile = profile.dropna(subset=["Symbol"])

# Optional: focus on protein-coding genes
profile = profile[profile["GeneType"] == "protein-coding"]

# Sort by most abnormal
profile = profile.sort_values("abs_score", ascending=False)

# Save top 200
top200 = profile.head(200)

top200.to_csv("profiles/patient_001_profile.csv", index=False)

print("\nTop 20 abnormal protein-coding genes for", patient)
print(top200[["Symbol", "GeneID", "score", "Description"]].head(20))

print("\nSaved:")
print("profiles/patient_001_profile.csv")