import pandas as pd
import numpy as np

EXPRESSION_FILE = "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
ANNOTATION_FILE = "data/raw/Human.GRCh38.p13.annot.tsv"

# Load files
expr = pd.read_csv(EXPRESSION_FILE, sep="\t")
annot = pd.read_csv(ANNOTATION_FILE, sep="\t", low_memory=False)

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
# Build profiles for all ALS patients

for i, patient in enumerate(als_samples, start=1):

    patient_value = expr[patient]

    score = np.log2(
        (patient_value + 1)
        /
        (control_mean + 1)
    )

    profile = pd.DataFrame({
        "GeneID": expr["GeneID"],
        "Symbol": expr["Symbol"],
        "Description": expr["Description"],
        "GeneType": expr["GeneType"],
        "score": score,
        "abs_score": score.abs()
    })

    profile = profile.dropna(subset=["Symbol"])
    profile = profile[profile["GeneType"] == "protein-coding"]

    profile = profile.sort_values(
        "abs_score",
        ascending=False
    )

    filename = (
        f"profiles/patient_{i:03d}_{patient}_profile.csv"
    )

    profile.head(200).to_csv(
        filename,
        index=False
    )

    print(f"Saved {filename}")

print("\nDone.")
print(f"Created {len(als_samples)} patient profiles.")