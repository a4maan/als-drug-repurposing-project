import glob
import os
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# -----------------------------
# Load all patient profiles
# -----------------------------

profile_files = glob.glob("profiles/patient_*_profile.csv")

rows = []

for file in profile_files:
    patient_id = os.path.basename(file).split("_")[2]

    df = pd.read_csv(file)

    for _, row in df.iterrows():
        rows.append({
            "patient": patient_id,
            "gene": row["Symbol"],
            "score": row["score"]
        })

long_df = pd.DataFrame(rows)

# -----------------------------
# Patient x Gene matrix
# -----------------------------

matrix = long_df.pivot_table(
    index="patient",
    columns="gene",
    values="score",
    fill_value=0
)

print("Patient-gene matrix shape:", matrix.shape)

# -----------------------------
# Scale data
# -----------------------------

X = StandardScaler().fit_transform(matrix)

# -----------------------------
# PCA
# -----------------------------

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X)

pca_df = pd.DataFrame({
    "patient": matrix.index,
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1]
})

# -----------------------------
# KMeans clustering
# -----------------------------

kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
pca_df["cluster"] = kmeans.fit_predict(X)

# -----------------------------
# Save results
# -----------------------------

pca_df.to_csv("results/phase2_patient_clusters.csv", index=False)

# -----------------------------
# Plot
# -----------------------------

plt.figure(figsize=(8, 6))

for cluster in sorted(pca_df["cluster"].unique()):
    subset = pca_df[pca_df["cluster"] == cluster]
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=f"Cluster {cluster}"
    )

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("ALS Patient Molecular Subtypes")
plt.legend()
plt.tight_layout()

plt.savefig("results/phase2_pca_clusters.png", dpi=300)

print("Saved:")
print("results/phase2_patient_clusters.csv")
print("results/phase2_pca_clusters.png")