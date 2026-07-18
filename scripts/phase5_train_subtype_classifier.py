import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score

# -----------------------------
# Load expression matrix
# -----------------------------

expr = pd.read_csv(
    "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",
    sep="\t"
)

# -----------------------------
# Load cluster assignments
# -----------------------------

clusters = pd.read_csv(
    "results/phase2_patient_clusters.csv"
)

# -----------------------------
# Build feature matrix
# -----------------------------

sample_cols = expr.columns[1:]

# Transpose so rows = patients, columns = genes
X = expr[sample_cols].T

# Gene IDs become feature names
X.columns = expr["GeneID"].astype(str)

# Patient IDs become a column
X.index.name = "patient"
X.reset_index(inplace=True)

# -----------------------------
# Merge expression data with subtype labels
# -----------------------------

data = X.merge(
    clusters[["patient", "cluster"]],
    on="patient"
)

print("Patients:", len(data))

# -----------------------------
# Features and labels
# -----------------------------

y = data["cluster"]

X = data.drop(
    columns=["patient", "cluster"]
)

# Fix sklearn column-name issue
X.columns = X.columns.astype(str)

# -----------------------------
# Train/test split
# -----------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.25,
    random_state=42,
    stratify=y
)

# -----------------------------
# Train Random Forest classifier
# -----------------------------

model = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced"
)

model.fit(X_train, y_train)

# -----------------------------
# Evaluate
# -----------------------------

pred = model.predict(X_test)

print("\nAccuracy:")
print(accuracy_score(y_test, pred))

print("\nClassification Report:")
print(classification_report(y_test, pred))

# -----------------------------
# Feature importance
# -----------------------------

importance = pd.DataFrame({
    "GeneID": X.columns,
    "importance": model.feature_importances_
})

importance = importance.sort_values(
    "importance",
    ascending=False
)

importance.to_csv(
    "results/subtype_classifier_gene_importance.csv",
    index=False
)

print("\nTop 30 predictive genes:")
print(importance.head(30))

print("\nSaved:")
print("results/subtype_classifier_gene_importance.csv")