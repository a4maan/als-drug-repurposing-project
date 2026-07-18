import pandas as pd

df = pd.read_csv(
    "results/phase2_patient_clusters.csv"
)

print(df["cluster"].value_counts())