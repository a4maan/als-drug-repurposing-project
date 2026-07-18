import pandas as pd

# STRING info
info = pd.read_csv(
    "data/raw/string/9606.protein.info.v12.0.txt",
    sep="\t"
)

# Hub proteins
hubs = pd.read_csv(
    "results/cluster1_hub_genes.csv"
)

# Merge
decoded = hubs.merge(
    info[
        ["#string_protein_id", "preferred_name"]
    ],
    left_on="protein",
    right_on="#string_protein_id",
    how="left"
)

decoded = decoded[
    ["preferred_name", "centrality"]
]

decoded = decoded.rename(
    columns={
        "preferred_name": "Gene"
    }
)

print(decoded.head(30))

decoded.to_csv(
    "results/cluster1_hub_genes_named.csv",
    index=False
)