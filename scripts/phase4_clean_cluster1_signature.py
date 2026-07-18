import pandas as pd

sig = pd.read_csv(
    "results/cluster1_signature_genes.csv"
)

REMOVE_PREFIXES = [
    "RPL",
    "RPS",
    "EEF",
    "MT-"
]

REMOVE_EXACT = {
    "COX1",
    "COX2",
    "COX3",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6",
    "ATP6",
    "ATP8",
    "CYTB"
}

def keep_gene(gene):

    gene = str(gene)

    for prefix in REMOVE_PREFIXES:
        if gene.startswith(prefix):
            return False

    if gene in REMOVE_EXACT:
        return False

    return True

clean = sig[
    sig["Gene"].apply(keep_gene)
].copy()

clean.to_csv(
    "results/cluster1_signature_clean.csv",
    index=False
)

print("Original genes:", len(sig))
print("Remaining genes:", len(clean))

print("\nTop 30 cleaned genes:")
print(
    clean[
        ["Gene", "Difference"]
    ].head(30)
)