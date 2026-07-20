"""
Phase 2b: Confounder analysis for the discovered subtypes (Flaw 2).

The original subtypes were dominated by globin / ribosomal / mitochondrial
genes -- families that track blood contamination and RNA quality, not ALS
biology. This script quantifies that directly: for each ALS sample it computes
the fraction of total FPKM contributed by each confounder family, then tests
whether the two clusters differ on those fractions (Mann-Whitney U).

A large, significant difference means the "subtype" split is (at least partly)
a technical artifact. Running the pipeline with EXCLUDE_CONFOUNDERS=1 should
shrink these differences.

Outputs:
    results/phase2b_sample_confounder_fractions.csv
    results/phase2b_confounder_report.csv
"""

import sys
import pandas as pd
from scipy.stats import mannwhitneyu

from confounders import sample_confounder_fractions

EXPRESSION_FILE = "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
ANNOTATION_FILE = "data/raw/Human.GRCh38.p13.annot.tsv"
CLUSTERS_FILE = "results/phase2_patient_clusters.csv"

# Tag results so we can compare the confounder-included vs -excluded runs.
LABEL = sys.argv[1] if len(sys.argv) > 1 else "default"


def cohens_d(a, b):
    import numpy as np
    a, b = np.asarray(a), np.asarray(b)
    na, nb = len(a), len(b)
    pooled = ((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2)
    if pooled == 0:
        return 0.0
    return (a.mean() - b.mean()) / (pooled ** 0.5)


def main():
    expr = pd.read_csv(EXPRESSION_FILE, sep="\t")
    annot = pd.read_csv(ANNOTATION_FILE, sep="\t", low_memory=False)
    expr = expr.merge(annot[["GeneID", "Symbol"]], on="GeneID", how="left")

    clusters = pd.read_csv(CLUSTERS_FILE)  # columns: patient (gsm), cluster, ...
    als_samples = [c for c in clusters["patient"] if c in expr.columns]

    frac = sample_confounder_fractions(expr, als_samples)
    frac = frac.merge(
        clusters[["patient", "cluster"]],
        left_on="gsm", right_on="patient", how="left",
    ).drop(columns="patient")

    frac.to_csv("results/phase2b_sample_confounder_fractions.csv", index=False)

    print(f"[{LABEL}] ALS samples analysed:", len(frac))
    sizes = frac["cluster"].value_counts()
    print("Cluster sizes:\n", sizes.to_string())

    # A confounder test only makes sense with two non-trivial groups. If the
    # clustering collapsed (e.g. an 85/1 split after correction) say so: that
    # degeneracy is itself the finding -- no balanced subtype structure remains.
    degenerate = sizes.min() < 3 or len(sizes) < 2
    if degenerate:
        print("\nNOTE: clustering is degenerate (a cluster has < 3 members).")
        print("      No balanced two-subtype structure to test -- reporting means only.")

    rows = []
    for col in ["globin_frac", "ribosomal_frac", "mitochondrial_frac"]:
        g0 = frac[frac["cluster"] == 0][col].dropna()
        g1 = frac[frac["cluster"] == 1][col].dropna()
        if len(g0) < 3 or len(g1) < 3:
            u, p, d = float("nan"), float("nan"), float("nan")
        else:
            u, p = mannwhitneyu(g0, g1, alternative="two-sided")
            d = cohens_d(g0, g1)
        rows.append({
            "run": LABEL,
            "family": col,
            "cluster0_mean": g0.mean(),
            "cluster1_mean": g1.mean(),
            "fold": (g0.mean() / g1.mean()) if g1.mean() else float("nan"),
            "cohens_d": d,
            "mannwhitney_p": p,
            "significant_0.05": bool(p < 0.05),  # NaN -> False
        })

    report = pd.DataFrame(rows)
    pd.set_option("display.float_format", lambda x: f"{x:.4g}")
    print("\n=== Confounder association with cluster label ===")
    print(report.to_string(index=False))

    n_sig = int(report["significant_0.05"].sum())
    if degenerate:
        print(f"\nVERDICT [{LABEL}]: clustering collapsed to a degenerate split; "
              "the balanced two-subtype structure did NOT survive correction.")
    else:
        print(f"\nVERDICT [{LABEL}]: {n_sig}/3 confounder families differ "
              f"significantly (p<0.05) between the two subtypes.")
        if n_sig:
            print("  -> The subtype split is associated with technical/QC signal.")
        else:
            print("  -> No significant technical association detected.")

    report.to_csv(f"results/phase2b_confounder_report_{LABEL}.csv", index=False)
    print(f"\nSaved results/phase2b_confounder_report_{LABEL}.csv")


if __name__ == "__main__":
    main()
