"""
Phase 1: Build a per-patient expression profile (log2 fold-change vs the control
baseline) for every ALS patient.

Fixes applied vs the original version:

  * Flaw 1 (label leakage): ALS vs control samples are now assigned from the
    GEO metadata table (results/sample_metadata.csv), joined by GSM accession,
    instead of the positional guess `sample_cols[:96]` / `[96:]`. Because 12
    samples are absent from NCBI's FPKM matrix, the positional split silently
    mislabeled 10 control subjects as ALS patients.

  * Flaw 2 (technical-artifact subtypes): with EXCLUDE_CONFOUNDERS=True, globin,
    ribosomal and mitochondrial gene families are dropped AND each sample's
    library is re-normalized over the remaining genes BEFORE selecting each
    patient's top genes. Simply dropping the genes is not enough: globin alone
    is ~60-80% of the FPKM in this peripheral-blood data, so it compositionally
    distorts every other gene's FPKM. Re-normalizing to a common library size
    over the non-confounder genes removes that distortion, so the downstream
    subtypes reflect biology rather than blood contamination / RNA quality.

Set EXCLUDE_CONFOUNDERS via the environment variable EXCLUDE_CONFOUNDERS
(1/true to enable) so the pipeline can be run both ways for comparison.
"""

import os
import glob
import pandas as pd
import numpy as np

from confounders import is_confounder

EXPRESSION_FILE = "data/raw/GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv"
ANNOTATION_FILE = "data/raw/Human.GRCh38.p13.annot.tsv"
METADATA_FILE = "results/sample_metadata.csv"

EXCLUDE_CONFOUNDERS = os.environ.get("EXCLUDE_CONFOUNDERS", "0").lower() in (
    "1", "true", "yes"
)


def main():
    expr = pd.read_csv(EXPRESSION_FILE, sep="\t")
    annot = pd.read_csv(ANNOTATION_FILE, sep="\t", low_memory=False)

    expr = expr.merge(
        annot[["GeneID", "Symbol", "Description", "GeneType"]],
        on="GeneID",
        how="left",
    )

    # ---- Flaw 1 fix: labels from metadata, not column position ----
    meta = pd.read_csv(METADATA_FILE)
    group_of = dict(zip(meta["gsm"], meta["group"]))

    sample_cols = [c for c in expr.columns if c.startswith("GSM")]
    als_samples = [c for c in sample_cols if group_of.get(c) == "ALS"]
    control_samples = [c for c in sample_cols if group_of.get(c) == "control"]
    unknown = [c for c in sample_cols if c not in group_of]

    print("Samples in expression matrix:", len(sample_cols))
    print("ALS samples (from metadata):", len(als_samples))
    print("Control samples (from metadata):", len(control_samples))
    if unknown:
        print("WARNING: samples with no metadata label (skipped):", unknown)

    print("Exclude confounder gene families:", EXCLUDE_CONFOUNDERS)

    # ---- Flaw 2 fix: deplete confounder families, then RE-NORMALIZE ----
    # Dropping globin/ribosomal/mito genes is not enough on its own: because
    # they make up the majority of the library, they inflate the FPKM
    # denominator and distort every other gene. After removing them we rescale
    # each sample to a common library size (TPM-like) over the remaining genes.
    if EXCLUDE_CONFOUNDERS:
        conf_mask = expr["Symbol"].map(is_confounder).fillna(False)
        n_conf = int(conf_mask.sum())
        expr = expr[~conf_mask].copy()
        lib = expr[sample_cols].sum(axis=0)
        expr[sample_cols] = expr[sample_cols].div(lib, axis=1) * 1e6
        print(f"Dropped {n_conf} confounder genes and re-normalized "
              f"{len(sample_cols)} samples to 1e6 over the remaining genes.")

    # Control baseline built from TRUE controls only.
    control_mean = expr[control_samples].mean(axis=1)

    # Clean out stale profiles so mislabeled runs don't linger.
    for old in glob.glob("profiles/patient_*_profile.csv"):
        os.remove(old)

    for i, patient in enumerate(als_samples, start=1):
        patient_value = expr[patient]

        score = np.log2((patient_value + 1) / (control_mean + 1))

        profile = pd.DataFrame({
            "GeneID": expr["GeneID"],
            "Symbol": expr["Symbol"],
            "Description": expr["Description"],
            "GeneType": expr["GeneType"],
            "score": score,
            "abs_score": score.abs(),
        })

        profile = profile.dropna(subset=["Symbol"])
        profile = profile[profile["GeneType"] == "protein-coding"]

        # ---- Flaw 2 fix: drop technical-confounder gene families ----
        if EXCLUDE_CONFOUNDERS:
            profile = profile[~profile["Symbol"].map(is_confounder)]

        profile = profile.sort_values("abs_score", ascending=False)

        filename = f"profiles/patient_{i:03d}_{patient}_profile.csv"
        profile.head(200).to_csv(filename, index=False)

    print(f"\nDone. Created {len(als_samples)} patient profiles.")


if __name__ == "__main__":
    main()
