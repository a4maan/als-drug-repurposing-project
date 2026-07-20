# Methodological fixes: Flaw 1 (label leakage) and Flaw 2 (technical-artifact subtypes)

This document records two validity-threatening flaws in the original pipeline,
the fixes applied, and the test evidence produced by re-running the pipeline on
the real GSE234297 data.

Reproduce the inputs with `python3 scripts/download_data.py` (fetches the GEO
series matrix, the NCBI FPKM matrix, and the gene annotation into `data/raw/`).

---

## Flaw 1 — ALS/control labels were assigned by column position

**Problem.** `phase1_build_profiles.py` and `phase6_average_als_ranking.py`
split samples with `sample_cols[:96]` = ALS and `sample_cols[96:]` = control.
This assumes the matrix contains all 96 ALS samples first, in order.

**Reality.** The GEO series matrix defines **96 sALS + 48 healthy controls**
(144 samples). NCBI's FPKM matrix, however, contains only **132** of them —
**12 samples are missing (10 ALS, 2 control)**. The positional split therefore
put **10 control subjects into the "ALS" group** and treated them as patients
throughout profiling, clustering, and signature/drug ranking.

| Assumed (positional) | Truth (metadata) |
|---|---|
| first 96 columns = ALS | 86 ALS + **10 controls mislabeled** |
| remaining 36 = control | 36 controls |
| **Total mislabeled** | **10 / 132 samples** |

**Fix.**
- New `scripts/phase0_parse_metadata.py` parses the series matrix into
  `results/sample_metadata.csv` (`gsm, title, group, disease_state, tissue`).
  Title-derived group (Case/Control) agrees 100% with the disease-state field
  (sALS / Healthy control).
- `phase1` and `phase6` now select ALS vs control **by joining on GSM
  accession**, and build the control baseline from **true controls only**.
  Result: 86 ALS / 46 control, correctly labeled.

---

## Flaw 2 — the "molecular subtypes" were a technical/QC artifact

**Problem.** The original driver genes were dominated by ribosomal,
mitochondrial and hemoglobin genes — families that track **RNA quality and
blood contamination**, not ALS biology.

**Reality (this is peripheral blood).** Globin genes alone are **63–84 % of
total FPKM** per sample. `scripts/phase2b_confounder_analysis.py` computes each
sample's globin/ribosomal/mitochondrial fraction and tests it against the
cluster label.

With **labels fixed but confounders still present** (`EXCLUDE_CONFOUNDERS=0`),
the 66/20 subtype split is strongly confounded:

| family | cluster0 mean | cluster1 mean | Cohen's d | Mann–Whitney p |
|---|---|---|---|---|
| globin | 0.785 | 0.626 | **1.54** | 5.9e-07 |
| ribosomal | 0.017 | 0.028 | -1.12 | 1.1e-05 |
| mitochondrial | 0.020 | 0.033 | -1.03 | 6.6e-05 |

**3 / 3 families differ significantly** — the split is essentially a QC axis.
(Its silhouette is only 0.030: weak cluster structure to begin with.)

**Fix — and why dropping the genes is not enough.** Simply removing the
confounder genes (what the old `phase4_clean` script did *after* the fact) still
left all 3/3 families significant (globin d=1.44, p=2e-9), because globin
dominating the library **compositionally distorts every other gene's FPKM**.
The correct fix, now in `phase1` under `EXCLUDE_CONFOUNDERS=1`, is to **drop the
confounder families AND re-normalize each sample's library over the remaining
genes** before computing fold changes.

**Result of the corrected run (`EXCLUDE_CONFOUNDERS=1`).** The technical
association disappears — and so does the balanced two-way split: KMeans
collapses to **85 / 1**, and silhouette across k = 2..6 only ever peels off
single outliers (0.13, 0.12, 0.04, …). In other words:

> The two "molecular subtypes" reported by this pipeline **do not survive
> confounder correction**. They were driven by blood-contamination / RNA-quality
> differences between samples, not by ALS biology.

---

## Consequence for downstream results

Because the subtype split was artifactual, the subtype-specific outputs that
depend on it (`cluster*_signature_genes`, `cluster1_*_drug_rankings`,
`cluster1_hub_genes`, the Phase-5 classifier) are **not supported as biology**.
Recovering genuine ALS blood subtypes (which the source study reports) requires
the further fixes flagged in the methodology review — start from raw counts with
proper normalization/DE (DESeq2/limma-voom), choose k with silhouette/consensus
clustering, and validate on an independent cohort — not the FPKM + top-200 +
KMeans approach.

## How to run

```bash
python3 scripts/download_data.py            # fetch data/raw/ inputs
python3 scripts/phase0_parse_metadata.py    # -> results/sample_metadata.csv

# corrected, biology-focused run
EXCLUDE_CONFOUNDERS=1 python3 scripts/phase1_build_profiles.py
python3 scripts/phase2_cluster_patients.py
python3 scripts/phase2b_confounder_analysis.py confounders_removed

# comparison run that leaves the confounders in
EXCLUDE_CONFOUNDERS=0 python3 scripts/phase1_build_profiles.py
python3 scripts/phase2_cluster_patients.py
python3 scripts/phase2b_confounder_analysis.py labels_fixed
```
