"""
Shared definitions of technical-confounder gene families and helpers.

Flaw 2 in the original pipeline: the "molecular subtypes" were dominated by
genes that track sample QUALITY and CELL COMPOSITION rather than ALS biology:

  * Hemoglobin / globin genes  -> red-blood-cell / globin contamination of
    peripheral-blood RNA (a well known confounder in blood RNA-seq).
  * Ribosomal protein genes    -> RNA integrity / degradation and library-prep
    effects; they are also universal STRING hubs, which is why they trivially
    dominated the network centrality step.
  * Mitochondrially-encoded genes -> mitochondrial content / RNA quality.

These modules provide the gene sets and a routine to compute each sample's
fractional expression in every family, so we can (a) test whether clusters are
driven by them and (b) optionally deplete them before subtype discovery.
"""

import re
import pandas as pd

# Hemoglobin / globin subunits (alpha, beta, delta, epsilon, gamma, mu, theta, zeta)
GLOBIN_GENES = {
    "HBA1", "HBA2", "HBB", "HBD", "HBE1",
    "HBG1", "HBG2", "HBM", "HBQ1", "HBZ",
}

# Mitochondrially-encoded genes. This NCBI annotation uses LEGACY symbols
# (COX1, ND1, ...), not the HGNC "MT-" prefix.
MITO_GENES = {
    "COX1", "COX2", "COX3",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "ATP6", "ATP8", "CYTB",
}


def is_ribosomal(symbol):
    """Cytosolic (RPL/RPS, RPLP*, RPSA) and mitochondrial (MRPL/MRPS) ribosomal proteins."""
    s = str(symbol)
    return bool(
        re.match(r"^RP[LS]\d", s)
        or re.match(r"^MRP[LS]\d", s)
        or s in {"RPLP0", "RPLP1", "RPLP2", "RPSA"}
    )


def is_confounder(symbol):
    """True if the gene belongs to any technical-confounder family."""
    s = str(symbol)
    return s in GLOBIN_GENES or s in MITO_GENES or is_ribosomal(s)


def family_of(symbol):
    s = str(symbol)
    if s in GLOBIN_GENES:
        return "globin"
    if s in MITO_GENES:
        return "mitochondrial"
    if is_ribosomal(s):
        return "ribosomal"
    return "other"


def sample_confounder_fractions(expr, sample_cols, symbol_col="Symbol"):
    """
    For each sample, compute the fraction of total FPKM contributed by each
    confounder family. Returns a DataFrame indexed by sample (gsm) with columns
    globin_frac, ribosomal_frac, mitochondrial_frac.

    `expr` must contain `symbol_col` plus the sample FPKM columns.
    """
    fam = expr[symbol_col].map(family_of)
    totals = expr[sample_cols].sum(axis=0)  # per-sample total FPKM
    out = {}
    for family in ["globin", "ribosomal", "mitochondrial"]:
        mask = (fam == family).values
        fam_sum = expr.loc[mask, sample_cols].sum(axis=0)
        out[f"{family}_frac"] = fam_sum / totals
    frac = pd.DataFrame(out)
    frac.index.name = "gsm"
    return frac.reset_index()
