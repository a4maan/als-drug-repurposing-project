"""
Phase 0: Parse ground-truth sample metadata from the GEO series matrix.

This fixes Flaw 1 (labels were assigned positionally with
`sample_cols[:96]` / `[96:]`). Instead of assuming column order, we read the
actual per-sample labels GEO ships in the series matrix and write a clean
metadata table that every downstream phase can join against by GSM accession.

Output: results/sample_metadata.csv with columns
    gsm, title, group ("ALS"/"control"), disease_state, tissue
"""

import re
import pandas as pd

SERIES_MATRIX = "data/raw/GSE234297_series_matrix.txt"
OUT = "results/sample_metadata.csv"


def parse_row(line):
    """Split a '!Key\t"v1"\t"v2"...' series-matrix row into a list of values."""
    parts = line.rstrip("\n").split("\t")
    # parts[0] is the key; the rest are quoted values
    return [p.strip().strip('"') for p in parts[1:]]


def main():
    titles = accessions = None
    disease = tissue = None

    with open(SERIES_MATRIX) as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = parse_row(line)
            elif line.startswith("!Sample_geo_accession"):
                accessions = parse_row(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = parse_row(line)
                # There are multiple characteristic rows; route by field name.
                field = vals[0].split(":", 1)[0].strip().lower() if vals else ""
                if field == "disease state":
                    disease = [v.split(":", 1)[1].strip() for v in vals]
                elif field == "tissue":
                    tissue = [v.split(":", 1)[1].strip() for v in vals]

    assert accessions and titles, "Could not parse series matrix"

    n = len(accessions)
    disease = disease or [""] * n
    tissue = tissue or [""] * n

    def group_of(title, dstate):
        # Ground truth comes from the sample title (Case/Control) and is
        # corroborated by the disease-state characteristic (sALS vs control).
        t = title.lower()
        if "control" in t or "ctrl" in t:
            return "control"
        if "case" in t:
            return "ALS"
        # Fall back to disease state if the title is unexpected.
        return "ALS" if "als" in dstate.lower() else "control"

    meta = pd.DataFrame({
        "gsm": accessions,
        "title": titles,
        "disease_state": disease,
        "tissue": tissue,
    })
    meta["group"] = [
        group_of(t, d) for t, d in zip(meta["title"], meta["disease_state"])
    ]

    meta.to_csv(OUT, index=False)

    print(f"Parsed {len(meta)} samples from {SERIES_MATRIX}")
    print(meta["group"].value_counts())
    print("\nCross-check group vs disease_state:")
    print(pd.crosstab(meta["group"], meta["disease_state"]))
    print(f"\nSaved {OUT}")


if __name__ == "__main__":
    main()
