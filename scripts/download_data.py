"""
Download the GEO inputs needed by phases 0-2 and 6 into data/raw/.

Files fetched:
  * GSE234297_series_matrix.txt              (sample metadata: case/control)
  * GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv   (expression matrix)
  * Human.GRCh38.p13.annot.tsv               (GeneID -> Symbol/GeneType)

STRING and DGIdb inputs (used by phases 3-4) are not fetched here; download
them separately into data/raw/string/ and data/raw/dgidb/.

Run from the repository root:  python3 scripts/download_data.py
"""

import gzip
import os
import shutil
import urllib.request

RAW = "data/raw"

FILES = {
    "GSE234297_series_matrix.txt":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234297/matrix/"
        "GSE234297_series_matrix.txt.gz",
    "GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv":
        "https://www.ncbi.nlm.nih.gov/geo/download/"
        "?type=rnaseq_counts&acc=GSE234297&format=file"
        "&file=GSE234297_norm_counts_FPKM_GRCh38.p13_NCBI.tsv.gz",
    "Human.GRCh38.p13.annot.tsv":
        "https://www.ncbi.nlm.nih.gov/geo/download/"
        "?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz",
}


def fetch(dest, url):
    if os.path.exists(dest):
        print(f"exists, skipping: {dest}")
        return
    gz = dest + ".gz"
    print(f"downloading {url}")
    with urllib.request.urlopen(url) as r, open(gz, "wb") as fh:
        shutil.copyfileobj(r, fh)
    with gzip.open(gz, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    os.remove(gz)
    print(f"  -> {dest}")


def main():
    os.makedirs(RAW, exist_ok=True)
    for name, url in FILES.items():
        fetch(os.path.join(RAW, name), url)
    print("\nDone.")


if __name__ == "__main__":
    main()
