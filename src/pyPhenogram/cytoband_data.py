"""Download, cache, and serve UCSC hg19 cytogenetic band data."""

import gzip
import io
from pathlib import Path

import pandas as pd
import requests

UCSC_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz"
)
CACHE_DIR = Path.home() / ".cache" / "pyPhenogram"
CACHE_FILE = CACHE_DIR / "cytoBand_hg19.tsv"

# Giemsa stain → display colour
BAND_COLORS = {
    "gneg":    "#FFFFFF",   # light (euchromatin)
    "gpos25":  "#C8C8C8",
    "gpos50":  "#A0A0A0",
    "gpos75":  "#686868",
    "gpos100": "#404040",   # dark (heterochromatin)
    "acen":    "#CC0000",   # centromere
    "gvar":    "#E0E0E0",   # variable heterochromatin
    "stalk":   "#7090C0",   # satellite stalk
}


def load_cytobands():
    """
    Return a DataFrame with columns:
        chrom (int), start (int), end (int), name (str), stain (str), color (str)

    Data is downloaded from UCSC on the first call and cached at
    ~/.cache/pyPhenogram/cytoBand_hg19.tsv.
    """
    if not CACHE_FILE.exists():
        _download_cytobands()

    df = pd.read_csv(CACHE_FILE, sep="\t")
    return df


def _download_cytobands():
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Downloading cytogenetic band data from UCSC → {CACHE_FILE}")

    response = requests.get(UCSC_URL, timeout=60)
    response.raise_for_status()

    with gzip.open(io.BytesIO(response.content)) as gz:
        raw = gz.read().decode()

    rows = []
    for line in raw.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        chrom_str, start, end, name, stain = parts[:5]
        if not chrom_str.startswith("chr"):
            continue
        chrom_num = chrom_str.replace("chr", "")
        # Only autosomes 1-22
        try:
            chrom_int = int(chrom_num)
        except ValueError:
            continue
        if not (1 <= chrom_int <= 22):
            continue
        rows.append(
            {
                "chrom": chrom_int,
                "start": int(start),
                "end": int(end),
                "name": name,
                "stain": stain,
                "color": BAND_COLORS.get(stain, "#CCCCCC"),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(CACHE_FILE, sep="\t", index=False)
    print(f"  Saved {len(df)} bands.")
    return df
