"""Parse GWAS hit files for chromosome idiogram plotting."""

import pandas as pd


def load_hits(filepath, sep="\t"):
    """
    Load GWAS hits from a tab-separated file.

    Expected columns (case-insensitive): CHR, POS.
    Any additional columns (e.g. GENENAME) are preserved for coloring.

    Returns a DataFrame with CHR as int and POS as int,
    limited to autosomes 1-22.
    """
    df = pd.read_csv(filepath, sep=sep)

    # Normalise column names so callers can use uppercase keys
    df.columns = [c.strip() for c in df.columns]

    # Find CHR and POS columns (flexible matching)
    chr_col = _find_column(df, ["CHR", "CHROM", "CHROMOSOME"])
    pos_col = _find_column(df, ["POS", "POSITION", "BP", "START"])

    # Rename to standard names
    df = df.rename(columns={chr_col: "CHR", pos_col: "POS"})

    # Strip "chr" prefix if present (e.g. "chr10" -> 10)
    df["CHR"] = df["CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce")
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    # Drop rows with missing data and restrict to autosomes 1-22.
    # Also drop POS == 0, which arises when the input file has long lines
    # that are split across multiple physical lines (the continuation row
    # maps column SNPID → blank, CHR → some unrelated number, POS → 0).
    df = df.dropna(subset=["CHR", "POS"])
    df = df[df["POS"] > 0]
    df = df[df["CHR"].between(1, 22)].copy()
    df["CHR"] = df["CHR"].astype(int)
    df["POS"] = df["POS"].astype(int)

    return df


def _find_column(df, candidates):
    """Return the first column name in *candidates* found in *df* (case-insensitive)."""
    lower_map = {c.lower(): c for c in df.columns}
    for name in candidates:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    raise ValueError(
        f"Could not find any of columns {candidates} in file. "
        f"Available columns: {list(df.columns)}"
    )
