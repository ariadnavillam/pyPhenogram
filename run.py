#!/usr/bin/env python
"""
Run the pyPhenogram chromosome idiogram plotter.

Usage
-----
    uv run run.py -i data/input.txt
    uv run run.py -i data/input.txt -o plots/phenogram --color-col GENENAME
    uv run run.py -i data/input.txt --cmap tab10
"""

import argparse

from pyPhenogram.parser import load_hits
from pyPhenogram.plotter import P_GENOME_WIDE, plot_all_chromosomes


def parse_args():
    p = argparse.ArgumentParser(description="Plot GWAS hits on chromosome idiograms.")
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Tab-separated input file with at least CHR and POS columns.",
    )
    p.add_argument(
        "-o",
        "--output",
        default="plots/phenogram",
        help="Output file path without extension (default: plots/phenogram). "
        "The .svg extension is added automatically.",
    )
    p.add_argument(
        "--color-col",
        default="GENENAME",
        help="Column used to colour the dots (default: GENENAME).",
    )
    p.add_argument(
        "--cmap",
        default=None,
        help="Matplotlib colormap name for dot colours (e.g. tab10, viridis). "
        "Requires matplotlib. Overrides the default HSV colour scheme.",
    )
    p.add_argument(
        "--p-col",
        default=None,
        metavar="COLUMN",
        help="Column containing p-values (e.g. P).  When supplied, markers "
        "with p < the genome-wide threshold are drawn as diamonds; "
        "all others as circles.",
    )
    p.add_argument(
        "--p-gws",
        type=float,
        default=P_GENOME_WIDE,
        metavar="THRESHOLD",
        help=f"Genome-wide significance threshold (default {P_GENOME_WIDE}). "
        "Hits below this p-value are shown as diamonds.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    print(f"Loading hits from: {args.input}")
    df = load_hits(args.input)
    print(f"  {len(df)} hit(s) on chromosome(s): {sorted(df['CHR'].unique())}")

    # Fall back to a neutral label when the colour column is missing
    if args.color_col not in df.columns:
        print(
            f"  Column '{args.color_col}' not found – "
            "dots will all use the same colour."
        )
        df["_label"] = "hit"
        color_col = "_label"
    else:
        color_col = args.color_col

    print(f"\nPlotting all chromosomes → {args.output}.svg")
    svg_path, color_map = plot_all_chromosomes(
        df,
        color_col,
        args.output,
        cmap=args.cmap,
        p_col=args.p_col,
        p_genome_wide=args.p_gws,
    )

    print(f"\nColour legend ({color_col}):")
    for label, color in color_map.items():
        print(f"  {label}: {color}")

    print(f"\nDone – plot saved to {svg_path}")


if __name__ == "__main__":
    main()
