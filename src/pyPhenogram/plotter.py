"""Draw all 22 chromosome idiograms with Giemsa bands and GWAS hit markers.

Outputs a pure SVG file using only Python's standard library xml.etree module.
No Plotly or kaleido required.
"""

import colorsys
import xml.etree.ElementTree as ET
from pathlib import Path

from .cytoband_data import load_cytobands

# ---------------------------------------------------------------------------
# Layout constants (pixels)
# ---------------------------------------------------------------------------
MARGIN_L = 35
MARGIN_T = 30
MARGIN_B = 25
CHR_SPACING = 68  # horizontal distance between chromosome centres
CHR_HW = 7  # half-width of chromosome body
NARROW_HW = 2.5  # half-width at centromere midpoint (constriction)
MAX_H = 340  # height of chr1 (longest chromosome)
CAP_H = 7  # bezier cap protrusion above / below chromosome body
ROW_GAP = 50  # gap between row 0 labels and row 1 top caps
LABEL_H = 18  # space reserved below each chromosome for its number
DOT_OFFSET = 13      # horizontal distance from chr right edge to dot centre
DOT_R = 3.5          # radius of GWAS hit dot
ELBOW_X_OFFSET = 5   # horizontal stub before the elbow turn

CHR_COLOR = "#DCDCDC"
OUTLINE_COLOR = "#555555"

# Row bottom y-coordinates (SVG: y increases downward)
# "Bottom" = q-arm tip of the chromosome body
ROW0_BOTTOM = MARGIN_T + CAP_H + MAX_H
ROW1_BOTTOM = ROW0_BOTTOM + CAP_H + LABEL_H + ROW_GAP + CAP_H + MAX_H

ROWS = {
    0: list(range(1, 12)),
    1: list(range(12, 23)),
}


# ---------------------------------------------------------------------------
# Colour assignment
# ---------------------------------------------------------------------------


def assign_colors(labels, cmap=None):
    """Return {label: hex_colour} with one distinct colour per unique label.

    If *cmap* is None, colours are evenly spaced around the HSV hue wheel.
    If *cmap* is a matplotlib colormap name (e.g. "tab10", "viridis"),
    N colours are sampled evenly from that colormap where N = unique labels.
    Requires matplotlib only when *cmap* is provided.
    """
    unique = sorted(set(labels))
    n = len(unique)

    if cmap is not None:
        try:
            import matplotlib.colors as mcolors
            import matplotlib
        except ImportError:
            raise ImportError(
                "matplotlib is required to use --cmap. "
                "Install it with: uv add matplotlib"
            )
        colormap = matplotlib.colormaps[cmap]
        return {
            label: mcolors.to_hex(colormap(i / max(n - 1, 1)))
            for i, label in enumerate(unique)
        }

    # Default: evenly spaced HSV hues
    colors = {}
    for i, label in enumerate(unique):
        hue = i / n
        r, g, b = colorsys.hsv_to_rgb(hue, 0.75, 0.85)
        colors[label] = f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
    return colors


# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------


def _f(n):
    """Format a float to 2 decimal places for SVG attributes."""
    return f"{n:.2f}"


def _chr_params(chrom_num, chr_sizes):
    """Return (cx, y_bottom, chr_h) in pixels for one chromosome."""
    max_size = chr_sizes[1]
    chr_h = chr_sizes[chrom_num] / max_size * MAX_H
    row = 0 if chrom_num <= 11 else 1
    col = ROWS[row].index(chrom_num)
    cx = MARGIN_L + col * CHR_SPACING
    y_bottom = ROW0_BOTTOM if row == 0 else ROW1_BOTTOM
    return cx, y_bottom, chr_h


def _bp_to_y(bp, y_bottom, chr_h, chr_size):
    """Convert a genomic bp position to SVG y (y increases downward).

    bp=0       → top of chromosome body (visually the p-arm tip)
    bp=chr_size → y_bottom (visually the q-arm tip)
    """
    return y_bottom - (1.0 - bp / chr_size) * chr_h


def _stack_positions(ys, min_sep):
    """Spread y positions so no two adjacent values are closer than *min_sep*.

    Sorting order is preserved: values that are close together are pushed
    downward (increasing y) to avoid overlap.  Returns a list in the same
    order as the input.
    """
    if not ys:
        return []
    # Sort keeping original indices
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    stacked = [ys[i] for i in order]
    # Greedy forward pass: push each value down if it collides with the previous
    for i in range(1, len(stacked)):
        if stacked[i] - stacked[i - 1] < min_sep:
            stacked[i] = stacked[i - 1] + min_sep
    # Map back to original order
    result = [0.0] * len(ys)
    for new_pos, orig_i in enumerate(order):
        result[orig_i] = stacked[new_pos]
    return result


# ---------------------------------------------------------------------------
# SVG path builders
# ---------------------------------------------------------------------------


def _outline_path(cx, y_bottom, chr_h, cent_yt, cent_ym, cent_yb):
    """SVG path for the chromosome outline: bezier caps + centromere V."""
    hw, nw, cap = CHR_HW, NARROW_HW, CAP_H
    yt = y_bottom - chr_h  # top of chromosome body (smaller y in SVG)
    return (
        f"M {_f(cx-hw)},{_f(yt)} "
        f"Q {_f(cx)},{_f(yt-cap)} {_f(cx+hw)},{_f(yt)} "  # top cap
        f"L {_f(cx+hw)},{_f(cent_yt)} "
        f"L {_f(cx+nw)},{_f(cent_ym)} "  # centromere right V
        f"L {_f(cx+hw)},{_f(cent_yb)} "
        f"L {_f(cx+hw)},{_f(y_bottom)} "
        f"Q {_f(cx)},{_f(y_bottom+cap)} {_f(cx-hw)},{_f(y_bottom)} "  # bottom cap
        f"L {_f(cx-hw)},{_f(cent_yb)} "
        f"L {_f(cx-nw)},{_f(cent_ym)} "  # centromere left V
        f"L {_f(cx-hw)},{_f(cent_yt)} "
        f"Z"
    )


# ---------------------------------------------------------------------------
# Per-chromosome drawing
# ---------------------------------------------------------------------------


def _add_chromosome(g, chrom_num, bands_chr, hits, color_map, chr_sizes):
    """Append SVG elements for one chromosome to the group element *g*."""
    cx, y_bottom, chr_h = _chr_params(chrom_num, chr_sizes)
    chr_size = chr_sizes[chrom_num]
    y_top = y_bottom - chr_h

    # Centromere extent from acen bands
    acen = bands_chr[bands_chr["stain"] == "acen"]
    if len(acen) > 0:
        cent_start = int(acen["start"].min())
        cent_end = int(acen["end"].max())
    else:
        mid = chr_size // 2
        cent_start, cent_end = mid - 1_000_000, mid + 1_000_000

    cent_yt = _bp_to_y(cent_start, y_bottom, chr_h, chr_size)
    cent_yb = _bp_to_y(cent_end, y_bottom, chr_h, chr_size)
    cent_ym = (cent_yt + cent_yb) / 2

    outline_d = _outline_path(cx, y_bottom, chr_h, cent_yt, cent_ym, cent_yb)

    # 1. Gray chromosome background (filled outline path)
    bg = ET.SubElement(g, "path")
    bg.set("d", outline_d)
    bg.set("fill", CHR_COLOR)
    bg.set("stroke", "none")

    # 2. Cytogenetic band rectangles
    for _, band in bands_chr.iterrows():
        by0 = _bp_to_y(band["start"], y_bottom, chr_h, chr_size)
        by1 = _bp_to_y(band["end"], y_bottom, chr_h, chr_size)
        # clamp to body rectangle
        by0 = max(by0, y_top)
        by1 = min(by1, y_bottom)
        h = by1 - by0
        if h <= 0:
            continue
        rect = ET.SubElement(g, "rect")
        rect.set("x", _f(cx - CHR_HW))
        rect.set("y", _f(by0))
        rect.set("width", _f(CHR_HW * 2))
        rect.set("height", _f(h))
        rect.set("fill", band["color"])
        rect.set("stroke", "none")

    # 3. White wedges that fill the centromere constriction notch
    for sign in (+1, -1):
        corner_x = cx + sign * CHR_HW
        notch_x = cx + sign * NARROW_HW
        bow = ET.SubElement(g, "path")
        bow.set(
            "d",
            (
                f"M {_f(corner_x)},{_f(cent_yt)} "
                f"L {_f(notch_x)},{_f(cent_ym)} "
                f"L {_f(corner_x)},{_f(cent_yb)} "
                f"Z"
            ),
        )
        bow.set("fill", "white")
        bow.set("stroke", "none")

    # 4. Chromosome outline (stroke only, on top)
    outline = ET.SubElement(g, "path")
    outline.set("d", outline_d)
    outline.set("fill", "none")
    outline.set("stroke", OUTLINE_COLOR)
    outline.set("stroke-width", "0.8")

    # 5. Chromosome number label
    label_y = y_bottom + CAP_H + LABEL_H - 4
    lbl = ET.SubElement(g, "text")
    lbl.set("x", _f(cx))
    lbl.set("y", _f(label_y))
    lbl.set("text-anchor", "middle")
    lbl.set("font-family", "sans-serif")
    lbl.set("font-size", "9")
    lbl.set("fill", "#333333")
    lbl.text = str(chrom_num)

    # 6. GWAS hit dots with elbow connectors
    # Dots are stacked to avoid overlap; elbow lines connect each dot back
    # to its true position on the chromosome.
    dot_x   = cx + CHR_HW + DOT_OFFSET
    elbow_x = cx + CHR_HW + ELBOW_X_OFFSET
    min_sep = DOT_R * 2 + 2   # minimum pixel gap between dot centres

    if hits:
        y_chrs = [_bp_to_y(pos, y_bottom, chr_h, chr_size) for pos, _ in hits]
        y_dots = _stack_positions(y_chrs, min_sep)

        for (pos_bp, label), y_chr, y_dot in zip(hits, y_chrs, y_dots):
            color = color_map.get(label, "#333333")

            # Elbow path: chr edge → short stub → vertical → dot
            connector = ET.SubElement(g, "path")
            connector.set("d", (
                f"M {_f(cx + CHR_HW)},{_f(y_chr)} "
                f"L {_f(elbow_x)},{_f(y_chr)} "
                f"L {_f(elbow_x)},{_f(y_dot)} "
                f"L {_f(dot_x)},{_f(y_dot)}"
            ))
            connector.set("fill", "none")
            connector.set("stroke", color)
            connector.set("stroke-width", "0.9")

            dot = ET.SubElement(g, "circle")
            dot.set("cx", _f(dot_x))
            dot.set("cy", _f(y_dot))
            dot.set("r",  _f(DOT_R))
            dot.set("fill",   color)
            dot.set("stroke", "white")
            dot.set("stroke-width", "0.7")


# ---------------------------------------------------------------------------
# Legend
# ---------------------------------------------------------------------------


def _add_legend(svg, color_map, color_col, x, y_start):
    """Append a colour legend to the SVG."""
    g = ET.SubElement(svg, "g")

    title = ET.SubElement(g, "text")
    title.set("x", _f(x))
    title.set("y", _f(y_start))
    title.set("font-family", "sans-serif")
    title.set("font-size", "10")
    title.set("font-weight", "bold")
    title.set("fill", "#333333")
    title.text = color_col

    for i, (label, color) in enumerate(sorted(color_map.items())):
        y = y_start + 18 + i * 18
        dot = ET.SubElement(g, "circle")
        dot.set("cx", _f(x + 6))
        dot.set("cy", _f(y - 4))
        dot.set("r", "5")
        dot.set("fill", color)

        txt = ET.SubElement(g, "text")
        txt.set("x", _f(x + 16))
        txt.set("y", _f(y))
        txt.set("font-family", "sans-serif")
        txt.set("font-size", "10")
        txt.set("fill", "#333333")
        txt.text = str(label)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def plot_all_chromosomes(df, color_col, output_path, cmap=None):
    """
    Draw all 22 chromosomes with Giemsa bands and GWAS hit markers to SVG.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns 'CHR' (int), 'POS' (int), and *color_col*.
    color_col : str
        Column whose values determine dot colours.
    output_path : str or Path
        Destination path (extension will be set to .svg).
    cmap : str or None
        Matplotlib colormap name (e.g. "tab10", "viridis").  When None the
        default HSV hue-wheel spacing is used.

    Returns
    -------
    Path
        Path to the saved SVG file.
    dict
        Color map used: { label -> colour }.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    bands = load_cytobands()
    chr_sizes = {c: int(bands[bands["chrom"] == c]["end"].max()) for c in range(1, 23)}
    color_map = assign_colors(df[color_col].tolist(), cmap=cmap)

    # Figure dimensions
    x_last_chr = MARGIN_L + 10 * CHR_SPACING
    legend_x = x_last_chr + CHR_HW + DOT_OFFSET + DOT_R + 20
    fig_w = legend_x + 150
    fig_h = ROW1_BOTTOM + CAP_H + LABEL_H + MARGIN_B

    svg = ET.Element("svg")
    svg.set("xmlns", "http://www.w3.org/2000/svg")
    svg.set("width", _f(fig_w))
    svg.set("height", _f(fig_h))

    # White background
    bg = ET.SubElement(svg, "rect")
    bg.set("x", "0")
    bg.set("y", "0")
    bg.set("width", _f(fig_w))
    bg.set("height", _f(fig_h))
    bg.set("fill", "white")

    g_chroms = ET.SubElement(svg, "g")
    g_chroms.set("id", "chromosomes")

    for chrom_num in range(1, 23):
        bands_chr = bands[bands["chrom"] == chrom_num]
        chrom_df = df[df["CHR"] == chrom_num]
        hits = [(row["POS"], row[color_col]) for _, row in chrom_df.iterrows()]
        _add_chromosome(g_chroms, chrom_num, bands_chr, hits, color_map, chr_sizes)

    _add_legend(svg, color_map, color_col, legend_x, MARGIN_T + 10)

    svg_path = output_path.with_suffix(".svg")
    tree = ET.ElementTree(svg)
    with open(svg_path, "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write(ET.tostring(svg, encoding="unicode"))

    print(f"  Saved SVG → {svg_path}")
    return svg_path, color_map
