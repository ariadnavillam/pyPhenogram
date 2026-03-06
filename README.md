# pyPhenogram

A Python tool for plotting chromosome ideograms annotated with GWAS hits. Outputs a pure SVG file with all 22 autosomes drawn to scale with coloured hit markers.

---

## Features

- All 22 autosomes drawn to scale (GRCh37/hg19)
- Giemsa cytogenetic banding (gneg, gpos, acen, etc.)
- Hit markers coloured by any column (e.g. gene name, trait, study)
- **Shape encodes significance**: diamond for genome-wide significant (p < 5x10-8), circle for suggestive
- Automatic vertical stacking when hits overlap on the same chromosome

---

## Installation

Requires Python ≥ 3.12. Install dependencies with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

---

## Input file format

A tab-separated file with at least two columns: chromosome and position. Column names are matched case-insensitively. Additional columns are carried through and can be used for colouring.

| Column | Accepted names | Required |
|---|---|---|
| Chromosome | `CHR`, `CHROM`, `CHROMOSOME` | Yes |
| Position (bp) | `POS`, `POSITION`, `BP`, `START` | Yes |
| P-value | any name (pass via `--p-col`) | No |
| Colour label | any name (pass via `--color-col`) | No |

A `chr` prefix on chromosome values is stripped automatically (e.g. `chr10` → `10`). Only autosomes 1–22 are plotted.


## Usage

### Command line

```bash
uv run run.py -i <input> [options]
```

| Option | Default | Description |
|---|---|---|
| `-i`, `--input` | *(required)* | Path to the input file |
| `-o`, `--output` | `plots/phenogram` | Output path without extension (`.svg` is added) |
| `--color-col` | `GENENAME` | Column whose values determine marker colour |
| `--cmap` | *(HSV wheel)* | Matplotlib colormap name (e.g. `tab10`, `viridis`) |
| `--p-col` | *(none)* | Column containing p-values — enables shape encoding |
| `--p-gws` | `5e-8` | Genome-wide significance threshold for diamond markers |

### Examples

**Basic plot — all markers as circles, coloured by gene:**

```bash
uv run run.py -i data/input.txt
```

**Enable shape encoding (diamonds for GWS hits, circles for suggestive):**

```bash
uv run run.py -i data/input.txt --p-col P
```

**Custom output path and colour column:**

```bash
uv run run.py -i data/input.txt -o plots/my_study --color-col TRAIT --p-col P
```

**Use a matplotlib colormap:**

```bash
uv run run.py -i data/input.txt --cmap tab10 --p-col P
```

**Override the genome-wide significance threshold:**

```bash
uv run run.py -i data/input.txt --p-col P --p-gws 1e-7
```

---

## Python API

```python
from pyPhenogram.parser import load_hits
from pyPhenogram.plotter import plot_all_chromosomes

df = load_hits("data/input.txt")

svg_path, color_map = plot_all_chromosomes(
    df,
    color_col="GENENAME",
    output_path="plots/phenogram",
    cmap=None,          # or e.g. "tab10"
    p_col="P",          # None → all circles
    p_genome_wide=5e-8, # threshold for diamond markers
)
```

---

## Output

An SVG file is written to the specified output path. The figure contains:

- **Chromosomes**: two rows (chr 1–11 top, chr 12–22 bottom), heights proportional to GRCh37 sizes
- **Cytogenetic bands**: Giemsa staining colours
- **Hit markers**: coloured by the chosen label column; shape indicates significance when `--p-col` is provided
  - ◆ Diamond — genome-wide significant (p < 5×10⁻⁸ by default)
  - ● Circle — suggestive
- **Colour legend**: one entry per unique label
- **Shape legend**: shown when `--p-col` is used

---

## Requirements

| Package | Version |
|---|---|
| Python | ≥ 3.12 |
| pandas | ≥ 2.0 |
| requests | ≥ 2.0 |
| matplotlib | ≥ 3.10 |
