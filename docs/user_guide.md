# geneview User Guide

**geneview** is a Python library for genomics data visualization, built on matplotlib and pandas. It provides publication-ready plots for GWAS, population genetics, set comparisons, chromosome karyotypes, and genome track browser-style visualizations.

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [GWAS Plots](#gwas-plots)
4. [Venn Diagrams](#venn-diagrams)
5. [Population Structure (ADMIXTURE)](#population-structure-admixture)
6. [Karyotype Plots](#karyotype-plots)
7. [Color Palettes](#color-palettes)
8. [Genome Tracks](#genome-tracks)
9. [Command-Line Interface (CLI)](#command-line-interface-cli)
10. [Utilities](#utilities)
11. [Saving Figures](#saving-figures)

---

## Installation

```bash
pip install geneview

# With optional genome tracks dependencies (for BigWig, BAM support):
pip install geneview[genometracks]
```

**Requirements:** Python >= 3.8, matplotlib, pandas, numpy, scipy.

---

## Quick Start

```python
import geneview as gv

# Load a built-in GWAS dataset
df = gv.load_dataset("gwas")

# Manhattan plot
gv.manhattanplot(data=df[["#CHROM", "POS", "P"]])

# Q-Q plot
gv.qqplot(df["P"])
```

---

## GWAS Plots

geneview provides two core plots for genome-wide association studies: the **Manhattan plot** and the **Q-Q plot**.

### Manhattan Plot

```python
import matplotlib.pyplot as plt
import geneview as gv

df = gv.utils.load_dataset("gwas")
gv.manhattanplot(
    data=df[["#CHROM", "POS", "P"]],
    hline_kws={"linestyle": "--", "lw": 1.3},
    xlabel="Chromosome",
    ylabel="-Log10(P-value)",
    xticklabel_kws={"rotation": 45},
)
plt.tight_layout()
plt.savefig("manhattan.png", dpi=300)
plt.show()
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `data` | DataFrame with columns `#CHROM`, `POS`, `P` (or specify with `chrom`, `pos`, `pv`) |
| `suggestiveline` | P-value threshold for suggestive line (default `1e-5`) |
| `genomewideline` | P-value threshold for genome-wide significance (default `5e-8`) |
| `sign_marker_p` | P-value threshold for highlighting significant SNPs |
| `sign_marker_color` | Color for significant markers (default `"r"`) |
| `is_annotate_topsnp` | Annotate the top SNP per chromosome |
| `xtick_label_set` | Set of chromosome labels to show on x-axis |
| `color` | Cycle colors for alternating chromosomes (comma-separated) |

### Q-Q Plot

```python
import matplotlib.pyplot as plt
import geneview as gv

df = gv.utils.load_dataset("gwas")
gv.qqplot(df["P"])
plt.tight_layout()
plt.savefig("qq.png", dpi=300)
plt.show()
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `data` | Array-like of P-values |
| `ax` | matplotlib Axes to draw on |
| `color` | Point color |
| `title` | Plot title |

### Q-Q Normal Plot

```python
gv.qqnorm(data_array)
```

Plots sample quantiles against theoretical normal quantiles.

---

## Venn Diagrams

Visualize set intersections for 2 to 6 sets.

```python
from itertools import chain, islice
from string import ascii_uppercase
from numpy.random import choice
import matplotlib.pyplot as plt
from geneview import venn

_, ax = plt.subplots(figsize=(8, 6))
dataset_dict = {
    name: set(choice(1000, 700, replace=False))
    for name in islice(ascii_uppercase, 3)  # 3 sets: A, B, C
}
venn(
    dataset_dict,
    fmt="{percentage:.1f}%",
    palette="cool",
    fontsize=12,
    legend_use_petal_color=True,
    ax=ax,
)
plt.savefig("venn.png", dpi=300)
plt.show()
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `data` | Dict mapping set names to sets of items |
| `fmt` | Label format: `"{size}"`, `"{percentage:.1f}%"`, or `"{logic}"` |
| `palette` | Color palette name or list of colors |
| `fontsize` | Label font size |
| `legend_use_petal_color` | Use petal colors in legend |

Supports 2, 3, 4, 5, and 6-set Venn diagrams.

---

## Population Structure (ADMIXTURE)

Visualize ADMIXTURE Q-matrix output as a stacked bar chart.

```python
import matplotlib.pyplot as plt
import geneview as gv

q_file = gv.load_dataset("admixture_output.Q")
pop_file = gv.load_dataset("admixture_population.info")

pop_order = ["KHV", "CDX", "CHS", "CHB", "JPT", "BEB", "STU", "ITU",
             "GIH", "PJL", "FIN", "CEU", "GBR", "IBS", "TSI", "PEL",
             "PUR", "MXL", "CLM", "ASW", "ACB", "GWD", "MSL", "YRI",
             "ESN", "LWK"]

fig, ax = plt.subplots(1, 1, figsize=(14, 2), facecolor="w",
                        constrained_layout=True)
gv.admixtureplot(
    data=q_file,
    population_info=pop_file,
    edgewidth=2.0,
    group_order=pop_order,
    shuffle_popsample_kws={"frac": 0.5},
    palette="Set1",
    xticklabel_kws={"rotation": 45},
    ylabel_kws={"rotation": 0, "ha": "right"},
    set_xticklabel_top=True,
    ax=ax,
)
plt.savefig("admixture.png", dpi=300)
plt.show()
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `data` | Path to ADMIXTURE `.Q` file or matrix |
| `population_info` | Path to population labels file |
| `group_order` | List of population names in desired order |
| `palette` | Color palette for ancestral components |
| `edgewidth` | Width of bar edges |

---

## Karyotype Plots

Display chromosome karyotypes with band annotations.

```python
import matplotlib.pyplot as plt
import geneview as gv

# Load karyotype data (or provide a file path / URL)
gv.karyoplot(data="path/to/karyotype.txt")
plt.savefig("karyotype.png", dpi=300)
plt.show()
```

**Key parameters:**

| Parameter | Description |
|-----------|-------------|
| `data` | Karyotype data (file path, URL, or list) |
| `CHR` | Plot a specific chromosome only |
| `width` | Chromosome width in plot |
| `alpha` | Opacity |
| `color4none` | Color for undefined bands |

---

## Color Palettes

geneview includes built-in color palettes for genomics figures.

```python
from geneview.palette import generate_colors_palette, circos, xkcd_rgb

# Generate a custom palette
colors = generate_colors_palette(n=10, palette="Set1")

# Circos-style colors
circos_colors = circos  # dict of named colors

# XKCD color survey colors
xkcd = xkcd_rgb  # dict of 949 named colors
```

Common palette names: `"Set1"`, `"Set2"`, `"Set3"`, `"pastel"`, `"deep"`, `"muted"`, `"bright"`, `"dark"`, `"colorblind"`.

---

## Genome Tracks

The genome tracks module provides Gviz-style track browser visualizations in Python. This is a powerful system for displaying genomic features along a shared coordinate axis.

### Overview

```python
from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack,
    DataTrack, HighlightTrack, GenomicInterval, plot_tracks,
)
```

| Track Type | Purpose |
|------------|---------|
| `GenomeAxisTrack` | Genomic coordinate ruler with auto-formatted labels |
| `IdeogramTrack` | Chromosome ideogram with cytobands (auto-loads from geneview-data) |
| `AnnotationTrack` | Generic genomic ranges (boxes, ellipses, arrows) |
| `GeneRegionTrack` | Gene models (exons, UTRs, introns) |
| `DataTrack` | Numeric data (line, histogram, heatmap, etc.) |
| `HighlightTrack` | Cross-track highlight regions |

### Quick Example

```python
import pandas as pd
from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, IdeogramTrack,
    GenomicInterval, plot_tracks,
)
import matplotlib.pyplot as plt

# Create annotation data
data = pd.DataFrame({
    "chrom": ["chr7"] * 4,
    "start": [2000000, 2070000, 2100000, 2160000],
    "end":   [2050000, 2130000, 2150000, 2170000],
    "strand": ["-", "+", "-", "-"],
    "name": ["feat1", "feat2", "feat3", "feat4"],
})

# Build tracks
itrack = IdeogramTrack(chromosome="chr7")  # Auto-loads hg38 karyotype
gtrack = GenomeAxisTrack()
atrack = AnnotationTrack(data, name="Features")

# Plot
region = GenomicInterval("chr7", 1900000, 2200000)
axes = plot_tracks([itrack, gtrack, atrack], region=region, figsize=(12, 5))
fig = axes[0].figure
fig.savefig("genome_tracks.png", dpi=300, bbox_inches="tight")
plt.show()
```

For a comprehensive guide with all track types, display parameters, file I/O, and advanced usage, see the [Genome Tracks Guide](./genome_tracks_guide.md).

---

## Command-Line Interface (CLI)

geneview can also be used as a command-line tool:

```bash
# Install
pip install geneview

# Manhattan plot from PLINK output
geneview manhattan --input gwas.assoc --output manhattan.png

# Q-Q plot
geneview qq --input gwas.assoc --output qq.png

# Venn diagram from gene list files
geneview venn --input list1.txt list2.txt list3.txt --output venn.png

# Admixture plot
geneview admixture --input output.Q --pop-file pop.info --output admixture.png
```

Run `geneview --help` for a full list of subcommands and options.

---

## Utilities

### Loading Built-in Datasets

```python
import geneview as gv

# List available datasets
names = gv.get_dataset_names()

# Load a dataset
df = gv.load_dataset("gwas")
```

### Text Adjustment

```python
from geneview.utils import adjust_text

# Avoid text label overlap (used internally by manhattanplot)
adjust_text(texts, ax=ax)
```

### Numeric Check

```python
from geneview.utils import is_numeric

is_numeric(42)      # True
is_numeric("hello") # False
```

---

## Saving Figures

All geneview plotting functions return matplotlib Axes objects. Use standard matplotlib methods to save figures:

```python
import matplotlib.pyplot as plt
import geneview as gv

df = gv.utils.load_dataset("gwas")
gv.manhattanplot(df[["#CHROM", "POS", "P"]])

# Get the current figure
fig = plt.gcf()

# Save in various formats
fig.savefig("manhattan.png", dpi=300, bbox_inches="tight")
fig.savefig("manhattan.pdf", bbox_inches="tight")
fig.savefig("manhattan.svg", bbox_inches="tight")
```

For genome tracks, get the figure from the axes list:

```python
axes = plot_tracks([gtrack, atrack], region=region)
fig = axes[0].figure
fig.savefig("tracks.png", dpi=300, bbox_inches="tight")
```

**Tips:**
- Use `bbox_inches="tight"` to remove extra whitespace
- Set `dpi=300` for publication-quality PNGs
- PDF/SVG formats preserve vector graphics for journal submission

---

## API Reference

### Top-level Functions

| Function | Module | Description |
|----------|--------|-------------|
| `manhattanplot()` | `geneview.gwas` | Manhattan plot for GWAS |
| `qqplot()` | `geneview.gwas` | Q-Q plot for P-values |
| `qqnorm()` | `geneview.gwas` | Q-Q normal plot |
| `venn()` | `geneview.baseplot` | Venn diagram (2-6 sets) |
| `admixtureplot()` | `geneview.popgene` | ADMIXTURE bar plot |
| `karyoplot()` | `geneview.karyotype` | Chromosome karyotype |
| `plot_tracks()` | `geneview.genometracks` | Genome track browser |

### Track Classes

| Class | Description |
|-------|-------------|
| `GenomeAxisTrack` | Coordinate axis |
| `IdeogramTrack` | Chromosome ideogram (auto-loads karyotype) |
| `AnnotationTrack` | Genomic range annotations |
| `GeneRegionTrack` | Gene model tracks |
| `DataTrack` | Numeric data tracks |
| `HighlightTrack` | Cross-track highlights |
| `GenomicInterval` | Genomic region dataclass |

---

## Examples

Full example scripts are available in the [`examples/scripts/`](../examples/scripts/) directory:

| Script | Description |
|--------|-------------|
| `manhattan.py` | Manhattan plot |
| `qq.py` | Q-Q plot |
| `venn.py` | Venn diagrams (2-6 sets) |
| `admixture.py` | ADMIXTURE plot |
| `genome_tracks_basic.py` | Basic genome tracks |
| `genome_tracks_advanced.py` | IdeogramTrack and advanced features |
| `genome_tracks_gene_region.py` | Gene region tracks |
| `genome_tracks_data.py` | Data track plot types |
| `genome_tracks_highlight.py` | Highlight regions |
| `genome_tracks_comprehensive.py` | Full showcase |

Run any example:

```bash
python examples/scripts/genome_tracks_basic.py
```
