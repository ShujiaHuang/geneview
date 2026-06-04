# geneview User Guide

**geneview** is a Python library for genomics data visualization, built on matplotlib and pandas. It provides publication-ready plots for GWAS, population genetics, set comparisons, chromosome karyotypes, and genome track browser-style visualizations.

---

## Table of Contents

- [geneview User Guide](#geneview-user-guide)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
  - [Quick Start](#quick-start)
  - [GWAS Plots](#gwas-plots)
    - [Manhattan Plot](#manhattan-plot)
    - [Q-Q Plot](#q-q-plot)
    - [Q-Q Normal Plot](#q-q-normal-plot)
  - [Venn Diagrams](#venn-diagrams)
  - [Population Structure (ADMIXTURE)](#population-structure-admixture)
  - [Karyotype Plots](#karyotype-plots)
  - [Color Palettes](#color-palettes)
  - [Genome Tracks](#genome-tracks)
    - [Overview](#overview)
    - [Quick Example](#quick-example)
    - [Display Parameters (Gviz-Compatible API)](#display-parameters-gviz-compatible-api)
    - [Enhanced GenomeAxisTrack](#enhanced-genomeaxistrack)
    - [Embedding Tracks in Existing Figures](#embedding-tracks-in-existing-figures)
  - [Plot Styles](#plot-styles)
    - [Available Styles](#available-styles)
    - [Applying a Style Globally](#applying-a-style-globally)
    - [Using a Style as a Context Manager](#using-a-style-as-a-context-manager)
    - [Passing a Style to a Plot Function](#passing-a-style-to-a-plot-function)
    - [Style Comparison Table](#style-comparison-table)
    - [Creating a Custom Style](#creating-a-custom-style)
  - [Command-Line Interface (CLI)](#command-line-interface-cli)
    - [Applying a Plot Style via CLI](#applying-a-plot-style-via-cli)
  - [Utilities](#utilities)
    - [Loading Built-in Datasets](#loading-built-in-datasets)
    - [Text Adjustment](#text-adjustment)
    - [Numeric Check](#numeric-check)
  - [Saving Figures](#saving-figures)
  - [API Reference](#api-reference)
    - [Top-level Functions](#top-level-functions)
    - [Track Classes](#track-classes)
    - [File I/O](#file-io)
  - [Examples](#examples)

---

## Installation

```bash
pip install geneview

# With optional genome tracks dependencies (for BigWig, BAM, and CRAM support):
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
    available_display_params,
)
```

| Track Type | Purpose |
|------------|--------|
| `GenomeAxisTrack` | Genomic coordinate ruler with auto-formatted labels |
| `IdeogramTrack` | Chromosome ideogram with cytobands (auto-loads from geneview-data) |
| `AnnotationTrack` | Generic genomic ranges (boxes, ellipses, arrows) |
| `GeneRegionTrack` | Gene models (CDS/UTR/introns, UCSC style, strand coloring, intron chevron arrows) |
| `DataTrack` | Numeric data (line, histogram, heatmap, etc.) |
| `HighlightTrack` | Cross-track highlight regions |
| `OverlayTrack` | Overlay multiple tracks on the same panel |

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

### Display Parameters (Gviz-Compatible API)

All track constructors accept display parameters as keyword arguments, mirroring Gviz's `DisplayPars` system. You can also use Gviz-style dot-notation aliases:

```python
from geneview.genometracks import available_display_params

# Query default display parameters
base_params = available_display_params()
axis_params = available_display_params("GenomeAxisTrack")

# Constructor kwargs as display params (Gviz convention)
atrack = AnnotationTrack(data, name="Styled",
                         col="darkgreen", fill="lightgreen",
                         fontsize=14, alpha=0.8)

# Dot-notation aliases (Gviz "col.title" → Python col_title)
atrack2 = AnnotationTrack(data, name="Alias",
                          **{"col.title": "navy", "background.panel": "#FFF8DC"})
# Equivalent to:
atrack2 = AnnotationTrack(data, name="Alias",
                          col_title="navy", background_panel="#FFF8DC")

# Modify display params after creation
atrack.set_param("col", "red")
atrack.set_params({"alpha": 0.5, "fontsize": 16})
print(atrack.get_param("col"))  # 'red'
```

### Enhanced GenomeAxisTrack

The axis track supports explicit tick positions and forced label units:

```python
# Force labels in Mb (exponent=6) instead of auto-detected kb
gtrack_mb = GenomeAxisTrack(exponent=6)

# Place ticks at explicit positions
gtrack_custom = GenomeAxisTrack(
    ticks_at=[26_500_000, 26_550_000, 26_600_000, 26_650_000]
)
```

### Embedding Tracks in Existing Figures

Use `panel_only=True` to render tracks into pre-existing matplotlib axes (e.g. GridSpec layouts):

```python
fig, axes = plt.subplots(2, 1, figsize=(14, 6))
plot_tracks([gtrack, data_track], region=region, panel_only=True, ax=axes[0])
plot_tracks([gtrack, annotation], region=region, panel_only=True, ax=axes[1])
```

Additional margin controls (`margin`, `inner_margin`) mirror Gviz's pixel-based spacing.

For a comprehensive guide with all track types, display parameters, file I/O, and advanced usage, see the [Genome Tracks Guide](./genome_tracks_guide.md).

---

## Plot Styles

geneview ships with a built-in style system that lets you produce figures compliant with the requirements of major scientific journals — **Nature**, **Science**, and **Cell** — with a single function call. Each style configures fonts, font sizes, figure dimensions, tick marks, spine visibility, colour palettes, and export settings (e.g. TrueType font embedding, DPI) so you do not have to set them manually.

### Available Styles

| Style name | Description | Key properties |
| :------------: | :-------------: | :----------------: |
| `"geneview"` | Default style — clean, readable genomics figures for exploration and general publication. | Arial / Lucida Sans; 10–12 pt text; top/right spines hidden; 300 dpi; legacy geneview colour palette. |
| `"nature"` | Nature Research Figure Guide compliant. | Arial / Helvetica; 5–7 pt text; no gridlines; Wong colour-blind-safe palette; 450 dpi raster export. |
| `"science"` | AAAS *Science* guidelines. | Arial / Helvetica; 6–10 pt text; Okabe–Ito palette; 600 dpi line-art export. |
| `"cell"` | Cell Press illustration guidelines. | Arial / Helvetica; 6–8 pt text; accessible colour palette; 300 dpi. |

List all registered styles at any time:

```python
from geneview.plotstyle import list_styles
print(list_styles())
# ['cell', 'geneview', 'nature', 'science']
```

### Applying a Style Globally

Call `apply_style()` once at the top of your script; all subsequent plots will use that style:

```python
import geneview as gv
from geneview.plotstyle import apply_style

apply_style("nature")  # all plots from here on use Nature style

df = gv.load_dataset("gwas")
ax = gv.manhattanplot(data=df[["#CHROM", "POS", "P"]])
```

To revert to the default:

```python
apply_style("geneview")
```

### Using a Style as a Context Manager

Use `use_style()` as a `with` block to apply a style temporarily. The previous style is automatically restored when the block exits:

```python
from geneview.plotstyle import use_style

# Default geneview style is active here
with use_style("science"):
    ax = gv.qqplot(data=df["P"])       # drawn with Science style
    plt.savefig("qq_science.pdf")

# geneview default is restored here
ax = gv.manhattanplot(data=df[["#CHROM", "POS", "P"]])
```

### Passing a Style to a Plot Function

Every major plotting function accepts an optional `style=` parameter, so you can apply a style on a per-call basis without touching global settings:

```python
# Manhattan plot in Nature style
ax = gv.manhattanplot(data=df, style="nature")

# Q-Q plot in Science style
ax = gv.qqplot(data=df["P"], style="science")

# Q-Q normal plot in Cell style
ax = gv.qqnorm(data=df["P"].values, style="cell")

# Venn diagram in Nature style
ax = gv.venn(data={"A": {1,2,3}, "B": {2,3,4}, "C": {3,4,5}}, style="nature")

# Admixture plot in Cell style
ax = gv.admixtureplot(data=admixture_dict, style="cell")

# Genome tracks in Nature style
from geneview.genometracks import plot_tracks, GenomeAxisTrack, IdeogramTrack, GenomicInterval
region = GenomicInterval("chr7", 20_000_000, 60_000_000)
axes = plot_tracks([IdeogramTrack(chromosome="chr7"), GenomeAxisTrack()], region=region, style="nature")
```

Pass `style=None` (the default) to use whatever style is currently active globally.

### Style Comparison Table

The table below summarises the key differences between the built-in styles:

| Parameter | `geneview` | `nature` | `science` | `cell` |
|-----------|:--------:|:------:|:--------:|:----:|
| Font family | Arial, Lucida Sans, DejaVu Sans | Arial, Helvetica | Arial, Helvetica | Arial, Helvetica |
| Title font size (pt) | 12 | 7 | 10 | 8 |
| Label font size (pt) | 10 | 7 | 8 | 7 |
| Tick font size (pt) | 10 | 5.5 | 6 | 6 |
| Figure width (in) | 9.0 | 3.50 | 2.36 | 3.35 |
| Top spine visible | No | No | No | No |
| Right spine visible | No | No | No | No |
| Grid lines | No | No | No | No |
| Raster export DPI | 300 | 450 | 600 | 300 |
| Colour palette | geneview legacy | Wong (8 colours) | Okabe–Ito (8 colours) | Cell accessible (8 colours) |
| `pdf.fonttype` | 42 | 42 | 42 | 42 |

### Creating a Custom Style

You can define and register your own style:

```python
from geneview.plotstyle import PlotStyle, register_style, apply_style

my_style = PlotStyle(
    name="my_journal",
    description="Custom style for my lab's journal",
    font_size_title=9.0,
    font_size_label=8.0,
    font_size_tick=6.5,
    figure_figsize=(4.0, 3.0),
    axes_linewidth=0.5,
    color_palette=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"],
)

register_style(my_style)
apply_style("my_journal")
```

Once registered, a custom style appears in `list_styles()` and can be used anywhere:

```python
ax = gv.manhattanplot(data=df, style="my_journal")
```

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

# Genome tracks from BED, GFF, and bedGraph files
geneview tracks --region chr7:26490000-26720000 \
    --ideogram \
    -a cpg_islands.bed \
    -g gene_models.gtf \
    -d coverage.bedgraph \
    -o genome_tracks.png
```

### Applying a Plot Style via CLI

Every subcommand accepts the `--style` flag, which applies a built-in journal-compliant plot style to the figure:

```bash
# Manhattan plot in Nature style
geneview manhattan -i gwas.assoc -o manhattan.png --style nature

# Q-Q plot in Science style
geneview qq -i gwas.assoc -o qq.png --style science

# Venn diagram in Cell style
geneview venn -i list1.txt list2.txt -o venn.png --style cell

# Admixture plot in Nature style
geneview admixture -i output.Q -p pop.info -o admixture.png --style nature

# Genome tracks in Science style
geneview tracks --region chr7:26490000-26720000 \
    --ideogram \
    -a cpg_islands.bed \
    -g gene_models.gtf \
    -d coverage.bedgraph \
    --style science \
    -o genome_tracks.png
```

The available styles are: `geneview` (default), `nature`, `science`, and `cell`.

Run `geneview --help` for a full list of subcommands and options.

```text
subcommands:
  manhattan    Create a Manhattan plot from GWAS association results.
  qq           Create a Q-Q plot from GWAS association results.
  venn         Create a Venn diagram from 2-6 input files.
  admixture    Create an Admixture plot from ADMIXTURE .Q output.
  tracks       Create a genome track plot from BED, GFF, or bedGraph files.
```

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
| `apply_style()` | `geneview.plotstyle` | Apply a plot style globally |
| `use_style()` | `geneview.plotstyle` | Context manager for temporary style |
| `list_styles()` | `geneview.plotstyle` | List all registered styles |
| `get_style()` | `geneview.plotstyle` | Retrieve a PlotStyle by name |
| `register_style()` | `geneview.plotstyle` | Register a custom PlotStyle |

### Track Classes

| Class | Description |
|-------|-------------|
| `GenomeAxisTrack` | Coordinate axis |
| `IdeogramTrack` | Chromosome ideogram (auto-loads karyotype) |
| `AnnotationTrack` | Genomic range annotations |
| `GeneRegionTrack` | Gene models (CDS/UTR/introns, UCSC style with strand coloring) |
| `DataTrack` | Numeric data tracks |
| `HighlightTrack` | Cross-track highlights |
| `OverlayTrack` | Overlay multiple tracks |
| `GenomicInterval` | Genomic region dataclass |

### File I/O

| Function | Description |
|----------|-------------|
| `read_bed` | Read BED files |
| `read_gff` | Read GFF/GTF files |
| `read_bedgraph` | Read bedGraph files |
| `read_bigwig` | Read BigWig files (requires pyBigWig) |
| `read_bam_coverage` | Compute coverage from BAM files (requires pysam) |
| `read_cram_coverage` | Compute coverage from CRAM files (requires pysam + reference FASTA) |
| `read_auto` | Auto-detect format from file extension |

---

## Examples

Full example scripts are available in the [`examples/scripts/`](../examples/scripts/) directory:

| Script | Description |
|--------|-------------|
| `manhattan.py` | Manhattan plot |
| `qq.py` | Q-Q plot |
| `venn.py` | Venn diagrams (2-6 sets) |
| `admixture.py` | ADMIXTURE plot |
| `plotstyle.py` | Plot styles (Nature, Science, Cell) |
| `genome_tracks_basic.py` | Basic genome tracks |
| `genome_tracks_advanced.py` | IdeogramTrack and advanced features |
| `genome_tracks_gene_region.py` | Gene region tracks |
| `genome_tracks_data.py` | Data track plot types |
| `genome_tracks_highlight.py` | Highlight regions |
| `genome_tracks_comprehensive.py` | Full showcase |
| `genome_tracks_bam_cram.py` | BAM/CRAM coverage as DataTrack |
| `genome_tracks_io.py` | File I/O with real test data (BED, bedGraph, BigWig, BAM, GTF, GFF3) |
| `genome_tracks_multi_sample.py` | Multi-sample heatmap and overlay |
| `genome_tracks_sequence.py` | SequenceTrack (nucleotide display) |
| `genome_tracks_alignments.py` | AlignmentsTrack (BAM/CRAM pileup, sashimi) |
| `genome_tracks_details.py` | DetailsAnnotationTrack (detail panels) |
| `genome_tracks_data_extended.py` | Extended DataTrack types (average, confint, smooth, horizon, grid, regression) |
| `genome_tracks_export.py` | Exporting tracks to BED, GFF, bedGraph, WIG |
| `genome_tracks_annotation_enhanced.py` | AnnotationTrack enhancements (just_group, overplotting, merge_groups, from_bam) |
| `genome_tracks_gene_region_enhanced.py` | GeneRegionTrack enhancements (exon_annotation, gene_symbols, shortest) |
| `genome_tracks_schemes.py` | Color schemes and plot_tracks params (cex, add, ylim, scheme) |

Run any example:

```bash
python examples/scripts/genome_tracks_basic.py
```
