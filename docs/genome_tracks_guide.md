# Genome Tracks Guide

The genome tracks module in geneview provides a powerful system for visualizing genomic features along a shared coordinate axis. Inspired by the R/Bioconductor [Gviz](https://bioconductor.org/packages/Gviz) package, it has been redesigned for Python with matplotlib rendering and pandas DataFrames for data handling.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Basic Features](#basic-features)
3. [Display Parameters](#display-parameters)
4. [GenomeAxisTrack](#genomeaxistrack)
5. [AnnotationTrack](#annotationtrack)
6. [GeneRegionTrack](#generegiontrack)
7. [DataTrack](#datatrack)
8. [HighlightTrack](#highlighttrack)
9. [OverlayTrack](#overlaytrack)
10. [File I/O](#file-io)
11. [Advanced: plot_tracks()](#advanced-plot_tracks)
12. [Display Parameters Reference](#display-parameters-reference)
13. [Complete Example](#complete-example)

---

## Introduction

### Core Concepts

The fundamental concept is similar to genome browsers: individual types of genomic features are represented by separate **track objects**. All tracks share a common genomic coordinate system, so alignment is handled automatically.

```
GenomeAxisTrack  ─── coordinate ruler
AnnotationTrack  ─── generic ranges (boxes/arrows)
GeneRegionTrack  ─── gene models (exons/UTRs/introns)
DataTrack        ─── numeric data (line/histogram/heatmap)
HighlightTrack   ─── cross-track highlights
OverlayTrack     ─── overlay multiple tracks on same axes
```

### Track Hierarchy

```
Track (abstract base)
  ├── GenomeAxisTrack
  └── RangeTrack (has a DataFrame)
        ├── StackedTrack (overlapping features stacked)
        │     ├── AnnotationTrack
        │     └── GeneRegionTrack
        └── NumericTrack (numeric value columns)
              └── DataTrack
        └── HighlightTrack (wrapper/container)
  └── OverlayTrack (container, overlays tracks)
```

### Design Philosophy

- **Hybrid API**: Lightweight Track classes for data + styling, with `plot_tracks()` as the function-based entry point (consistent with geneview)
- **matplotlib GridSpec**: Vertical track stacking via `matplotlib.gridspec`
- **GenomicInterval**: A simple dataclass replacing R's GRanges
- **Display parameters**: Python dicts with sensible defaults

---

## Basic Features

### Your First Track Plot

```python
import pandas as pd
from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GenomicInterval, plot_tracks,
)
import matplotlib.pyplot as plt

# Create some genomic features
data = pd.DataFrame({
    "chrom": ["chr7"] * 7,
    "start": [26500100, 26520000, 26560200, 26580000, 26610500, 26660000, 26700800],
    "end":   [26500800, 26520500, 26560900, 26580350, 26611200, 26660600, 26701500],
    "strand": ["+"] * 4 + ["-"] * 3,
    "name": ["CpG1", "CpG2", "CpG3", "CpG4", "CpG5", "CpG6", "CpG7"],
})

# Build tracks
atrack = AnnotationTrack(data, name="CpG Islands")
gtrack = GenomeAxisTrack()

# Define the region to display
region = GenomicInterval("chr7", 26500000, 26800000)

# Plot
axes = plot_tracks([gtrack, atrack], region=region)
fig = axes[0].figure
plt.show()
```

The result shows a title region with the track name on the left, and the data region on the right — similar to the UCSC Genome Browser layout.

### Adding a Genome Axis

The `GenomeAxisTrack` provides a coordinate reference. It is always relative to the other tracks being plotted.

```python
gtrack = GenomeAxisTrack()
axes = plot_tracks([gtrack, atrack], region=region)
```

### Auto Region Detection

If you don't specify a `region`, `plot_tracks()` will auto-derive it from the track data:

```python
# Region is automatically set to span all features
axes = plot_tracks([atrack])
```

### Zooming with extend_left / extend_right

Rather than specifying exact coordinates, you can extend or shrink the view relative to the data range:

```python
# Extend 100kb on each side
axes = plot_tracks([gtrack, atrack], extend_left=100000, extend_right=100000)

# Or use fractional values (0.1 = extend by 10% of the data range)
axes = plot_tracks([gtrack, atrack], extend_left=0.1, extend_right=0.1)
```

---

## Display Parameters

Every track has a `display_params` dictionary controlling its appearance. Parameters can be set at construction time or modified later.

### Setting Parameters at Construction

```python
atrack = AnnotationTrack(data, name="Features", display_params={
    "fill": "#3C5488",
    "col": "#333333",
    "fontsize": 9,
})
```

### Getting and Setting Parameters

```python
# Get a parameter
color = atrack.get_param("fill")

# Set a parameter
atrack.set_param("fill", "#E64B35")

# Set multiple parameters at once
atrack.set_params({"fill": "#E64B35", "col": "white"})
```

### Common Display Parameters (All Tracks)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 1.0 | Opacity of plot elements |
| `background_panel` | `"white"` | Background color of data panel |
| `background_title` | `"#E8E8E8"` | Background color of title panel |
| `col` | `"#3C5488"` | Line/border color |
| `fill` | `"#5B8DB8"` | Fill color |
| `col_border` | `"#333333"` | Border color |
| `fontsize` | 10 | Base font size |
| `fontsize_title` | 10 | Title font size |
| `lwd` | 1.0 | Line width |
| `show_title` | True | Show the title panel |
| `min_width` | 1 | Minimum feature width in pixels |
| `min_distance` | 1 | Minimum distance for stacking |
| `collapse` | True | Collapse overlapping features |

---

## GenomeAxisTrack

Displays a genomic coordinate ruler with tick marks and auto-formatted labels (bp/kb/Mb/Gb).

### Constructor

```python
GenomeAxisTrack(
    ranges=None,           # DataFrame of regions to highlight
    name="Axis",
    height=0.3,            # Relative height (small)
    scale=None,            # None = full axis; float = scale bar
    label_pos="alternating",  # "above", "below", "alternating"
    little_ticks=False,    # Add minor tick marks
    add53=False,           # Draw 5'->3' direction arrow
    add35=False,           # Draw 3'->5' direction arrow
    display_params=None,
)
```

### Full Axis

```python
gtrack = GenomeAxisTrack()
axes = plot_tracks([gtrack], region=GenomicInterval("chr7", 2000000, 2200000))
```

### Scale Bar Mode

Show a simple scale indicator instead of a full axis:

```python
# Fractional: use 30% of the plot width for the scale
gtrack = GenomeAxisTrack(scale=0.3, label_pos="below")

# Absolute: show a 50kb scale bar
gtrack = GenomeAxisTrack(scale=50000, label_pos="below")
```

### Label Positioning

```python
# All labels below the axis
gtrack = GenomeAxisTrack(label_pos="below")

# All labels above
gtrack = GenomeAxisTrack(label_pos="above")

# Alternating (default)
gtrack = GenomeAxisTrack(label_pos="alternating")
```

### Minor Ticks

```python
gtrack = GenomeAxisTrack(little_ticks=True)
```

### Direction Indicators

Show 5'→3' or 3'→5' direction arrows on the axis:

```python
# Draw 5'->3' arrow
gtrack = GenomeAxisTrack(add53=True)

# Draw 3'->5' arrow
gtrack = GenomeAxisTrack(add35=True)

# Both directions
gtrack = GenomeAxisTrack(add53=True, add35=True)
```

### Highlighting Regions on the Axis

```python
import pandas as pd
ranges = pd.DataFrame({
    "chrom": ["chr7", "chr7"],
    "start": [2050000, 2100000],
    "end":   [2070000, 2150000],
})
gtrack = GenomeAxisTrack(ranges=ranges)
axes = plot_tracks([gtrack], region=GenomicInterval("chr7", 2000000, 2200000))
```

### Display Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#555555"` | Axis line and tick color |
| `fontcolor` | `"#555555"` | Label text color |
| `fontsize` | 9 | Label font size |
| `lwd` | 2 | Axis line width |
| `col_range` | `"#8B8378"` | Highlighted region border |
| `fill_range` | `"#CDC8B1"` | Highlighted region fill |
| `col_53` | `"#555555"` | 5'→3' arrow color |
| `col_35` | `"#555555"` | 3'→5' arrow color |

---

## AnnotationTrack

Displays genomic ranges as shapes (boxes, ellipses, arrows) on a track. Supports grouping, strand direction, feature coloring, and stacking.

### Constructor

```python
AnnotationTrack(
    data,                  # DataFrame or file path (BED/GFF)
    stacking="squish",     # "squish", "pack", "dense", "full", "hide"
    shape="box",           # "box", "ellipse", "arrow", "fixedArrow", "smallArrow"
    show_label=False,      # Show feature names
    label_pos="above",     # "above" or "below"
    group_annotation=None, # Group features: "id", "group", "feature"
    name="Annotation",
    height=1.0,
    display_params=None,
)
```

### Basic Usage

```python
import pandas as pd
from geneview.genometracks import AnnotationTrack, GenomicInterval, plot_tracks

data = pd.DataFrame({
    "chrom": ["chr7"] * 4,
    "start": [2000000, 2070000, 2100000, 2160000],
    "end":   [2050000, 2130000, 2150000, 2170000],
    "strand": ["-", "+", "-", "-"],
    "name": ["feat1", "feat2", "feat3", "feat4"],
    "feature": ["exon", "CDS", "exon", "CDS"],
})
atrack = AnnotationTrack(data, name="Features")
axes = plot_tracks([atrack], region=GenomicInterval("chr7", 1900000, 2200000))
```

### Shapes

```python
# Boxes (default)
atrack = AnnotationTrack(data, shape="box")

# Ellipses
atrack = AnnotationTrack(data, shape="ellipse")

# Arrows (strand-aware)
atrack = AnnotationTrack(data, shape="arrow")

# Fixed-width arrows (head width controlled by arrowHeadWidth param)
atrack = AnnotationTrack(data, shape="fixedArrow",
    display_params={"arrowHeadWidth": 0.4})

# Small arrows (50% head size)
atrack = AnnotationTrack(data, shape="smallArrow")
```

### Stacking Modes

When features overlap, stacking determines their vertical arrangement:

| Mode | Behavior |
|------|----------|
| `squish` | Stack overlapping features in separate rows (default) |
| `pack` | Similar to squish but more compact |
| `dense` | Merge overlapping features into one row |
| `full` | Each feature gets its own row |
| `hide` | Don't draw features |

```python
atrack = AnnotationTrack(data, stacking="squish")
```

### Feature Coloring

Different feature types can have different colors:

```python
atrack = AnnotationTrack(data, display_params={
    "fill": "#3C5488",          # Default fill
    "feature_color": True,      # Auto-color by feature type
})
```

### Labels

```python
atrack = AnnotationTrack(data, show_label=True, label_pos="above")
```

### Group Annotation

Draw connecting lines and labels between grouped features:

```python
# Group by 'feature' column and draw connecting lines
atrack = AnnotationTrack(data, group_annotation="feature",
    display_params={"col_line": "#888888"})

# Group by 'id' column
atrack = AnnotationTrack(data, group_annotation="id")

# Group by custom 'group' column
atrack = AnnotationTrack(data, group_annotation="group")
```

### Loading from Files

```python
from geneview.genometracks import read_bed, AnnotationTrack

# Load directly from BED file
bed_data = read_bed("features.bed")
atrack = AnnotationTrack(bed_data, name="BED Features")

# Or pass file path directly
atrack = AnnotationTrack("features.bed", name="BED Features")
```

---

## GeneRegionTrack

Displays gene models with exons as thick boxes, UTRs as thin boxes, and introns as connecting lines.

### Constructor

```python
GeneRegionTrack(
    data,                       # DataFrame or GFF/GTF file path
    stacking="squish",
    collapse_transcripts=False, # False, True/"gene", "longest", "shortest", "meta"
    show_id="gene",             # "gene", "transcript", "exon", or None
    thin_box_features=None,     # Set of feature types drawn as thin boxes
    name="GeneRegion",
    height=1.5,
    display_params=None,
)
```

### Basic Usage

```python
from geneview.genometracks import GeneRegionTrack, read_gff, GenomicInterval, plot_tracks

gene_data = read_gff("gene_models.gtf")
grtrack = GeneRegionTrack(gene_data, name="Gene Models")
region = GenomicInterval("chr7", 26490000, 26720000)
axes = plot_tracks([grtrack], region=region)
```

### Drawing Details

- **Thick boxes**: CDS (coding sequence) features
- **Thin boxes**: UTR features (5' UTR, 3' UTR, and other non-coding features)
- **Lines**: Introns connecting exons
- **Chevrons**: Small arrows inside exons indicating strand direction

### Transcript Collapsing

```python
# Show all transcripts (default)
grtrack = GeneRegionTrack(data, collapse_transcripts=False)

# Collapse to gene level (merge all exons per gene)
grtrack = GeneRegionTrack(data, collapse_transcripts="gene")
# or: collapse_transcripts=True

# Show only the longest transcript per gene
grtrack = GeneRegionTrack(data, collapse_transcripts="longest")

# Show only the shortest transcript per gene
grtrack = GeneRegionTrack(data, collapse_transcripts="shortest")

# Meta-transcript: union of all exon positions per gene
grtrack = GeneRegionTrack(data, collapse_transcripts="meta")
```

### Label Options

```python
# Show gene names (default)
grtrack = GeneRegionTrack(data, show_id="gene")

# Show transcript IDs
grtrack = GeneRegionTrack(data, show_id="transcript")

# No labels
grtrack = GeneRegionTrack(data, show_id=None)
```

### Display Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fill` | `"#3C5488"` | CDS exon fill color |
| `fill_utr` | `"#8DB4E2"` | UTR fill color |
| `col_intron` | `"#888888"` | Intron line color |
| `fontsize` | 8 | Label font size |
| `fontcolor` | `"#333333"` | Label text color |
| `lwd` | 0.8 | Line width |

---

## DataTrack

Visualizes numeric genomic data along coordinates. Supports multiple plot types and multi-sample data.

### Constructor

```python
DataTrack(
    data,              # DataFrame with chrom/start/end + value columns
    type="line",       # Plot type (see table below)
    ylim=None,         # Custom y-axis limits
    groups=None,       # Sample grouping factor
    transformation=None,  # Callable applied to values (e.g., np.log2)
    window=None,       # Windowing: int, "auto", or "fixed"
    window_size=None,  # Bin size for window="fixed"
    aggregation="mean", # Aggregation: "mean", "median", "sum", "min", "max"
    legend=False,      # Show legend for groups
    name="Data",
    height=1.5,
    display_params=None,
)
```

### Plot Types

| Type | Description |
|------|-------------|
| `"line"` | Connected line plot |
| `"histogram"` | Bar chart (bar width = range width) |
| `"polygon"` | Filled area plot relative to baseline |
| `"points"` | Scatter/dot plot |
| `"mountain"` | Mountain-type plot relative to baseline |
| `"gradient"` | False-color image of summarized values |
| `"heatmap"` | False-color image per sample |
| `"boxplot"` | Box-and-whisker plot |
| `"b"` | Combined line + points |
| `"s"` | Stair steps (horizontal first, then vertical) |
| `"S"` | Stair steps (vertical first, then horizontal) |

### Line Plot

```python
import pandas as pd
import numpy as np
from geneview.genometracks import DataTrack, GenomicInterval, plot_tracks

data = pd.DataFrame({
    "chrom": ["chr7"] * 100,
    "start": np.arange(26500000, 26800000, 3000),
    "end": np.arange(26503000, 26803000, 3000),
    "value": np.random.randn(100).cumsum(),
})
dtrack = DataTrack(data, type="line", name="Signal")
axes = plot_tracks([dtrack], region=GenomicInterval("chr7", 26500000, 26800000))
```

### Histogram

```python
dtrack = DataTrack(data, type="histogram", name="Coverage")
```

### Polygon (Area Fill)

```python
dtrack = DataTrack(data, type="polygon", name="Area")
```

### Points

```python
dtrack = DataTrack(data, type="points", name="Scatter")
```

### Mountain

```python
dtrack = DataTrack(data, type="mountain", name="Mountain")
```

### Gradient

```python
dtrack = DataTrack(data, type="gradient", name="Gradient")
```

### Heatmap (Multi-sample)

The heatmap type works best with multiple value columns:

```python
data = pd.DataFrame({
    "chrom": ["chr7"] * 50,
    "start": np.arange(26500000, 26800000, 6000),
    "end": np.arange(26506000, 26806000, 6000),
    "sample_A": np.random.randn(50).cumsum() / 3,
    "sample_B": np.random.randn(50).cumsum() / 3 + 2,
    "sample_C": np.random.randn(50).cumsum() / 3 - 1,
})
dtrack = DataTrack(data, type="heatmap", name="Heatmap")
```

### Custom Y-Axis Limits

```python
dtrack = DataTrack(data, type="line", ylim=(-10, 10))
```

### Combined Line + Points

```python
dtrack = DataTrack(data, type="b", name="Signal+Points")
```

### Stair Steps

```python
# Horizontal first, then vertical (post-step)
dtrack = DataTrack(data, type="s", name="Steps")

# Vertical first, then horizontal (pre-step)
dtrack = DataTrack(data, type="S", name="Reverse Steps")
```

### Transformation

Apply a function to all values before plotting:

```python
import numpy as np

# Log-transform values
dtrack = DataTrack(data, type="line", transformation=np.log2)

# Custom transformation
dtrack = DataTrack(data, type="line", transformation=lambda x: x**2)
```

### Windowing and Smoothing

Reduce noise by binning data into windows:

```python
# Bin into 50 equal windows across the range
dtrack = DataTrack(data, type="line", window=50)

# Running window of 10 bins
dtrack = DataTrack(data, type="line", window=-10)

# Auto-detect window size based on data density
dtrack = DataTrack(data, type="line", window="auto")

# Fixed-size windows of 20 bins
dtrack = DataTrack(data, type="line", window="fixed", window_size=20)
```

### Aggregation

Control how values within windows are combined:

```python
dtrack = DataTrack(data, window=50, aggregation="median")
dtrack = DataTrack(data, window=50, aggregation="sum")
dtrack = DataTrack(data, window=50, aggregation="max")
```

### Legend

Show a legend for grouped data:

```python
dtrack = DataTrack(data, type="line", groups=["A", "A", "B", "B"],
                   legend=True,
                   display_params={"group_colors": ["#3C5488", "#E64B35"]})
```

### Grid Lines

Draw horizontal grid lines:

```python
dtrack = DataTrack(data, type="line",
                   display_params={"grid": True, "col_grid": "#DDDDDD"})
```

### Display Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fill` | `"#5B8DB8"` | Fill color for histogram/polygon/mountain |
| `col` | `"#3C5488"` | Line color |
| `baseline` | 0 | Y-position of baseline for mountain/polygon |
| `ncolor` | 100 | Number of colors for gradient |

---

## HighlightTrack

A wrapper/container track that adds colored highlight regions across multiple tracks.

### Constructor

```python
HighlightTrack(
    regions,           # DataFrame or list of GenomicInterval
    track_list=None,   # List of tracks to highlight
    fill="#FFE3E6",    # Highlight fill color (str or list of str)
    col=None,          # Highlight border color (str or list of str)
    alpha=0.3,         # Transparency
    name="Highlight",
    height=1.0,
    display_params=None,
)
```

### Basic Usage

```python
import pandas as pd
from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, DataTrack, HighlightTrack,
    GenomicInterval, plot_tracks,
)

# Create tracks
gtrack = GenomeAxisTrack()
atrack = AnnotationTrack(data, name="Features")
dtrack = DataTrack(cov_data, type="line", name="Signal")

# Define highlight regions
highlights = pd.DataFrame({
    "chrom": ["chr7", "chr7"],
    "start": [26520000, 26640000],
    "end":   [26540000, 26660000],
})

# Wrap tracks in HighlightTrack
ht = HighlightTrack(
    regions=highlights,
    track_list=[atrack, dtrack],
    fill="#FFFF99",
    alpha=0.35,
)

# Plot -- HighlightTrack expands to show wrapped tracks + highlights
axes = plot_tracks([gtrack, ht], region=region)
```

### Multiple Highlight Collections

```python
ht1 = HighlightTrack(regions=regions_a, track_list=[atrack, grtrack],
                     fill="#90EE90", alpha=0.3)
ht2 = HighlightTrack(regions=regions_b, track_list=[dtrack],
                     fill="#FFB6C1", alpha=0.3)

axes = plot_tracks([gtrack, ht1, grtrack, ht2], region=region)
```

### Using GenomicInterval Objects

```python
from geneview.genometracks import GenomicInterval

highlights = [
    GenomicInterval("chr7", 26520000, 26540000),
    GenomicInterval("chr7", 26640000, 26660000),
]
ht = HighlightTrack(regions=highlights, track_list=[atrack])
```

### Per-Region Colors

Use a list of colors to assign different colors to each highlight region:

```python
ht = HighlightTrack(
    regions=highlights,
    track_list=[atrack],
    fill=["#FF9999", "#99FF99", "#9999FF"],  # one color per region
    col=["#CC0000", "#00CC00", "#0000CC"],
    alpha=0.3,
)
```

Colors are cycled if there are more regions than colors.

### Highlight Targeting

Highlights are only drawn on panels for tracks listed in `track_list`. Other tracks in the plot are not affected:

```python
# Only atrack and grtrack get highlighted; dtrack is unaffected
ht = HighlightTrack(regions=regions, track_list=[atrack, grtrack])
axes = plot_tracks([gtrack, ht, dtrack], region=region)
```

---

## OverlayTrack

Overlays multiple tracks on the same plot area, useful for comparing signals from different DataTracks.

### Constructor

```python
OverlayTrack(
    track_list,        # List of Track objects to overlay
    name=None,         # Track name (defaults to first track's name)
    display_params=None,
)
```

### Basic Usage

```python
from geneview.genometracks import OverlayTrack

# Create two data tracks with different signals
dtrack1 = DataTrack(data1, type="line", name="Sample A",
    display_params={"col": "#3C5488"})
dtrack2 = DataTrack(data2, type="line", name="Sample B",
    display_params={"col": "#E64B35", "alpha": 0.7})

# Overlay them on the same axes
otrack = OverlayTrack([dtrack1, dtrack2], name="Comparison")
axes = plot_tracks([gtrack, otrack], region=region)
```

---

## File I/O

The genome tracks module includes utilities for reading common genomic file formats.

### BED Files

```python
from geneview.genometracks import read_bed

df = read_bed("features.bed")
# Returns DataFrame with columns: chrom, start, end, name, score, strand, ...
```

### GFF/GTF Files

```python
from geneview.genometracks import read_gff

# Works with both GFF3 and GTF formats
df = read_gff("annotations.gtf")
# Returns DataFrame with: chrom, source, feature, start, end, score, strand,
#   frame, gene_id, transcript_id, gene_name, ...
```

GFF3 uses `key=value` attributes; GTF uses `key "value"` attributes. Both are handled automatically.

### bedGraph Files

```python
from geneview.genometracks import read_bedgraph

df = read_bedgraph("signal.bedgraph")
# Returns DataFrame with: chrom, start, end, value
```

### BigWig Files (optional)

Requires `pyBigWig`:

```python
from geneview.genometracks import read_bigwig

df = read_bigwig("signal.bw", region=GenomicInterval("chr7", 26500000, 26800000))
```

### BAM Coverage (optional)

Requires `pysam`:

```python
from geneview.genometracks import read_bam_coverage

df = read_bam_coverage("alignments.bam",
                        region=GenomicInterval("chr7", 26500000, 26800000))
```

### Auto-detect Format

```python
from geneview.genometracks import read_auto

# Detects format from file extension
df = read_auto("features.bed")    # calls read_bed
df = read_auto("annotations.gtf") # calls read_gff
df = read_auto("signal.bedgraph") # calls read_bedgraph
```

### Using File Paths with Track Constructors

All track constructors accept file paths directly:

```python
atrack = AnnotationTrack("features.bed", name="BED Features")
grtrack = GeneRegionTrack("gene_models.gtf", name="Genes")
dtrack = DataTrack("signal.bedgraph", type="line", name="Signal")
```

---

## Advanced: plot_tracks()

### Function Signature

```python
plot_tracks(
    track_list,          # List of Track objects
    region=None,         # GenomicInterval (auto-detected if None)
    sizes=None,          # Relative heights for each track
    title=None,          # Main figure title
    figsize=None,        # Figure size (width, height) in inches
    extend_left=0,       # Extend left boundary (int bp or float fraction)
    extend_right=0,      # Extend right boundary
    show_title=True,     # Show track title panels
    reverse_strand=False, # Flip x-axis (3' on left)
    ax=None,             # Existing axes to plot into (single track only)
    **kwargs,            # Additional kwargs
)
# Returns: list of matplotlib Axes
```

### Custom Track Sizes

Control the relative height of each track:

```python
axes = plot_tracks(
    [gtrack, atrack, grtrack, dtrack],
    region=region,
    sizes=[0.2, 1.0, 1.5, 2.0],  # relative heights
)
```

### Custom Figure Size

```python
axes = plot_tracks([gtrack, atrack], region=region, figsize=(16, 6))
```

### Main Title

```python
axes = plot_tracks([gtrack, atrack], region=region, title="My Genomic Region")
```

### Region Extension

```python
# Extend by absolute base pairs
axes = plot_tracks([gtrack, atrack], extend_left=50000, extend_right=50000)

# Extend by fraction of data range (float between -1 and 1)
axes = plot_tracks([gtrack, atrack], extend_left=0.1, extend_right=0.1)
```

### Toggle Title Panel

```python
# Hide track title panels for a more compact layout
axes = plot_tracks([gtrack, atrack], show_title=False)
```

### Reverse Strand

Flip the x-axis so 3' is on the left:

```python
axes = plot_tracks([gtrack, atrack], reverse_strand=True)
```

### Shared X-axis

All tracks automatically share the same genomic x-axis. The last non-axis track shows coordinate labels; intermediate tracks hide their x-axis labels for a clean appearance.

### Getting the Figure

```python
axes = plot_tracks([gtrack, atrack], region=region)
fig = axes[0].figure  # Get the matplotlib Figure
fig.savefig("output.png", dpi=300, bbox_inches="tight")
```

---

## Display Parameters Reference

### GenomeAxisTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#555555"` | Axis line and tick color |
| `col_range` | `"#8B8378"` | Highlight region border |
| `fill_range` | `"#CDC8B1"` | Highlight region fill |
| `fontcolor` | `"#555555"` | Label text color |
| `fontsize` | 9 | Label font size |
| `lwd` | 2 | Axis line width |

### AnnotationTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#3C5488"` | Border color |
| `fill` | `"#5B8DB8"` | Fill color |
| `col_border` | `"#333333"` | Feature border color (uses `col` value) |
| `fontsize` | 10 | Label font size |
| `lwd` | 1.0 | Border line width |
| `min_width` | 1 | Minimum feature width (px) |
| `min_distance` | 1 | Minimum gap for stacking |
| `arrowHeadWidth` | 0.3 | Head width for fixedArrow shape |

### GeneRegionTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#3C5488"` | Exon border color |
| `fill` | `"#3C5488"` | CDS exon fill color |
| `fill_utr` | `"#8DB4E2"` | UTR fill color |
| `col_intron` | `"#888888"` | Intron line color |
| `fontsize` | 8 | Label font size |
| `fontcolor` | `"#333333"` | Label text color |
| `lwd` | 0.8 | Line width |

### DataTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#3C5488"` | Line color |
| `fill` | `"#5B8DB8"` | Fill color |
| `baseline` | 0 | Baseline y-position |
| `ncolor` | 100 | Gradient color count |
| `grid` | False | Draw horizontal grid lines |
| `col_grid` | `"#DDDDDD"` | Grid line color |

### Common (All Tracks)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 1.0 | Opacity |
| `background_panel` | `"white"` | Data panel background |
| `background_title` | `"#E8E8E8"` | Title panel background |
| `show_title` | True | Show title panel |
| `fontsize_title` | 10 | Title font size |
| `fontface_title` | `"bold"` | Title font weight |
| `col_title` | `"#333333"` | Title text color |
| `grid` | False | Show grid lines |
| `col_grid` | `"#DDDDDD"` | Grid line color |
| `frame` | False | Draw frame around panel |

---

## Complete Example

A full reproducible example combining all track types:

```python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
    HighlightTrack, OverlayTrack, GenomicInterval, plot_tracks,
    read_bed, read_gff, read_bedgraph,
)

# --- Load data ---
DATA_DIR = "examples/data/genome_tracks"
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
cov_data = read_bedgraph(os.path.join(DATA_DIR, "coverage.bedgraph"))
ann_data = read_bed(os.path.join(DATA_DIR, "annotations.bed"))

# --- Define region ---
region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# --- Create tracks ---
gtrack = GenomeAxisTrack(little_ticks=True)

atrack_cpg = AnnotationTrack(cpg_data, name="CpG Islands",
    display_params={"fill": "#3C5488"})

atrack_ann = AnnotationTrack(ann_data, name="Regulatory",
    shape="ellipse",
    display_params={"fill": "#E64B35"})

grtrack = GeneRegionTrack(gene_data, name="Gene Models",
    collapse_transcripts="longest")

dtrack = DataTrack(cov_data, type="histogram", name="Coverage",
    display_params={"fill": "#4DBBD5", "col": "#4DBBD5"})

# --- Highlight regions ---
ht = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7", "chr7"],
        "start": [26_505_000, 26_600_000],
        "end":   [26_535_000, 26_665_000],
    }),
    track_list=[atrack_cpg, atrack_ann, grtrack, dtrack],
    fill="#FFF3BF", alpha=0.3,
)

# --- Plot ---
axes = plot_tracks(
    [gtrack, ht],
    region=region,
    sizes=[0.2, 0.8, 0.8, 1.5, 1.5],
    title="Gene Region Overview",
    figsize=(16, 10),
    extend_left=0.05,
    extend_right=0.05,
)

fig = axes[0].figure
fig.savefig("genome_tracks_complete.png", dpi=300, bbox_inches="tight")
plt.show()
```

This produces a multi-panel figure showing:
- A genome axis with minor ticks at the top
- CpG island annotations (blue boxes)
- Regulatory features (red ellipses)
- Gene models with exons, UTRs, and introns
- Coverage histogram
- Yellow highlighted regions spanning all tracks
