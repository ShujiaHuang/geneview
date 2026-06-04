# Genome Tracks Guide

The genome tracks module in geneview provides a powerful system for visualizing genomic features along a shared coordinate axis. Inspired by the R/Bioconductor [Gviz](https://bioconductor.org/packages/Gviz) package, it has been redesigned for Python with matplotlib rendering and pandas DataFrames for data handling.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Basic Features](#basic-features)
3. [Display Parameters](#display-parameters)
4. [GenomeAxisTrack](#genomeaxistrack)
5. [IdeogramTrack](#ideogramtrack)
6. [AnnotationTrack](#annotationtrack)
7. [DetailsAnnotationTrack](#detailsannotationtrack)
8. [GeneRegionTrack](#generegiontrack)
9. [DataTrack](#datatrack)
10. [SequenceTrack](#sequencetrack)
11. [AlignmentsTrack](#alignmentstrack)
12. [HighlightTrack](#highlighttrack)
13. [OverlayTrack](#overlaytrack)
14. [File I/O](#file-io)
15. [Exporting Tracks](#exporting-tracks)
16. [Color Schemes](#color-schemes)
17. [Advanced: plot_tracks()](#advanced-plot_tracks)
18. [Display Parameters Reference](#display-parameters-reference)
19. [Complete Example](#complete-example)

---

## Introduction

### Core Concepts

The fundamental concept is similar to genome browsers: individual types of genomic features are represented by separate **track objects**. All tracks share a common genomic coordinate system, so alignment is handled automatically.

```
GenomeAxisTrack  ─── coordinate ruler
IdeogramTrack    ─── chromosome ideogram (cytobands)
AnnotationTrack  ─── generic ranges (boxes/arrows)
  └── DetailsAnnotationTrack  ─── annotations with detail panels
GeneRegionTrack  ─── gene models (exons/UTRs/introns)
DataTrack        ─── numeric data (line/histogram/heatmap)
SequenceTrack    ─── nucleotide sequence display
AlignmentsTrack  ─── BAM/CRAM read alignments
HighlightTrack   ─── cross-track highlights
OverlayTrack     ─── overlay multiple tracks on same axes
```

### Track Hierarchy

```
Track (abstract base)
  ├── GenomeAxisTrack
  ├── IdeogramTrack (chromosome schematic)
  ├── SequenceTrack (nucleotide sequence)
  └── RangeTrack (has a DataFrame)
        ├── StackedTrack (overlapping features stacked)
        │     ├── AnnotationTrack
        │     │     └── DetailsAnnotationTrack
        │     ├── GeneRegionTrack
        │     └── AlignmentsTrack (BAM/CRAM reads)
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

Every track has a `display_params` dictionary controlling its appearance. Parameters can be set at construction time or modified later. Following the Gviz convention, all extra keyword arguments passed to a track constructor are automatically treated as display parameters.

### Setting Parameters at Construction

```python
# Using display_params dict
atrack = AnnotationTrack(data, name="Features", display_params={
    "fill": "lightblue",
    "col": "#333333",
    "fontsize": 9,
})

# Using keyword arguments (Gviz-style, recommended)
atrack = AnnotationTrack(data, name="Features",
    fill="lightblue", col="#333333", fontsize=9)
```

### Getting and Setting Parameters

```python
# Get a parameter
color = atrack.get_param("fill")

# Set a parameter
atrack.set_param("fill", "#DC0000")

# Set multiple parameters at once
atrack.set_params({"fill": "#DC0000", "col": "white"})

# Get all display parameters
all_params = atrack.display_params()
```

### Alias System

Display parameters support both snake_case (Python) and dot-notation (Gviz-style) names:

```python
# These are equivalent:
atrack.set_param("col_title", "white")
atrack.set_param("col.title", "white")

# Query available defaults for any track class:
from geneview.genometracks import available_display_params
base_params = available_display_params()
axis_params = available_display_params("GenomeAxisTrack")
```

### Common Display Parameters (All Tracks)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 1.0 | Opacity of plot elements |
| `background_panel` | `"transparent"` | Background color of data panel |
| `background_title` | `"#D3D3D3"` | Background color of title panel |
| `col` | `"#0080FF"` | Line/border color |
| `fill` | `"lightgray"` | Fill color |
| `col_border` | `"transparent"` | Border color |
| `fontcolor` | `"black"` | Text color |
| `fontsize` | 12 | Base font size |
| `fontsize_title` | 12 | Title font size |
| `lwd` | 1.0 | Line width |
| `lty` | `"solid"` | Line type |
| `show_title` | True | Show the title panel |
| `show_axis` | True | Show axis elements |
| `cex` | 1.0 | Font expansion factor |
| `min_width` | 1 | Minimum feature width in pixels |
| `min_distance` | 1 | Minimum distance for stacking |
| `collapse` | True | Collapse overlapping features |
| `rotation` | 0 | Text rotation in degrees |
| `rotation_title` | 90 | Title text rotation |

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
    exponent=None,         # Force label exponent (3=kb, 6=Mb, 9=Gb)
    ticks_at=None,         # Explicit tick positions (list of bp)
    display_params=None,
    **kwargs,              # Additional display parameters
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

## IdeogramTrack

Displays a chromosome ideogram with cytogenetic bands (cytobands), centromere, and optional region highlighting. Ported from Gviz's IdeogramTrack.

### Constructor

```python
IdeogramTrack(
    bands,                  # DataFrame with cytoband data or file path
    chromosome=None,        # Chromosome to display (e.g., "chr7")
    show_id=True,           # Show chromosome name
    show_band_id=False,     # Show band labels (p11.1, q21, etc.)
    centromere_shape="triangle",  # "triangle" or "circle"
    outline=False,          # Draw outline around ideogram
    name=None,              # Track name (auto-generated if None)
    height=0.5,             # Relative height (compact)
    display_params=None,
)
```

### Basic Usage

```python
import pandas as pd
from geneview.genometracks import IdeogramTrack, GenomicInterval, plot_tracks

# Cytoband data: chrom, start, end, band name, stain type
bands = pd.DataFrame({
    "chrom": ["chr7"] * 5,
    "chromStart": [0, 2000000, 5000000, 8000000, 10000000],
    "chromEnd": [2000000, 5000000, 8000000, 10000000, 12000000],
    "name": ["p25", "p21", "p11", "q11", "q21"],
    "gieStain": ["gneg", "gpos25", "acen", "acen", "gpos50"],
})

itrack = IdeogramTrack(bands, chromosome="chr7")
region = GenomicInterval("chr7", 3000000, 9000000)
axes = plot_tracks([itrack], region=region)
```

### Cytoband Data Format

The `bands` DataFrame must have these columns:
- `chrom`: Chromosome name (e.g., "chr7")
- `chromStart`: Start position
- `chromEnd`: End position
- `name`: Band name (e.g., "p21.1", "q31.2")
- `gieStain`: Stain type (determines color)

Supported `gieStain` values and their colors:

| Stain | Color | Description |
|-------|-------|-------------|
| `gneg` | White | Negative (light) band |
| `gpos25` | Light gray | 25% positive |
| `gpos33` | Light gray | 33% positive |
| `gpos50` | Medium gray | 50% positive |
| `gpos66` | Dark gray | 66% positive |
| `gpos75` | Dark gray | 75% positive |
| `gpos100` | Black | 100% positive (darkest) |
| `acen` | Red | Centromere region |
| `stalk` | Blue-gray | Acrocentric stalk |
| `gvar` | Light gray | Variable region |

### Centromere Shapes

```python
# Triangle centromere (default, pinched waist)
itrack = IdeogramTrack(bands, chromosome="chr7", centromere_shape="triangle")

# Circle centromere
itrack = IdeogramTrack(bands, chromosome="chr7", centromere_shape="circle")
```

### Band Labels

```python
# Show band names (p21, q31, etc.)
itrack = IdeogramTrack(bands, chromosome="chr7", show_band_id=True)
```

### Region Highlighting

When plotted with a `region`, the IdeogramTrack automatically draws a highlight box showing the current view:

```python
# The current region is highlighted with a pink box
region = GenomicInterval("chr7", 5000000, 8000000)
axes = plot_tracks([itrack], region=region)
```

The highlight box colors can be customized:

```python
itrack = IdeogramTrack(bands, chromosome="chr7", display_params={
    "col": "blue",        # Highlight box border
    "fill": "#FFFFCC",    # Highlight box fill
})
```

### Loading from File

```python
# Load cytoband data from a tab-separated file
itrack = IdeogramTrack("cytobands.txt", chromosome="chr7")
```

### Display Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"red"` | Highlight box border color |
| `fill` | `"#FFE3E6"` | Highlight box fill color (pink) |
| `lwd` | 1.0 | Highlight box line width |
| `fontsize` | 10 | Band label font size |
| `fontcolor` | `"#808080"` | Band label color |
| `show_title` | False | Hide title panel by default |

### Integration with Other Tracks

The IdeogramTrack shows the full chromosome while other tracks zoom into a region. The current region is highlighted on the ideogram:

```python
# IdeogramTrack shows full chromosome
itrack = IdeogramTrack(bands, chromosome="chr7")

# Other tracks show the zoomed region
gtrack = GenomeAxisTrack()
atrack = AnnotationTrack(data, name="Features")

# All tracks plotted together
region = GenomicInterval("chr7", 5000000, 8000000)
axes = plot_tracks([itrack, gtrack, atrack], region=region)
```

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
    "fill": "lightblue",          # Default fill
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

### Group Label Justification

Control where group labels appear relative to the feature row:

```python
atrack = AnnotationTrack(data, group_annotation="feature",
    just_group="left")    # left, right, above, below
```

### Overplotting Indicator

Color-code features by overlap depth when collapsing:

```python
atrack = AnnotationTrack(data, show_overplotting=True)
```

### Merging Groups

When `collapse=True`, merge overlapping items within the same group:

```python
atrack = AnnotationTrack(data, merge_groups=True, stacking="squish")
```

### Connector Line Color

Customize the line connecting grouped features:

```python
atrack = AnnotationTrack(data, group_annotation="id",
    col_line="#CC0000")
```

### Creating from BAM Files

Create an AnnotationTrack from aligned reads in a BAM file (requires `pysam`):

```python
atrack = AnnotationTrack.from_bam(
    "alignments.bam",
    region=GenomicInterval("chr7", 26500000, 26800000),
    name="Reads",
)
```

Each aligned read becomes a separate annotation feature.

---

## DetailsAnnotationTrack

An extended AnnotationTrack that draws detail panels below selected features, connected by lines. Useful for showing per-feature zoom-in views or custom sub-plots.

### Constructor

```python
DetailsAnnotationTrack(
    data,                      # DataFrame or file path
    fun=None,                  # Callable(identifier, ax, data_row) for custom detail
    select_fun=None,           # Callable(data_row) -> bool; which features get details
    details_size=0.3,          # Fraction of track height for detail panels
    details_connector_col="#888888",
    details_connector_lty="-",
    details_connector_lwd=0.5,
    details_border_col="#CCCCCC",
    details_border_fill="#FFFFEE",
    details_min_width=20,      # Minimum pixel width for a detail panel
    details_ratio=1.0,         # Width ratio relative to feature
    group_details=False,       # Group detail panels by feature group
    name="DetailsAnnotation",
    height=2.0,
    **kwargs,                  # Passed to AnnotationTrack
)
```

### Basic Usage

```python
from geneview.genometracks import DetailsAnnotationTrack, GenomicInterval, plot_tracks

data = pd.DataFrame({
    "chrom": ["chr7"] * 3,
    "start": [2000000, 2050000, 2100000],
    "end":   [2020000, 2070000, 2120000],
    "name": ["geneA", "geneB", "geneC"],
})

# Default detail panels (shows feature info text)
dtrack = DetailsAnnotationTrack(data, name="Details")
axes = plot_tracks([dtrack], region=GenomicInterval("chr7", 1950000, 2150000))
```

### Custom Detail Function

```python
def my_detail(identifier, ax, data_row):
    ax.text(0.5, 0.5, f"Detail for {identifier}",
            transform=ax.transAxes, ha="center", va="center")
    ax.set_title(identifier, fontsize=6)

dtrack = DetailsAnnotationTrack(data, fun=my_detail, details_size=0.4)
```

### Selective Details

Only show detail panels for specific features:

```python
dtrack = DetailsAnnotationTrack(
    data,
    select_fun=lambda row: row.get("name", "").startswith("geneA"),
)
```

---

## GeneRegionTrack

Displays gene models with exons as thick boxes, UTRs as thin boxes, and introns as connecting lines.

### GeneRegionTrack Constructor

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

### Basic Usage of genometracks

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

### Display Parameters of GeneRegionTrack

| Parameter | Default | Description |
| :-----------: | :---------: | :-------------: |
| `fill` | `"orange"` | CDS exon fill color |
| `fill_utr` | `"#FFD699"` | UTR fill color |
| `col_intron` | `"#808080"` | Intron line color |
| `fontsize` | 8 | Label font size |
| `fontcolor` | `"#333333"` | Label text color |
| `lwd` | 0.8 | Line width |

### Exon Annotation Labels

Label individual exons with different information modes:

```python
# Show exon numbers on each exon box
grtrack = GeneRegionTrack(data, exon_annotation="exon")

# Show gene symbol on each exon
grtrack = GeneRegionTrack(data, exon_annotation="symbol")

# Show gene ID
grtrack = GeneRegionTrack(data, exon_annotation="gene")

# Show transcript ID
grtrack = GeneRegionTrack(data, exon_annotation="transcript")

# Show feature type
grtrack = GeneRegionTrack(data, exon_annotation="feature")
```

### Gene Symbols

Use gene_name (symbol) instead of gene_id for transcript labels:

```python
grtrack = GeneRegionTrack(data, gene_symbols=True)
```

### Transcript Annotation Alias

`transcript_annotation` is an alias for `show_id` with additional `"symbol"` mode:

```python
# Equivalent to show_id="transcript"
grtrack = GeneRegionTrack(data, transcript_annotation="transcript")

# Show gene symbol instead of ID
grtrack = GeneRegionTrack(data, transcript_annotation="symbol")
```

---

## DataTrack

Visualizes numeric genomic data along coordinates. Supports multiple plot types and multi-sample data.

### DataTrack Constructor

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
    show_sample_names=False,  # Show sample names on heatmap (Gviz default: FALSE)
    separator=0,       # Separator line width between heatmap rows (0=none)
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
| `"a"` | Average (mean across value columns per position) |
| `"confint"` | Confidence interval (mean ± 1.96×SE band) |
| `"smooth"` | LOWESS/loess smoothed line |
| `"horizon"` | Horizon plot (positive/negative bands from baseline) |
| `"g"` | Horizontal grid lines (standalone) |
| `"r"` | Linear regression line |

Multiple types can be combined by passing a list:

```python
dtrack = DataTrack(data, type=["boxplot", "a", "g"])
```

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

The heatmap type works best with multiple value columns. It uses a sequential blue gradient (matching Gviz's `colorRampPalette(brewer.pal(9, "Blues"))`):

```python
data = pd.DataFrame({
    "chrom": ["chr7"] * 50,
    "start": np.arange(26500000, 26800000, 6000),
    "end": np.arange(26506000, 26806000, 6000),
    "sample_A": np.random.randn(50).cumsum() / 3,
    "sample_B": np.random.randn(50).cumsum() / 3 + 2,
    "sample_C": np.random.randn(50).cumsum() / 3 - 1,
})

# Basic heatmap
dtrack = DataTrack(data, type="heatmap", name="Heatmap")

# With sample names and row separators (matching Gviz style)
dtrack = DataTrack(
    data,
    type="heatmap",
    name="Heatmap",
    show_sample_names=True,  # Show sample names on y-axis
    separator=2,             # Add white lines between rows
)
```

The heatmap uses a **sequential blue gradient** (white to dark blue) by default, matching Gviz's default color scheme. The gradient ranges from white (low values) to dark blue `#08306B` (high values).

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
                   display_params={"group_colors": ["#0080FF", "#DC0000"]})
```

### Grid Lines

Draw horizontal grid lines:

```python
dtrack = DataTrack(data, type="line",
                   display_params={"grid": True, "col_grid": "#DDDDDD"})
```

### Average (type "a")

Compute and draw the mean across all value columns at each position:

```python
# Multi-sample data
data = pd.DataFrame({
    "chrom": ["chr7"] * 50,
    "start": np.arange(1000, 2000, 20),
    "end":   np.arange(1020, 2020, 20),
    "sample_A": rng.randn(50).cumsum(),
    "sample_B": rng.randn(50).cumsum() + 2,
    "sample_C": rng.randn(50).cumsum() - 1,
})
dtrack = DataTrack(data, type="a", name="Average", col="#3C5488")
```

### Confidence Interval (type "confint")

Draw mean ± 1.96×SE as a filled band:

```python
dtrack = DataTrack(data, type="confint", name="95% CI", col="#00A087")
```

### Smooth / Loess (type "smooth")

Apply LOWESS smoothing (or rolling mean fallback):

```python
dtrack = DataTrack(data, type="smooth", smooth_span=0.3, name="Smooth")
```

The `smooth_span` parameter controls the fraction of data used for local regression (default 0.3). Requires `statsmodels` for LOWESS; falls back to a rolling mean if not available.

### Horizon Plot (type "horizon")

Split data into positive and negative bands relative to a baseline:

```python
dtrack = DataTrack(data, type="horizon", name="Horizon", baseline=0)
```

Positive values are drawn above the baseline in blue; negative values below in red.

### Grid Type (type "g")

Standalone horizontal grid lines:

```python
dtrack = DataTrack(data, type="g", name="Grid")
```

### Regression (type "r")

Fit and draw a linear regression line:

```python
dtrack = DataTrack(data, type="r", name="Regression")
```

### Composite Plot Types

Combine multiple plot types by passing a list:

```python
# Boxplot + average line + grid
dtrack = DataTrack(data, type=["boxplot", "a", "g"], name="Composite")

# Line + confidence interval
dtrack = DataTrack(data, type=["line", "confint"], name="Signal + CI")
```

Each type is drawn in order on the same axes.

### Display Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fill` | `"#808080"` | Fill color for histogram/polygon/mountain |
| `col` | `"#0080FF"` | Line color |
| `baseline` | 0 | Y-position of baseline for mountain/polygon |
| `ncolor` | 100 | Number of colors for gradient |

---

## SequenceTrack

Displays nucleotide sequence as colored letters, boxes, or a continuous line depending on zoom level.

### Constructor

```python
SequenceTrack(
    sequence=None,        # str, DataFrame with chrom/start/end/seq, or None
    fasta_path=None,      # Path to indexed FASTA file (.fa + .fai)
    twobit_path=None,     # Path to 2bit genome file
    chromosome=None,      # Chromosome name
    complement=False,     # Show reverse complement
    add53=False,          # Draw 5'→3' direction arrow
    fontcolor=None,       # Dict of nucleotide -> color
    noLetters=False,      # Force box mode even at high zoom
    cex=1.0,              # Font size scaling factor
    name="Sequence",
    height=0.5,
)
```

### Zoom-Level Rendering

The display automatically adapts to the visible region width:

| Region Width | Display Mode |
|--------------|-------------|
| < 200 bp | Individual colored letters (A, C, G, T) |
| 200–2000 bp | Colored boxes per nucleotide |
| > 2000 bp | Single continuous colored line |

Default nucleotide colors: A=green (`#009E73`), C=blue (`#0072B2`), G=yellow (`#E69F00`), T=red (`#D55E00`).

### String Sequence

```python
from geneview.genometracks import SequenceTrack, GenomicInterval, plot_tracks

seq_track = SequenceTrack(
    sequence="ATCGATCGATCGATCG" * 5,
    chromosome="chr7",
    name="Sequence",
)
axes = plot_tracks([seq_track], region=GenomicInterval("chr7", 0, 80))
```

### Loading from FASTA

```python
seq_track = SequenceTrack(
    fasta_path="hg38.fa",
    chromosome="chr7",
    name="Sequence",
)
axes = plot_tracks([seq_track], region=GenomicInterval("chr7", 26500000, 26500100))
```

### Reverse Complement

```python
seq_track = SequenceTrack(sequence="ATCGATCG", complement=True, add53=True)
```

### Custom Nucleotide Colors

```python
seq_track = SequenceTrack(
    sequence="ATCGATCG",
    fontcolor={"A": "#FF0000", "T": "#00FF00", "C": "#0000FF", "G": "#FFFF00"},
)
```

---

## AlignmentsTrack

Visualizes aligned reads from BAM or CRAM files with coverage, pileup, and sashimi plot modes. Requires `pysam`.

### Constructor

```python
AlignmentsTrack(
    filepath,              # Path to BAM or CRAM file
    is_paired=False,       # Paired-end reads
    show_mismatches=False, # Color mismatched bases vs reference
    show_indels=False,     # Show insertions/deletions
    reference=None,        # FASTA path for mismatch detection
    type="coverage",       # "coverage", "pileup", "sashimi", or list
    coverage_height=0.3,   # Relative height of coverage panel
    sashimi_height=0.5,    # Relative height of sashimi panel
    sashimi_score=1,       # Minimum reads for sashimi arcs
    sashimi_filter=None,   # DataFrame of junction regions to show
    reverse_stacking=False, # Reverse stacking order
    col_mates=None,        # Color for mate connectors
    col_gap="#888888",     # Color for gap/junction lines
    col_deletion="#FF0000", # Deletion marker color
    col_insertion="#0000FF", # Insertion marker color
    fill_coverage="#B3CDE3", # Coverage fill color
    fill_reads="#C8C8C8",   # Read fill color
    alpha_reads=0.7,         # Read transparency
    name="Alignments",
    height=2.0,
)
```

### Plot Modes

| Mode | Description |
|------|-------------|
| `coverage` | Coverage depth histogram |
| `pileup` | Individual reads with CIGAR-aware rendering |
| `sashimi` | Arc plot showing splice junctions with read counts |

Multiple modes can be combined:

```python
# Coverage + pileup
aln = AlignmentsTrack("sample.bam", type=["coverage", "pileup"])

# Coverage + sashimi
aln = AlignmentsTrack("sample.bam", type=["coverage", "sashimi"])
```

### Basic Usage

```python
from geneview.genometracks import (
    AlignmentsTrack, GeneRegionTrack, GenomeAxisTrack,
    GenomicInterval, plot_tracks,
)

region = GenomicInterval("chr7", 26500000, 26550000)

aln = AlignmentsTrack(
    "alignments.bam",
    type="coverage",
    name="Coverage",
)
grtrack = GeneRegionTrack("gene_models.gtf", name="Genes")

axes = plot_tracks(
    [GenomeAxisTrack(), aln, grtrack],
    region=region,
)
```

### Mismatch Display

Color-code mismatched bases against a reference sequence:

```python
aln = AlignmentsTrack(
    "alignments.bam",
    type="pileup",
    show_mismatches=True,
    reference="hg38.fa",
)
```

### Indel Display

Show insertions as vertical bars and deletions as bridging lines:

```python
aln = AlignmentsTrack(
    "alignments.bam",
    type="pileup",
    show_indels=True,
)
```

### Paired-End Reads

Draw connectors between read mates:

```python
aln = AlignmentsTrack(
    "paired.bam",
    type="pileup",
    is_paired=True,
    col_mates="#CC0000",
)
```

### Sashimi Plots

Visualize splice junctions as arcs with read count labels:

```python
aln = AlignmentsTrack(
    "rnaseq.bam",
    type=["coverage", "sashimi"],
    sashimi_score=5,  # Only show junctions with >= 5 reads
)
```

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
    display_params={"col": "#0080FF"})
dtrack2 = DataTrack(data2, type="line", name="Sample B",
    display_params={"col": "#DC0000", "alpha": 0.7})

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
from geneview.genometracks import read_bam_coverage, GenomicInterval

df = read_bam_coverage("alignments.bam",
                        region=GenomicInterval("chr7", 26500000, 26800000))
```

### CRAM Coverage (optional)

CRAM is the reference-based compressed counterpart of BAM. Requires `pysam` and typically a reference FASTA file (indexed with `samtools faidx`):

```python
from geneview.genometracks import read_cram_coverage, GenomicInterval

region = GenomicInterval("chr7", 26500000, 26800000)
df = read_cram_coverage("alignments.cram", region=region, reference="hg38.fa")
```

The `reference` parameter can be omitted if the reference path is embedded in the CRAM header and accessible on the local filesystem.

You can also feed BAM or CRAM coverage directly into a DataTrack:

```python
from geneview.genometracks import DataTrack, read_bam_coverage, read_cram_coverage, plot_tracks

region = GenomicInterval("chr7", 26500000, 26800000)
bam_cov = read_bam_coverage("sample.bam", region=region)
cram_cov = read_cram_coverage("sample.cram", region=region, reference="hg38.fa")

bam_track = DataTrack(bam_cov, type="histogram", name="BAM Coverage")
cram_track = DataTrack(cram_cov, type="histogram", name="CRAM Coverage")

axes = plot_tracks([bam_track, cram_track], region=region)
```

### Auto-detect Format

```python
from geneview.genometracks import read_auto

# Detects format from file extension
df = read_auto("features.bed")    # calls read_bed
df = read_auto("annotations.gtf") # calls read_gff
df = read_auto("signal.bedgraph") # calls read_bedgraph
df = read_auto("signal.bw")       # calls read_bigwig
df = read_auto("alignments.bam")  # calls read_bam_coverage
df = read_auto("alignments.cram") # calls read_cram_coverage
df = read_auto("signal.wig")      # calls read_wig
seq = read_auto("genome.fa")      # calls read_fasta
seq = read_auto("genome.2bit")    # calls read_2bit
```

### WIG Files

Parse WIG format files (both fixedStep and variableStep):

```python
from geneview.genometracks import read_wig

df = read_wig("signal.wig")
# Returns DataFrame with: chrom, start, end, value
```

### FASTA Files (optional)

Read nucleotide sequences from indexed FASTA files. Requires `pysam`:

```python
from geneview.genometracks import read_fasta, GenomicInterval

seq = read_fasta("hg38.fa", region=GenomicInterval("chr7", 26500000, 26500100))
# Returns a nucleotide string
```

### 2bit Files (optional)

Read nucleotide sequences from 2bit genome files. Requires `py2bit`:

```python
from geneview.genometracks import read_2bit, GenomicInterval

seq = read_2bit("hg38.2bit", region=GenomicInterval("chr7", 26500000, 26500100))
# Returns a nucleotide string
```

### Using File Paths with Track Constructors

All track constructors accept file paths directly:

```python
atrack = AnnotationTrack("features.bed", name="BED Features")
grtrack = GeneRegionTrack("gene_models.gtf", name="Genes")
dtrack = DataTrack("signal.bedgraph", type="line", name="Signal")
```

### Loading Included Test Data

The `examples/data/genome_tracks/` directory ships real test files that you can use
to try each reader:

| File | Format | Region | Description |
|------|--------|--------|-------------|
| `test.bed` | BED9 | chr7:127.47M | 9 features with RGB colors |
| `test.bedGraph` | bedGraph | chr19:49.3M | 9 bins, signed values (-1 to +1) |
| `test.bw` | BigWig | chr19:49.3M | BigWig coverage |
| `test.bam` | BAM | chr1:189.89M | 3118 reads (hg19) |
| `test.gtf` | GTF | chr1:67M | SGIP1 gene model |
| `test.gff3` | GFF3 | chr1:67M | Same gene, GFF3 format |
| `cpg_islands.bed` | BED6 | chr7:26.5M | 9 synthetic CpG islands |
| `gene_models.gtf` | GTF | chr7:26.5M | 3 synthetic gene models |
| `coverage.bedgraph` | bedGraph | chr7:26.5M | 100-bin synthetic coverage |
| `annotations.bed` | BED6 | chr7:26.5M | 8 synthetic regulatory features |
| `multi_sample.tsv` | TSV | chr7:26.5M | 6 samples (3 ctrl + 3 treat) |

```python
from geneview.genometracks import (
    read_bed, read_bedgraph, read_gff, read_bigwig,
    read_bam_coverage, read_auto, GenomicInterval,
)

# Real test data examples
bed_data  = read_bed("examples/data/genome_tracks/test.bed")
bg_data   = read_bedgraph("examples/data/genome_tracks/test.bedGraph")
bigwig    = read_bigwig("examples/data/genome_tracks/test.bw")       # requires pyBigWig
bam_cov   = read_bam_coverage(                                        # requires pysam
    "examples/data/genome_tracks/test.bam",
    region=GenomicInterval("chr1", 189_891_000, 189_900_000),
)
gtf_data  = read_gff("examples/data/genome_tracks/test.gtf")
gff3_data = read_gff("examples/data/genome_tracks/test.gff3")

# Auto-detect format
auto_bed  = read_auto("examples/data/genome_tracks/test.bed")
auto_bw   = read_auto("examples/data/genome_tracks/test.bw")
```

---

## Exporting Tracks

Export track data to common genomic file formats using `export_tracks`:

```python
from geneview.genometracks import export_tracks

# Export AnnotationTrack to BED
export_tracks(atrack, "features.bed", format="bed")

# Export DataTrack to bedGraph
export_tracks(dtrack, "signal.bedgraph", format="bedgraph")

# Export DataTrack to WIG
export_tracks(dtrack, "signal.wig", format="wig")

# Export GeneRegionTrack to GFF
export_tracks(grtrack, "genes.gff", format="gff")
```

| Format | Output Columns |
|--------|---------------|
| `bed` | chrom, start, end, name, score, strand |
| `gff` | chrom, source, feature, start (1-based), end, score, strand, frame, attributes |
| `bedgraph` | chrom, start, end, value |
| `wig` | fixedStep or variableStep WIG format |

---

## Color Schemes

Apply predefined color schemes to tracks using `apply_scheme` or the `scheme` parameter in `plot_tracks`:

```python
from geneview.genometracks import apply_scheme, plot_tracks

# Apply a scheme to a single track
apply_scheme(grtrack, "genes")

# Or apply via plot_tracks
axes = plot_tracks([grtrack], region=region, scheme="transcripts")
```

### Available Schemes

| Scheme | Description |
|--------|-------------|
| `"default"` | Light gray fill, dark gray border |
| `"genes"` | Distinct color per gene |
| `"transcripts"` | Distinct color per transcript |

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
    cex=None,            # Global font expansion factor
    add=False,           # Plot into existing axes (no new figure)
    ylim=None,           # Global y-axis limits for DataTrack panels
    scheme=None,         # Apply named color scheme ("default", "genes", "transcripts")
    panel_only=False,    # Data panels only, no title/new figure (implies add=True)
    margin=6,            # Pixel margin around the plot
    inner_margin=3,      # Inner margin between title and data panels (px)
    **kwargs,            # Additional display parameters for all tracks
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

### Font Scaling (cex)

Globally expand or shrink all text:

```python
axes = plot_tracks([gtrack, grtrack], cex=1.5)  # 150% font size
```

### Plot Into Existing Axes (add)

Add a track to an already-existing figure:

```python
fig, ax = plt.subplots(figsize=(10, 2))
axes = plot_tracks([dtrack], region=region, add=True, ax=ax)
```

### Global Y-axis Limits (ylim)

Set the same y-axis limits for all DataTrack panels:

```python
axes = plot_tracks([dtrack1, dtrack2], region=region, ylim=(-5, 10))
```

### Color Scheme

Apply a named color scheme to all annotation and gene tracks:

```python
axes = plot_tracks([grtrack, atrack], region=region, scheme="genes")
```

### Panel-Only Mode (Embedding)

Draw only the data panels without title panels or creating a new figure. Useful for embedding tracks in larger matplotlib layouts:

```python
fig, axes = plt.subplots(2, 1, figsize=(12, 6))
plot_tracks([gtrack], region=region, panel_only=True, ax=axes[0])
plot_tracks([dtrack], region=region, panel_only=True, ax=axes[1])
```

### Margin Controls

Control spacing around and between panels:

```python
# Adjust margins (in pixels, Gviz-compatible)
axes = plot_tracks([gtrack, atrack], region=region,
    margin=10,          # outer margin
    inner_margin=5,     # space between title and data panels
)
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

### IdeogramTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"red"` | Highlight box border color |
| `fill` | `"#FFE3E6"` | Highlight box fill color (pink) |
| `lwd` | 1.0 | Highlight box line width |
| `fontsize` | 10 | Band label font size |
| `fontcolor` | `"#808080"` | Band label color |
| `show_title` | False | Hide title panel by default |

### AnnotationTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"transparent"` | Border color (edges match fill, Gviz default) |
| `fill` | `"lightblue"` | Fill color |
| `col_border` | `"#333333"` | Feature border color (uses `col` value) |
| `fontsize` | 10 | Label font size |
| `lwd` | 1.0 | Border line width |
| `min_width` | 1 | Minimum feature width (px) |
| `min_distance` | 1 | Minimum gap for stacking |
| `arrowHeadWidth` | 0.3 | Head width for fixedArrow shape |

### GeneRegionTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"orange"` | Exon border color |
| `fill` | `"orange"` | CDS exon fill color |
| `fill_utr` | `"#FFD699"` | UTR fill color |
| `col_intron` | `"#888888"` | Intron line color |
| `fontsize` | 8 | Label font size |
| `fontcolor` | `"#333333"` | Label text color |
| `lwd` | 0.8 | Line width |

### DataTrack

| Parameter | Default | Description |
|-----------|---------|-------------|
| `col` | `"#0080FF"` | Line color |
| `fill` | `"#808080"` | Fill color |
| `col_histogram` | `"#808080"` | Histogram bar color |
| `col_mountain` | None (uses `col`) | Mountain plot line color |
| `fill_mountain` | None (uses `fill`) | Mountain plot fill color |
| `col_baseline` | `"#888888"` | Baseline line color |
| `alpha_confint` | 0.3 | Confidence interval transparency |
| `cex_legend` | 1.0 | Legend font size factor |
| `baseline` | 0 | Baseline y-position |
| `ncolor` | 100 | Gradient color count |
| `show_sample_names` | False | Show sample names on heatmap y-axis (Gviz default: FALSE) |
| `separator` | 0 | Separator line width between heatmap rows (0=none) |
| `grid` | False | Draw horizontal grid lines |
| `col_grid` | `"#DDDDDD"` | Grid line color |

**Note:** The heatmap type uses a sequential blue gradient (white to dark blue `#08306B`), matching Gviz's `colorRampPalette(brewer.pal(9, "Blues"))`.

### Common (All Tracks)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `alpha` | 1.0 | Opacity |
| `background_panel` | `"transparent"` | Data panel background |
| `background_title` | `"#D3D3D3"` | Title panel background |
| `show_title` | True | Show title panel |
| `fontsize_title` | 12 | Title font size |
| `fontface_title` | `"bold"` | Title font weight |
| `col_title` | `"white"` | Title text color |
| `grid` | False | Show grid lines |
| `col_grid` | `"#808080"` | Grid line color |
| `frame` | False | Draw frame around panel |
| `fontcolor` | `"black"` | Text color |
| `fontsize` | 12 | Base font size |
| `lty` | `"solid"` | Line type |

---

## BiomartGeneRegionTrack and UcscTrack (Stubs)

These classes provide Gviz-compatible API signatures but raise `NotImplementedError` when drawn, since they require external server connections.

### BiomartGeneRegionTrack

```python
from geneview.genometracks import BiomartGeneRegionTrack

# Accepts Gviz-compatible parameters
bm = BiomartGeneRegionTrack(
    genome="hg38",
    chromosome="chr7",
    start=26500000,
    end=26800000,
    symbol="EGFR",
)
# bm.draw() raises NotImplementedError
```

**Alternative:** Use `GeneRegionTrack` with a local GTF/GFF file:

```python
from geneview.genometracks import GeneRegionTrack, read_gff
gene_data = read_gff("Homo_sapiens.gtf")
grtrack = GeneRegionTrack(gene_data, name="Genes")
```

### UcscTrack

```python
from geneview.genometracks import UcscTrack

# Accepts Gviz-compatible parameters
ucsc = UcscTrack(
    genome="hg38",
    chromosome="chr7",
    track="knownGene",
    table="knownGene",
)
# ucsc.draw() raises NotImplementedError
```

**Alternative:** Download the file from UCSC and use the appropriate `read_*` function:

```python
from geneview.genometracks import read_bed, AnnotationTrack
df = read_bed("ucsc_known_genes.bed")
atrack = AnnotationTrack(df, name="UCSC Genes")
```

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
    display_params={"fill": "lightblue"})

atrack_ann = AnnotationTrack(ann_data, name="Regulatory",
    shape="ellipse",
    display_params={"fill": "#DC0000"})

grtrack = GeneRegionTrack(gene_data, name="Gene Models",
    collapse_transcripts="longest")

dtrack = DataTrack(cov_data, type="histogram", name="Coverage",
    display_params={"col": "#0080FF"})

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
