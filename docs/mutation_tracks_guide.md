# Mutation Tracks Guide

The **mutation tracks** module in geneview provides lollipop-style and dandelion-style visualizations for mutations, variants, and methylation sites along protein features. Ported from the R/Bioconductor [trackViewer](https://bioconductor.org/packages/trackViewer) package, these tracks are designed for publication-ready figures showing variant annotations with per-site customization.

---

## Table of Contents

- [Mutation Tracks Guide](#mutation-tracks-guide)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
    - [Core Concepts](#core-concepts)
    - [Data Format](#data-format)
  - [LolliplotTrack](#lolliplottrack)
    - [Basic Usage](#basic-usage)
    - [Shape Types](#shape-types)
    - [Per-SNP Customization](#per-snp-customization)
    - [Tanghulu Stacking](#tanghulu-stacking)
    - [Node Labels](#node-labels)
    - [Caterpillar Layout](#caterpillar-layout)
    - [Multiple Shapes in Tanghulu](#multiple-shapes-in-tanghulu)
    - [Legend](#legend)
    - [Aligned Labels (jitter="label")](#aligned-labels-jitterlabel)
    - [Custom Y-axis and ylab](#custom-y-axis-and-ylab)
    - [pie.stack Type](#piestack-type)
    - [Coordinate Rescaling](#coordinate-rescaling)
    - [Multi-layer Features](#multi-layer-features)
    - [Shape Types per SNP](#shape-types-per-snp)
  - [DandelionTrack](#dandeliontrack)
    - [DandelionTrack Basic Usage](#dandeliontrack-basic-usage)
    - [DandelionTrack Shape Types](#dandeliontrack-shape-types)
    - [Clustering Parameters](#clustering-parameters)
  - [Convenience Functions](#convenience-functions)
    - [lolliplot()](#lolliplot)
    - [dandelion\_plot()](#dandelion_plot)
  - [Composing with Other Tracks](#composing-with-other-tracks)
  - [API Reference](#api-reference)
    - [LolliplotTrack](#lolliplottrack-1)
    - [DandelionTrack](#dandeliontrack-1)
    - [Per-SNP DataFrame Columns](#per-snp-dataframe-columns)
  - [Complete Example](#complete-example)

---

## Introduction

### Core Concepts

```
LolliplotTrack   ─── lollipop-style variant/mutation plot (stems + shapes)
DandelionTrack   ─── clustered variant plot (stems fanning out from clusters)
```

Both tracks accept a **variant DataFrame** (`snp_data`) and an optional **feature DataFrame** (`features`) representing protein domains, exons, or other annotations along a horizontal baseline.

### Data Format

**snp_data** (variant DataFrame):

| Column | Required | Description |
|--------|----------|-------------|
| `chrom` | Yes | Chromosome name |
| `start` | Yes | Genomic position (1-based) |
| `score` | No | Numeric score (height). Default 1 |
| `label` | No | External text label |
| `fill` | No | Shape fill color |
| `border` | No | Shape edge color |
| `alpha` | No | Opacity (0-1) |

**features** (protein domain DataFrame):

| Column | Required | Description |
|--------|----------|-------------|
| `chrom` | Yes | Chromosome name |
| `start` | Yes | Feature start position |
| `end` | Yes | Feature end position |
| `name` | No | Feature label |
| `fill` | No | Rectangle fill color |
| `height` | No | Rectangle height (axes fraction) |

---

## LolliplotTrack

### Basic Usage

```python
import pandas as pd
import numpy as np
from geneview.genometracks import LolliplotTrack, GenomicInterval, plot_tracks

# Variant data
snp_data = pd.DataFrame({
    "chrom": ["chr1"] * 8,
    "start": [10, 100, 200, 400, 500, 700, 900, 1200],
    "score": [3, 5, 2, 7, 1, 4, 6, 3],
    "label": ["rs10", "rs100", "rs200", "rs400",
              "rs500", "rs700", "rs900", "rs1200"],
})

# Protein domain features
features = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 401, 801],
    "end": [150, 700, 1300],
    "name": ["Kinase", "SH2", "DNA-binding"],
    "fill": ["#FF8833", "#51C6E6", "#DFA32D"],
})

# Plot
track = LolliplotTrack(snp_data, features=features)
region = GenomicInterval("chr1", 0, 1400)
axes = plot_tracks([track], region=region, figsize=(12, 4))
axes[0].figure.savefig("lolliplot_basic.png", dpi=300, bbox_inches="tight")
```

**Constructor parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `snp_data` | — | Variant DataFrame |
| `features` | `None` | Protein domain DataFrame |
| `type` | `"circle"` | Shape type (see below) |
| `cex` | `1.0` | Global shape size multiplier |
| `dashline_col` | `"#CCCCCC"` | Guide line color |
| `show_yaxis` | `True` | Show y-axis scale |
| `label_on_feature` | `False` | Draw labels on feature rectangles |
| `lollipop_style_switch_limit` | `10` | Max score for Tanghulu stacking |
| `legend` | `None` | Legend specification dict |
| `jitter` | `None` | `"label"` for aligned labels |
| `yaxis` | `None` | Custom y-axis tick positions |
| `ylab` | `None` | Y-axis label text |
| `rescale` | `None` | Coordinate remapping rules |

### Shape Types

| Type | Description |
|------|-------------|
| `"circle"` | Filled circles (default). Supports Tanghulu stacking for integer scores |
| `"pie"` | Pie charts at each position (requires `pie_values` and `pie_colors` columns) |
| `"pin"` | Map-pin / teardrop shapes |
| `"flag"` | Rectangular banners with text inside |
| `"pie.stack"` | Stacked pie charts at grouped positions |

```python
# Pie type with per-site pie values
snp_pie = snp_data.copy()
snp_pie["score"] = 1  # reset score for pie type
snp_pie["pie_values"] = [[3, 7], [5, 5], [2, 8], [6, 4],
                          [1, 9], [4, 6], [7, 3], [5, 5]]
snp_pie["pie_colors"] = [["#87CEFA", "#98CE31"]] * 8

track = LolliplotTrack(snp_pie, features=features, type="pie")
```

### Per-SNP Customization

All of the following columns in `snp_data` are **optional** and override global defaults per row:

| Column | Purpose | Default |
| :--------: | :---------: | :---------: |
| `cex` | Per-SNP node size multiplier | Global `cex` |
| `node_label` | Text rendered **inside** the shape | None |
| `node_label_color` | Color of node label text | `"white"` |
| `node_label_size` | Font size of node label | `6 * cex` |
| `label_rotation` | External label rotation (degrees) | `90` |
| `label_color` | External label color | `"black"` |
| `dashline_col` | Per-SNP guide line color | Global `dashline_col` |
| `side` | `"top"` or `"bottom"` (caterpillar) | `"top"` |
| `shape` | Sub-shape name (see below) | `"circle"` |
| `stack_factor` | Grouping for `pie.stack` type | None |

```python
# Per-SNP customization example
snp_custom = snp_data.copy()
snp_custom["cex"] = [0.5, 1.0, 1.5, 2.0, 0.8, 1.2, 0.6, 1.8]
snp_custom["label_rotation"] = [0, 45, 90, 45, 0, 90, 45, 0]
snp_custom["label_color"] = ["red", "blue", "green", "orange",
                              "purple", "red", "blue", "green"]
snp_custom["dashline_col"] = [
    "#FF0000", "#0000FF", "#00FF00", "#FF8800", "#8800FF",
    "#FF0000", "#0000FF", "#00FF00"
]

track = LolliplotTrack(snp_custom, features=features)
```

### Tanghulu Stacking

When `type="circle"`, scores are integers, and `max(score) <= lollipop_style_switch_limit`, variants are drawn as **stacked circles** (Tanghulu style) — one circle per score unit.

```python
# Tanghulu: integer scores <= limit → stacked circles
snp_tanghulu = snp_data.copy()
snp_tanghulu["score"] = [2, 3, 1, 4, 1, 2, 5, 3]

track = LolliplotTrack(snp_tanghulu, features=features,
                       lollipop_style_switch_limit=10)
```

When `max(score) > limit`, the track switches to **proportional-height** mode (single circle whose y-position encodes the score).

### Node Labels

Render text **inside** each shape:

```python
snp_nl = snp_data.copy()
snp_nl["score"] = [2, 3, 1, 4, 1, 2, 5, 3]
snp_nl["node_label"] = ["A", "B", "C", "D", "E", "F", "G", "H"]
snp_nl["node_label_color"] = "white"
snp_nl["node_label_size"] = 5

track = LolliplotTrack(snp_nl, features=features)
```

### Caterpillar Layout

Place variants above and below the feature baseline using the `side` column:

```python
snp_cat = snp_data.copy()
snp_cat["score"] = [2, 3, 1, 4, 1, 2, 5, 3]
snp_cat["side"] = ["top", "bottom", "top", "bottom",
                    "top", "bottom", "top", "bottom"]

track = LolliplotTrack(snp_cat, features=features)
```

When `side` is detected, the feature baseline is centered vertically (at ~0.45) to leave room on both sides.

### Multiple Shapes in Tanghulu

Use list-valued `shape` and `fill` columns to alternate shapes/colors in stacked circles:

```python
snp_multi = pd.DataFrame({
    "chrom": ["chr1"] * 3,
    "start": [100, 400, 800],
    "score": [3, 4, 2],
    "fill": [["#FF0000", "#0000FF", "#00FF00"],
             ["#FF00FF", "#FFFF00"],
             ["#0080FF", "#E69F00"]],
    "shape": [["circle", "square", "diamond"],
              ["triangle_point_up", "circle"],
              ["square", "diamond"]],
})

track = LolliplotTrack(snp_multi, features=features)
```

### Legend

Pass a `legend` dict with `labels` and `fill` keys:

```python
legend = {
    "labels": ["Wild Type", "Mutant", "Unknown"],
    "fill": ["#87CEFA", "#98CE31", "#CCCCCC"],
}

track = LolliplotTrack(snp_data, features=features, legend=legend)
```

Optional legend keys: `ncol`, `loc`, `fontsize`.

### Aligned Labels (jitter="label")

When `jitter="label"`, all external labels are placed at a **uniform y-level at the top** of the plot area, connected by dashed lines from each shape. Labels that are too close horizontally are **automatically staggered** to different rows to avoid overlap — mirroring trackViewer's behavior.

```python
snp_jit = snp_data.copy()
snp_jit["dashline_col"] = ["#FF0000", "#0000FF", "#00FF00", "#FF8800",
                            "#8800FF", "#FF0000", "#0000FF", "#00FF00"]

track = LolliplotTrack(snp_jit, features=features, jitter="label")
```

**How anti-overlap works:**

1. All labels start at a base y-row above the tallest shape
2. Label widths are estimated (accounting for rotation angle)
3. When two labels overlap horizontally, the later one is shifted to the next row up
4. Dashed connector lines link each shape to its (possibly elevated) label

### Custom Y-axis and ylab

Control the y-axis ticks and label:

```python
track = LolliplotTrack(
    snp_data, features=features,
    yaxis=[0, 3, 7],       # custom tick positions
    ylab="# evidences",     # y-axis label
)
```

### pie.stack Type

Draw stacked pie charts at grouped positions. Requires `pie_values` and optionally `stack_factor`:

```python
snp_ps = pd.DataFrame({
    "chrom": ["chr1"] * 6,
    "start": [100, 100, 100, 400, 400, 400],
    "score": [1, 1, 1, 1, 1, 1],
    "pie_values": [[3, 7], [5, 5], [2, 8], [6, 4], [1, 9], [4, 6]],
    "pie_colors": [["#0080FF", "#E69F00"]] * 6,
    "stack_factor": ["A", "B", "C", "A", "B", "C"],
    "label": ["p1", "p2", "p3", "p4", "p5", "p6"],
})

track = LolliplotTrack(snp_ps, type="pie.stack")
```

### Coordinate Rescaling

Remap genomic positions using `rescale` — a list of `(from_start, from_end, to_start, to_end)` tuples:

```python
rescale_map = [
    (0, 500, 0, 250),       # compress first half
    (500, 1400, 300, 900),   # expand second half
]

track = LolliplotTrack(snp_data, features=features, rescale=rescale_map)
```

Positions are linearly interpolated within each mapping range; positions outside all ranges are unchanged.

### Multi-layer Features

Use `feature_layer_id` to draw features on separate baselines:

```python
features_ml = pd.DataFrame({
    "chrom": ["chr1"] * 5,
    "start": [1, 401, 801, 200, 600],
    "end": [150, 700, 1300, 350, 750],
    "name": ["Kinase", "SH2", "DNA-bind", "Motif X", "Motif Y"],
    "fill": ["#FF8833", "#51C6E6", "#DFA32D", "#CC79A7", "#009E73"],
    "feature_layer_id": [1, 1, 1, 2, 2],  # two separate layers
})

track = LolliplotTrack(snp_data, features=features_ml)
```

### Shape Types per SNP

Use the `shape` column to assign different shapes per variant:

| Shape name | Description |
| :------------: | :-------------: |
| `"circle"` | Filled circle (default) |
| `"square"` | Filled square |
| `"diamond"` | Rotated square |
| `"triangle_point_up"` | Upward-pointing triangle |
| `"triangle_point_down"` | Downward-pointing triangle |

```python
snp_shapes = snp_data.copy()
snp_shapes["shape"] = ["circle", "square", "diamond",
                        "triangle_point_up", "triangle_point_down",
                        "circle", "square", "diamond"]
snp_shapes["fill"] = ["#0080FF", "#E69F00", "#009E73", "#D55E00",
                       "#CC79A7", "#56B4E9", "#F0E442", "#FF8833"]

track = LolliplotTrack(snp_shapes, features=features)
```

---

## DandelionTrack

### DandelionTrack Basic Usage

Nearby variants are clustered into groups. Each cluster is drawn as a dandelion: a central stem with individual variants fanning out from the top.

```python
from geneview.genometracks import DandelionTrack

# Dense variant data
np.random.seed(42)
positions = np.sort(np.random.randint(10, 1400, 30))
dense_snps = pd.DataFrame({
    "chrom": ["chr1"] * 30,
    "start": positions,
    "score": np.random.randint(1, 8, 30),
    "label": [f"var{i}" for i in range(30)],
})

features = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 501, 1001],
    "end": [120, 900, 1405],
    "name": ["Domain A", "Kinase", "DNA Binding"],
    "fill": ["#FF8833", "#51C6E6", "#DFA32D"],
})

track = DandelionTrack(dense_snps, features=features)
region = GenomicInterval("chr1", 0, 1500)
axes = plot_tracks([track], region=region, figsize=(12, 4))
```

### DandelionTrack Shape Types

| Type | Description |
| :------: | :------------: |
| `"fan"` | Filled sector (default). Angular span encodes score |
| `"circle"` | Filled circles |
| `"pie"` | Pie charts (requires `pie_values`, `pie_colors`) |
| `"pin"` | Map-pin / teardrop shapes |

### Clustering Parameters

| Parameter | Default | Description |
| :---------: | :-------: | :-----------: |
| `maxgaps` | `1/50` | Max gap between clustered variants (fraction of region width) |
| `height_method` | `len` | Function `f(scores) -> float` for stem height |
| `cex` | `1.0` | Shape size multiplier |
| `show_yaxis` | `False` | Show y-axis |
| `rescale` | `None` | Coordinate remapping (same as LolliplotTrack) |

DandelionTrack also supports per-SNP `cex`, `label_rotation`, `label_color`, and `dashline_col` columns.

---

## Convenience Functions

### lolliplot()

```python
from geneview.genometracks import lolliplot

ax = lolliplot(
    snp_data,
    features=features,
    type="circle",         # or "pie", "pin", "flag", "pie.stack"
    region=None,           # auto-detect from data
    figsize=(12, 4),
    title="My Lolliplot",
    legend=None,           # legend dict
    jitter=None,           # "label" for aligned labels
    yaxis=None,            # custom y-axis ticks
    ylab=None,             # y-axis label
    rescale=None,          # coordinate remapping
    cex=1.0,               # global shape size
)
```

Returns a matplotlib Axes. When `ax=` is provided, draws into an existing axes instead of creating a new figure.

### dandelion_plot()

```python
from geneview.genometracks import dandelion_plot

ax = dandelion_plot(
    dense_snps,
    features=features,
    type="fan",            # or "circle", "pie", "pin"
    region=None,
    figsize=(12, 4),
    title="Dandelion Plot",
)
```

---

## Composing with Other Tracks

Both mutation tracks can be combined with `GenomeAxisTrack`, `AnnotationTrack`, and other track types in `plot_tracks()`:

```python
from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack,
    LolliplotTrack, plot_tracks, GenomicInterval,
)

region = GenomicInterval("chr1", 0, 1500)
tracks = [
    GenomeAxisTrack(),
    AnnotationTrack(features, name="Domains"),
    LolliplotTrack(snp_data, features=features, name="Mutations"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 7))
```

Multiple LolliplotTracks can be stacked for multi-sample comparison:

```python
tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(sample_a, features=features, name="Sample A"),
    LolliplotTrack(sample_b, features=features, name="Sample B"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 8))
```

---

## API Reference

### LolliplotTrack

```python
LolliplotTrack(
    snp_data,                          # pd.DataFrame (required)
    features=None,                     # pd.DataFrame
    type="circle",                     # str: circle|pie|pin|flag|pie.stack
    cex=1.0,                           # float: global shape size
    dashline_col="#CCCCCC",            # str: guide line color
    show_yaxis=True,                   # bool: show y-axis
    label_on_feature=False,            # bool: labels on rectangles
    lollipop_style_switch_limit=10,    # int: Tanghulu threshold
    legend=None,                       # dict: {labels, fill, ...}
    jitter=None,                       # str: "label" for aligned
    yaxis=None,                        # list: custom tick positions
    ylab=None,                         # str: y-axis label
    rescale=None,                      # list: [(from, to, from, to), ...]
    name="Lolliplot",                  # str: track name
    height=3.0,                        # float: relative height
)
```

### DandelionTrack

```python
DandelionTrack(
    snp_data,                          # pd.DataFrame (required)
    features=None,                     # pd.DataFrame
    type="fan",                        # str: fan|circle|pie|pin
    maxgaps=1/50,                      # float: clustering threshold
    height_method=None,                # callable: f(scores) -> float
    cex=1.0,                           # float: shape size
    show_yaxis=False,                  # bool: show y-axis
    label_on_feature=False,            # bool: labels on rectangles
    rescale=None,                      # list: coordinate remapping
    name="Dandelion",                  # str: track name
    height=3.0,                        # float: relative height
)
```

### Per-SNP DataFrame Columns

All optional; override constructor defaults per row:

`cex`, `node_label`, `node_label_color`, `node_label_size`,
`label_rotation`, `label_color`, `dashline_col`, `side`,
`shape`, `stack_factor`, `fill`, `border`, `alpha`,
`pie_values`, `pie_colors`, `label`

---

## Complete Example

```python
import numpy as np
import pandas as pd
from geneview.genometracks import (
    GenomeAxisTrack, LolliplotTrack, GenomicInterval, plot_tracks,
)

# --- Generate sample data ---
np.random.seed(42)
SNP = [10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402]
palette = ["#0080FF", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]

snp_data = pd.DataFrame({
    "chrom": ["chr1"] * len(SNP),
    "start": SNP,
    "label": [f"snp{s}" for s in SNP],
    "score": np.random.randint(1, 8, len(SNP)),
    "fill": [palette[i % len(palette)] for i in range(len(SNP))],
})

features = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 501, 1001],
    "end": [120, 900, 1405],
    "name": ["Domain A", "Kinase", "DNA Binding"],
    "fill": ["#FF8833", "#51C6E6", "#DFA32D"],
    "height": [0.04, 0.06, 0.08],
})

# --- Build tracks ---
region = GenomicInterval("chr1", 0, 1500)
tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(
        snp_data, features=features,
        jitter="label",
        legend={"labels": ["Type A", "Type B"], "fill": ["#0080FF", "#E69F00"]},
        yaxis=[0, 4, 8], ylab="Score",
        name="Mutations",
    ),
]

# --- Plot ---
axes = plot_tracks(tracks, region=region, figsize=(14, 6),
                   title="Comprehensive Lolliplot")
fig = axes[0].figure
fig.savefig("lolliplot_complete.png", dpi=300, bbox_inches="tight")
```

This produces a publication-quality figure with:

- Genome axis at the top
- Aligned labels with anti-overlap and dashed connector lines
- Colored shapes with per-SNP customization
- Y-axis with custom ticks and label
- Protein domain features on the baseline
- Legend in the upper right corner
