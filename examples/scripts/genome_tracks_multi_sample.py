"""Multi-sample data -- heatmap, overlay, and grouped comparisons.

Uses the ``multi_sample.tsv`` file (6 samples: 3 control + 3 treatment)
to demonstrate DataTrack with heatmap, overlay, and per-sample line plots.

Run:  python examples/scripts/genome_tracks_multi_sample.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack,
    DataTrack,
    OverlayTrack,
    GenomicInterval,
    plot_tracks,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

TSV_FILE = os.path.join(DATA_DIR, "multi_sample.tsv")
region = GenomicInterval("chr7", 26_500_000, 26_800_000)
gtrack = GenomeAxisTrack()

# ---------------------------------------------------------------------------
# Load multi-sample data
# ---------------------------------------------------------------------------
multi_data = pd.read_csv(TSV_FILE, sep="\t")
# Normalize column names to lowercase (matching geneview conventions)
multi_data.columns = [c.lower().strip() for c in multi_data.columns]
print(f"[INFO] Loaded multi_sample.tsv: {multi_data.shape[0]} bins x "
      f"{multi_data.shape[1]} columns")
print(f"       Columns: {list(multi_data.columns)}")

# Column groups
ctrl_cols = ["sample_ctrl_1", "sample_ctrl_2", "sample_ctrl_3"]
treat_cols = ["sample_treat_1", "sample_treat_2", "sample_treat_3"]
all_sample_cols = ctrl_cols + treat_cols

figs = []
labels = []


# ---------------------------------------------------------------------------
# 1. Heatmap with all 6 samples
# ---------------------------------------------------------------------------
dtrack_heat = DataTrack(
    multi_data, type="heatmap", name="All Samples (heatmap)",
    show_sample_names=True, separator=3,
)
axes = plot_tracks(
    [gtrack, dtrack_heat], region=region,
    title="Multi-sample heatmap (6 samples, ctrl | treat)",
    figsize=(14, 5),
)
figs.append(axes[0].figure)
labels.append("multi_heatmap")


# ---------------------------------------------------------------------------
# 2. Heatmap — control samples only
# ---------------------------------------------------------------------------
ctrl_data = multi_data[["chrom", "start", "end"] + ctrl_cols].copy()
dtrack_ctrl = DataTrack(
    ctrl_data, type="heatmap", name="Control (heatmap)",
    show_sample_names=True,
)
axes = plot_tracks(
    [gtrack, dtrack_ctrl], region=region,
    title="Control samples heatmap", figsize=(14, 4),
)
figs.append(axes[0].figure)
labels.append("multi_heatmap_ctrl")


# ---------------------------------------------------------------------------
# 3. Overlay — control vs treatment mean
# ---------------------------------------------------------------------------
# Compute mean signal per group
overlay_data = multi_data[["chrom", "start", "end"]].copy()
overlay_data["ctrl_mean"] = multi_data[ctrl_cols].mean(axis=1)
overlay_data["treat_mean"] = multi_data[treat_cols].mean(axis=1)

dtrack_overlay = DataTrack(
    overlay_data, type="line", name="Ctrl vs Treat (mean)",
)
axes = plot_tracks(
    [gtrack, dtrack_overlay], region=region,
    title="Overlay: control mean vs treatment mean", figsize=(14, 4),
)
figs.append(axes[0].figure)
labels.append("multi_overlay_mean")


# ---------------------------------------------------------------------------
# 4. Per-sample line overlay (all 6 samples)
# ---------------------------------------------------------------------------
dtrack_all = DataTrack(
    multi_data, type="line", name="All samples (line)",
)
axes = plot_tracks(
    [gtrack, dtrack_all], region=region,
    title="All 6 samples — line overlay", figsize=(14, 5),
)
figs.append(axes[0].figure)
labels.append("multi_all_lines")


# ---------------------------------------------------------------------------
# 5. Heatmap + mean line — stacked tracks
# ---------------------------------------------------------------------------
dtrack_heat2 = DataTrack(
    multi_data, type="heatmap", name="Heatmap",
    show_sample_names=True, separator=3,
)
dtrack_mean = DataTrack(
    overlay_data, type="line", name="Mean signal",
)
axes = plot_tracks(
    [gtrack, dtrack_heat2, dtrack_mean],
    region=region,
    sizes=[0.2, 1.5, 1.0],
    title="Heatmap + mean signal (stacked)",
    figsize=(14, 7),
)
figs.append(axes[0].figure)
labels.append("multi_heatmap_mean")


# ---------------------------------------------------------------------------
# 6. Polygon overlay — control vs treatment
# ---------------------------------------------------------------------------
dtrack_poly = DataTrack(
    overlay_data, type="polygon", name="Ctrl vs Treat (polygon)",
)
axes = plot_tracks(
    [gtrack, dtrack_poly], region=region,
    title="Polygon overlay: control vs treatment", figsize=(14, 4),
)
figs.append(axes[0].figure)
labels.append("multi_polygon_overlay")


# ---------------------------------------------------------------------------
# Save all figures
# ---------------------------------------------------------------------------
for label, fig in zip(labels, figs):
    path = os.path.join(OUT_DIR, f"genome_tracks_{label}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {path}")

plt.show()
print(f"\n[DONE] Generated {len(figs)} multi-sample figures.")
