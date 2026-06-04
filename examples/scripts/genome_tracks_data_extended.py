"""
Extended DataTrack Plot Types Example
=====================================

Demonstrates the new DataTrack plot types: average (a), confint,
smooth, horizon, grid (g), and regression (r), plus composite types.
"""
import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, DataTrack, GenomicInterval,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")

# Create multi-sample data
rng = np.random.RandomState(42)
n = 50
starts = np.linspace(1000, 2000, n, dtype=int)
data = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts,
    "end": starts + 20,
    "sample_A": rng.randn(n).cumsum(),
    "sample_B": rng.randn(n).cumsum() + 2,
    "sample_C": rng.randn(n).cumsum() - 1,
})

region = GenomicInterval("chr7", 1000, 2050)

# ---------------------------------------------------------------------------
# 1. Average (type "a")
# ---------------------------------------------------------------------------
dt_avg = DataTrack(data, type="a", name="Average", col="#3C5488")
axes = plot_tracks([GenomeAxisTrack(), dt_avg], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_average.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Confidence Interval (type "confint")
# ---------------------------------------------------------------------------
dt_ci = DataTrack(data, type="confint", name="Conf. Interval", col="#00A087")
axes = plot_tracks([GenomeAxisTrack(), dt_ci], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_confint.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Smooth/Loess (type "smooth")
# ---------------------------------------------------------------------------
dt_smooth = DataTrack(data, type="smooth", name="Smooth", smooth_span=0.3)
axes = plot_tracks([GenomeAxisTrack(), dt_smooth], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_smooth.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Horizon plot (type "horizon")
# ---------------------------------------------------------------------------
# Center data around zero for horizon plot
data_hz = data.copy()
mean_vals = data_hz[["sample_A", "sample_B", "sample_C"]].mean(axis=1)
for col in ["sample_A", "sample_B", "sample_C"]:
    data_hz[col] = data_hz[col] - mean_vals.mean()

dt_horizon = DataTrack(data_hz, type="horizon", name="Horizon", baseline=0)
axes = plot_tracks([GenomeAxisTrack(), dt_horizon], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_horizon.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. Grid type (type "g")
# ---------------------------------------------------------------------------
dt_grid = DataTrack(data, type="g", name="Grid")
axes = plot_tracks([GenomeAxisTrack(), dt_grid], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_grid_type.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 6. Regression line (type "r")
# ---------------------------------------------------------------------------
dt_reg = DataTrack(data, type="r", name="Regression")
axes = plot_tracks([GenomeAxisTrack(), dt_reg], region=region, figsize=(10, 3))
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_regression.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 7. Composite: boxplot + average + grid
# ---------------------------------------------------------------------------
dt_composite = DataTrack(
    data, type=["boxplot", "a", "g"],
    name="Composite",
)
axes = plot_tracks(
    [GenomeAxisTrack(), dt_composite],
    region=region,
    figsize=(10, 4),
    title="Composite: boxplot + average + grid",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_data_composite.png"), dpi=150, bbox_inches="tight")
plt.close("all")

print("All extended DataTrack examples saved.")
