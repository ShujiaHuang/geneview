"""DataTrack -- visualize numeric genomic data with multiple plot types.

Demonstrates DataTrack with line, histogram, polygon, heatmap, points,
mountain, and gradient plot types. Also shows multi-sample data and
custom y-axis limits.

Run:  python examples/scripts/genome_tracks_data.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, DataTrack, GenomicInterval, plot_tracks, read_bedgraph,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
BEDGRAPH_FILE = os.path.join(DATA_DIR, "coverage.bedgraph")

# Load coverage data
cov_data = read_bedgraph(BEDGRAPH_FILE)
region = GenomicInterval("chr7", 26_500_000, 26_800_000)
gtrack = GenomeAxisTrack()

# ---------------------------------------------------------------------------
# 1. Single DataTrack with different plot types
# ---------------------------------------------------------------------------
plot_types = ["line", "histogram", "polygon", "points", "mountain", "gradient"]
figs = []

for ptype in plot_types:
    dtrack = DataTrack(cov_data, type=ptype, name=f"{ptype.title()}")
    axes = plot_tracks(
        [gtrack, dtrack], region=region,
        title=f"DataTrack: {ptype}", figsize=(12, 4),
    )
    figs.append(axes[0].figure)

# ---------------------------------------------------------------------------
# 2. Heatmap (needs multi-column data)
# ---------------------------------------------------------------------------
# Create multi-column data programmatically
rng = np.random.RandomState(99)
n = 50
starts = np.linspace(26_500_000, 26_790_000, n, dtype=int)
ends = starts + 6000
heat_data = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts,
    "end": ends,
    "sample_A": rng.randn(n).cumsum() / 3,
    "sample_B": rng.randn(n).cumsum() / 3 + 2,
    "sample_C": rng.randn(n).cumsum() / 3 - 1,
    "sample_D": rng.randn(n).cumsum() / 3 + 1,
})
dtrack_heat = DataTrack(heat_data, type="heatmap", name="Heatmap")
axes_heat = plot_tracks(
    [gtrack, dtrack_heat], region=region,
    title="DataTrack: heatmap (4 samples)", figsize=(12, 5),
)
figs.append(axes_heat[0].figure)

# ---------------------------------------------------------------------------
# 3. Custom y-axis limits
# ---------------------------------------------------------------------------
dtrack_ylim = DataTrack(cov_data, type="line", name="Custom YLim",
                        ylim=(-5, 5))
axes_ylim = plot_tracks(
    [gtrack, dtrack_ylim], region=region,
    title="DataTrack: custom ylim=(-5, 5)", figsize=(12, 4),
)
figs.append(axes_ylim[0].figure)

# ---------------------------------------------------------------------------
# Save all figures
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(out_dir, exist_ok=True)
labels = plot_types + ["heatmap", "custom_ylim"]
for label, fig in zip(labels, figs):
    path = os.path.join(out_dir, f"genome_tracks_data_{label}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {path}")

plt.show()
