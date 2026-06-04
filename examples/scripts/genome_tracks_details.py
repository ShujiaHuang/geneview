"""
DetailsAnnotationTrack Example
===============================

Demonstrates the DetailsAnnotationTrack which extends AnnotationTrack
with detail panels below each feature.
"""
import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, DetailsAnnotationTrack, GenomicInterval,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")

# Create annotation data
data = pd.DataFrame({
    "chrom": ["chr7"] * 4,
    "start": [1000, 2000, 3000, 4000],
    "end": [1500, 2800, 3600, 4400],
    "strand": ["+", "-", "+", "-"],
    "name": ["geneA", "geneB", "geneC", "geneD"],
})

# ---------------------------------------------------------------------------
# 1. Default detail panels
# ---------------------------------------------------------------------------
track = DetailsAnnotationTrack(data, name="Details (default)")
region = GenomicInterval("chr7", 800, 4600)

axes = plot_tracks(
    [GenomeAxisTrack(), track],
    region=region,
    figsize=(12, 4),
    title="DetailsAnnotationTrack - Default Details",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_details_default.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Custom detail function
# ---------------------------------------------------------------------------
def my_detail_fun(ax, identifier, region, **kwargs):
    """Draw a simple text label in the detail panel."""
    ax.text(0.5, 0.5, f"Detail: {identifier}",
            ha="center", va="center", fontsize=6,
            transform=ax.transAxes)

track_custom = DetailsAnnotationTrack(
    data,
    fun=my_detail_fun,
    details_connector_col="blue",
    details_connector_lty="-",
    name="Details (custom)",
)

axes = plot_tracks(
    [GenomeAxisTrack(), track_custom],
    region=region,
    figsize=(12, 4),
    title="DetailsAnnotationTrack - Custom Function",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_details_custom.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Selective details (only certain features)
# ---------------------------------------------------------------------------
track_selective = DetailsAnnotationTrack(
    data,
    select_fun=lambda row: row["name"] in ("geneA", "geneC"),
    details_min_width=400,
    name="Selective Details",
)

axes = plot_tracks(
    [GenomeAxisTrack(), track_selective],
    region=region,
    figsize=(12, 4),
    title="DetailsAnnotationTrack - Selective",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_details_selective.png"), dpi=150, bbox_inches="tight")
plt.close("all")

print("All DetailsAnnotationTrack examples saved.")
