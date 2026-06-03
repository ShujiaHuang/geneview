"""HighlightTrack -- cross-track highlighting of genomic regions.

Demonstrates wrapping multiple tracks in a HighlightTrack to draw
semi-transparent colored regions across all contained tracks.

Run:  python examples/scripts/genome_tracks_highlight.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
    HighlightTrack, GenomicInterval, plot_tracks, read_bed, read_gff,
    read_bedgraph,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

# Load data
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
cov_data = read_bedgraph(os.path.join(DATA_DIR, "coverage.bedgraph"))

region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# Create individual tracks
gtrack = GenomeAxisTrack()
atrack = AnnotationTrack(cpg_data, name="CpG Islands")
grtrack = GeneRegionTrack(gene_data, name="Gene Models")
dtrack = DataTrack(cov_data, type="line", name="Coverage")

# ---------------------------------------------------------------------------
# 1. Single highlight region (yellow)
# ---------------------------------------------------------------------------
highlight_regions = pd.DataFrame({
    "chrom": ["chr7", "chr7"],
    "start": [26_520_000, 26_640_000],
    "end":   [26_540_000, 26_660_000],
})
ht = HighlightTrack(
    regions=highlight_regions,
    track_list=[atrack, grtrack, dtrack],
    fill="#FFFF99",
    alpha=0.35,
)

axes1 = plot_tracks(
    [gtrack, ht],
    region=region,
    title="Highlight: Two Regions (yellow)",
    figsize=(14, 7),
)
fig1 = axes1[0].figure

# ---------------------------------------------------------------------------
# 2. Multiple highlight tracks with different colors
# ---------------------------------------------------------------------------
ht_green = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7"], "start": [26_500_000], "end": [26_530_000],
    }),
    track_list=[atrack, grtrack],
    fill="#90EE90", alpha=0.3,
)
ht_red = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7"], "start": [26_600_000], "end": [26_665_000],
    }),
    track_list=[dtrack],
    fill="#FFB6C1", alpha=0.3,
)

axes2 = plot_tracks(
    [gtrack, ht_green, grtrack, ht_red],
    region=region,
    title="Multiple Highlight Tracks (green + pink)",
    figsize=(14, 7),
)
fig2 = axes2[0].figure

# ---------------------------------------------------------------------------
# 3. Per-region colors (list of fill/col)
# ---------------------------------------------------------------------------
ht_perregion = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7", "chr7", "chr7"],
        "start": [26_500_000, 26_560_000, 26_640_000],
        "end":   [26_530_000, 26_590_000, 26_670_000],
    }),
    track_list=[atrack, grtrack, dtrack],
    fill=["#FF9999", "#99FF99", "#9999FF"],
    col=["#CC0000", "#00CC00", "#0000CC"],
    alpha=0.3,
)

axes3 = plot_tracks(
    [gtrack, ht_perregion],
    region=region,
    title="Per-Region Colors (red / green / blue)",
    figsize=(14, 7),
)
fig3 = axes3[0].figure

# ---------------------------------------------------------------------------
# Save figures
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(out_dir, exist_ok=True)
for i, fig in enumerate([fig1, fig2, fig3], 1):
    path = os.path.join(out_dir, f"genome_tracks_highlight_{i}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {path}")

plt.show()
