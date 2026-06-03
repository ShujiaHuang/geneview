"""Comprehensive genome tracks showcase -- all track types in one plot.

Combines GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
and HighlightTrack into a single multi-panel figure with custom sizes,
title, and region extension.

Run:  python examples/scripts/genome_tracks_comprehensive.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
    HighlightTrack, GenomicInterval, plot_tracks,
    read_bed, read_gff, read_bedgraph,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
cov_data = read_bedgraph(os.path.join(DATA_DIR, "coverage.bedgraph"))
ann_data = read_bed(os.path.join(DATA_DIR, "annotations.bed"))

# ---------------------------------------------------------------------------
# Define region
# ---------------------------------------------------------------------------
region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# Create tracks
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Highlight two regions of interest
# ---------------------------------------------------------------------------
ht = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7", "chr7"],
        "start": [26_505_000, 26_600_000],
        "end":   [26_535_000, 26_665_000],
    }),
    track_list=[atrack_cpg, atrack_ann, grtrack, dtrack],
    fill="#FFF3BF",
    alpha=0.3,
)

# ---------------------------------------------------------------------------
# Plot everything
# ---------------------------------------------------------------------------
axes = plot_tracks(
    [gtrack, ht],
    region=region,
    sizes=[0.2, 0.8, 0.8, 1.5, 1.5],  # relative heights (gtrack is auto)
    title="Comprehensive Genome Tracks",
    figsize=(16, 10),
    extend_left=0.05,
    extend_right=0.05,
)
fig = axes[0].figure
plt.suptitle("geneview Genome Tracks — Full Showcase",
             fontsize=14, fontweight="bold", y=1.01)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "genome_tracks_comprehensive.png")
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"[INFO] Saved comprehensive figure to {out_path}")
plt.show()
