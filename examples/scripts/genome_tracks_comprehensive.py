"""Comprehensive genome tracks showcase -- all track types in one plot.

Combines IdeogramTrack, GenomeAxisTrack, AnnotationTrack, GeneRegionTrack,
DataTrack, HighlightTrack, and OverlayTrack into a single multi-panel figure
with direction indicators, per-region highlight colors, and fractional
region extension.

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
    HighlightTrack, OverlayTrack, IdeogramTrack, GenomicInterval, plot_tracks,
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
gtrack = GenomeAxisTrack(little_ticks=True, add53=True)

# IdeogramTrack with synthetic cytoband data for chr7
bands_chr7 = pd.DataFrame({
    "chrom": ["chr7"] * 10,
    "chromStart": [0, 7100000, 19000000, 27400000, 40800000,
                   58800000, 62900000, 83200000, 103400000, 127500000],
    "chromEnd": [7100000, 19000000, 27400000, 40800000, 58800000,
                 62900000, 83200000, 103400000, 127500000, 159138663],
    "name": ["p22.3", "p21.3", "p21.1", "p15.3", "p14.3",
             "acen", "p13.1", "p12.1", "q11.23", "q21.11"],
    "gieStain": ["gneg", "gpos75", "gpos25", "gneg", "gpos75",
                 "acen", "gneg", "gneg", "gpos50", "gpos100"],
})
ideo = IdeogramTrack(bands_chr7, chromosome="chr7")

atrack_cpg = AnnotationTrack(cpg_data, name="CpG Islands")
atrack_ann = AnnotationTrack(ann_data, name="Regulatory",
                             shape="ellipse",
                             display_params={"fill": "#DC0000"})
grtrack = GeneRegionTrack(gene_data, name="Gene Models",
                          collapse_transcripts="longest")
dtrack = DataTrack(cov_data, type="histogram", name="Coverage",
                   display_params={"grid": True})

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
    fill=["#FFF3BF", "#FFD6E0"],
    alpha=0.3,
)

# ---------------------------------------------------------------------------
# Plot everything
# ---------------------------------------------------------------------------
axes = plot_tracks(
    [ideo, gtrack, ht],
    region=region,
    sizes=[0.4, 0.2, 0.8, 0.8, 1.5, 1.5],  # relative heights
    title="Comprehensive Genome Tracks",
    figsize=(16, 12),
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
