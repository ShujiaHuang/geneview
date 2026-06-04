"""BAM and CRAM coverage -- visualize alignment depth as DataTracks.

This script demonstrates how to compute coverage from BAM and CRAM alignment
files and display them as genome track DataTracks. Because BAM/CRAM binary
files are not shipped with the repository, the script uses **synthetic
coverage data** that mimics the output of ``read_bam_coverage`` /
``read_cram_coverage`` so the figure can be generated without any external
data.

The real API calls (commented out) are shown alongside the synthetic
workflow so you can easily adapt this to your own alignments.

Real usage (requires indexed BAM/CRAM files)::

    from geneview.genometracks import read_bam_coverage, read_cram_coverage

    bam_cov = read_bam_coverage("sample.bam", region=region)
    cram_cov = read_cram_coverage("sample.cram", region=region, reference="hg38.fa")

Run:  python examples/scripts/genome_tracks_bam_cram.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, DataTrack, HighlightTrack,
    GenomicInterval, plot_tracks, read_bed,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Define the region of interest
# ---------------------------------------------------------------------------
region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# Generate synthetic coverage data (mimics read_bam_coverage / read_cram_coverage output)
# ---------------------------------------------------------------------------
# In real usage, you would replace these blocks with:
#
#   bam_cov = read_bam_coverage("sample.bam", region=region, bins=500)
#   cram_cov = read_cram_coverage("sample.cram", region=region,
#                                  reference="hg38.fa", bins=500)
#
rng = np.random.RandomState(42)
n_bins = 230  # matches typical bedGraph resolution in this region
bin_size = (region.end - region.start) // n_bins

starts = np.arange(region.start, region.end - bin_size + 1, bin_size)
ends = starts + bin_size

# Simulate two samples with different coverage profiles
bam_values = np.abs(rng.poisson(lam=30, size=n_bins) + rng.randn(n_bins).cumsum() * 0.3)
cram_values = np.abs(rng.poisson(lam=28, size=n_bins) + rng.randn(n_bins).cumsum() * 0.25 + 2)

bam_cov = pd.DataFrame({
    "chrom": [region.chrom] * n_bins,
    "start": starts[:n_bins],
    "end": ends[:n_bins],
    "value": bam_values[:n_bins],
})
cram_cov = pd.DataFrame({
    "chrom": [region.chrom] * n_bins,
    "start": starts[:n_bins],
    "end": ends[:n_bins],
    "value": cram_values[:n_bins],
})

# ---------------------------------------------------------------------------
# 1. Side-by-side BAM vs CRAM coverage (histogram)
# ---------------------------------------------------------------------------
gtrack = GenomeAxisTrack(little_ticks=True)

bam_track = DataTrack(bam_cov, type="histogram", name="BAM Coverage",
                      display_params={"col": "#2196F3", "fill": "#90CAF9"})
cram_track = DataTrack(cram_cov, type="histogram", name="CRAM Coverage",
                       display_params={"col": "#4CAF50", "fill": "#A5D6A7"})

axes1 = plot_tracks(
    [gtrack, bam_track, cram_track],
    region=region,
    title="BAM vs CRAM Coverage (histogram)",
    figsize=(14, 6),
)
fig1 = axes1[0].figure
fig1.savefig(os.path.join(OUT_DIR, "genome_tracks_bam_cram_histogram.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_bam_cram_histogram.png")

# ---------------------------------------------------------------------------
# 2. Line plot overlay of BAM and CRAM signals
# ---------------------------------------------------------------------------
bam_line = DataTrack(bam_cov, type="line", name="BAM (line)",
                     display_params={"col": "#2196F3", "lwd": 1.5})
cram_line = DataTrack(cram_cov, type="line", name="CRAM (line)",
                      display_params={"col": "#4CAF50", "lwd": 1.5, "alpha": 0.8})

axes2 = plot_tracks(
    [gtrack, bam_line, cram_line],
    region=region,
    title="BAM vs CRAM Coverage (line overlay)",
    figsize=(14, 5),
)
fig2 = axes2[0].figure
fig2.savefig(os.path.join(OUT_DIR, "genome_tracks_bam_cram_line.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_bam_cram_line.png")

# ---------------------------------------------------------------------------
# 3. Combined view: annotation + BAM coverage + highlights
# ---------------------------------------------------------------------------
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
atrack = AnnotationTrack(cpg_data, name="CpG Islands")

bam_hist = DataTrack(bam_cov, type="histogram", name="BAM Coverage",
                     display_params={"col": "#2196F3", "fill": "#BBDEFB",
                                     "grid": True})

ht = HighlightTrack(
    regions=pd.DataFrame({
        "chrom": ["chr7", "chr7"],
        "start": [26_505_000, 26_600_000],
        "end":   [26_535_000, 26_665_000],
    }),
    track_list=[atrack, bam_hist],
    fill="#FFF3BF", alpha=0.3,
)

axes3 = plot_tracks(
    [gtrack, ht],
    region=region,
    title="Annotation + BAM Coverage with Highlighted Regions",
    figsize=(14, 7),
    extend_left=0.03, extend_right=0.03,
)
fig3 = axes3[0].figure
fig3.savefig(os.path.join(OUT_DIR, "genome_tracks_bam_cram_combined.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_bam_cram_combined.png")

plt.show()
print("[INFO] All BAM/CRAM example figures saved.")
