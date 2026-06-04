#!/usr/bin/env python
"""
Custom drawing callbacks — pre/post rendering on genome tracks.

Mirrors genomeview's pre-renderer/post-renderer pattern (Cell 14 of the
examples notebook).  In genometracks, since ``plot_tracks()`` returns
standard matplotlib Axes, users can draw custom annotations before or
after rendering tracks — no special callback mechanism is needed.

Examples:
  1. Draw a shaded region overlay before track rendering (pre-render)
  2. Add text annotations after track rendering (post-render)
  3. Combined: pre-render highlight + post-render text + custom markers

Run:  python examples/scripts/genome_tracks_custom_drawing.py
"""

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "genome_tracks")
FIG_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Data paths
BAM_CHR1 = os.path.join(DATA_DIR, "illumina.chr1.bam")
ILLUMINA_BAM = os.path.join(DATA_DIR, "illumina.chr14.bam")
BED_FILE = os.path.join(DATA_DIR, "chr1_200mb.refseq.sorted.bed.gz")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    AnnotationTrack,
    GenomicInterval,
    read_bed,
    save_figure,
)

# =========================================================================
# 1. Pre-render: shaded region overlay before drawing tracks
# =========================================================================
print("[1] Custom drawing — pre-render shaded region")

chrom = "14"
start = 66901400
region = GenomicInterval(chrom, start, start + 5000)

# Create tracks
tracks = [
    GenomeAxisTrack(),
    BAMCoverageTrack(filepath=ILLUMINA_BAM, name="Coverage"),
]

# Render the tracks first
axes = plot_tracks(tracks, region=region, figsize=(14, 6))

# Now draw custom overlays (post-render on the returned axes)
# Highlight a specific region with a green shading
highlight_start = start + 1500
highlight_end = start + 3500

for ax in axes:
    ax.axvspan(highlight_start, highlight_end, alpha=0.15, color="green",
               zorder=0)

# Add a text label above the highlighted region
axes[0].text(
    (highlight_start + highlight_end) / 2, 1.05,
    "Region of Interest",
    ha="center", va="bottom", fontsize=9, color="green",
    fontweight="bold", transform=axes[0].transData,
)

save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_custom_prerender.png"))
plt.close("all")
print("  → genome_tracks_custom_prerender.png")

# =========================================================================
# 2. Post-render: text annotations and markers after drawing
# =========================================================================
print("[2] Custom drawing — post-render annotations")

region_zoom = GenomicInterval(chrom, start, start + 1500)

tracks_zoom = [
    GenomeAxisTrack(),
    AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Illumina Reads",
        color_fn=lambda r: "lightgray",
    ),
]

axes = plot_tracks(tracks_zoom, region=region_zoom, figsize=(14, 7))

# Add post-render annotations on the pileup axes (last axis)
pileup_ax = axes[-1]

# Draw an arrow pointing to a specific position
arrow_pos = start + 750
pileup_ax.annotate(
    "SNP",
    xy=(arrow_pos, 0.5), xytext=(arrow_pos + 200, 0.8),
    fontsize=8, color="red", fontweight="bold",
    arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
    xycoords="data",
)

# Add a vertical dashed line at a position of interest
pileup_ax.axvline(start + 750, color="red", linestyle="--", linewidth=1,
                  alpha=0.7, zorder=10)

# Add a text box with summary statistics
pileup_ax.text(
    0.02, 0.02,
    f"Region: {chrom}:{start}-{start+1500}\nCustom annotation added",
    transform=pileup_ax.transAxes,
    fontsize=7, color="#333333",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
)

save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_custom_postrender.png"))
plt.close("all")
print("  → genome_tracks_custom_postrender.png")

# =========================================================================
# 3. Combined: pre-render + post-render with multi-track
# =========================================================================
print("[3] Combined custom drawing — pre + post render")

tracks_combo = [
    GenomeAxisTrack(),
    BAMCoverageTrack(filepath=ILLUMINA_BAM, name="Coverage"),
    AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Reads",
        color_fn=lambda r: "lightgray",
    ),
]

axes = plot_tracks(tracks_combo, region=region, figsize=(14, 10))

# Pre-render style: add background shading to all axes
region_a_start = start + 500
region_a_end = start + 1500
region_b_start = start + 3000
region_b_end = start + 4000

for ax in axes:
    # Region A shading
    ax.axvspan(region_a_start, region_a_end, alpha=0.12, color="blue", zorder=0)
    # Region B shading
    ax.axvspan(region_b_start, region_b_end, alpha=0.12, color="orange", zorder=0)

# Post-render: add labels on the coverage axis
cov_ax = axes[1]
cov_ax.text(
    (region_a_start + region_a_end) / 2, 0.95,
    "Region A", ha="center", va="top", fontsize=8, color="blue",
    fontweight="bold",
)
cov_ax.text(
    (region_b_start + region_b_end) / 2, 0.95,
    "Region B", ha="center", va="top", fontsize=8, color="darkorange",
    fontweight="bold",
)

# Post-render: add custom markers on the pileup axis
pileup_ax = axes[-1]
for pos in [start + 1000, start + 2000, start + 3500]:
    pileup_ax.plot(pos, 0.5, marker="v", markersize=8, color="red",
                   zorder=10, markeredgecolor="darkred")

# Add overall title
axes[0].set_title(
    "Custom Drawing: Pre-render Highlights + Post-render Annotations",
    fontsize=11, loc="left",
)

save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_custom_combined.png"))
plt.close("all")
print("  → genome_tracks_custom_combined.png")

# =========================================================================
# 4. Modular layout with custom annotations (genomeview Cell 12 style)
# =========================================================================
print("[4] Modular side-by-side views with custom annotations")

from geneview.genometracks import plot_tracks_grid

region1 = GenomicInterval(chrom, 66901400, 66901400 + 1000)
region2 = GenomicInterval(chrom, 66901400 + 5000, 66901400 + 6000)

view1_tracks = [
    GenomeAxisTrack(),
    AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Region 1",
        color_fn=lambda r: "lightgray",
    ),
]
view2_tracks = [
    GenomeAxisTrack(),
    AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Region 2",
        color_fn=lambda r: "lightgray",
    ),
]

all_axes = plot_tracks_grid(
    [view1_tracks, view2_tracks],
    regions=[region1, region2],
    columns=2,
    figsize=(16, 6),
    title="Side-by-Side Modular Views (genomeview Cell 12 style)",
)

# Add custom post-render annotations on each view
for view_axes, label in zip(all_axes, ["Region 1", "Region 2"]):
    ax = view_axes[-1]
    ax.text(
        0.98, 0.98, label,
        transform=ax.transAxes, fontsize=9,
        ha="right", va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

save_figure(
    [ax for view_axes in all_axes for ax in view_axes],
    os.path.join(FIG_DIR, "genome_tracks_custom_modular.png"),
)
plt.close("all")
print("  → genome_tracks_custom_modular.png")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All custom drawing figures generated in:", FIG_DIR)
