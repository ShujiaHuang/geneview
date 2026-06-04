#!/usr/bin/env python
"""
Custom read coloring via ``color_fn`` — mirrors genomeview's color_fn demo.

Demonstrates the AlignmentsTrack ``color_fn`` parameter, equivalent to
genomeview's ``track.color_fn = lambda x: ...`` pattern.

Examples:
  1. Color by insert size (paired-end data)
  2. Uniform gray coloring (useful as a backdrop for variant tracks)
  3. Color by mapping quality
  4. Side-by-side comparison with default coloring

Run:  python examples/scripts/genome_tracks_custom_colors.py
"""

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "genome_tracks")
FIG_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Data paths
ILLUMINA_BAM = os.path.join(DATA_DIR, "illumina.chr14.bam")
PACBIO_BAM = os.path.join(DATA_DIR, "pacbio.chr14.bam")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    GenomicInterval,
    save_figure,
)

# Check data availability
if not os.path.exists(ILLUMINA_BAM):
    print(f"[SKIP] Missing {ILLUMINA_BAM}")
    sys.exit(0)

# =========================================================================
# 1. Color by insert size (genomeview Cell 8 equivalent)
# =========================================================================
print("[1] color_fn — color by insert size (paired-end)")

chrom = "14"
start = 66901400
region = GenomicInterval(chrom, start, start + 10000)


def color_by_insert_size(read):
    """Color reads by insert size, matching genomeview's example.

    - Red: insert < 100bp or > 1500bp (abnormal)
    - Blue: insert > 550bp (large)
    - Green: normal insert
    """
    try:
        isize = abs(read.template_length)  # pysam equivalent of isize
    except (AttributeError, TypeError):
        return "#BDBDBD"
    if isize < 100 or isize > 1500:
        return "red"
    if isize > 550:
        return "blue"
    return "green"


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

# Default coloring
track_default = AlignmentsTrack(
    filepath=ILLUMINA_BAM, is_paired=True,
    type="pileup", name="Default",
)
track_default.draw(ax1, region)
ax1.set_title("Paired-end Pileup — Default Coloring", fontsize=10, loc="left")

# Insert-size coloring
track_insert = AlignmentsTrack(
    filepath=ILLUMINA_BAM, is_paired=True,
    type="pileup", name="Insert Size",
    color_fn=color_by_insert_size,
)
track_insert.draw(ax2, region)
ax2.set_title("Paired-end Pileup — color_fn by Insert Size", fontsize=10, loc="left")

# Add a legend for the custom colors
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="green", label="Normal (100-550bp)"),
    Patch(facecolor="blue", label="Large (550-1500bp)"),
    Patch(facecolor="red", label="Abnormal (<100 or >1500bp)"),
]
ax2.legend(handles=legend_elements, loc="upper right", fontsize=7,
           framealpha=0.8, title="Insert Size")

plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_color_by_insert.png"))
plt.close(fig)
print("  → genome_tracks_color_by_insert.png")

# =========================================================================
# 2. Uniform gray coloring (genomeview VCF example style)
# =========================================================================
print("[2] color_fn — uniform gray backdrop for variant display")

region_zoom = GenomicInterval(chrom, start, start + 1500)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

track_gray_illumina = AlignmentsTrack(
    filepath=ILLUMINA_BAM, is_paired=True,
    type="pileup", name="Illumina (gray)",
    color_fn=lambda read: "lightgray",
)
track_gray_illumina.draw(ax1, region_zoom)
ax1.set_title("Illumina — color_fn = lightgray", fontsize=10, loc="left")

if os.path.exists(PACBIO_BAM):
    track_gray_pacbio = AlignmentsTrack(
        filepath=PACBIO_BAM,
        type="pileup", name="PacBio (gray)",
        min_indel_size=10,
        color_fn=lambda read: "lightgray",
    )
    track_gray_pacbio.draw(ax2, region_zoom)
    ax2.set_title("PacBio — color_fn = lightgray, min_indel_size=10",
                  fontsize=10, loc="left")
else:
    ax2.text(0.5, 0.5, "PacBio BAM not found", ha="center", va="center",
             transform=ax2.transAxes)

plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_color_gray.png"))
plt.close(fig)
print("  → genome_tracks_color_gray.png")

# =========================================================================
# 3. Color by mapping quality
# =========================================================================
print("[3] color_fn — color by mapping quality")


def color_by_mapq(read):
    """Color reads by mapping quality (MAPQ).

    - Dark blue: high MAPQ (≥40)
    - Light blue: medium MAPQ (20-39)
    - Orange: low MAPQ (<20)
    """
    mapq = read.mapping_quality
    if mapq >= 40:
        return "#1565C0"  # dark blue
    elif mapq >= 20:
        return "#64B5F6"  # light blue
    else:
        return "#FF9800"  # orange


fig, ax = plt.subplots(figsize=(14, 5))

track_mapq = AlignmentsTrack(
    filepath=ILLUMINA_BAM, is_paired=True,
    type="pileup", name="MAPQ Coloring",
    color_fn=color_by_mapq,
)
track_mapq.draw(ax, region_zoom)
ax.set_title("Pileup — color_fn by Mapping Quality (MAPQ)",
             fontsize=10, loc="left")

legend_elements = [
    Patch(facecolor="#1565C0", label="MAPQ ≥ 40"),
    Patch(facecolor="#64B5F6", label="MAPQ 20–39"),
    Patch(facecolor="#FF9800", label="MAPQ < 20"),
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=7,
          framealpha=0.8, title="Mapping Quality")

plt.tight_layout()
save_figure([ax], os.path.join(FIG_DIR, "genome_tracks_color_by_mapq.png"))
plt.close(fig)
print("  → genome_tracks_color_by_mapq.png")

# =========================================================================
# 4. Combined: coverage + pileup with color_fn
# =========================================================================
print("[4] Combined coverage + pileup with color_fn")

axes = plot_tracks(
    [
        GenomeAxisTrack(),
        BAMCoverageTrack(filepath=ILLUMINA_BAM, name="Illumina Coverage"),
        AlignmentsTrack(
            filepath=ILLUMINA_BAM, is_paired=True,
            type="pileup", name="Insert Size",
            color_fn=color_by_insert_size,
        ),
    ],
    region=region,
    figsize=(14, 8),
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_color_fn_combined.png"))
plt.close("all")
print("  → genome_tracks_color_fn_combined.png")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All custom color figures generated in:", FIG_DIR)
