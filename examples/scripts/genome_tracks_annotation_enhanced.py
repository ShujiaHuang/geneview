"""
AnnotationTrack Enhancements Example
=====================================

Demonstrates new AnnotationTrack features: just_group, show_overplotting,
merge_groups, col_line, and from_bam classmethod.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, AnnotationTrack, GenomicInterval,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

# Create annotation data with groups
data = pd.DataFrame({
    "chrom": ["chr7"] * 6,
    "start": [1000, 1200, 2000, 2100, 3000, 3200],
    "end":   [1500, 1600, 2400, 2800, 3500, 3700],
    "strand": ["+", "+", "-", "+", "+", "-"],
    "name": ["f1", "f2", "f3", "f4", "f5", "f6"],
    "group": ["A", "A", "A", "B", "B", "B"],
})

region = GenomicInterval("chr7", 800, 4000)

# ---------------------------------------------------------------------------
# 1. Group labels with justification
# ---------------------------------------------------------------------------
track_just = AnnotationTrack(
    data, group_annotation="group",
    just_group="above",
    name="just_group=above",
    display_params={"col_line": "#888888"},
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_just],
    region=region, figsize=(12, 3),
    title="AnnotationTrack - Group Label Justification",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_annotation_just_group.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Custom connector line color
# ---------------------------------------------------------------------------
track_col = AnnotationTrack(
    data, group_annotation="group",
    col_line="#CC0000",
    name="col_line=red",
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_col],
    region=region, figsize=(12, 3),
    title="AnnotationTrack - Custom Connector Lines",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_annotation_col_line.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Overplotting indicator
# ---------------------------------------------------------------------------
overlap_data = pd.DataFrame({
    "chrom": ["chr7"] * 8,
    "start": [1000, 1050, 1100, 1150, 2000, 2050, 2100, 2150],
    "end":   [1400, 1450, 1500, 1550, 2400, 2450, 2500, 2550],
    "name": [f"ov_{i}" for i in range(8)],
})

track_over = AnnotationTrack(
    overlap_data,
    show_overplotting=True,
    name="Overplotting",
    stacking="squish",
)
region_over = GenomicInterval("chr7", 900, 2600)
axes = plot_tracks(
    [GenomeAxisTrack(), track_over],
    region=region_over, figsize=(12, 3),
    title="AnnotationTrack - Overplotting Indicator",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_annotation_overplotting.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Merge groups
# ---------------------------------------------------------------------------
track_merge = AnnotationTrack(
    data, group_annotation="group",
    merge_groups=True,
    name="merge_groups",
    stacking="squish",
    display_params={"col_line": "#666666"},
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_merge],
    region=region, figsize=(12, 3),
    title="AnnotationTrack - Merge Groups",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_annotation_merge_groups.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. from_bam classmethod (requires pysam)
# ---------------------------------------------------------------------------
try:
    import pysam
    bam_path = os.path.join(DATA_DIR, "gapped.bam")
    bam_region = GenomicInterval("chr12", 2966840, 2966920)

    track_from_bam = AnnotationTrack.from_bam(
        bam_path, region=bam_region, name="Reads from BAM",
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_from_bam],
        region=bam_region, figsize=(12, 4),
        title="AnnotationTrack.from_bam()",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_annotation_from_bam.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")
    print("AnnotationTrack.from_bam example saved.")
except ImportError:
    print("from_bam example skipped: pysam not installed.")

print("All AnnotationTrack enhancement examples saved.")
