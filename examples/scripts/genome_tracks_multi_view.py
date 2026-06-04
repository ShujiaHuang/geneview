"""Multi-View Layouts Example
===============================

Demonstrates the new multi-view layout functions ported from genomeview:

    - ``plot_tracks_grid`` — render independent views side-by-side (grid)
    - ``plot_tracks_multi`` — stack sections from different regions vertically

Useful for comparing loci, samples, or conditions in a single figure.

Run:  python examples/scripts/genome_tracks_multi_view.py
"""
import os
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack,
    AnnotationTrack,
    DataTrack,
    GenomicInterval,
    plot_tracks,
    plot_tracks_grid,
    plot_tracks_multi,
    read_bed,
    read_bedgraph,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
os.makedirs(FIG_DIR, exist_ok=True)

gtrack = GenomeAxisTrack()

# ---------------------------------------------------------------------------
# Load shared annotation data
# ---------------------------------------------------------------------------
BED_FILE = os.path.join(DATA_DIR, "cpg_islands.bed")
cpg_data = read_bed(BED_FILE)


# ---------------------------------------------------------------------------
# 1. plot_tracks_grid — compare two regions side-by-side
# ---------------------------------------------------------------------------
region_left = GenomicInterval("chr7", 26_500_000, 26_800_000)
region_right = GenomicInterval("chr7", 27_000_000, 27_300_000)

view_left = [
    gtrack,
    AnnotationTrack(cpg_data, name="CpG Islands", stacking="squish"),
]
view_right = [
    GenomeAxisTrack(),
    AnnotationTrack(cpg_data, name="CpG Islands", stacking="squish"),
]

axes = plot_tracks_grid(
    [view_left, view_right],
    regions=[region_left, region_right],
    columns=2,
    figsize=(20, 6),
    title="Side-by-Side Comparison: Two Regions on chr7",
)
fig = plt.gcf()
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_multi_view_grid.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: genome_tracks_multi_view_grid.png")


# ---------------------------------------------------------------------------
# 2. plot_tracks_grid — 2×2 grid with four views
# ---------------------------------------------------------------------------
regions_4 = [
    GenomicInterval("chr7", 26_000_000, 26_300_000),
    GenomicInterval("chr7", 26_500_000, 26_800_000),
    GenomicInterval("chr7", 27_000_000, 27_300_000),
    GenomicInterval("chr7", 27_500_000, 27_800_000),
]

views_4 = []
for i, reg in enumerate(regions_4):
    view = [
        GenomeAxisTrack(),
        AnnotationTrack(cpg_data, name=f"Region {i+1}", stacking="squish"),
    ]
    views_4.append(view)

axes = plot_tracks_grid(
    views_4,
    regions=regions_4,
    columns=2,
    figsize=(20, 12),
    title="2×2 Grid: Four chr7 Regions",
)
fig = plt.gcf()
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_multi_view_grid2x2.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: genome_tracks_multi_view_grid2x2.png")


# ---------------------------------------------------------------------------
# 3. plot_tracks_multi — stacked sections from different regions
# ---------------------------------------------------------------------------
section1_tracks = [
    gtrack,
    AnnotationTrack(cpg_data, name="Region A", stacking="squish"),
]
section2_tracks = [
    GenomeAxisTrack(),
    AnnotationTrack(cpg_data, name="Region B", stacking="squish"),
]

sections = [
    (section1_tracks, GenomicInterval("chr7", 26_500_000, 26_800_000)),
    (section2_tracks, GenomicInterval("chr7", 27_000_000, 27_300_000)),
]

axes = plot_tracks_multi(
    sections,
    title="Multi-Region Stacked View",
)
fig = plt.gcf()
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_multi_view_stacked.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: genome_tracks_multi_view_stacked.png")


# ---------------------------------------------------------------------------
# 4. plot_tracks_grid with DataTrack and AnnotationTrack combined
# ---------------------------------------------------------------------------
BEDGRAPH_FILE = os.path.join(DATA_DIR, "coverage.bedgraph")
try:
    bg_data = read_bedgraph(BEDGRAPH_FILE)
    has_bedgraph = True
except Exception:
    has_bedgraph = False

if has_bedgraph:
    region_d = GenomicInterval("chr19", 49_300_000, 49_500_000)

    view_data = [
        gtrack,
        DataTrack(bg_data, name="Coverage", type="line"),
        AnnotationTrack(cpg_data, name="CpG", stacking="squish"),
    ]
    view_annot = [
        GenomeAxisTrack(),
        AnnotationTrack(cpg_data, name="CpG Full", stacking="squish"),
    ]

    axes = plot_tracks_grid(
        [view_data, view_annot],
        regions=[region_d, GenomicInterval("chr7", 26_500_000, 26_800_000)],
        columns=2,
        figsize=(20, 8),
        title="Mixed Track Types in Grid Layout",
    )
    fig = plt.gcf()
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_multi_view_mixed.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_multi_view_mixed.png")

print("[INFO] Multi-view examples complete.")
