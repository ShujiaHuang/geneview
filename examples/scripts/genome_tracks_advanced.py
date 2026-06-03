"""Advanced genome tracks features -- OverlayTrack, direction indicators, etc.

Showcases newly added features:
- OverlayTrack: overlay multiple DataTracks on the same axes
- GenomeAxisTrack direction indicators (add53/add35)
- AnnotationTrack fixedArrow and smallArrow shapes
- AnnotationTrack group_annotation with connecting lines
- plot_tracks show_title=False and reverse_strand=True
- Fractional extend_left/extend_right

Run:  python examples/scripts/genome_tracks_advanced.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
    OverlayTrack, HighlightTrack, GenomicInterval, plot_tracks,
    read_bed, read_gff, read_bedgraph,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Load data
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
cov_data = read_bedgraph(os.path.join(DATA_DIR, "coverage.bedgraph"))
ann_data = read_bed(os.path.join(DATA_DIR, "annotations.bed"))

region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# 1. OverlayTrack — compare two signals on the same axes
# ---------------------------------------------------------------------------
rng = np.random.RandomState(42)
n = len(cov_data)
starts = cov_data["start"].values
ends = cov_data["end"].values

overlay_data_a = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts, "end": ends,
    "signal_A": rng.randn(n).cumsum() / 2,
})
overlay_data_b = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts, "end": ends,
    "signal_B": rng.randn(n).cumsum() / 2 + 3,
})

dtrack_a = DataTrack(overlay_data_a, type="line", name="Sample A",
                     display_params={"col": "#3C5488", "alpha": 0.9})
dtrack_b = DataTrack(overlay_data_b, type="line", name="Sample B",
                     display_params={"col": "#E64B35", "alpha": 0.7})

otrack = OverlayTrack([dtrack_a, dtrack_b], name="Overlay A vs B")

gtrack = GenomeAxisTrack()
axes1 = plot_tracks(
    [gtrack, otrack], region=region,
    title="OverlayTrack: Two Signals on Same Axes", figsize=(14, 5),
)
fig1 = axes1[0].figure
fig1.savefig(os.path.join(OUT_DIR, "genome_tracks_overlay.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_overlay.png")

# ---------------------------------------------------------------------------
# 2. Direction indicators (add53 / add35)
# ---------------------------------------------------------------------------
gtrack_53 = GenomeAxisTrack(little_ticks=True, add53=True)
atrack = AnnotationTrack(cpg_data, name="CpG Islands")

axes2 = plot_tracks(
    [gtrack_53, atrack], region=region,
    title="GenomeAxisTrack with 5'->3' Direction Arrow", figsize=(14, 4),
)
fig2 = axes2[0].figure
fig2.savefig(os.path.join(OUT_DIR, "genome_tracks_direction_53.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_direction_53.png")

# ---------------------------------------------------------------------------
# 3. AnnotationTrack shapes: fixedArrow and smallArrow
# ---------------------------------------------------------------------------
# Create data with strand info
arrow_data = pd.DataFrame({
    "chrom": ["chr7"] * 6,
    "start": [26_500_000, 26_530_000, 26_560_000,
              26_600_000, 26_640_000, 26_680_000],
    "end":   [26_520_000, 26_550_000, 26_580_000,
              26_620_000, 26_660_000, 26_700_000],
    "strand": ["+", "+", "-", "+", "-", "+"],
    "name": ["f1", "f2", "f3", "f4", "f5", "f6"],
})

fig_shapes, axes_shapes = plt.subplots(3, 1, figsize=(14, 8))
plt.close(fig_shapes)

for idx, shape in enumerate(["arrow", "fixedArrow", "smallArrow"]):
    at = AnnotationTrack(arrow_data, shape=shape, name=f"Shape: {shape}",
                         show_label=True,
                         display_params={"fill": "#3C5488"})
    ax_list = plot_tracks([at], region=region, figsize=(14, 2.5))
    fig_s = ax_list[0].figure
    fig_s.savefig(os.path.join(OUT_DIR, f"genome_tracks_shape_{shape}.png"),
                  dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved genome_tracks_shape_{shape}.png")

# ---------------------------------------------------------------------------
# 4. GeneRegionTrack meta collapse
# ---------------------------------------------------------------------------
gene_region = GenomicInterval("chr7", 26_490_000, 26_720_000)

grtrack_all = GeneRegionTrack(gene_data, name="All Transcripts")
grtrack_meta = GeneRegionTrack(gene_data, name="Meta-transcript",
                               collapse_transcripts="meta")

axes4 = plot_tracks(
    [GenomeAxisTrack(), grtrack_all, grtrack_meta],
    region=gene_region,
    title="All Transcripts vs Meta-transcript",
    figsize=(14, 6),
)
fig4 = axes4[0].figure
fig4.savefig(os.path.join(OUT_DIR, "genome_tracks_meta_collapse.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_meta_collapse.png")

# ---------------------------------------------------------------------------
# 5. show_title=False and reverse_strand=True
# ---------------------------------------------------------------------------
dtrack = DataTrack(cov_data, type="histogram", name="Coverage",
                   display_params={"fill": "#4DBBD5", "col": "#4DBBD5"})

# Compact layout without title panels
axes5a = plot_tracks(
    [GenomeAxisTrack(), atrack, dtrack],
    region=region,
    show_title=False,
    title="Compact Layout (no title panels)",
    figsize=(14, 5),
)
fig5a = axes5a[0].figure
fig5a.savefig(os.path.join(OUT_DIR, "genome_tracks_no_title.png"),
              dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_no_title.png")

# Reverse strand
axes5b = plot_tracks(
    [GenomeAxisTrack(add35=True), atrack, dtrack],
    region=region,
    reverse_strand=True,
    title="Reverse Strand (3' on left)",
    figsize=(14, 5),
)
fig5b = axes5b[0].figure
fig5b.savefig(os.path.join(OUT_DIR, "genome_tracks_reverse_strand.png"),
              dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_reverse_strand.png")

plt.show()
print("[INFO] All advanced example figures saved.")
