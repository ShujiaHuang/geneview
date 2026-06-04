"""Advanced genome tracks features -- OverlayTrack, IdeogramTrack, direction indicators, etc.

Showcases features including:
- OverlayTrack: overlay multiple DataTracks on the same axes
- IdeogramTrack: chromosome ideogram with cytoband coloring (Gviz-style)
- GenomeAxisTrack direction indicators (add53/add35)
- AnnotationTrack fixedArrow and smallArrow shapes
- AnnotationTrack group_annotation with connecting lines
- plot_tracks show_title=False and reverse_strand=True
- Fractional extend_left/extend_right
- GenomeAxisTrack exponent and ticks_at for axis label control
- HighlightTrack in_background for z-order control
- Constructor kwargs as display parameters (Gviz convention)

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
    OverlayTrack, HighlightTrack, IdeogramTrack, GenomicInterval, plot_tracks,
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
                     col="#0080FF", alpha=0.9)
dtrack_b = DataTrack(overlay_data_b, type="line", name="Sample B",
                     col="#DC0000", alpha=0.7)

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
                         show_label=True)
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
dtrack = DataTrack(cov_data, type="histogram", name="Coverage")

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

# ---------------------------------------------------------------------------
# 6. IdeogramTrack — chromosome ideogram (Gviz-style)
# ---------------------------------------------------------------------------
# Create cytoband data for chr7 (synthetic data mimicking UCSC format)
bands_chr7 = pd.DataFrame({
    "chrom": ["chr7"] * 18,
    "chromStart": [
        0, 7100000, 11500000, 19000000, 24900000, 27400000,
        31600000, 35800000, 40800000, 45600000,
        58800000, 61000000, 62900000, 72700000, 83200000,
        92100000, 103400000, 127500000,
    ],
    "chromEnd": [
        7100000, 11500000, 19000000, 24900000, 27400000, 31600000,
        35800000, 40800000, 45600000, 58800000,
        61000000, 62900000, 72700000, 83200000, 92100000,
        103400000, 127500000, 159138663,
    ],
    "name": [
        "p22.3", "p22.2", "p22.1", "p21.3", "p21.2", "p21.1",
        "p15.3", "p15.2", "p15.1", "p14.3",
        "p14.1", "p13.3", "p13.1", "p12.3", "p12.1",
        "q11.21", "q11.23", "q21.11",
    ],
    "gieStain": [
        "gneg", "gpos50", "gneg", "gpos75", "gneg", "gpos25",
        "gneg", "gpos50", "gneg", "gpos75",
        "acen", "acen", "gneg", "gpos75", "gneg",
        "gpos50", "gneg", "gpos100",
    ],
})

ideo = IdeogramTrack(bands_chr7, chromosome="chr7", show_band_id=True)
axes6 = plot_tracks(
    [ideo, gtrack, atrack, dtrack],
    region=region,
    title="IdeogramTrack + Annotation + Coverage",
    figsize=(14, 6),
)
fig6 = axes6[0].figure
fig6.savefig(os.path.join(OUT_DIR, "genome_tracks_ideogram.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_ideogram.png")

# Ideogram with circle centromere
ideo_circle = IdeogramTrack(
    bands_chr7, chromosome="chr7",
    centromere_shape="circle", show_band_id=True,
    name="chr7",
)
axes6b = plot_tracks(
    [ideo_circle, gtrack],
    region=region,
    title="IdeogramTrack (circle centromere)",
    figsize=(14, 3),
)
fig6b = axes6b[0].figure
fig6b.savefig(os.path.join(OUT_DIR, "genome_tracks_ideogram_circle.png"),
              dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_ideogram_circle.png")

# ---------------------------------------------------------------------------
# 7. GenomeAxisTrack with explicit ticks_at and forced exponent
# ---------------------------------------------------------------------------
tick_positions = [26_500_000, 26_550_000, 26_600_000, 26_650_000, 26_700_000]
gtrack_exp = GenomeAxisTrack(exponent=6, name="Forced Mb")
gtrack_ticks = GenomeAxisTrack(ticks_at=tick_positions, name="Custom Ticks")

axes7 = plot_tracks(
    [GenomeAxisTrack(name="Default"), gtrack_exp, gtrack_ticks],
    region=region,
    title="GenomeAxisTrack: exponent and ticks_at",
    figsize=(14, 5),
)
fig7 = axes7[0].figure
fig7.savefig(os.path.join(OUT_DIR, "genome_tracks_axis_enhanced.png"),
             dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_axis_enhanced.png")

# ---------------------------------------------------------------------------
# 8. HighlightTrack with in_background z-order control
# ---------------------------------------------------------------------------
hl_regions = pd.DataFrame({
    "chrom": ["chr7", "chr7"],
    "start": [26_520_000, 26_620_000],
    "end":   [26_560_000, 26_660_000],
})

atrack_hl = AnnotationTrack(cpg_data, name="CpG Islands")
htrack_bg = HighlightTrack(
    [atrack_hl], regions=hl_regions,
    fill="yellow", alpha=0.3, in_background=True, name="Background",
)
htrack_fg = HighlightTrack(
    [AnnotationTrack(cpg_data, name="CpG Islands")],
    regions=hl_regions,
    fill="yellow", alpha=0.3, in_background=False, name="Foreground",
)

axes8a = plot_tracks(
    [GenomeAxisTrack(), htrack_bg], region=region,
    title="HighlightTrack: in_background=True (behind)", figsize=(14, 4),
)
fig8a = axes8a[0].figure
fig8a.savefig(os.path.join(OUT_DIR, "genome_tracks_highlight_background.png"),
              dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_highlight_background.png")

axes8b = plot_tracks(
    [GenomeAxisTrack(), htrack_fg], region=region,
    title="HighlightTrack: in_background=False (foreground)", figsize=(14, 4),
)
fig8b = axes8b[0].figure
fig8b.savefig(os.path.join(OUT_DIR, "genome_tracks_highlight_foreground.png"),
              dpi=150, bbox_inches="tight")
print(f"[INFO] Saved genome_tracks_highlight_foreground.png")

plt.show()
print("[INFO] All advanced example figures saved.")
