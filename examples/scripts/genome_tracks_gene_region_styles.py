"""
GeneRegionTrack Drawing Styles
================================

Demonstrates the four gene drawing styles ported from pyGenomeTracks:

- **UCSC** (default): thick CDS blocks, thin UTR blocks, intron chevron arrows
- **flybase**: backbone line, last exon as filled directional arrow
- **tssarrow**: vertical line + arrow at TSS, half-height exon boxes, L-shaped introns
- **exonarrows**: full-height exon boxes with directional arrows inside, filled intron connectors

Also shows customisation via display params: ``color_utr``, ``color_backbone``,
``height_utr``, ``height_intron``, ``arrow_interval``, ``arrowhead_fraction``.

Run:  python examples/scripts/genome_tracks_gene_region_styles.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, GeneRegionTrack,
    GenomicInterval, read_gff,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
os.makedirs(FIG_DIR, exist_ok=True)

# Load gene model data
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
region = GenomicInterval("chr7", 26_490_000, 26_720_000)
gtrack = GenomeAxisTrack()

# ---------------------------------------------------------------------------
# 1. Side-by-side comparison of all 4 styles (longest transcript)
# ---------------------------------------------------------------------------
styles = ["UCSC", "flybase", "tssarrow", "exonarrows"]
fig, axes = plt.subplots(len(styles) + 1, 1, figsize=(14, 2.8 * (len(styles) + 1)))
plt.close(fig)  # we'll fill axes manually

for i, style in enumerate(styles):
    grtrack = GeneRegionTrack(
        gene_data, name=style,
        style=style,
        collapse_transcripts="longest",
    )
    axes_list = plot_tracks(
        [gtrack, grtrack],
        region=region, figsize=(14, 3),
        title=f"GeneRegionTrack  —  style='{style}'",
    )
    fig_style = axes_list[0].figure
    fig_style.savefig(
        os.path.join(FIG_DIR, f"genome_tracks_gene_region_style_{style.lower()}.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close(fig_style)
    print(f"[INFO] Saved genome_tracks_gene_region_style_{style.lower()}.png")

# ---------------------------------------------------------------------------
# 2. All four styles in a single combined figure (4 panels)
# ---------------------------------------------------------------------------
fig_combined, axs = plt.subplots(4, 1, figsize=(14, 10))
plt.close(fig_combined)

for i, style in enumerate(styles):
    grtrack = GeneRegionTrack(
        gene_data, name=style,
        style=style,
        collapse_transcripts="longest",
    )
    plot_tracks(
        [grtrack], region=region,
        ax=axs[i], add=True,
        title=f"style='{style}'",
    )

axs[0].set_title("GeneRegionTrack Drawing Styles (longest transcript)", fontsize=12)
fig_combined.tight_layout()
fig_combined.savefig(
    os.path.join(FIG_DIR, "genome_tracks_gene_region_styles_combined.png"),
    dpi=150, bbox_inches="tight",
)
plt.close(fig_combined)
print("[INFO] Saved genome_tracks_gene_region_styles_combined.png")

# ---------------------------------------------------------------------------
# 3. Flybase style with custom UTR colour and height
# ---------------------------------------------------------------------------
grtrack_flybase_custom = GeneRegionTrack(
    gene_data, name="Flybase (custom)",
    style="flybase",
    collapse_transcripts="longest",
    display_params={
        "color_utr": "#AAAAAA",
        "height_utr": 0.5,
        "color_backbone": "#555555",
    },
)
axes_fb = plot_tracks(
    [GenomeAxisTrack(), grtrack_flybase_custom],
    region=region, figsize=(14, 4),
    title="Flybase style  —  color_utr=grey, height_utr=0.5, color_backbone=#555",
)
fig_fb = axes_fb[0].figure
fig_fb.savefig(
    os.path.join(FIG_DIR, "genome_tracks_gene_region_style_flybase_custom.png"),
    dpi=150, bbox_inches="tight",
)
plt.close(fig_fb)
print("[INFO] Saved genome_tracks_gene_region_style_flybase_custom.png")

# ---------------------------------------------------------------------------
# 4. Exonarrows style with custom intron height and arrow interval
# ---------------------------------------------------------------------------
grtrack_ea_custom = GeneRegionTrack(
    gene_data, name="Exonarrows (custom)",
    style="exonarrows",
    collapse_transcripts="longest",
    display_params={
        "height_intron": 0.3,
        "arrow_interval": 3,
        "color_arrow": "white",
    },
)
axes_ea = plot_tracks(
    [GenomeAxisTrack(), grtrack_ea_custom],
    region=region, figsize=(14, 4),
    title="Exonarrows style  —  height_intron=0.3, arrow_interval=3, color_arrow=white",
)
fig_ea = axes_ea[0].figure
fig_ea.savefig(
    os.path.join(FIG_DIR, "genome_tracks_gene_region_style_exonarrows_custom.png"),
    dpi=150, bbox_inches="tight",
)
plt.close(fig_ea)
print("[INFO] Saved genome_tracks_gene_region_style_exonarrows_custom.png")

# ---------------------------------------------------------------------------
# 5. TSS arrow style with reverse strand gene
# ---------------------------------------------------------------------------
# Use a region known to have reverse-strand genes, or just show the style
grtrack_tss = GeneRegionTrack(
    gene_data, name="TSS Arrow",
    style="tssarrow",
    collapse_transcripts="longest",
    display_params={
        "color_utr": "#BBBBBB",
        "height_utr": 0.6,
    },
)
axes_tss = plot_tracks(
    [GenomeAxisTrack(), grtrack_tss],
    region=region, figsize=(14, 4),
    title="TSS Arrow style  —  color_utr=#BBB, height_utr=0.6",
)
fig_tss = axes_tss[0].figure
fig_tss.savefig(
    os.path.join(FIG_DIR, "genome_tracks_gene_region_style_tssarrow_custom.png"),
    dpi=150, bbox_inches="tight",
)
plt.close(fig_tss)
print("[INFO] Saved genome_tracks_gene_region_style_tssarrow_custom.png")

# ---------------------------------------------------------------------------
# 6. All transcripts (not collapsed) with flybase style
# ---------------------------------------------------------------------------
grtrack_flybase_all = GeneRegionTrack(
    gene_data, name="Flybase (all transcripts)",
    style="flybase",
    collapse_transcripts=False,
)
axes_fba = plot_tracks(
    [GenomeAxisTrack(), grtrack_flybase_all],
    region=region, figsize=(14, 5),
    title="Flybase style  —  all transcripts (collapse_transcripts=False)",
)
fig_fba = axes_fba[0].figure
fig_fba.savefig(
    os.path.join(FIG_DIR, "genome_tracks_gene_region_style_flybase_all.png"),
    dpi=150, bbox_inches="tight",
)
plt.close(fig_fba)
print("[INFO] Saved genome_tracks_gene_region_style_flybase_all.png")

print("\nAll GeneRegionTrack style examples saved.")
