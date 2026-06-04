"""GeneRegionTrack -- visualize gene models (exons, UTRs, introns).

Demonstrates loading a GTF file, creating a GeneRegionTrack, and exploring
transcript collapsing modes, drawing styles, and display options.

Run:  python examples/scripts/genome_tracks_gene_region.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, GeneRegionTrack, AnnotationTrack,
    GenomicInterval, plot_tracks, read_gff,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
GTF_FILE = os.path.join(DATA_DIR, "gene_models.gtf")

# Load gene model data from GTF
gene_data = read_gff(GTF_FILE)
region = GenomicInterval("chr7", 26_490_000, 26_720_000)
gtrack = GenomeAxisTrack()

# ---------------------------------------------------------------------------
# 1. All transcripts (default)
# ---------------------------------------------------------------------------
grtrack_all = GeneRegionTrack(gene_data, name="All Transcripts")
figs = []

axes1 = plot_tracks([gtrack, grtrack_all], region=region,
                    title="All Transcripts", figsize=(14, 5))
figs.append(axes1[0].figure)

# ---------------------------------------------------------------------------
# 2. Collapse to gene level
# ---------------------------------------------------------------------------
grtrack_gene = GeneRegionTrack(gene_data, name="Gene Level",
                               collapse_transcripts="gene")
axes2 = plot_tracks([gtrack, grtrack_gene], region=region,
                    title="Collapsed to Gene Level", figsize=(14, 4))
figs.append(axes2[0].figure)

# ---------------------------------------------------------------------------
# 3. Show only the longest transcript
# ---------------------------------------------------------------------------
grtrack_long = GeneRegionTrack(gene_data, name="Longest",
                               collapse_transcripts="longest")
axes3 = plot_tracks([gtrack, grtrack_long], region=region,
                    title="Longest Transcript Only", figsize=(14, 4))
figs.append(axes3[0].figure)

# ---------------------------------------------------------------------------
# 4. Meta-transcript (union of all exons)
# ---------------------------------------------------------------------------
grtrack_meta = GeneRegionTrack(gene_data, name="Meta-transcript",
                               collapse_transcripts="meta")
axes4 = plot_tracks([gtrack, grtrack_meta], region=region,
                    title="Meta-transcript (union of all exons)",
                    figsize=(14, 4))
figs.append(axes4[0].figure)

# ---------------------------------------------------------------------------
# 5. Combined: axis + gene region + annotation (CpG islands)
# ---------------------------------------------------------------------------
from geneview.genometracks import read_bed
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
atrack = AnnotationTrack(cpg_data, name="CpG Islands")

axes5 = plot_tracks(
    [gtrack, atrack, grtrack_all],
    region=region,
    title="Gene Models + CpG Islands",
    figsize=(14, 6),
)
figs.append(axes5[0].figure)

# ---------------------------------------------------------------------------
# 6. Drawing styles — 4-panel comparison (UCSC, flybase, tssarrow, exonarrows)
# ---------------------------------------------------------------------------
styles = ["UCSC", "flybase", "tssarrow", "exonarrows"]
fig6, axs6 = plt.subplots(4, 1, figsize=(14, 10))
plt.close(fig6)  # we'll fill axes manually then re-show

for i, style in enumerate(styles):
    grtrack_s = GeneRegionTrack(
        gene_data, name=style,
        style=style,
        collapse_transcripts="longest",
    )
    plot_tracks(
        [grtrack_s], region=region,
        ax=axs6[i], add=True,
        title=f"style='{style}'",
    )

axs6[0].set_title("GeneRegionTrack Drawing Styles (longest transcript)", fontsize=12)
fig6.tight_layout()
figs.append(fig6)

# ---------------------------------------------------------------------------
# 7. Style customisation — flybase + exonarrows with display_params
# ---------------------------------------------------------------------------
grtrack_fb = GeneRegionTrack(
    gene_data, name="Flybase (custom)",
    style="flybase",
    collapse_transcripts="longest",
    display_params={
        "color_utr": "#AAAAAA",
        "height_utr": 0.5,
        "color_backbone": "#555555",
    },
)
grtrack_ea = GeneRegionTrack(
    gene_data, name="Exonarrows (custom)",
    style="exonarrows",
    collapse_transcripts="longest",
    display_params={
        "height_intron": 0.3,
        "arrow_interval": 3,
        "color_arrow": "white",
    },
)
axes7 = plot_tracks(
    [gtrack, grtrack_fb, grtrack_ea],
    region=region,
    title="Customised Drawing Styles",
    figsize=(14, 7),
)
figs.append(axes7[0].figure)

# ---------------------------------------------------------------------------
# Save all figures
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(out_dir, exist_ok=True)
for i, fig in enumerate(figs, 1):
    path = os.path.join(out_dir, f"genome_tracks_gene_region_{i}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {path}")

plt.show()
