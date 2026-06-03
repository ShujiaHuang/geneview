"""GeneRegionTrack -- visualize gene models (exons, UTRs, introns).

Demonstrates loading a GTF file, creating a GeneRegionTrack, and exploring
transcript collapsing modes and display options.

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
# 4. Combined: axis + gene region + annotation (CpG islands)
# ---------------------------------------------------------------------------
from geneview.genometracks import read_bed
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
atrack = AnnotationTrack(cpg_data, name="CpG Islands")

axes4 = plot_tracks(
    [gtrack, atrack, grtrack_all],
    region=region,
    title="Gene Models + CpG Islands",
    figsize=(14, 6),
)
figs.append(axes4[0].figure)

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
