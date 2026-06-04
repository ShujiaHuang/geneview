"""
Color Schemes and plot_tracks Enhancements
==========================================

Demonstrates apply_scheme, and the new plot_tracks parameters:
cex (font scaling), add (existing axes), ylim (global y-limits),
and scheme (color scheme).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, GeneRegionTrack, AnnotationTrack,
    DataTrack, GenomicInterval, apply_scheme, read_gff, read_bed,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# 1. Color scheme: "genes" — distinct color per gene
# ---------------------------------------------------------------------------
grtrack = GeneRegionTrack(gene_data, name="Genes (scheme=genes)",
                          collapse_transcripts="longest")
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack],
    region=region, figsize=(14, 4),
    scheme="genes",
    title="Color Scheme: genes",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_scheme_genes.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Color scheme: "transcripts" — distinct color per transcript
# ---------------------------------------------------------------------------
grtrack2 = GeneRegionTrack(gene_data, name="Transcripts")
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack2],
    region=region, figsize=(14, 5),
    scheme="transcripts",
    title="Color Scheme: transcripts",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_scheme_transcripts.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Font scaling (cex)
# ---------------------------------------------------------------------------
atrack = AnnotationTrack(cpg_data, name="CpG Islands")
grtrack3 = GeneRegionTrack(gene_data, name="Genes", collapse_transcripts="longest")
axes = plot_tracks(
    [GenomeAxisTrack(), atrack, grtrack3],
    region=region, figsize=(14, 5),
    cex=1.5,
    title="Font Scaling: cex=1.5",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_cex.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Global y-axis limits (ylim) for DataTrack
# ---------------------------------------------------------------------------
rng = np.random.RandomState(42)
n = 50
starts = np.linspace(26_490_000, 26_720_000, n, dtype=int)
d1 = pd.DataFrame({"chrom": ["chr7"] * n, "start": starts,
                    "end": starts + 5000, "value": rng.randn(n).cumsum()})
d2 = pd.DataFrame({"chrom": ["chr7"] * n, "start": starts,
                    "end": starts + 5000, "value": rng.randn(n).cumsum() * 3})

dt1 = DataTrack(d1, type="line", name="Signal A")
dt2 = DataTrack(d2, type="line", name="Signal B")

axes = plot_tracks(
    [GenomeAxisTrack(), dt1, dt2],
    region=region, figsize=(14, 5),
    ylim=(-10, 10),
    title="Global ylim=(-10, 10)",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_ylim.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. Plot into existing axes (add)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 2))
dtrack = DataTrack(d1, type="histogram", name="Added Track")
axes = plot_tracks(
    [dtrack], region=region, add=True, ax=ax,
    title="Plot into existing axes (add=True)",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_add.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

print("All color scheme and plot_tracks enhancement examples saved.")
