"""
GeneRegionTrack Enhancements Example
=====================================

Demonstrates new GeneRegionTrack features: exon_annotation, gene_symbols,
transcript_annotation, and collapse_transcripts="shortest".
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, GeneRegionTrack, GenomicInterval, read_gff,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

# Load gene model data
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))
region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# 1. Exon annotation labels — show exon numbers
# ---------------------------------------------------------------------------
grtrack_exon = GeneRegionTrack(
    gene_data, name="Exon Labels",
    exon_annotation="exon",
    collapse_transcripts="longest",
)
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack_exon],
    region=region, figsize=(14, 4),
    title="GeneRegionTrack - Exon Annotation Labels",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_gene_region_exon_annotation.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Gene symbols for labels
# ---------------------------------------------------------------------------
grtrack_sym = GeneRegionTrack(
    gene_data, name="Gene Symbols",
    gene_symbols=True,
    collapse_transcripts="longest",
)
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack_sym],
    region=region, figsize=(14, 4),
    title="GeneRegionTrack - Gene Symbols",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_gene_region_gene_symbols.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Transcript annotation
# ---------------------------------------------------------------------------
grtrack_tx = GeneRegionTrack(
    gene_data, name="Transcript Annotation",
    transcript_annotation="transcript",
)
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack_tx],
    region=region, figsize=(14, 5),
    title="GeneRegionTrack - Transcript Annotation",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_gene_region_transcript_annotation.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Collapse to shortest transcript
# ---------------------------------------------------------------------------
grtrack_short = GeneRegionTrack(
    gene_data, name="Shortest",
    collapse_transcripts="shortest",
)
axes = plot_tracks(
    [GenomeAxisTrack(), grtrack_short],
    region=region, figsize=(14, 4),
    title="GeneRegionTrack - Shortest Transcript",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_gene_region_shortest.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")

print("All GeneRegionTrack enhancement examples saved.")
