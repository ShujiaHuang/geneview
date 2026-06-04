"""File I/O readers -- load real genomic data and visualize as genome tracks.

Demonstrates all available I/O readers with real test data shipped in
``examples/data/genome_tracks/``:

    - ``read_bed``       -> test.bed (chr7:127M, BED9 with RGB colors)
    - ``read_bedgraph``  -> test.bedGraph (chr19:49.3M, signed values)
    - ``read_bigwig``    -> test.bw (chr19:49.3M, BigWig)
    - ``read_bam_coverage`` -> test.bam (chr1:189.9M, hg19 alignments)
    - ``read_gff``       -> test.gtf and test.gff3 (chr1:67M, SGIP1 gene)
    - ``read_auto``      -> auto-detect format by extension

Optional dependencies (``pyBigWig``, ``pysam``) are handled gracefully --
the script falls back to synthetic data when they are unavailable.

Run:  python examples/scripts/genome_tracks_io.py
"""
import os
import warnings
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from geneview.genometracks import (
    GenomeAxisTrack,
    AnnotationTrack,
    GeneRegionTrack,
    DataTrack,
    GenomicInterval,
    plot_tracks,
    read_bed,
    read_bedgraph,
    read_gff,
    read_auto,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

gtrack = GenomeAxisTrack()
figs = []
labels = []


# ---------------------------------------------------------------------------
# 1. read_bed -- real BED9 file with RGB colors
# ---------------------------------------------------------------------------
bed_path = os.path.join(DATA_DIR, "test.bed")
bed_data = read_bed(bed_path)
print(f"[INFO] read_bed('{os.path.basename(bed_path)}') -> {len(bed_data)} features")

region_bed = GenomicInterval("chr7", 127_471_000, 127_482_000)
atrack = AnnotationTrack(bed_data, name="BED (test.bed)")
axes = plot_tracks(
    [gtrack, atrack], region=region_bed,
    title="read_bed: test.bed (chr7:127.47M)", figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_bed")


# ---------------------------------------------------------------------------
# 2. read_bedgraph -- real bedGraph with signed values
# ---------------------------------------------------------------------------
bg_path = os.path.join(DATA_DIR, "test.bedGraph")
bg_data = read_bedgraph(bg_path)
print(f"[INFO] read_bedgraph('{os.path.basename(bg_path)}') -> {len(bg_data)} bins")

region_bg = GenomicInterval("chr19", 49_302_000, 49_305_000)
dtrack_bg = DataTrack(bg_data, type="histogram", name="bedGraph (test.bedGraph)")
axes = plot_tracks(
    [gtrack, dtrack_bg], region=region_bg,
    title="read_bedgraph: test.bedGraph (chr19:49.3M)", figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_bedgraph")

# Line plot variant
dtrack_bg_line = DataTrack(bg_data, type="line", name="bedGraph line")
axes = plot_tracks(
    [gtrack, dtrack_bg_line], region=region_bg,
    title="read_bedgraph: line plot", figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_bedgraph_line")


# ---------------------------------------------------------------------------
# 3. read_bigwig -- real BigWig (requires pyBigWig)
# ---------------------------------------------------------------------------
bw_path = os.path.join(DATA_DIR, "test.bw")
try:
    from geneview.genometracks import read_bigwig

    bw_data = read_bigwig(bw_path)
    print(f"[INFO] read_bigwig('{os.path.basename(bw_path)}') -> {len(bw_data)} bins")

    dtrack_bw = DataTrack(bw_data, type="line", name="BigWig (test.bw)")
    axes = plot_tracks(
        [gtrack, dtrack_bw], region=region_bg,
        title="read_bigwig: test.bw (chr19:49.3M)", figsize=(12, 4),
    )
    figs.append(axes[0].figure)
    labels.append("io_bigwig")
except ImportError:
    print("[WARN] pyBigWig not installed -- skipping BigWig example.")
    print("       Install with:  pip install pyBigWig")


# ---------------------------------------------------------------------------
# 4. read_bam_coverage -- real BAM file (requires pysam)
# ---------------------------------------------------------------------------
bam_path = os.path.join(DATA_DIR, "test.bam")
try:
    from geneview.genometracks import read_bam_coverage

    region_bam = GenomicInterval("chr1", 189_891_000, 189_900_000)
    bam_data = read_bam_coverage(bam_path, region=region_bam, bins=200)
    print(f"[INFO] read_bam_coverage('{os.path.basename(bam_path)}') -> {len(bam_data)} bins")

    dtrack_bam = DataTrack(bam_data, type="histogram", name="BAM (test.bam)")
    axes = plot_tracks(
        [gtrack, dtrack_bam], region=region_bam,
        title="read_bam_coverage: test.bam (chr1:189.89M)", figsize=(12, 4),
    )
    figs.append(axes[0].figure)
    labels.append("io_bam")
except ImportError:
    print("[WARN] pysam not installed -- skipping BAM example.")
    print("       Install with:  pip install pysam")


# ---------------------------------------------------------------------------
# 5. read_gff -- real GTF and GFF3 files
# ---------------------------------------------------------------------------
gtf_path = os.path.join(DATA_DIR, "test.gtf")
gtf_data = read_gff(gtf_path)
print(f"[INFO] read_gff('{os.path.basename(gtf_path)}') -> {len(gtf_data)} features")

gff3_path = os.path.join(DATA_DIR, "test.gff3")
gff3_data = read_gff(gff3_path)
print(f"[INFO] read_gff('{os.path.basename(gff3_path)}') -> {len(gff3_data)} features")

region_gene = GenomicInterval("chr1", 66_999_000, 67_200_000)

# GTF -> GeneRegionTrack
grtrack_gtf = GeneRegionTrack(gtf_data, name="GTF (test.gtf)")
axes = plot_tracks(
    [gtrack, grtrack_gtf], region=region_gene,
    title="read_gff: test.gtf (chr1:67M, SGIP1)", figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_gtf")

# GFF3 -> GeneRegionTrack
grtrack_gff3 = GeneRegionTrack(gff3_data, name="GFF3 (test.gff3)")
axes = plot_tracks(
    [gtrack, grtrack_gff3], region=region_gene,
    title="read_gff: test.gff3 (chr1:67M, SGIP1)", figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_gff3")

# Side-by-side comparison
axes = plot_tracks(
    [gtrack, grtrack_gtf, grtrack_gff3],
    region=region_gene,
    title="GTF vs GFF3 — same gene, different formats",
    figsize=(12, 6),
)
figs.append(axes[0].figure)
labels.append("io_gtf_vs_gff3")


# ---------------------------------------------------------------------------
# 6. read_auto -- format auto-detection
# ---------------------------------------------------------------------------
print("\n[INFO] read_auto dispatch demo:")
auto_files = [
    ("test.bed", "chr7", 127_471_000, 127_482_000),
    ("test.bedGraph", "chr19", 49_302_000, 49_305_000),
    ("test.gtf", "chr1", 66_999_000, 67_200_000),
    ("test.gff3", "chr1", 66_999_000, 67_200_000),
]

for fname, chrom, start, end in auto_files:
    fpath = os.path.join(DATA_DIR, fname)
    df = read_auto(fpath)
    print(f"  read_auto('{fname}') -> {df.shape[0]} rows, {df.shape[1]} cols")


# ---------------------------------------------------------------------------
# 7. Combined: BED annotation + bedGraph data + GTF gene model
# ---------------------------------------------------------------------------
# Multi-panel figure combining different real data sources
atrack_bed = AnnotationTrack(bed_data, name="test.bed")
dtrack_bg2 = DataTrack(bg_data, type="polygon", name="test.bedGraph")

# Use bedGraph data for the combined panel since it's on chr19
axes = plot_tracks(
    [gtrack, dtrack_bg2],
    region=region_bg,
    title="Combined: bedGraph polygon",
    figsize=(12, 4),
)
figs.append(axes[0].figure)
labels.append("io_combined")


# ---------------------------------------------------------------------------
# Save all figures
# ---------------------------------------------------------------------------
for label, fig in zip(labels, figs):
    path = os.path.join(OUT_DIR, f"genome_tracks_{label}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {path}")

plt.show()
print(f"\n[DONE] Generated {len(figs)} figures from real test data.")
