"""Convenience API Example — visualize_files()
================================================

Demonstrates the ``visualize_files()`` helper that auto-detects file types
from extensions and builds the appropriate tracks automatically.

This is the quickest way to visualize a set of genomic files without
manually constructing each track type.

Run:  python examples/scripts/genome_tracks_convenience.py
"""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomicInterval,
    plot_tracks,
    visualize_files,
    match_chrom_format,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
os.makedirs(FIG_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Basic usage — pass a list of files
# ---------------------------------------------------------------------------
region = GenomicInterval("chr7", 26_500_000, 26_800_000)
bed_file = os.path.join(DATA_DIR, "cpg_islands.bed")

tracks = visualize_files(
    [bed_file],
    region=region,
)

axes = plot_tracks(tracks, region=region, figsize=(14, 4),
                   title="visualize_files() — BED File")
fig = axes[0].figure
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_convenience_bed.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: genome_tracks_convenience_bed.png")


# ---------------------------------------------------------------------------
# 2. Multiple file types at once (BED + BedGraph)
# ---------------------------------------------------------------------------
bedgraph_file = os.path.join(DATA_DIR, "coverage.bedgraph")
region_bg = GenomicInterval("chr19", 49_300_000, 49_500_000)

try:
    tracks = visualize_files(
        [bedgraph_file, bed_file],
        region=region_bg,
    )
    axes = plot_tracks(tracks, region=region_bg, figsize=(14, 6),
                       title="visualize_files() — BedGraph + BED")
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_convenience_multi.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_convenience_multi.png")
except Exception as e:
    print(f"[SKIP] BedGraph example: {e}")


# ---------------------------------------------------------------------------
# 3. BAM file with auto-detection of paired-end status
# ---------------------------------------------------------------------------
try:
    import pysam
    has_pysam = True
except ImportError:
    has_pysam = False

if has_pysam:
    bam_file = os.path.join(DATA_DIR, "test.bam")
    region_bam = GenomicInterval("chr1", 189_892_000, 189_895_000)

    tracks = visualize_files(
        [bam_file],
        region=region_bam,
    )
    axes = plot_tracks(tracks, region=region_bam, figsize=(14, 5),
                       title="visualize_files() — BAM (auto-detect PE)")
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_convenience_bam.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_convenience_bam.png")

    # -----------------------------------------------------------------------
    # 4. Dict input — use descriptive track names
    # -----------------------------------------------------------------------
    tracks = visualize_files(
        {"Control BAM": bam_file, "CpG Islands": bed_file},
        region=region_bam,
    )
    axes = plot_tracks(tracks, region=region_bam, figsize=(14, 6),
                       title="visualize_files() — Dict with Custom Names")
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_convenience_dict.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_convenience_dict.png")


# ---------------------------------------------------------------------------
# 5. Utility: match_chrom_format — handle chr-prefix mismatches
# ---------------------------------------------------------------------------
print("\n--- match_chrom_format demo ---")
print(f"  match_chrom_format('chr1', ['1', '2']) => "
      f"'{match_chrom_format('chr1', ['1', '2'])}'")
print(f"  match_chrom_format('1', ['chr1', 'chr2']) => "
      f"'{match_chrom_format('1', ['chr1', 'chr2'])}'")
print(f"  match_chrom_format('chrM', ['MT']) => "
      f"'{match_chrom_format('chrM', ['MT'])}'")
print(f"  match_chrom_format('MT', ['chrM']) => "
      f"'{match_chrom_format('MT', ['chrM'])}'")

print("\n[INFO] Convenience examples complete.")
