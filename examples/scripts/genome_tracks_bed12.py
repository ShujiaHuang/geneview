"""BED12 Transcript Rendering Example
======================================

Demonstrates BED12-aware rendering in AnnotationTrack:

    - Thick (CDS) vs thin (UTR) exon boxes
    - Intron connecting lines with strand-direction chevrons
    - Automatic fallback for features without CDS annotation

Uses both a synthetic BED12 DataFrame and real annotation files.

Run:  python examples/scripts/genome_tracks_bed12.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack,
    AnnotationTrack,
    GenomicInterval,
    plot_tracks,
    read_bed,
    read_gff,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
os.makedirs(FIG_DIR, exist_ok=True)

gtrack = GenomeAxisTrack()


# ---------------------------------------------------------------------------
# 1. Synthetic BED12 transcripts — two genes with UTR/CDS structure
# ---------------------------------------------------------------------------
bed12_df = pd.DataFrame([
    {
        "chrom": "chr1", "start": 10_000, "end": 30_000,
        "name": "GeneA (+)", "score": 0, "strand": "+",
        "thick_start": 12_000, "thick_end": 27_000,
        "block_count": 4,
        "block_sizes": "1500,2000,3000,1500",
        "block_starts": "0,5000,12000,18500",
    },
    {
        "chrom": "chr1", "start": 35_000, "end": 55_000,
        "name": "GeneB (-)", "score": 0, "strand": "-",
        "thick_start": 37_000, "thick_end": 52_000,
        "block_count": 3,
        "block_sizes": "2000,4000,2500",
        "block_starts": "0,7000,17500",
    },
    {
        "chrom": "chr1", "start": 60_000, "end": 72_000,
        "name": "GeneC (no CDS)", "score": 0, "strand": "+",
        "thick_start": np.nan, "thick_end": np.nan,
        "block_count": 2,
        "block_sizes": "3000,3000",
        "block_starts": "0,9000",
    },
])

region = GenomicInterval("chr1", 8_000, 75_000)

atrack = AnnotationTrack(bed12_df, name="BED12 Transcripts", stacking="squish")
axes = plot_tracks(
    [gtrack, atrack],
    region=region,
    figsize=(14, 4),
    title="BED12 Transcript Rendering",
)
fig = axes[0].figure
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_bed12.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: {os.path.join(FIG_DIR, 'genome_tracks_bed12.png')}")


# ---------------------------------------------------------------------------
# 2. Side-by-side: BED12 with feature labels
# ---------------------------------------------------------------------------
atrack_labeled = AnnotationTrack(
    bed12_df, name="BED12 + Labels",
    stacking="squish", show_label=True, label_pos="above",
)
axes = plot_tracks(
    [gtrack, atrack_labeled],
    region=region,
    figsize=(14, 5),
    title="BED12 with Feature Labels",
)
fig = axes[0].figure
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_bed12_labels.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: {os.path.join(FIG_DIR, 'genome_tracks_bed12_labels.png')}")


# ---------------------------------------------------------------------------
# 3. BED12 with feature_filter — show only CDS-containing transcripts
# ---------------------------------------------------------------------------
def has_cds(row):
    """Keep only transcripts that have a CDS region."""
    return pd.notna(row.get("thick_start")) and pd.notna(row.get("thick_end"))


atrack_filtered = AnnotationTrack(
    bed12_df,
    name="CDS Only (feature_filter)",
    stacking="squish",
    feature_filter=has_cds,
)
axes = plot_tracks(
    [gtrack, atrack_filtered],
    region=region,
    figsize=(14, 4),
    title="BED12 with feature_filter (CDS only)",
)
fig = axes[0].figure
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_bed12_filtered.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: {os.path.join(FIG_DIR, 'genome_tracks_bed12_filtered.png')}")


# ---------------------------------------------------------------------------
# 4. Comparison: regular BED vs BED12 rendering
# ---------------------------------------------------------------------------
# Regular BED (no thick columns) — all features drawn the same
bed_simple = bed12_df[["chrom", "start", "end", "name", "score", "strand"]].copy()
atrack_simple = AnnotationTrack(bed_simple, name="Simple BED (no CDS)")

atrack_full = AnnotationTrack(bed12_df, name="BED12 (CDS + UTR)")

axes = plot_tracks(
    [gtrack, atrack_simple, atrack_full],
    region=GenomicInterval("chr1", 8_000, 32_000),
    figsize=(14, 6),
    title="BED vs BED12 Comparison",
)
fig = axes[0].figure
fig.savefig(
    os.path.join(FIG_DIR, "genome_tracks_bed12_comparison.png"),
    dpi=150, bbox_inches="tight",
)
plt.close("all")
print(f"[INFO] Saved: {os.path.join(FIG_DIR, 'genome_tracks_bed12_comparison.png')}")

print("[INFO] BED12 examples complete.")
