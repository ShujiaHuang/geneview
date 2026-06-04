#!/usr/bin/env python
"""
Basic BAM visualization — mirrors genomeview's Cells 4/6.

Shows PacBio and Illumina read alignments side by side using default
pileup rendering, equivalent to genomeview's basic ``visualize_data()``
example.

Examples:
  1. PacBio pileup (long reads, chr14:66901400)
  2. Illumina paired-end pileup (short reads)
  3. Combined: axis + pacbio + illumina in one figure
  4. Illumina coverage + pileup

Run:  python examples/scripts/genome_tracks_bam_basic.py
"""

import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "genome_tracks")
FIG_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Data paths (genomeview uses data/pacbio.chr14.bam, data/illumina.chr14.bam)
PACBIO_BAM = os.path.join(DATA_DIR, "pacbio.chr14.bam")
ILLUMINA_BAM = os.path.join(DATA_DIR, "illumina.chr14.bam")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    GenomicInterval,
    save_figure,
)

# Check data availability
if not os.path.exists(PACBIO_BAM) or not os.path.exists(ILLUMINA_BAM):
    missing = []
    if not os.path.exists(PACBIO_BAM):
        missing.append(PACBIO_BAM)
    if not os.path.exists(ILLUMINA_BAM):
        missing.append(ILLUMINA_BAM)
    print(f"[SKIP] Missing data files: {missing}")
    sys.exit(0)

# =========================================================================
# Genomeview parameters: chrom = "14", start = 66901400
# =========================================================================
chrom = "14"
start = 66901400
region = GenomicInterval(chrom, start, start + 10000)

# =========================================================================
# 1. PacBio pileup — long reads (genomeview Cell 4/6, pacbio panel)
# =========================================================================
print("[1] Basic PacBio pileup — long reads on chr14")

axes = plot_tracks(
    [
        GenomeAxisTrack(),
        AlignmentsTrack(
            filepath=PACBIO_BAM,
            type="pileup",
            name="pacbio",
        ),
    ],
    region=region,
    figsize=(14, 8),
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_bam_basic_pacbio.png"))
plt.close("all")
print("  → genome_tracks_bam_basic_pacbio.png")

# =========================================================================
# 2. Illumina paired-end pileup (genomeview Cell 4/6, illumina panel)
# =========================================================================
print("[2] Basic Illumina paired-end pileup — short reads on chr14")

axes = plot_tracks(
    [
        GenomeAxisTrack(),
        AlignmentsTrack(
            filepath=ILLUMINA_BAM,
            is_paired=True,
            type="pileup",
            name="illumina",
        ),
    ],
    region=region,
    figsize=(14, 8),
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_bam_basic_illumina.png"))
plt.close("all")
print("  → genome_tracks_bam_basic_illumina.png")

# =========================================================================
# 3. Combined: PacBio + Illumina in one figure (genomeview Cell 4/6)
#    genomeview.visualize_data(track_info, chrom, start, start+10000, ...)
# =========================================================================
print("[3] Combined PacBio + Illumina pileup — genomeview Cell 4/6 equivalent")

axes = plot_tracks(
    [
        GenomeAxisTrack(),
        AlignmentsTrack(
            filepath=PACBIO_BAM,
            type="pileup",
            name="pacbio",
        ),
        AlignmentsTrack(
            filepath=ILLUMINA_BAM,
            is_paired=True,
            type="pileup",
            name="illumina",
        ),
    ],
    region=region,
    figsize=(14, 12),
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_bam_basic_combined.png"))
plt.close("all")
print("  → genome_tracks_bam_basic_combined.png")

# =========================================================================
# 4. Illumina coverage + pileup (enhanced version with coverage overlay)
# =========================================================================
print("[4] Illumina coverage histogram + pileup")

axes = plot_tracks(
    [
        GenomeAxisTrack(),
        BAMCoverageTrack(filepath=ILLUMINA_BAM, name="Illumina Coverage"),
        AlignmentsTrack(
            filepath=ILLUMINA_BAM,
            is_paired=True,
            type="pileup",
            name="illumina",
        ),
    ],
    region=region,
    figsize=(14, 10),
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_bam_basic_coverage.png"))
plt.close("all")
print("  → genome_tracks_bam_basic_coverage.png")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All basic BAM figures generated in:", FIG_DIR)
