"""
AlignmentsTrack Example
=======================

Demonstrates the AlignmentsTrack for displaying BAM/CRAM read alignments
using real test data: coverage, pileup, sashimi, and combined views.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, AlignmentsTrack, GenomicInterval,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

# ---------------------------------------------------------------------------
# 1. Coverage view — gapped.bam (chr12:2966851-2966934, 50M reads)
# ---------------------------------------------------------------------------
try:
    import pysam
    has_pysam = True
except ImportError:
    has_pysam = False

if has_pysam:
    # Coverage histogram from gapped reads
    gapped_bam = os.path.join(DATA_DIR, "gapped.bam")
    region_cov = GenomicInterval("chr12", 2966800, 2966950)

    track_cov = AlignmentsTrack(
        filepath=gapped_bam,
        type="coverage",
        name="Coverage",
        fill_coverage="#5B8DB8",
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_cov],
        region=region_cov,
        figsize=(12, 4),
        title="AlignmentsTrack - Coverage",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_alignments_coverage.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # -----------------------------------------------------------------------
    # 2. Pileup view — individual reads
    # -----------------------------------------------------------------------
    track_pileup = AlignmentsTrack(
        filepath=gapped_bam,
        type="pileup",
        name="Pileup",
        fill_reads="#C8C8C8",
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_pileup],
        region=region_cov,
        figsize=(12, 6),
        title="AlignmentsTrack - Pileup",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_alignments_pileup.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # -----------------------------------------------------------------------
    # 3. Combined coverage + pileup
    # -----------------------------------------------------------------------
    track_combined = AlignmentsTrack(
        filepath=gapped_bam,
        type=["coverage", "pileup"],
        name="Coverage + Pileup",
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_combined],
        region=region_cov,
        figsize=(12, 6),
        title="AlignmentsTrack - Combined",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_alignments_combined.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # -----------------------------------------------------------------------
    # 4. Splice junction reads — test.bam (chr1, has N cigar op)
    # -----------------------------------------------------------------------
    test_bam = os.path.join(DATA_DIR, "test.bam")
    region_sashimi = GenomicInterval("chr1", 189891400, 189897800)

    track_sashimi = AlignmentsTrack(
        filepath=test_bam,
        type=["coverage", "sashimi"],
        name="Sashimi",
        sashimi_score=1,
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_sashimi],
        region=region_sashimi,
        figsize=(14, 6),
        title="AlignmentsTrack - Sashimi (Splice Junctions)",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_alignments_sashimi.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    # -----------------------------------------------------------------------
    # 5. Indel display — indels.bam (chr2, has I/D operations)
    # -----------------------------------------------------------------------
    indels_bam = os.path.join(DATA_DIR, "indels.bam")
    region_indels = GenomicInterval("chr2", 126389500, 126392000)

    track_indels = AlignmentsTrack(
        filepath=indels_bam,
        type="pileup",
        show_indels=True,
        name="Indels",
    )
    axes = plot_tracks(
        [GenomeAxisTrack(), track_indels],
        region=region_indels,
        figsize=(14, 6),
        title="AlignmentsTrack - Indels",
    )
    plt.savefig(os.path.join(FIG_DIR, "genome_tracks_alignments_indels.png"),
                dpi=150, bbox_inches="tight")
    plt.close("all")

    print("All AlignmentsTrack examples saved.")
else:
    print("AlignmentsTrack examples skipped: pysam not installed.")
