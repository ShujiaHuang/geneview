"""Grouped Alignments & Quick Consensus Example
=================================================

Demonstrates new features ported from genomeview:

    - ``GroupedAlignmentsTrack`` — split BAM reads into groups
    - ``get_group_by_tag_fn`` — group reads by BAM tag value
    - ``quick_consensus`` — filter pileup mismatches by consensus threshold
      (essential for noisy long-read data)
    - ``read_filter`` — custom callable for filtering reads before display

Uses ``quick_consensus_test.bam`` (PacBio) and ``test.bam`` (Illumina PE).

Run:  python examples/scripts/genome_tracks_grouped.py
"""
import os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    GroupedAlignmentsTrack,
    GenomicInterval,
    get_group_by_tag_fn,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
os.makedirs(FIG_DIR, exist_ok=True)

try:
    import pysam
    has_pysam = True
except ImportError:
    has_pysam = False
    print("[WARN] pysam not available; skipping BAM examples.")

if has_pysam:
    gtrack = GenomeAxisTrack()

    # -----------------------------------------------------------------------
    # 1. Quick consensus filtering — PacBio long reads
    #    quick_consensus_test.bam contains ~282 PacBio reads on chr4:96.5M
    #    Without consensus filtering, pileup is overwhelmed by errors.
    # -----------------------------------------------------------------------
    pacbio_bam = os.path.join(DATA_DIR, "quick_consensus_test.bam")

    # No reference passed — test.fa only covers chr1/chr2, not chr4.
    # Mismatches are still shown in pileup as coloured segments (no base
    # identity without a reference, but structural variants visible).

    # Without quick consensus — all mismatches shown
    track_no_qc = AlignmentsTrack(
        filepath=pacbio_bam,
        type="pileup",
        name="No Consensus Filter",
        quick_consensus=False,
    )
    region_pb = GenomicInterval("4", 96_543_500, 96_544_000)

    axes = plot_tracks(
        [gtrack, track_no_qc],
        region=region_pb,
        figsize=(14, 8),
        title="PacBio Pileup — No Consensus Filter",
    )
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_grouped_no_consensus.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_grouped_no_consensus.png")

    # With quick consensus — only confident mismatches shown
    track_qc = AlignmentsTrack(
        filepath=pacbio_bam,
        type="pileup",
        name="Consensus Filter (0.2)",
        quick_consensus=True,
        consensus_threshold=0.2,
    )

    axes = plot_tracks(
        [gtrack, track_qc],
        region=region_pb,
        figsize=(14, 8),
        title="PacBio Pileup — Quick Consensus (threshold=0.2)",
    )
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_grouped_consensus.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_grouped_consensus.png")

    # -----------------------------------------------------------------------
    # 2. Grouped alignments — group test.bam reads by mapping quality
    #    Uses a custom keyfn (not a tag-based one, since test.bam has no
    #    HP/RG tags).
    # -----------------------------------------------------------------------
    test_bam = os.path.join(DATA_DIR, "test.bam")

    def mapq_group(read):
        """Group reads into high/low mapping quality bins."""
        if read.mapping_quality >= 50:
            return "High MAPQ (>=50)"
        else:
            return "Low MAPQ (<50)"

    region_test = GenomicInterval("chr1", 189_892_000, 189_895_000)

    grouped_track = GroupedAlignmentsTrack(
        filepath=test_bam,
        keyfn=mapq_group,
        type="coverage",
        name="Grouped by MAPQ",
        height=4.0,
    )

    axes = plot_tracks(
        [gtrack, grouped_track],
        region=region_test,
        figsize=(14, 6),
        title="Grouped Alignments — Reads by Mapping Quality",
    )
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_grouped_mapq.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_grouped_mapq.png")

    # -----------------------------------------------------------------------
    # 3. Custom read_filter — show only primary alignments
    # -----------------------------------------------------------------------
    def primary_only(read):
        """Keep only primary alignments (flag & 0x900 == 0)."""
        return not read.is_secondary and not read.is_supplementary

    track_primary = AlignmentsTrack(
        filepath=test_bam,
        type="coverage",
        name="Primary Only (read_filter)",
        read_filter=primary_only,
    )
    track_all = AlignmentsTrack(
        filepath=test_bam,
        type="coverage",
        name="All Reads",
    )

    axes = plot_tracks(
        [gtrack, track_all, track_primary],
        region=region_test,
        figsize=(14, 8),
        title="read_filter — Primary Alignments Only",
    )
    fig = axes[0].figure
    fig.savefig(
        os.path.join(FIG_DIR, "genome_tracks_grouped_read_filter.png"),
        dpi=150, bbox_inches="tight",
    )
    plt.close("all")
    print(f"[INFO] Saved: genome_tracks_grouped_read_filter.png")

    # -----------------------------------------------------------------------
    # 4. Group by BAM tag (HP) — 10x Genomics data (genomeview Cell 16)
    #    Mirrors genomeview's GroupedBAMTrack with get_group_by_tag_fn("HP")
    # -----------------------------------------------------------------------
    tenx_bam = os.path.join(DATA_DIR, "10x.chr14.bam")

    if os.path.exists(tenx_bam):
        region_10x = GenomicInterval("14", 66901400, 66901400 + 450)

        grouped_hp = GroupedAlignmentsTrack(
            filepath=tenx_bam,
            keyfn=get_group_by_tag_fn("HP"),
            is_paired=True,
            type="pileup",
            name="10x Grouped by HP",
            category_label_fn=lambda x: "haplotype_{}".format(x),
            height=5.0,
        )

        axes = plot_tracks(
            [gtrack, grouped_hp],
            region=region_10x,
            figsize=(14, 10),
            title="10x Genomics — Grouped by HP Tag (genomeview Cell 16)",
        )
        fig = axes[0].figure
        fig.savefig(
            os.path.join(FIG_DIR, "genome_tracks_grouped_hp_tag.png"),
            dpi=150, bbox_inches="tight",
        )
        plt.close("all")
        print(f"[INFO] Saved: genome_tracks_grouped_hp_tag.png")
    else:
        print(f"[SKIP] 10x BAM not found: {tenx_bam}")

print("[INFO] Grouped alignments examples complete.")
