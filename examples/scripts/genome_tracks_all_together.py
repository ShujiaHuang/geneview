#!/usr/bin/env python
"""
"All Together" combined figure — mirrors genomeview's Cell 18/19 (fig1.svg).

Shows PacBio + Illumina + BED annotations in a single multi-panel figure,
with custom indel display settings (min_indel_size, min_insertion_label_size).
This is the "showcase" figure from the genomeview examples notebook.

Examples:
  1. PacBio + Illumina + BED on chr1 (genomeview Cell 18 equivalent)
  2. Same with VCF variants overlaid
  3. VCF + PacBio + Illumina + BED on chr14 (genomeview Cell 24 fig1)

Run:  python examples/scripts/genome_tracks_all_together.py
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

# Data paths (genomeview Cell 18 uses chr1 data)
PACBIO_CHR1 = os.path.join(DATA_DIR, "pacbio.chr1.bam")
ILLUMINA_CHR1 = os.path.join(DATA_DIR, "illumina.chr1.bam")
BED_CHR1 = os.path.join(DATA_DIR, "chr1_200mb.refseq.sorted.bed.gz")

# chr14 data (for VCF combined example)
PACBIO_CHR14 = os.path.join(DATA_DIR, "pacbio.chr14.bam")
ILLUMINA_CHR14 = os.path.join(DATA_DIR, "illumina.chr14.bam")
VCF_CHR14 = os.path.join(DATA_DIR, "hg002.chr14.vcf.gz")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    AnnotationTrack,
    VCFTrack,
    GenomicInterval,
    read_bed,
    save_figure,
)

# =========================================================================
# 1. "All Together" — PacBio + Illumina + BED on chr1 (genomeview Cell 18)
# =========================================================================
print("[1] All Together — PacBio + Illumina + BED on chr1 (genomeview fig1)")

chrom = "chr1"
start = 224375300
end = start + 20000
region_chr1 = GenomicInterval(chrom, start, end)

tracks = [GenomeAxisTrack()]

# BED annotations
if os.path.exists(BED_CHR1):
    try:
        bed_data = read_bed(BED_CHR1)
        bed_sub = bed_data[
            (bed_data["chrom"].isin(["chr1", "1"])) &
            (bed_data["end"] > start) &
            (bed_data["start"] < end)
        ]
        if len(bed_sub) > 0:
            tracks.append(AnnotationTrack(bed_sub, name="bed"))
    except Exception:
        pass

# PacBio pileup with custom indel display
if os.path.exists(PACBIO_CHR1):
    tracks.append(AlignmentsTrack(
        filepath=PACBIO_CHR1,
        type="pileup",
        name="pacbio",
        min_indel_size=50,           # genomeview: doc.get_tracks("pacbio")[0].min_indel_size = 50
        min_insertion_label_size=100,  # genomeview: min_insertion_label_size = 100
    ))

# Illumina paired-end pileup
if os.path.exists(ILLUMINA_CHR1):
    tracks.append(AlignmentsTrack(
        filepath=ILLUMINA_CHR1,
        is_paired=True,
        type="pileup",
        name="illumina",
    ))

if len(tracks) > 1:
    axes = plot_tracks(tracks, region=region_chr1, figsize=(14, 14))
    save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_all_together.png"))
    plt.close("all")
    print("  → genome_tracks_all_together.png")
else:
    print("  [SKIP] No data files found for chr1")

# =========================================================================
# 2. VCF + PacBio + Illumina on chr14 (genomeview Cell 24 fig1)
#    Mirrors the VCFTrack example: variants + reads in one view
# =========================================================================
print("[2] VCF + PacBio + Illumina on chr14 (genomeview Cell 24 fig1)")

chrom14 = "14"
start14 = 66903600
end14 = start14 + 1500
region_chr14 = GenomicInterval(chrom14, start14, end14)

tracks_vcf = [GenomeAxisTrack()]

if os.path.exists(VCF_CHR14):
    tracks_vcf.append(VCFTrack(
        VCF_CHR14,
        name="HG002 small variants",
    ))

if os.path.exists(ILLUMINA_CHR14):
    tracks_vcf.append(AlignmentsTrack(
        filepath=ILLUMINA_CHR14,
        type="pileup",
        name="HG002 illumina",
        color_fn=lambda x: "lightgray",
    ))

if os.path.exists(PACBIO_CHR14):
    tracks_vcf.append(AlignmentsTrack(
        filepath=PACBIO_CHR14,
        type="pileup",
        name="HG002 pacbio",
        min_indel_size=10,
        color_fn=lambda x: "lightgray",
    ))

if len(tracks_vcf) > 1:
    axes = plot_tracks(tracks_vcf, region=region_chr14, figsize=(14, 12))
    save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_all_together_vcf.png"))
    plt.close("all")
    print("  → genome_tracks_all_together_vcf.png")
else:
    print("  [SKIP] Missing chr14 data files")

# =========================================================================
# 3. Full showcase: VCF + Coverage + BED + PacBio + Illumina on chr14
# =========================================================================
print("[3] Full showcase — VCF + Coverage + BED + Reads on chr14")

region_wide = GenomicInterval("14", 66901400, 66901400 + 5000)

tracks_full = [GenomeAxisTrack()]

if os.path.exists(BED_CHR1):
    try:
        bed_wide = bed_data[
            (bed_data["chrom"].isin(["14", "chr14"])) &
            (bed_data["end"] > region_wide.start) &
            (bed_data["start"] < region_wide.end)
        ]
        if len(bed_wide) > 0:
            tracks_full.append(AnnotationTrack(bed_wide, name="Genes"))
    except Exception:
        pass

if os.path.exists(VCF_CHR14):
    tracks_full.append(VCFTrack(VCF_CHR14, name="HG002 Variants"))

if os.path.exists(ILLUMINA_CHR14):
    tracks_full.append(BAMCoverageTrack(
        filepath=ILLUMINA_CHR14, name="Illumina Coverage"
    ))
    tracks_full.append(AlignmentsTrack(
        filepath=ILLUMINA_CHR14,
        is_paired=True,
        type="pileup",
        name="Illumina Reads",
        color_fn=lambda r: "lightgray",
    ))

if os.path.exists(PACBIO_CHR14):
    tracks_full.append(AlignmentsTrack(
        filepath=PACBIO_CHR14,
        type="pileup",
        name="PacBio Reads",
        min_indel_size=10,
        color_fn=lambda r: "lightgray",
    ))

if len(tracks_full) > 2:
    axes = plot_tracks(tracks_full, region=region_wide, figsize=(14, 18))
    save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_all_together_full.png"))
    plt.close("all")
    print("  → genome_tracks_all_together_full.png")
else:
    print("  [SKIP] Insufficient data for full showcase")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All 'all together' figures generated in:", FIG_DIR)
