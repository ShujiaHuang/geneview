#!/usr/bin/env python
"""
VCFTrack — VCF variant visualization with BAM alignments.

Mirrors the genomeview examples notebook's VCFTrack extension demo
(Cell 24), showing how to display SNPs from a VCF file alongside
read alignments.

Examples:
  1. Basic VCFTrack with default alt-allele coloring
  2. Custom color_fn (by variant quality)
  3. VCF + BAM pileup (genomeview fig1 equivalent)
  4. VCF with coverage and annotation tracks

Run:  python examples/scripts/genome_tracks_vcf.py
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

# Data paths
VCF_FILE = os.path.join(DATA_DIR, "hg002.chr14.vcf.gz")
ILLUMINA_BAM = os.path.join(DATA_DIR, "illumina.chr14.bam")
PACBIO_BAM = os.path.join(DATA_DIR, "pacbio.chr14.bam")
BED_FILE = os.path.join(DATA_DIR, "chr1_200mb.refseq.sorted.bed.gz")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    VCFTrack,
    AnnotationTrack,
    GenomicInterval,
    read_bed,
    save_figure,
)

# Check data availability
if not os.path.exists(VCF_FILE):
    print(f"[SKIP] Missing {VCF_FILE}")
    sys.exit(0)

# =========================================================================
# 1. Basic VCFTrack with default coloring (by alt allele)
# =========================================================================
print("[1] VCFTrack — default alt-allele coloring")

chrom = "14"
start = 66903600
end = start + 1500
region = GenomicInterval(chrom, start, end)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

axis = GenomeAxisTrack()
axis.draw(ax1, region)
ax1.set_title(f"Genome Axis — {chrom}:{start}-{end}", fontsize=10, loc="left")

vcf_track = VCFTrack(VCF_FILE, name="HG002 Variants")
vcf_track.draw(ax2, region)
ax2.set_title("VCFTrack — colored by alt allele", fontsize=10, loc="left")

# Add legend for nucleotide colors
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#0072B2", label="A"),
    Patch(facecolor="#E69F00", label="C"),
    Patch(facecolor="#009E73", label="G"),
    Patch(facecolor="#D55E00", label="T"),
]
ax2.legend(handles=legend_elements, loc="upper right", fontsize=7,
           framealpha=0.8, title="Alt Allele", ncol=4)

plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_vcf_basic.png"))
plt.close(fig)
print("  → genome_tracks_vcf_basic.png")

# =========================================================================
# 2. VCFTrack with custom color_fn (by variant quality)
# =========================================================================
print("[2] VCFTrack — color by variant quality")


def color_by_quality(variant):
    """Color variants by QUAL score."""
    qual = variant.qual
    if qual is None:
        return "#999999"
    if qual >= 50:
        return "#2E7D32"   # dark green = high quality
    elif qual >= 20:
        return "#FFC107"   # amber = medium quality
    else:
        return "#C62828"   # dark red = low quality


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

vcf_default = VCFTrack(VCF_FILE, name="Default (alt allele)")
vcf_default.draw(ax1, region)
ax1.set_title("VCFTrack — Default (alt allele)", fontsize=10, loc="left")

vcf_qual = VCFTrack(VCF_FILE, color_fn=color_by_quality, name="By Quality")
vcf_qual.draw(ax2, region)
ax2.set_title("VCFTrack — color_fn by QUAL score", fontsize=10, loc="left")

legend_qual = [
    Patch(facecolor="#2E7D32", label="QUAL ≥ 50"),
    Patch(facecolor="#FFC107", label="QUAL 20-49"),
    Patch(facecolor="#C62828", label="QUAL < 20"),
]
ax2.legend(handles=legend_qual, loc="upper right", fontsize=7,
           framealpha=0.8, title="Variant Quality")

plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_vcf_quality.png"))
plt.close(fig)
print("  → genome_tracks_vcf_quality.png")

# =========================================================================
# 3. VCF + BAM pileup (genomeview "fig1" equivalent)
# =========================================================================
print("[3] VCF + BAM pileup — genomeview fig1 equivalent")

tracks = []
tracks.append(GenomeAxisTrack())

if os.path.exists(VCF_FILE):
    tracks.append(VCFTrack(
        VCF_FILE, name="HG002 Variants",
        color_fn=lambda v: {
            "A": "#0072B2", "C": "#E69F00",
            "G": "#009E73", "T": "#D55E00",
        }.get(str(v.alts[0]).upper(), "#999999") if v.alts else "#999999",
    ))

if os.path.exists(ILLUMINA_BAM):
    tracks.append(AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Illumina",
        color_fn=lambda read: "lightgray",
    ))

if os.path.exists(PACBIO_BAM):
    tracks.append(AlignmentsTrack(
        filepath=PACBIO_BAM,
        type="pileup", name="PacBio",
        min_indel_size=10,
        color_fn=lambda read: "lightgray",
    ))

axes = plot_tracks(tracks, region=region, figsize=(14, 12))
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_vcf_with_bam.png"))
plt.close("all")
print("  → genome_tracks_vcf_with_bam.png")

# =========================================================================
# 4. VCF with wider region + coverage + BED annotations
# =========================================================================
print("[4] VCF with coverage track and gene annotations")

region_wide = GenomicInterval(chrom, 66901400, 66901400 + 5000)

tracks_wide = [GenomeAxisTrack()]

if os.path.exists(BED_FILE):
    try:
        bed_data = read_bed(BED_FILE)
        bed_sub = bed_data[
            (bed_data["chrom"].isin(["14", "chr14"])) &
            (bed_data["end"] > region_wide.start) &
            (bed_data["start"] < region_wide.end)
        ]
        if len(bed_sub) > 0:
            tracks_wide.append(AnnotationTrack(bed_sub, name="Genes"))
    except Exception:
        pass

tracks_wide.append(VCFTrack(VCF_FILE, name="HG002 SNPs"))

if os.path.exists(ILLUMINA_BAM):
    tracks_wide.append(
        BAMCoverageTrack(filepath=ILLUMINA_BAM, name="Illumina Coverage")
    )
    tracks_wide.append(AlignmentsTrack(
        filepath=ILLUMINA_BAM, is_paired=True,
        type="pileup", name="Illumina Reads",
        color_fn=lambda r: "lightgray",
    ))

axes = plot_tracks(tracks_wide, region=region_wide, figsize=(14, 14))
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_vcf_comprehensive.png"))
plt.close("all")
print("  → genome_tracks_vcf_comprehensive.png")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All VCF figures generated in:", FIG_DIR)
