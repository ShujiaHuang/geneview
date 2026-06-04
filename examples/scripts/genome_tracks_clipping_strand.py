"""
Clipping, Strand Coloring, and Indel Filtering Examples
=======================================================

Demonstrates new AlignmentsTrack features ported from genomeview:
  - show_clipping / col_clipping
  - color_by_strand / fill_reads_fwd / fill_reads_rev
  - min_indel_size
  - show_insertion_labels
  - include_secondary
  - save_figure()
  - reverse_comp()
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, AlignmentsTrack, GenomicInterval,
    save_figure, reverse_comp,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")

try:
    import pysam
    has_pysam = True
except ImportError:
    has_pysam = False

if not has_pysam:
    print("Skipped: pysam not installed.")
    exit()

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------
test_bam = os.path.join(DATA_DIR, "test.bam")
test_fa = os.path.join(DATA_DIR, "test.fa")
indels_bam = os.path.join(DATA_DIR, "indels.bam")

region_chr1 = GenomicInterval("chr1", 189_892_000, 189_896_000)
region_indels = GenomicInterval("chr2", 126_389_500, 126_392_000)

# ---------------------------------------------------------------------------
# 1. Default pileup vs. color_by_strand
# ---------------------------------------------------------------------------
print("1. color_by_strand example...")
track_default = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    name="Default Coloring",
)
track_strand = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    color_by_strand=True,
    name="Color by Strand",
)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
track_default.draw(ax1, region_chr1)
track_strand.draw(ax2, region_chr1)
ax1.set_title("Default Coloring", fontsize=10, loc="left")
ax2.set_title("Color by Strand (fwd=red, rev=blue)", fontsize=10, loc="left")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_strand_coloring.png"),
            dpi=150, bbox_inches="tight")
plt.close("all")
print("  → genome_tracks_strand_coloring.png")

# ---------------------------------------------------------------------------
# 2. Custom strand colors
# ---------------------------------------------------------------------------
print("2. custom strand colors...")
track_custom = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    color_by_strand=True,
    fill_reads_fwd="#FF6B6B",
    fill_reads_rev="#4ECDC4",
    name="Custom Strand Colors",
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_custom],
    region=region_chr1,
    figsize=(14, 5),
    title="Custom Strand Colors (fwd=#FF6B6B, rev=#4ECDC4)",
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_strand_custom_colors.png"))
plt.close("all")
print("  → genome_tracks_strand_custom_colors.png")

# ---------------------------------------------------------------------------
# 3. Show clipping visualization
# ---------------------------------------------------------------------------
print("3. show_clipping example...")
track_clipping = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    show_clipping=True,
    col_clipping="magenta",
    name="Clipping (magenta)",
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_clipping],
    region=region_chr1,
    figsize=(14, 5),
    title="Soft Clipping Visualization",
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_clipping.png"))
plt.close("all")
print("  → genome_tracks_clipping.png")

# ---------------------------------------------------------------------------
# 4. min_indel_size + show_insertion_labels
# ---------------------------------------------------------------------------
print("4. min_indel_size and insertion labels...")
track_all_indels = AlignmentsTrack(
    filepath=indels_bam, type="pileup",
    show_indels=True,
    name="All Indels",
)
track_filtered = AlignmentsTrack(
    filepath=indels_bam, type="pileup",
    show_indels=True,
    min_indel_size=5,
    show_insertion_labels=True,
    name="Filtered (≥5bp)",
)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
track_all_indels.draw(ax1, region_indels)
track_filtered.draw(ax2, region_indels)
ax1.set_title("All Indels", fontsize=10, loc="left")
ax2.set_title("Indels ≥5bp with Insertion Labels", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_indel_filtering.png"))
plt.close("all")
print("  → genome_tracks_indel_filtering.png")

# ---------------------------------------------------------------------------
# 5. include_secondary filtering
# ---------------------------------------------------------------------------
print("5. include_secondary filtering...")
track_with_sec = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    include_secondary=True,
    name="With Secondary",
)
track_no_sec = AlignmentsTrack(
    filepath=test_bam, type="pileup",
    include_secondary=False,
    name="Primary Only",
)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
track_with_sec.draw(ax1, region_chr1)
track_no_sec.draw(ax2, region_chr1)
ax1.set_title("Include Secondary/Supplementary", fontsize=10, loc="left")
ax2.set_title("Primary Alignments Only", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_secondary_filter.png"))
plt.close("all")
print("  → genome_tracks_secondary_filter.png")

# ---------------------------------------------------------------------------
# 6. Combined: clipping + strand coloring + filtered indels
# ---------------------------------------------------------------------------
print("6. combined features...")
track_combined = AlignmentsTrack(
    filepath=indels_bam, type="pileup",
    show_clipping=True, col_clipping="cyan",
    color_by_strand=True,
    show_indels=True, min_indel_size=3,
    show_insertion_labels=True,
    include_secondary=False,
    name="Combined Features",
)
axes = plot_tracks(
    [GenomeAxisTrack(), track_combined],
    region=region_indels,
    figsize=(14, 5),
    title="Combined: Clipping + Strand Color + Indel Filter + Labels",
)
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_combined_features.png"))
plt.close("all")
print("  → genome_tracks_combined_features.png")

# ---------------------------------------------------------------------------
# 7. reverse_comp() utility
# ---------------------------------------------------------------------------
print("7. reverse_comp() utility demo...")
seq = "ATCGATCG"
rc = reverse_comp(seq)
print(f"  {seq} → reverse complement → {rc}")

# Demonstrate in a SequenceTrack-like context
from geneview.genometracks import SequenceTrack
region_seq = GenomicInterval("chr1", 0, 50)
# Create a short sequence to show complement
seq_data = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
track_fwd = SequenceTrack(sequence=seq_data, chromosome="chr1", name="Forward")
track_rev = SequenceTrack(sequence=seq_data, chromosome="chr1",
                          complement=True, name="Reverse Comp")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 3), sharex=True)
track_fwd.draw(ax1, region_seq)
track_rev.draw(ax2, region_seq)
ax1.set_title("Forward 5'→3'", fontsize=10, loc="left")
ax2.set_title(f"Reverse Complement ({seq_data[:16]}... → {reverse_comp(seq_data[:16])}...)",
              fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_reverse_comp_demo.png"))
plt.close("all")
print("  → genome_tracks_reverse_comp_demo.png")

print("\nAll clipping/strand/export examples saved.")
