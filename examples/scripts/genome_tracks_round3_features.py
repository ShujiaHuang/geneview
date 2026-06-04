#!/usr/bin/env python
"""
Demonstrate genomeview features ported to genometracks (round 3).

New features covered:
  1. BAMCoverageTrack — standalone per-base coverage line/fill
  2. overlap_color — highlight overlapping paired-end read mates
  3. draw_read_labels — show read query names on pileup
  4. min_insertion_label_size — auto-label large insertions
  5. get_ticks() — standalone axis tick utility
  6. find_tracks() — retrieve tracks by name or type
  7. read_bigbed — BigBed I/O support (stub — no test data)
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "genome_tracks")
FIG_DIR = os.path.join(SCRIPT_DIR, "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Test BAM files
BAM_CHR1 = os.path.join(DATA_DIR, "test.bam")
BAM_INDELS = os.path.join(DATA_DIR, "indels.bam")

# Reference FASTA
REF_FA = os.path.join(DATA_DIR, "test.fa")

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AnnotationTrack,
    AlignmentsTrack,
    BAMCoverageTrack,
    GenomicInterval,
    get_ticks,
    find_tracks,
    save_figure,
)

# =========================================================================
# 1. BAMCoverageTrack — line and fill modes
# =========================================================================
print("[1] BAMCoverageTrack — line vs fill modes")

region_chr1 = GenomicInterval("chr1", 189_860_000, 189_870_000)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

cov_line = BAMCoverageTrack(filepath=BAM_CHR1, type="line", col="#2196F3",
                             name="Coverage (line)")
cov_fill = BAMCoverageTrack(filepath=BAM_CHR1, type="fill", col="#4CAF50",
                             name="Coverage (fill)")

cov_line.draw(ax1, region_chr1)
cov_fill.draw(ax2, region_chr1)

ax1.set_title("BAMCoverageTrack — Line Mode", fontsize=10, loc="left")
ax2.set_title("BAMCoverageTrack — Fill Mode", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_bam_coverage.png"))
plt.close(fig)
print("  → genome_tracks_bam_coverage.png")

# =========================================================================
# 2. BAMCoverageTrack with transformation (log1p)
# =========================================================================
print("[2] BAMCoverageTrack — log-transformed coverage")

fig, ax = plt.subplots(figsize=(14, 3))
cov_log = BAMCoverageTrack(filepath=BAM_CHR1, type="fill", col="#FF9800",
                            transformation=np.log1p,
                            name="log1p Coverage")
cov_log.draw(ax, region_chr1)
ax.set_title("BAMCoverageTrack — log1p Transformed", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax], os.path.join(FIG_DIR, "genome_tracks_bam_coverage_log.png"))
plt.close(fig)
print("  → genome_tracks_bam_coverage_log.png")

# =========================================================================
# 3. BAMCoverageTrack in multi-track plot
# =========================================================================
print("[3] BAMCoverageTrack in multi-track plot with GenomeAxisTrack")

axis = GenomeAxisTrack()
cov_track = BAMCoverageTrack(filepath=BAM_CHR1, type="fill", col="#5B8DB8",
                              name="BAM Coverage")

axes = plot_tracks([axis, cov_track], region=region_chr1, figsize=(14, 4))
save_figure(axes, os.path.join(FIG_DIR, "genome_tracks_bam_coverage_multi.png"))
plt.close("all")
print("  → genome_tracks_bam_coverage_multi.png")

# =========================================================================
# 4. draw_read_labels — show read names on pileup
# =========================================================================
print("[4] draw_read_labels — read names on pileup")

region_zoom = GenomicInterval("chr1", 189_862_000, 189_863_000)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

track_no_labels = AlignmentsTrack(
    filepath=BAM_CHR1, type="pileup", name="No Labels",
    reference=REF_FA,
)
track_labels = AlignmentsTrack(
    filepath=BAM_CHR1, type="pileup", name="With Read Labels",
    reference=REF_FA, draw_read_labels=True,
)

track_no_labels.draw(ax1, region_zoom)
track_labels.draw(ax2, region_zoom)

ax1.set_title("Pileup — Default (no labels)", fontsize=10, loc="left")
ax2.set_title("Pileup — draw_read_labels=True", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_read_labels.png"))
plt.close(fig)
print("  → genome_tracks_read_labels.png")

# =========================================================================
# 5. overlap_color — paired-end overlap highlighting
# =========================================================================
print("[5] overlap_color — paired-end overlap highlighting")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

track_default = AlignmentsTrack(
    filepath=BAM_CHR1, type="pileup", name="Default",
    is_paired=True, reference=REF_FA,
)
track_overlap = AlignmentsTrack(
    filepath=BAM_CHR1, type="pileup", name="Overlap Highlight",
    is_paired=True, reference=REF_FA, overlap_color="lime",
)

track_default.draw(ax1, region_zoom)
track_overlap.draw(ax2, region_zoom)

ax1.set_title("Paired-end — Default", fontsize=10, loc="left")
ax2.set_title("Paired-end — overlap_color='lime'", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_overlap_color.png"))
plt.close(fig)
print("  → genome_tracks_overlap_color.png")

# =========================================================================
# 6. min_insertion_label_size — auto-labelling large insertions
# =========================================================================
print("[6] min_insertion_label_size — auto-labelling large insertions")

region_indels = GenomicInterval("chr2", 126_455_000, 126_460_000)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

track_small_thresh = AlignmentsTrack(
    filepath=BAM_INDELS, type="pileup", name="min_insertion_label_size=5 (default)",
    reference=REF_FA, show_indels=True, min_insertion_label_size=5,
)
track_large_thresh = AlignmentsTrack(
    filepath=BAM_INDELS, type="pileup", name="min_insertion_label_size=20",
    reference=REF_FA, show_indels=True, min_insertion_label_size=20,
)

track_small_thresh.draw(ax1, region_indels)
track_large_thresh.draw(ax2, region_indels)

ax1.set_title("min_insertion_label_size=5 (default)", fontsize=10, loc="left")
ax2.set_title("min_insertion_label_size=20", fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax1, ax2], os.path.join(FIG_DIR, "genome_tracks_min_insertion_label.png"))
plt.close(fig)
print("  → genome_tracks_min_insertion_label.png")

# =========================================================================
# 7. Combined features — all new options together
# =========================================================================
print("[7] Combined features — labels + overlap + auto-insertion labels")

track_combined = AlignmentsTrack(
    filepath=BAM_CHR1,
    type="pileup",
    is_paired=True,
    reference=REF_FA,
    draw_read_labels=True,
    overlap_color="lime",
    min_insertion_label_size=3,
    color_by_strand=True,
    show_clipping=True,
    name="All Features Combined",
)

fig, ax = plt.subplots(figsize=(14, 5))
track_combined.draw(ax, region_zoom)
ax.set_title("Combined: read_labels + overlap_color + strand + clipping",
             fontsize=10, loc="left")
plt.tight_layout()
save_figure([ax], os.path.join(FIG_DIR, "genome_tracks_round3_combined.png"))
plt.close(fig)
print("  → genome_tracks_round3_combined.png")

# =========================================================================
# 8. get_ticks() — standalone tick computation demo
# =========================================================================
print("[8] get_ticks() — standalone tick utility")

fig, axes = plt.subplots(3, 1, figsize=(12, 6))

ranges = [
    ("chr1:0-500 (bp)", 0, 500),
    ("chr7:26M-27M (kb)", 26_000_000, 27_000_000),
    ("chr1:1-250M (Mb)", 1, 250_000_000),
]

for ax, (title, start, end) in zip(axes, ranges):
    ticks = get_ticks(start, end, target_n_labels=8)
    positions = [pos for pos, _ in ticks]
    labels = [lbl for _, lbl in ticks]

    ax.plot([start, end], [0.5, 0.5], "k-", linewidth=2)
    for pos, lbl in ticks:
        ax.plot([pos, pos], [0.3, 0.7], "k-", linewidth=1)
        ax.text(pos, 0.1, lbl, ha="center", va="top", fontsize=7)
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=9, loc="left")
    ax.axis("off")

plt.suptitle("get_ticks() — Auto-Scaling Genomic Tick Labels", fontsize=12)
plt.tight_layout()
save_figure(list(axes), os.path.join(FIG_DIR, "genome_tracks_get_ticks.png"))
plt.close(fig)
print("  → genome_tracks_get_ticks.png")

# =========================================================================
# 9. find_tracks() — track retrieval demo
# =========================================================================
print("[9] find_tracks() — track retrieval by name and type")

import pandas as pd

ann1 = AnnotationTrack(
    pd.DataFrame({"chrom": ["chr1"] * 3, "start": [100, 500, 900],
                   "end": [200, 700, 1000], "name": ["A", "B", "C"]}),
    name="Features",
)
ann2 = AnnotationTrack(
    pd.DataFrame({"chrom": ["chr1"], "start": [300], "end": [400],
                   "name": ["D"]}),
    name="Genes",
)
axis = GenomeAxisTrack(name="Axis")

track_list = [axis, ann1, ann2]

# Find by name
by_name = find_tracks(track_list, name="Genes")
print(f"  find_tracks(name='Genes') → {by_name}")

# Find by type
by_type = find_tracks(track_list, track_type=GenomeAxisTrack)
print(f"  find_tracks(track_type=GenomeAxisTrack) → {by_type}")

# Find all
all_tracks = find_tracks(track_list)
print(f"  find_tracks() → {len(all_tracks)} tracks")

# =========================================================================
# Summary
# =========================================================================
print("\n✓ All round-3 figures generated in:", FIG_DIR)
