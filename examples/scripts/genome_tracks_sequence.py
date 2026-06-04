"""
SequenceTrack Example
=====================

Demonstrates the SequenceTrack for displaying nucleotide sequences,
including zoom-level rendering, complement mode, and custom colors.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, SequenceTrack, GenomicInterval,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")

# ---------------------------------------------------------------------------
# 1. Basic sequence with colored letters (small window)
# ---------------------------------------------------------------------------
seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG" * 3  # ~108 bp
track = SequenceTrack(sequence=seq, name="Sequence (letters)")
region = GenomicInterval("chr1", 0, len(seq))

axes = plot_tracks(
    [GenomeAxisTrack(), track],
    region=region,
    figsize=(12, 3),
    title="SequenceTrack - Letters",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_sequence_letters.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Medium zoom - colored boxes
# ---------------------------------------------------------------------------
seq_medium = "ATCGATCG" * 100  # ~800 bp
track_medium = SequenceTrack(sequence=seq_medium, name="Sequence (boxes)")
region_medium = GenomicInterval("chr1", 0, len(seq_medium))

axes = plot_tracks(
    [GenomeAxisTrack(), track_medium],
    region=region_medium,
    figsize=(12, 3),
    title="SequenceTrack - Boxes",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_sequence_boxes.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Large zoom - colored line
# ---------------------------------------------------------------------------
seq_large = "ATCGATCG" * 500  # ~4000 bp
track_large = SequenceTrack(sequence=seq_large, name="Sequence (line)")
region_large = GenomicInterval("chr1", 0, len(seq_large))

axes = plot_tracks(
    [GenomeAxisTrack(), track_large],
    region=region_large,
    figsize=(12, 3),
    title="SequenceTrack - Line",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_sequence_line.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Complement and 5'->3' arrow
# ---------------------------------------------------------------------------
seq_small = "AAAAATTTTT" * 5
track_comp = SequenceTrack(
    sequence=seq_small, complement=True, add53=True,
    name="Complement + 5'->3'",
)
region_small = GenomicInterval("chr1", 0, len(seq_small))

axes = plot_tracks(
    [GenomeAxisTrack(), track_comp],
    region=region_small,
    figsize=(12, 3),
    title="SequenceTrack - Complement & Direction",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_sequence_complement.png"), dpi=150, bbox_inches="tight")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. Custom nucleotide colors
# ---------------------------------------------------------------------------
custom_colors = {
    "A": "#FF0000",  # red
    "C": "#00FF00",  # green
    "G": "#0000FF",  # blue
    "T": "#FFFF00",  # yellow
}
track_custom = SequenceTrack(
    sequence=seq, fontcolor=custom_colors, name="Custom Colors",
)

axes = plot_tracks(
    [GenomeAxisTrack(), track_custom],
    region=region,
    figsize=(12, 3),
    title="SequenceTrack - Custom Colors",
)
plt.savefig(os.path.join(FIG_DIR, "genome_tracks_sequence_custom.png"), dpi=150, bbox_inches="tight")
plt.close("all")

print("All SequenceTrack examples saved.")
