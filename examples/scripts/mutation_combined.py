"""Mutation tracks combined with other genome tracks.

Demonstrates LolliplotTrack and DandelionTrack used inside plot_tracks()
alongside GenomeAxisTrack and AnnotationTrack — the key feature of moving
mutation plots into the genometracks system.

Run:  python examples/scripts/mutation_combined.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from geneview.genometracks import (
    plot_tracks,
    GenomeAxisTrack,
    AnnotationTrack,
    LolliplotTrack,
    DandelionTrack,
    GenomicInterval,
)

# ---------------------------------------------------------------------------
# Locate example data
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "mutation")
SNP_FILE = os.path.join(DATA_DIR, "snps_example.tsv")
FEAT_FILE = os.path.join(DATA_DIR, "features_example.tsv")

snps = pd.read_csv(SNP_FILE, sep="\t")
feats = pd.read_csv(FEAT_FILE, sep="\t")

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

region = GenomicInterval("chr1", 0, 1500)

# ---------------------------------------------------------------------------
# 1. LolliplotTrack + GenomeAxisTrack
# ---------------------------------------------------------------------------
tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(snps, features=feats, type="circle", name="Lollipop"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 5),
                   title="LolliplotTrack + GenomeAxis")
plt.savefig(os.path.join(OUT_DIR, "mutation_combined_lollipop_axis.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_combined_lollipop_axis.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. AnnotationTrack (domains) + LolliplotTrack
# ---------------------------------------------------------------------------
tracks = [
    AnnotationTrack(feats, name="Domains"),
    LolliplotTrack(snps, features=feats, type="circle", name="Lollipop"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 6),
                   title="Annotation + LolliplotTrack")
plt.savefig(os.path.join(OUT_DIR, "mutation_combined_annotation_lollipop.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_combined_annotation_lollipop.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Full stack: GenomeAxis + Annotation + LolliplotTrack + DandelionTrack
# ---------------------------------------------------------------------------
tracks = [
    GenomeAxisTrack(),
    AnnotationTrack(feats, name="Domains"),
    LolliplotTrack(snps, features=feats, type="circle", name="Lollipop"),
    DandelionTrack(snps, features=feats, type="fan", name="Dandelion"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 9),
                   title="All Mutation Tracks Combined")
plt.savefig(os.path.join(OUT_DIR, "mutation_combined_all_tracks.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_combined_all_tracks.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. LolliplotTrack with different shape types in separate tracks
# ---------------------------------------------------------------------------
tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(snps, features=feats, type="circle", name="Circles"),
    LolliplotTrack(snps, features=feats, type="pin", name="Pins"),
    LolliplotTrack(snps, features=feats, type="flag", name="Flags"),
]
axes = plot_tracks(tracks, region=region, figsize=(12, 9),
                   title="LolliplotTrack — Multiple Shape Types")
plt.savefig(os.path.join(OUT_DIR, "mutation_combined_multi_shape.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_combined_multi_shape.png")
plt.close("all")

print("[INFO] All combined figures generated successfully.")
