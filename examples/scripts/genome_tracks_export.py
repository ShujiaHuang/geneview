"""
export_tracks Example
=====================

Demonstrates exporting track data to BED, GFF, bedGraph, and WIG formats.
"""
import os
import tempfile
import pandas as pd
import numpy as np

from geneview.genometracks import (
    AnnotationTrack, DataTrack, GenomicInterval, export_tracks,
)

# ---------------------------------------------------------------------------
# 1. Export AnnotationTrack to BED
# ---------------------------------------------------------------------------
ann_data = pd.DataFrame({
    "chrom": ["chr7"] * 3,
    "start": [1000, 2000, 3000],
    "end": [1500, 2800, 3600],
    "name": ["geneA", "geneB", "geneC"],
    "strand": ["+", "-", "+"],
})
ann_track = AnnotationTrack(ann_data, name="Genes")

bed_path = os.path.join(tempfile.gettempdir(), "example_export.bed")
export_tracks(ann_track, bed_path, fmt="bed")
print(f"Exported BED: {bed_path}")
with open(bed_path) as f:
    print(f.read())

# ---------------------------------------------------------------------------
# 2. Export AnnotationTrack to GFF
# ---------------------------------------------------------------------------
gff_path = os.path.join(tempfile.gettempdir(), "example_export.gff")
export_tracks(ann_track, gff_path, fmt="gff", source="geneview_example")
print(f"Exported GFF: {gff_path}")
with open(gff_path) as f:
    print(f.read())

# ---------------------------------------------------------------------------
# 3. Export DataTrack to bedGraph
# ---------------------------------------------------------------------------
rng = np.random.RandomState(42)
n = 20
starts = np.linspace(1000, 2000, n, dtype=int)
data_values = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts,
    "end": starts + 50,
    "value": rng.randn(n).cumsum(),
})
data_track = DataTrack(data_values, type="line")

bdg_path = os.path.join(tempfile.gettempdir(), "example_export.bdg")
export_tracks(data_track, bdg_path, fmt="bedgraph")
print(f"Exported bedGraph: {bdg_path}")
with open(bdg_path) as f:
    print(f.read()[:500])

# ---------------------------------------------------------------------------
# 4. Export DataTrack to WIG
# ---------------------------------------------------------------------------
wig_path = os.path.join(tempfile.gettempdir(), "example_export.wig")
export_tracks(data_track, wig_path, fmt="wig")
print(f"Exported WIG: {wig_path}")
with open(wig_path) as f:
    print(f.read()[:500])

# ---------------------------------------------------------------------------
# 5. Export a raw DataFrame
# ---------------------------------------------------------------------------
df = pd.DataFrame({
    "chrom": ["chr1", "chr1"],
    "start": [100, 500],
    "end": [200, 600],
    "value": [1.5, 2.3],
})
raw_path = os.path.join(tempfile.gettempdir(), "example_export_raw.bdg")
export_tracks(df, raw_path, fmt="bedgraph")
print(f"Exported raw DataFrame: {raw_path}")

print("All export_tracks examples complete.")

# ---------------------------------------------------------------------------
# 6. save_figure() — convenience figure export
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from geneview.genometracks import (
    GenomeAxisTrack, GenomicInterval, plot_tracks, save_figure,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")

# Create a simple track figure
ann_track2 = AnnotationTrack(ann_data, name="Export Demo")
region = GenomicInterval("chr7", 500, 4000)
axes = plot_tracks([GenomeAxisTrack(), ann_track2], region=region,
                   figsize=(12, 4))

# Save using save_figure (auto-detects format from extension)
png_path = os.path.join(FIG_DIR, "genome_tracks_save_figure.png")
save_figure(axes, png_path, dpi=150)
print(f"\nSaved PNG: {png_path}")

# save_figure also supports PDF, SVG, EPS
pdf_path = os.path.join(tempfile.gettempdir(), "track_figure.pdf")
save_figure(axes, pdf_path)
print(f"Saved PDF: {pdf_path}")

svg_path = os.path.join(tempfile.gettempdir(), "track_figure.svg")
save_figure(axes, svg_path)
print(f"Saved SVG: {svg_path}")

# You can also override the format explicitly
save_figure(axes, os.path.join(tempfile.gettempdir(), "track_eps"), fmt="eps")
print("Saved EPS (explicit fmt)")

plt.close("all")
print("\nAll save_figure examples complete.")
