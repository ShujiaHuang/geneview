#!/usr/bin/env python
"""
Palette comparison — the *same* figure rendered with different color palettes.

Multi-series data is where the color palette matters most: geneview needs to
tell several overlaid series apart at a glance.  This example draws one
identical multi-series ``DataTrack`` and renders it under each built-in style
so the palettes can be compared directly:

  * ``geneview`` (default) — no categorical palette, so every series shares a
    single color (hard to distinguish; motivates the journal palettes).
  * ``nature``  — Wong colorblind-safe palette.
  * ``science`` — Okabe-Ito colorblind-safe palette.
  * ``cell``    — Cell colorblind-safe palette.

Outputs:
  * genome_tracks_palette_<style>.png  — one figure per palette
  * genome_tracks_palette_comparison.png — 2x2 side-by-side montage (a/b/c/d)

Run:  python examples/scripts/genome_tracks_palette_comparison.py
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.plotstyle import get_style
from geneview.genometracks import (
    GenomeAxisTrack, DataTrack, GenomicInterval, plot_tracks,
    save_figure, add_panel_labels,
)

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# The palettes we compare.  Each entry: (style, human-readable palette name).
STYLES = [
    ("geneview", "default (single color)"),
    ("nature", "Wong"),
    ("science", "Okabe-Ito"),
    ("cell", "Cell"),
]

# ---------------------------------------------------------------------------
# One shared dataset — four overlaid series (identical across every palette).
# ---------------------------------------------------------------------------
rng = np.random.RandomState(2024)
n = 120
starts = np.linspace(26_500_000, 26_795_000, n, dtype=int)
data = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": starts,
    "end": starts + 2500,
    "Sample A": rng.randn(n).cumsum() / 4 + 2,
    "Sample B": rng.randn(n).cumsum() / 4,
    "Sample C": rng.randn(n).cumsum() / 4 - 2,
    "Sample D": rng.randn(n).cumsum() / 4 + 4,
})
region = GenomicInterval("chr7", 26_500_000, 26_800_000)


# ---------------------------------------------------------------------------
# 1. Render the same figure once per palette.
# ---------------------------------------------------------------------------
paths = []
for style_name, palette_name in STYLES:
    axes = plot_tracks(
        [
            GenomeAxisTrack(),
            DataTrack(data, type="line", name="4 series", lwd=1.6),
        ],
        region=region,
        figsize=(9, 3.4),
        title=f"{style_name} — {palette_name}",
        style=style_name,
    )
    # Series carry labels from the DataFrame columns; add a legend so readers
    # can map each color to its sample (axes[1] is the DataTrack panel).
    axes[1].legend(loc="upper right", fontsize=6, ncol=4, framealpha=0.8)
    out = save_figure(
        axes,
        os.path.join(FIG_DIR, f"genome_tracks_palette_{style_name}.png"),
    )
    paths.append(out)
    plt.close("all")
    print(f"  -> {os.path.basename(out)}")


# ---------------------------------------------------------------------------
# 2. Assemble a single 2x2 montage so the palettes sit side-by-side.
# ---------------------------------------------------------------------------
fig, axgrid = plt.subplots(2, 2, figsize=(15, 8))
for ax, (style_name, palette_name), path in zip(axgrid.flat, STYLES, paths):
    ax.imshow(plt.imread(path))
    ax.set_axis_off()
    # Swatch row showing the palette actually used for the series.
    palette = get_style(style_name).tracks_categorical_palette or ["#0080FF"]
    for i, hexc in enumerate(palette[:8]):
        ax.add_patch(plt.Rectangle(
            (0.02 + i * 0.035, -0.06), 0.03, 0.035,
            transform=ax.transAxes, clip_on=False,
            facecolor=hexc, edgecolor="none",
        ))

# a / b / c / d panel labels (respects the montage figure's own default style).
add_panel_labels(list(axgrid.flat), x_offset=4.0, y_offset=-4.0)

fig.suptitle("Same multi-series figure under different color palettes",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.97))
montage = os.path.join(FIG_DIR, "genome_tracks_palette_comparison.png")
fig.savefig(montage, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  -> {os.path.basename(montage)}")

print("\n\u2713 Palette comparison figures generated in:", FIG_DIR)
