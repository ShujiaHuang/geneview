"""Display Parameters & Gviz-compatible API Example
===================================================

Demonstrates the fine-tuned features inspired by Gviz:
- Constructor **kwargs as display parameters (Gviz convention)
- Dot-notation aliases ("col.title" → col_title)
- available_display_params() to query class-specific defaults
- GenomeAxisTrack with explicit ticks_at and exponent
- HighlightTrack in_background for z-order control
- reverse_stacking for AnnotationTrack y-order
- plot_tracks panel_only, margin, and inner_margin

Run:  python examples/scripts/genome_tracks_display_params.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    GenomeAxisTrack, AnnotationTrack, GeneRegionTrack, DataTrack,
    HighlightTrack, GenomicInterval, plot_tracks, available_display_params,
    read_bed, read_bedgraph, read_gff,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "genome_tracks")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Load data
cpg_data = read_bed(os.path.join(DATA_DIR, "cpg_islands.bed"))
cov_data = read_bedgraph(os.path.join(DATA_DIR, "coverage.bedgraph"))
gene_data = read_gff(os.path.join(DATA_DIR, "gene_models.gtf"))

region = GenomicInterval("chr7", 26_490_000, 26_720_000)

# ---------------------------------------------------------------------------
# 1. Constructor **kwargs as display parameters (Gviz convention)
# ---------------------------------------------------------------------------
# Instead of passing display_params={"col": "red", "fontsize": 14},
# you can pass them directly as keyword arguments:
atrack_styled = AnnotationTrack(
    cpg_data, name="CpG (styled)",
    col="darkgreen", fill="lightgreen", fontsize=14, alpha=0.8,
)

gtrack = GenomeAxisTrack()
axes1 = plot_tracks(
    [gtrack, atrack_styled], region=region,
    title="Constructor kwargs as Display Parameters",
    figsize=(14, 4),
)
fig1 = axes1[0].figure
fig1.savefig(os.path.join(OUT_DIR, "genome_tracks_kwargs_api.png"),
             dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_kwargs_api.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Dot-notation aliases (Gviz "col.title" → Python col_title)
# ---------------------------------------------------------------------------
# You can use either form:
atrack_alias1 = AnnotationTrack(
    cpg_data, name="Dot Alias",
    **{"col.title": "navy", "background.panel": "#FFF8DC"},
)
# is equivalent to:
atrack_alias2 = AnnotationTrack(
    cpg_data, name="Snake Alias",
    col_title="navy", background_panel="#FFF8DC",
)

fig_aliases, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 5))
plt.close(fig_aliases)
plot_tracks([atrack_alias1], region=region, figsize=(14, 2.5), ax=ax1)
ax1.set_title('Using "col.title" (dot notation)', fontsize=10)
plot_tracks([atrack_alias2], region=region, figsize=(14, 2.5), ax=ax2)
ax2.set_title('Using col_title (snake_case)', fontsize=10)
fig_aliases.tight_layout()
fig_aliases.savefig(os.path.join(OUT_DIR, "genome_tracks_alias_comparison.png"),
                    dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_alias_comparison.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. available_display_params() — query class-specific defaults
# ---------------------------------------------------------------------------
base_defaults = available_display_params()
axis_defaults = available_display_params("GenomeAxisTrack")
data_defaults = available_display_params("DataTrack")

print("\n--- available_display_params() ---")
print(f"Base defaults (excerpt): alpha={base_defaults['alpha']}, "
      f"fontsize={base_defaults['fontsize']}, col={base_defaults['col']}")
print(f"GenomeAxisTrack-specific: col_range={axis_defaults.get('col_range', 'N/A')}")
print(f"DataTrack-specific: col_histogram={data_defaults.get('col_histogram', 'N/A')}")

# ---------------------------------------------------------------------------
# 4. GenomeAxisTrack: explicit ticks_at and exponent
# ---------------------------------------------------------------------------
# 4a. Force labels in Mb units (exponent=6)
gtrack_exp = GenomeAxisTrack(exponent=6, name="Forced Mb")

# 4b. Explicit tick positions
tick_positions = [26_500_000, 26_550_000, 26_600_000, 26_650_000, 26_700_000]
gtrack_ticks = GenomeAxisTrack(ticks_at=tick_positions, name="Custom Ticks")

# Default for comparison
gtrack_default = GenomeAxisTrack(name="Default")

axes4 = plot_tracks(
    [gtrack_default, gtrack_exp, gtrack_ticks],
    region=region,
    title="GenomeAxisTrack: exponent and ticks_at",
    figsize=(14, 5),
)
fig4 = axes4[0].figure
fig4.savefig(os.path.join(OUT_DIR, "genome_tracks_axis_enhanced.png"),
             dpi=150, bbox_inches="tight")
print("\n[INFO] Saved genome_tracks_axis_enhanced.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. HighlightTrack in_background (z-order control)
# ---------------------------------------------------------------------------
hl_regions = pd.DataFrame({
    "chrom": ["chr7", "chr7"],
    "start": [26_520_000, 26_620_000],
    "end":   [26_560_000, 26_660_000],
})

# Background highlighting (zorder=1) — behind track content
atrack_bg = AnnotationTrack(cpg_data, name="CpG Islands")
htrack_bg = HighlightTrack(
    [atrack_bg], regions=hl_regions,
    fill="yellow", alpha=0.3, in_background=True,
)

# Foreground highlighting (zorder=10) — on top of track content
atrack_fg = AnnotationTrack(cpg_data, name="CpG Islands")
htrack_fg = HighlightTrack(
    [atrack_fg], regions=hl_regions,
    fill="yellow", alpha=0.3, in_background=False,
)

axes5a = plot_tracks(
    [GenomeAxisTrack(), htrack_bg], region=region,
    title="HighlightTrack: in_background=True (behind)",
    figsize=(14, 4),
)
fig5a = axes5a[0].figure
fig5a.savefig(os.path.join(OUT_DIR, "genome_tracks_highlight_background.png"),
              dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_highlight_background.png")
plt.close("all")

axes5b = plot_tracks(
    [GenomeAxisTrack(), htrack_fg], region=region,
    title="HighlightTrack: in_background=False (foreground)",
    figsize=(14, 4),
)
fig5b = axes5b[0].figure
fig5b.savefig(os.path.join(OUT_DIR, "genome_tracks_highlight_foreground.png"),
              dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_highlight_foreground.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 6. reverse_stacking — reverse y-order of stacked features
# ---------------------------------------------------------------------------
# Normal stacking (top-to-bottom)
atrack_normal = AnnotationTrack(
    cpg_data, name="Normal Stacking", stacking="squish",
)
# Reversed stacking (bottom-to-top)
atrack_reversed = AnnotationTrack(
    cpg_data, name="Reversed Stacking", stacking="squish",
    reverse_stacking=True,
)

axes6a = plot_tracks(
    [GenomeAxisTrack(), atrack_normal], region=region,
    title="Normal Stacking Order", figsize=(14, 4),
)
fig6a = axes6a[0].figure
fig6a.savefig(os.path.join(OUT_DIR, "genome_tracks_stacking_normal.png"),
              dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_stacking_normal.png")
plt.close("all")

axes6b = plot_tracks(
    [GenomeAxisTrack(), atrack_reversed], region=region,
    title="Reversed Stacking Order", figsize=(14, 4),
)
fig6b = axes6b[0].figure
fig6b.savefig(os.path.join(OUT_DIR, "genome_tracks_stacking_reversed.png"),
              dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_stacking_reversed.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 7. DataTrack-specific parameters: col_mountain, fill_mountain, etc.
# ---------------------------------------------------------------------------
rng = np.random.RandomState(42)
n = len(cov_data)
signal_data = pd.DataFrame({
    "chrom": ["chr7"] * n,
    "start": cov_data["start"].values,
    "end": cov_data["end"].values,
    "signal": cov_data["value"].values + rng.randn(n) * 0.5,
})

# Mountain plot with custom colors
dtrack_mountain = DataTrack(
    signal_data, type="mountain", name="Mountain",
    col_mountain="darkblue", fill_mountain="#ADD8E6", lwd_mountain=2,
    col_baseline="gray",
)

# Line plot with custom baseline
dtrack_line = DataTrack(
    signal_data, type="line", name="Line",
    col="darkred", col_baseline="#CCCCCC",
)

axes7 = plot_tracks(
    [GenomeAxisTrack(), dtrack_mountain, dtrack_line],
    region=region,
    title="DataTrack-Specific Parameters (col_mountain, fill_mountain, col_baseline)",
    figsize=(14, 6),
)
fig7 = axes7[0].figure
fig7.savefig(os.path.join(OUT_DIR, "genome_tracks_data_specific_params.png"),
             dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_data_specific_params.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 8. plot_tracks panel_only, margin, and inner_margin
# ---------------------------------------------------------------------------
# panel_only=True: render into existing axes without creating a new figure
fig8, axes8 = plt.subplots(2, 1, figsize=(14, 6))

# Create tracks
gaxis = GenomeAxisTrack()
dtrack = DataTrack(cov_data, type="histogram", name="Coverage")
atrack = AnnotationTrack(cpg_data, name="CpG Islands")

# Plot into existing axes using panel_only
plot_tracks([gaxis, dtrack], region=region, panel_only=True, ax=axes8[0])
plot_tracks([gaxis, atrack], region=region, panel_only=True, ax=axes8[1])

fig8.suptitle("panel_only=True: Embedding Tracks in Existing GridSpec", fontsize=14)
fig8.tight_layout()
fig8.savefig(os.path.join(OUT_DIR, "genome_tracks_panel_only.png"),
             dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_panel_only.png")
plt.close("all")

# margin and inner_margin (pixel-based spacing, Gviz-compatible)
axes9 = plot_tracks(
    [GenomeAxisTrack(), dtrack, atrack],
    region=region,
    title="Custom Margins (margin=12, inner_margin=6)",
    figsize=(14, 6),
    margin=12,
    inner_margin=6,
)
fig9 = axes9[0].figure
fig9.savefig(os.path.join(OUT_DIR, "genome_tracks_custom_margins.png"),
             dpi=150, bbox_inches="tight")
print("[INFO] Saved genome_tracks_custom_margins.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 9. Improved __repr__ shows non-default params
# ---------------------------------------------------------------------------
print("\n--- Improved __repr__ ---")
t1 = AnnotationTrack(cpg_data, name="Default")
print(f"Default:  {t1}")

t2 = AnnotationTrack(cpg_data, name="Custom", col="red", fontsize=16, alpha=0.5)
print(f"Custom:   {t2}")

print("\n[INFO] All display parameter example figures saved.")
