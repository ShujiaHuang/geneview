"""Cell Press journal plot style.

Derived from the *Cell Press Illustration Guidelines*:
https://www.cell.com/pb/assets/raw/shared/figure-guidelines/GA_guide_full.pdf

Key requirements implemented here:

* Sans-serif font (Arial or Helvetica).
* Minimum text size 6 pt; recommended 7 pt for axis labels.
* Figure widths: single column = 85 mm (~3.35 in), 1.5 column = 114 mm
  (~4.49 in), double column = 174 mm (~6.85 in).
* Clean, minimal design – no background shading, no gratuitous 3-D.
* Accessible, high-contrast colour palette.
* TrueType font embedding (``pdf.fonttype = 42``).
* 300 dpi minimum for photographic images.

"""
from ._core import PlotStyle, register_style

# ---------------------------------------------------------------------------
# Cell Press accessible colour palette
# A curated set inspired by the Cell / Lancet colour guidelines, extended
# with Tol's bright qualitative palette for colour-blind safety.
# ---------------------------------------------------------------------------
_CELL_PALETTE = [
    "#0072B2",  # strong blue
    "#D55E00",  # vermilion
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

# ---------------------------------------------------------------------------
# Figure dimensions recommended by Cell Press
# Single column:  85 mm  ≈ 3.35 in
# 1.5 column:    114 mm  ≈ 4.49 in
# Double column: 174 mm  ≈ 6.85 in
# ---------------------------------------------------------------------------
_CELL_SINGLE_COL = (3.35, 2.50)  # width × height in inches

# ---------------------------------------------------------------------------
# Style definition
# ---------------------------------------------------------------------------
CELL_STYLE = PlotStyle(
    name="cell",
    description=(
        "Cell Press journal style – compact, accessible figures sized "
        "for single- or double-column layouts."
    ),

    # -- Fonts ---------------------------------------------------------------
    font_family="sans-serif",
    font_sans_serif=["Arial", "Helvetica", "DejaVu Sans"],
    # Cell recommends 7 pt for labels; 6 pt minimum for all text.
    font_size_title=8.0,
    font_size_label=7.0,
    font_size_tick=6.0,
    font_size_legend=6.0,

    # -- Figure --------------------------------------------------------------
    figure_dpi=300,
    figure_facecolor="white",
    figure_figsize=_CELL_SINGLE_COL,

    # -- Axes ----------------------------------------------------------------
    axes_linewidth=0.4,
    axes_labelpad=3.0,
    # Cell figures typically show only bottom + left axes.
    axes_spines_top=False,
    axes_spines_right=False,
    axes_spines_bottom=True,
    axes_spines_left=True,
    axes_grid=False,

    # -- Ticks ---------------------------------------------------------------
    xtick_major_size=2.5,
    ytick_major_size=2.5,
    xtick_major_width=0.4,
    ytick_major_width=0.4,
    xtick_minor_size=1.2,
    ytick_minor_size=1.2,
    xtick_minor_width=0.3,
    ytick_minor_width=0.3,
    xtick_direction="out",
    ytick_direction="out",

    # -- Lines / markers / scatter ------------------------------------------
    lines_linewidth=0.5,
    lines_markersize=2.5,
    scatter_alpha=0.85,

    # -- Colour palette ------------------------------------------------------
    color_palette=_CELL_PALETTE,

    # -- Export --------------------------------------------------------------
    pdf_fonttype=42,
    ps_fonttype=42,
    savefig_dpi=300,
    savefig_bbox="tight",

    # -- Escape hatch --------------------------------------------------------
    rc_params={
        "svg.fonttype": "none",
    },
)


# Register on import
register_style(CELL_STYLE)
