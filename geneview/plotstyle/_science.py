"""Science (AAAS) journal plot style.

Derived from the *Science* journal figure preparation guidelines:
https://www.science.org/content/page/preparing-initial-manuscript

Key requirements implemented here:

* Sans-serif font (Arial or Helvetica).
* Minimum text size 6 pt; labels 7–8 pt; titles up to 10 pt.
* Figure widths: single column = 2⅜ in (~60 mm), 1.5 column = 3¼ in
  (~82 mm), double column = 5 in (~127 mm).
* Clean, uncluttered design – no background shading, no 3-D effects.
* High-contrast, colour-blind-accessible palette (Okabe–Ito).
* TrueType font embedding for editable text in PDF/EPS exports.
* 300 dpi minimum for photographs; 600 dpi for line art.

"""
from ._core import PlotStyle, register_style

# ---------------------------------------------------------------------------
# Okabe–Ito colour-blind-safe palette
# https://jfly.uni-koeln.de/color/
# ---------------------------------------------------------------------------
_OKABE_ITO_PALETTE = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#CC79A7",  # reddish purple
    "#000000",  # black
]

# ---------------------------------------------------------------------------
# Figure dimensions recommended by Science
# Single column:  2 3/8 in  ≈  60 mm
# 1.5 column:     3 1/4 in  ≈  82 mm
# Double column:  5 in      ≈ 127 mm
# ---------------------------------------------------------------------------
_SCIENCE_SINGLE_COL = (2.36, 2.00)  # width × height in inches

# ---------------------------------------------------------------------------
# Style definition
# ---------------------------------------------------------------------------
SCIENCE_STYLE = PlotStyle(
    name="science",
    description=(
        "Science (AAAS) journal style – compact, high-contrast figures "
        "suitable for single- or double-column layouts."
    ),

    # -- Fonts ---------------------------------------------------------------
    font_family="sans-serif",
    font_sans_serif=["Arial", "Helvetica", "DejaVu Sans"],
    # Titles may be up to 10 pt; axis labels 7–8 pt; tick/legend ≥ 6 pt.
    font_size_title=10.0,
    font_size_label=8.0,
    font_size_tick=6.0,
    font_size_legend=6.0,

    # -- Figure --------------------------------------------------------------
    figure_dpi=300,
    figure_facecolor="white",
    figure_figsize=_SCIENCE_SINGLE_COL,

    # -- Axes ----------------------------------------------------------------
    axes_linewidth=0.4,
    axes_labelpad=3.0,
    # Clean design: only bottom + left spines.
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
    xtick_minor_size=1.0,
    ytick_minor_size=1.0,
    xtick_minor_width=0.3,
    ytick_minor_width=0.3,
    xtick_direction="out",
    ytick_direction="out",

    # -- Lines / markers / scatter ------------------------------------------
    lines_linewidth=0.5,
    lines_markersize=2.5,
    scatter_alpha=0.85,

    # -- Colour palette ------------------------------------------------------
    color_palette=_OKABE_ITO_PALETTE,

    # -- Export --------------------------------------------------------------
    pdf_fonttype=42,
    ps_fonttype=42,
    savefig_dpi=600,   # Science: 600 dpi for line art
    savefig_bbox="tight",

    # -- Escape hatch --------------------------------------------------------
    rc_params={
        "svg.fonttype": "none",
    },
)


# Register on import
register_style(SCIENCE_STYLE)
