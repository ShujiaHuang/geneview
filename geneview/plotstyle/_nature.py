"""Nature journal plot style.

Derived from the official *Nature Research Figure Guide*:
https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications

Key requirements implemented here:

* Sans-serif font (Arial / Helvetica), **not** outlined (``pdf.fonttype=42``).
* Maximum text size 7 pt, minimum 5 pt; panel labels 8 pt bold.
* No background gridlines, no drop shadows, no patterns.
* Accessible colour palette based on Bang Wong's *Points of View* palette
  (Nature Methods 8, 441, 2011) – safe for most forms of colour blindness.
* Axes must have tick marks and be labelled with units in parentheses.
* No coloured text – keylines + black labels only.
* RGB colour space; vector artwork export (PDF/EPS).
* Figure widths: single column ≈ 89 mm (3.5 in), 1.5 column ≈ 120 mm
  (4.7 in), double column ≈ 183 mm (7.2 in).

"""
from ._core import PlotStyle, register_style

# ---------------------------------------------------------------------------
# Wong colour-blind-safe palette (Nature Methods 8, 441, 2011)
# https://www.nature.com/articles/nmeth.1618
# ---------------------------------------------------------------------------
_WONG_PALETTE = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#CC79A7",  # reddish purple
]

# ---------------------------------------------------------------------------
# Figure dimensions recommended by Nature
# Single column:  ~89 mm  ≈ 3.50 in
# 1.5 column:    ~120 mm  ≈ 4.72 in
# Double column: ~183 mm  ≈ 7.20 in
# ---------------------------------------------------------------------------
_NATURE_SINGLE_COL = (3.50, 2.60)  # width × height in inches

# ---------------------------------------------------------------------------
# Style definition
# ---------------------------------------------------------------------------
NATURE_STYLE = PlotStyle(
    name="nature",
    description=(
        "Nature Research journal style – compact, accessible, vector-ready "
        "figures compliant with the Nature Research Figure Guide."
    ),

    # -- Fonts ---------------------------------------------------------------
    font_family="sans-serif",
    font_sans_serif=["Arial", "Helvetica", "DejaVu Sans"],
    # Max 7 pt body text; panel labels (a, b, c…) handled separately as
    # 8 pt bold by the user/annotation code.
    font_size_title=7.0,
    font_size_label=7.0,
    font_size_tick=5.5,
    font_size_legend=5.5,

    # -- Figure --------------------------------------------------------------
    figure_dpi=300,
    figure_facecolor="white",
    figure_figsize=_NATURE_SINGLE_COL,

    # -- Axes ----------------------------------------------------------------
    axes_linewidth=0.4,
    axes_labelpad=3.0,
    # Nature figures typically show only bottom + left spines.
    axes_spines_top=False,
    axes_spines_right=False,
    axes_spines_bottom=True,
    axes_spines_left=True,
    # "We avoid: background gridlines"
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
    lines_linewidth=0.4,
    lines_markersize=2.5,
    scatter_alpha=0.85,

    # -- Colour palette ------------------------------------------------------
    color_palette=_WONG_PALETTE,

    # -- Genome tracks -------------------------------------------------------
    # Nature: compact, minimal. White title panels with dark text.
    # Smaller fonts per the 5–7 pt rule.
    tracks_title_bg="white",
    tracks_title_color="#333333",
    tracks_title_fontsize=7.0,
    tracks_title_border="lightgray",
    tracks_axis_color="#666666",
    tracks_axis_linewidth=0.4,
    tracks_tick_fontsize=5.5,
    tracks_feature_fontsize=6.0,
    tracks_linewidth=0.4,
    tracks_figsize_width=7.2,   # Nature double-column ≈ 183 mm
    tracks_height_per_track=0.8,

    # -- Export --------------------------------------------------------------
    # TrueType embedding (type 42) – required so text stays editable.
    pdf_fonttype=42,
    ps_fonttype=42,
    savefig_dpi=450,   # Nature recommends 450 dpi for raster images
    savefig_bbox="tight",

    # -- Escape hatch --------------------------------------------------------
    # Nature asks for RGB colour space; matplotlib's default is already RGB.
    rc_params={
        "svg.fonttype": "none",   # keep SVG text editable
    },
)


# Register on import
register_style(NATURE_STYLE)
