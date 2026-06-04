"""Default 'geneview' plot style.

This style merges the historical hardcoded settings from
``geneview/__init__.py`` (rcParams lines) with the common defaults
observed across the plotting functions (manhattan, qq, venn, admixture).

It is designed to produce clean, publication-ready figures out of the
box while remaining fully customisable via journal-specific overrides.

"""
from ._core import PlotStyle, register_style

# ---------------------------------------------------------------------------
# Accessible colour palette
# ---------------------------------------------------------------------------
# A perceptually distinct, colour-blind-friendly palette suitable for
# general genomics visualisation.  The first two entries match the
# historical geneview manhattan colours (#3B5488 / #53BBD5) so that
# existing tutorials and examples remain visually consistent.
_DEFAULT_PALETTE = [
    "#3B5488",  # steel blue  (geneview legacy primary)
    "#53BBD5",  # sky blue    (geneview legacy secondary)
    "#E6553A",  # vermilion
    "#31A354",  # green
    "#756BB1",  # purple
    "#E6B422",  # gold
    "#A6761D",  # brown
    "#E377C2",  # pink
    "#17BECF",  # cyan
    "#9467BD",  # violet
]

# ---------------------------------------------------------------------------
# Style definition
# ---------------------------------------------------------------------------
DEFAULT_STYLE = PlotStyle(
    name="geneview",
    description=(
        "Default geneview style – clean, readable genomics figures "
        "suitable for exploration and general-purpose publication."
    ),

    # -- Fonts ---------------------------------------------------------------
    # Matches the historical rcParams from geneview/__init__.py:
    #   font.sans-serif = ["Arial","Lucida Sans","DejaVu Sans",
    #                      "Lucida Grande","Verdana"]
    font_family="sans-serif",
    font_sans_serif=[
        "Arial", "Lucida Sans", "DejaVu Sans", "Lucida Grande", "Verdana"
    ],
    # Larger than strict journal limits – comfortable for screen and poster
    # use.  Journal-specific styles override these with tighter constraints.
    font_size_title=12.0,
    font_size_label=10.0,
    font_size_tick=10.0,
    font_size_legend=9.0,

    # -- Figure --------------------------------------------------------------
    figure_dpi=300,
    figure_facecolor="white",
    # 9 × 3 inches matches the historical manhattanplot default.  Individual
    # plot functions (qq, venn, admixture) override this when they create
    # their own Figure/Axes, so this is only a fallback.
    figure_figsize=(9.0, 3.0),

    # -- Axes ----------------------------------------------------------------
    axes_linewidth=0.8,
    axes_labelpad=4.0,
    # Top / right spines hidden – matches manhattanplot and common
    # genomics plotting convention.
    axes_spines_top=False,
    axes_spines_right=False,
    axes_spines_bottom=True,
    axes_spines_left=True,
    axes_grid=False,

    # -- Ticks ---------------------------------------------------------------
    xtick_major_size=3.5,
    ytick_major_size=3.5,
    xtick_major_width=0.8,
    ytick_major_width=0.8,
    xtick_minor_size=2.0,
    ytick_minor_size=2.0,
    xtick_minor_width=0.5,
    ytick_minor_width=0.5,
    xtick_direction="out",
    ytick_direction="out",

    # -- Lines / markers / scatter ------------------------------------------
    lines_linewidth=1.0,
    lines_markersize=4.0,
    scatter_alpha=0.8,

    # -- Colour palette ------------------------------------------------------
    color_palette=_DEFAULT_PALETTE,

    # -- Genome tracks -------------------------------------------------------
    # Default style uses comfortable sizes and the historical Gviz-like
    # title panel (lightgray background, white bold text).
    tracks_title_bg="#D3D3D3",
    tracks_title_color="white",
    tracks_title_fontsize=10.0,
    tracks_title_border="none",
    tracks_axis_color="#A9A9A9",
    tracks_axis_linewidth=0.8,
    tracks_tick_fontsize=7.0,
    tracks_feature_fontsize=10.0,
    tracks_linewidth=1.0,
    tracks_figsize_width=12.0,
    tracks_height_per_track=1.2,

    # -- Export --------------------------------------------------------------
    # TrueType embedding (type 42) is required by almost all journals so
    # that text remains editable in Illustrator / Inkscape.
    pdf_fonttype=42,
    ps_fonttype=42,
    savefig_dpi=300,
    savefig_bbox="tight",
)


# Register on import
register_style(DEFAULT_STYLE)
