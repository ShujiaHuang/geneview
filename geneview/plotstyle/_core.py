"""Core style infrastructure for geneview plot styles.

Defines the PlotStyle dataclass, a global style registry, and helpers for
applying styles globally or as context managers.

"""
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Union

import matplotlib

# ---------------------------------------------------------------------------
# PlotStyle dataclass
# ---------------------------------------------------------------------------

@dataclass
class PlotStyle:
    """Complete style specification for a geneview plot.

    All dimension fields map (directly or via ``style_to_rcparams``) to
    matplotlib rcParams.  Spine visibility is handled separately because
    matplotlib has no rcParam for it; plot functions can read these flags
    from the active style and call ``ax.spines[...].set_visible()``.

    Parameters
    ----------
    name : str
        Short identifier used in the registry (e.g. ``"nature"``).
    description : str
        One-line human-readable description.
    font_family : str
        Matplotlib font family (``"sans-serif"``, ``"serif"``, etc.).
    font_sans_serif : list of str
        Ordered fallback list for the ``sans-serif`` family.
    font_size_title : float
        Font size in pt for axes titles.
    font_size_label : float
        Font size in pt for axes labels.
    font_size_tick : float
        Font size in pt for tick labels.
    font_size_legend : float
        Font size in pt for legend text.
    figure_dpi : int
        Default DPI for screen display.
    figure_facecolor : str
        Background colour of the figure canvas.
    figure_figsize : tuple of float
        Default ``(width, height)`` in inches for a single panel.
    axes_linewidth : float
        Width of the axes spines / box lines.
    axes_labelpad : float
        Padding between axis label and tick labels (pt).
    axes_spines_top : bool
        Whether the top spine is visible.
    axes_spines_right : bool
        Whether the right spine is visible.
    axes_spines_bottom : bool
        Whether the bottom spine is visible.
    axes_spines_left : bool
        Whether the left spine is visible.
    axes_grid : bool
        Whether to show background grid lines.
    xtick_major_size, ytick_major_size : float
        Length of major tick marks (pt).
    xtick_major_width, ytick_major_width : float
        Width of major tick marks (pt).
    xtick_minor_size, ytick_minor_size : float
        Length of minor tick marks (pt).
    xtick_minor_width, ytick_minor_width : float
        Width of minor tick marks (pt).
    xtick_direction, ytick_direction : str
        Tick direction (``"out"``, ``"in"``, or ``"inout"``).
    lines_linewidth : float
        Default line width for ``plot()`` calls.
    lines_markersize : float
        Default marker size.
    scatter_alpha : float
        Default alpha for scatter plots.
    color_palette : list
        Ordered list of hex colours forming the journal's accessible palette.
        An empty list means "use matplotlib's current colour cycle".
    pdf_fonttype : int
        Value for ``pdf.fonttype`` (42 = TrueType, required by most journals).
    ps_fonttype : int
        Value for ``ps.fonttype``.
    savefig_dpi : int
        DPI used when saving figures.
    savefig_bbox : str
        Bounding-box strategy for ``savefig`` (``"tight"`` or ``"standard"``).
    rc_params : dict
        Escape hatch: arbitrary extra rcParams merged last, overriding any
        field above.
    """

    name: str
    description: str = ""

    # -- Fonts ---------------------------------------------------------------
    font_family: str = "sans-serif"
    font_sans_serif: List[str] = field(
        default_factory=lambda: ["Arial", "Helvetica", "DejaVu Sans"]
    )
    font_size_title: float = 7.0
    font_size_label: float = 7.0
    font_size_tick: float = 6.0
    font_size_legend: float = 6.0

    # -- Figure --------------------------------------------------------------
    figure_dpi: int = 300
    figure_facecolor: str = "white"
    figure_figsize: tuple = (3.5, 2.5)

    # -- Axes ----------------------------------------------------------------
    axes_linewidth: float = 0.5
    axes_labelpad: float = 4.0
    axes_spines_top: bool = False
    axes_spines_right: bool = False
    axes_spines_bottom: bool = True
    axes_spines_left: bool = True
    axes_grid: bool = False

    # -- Ticks ---------------------------------------------------------------
    xtick_major_size: float = 3.0
    ytick_major_size: float = 3.0
    xtick_major_width: float = 0.5
    ytick_major_width: float = 0.5
    xtick_minor_size: float = 1.5
    ytick_minor_size: float = 1.5
    xtick_minor_width: float = 0.4
    ytick_minor_width: float = 0.4
    xtick_direction: str = "out"
    ytick_direction: str = "out"

    # -- Lines / markers / scatter ------------------------------------------
    lines_linewidth: float = 0.5
    lines_markersize: float = 3.0
    scatter_alpha: float = 0.8

    # -- Colour palette ------------------------------------------------------
    color_palette: List[str] = field(default_factory=list)

    # -- Export --------------------------------------------------------------
    pdf_fonttype: int = 42
    ps_fonttype: int = 42
    savefig_dpi: int = 300
    savefig_bbox: str = "tight"

    # -- Escape hatch --------------------------------------------------------
    rc_params: Dict = field(default_factory=dict)

    # --------------------------------------------------------------------- #
    # Derived helpers                                                        #
    # --------------------------------------------------------------------- #

    def to_rc_params(self) -> Dict:
        """Convert this style into a dict suitable for ``matplotlib.rcParams.update()``.

        Spine visibility flags are **not** included (they have no rcParam
        equivalent) and must be applied to each ``Axes`` object separately.
        """
        params: Dict = {
            # Fonts
            "font.family": self.font_family,
            "font.sans-serif": self.font_sans_serif,
            "axes.titlesize": self.font_size_title,
            "axes.labelsize": self.font_size_label,
            "xtick.labelsize": self.font_size_tick,
            "ytick.labelsize": self.font_size_tick,
            "legend.fontsize": self.font_size_legend,
            # Figure
            "figure.dpi": self.figure_dpi,
            "figure.facecolor": self.figure_facecolor,
            "figure.figsize": self.figure_figsize,
            # Axes
            "axes.linewidth": self.axes_linewidth,
            "axes.labelpad": self.axes_labelpad,
            "axes.grid": self.axes_grid,
            # Ticks
            "xtick.major.size": self.xtick_major_size,
            "ytick.major.size": self.ytick_major_size,
            "xtick.major.width": self.xtick_major_width,
            "ytick.major.width": self.ytick_major_width,
            "xtick.minor.size": self.xtick_minor_size,
            "ytick.minor.size": self.ytick_minor_size,
            "xtick.minor.width": self.xtick_minor_width,
            "ytick.minor.width": self.ytick_minor_width,
            "xtick.direction": self.xtick_direction,
            "ytick.direction": self.ytick_direction,
            # Lines / markers
            "lines.linewidth": self.lines_linewidth,
            "lines.markersize": self.lines_markersize,
            # Export / font embedding
            "pdf.fonttype": self.pdf_fonttype,
            "ps.fonttype": self.ps_fonttype,
            "savefig.dpi": self.savefig_dpi,
            "savefig.bbox": self.savefig_bbox,
        }

        # Colour cycle override
        if self.color_palette:
            params["axes.prop_cycle"] = matplotlib.cycler(color=self.color_palette)

        # Escape-hatch overrides last
        params.update(self.rc_params)
        return params

    def apply_to_axes(self, ax) -> None:
        """Apply spine visibility settings to a matplotlib ``Axes`` object.

        Call this after creating the axes to enforce the style's spine rules.
        """
        ax.spines["top"].set_visible(self.axes_spines_top)
        ax.spines["right"].set_visible(self.axes_spines_right)
        ax.spines["bottom"].set_visible(self.axes_spines_bottom)
        ax.spines["left"].set_visible(self.axes_spines_left)


# ---------------------------------------------------------------------------
# Style registry
# ---------------------------------------------------------------------------

_STYLES_REGISTRY: Dict[str, PlotStyle] = {}


def register_style(style: PlotStyle) -> None:
    """Register a ``PlotStyle`` instance in the global registry.

    Raises
    ------
    TypeError
        If *style* is not a ``PlotStyle`` instance.
    """
    if not isinstance(style, PlotStyle):
        raise TypeError(f"Expected a PlotStyle instance, got {type(style).__name__}")
    _STYLES_REGISTRY[style.name] = style


def get_style(name: str) -> PlotStyle:
    """Return the registered ``PlotStyle`` with the given *name*.

    Raises
    ------
    ValueError
        If no style with that name is registered.
    """
    if name not in _STYLES_REGISTRY:
        available = ", ".join(sorted(_STYLES_REGISTRY)) or "(none)"
        raise ValueError(
            f"Unknown plot style '{name}'. Available styles: {available}"
        )
    return _STYLES_REGISTRY[name]


def list_styles() -> List[str]:
    """Return a sorted list of all registered style names."""
    return sorted(_STYLES_REGISTRY)


# ---------------------------------------------------------------------------
# Apply / context-manager helpers
# ---------------------------------------------------------------------------

def _resolve_style(name_or_style: Union[str, PlotStyle]) -> PlotStyle:
    """Accept either a name string or a PlotStyle and return a PlotStyle."""
    if isinstance(name_or_style, PlotStyle):
        return name_or_style
    if isinstance(name_or_style, str):
        return get_style(name_or_style)
    raise TypeError(
        f"Expected a style name (str) or PlotStyle instance, "
        f"got {type(name_or_style).__name__}"
    )


def apply_style(name_or_style: Union[str, PlotStyle]) -> PlotStyle:
    """Apply a style **globally** by updating ``matplotlib.rcParams``.

    Parameters
    ----------
    name_or_style : str or PlotStyle
        A registered style name (e.g. ``"nature"``) or a ``PlotStyle`` object.

    Returns
    -------
    PlotStyle
        The resolved style that was applied.
    """
    style = _resolve_style(name_or_style)
    matplotlib.rcParams.update(style.to_rc_params())
    return style


@contextmanager
def use_style(name_or_style: Union[str, PlotStyle, None]):
    """Context manager that temporarily applies a style.

    On exit the previous ``matplotlib.rcParams`` are automatically restored
    (via ``matplotlib.rc_context``).

    Parameters
    ----------
    name_or_style : str, PlotStyle, or None
        A registered style name, a ``PlotStyle`` object, or ``None``.
        When ``None`` (the default), no style change is applied – this
        allows plot functions to write ``with use_style(style):`` without
        a separate ``if style is not None`` branch.

    Yields
    ------
    PlotStyle or None
        The resolved style, or ``None`` when no style was requested.

    Example
    -------
    >>> with use_style("nature"):
    ...     ax = manhattanplot(data=df)
    """
    if name_or_style is None:
        yield None
        return
    style = _resolve_style(name_or_style)
    with matplotlib.rc_context(style.to_rc_params()):
        yield style
