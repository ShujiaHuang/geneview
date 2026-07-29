"""Colour helpers for the geneview base plotting engine.

Plotting functions historically re-implemented two small colour chores:

* turning a comma-separated colour string (``"#3B5488,#53BBD5"``) or a plain
  colour string (``"rb"``) into an infinite iterator to cycle through, and
* turning a palette specification (colormap name / list / ``Colormap``) into a
  concrete list of ``n`` colours.

Both are collected here so that every module resolves colours the same way.
"""
from itertools import cycle
from typing import List

from matplotlib.colors import is_color_like

from ..palette import generate_colors_palette


def color_cycle(color):
    """Return an infinite iterator cycling through *color*.

    Mirrors the historical manhattan/qq behaviour, with one robustness guard:

    * a string containing commas is split into a list of colours
      (``"#3B5488,#53BBD5"`` -> ``["#3B5488", "#53BBD5"]``);
    * a comma-free string that matplotlib recognises as a *single* colour is
      cycled as one element (``"#FF0000"`` / ``"red"`` -> that colour, repeated)
      instead of being split character by character;
    * any other comma-free string is cycled **character by character**
      (``"rb"`` -> ``"r", "b", "r", ...``), matching ``itertools.cycle`` on a
      bare string;
    * a list / tuple is cycled element by element.

    Parameters
    ----------
    color : str, list, or tuple
        The colour specification to cycle through.

    Returns
    -------
    itertools.cycle
        An infinite iterator over the resolved colours.
    """
    if isinstance(color, str):
        if "," in color:
            color = [c.strip() for c in color.split(",")]
        elif is_color_like(color):
            # A single colour (e.g. "#FF0000", "red") must not be split into
            # characters; treat it as a one-element cycle.
            color = [color]
    return cycle(color)


def resolve_colors(palette, n_colors: int, alpha: float = 1.0) -> List:
    """Resolve a palette specification into a concrete list of colours.

    Thin wrapper around :func:`geneview.palette.generate_colors_palette` so
    that callers share a single palette-resolution path.

    Parameters
    ----------
    palette : str, list, or matplotlib.colors.Colormap
        A colormap name, an explicit list of colours, or a ``Colormap``.
    n_colors : int
        Number of colours to generate.
    alpha : float, optional, default: 1.0
        Alpha blending value applied to the generated colours.

    Returns
    -------
    list
        A list of colours (length may be shorter than ``n_colors`` if the
        underlying palette cannot supply enough distinct colours).
    """
    return generate_colors_palette(cmap=palette, n_colors=n_colors, alpha=alpha)
