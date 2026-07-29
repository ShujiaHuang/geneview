"""geneview base plotting engine (``geneview._core``).

This package holds the small, dependency-light primitives shared by every
geneview plotting function.  It sits *above* the infrastructure layer
(``plotstyle`` / ``palette`` / ``utils``) and *below* the individual plot
modules (``gwas`` / ``popgene`` / ``baseplot`` / ``karyotype`` / ...).

Public building blocks
----------------------
``styled_plot``
    Decorator that wraps a plotting function with the shared style + canvas
    lifecycle (enter ``use_style``, create/acquire the axes, enforce spine
    rules).  This is the recommended entry point for new plot functions.
``get_or_create_axes``
    The ``if ax is None: subplots(...)`` replacement, with style-aware default
    figure size and optional spine enforcement.
``color_cycle`` / ``resolve_colors``
    Shared colour-resolution helpers.

New plotting modules should build on these primitives instead of
re-implementing the boilerplate, so that the whole toolkit stays consistent
and picks up future engine improvements for free.
"""
from .canvas import get_or_create_axes
from .colors import color_cycle, resolve_colors
from .decorators import styled_plot

__all__ = [
    "styled_plot",
    "get_or_create_axes",
    "color_cycle",
    "resolve_colors",
]
