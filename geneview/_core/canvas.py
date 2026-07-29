"""Canvas helpers for the geneview base plotting engine.

This module centralises the *axes acquisition* boilerplate that used to be
duplicated across every plotting function (``if ax is None: subplots(...)``)
plus the follow-up step of enforcing the active :class:`~geneview.plotstyle.PlotStyle`
spine rules on the resulting axes.

The single public entry point is :func:`get_or_create_axes`.  It is consumed
directly by the :func:`geneview._core.decorators.styled_plot` decorator, and
may also be called by hand from any plotting function that needs finer control
(for example a function that draws several sub-axes).

"""
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt

from ..plotstyle import PlotStyle, get_active_style


def get_or_create_axes(
    ax=None,
    *,
    figsize: Optional[Tuple[float, float]] = None,
    style: Union[str, PlotStyle, None] = None,
    apply_spines: bool = False,
    use_gca: bool = False,
    **subplot_kws,
):
    """Return a matplotlib ``Axes``, creating one if the caller passed ``None``.

    This is the shared replacement for the ``if ax is None: _, ax = subplots(...)``
    idiom.  When it creates the axes it honours the ``figsize`` requested by the
    caller (falling back to the active style's ``figure_figsize`` when both the
    argument and the style are available), and it can optionally enforce the
    active style's spine-visibility rules on the axes.

    Parameters
    ----------
    ax : matplotlib Axes or None, optional
        An existing axes to draw on.  When ``None`` a new axes is created.
    figsize : tuple of float, optional
        ``(width, height)`` in inches for the newly created figure.  Ignored
        when ``ax`` is provided or when ``use_gca`` is True.  When ``None`` and
        a style is resolvable, the style's ``figure_figsize`` is used.
    style : str, PlotStyle, or None, optional
        Style used to resolve spine rules / default figure size.  When ``None``
        the currently active style (see :func:`geneview.plotstyle.get_active_style`)
        is used.  Callers usually leave this ``None`` and rely on being inside a
        ``use_style(...)`` context.
    apply_spines : bool, optional, default: False
        When True, apply the resolved style's spine-visibility flags to the
        axes via :meth:`PlotStyle.apply_to_axes`.
    use_gca : bool, optional, default: False
        When True and ``ax`` is ``None``, reuse matplotlib's current axes
        (``plt.gca()``) instead of creating a brand new figure.  Preserves the
        historical behaviour of functions such as ``karyoplot``.
    **subplot_kws
        Extra keyword arguments forwarded to ``matplotlib.pyplot.subplots``
        (e.g. ``facecolor``, ``edgecolor``, ``constrained_layout``).

    Returns
    -------
    ax : matplotlib Axes
        The existing or newly created axes.
    """
    resolved_style = _resolve(style)

    if ax is None:
        if use_gca:
            ax = plt.gca()
        else:
            if figsize is None and resolved_style is not None:
                figsize = resolved_style.figure_figsize
            _, ax = plt.subplots(figsize=figsize, **subplot_kws)

    if apply_spines and resolved_style is not None:
        resolved_style.apply_to_axes(ax)

    return ax


def _resolve(style: Union[str, PlotStyle, None]) -> Optional[PlotStyle]:
    """Turn a style argument into a PlotStyle, defaulting to the active one."""
    if style is None:
        return get_active_style()
    if isinstance(style, PlotStyle):
        return style
    # A style *name*: reuse the registry lookup exposed by plotstyle.
    from ..plotstyle import get_style
    return get_style(style)
