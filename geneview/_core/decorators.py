"""The ``@styled_plot`` decorator - heart of the geneview base plotting engine.

Every top-level plotting function in geneview shares the same wrapper logic:

1. read the ``style`` keyword and enter a ``use_style(style)`` context so the
   whole body renders under the requested (or active) style;
2. read the ``ax`` keyword and, when it is ``None``, create a figure/axes with
   a sensible default figure size;
3. optionally enforce the active style's spine-visibility rules on the axes.

``styled_plot`` factors that boilerplate out.  A decorated function is written
as if ``ax`` is always a real ``Axes`` and as if the style context is already
active - it only has to draw.  The public signature (including the ``ax`` and
``style`` keywords and the docstring) is preserved via :func:`functools.wraps`,
so nothing changes for callers.

Example
-------
::

    @styled_plot(figsize=(9, 3), subplot_kws={"facecolor": "w", "edgecolor": "k"})
    def manhattanplot(data, ..., ax=None, style=None, **kwargs):
        # ax is guaranteed non-None and the style context is active here
        ax.scatter(...)
        return ax
"""
import functools
import inspect
from typing import Optional, Tuple

import matplotlib.pyplot as plt

from ..plotstyle import use_style
from .canvas import get_or_create_axes


def styled_plot(
    func=None,
    *,
    figsize: Optional[Tuple[float, float]] = None,
    apply_spines: bool = True,
    use_gca: bool = False,
    subplot_kws: Optional[dict] = None,
):
    """Wrap a plotting function with geneview's shared style/canvas lifecycle.

    The wrapped function must expose (at least) an ``ax`` keyword argument;
    a ``style`` keyword argument is honoured when present.  At call time the
    decorator:

    * enters ``use_style(style)`` (a no-op when ``style is None``);
    * replaces ``ax=None`` with a freshly created axes (see
      :func:`geneview._core.canvas.get_or_create_axes`);
    * applies the active style's spine rules when ``apply_spines`` is True;
    * closes the auto-created figure if the body raises, so invalid input does
      not leak half-drawn figures.

    Parameters
    ----------
    func : callable, optional
        Supports both bare (``@styled_plot``) and parametrised
        (``@styled_plot(...)``) decorator syntax.
    figsize : tuple of float, optional
        Default ``(width, height)`` used when the function creates its own
        figure.  ``None`` falls back to the active style's ``figure_figsize``.
    apply_spines : bool, optional, default: True
        Whether to enforce the active style's spine-visibility flags on the
        axes.  Set ``False`` for plots that manage their own frame (e.g. an
        admixture bar plot that keeps all four spines, or a Venn diagram that
        turns the axis off entirely).
    use_gca : bool, optional, default: False
        When True, an omitted ``ax`` reuses ``plt.gca()`` instead of creating a
        new figure (preserves ``karyoplot``'s historical behaviour).
    subplot_kws : dict, optional
        Extra keyword arguments forwarded to ``matplotlib.pyplot.subplots`` when
        a new figure is created (e.g. ``facecolor``, ``constrained_layout``).

    Returns
    -------
    callable
        The wrapped plotting function with an unchanged public signature.
    """

    def decorator(fn):
        sig = inspect.signature(fn)
        if "ax" not in sig.parameters:
            raise TypeError(
                "@styled_plot can only decorate functions that expose an "
                "'ax' keyword argument; %r has none." % fn.__name__
            )

        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            style = bound.arguments.get("style", None)
            with use_style(style):
                created_fig = None
                incoming_ax = bound.arguments.get("ax", None)
                ax = get_or_create_axes(
                    incoming_ax,
                    figsize=figsize,
                    apply_spines=apply_spines,
                    use_gca=use_gca,
                    **(subplot_kws or {}),
                )
                if incoming_ax is None and not use_gca:
                    created_fig = ax.get_figure()
                bound.arguments["ax"] = ax

                try:
                    return fn(*bound.args, **bound.kwargs)
                except Exception:
                    # Do not leak the figure we just created on the error path.
                    if created_fig is not None:
                        plt.close(created_fig)
                    raise

        return wrapper

    # Support both @styled_plot and @styled_plot(...).
    if func is not None:
        return decorator(func)
    return decorator
