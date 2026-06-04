"""Plot styles for geneview.

This sub-package provides a registry of pre-defined plot styles that bring
geneview figures into compliance with the requirements of major scientific
journals (Nature, Science, Cell) or the geneview default.

Quick start
-----------

Apply a journal style globally::

    from geneview.plotstyle import apply_style
    apply_style("nature")      # all subsequent plots use Nature style

Apply a style temporarily via a context manager::

    from geneview.plotstyle import use_style
    with use_style("science"):
        ax = manhattanplot(data=df)

Pass a style directly to a plotting function::

    ax = manhattanplot(data=df, style="cell")

List all available styles::

    from geneview.plotstyle import list_styles
    print(list_styles())
    # ['cell', 'geneview', 'nature', 'science']

"""
# Core public API – must come first so that style modules can import it.
from ._core import (
    PlotStyle,
    apply_style,
    use_style,
    get_style,
    list_styles,
    register_style,
)

# Importing each style module triggers ``register_style(...)`` at module
# load time, populating the global registry.
from . import _default   # noqa: F401  (registers "geneview")
from . import _nature    # noqa: F401  (registers "nature")
from . import _science   # noqa: F401  (registers "science")
from . import _cell      # noqa: F401  (registers "cell")

__all__ = [
    "PlotStyle",
    "apply_style",
    "use_style",
    "get_style",
    "list_styles",
    "register_style",
]
