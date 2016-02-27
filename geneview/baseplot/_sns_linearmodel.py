"""Plotting functions for linear models (broadly construed)."""
from __future__ import division
import matplotlib.pyplot as plt

from ..palette import color_palette
from ._sns_axisgrid import PairGrid
from ._sns_distribution import kdeplot


def pairplot(data, hue=None, hue_order=None, palette=None,
             vars=None, x_vars=None, y_vars=None,
             kind="scatter", diag_kind="hist", markers=None,
             size=2.5, aspect=1, dropna=True,
             plot_kws=None, diag_kws=None, grid_kws=None):
    """Plot pairwise relationships in a dataset.

    By default, this function will create a grid of Axes such that each
    variable in ``data`` will by shared in the y-axis across a single row and
    in the x-axis across a single column. The diagonal Axes are treated
    differently, drawing a plot to show the univariate distribution of the data
    for the variable in that column.

    It is also possible to show a subset of variables or plot different
    variables on the rows and columns.

    This is a high-level interface for :class:`PairGrid` that is intended to
    make it easy to draw a few common styles. You should use :class`PairGrid`
    directly if you need more flexibility.

    Parameters
    ----------
    data : DataFrame
        Tidy (long-form) dataframe where each column is a variable and
        each row is an observation.
    hue : string (variable name), optional
        Variable in ``data`` to map plot aspects to different colors.
    hue_order : list of strings
        Order for the levels of the hue variable in the palette
    palette : dict or geneview color palette
        Set of colors for mapping the ``hue`` variable. If a dict, keys
        should be values  in the ``hue`` variable.
    vars : list of variable names, optional
        Variables within ``data`` to use, otherwise use every column with
        a numeric datatype.
    {x, y}_vars : lists of variable names, optional
        Variables within ``data`` to use separately for the rows and
        columns of the figure; i.e. to make a non-square plot.
    kind : {'scatter', 'reg'}, optional
        Kind of plot for the non-identity relationships.
    diag_kind : {'hist', 'kde'}, optional
        Kind of plot for the diagonal subplots.
    markers : single matplotlib marker code or list, optional
        Either the marker to use for all datapoints or a list of markers with
        a length the same as the number of levels in the hue variable so that
        differently colored points will also have different scatterplot
        markers.
    size : scalar, optional
        Height (in inches) of each facet.
    aspect : scalar, optional
        Aspect * size gives the width (in inches) of each facet.
    dropna : boolean, optional
        Drop missing values from the data before plotting.
    {plot, diag, grid}_kws : dicts, optional
        Dictionaries of keyword arguments.

    Returns
    -------
    grid : PairGrid
        Returns the underlying ``PairGrid`` instance for further tweaking.

    See Also
    --------
    PairGrid : Subplot grid for more flexible plotting of pairwise
               relationships.

    Examples
    --------

    Draw scatterplots for joint relationships and histograms for univariate
    distributions:

    .. plot::
        :context: close-figs

        >>> import geneview as gv; gv.setup(style="ticks", color_codes=True)
        >>> iris = gv.util.load_dataset("iris")
        >>> g = gv.pairplot(iris)

    Show different levels of a categorical variable by the color of plot
    elements:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, hue="species")

    Use a different color palette:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, hue="species", palette="husl")

    Use different markers for each level of the hue variable:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, hue="species", markers=["o", "s", "D"])

    Plot a subset of variables:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, vars=["sepal_width", "sepal_length"])

    Draw larger plots:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, size=3,
        ...                 vars=["sepal_width", "sepal_length"])

    Plot different variables in the rows and columns:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris,
        ...                 x_vars=["sepal_width", "sepal_length"],
        ...                 y_vars=["petal_width", "petal_length"])

    Use kernel density estimates for univariate plots:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, diag_kind="kde")

    Fit linear regression models to the scatter plots:

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, kind="reg")

    Pass keyword arguments down to the underlying functions (it may be easier
    to use :class:`PairGrid` directly):

    .. plot::
        :context: close-figs

        >>> g = gv.pairplot(iris, diag_kind="kde", markers="+",
        ...                 plot_kws=dict(s=50, edgecolor="b", linewidth=1),
        ...                 diag_kws=dict(shade=True))

    """
    if plot_kws is None:
        plot_kws = {}
    if diag_kws is None:
        diag_kws = {}
    if grid_kws is None:
        grid_kws = {}

    # Set up the PairGrid
    diag_sharey = diag_kind == "hist"
    grid = PairGrid(data, vars=vars, x_vars=x_vars, y_vars=y_vars, hue=hue,
                    hue_order=hue_order, palette=palette,
                    diag_sharey=diag_sharey,
                    size=size, aspect=aspect, dropna=dropna, **grid_kws)

    # Add the markers here as PairGrid has figured out how many levels of the
    # hue variable are needed and we don't want to duplicate that process
    if markers is not None:
        if grid.hue_names is None:
            n_markers = 1
        else:
            n_markers = len(grid.hue_names)
        if not isinstance(markers, list):
            markers = [markers] * n_markers
        if len(markers) != n_markers:
            raise ValueError(("markers must be a singeton or a list of markers"
                              " for each level of the hue variable"))
        grid.hue_kws = {"marker": markers}

    # Maybe plot on the diagonal
    if grid.square_grid:
        if diag_kind == "hist":
            grid.map_diag(plt.hist, **diag_kws)
        elif diag_kind == "kde":
            diag_kws["legend"] = False
            grid.map_diag(kdeplot, **diag_kws)

    # Maybe plot on the off-diagonals
    if grid.square_grid and diag_kind is not None:
        plotter = grid.map_offdiag
    else:
        plotter = grid.map

    if kind == "scatter":
        plot_kws.setdefault("edgecolor", "white")
        plotter(plt.scatter, **plot_kws)
    elif kind == "reg":
        plotter(regplot, **plot_kws)

    # Add a legend
    if hue is not None:
        grid.add_legend()

    return grid
