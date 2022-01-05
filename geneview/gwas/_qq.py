"""Functions for plotting qq-plot.

Copytright (c) Shujia Huang
Date: 2021-02-21

"""
import numpy as np
from scipy.stats import norm, chi2
from matplotlib.pyplot import subplots

from ..utils import is_numeric


def ppoints(n, a=0.5):
    """ 
    numpy analogue `R`'s `ppoints` function.

    ppoints() is used in qqplot and qqnorm to generate the set of 
    probabilities at which to evaluate the inverse distribution.

    see details at 
        http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html 

    Parameters
    ----------
    n : a number or array type
        If n is an array or list, ``n`` will re-assign to ``n=len(n)``
    
    a : float, optional, default: 0.5
        The offset fraction to be used; typically in (0, 1)
        Recommend to set: a = 3.0/8.0 if n <= 10 else 0.5 


    Returns
    -------
        A list of value calculates by this formular:
        (1:n - a)/(n + (1-2a)
        Typically it's an array with elements in (0, 1)


    Notes
    -----
    * Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. 
      Wadsworth & Brooks/Cole.
    * Blom, G. (1958) Statistical Estimates and Transformed Beta Variables. Wiley

    """

    if a < 0 or a > 1:
        msg = "`a` could just be any float value in (0, 1)."
        raise ValueError(msg)

    try:
        n = len(n)
    except TypeError:
        pass

    # ``n`` is a float value
    return (np.arange(n, dtype=float) + 1 - a) / (n + 1 - 2 * a)


def qqplot(data, other=None, logp=True, ax=None, marker="o", color=None, alpha=0.8, title=None,
           xlabel=None, ylabel=None, ablinecolor="r", **kwargs):
    """Creat Q-Q plot.
    **CAUSION: The x-axis(expected) is created from uniform distribution.**

    Parameters
    ----------
    data : list, 1d-array-like, or Series
        Data to be plotted

    other : list, 1d-array-like, Series or None, optional
        If provided, the sample quantiles of the `data` array-like object are 
        plotted against the sample quantiles of the `other` array-like object. 
        If not provided (default), the theoretical quantiles are used.

    logp : bool, default is 'True', optional
        If true, -log10 of the y_value(always be the p-value) is plotted.
        It isn't very useful to plot raw p-values in GWAS QQ plot.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    marker : matplotlib markers for scatter plot, default is "o", optional

    color : matplotlib color, optional
        The dots color in the plot

    alpha : float scalar, default is 0.8, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    title : string, or None, optional
        Set the title of the current plot.

    xlabel : string, or None, optional
        Set the x axis label of the current axis.
        CAUSION: The x axis will always be the expected value.

    ylabel : string, or None, optional
        Set the y axis label of the current axis.
        CAUSION: The y axis will always be the observed value.

    ablinecolor : matplotlib color, default is 'r' (red), optional
        Color for the abline in plot. if set ``ablinecolor=None`` 
        means do not plot the abline.

    kwargs : key, value pairings, optional
        Other keyword arguments are passed to ``plt.scatter()``
        (in matplotlib.pyplot).


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the plot.

    Notes
    -----
    1. The X axis will always be the expected values and Y axis always be
       observed values in the plot.
    2. This plot function is not just suit for GWAS QQ plot, it could
       also be used for creating QQ plot for other data, which format 
       are list-like ::
        [value1, value2, ...] (all the values should between 0 and 1)

    See Also
    --------
    qqnorm : Q-Q plot against the normal distribution.

    Examples
    --------

    Plot a basic QQ plot:

    .. plot::
        :context: close-figs

        >>> import pandas as pd
        >>> from geneview import qqplot
        >>> from geneview.utils import load_dataset
        >>> df = load_dataset("gwas")
        >>> ax = qqplot(data=df["P"])

    Plot a better QQ plot with title and label texts.

    .. plot::
        :context: close-figs

        >>> from matplotlib.pyplot import subplots
        >>> fig, ax = subplots(figsize=(6, 6), facecolor="w", edgecolor="k")  # define a plot
        >>> ax = qqplot(data=df["P"],
        ...             marker="o",
        ...             title="Test",
        ...             xlabel=r"Expected $-log_{10}{(P)}$",
        ...             ylabel=r"Observed $-log_{10}{(P)}$",
        ...             ax=ax)

    We could even create a QQ plot base on two different dataset:

    .. plot::
        :context: close-figs

        >>> import numpy as np
        >>> data1 = np.random.normal(size = 100)
        >>> data2 = np.random.normal(5.0, 1.0, size = 100)
        >>> ax = qqplot(data=data1, other=data2, logp=False, xlabel="Expected", ylabel="Observe")
    """
    if not all(map(is_numeric, data)):
        msg = 'Input must all be numeric in `data`.'
        raise ValueError(msg)

    if other is not None and not all(map(is_numeric, other)):
        msg = 'Input must all be numeric in `other`.'
        raise ValueError(msg)

    if other is not None and len(other) != len(data):
        msg = 'Input `data` and `other` must all be the same size.'
        raise ValueError(msg)

    if xlabel is None:
        xlabel = r"$Expected(-log_{10}{(P)})$" if other is None else r"$-log_{10}{(Value)} of 2nd Sample$"
    if ylabel is None:
        ylabel = r"$Observed(-log_{10}{(P)})$" if other is None else r"(-log_{10}{(Value)}) of 1st Sample$"

    data = np.array(data, dtype=float)

    # create observed and expected
    e = ppoints(len(data)) if other is None else sorted(other)

    if logp:
        o = -np.log10(sorted(data))
        e = -np.log10(e)
    else:
        o = np.array(sorted(data))
        e = np.array(e)

    if "marker" not in kwargs:
        kwargs["marker"] = marker
    ax = _do_plot(e, o, ax=ax, color=color, ablinecolor=ablinecolor, alpha=alpha, **kwargs)

    expected_median = chi2.ppf(0.5, 1)  # This value is equal to 0.4549364
    lambda_value = round(np.median(norm.ppf(data) ** 2) / expected_median, 3)
    if title:
        title += r"$(\lambda = %s)$" % lambda_value
    else:
        title = r"$\lambda = %s$" % lambda_value

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


def qqnorm(data, ax=None, xlabel="Expected normal distribution", ylabel="Observed distribution", 
           color=None, ablinecolor="r", alpha=0.8, **kwargs):
    """Creat Q-Q plot against the normal distribution values.
    *CAUSION: The x-axis(expected) is created from normal distribution.*

    Parameters
    ----------
    data : array-like. ``Series`` of ``pandas`` or 1d array-like
        The data will be normalization be this formula:
        (data - data.mean()) / data.std() before to be plotted.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    xlabel: string, or 'Expected', optional
        Set the x axis label of the current axis.
        CAUSION: The x axis will always be the expected value.

    ylabel: string, or 'Observed', optional
        Set the y axis label of the current axis.
        CAUSION: The y axis will always be the observed value.

    color : matplotlib color, optional
        The dots color in the plot

    ablinecolor: matplotlib color, default is 'r' (red), optional
        Color for the abline in plot. if set ``ablinecolor=None`` 
        means do not plot the abline.

    alpha : float scalar, default is 0.8, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    kwargs : key, value pairings
        Other keyword arguments are passed to ``plt.scatter()``
        (in matplotlib.pyplot).


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the plot.


    Notes
    -----
    1. The X axis will always be the expected values and Y axis always be
       observed values in the plot.


    See Also
    --------
    qqplot : Q-Q plot against other values or uniform distribution(default)


    Examples
    --------

    Plot a basic QQ norm plot:

    .. plot::
        :context: close-figs

        >>> import numpy as np
        >>> data = np.random.normal(size = 100)
        >>> ax = qqnorm(data=data)

    Plot a QQ norm plot with GOYA_preview data:

    .. plot::
        :context: close-figs

        >>> import pandas as pd
        >>> from geneview.utils import load_dataset
        >>> df = load_dataset("gwas")
        >>> ax = qqnorm(data=df["P"], xlabel="Expected value", ylabel="Observed value")
    """
    if not all(map(is_numeric, data)):
        msg = 'Input must all be numeric in `data`.'
        raise ValueError(msg)

    # Normalization the data to be in (mu=0.0, std=1.0) normal distribution
    obs = np.array(data, dtype=float)
    obs = (obs - obs.mean()) / obs.std()
    obs.sort()

    # create expected
    e = norm.ppf(ppoints(len(obs)))

    ax = _do_plot(e, obs, ax=ax, color=color, ablinecolor=ablinecolor,
                  alpha=alpha, **kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax


def _do_plot(x, y, ax=None, color=None, ablinecolor="r", alpha=0.8, **kwargs):
    """
    Boiler plate plotting function for the `qqplot` and `qqnorm`

    Parameters
    ----------
    x, y : array-like
        Data to be plotted

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    color : matplotlib color, optional
        The dots color in the plot
        
    ablinecolor: matplotlib color, default is 'r' (red), optional
        Color for the abline in plot. if set ``ablinecolor=None`` 
        means do not plot the abline.

    alpha : float scalar, default is 0.8, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    kwargs : key, value pairings
        Other keyword arguments are passed to ``plt.scatter()``
        (in matplotlib.pyplot).

    Returns
    -------
    ax : matplotlib Axes
        Axes object with the plot.
    """
    # Draw the plot and return the Axes
    if ax is None:
        # ax = plt.gca()
        _, ax = subplots(figsize=(5, 5), facecolor="w", edgecolor="k")

    # Get the color from the current color cycle
    if color is None:
        line, = ax.plot(0, x.mean())
        color = line.get_color()
        line.remove()

    # x is for expected; y is for observed value
    ax.scatter(x, y, c=color, alpha=alpha, edgecolors='none', **kwargs)
    ax.set_xlim(xmin=x.min(), xmax=1.05 * x.max())
    ax.set_ylim(ymin=y.min())

    if ablinecolor:
        # plot the y=x line by expected: uniform distribution data
        ax.plot([x.min(), ax.get_xlim()[1]], [x.min(), ax.get_xlim()[1]],
                color=ablinecolor, linestyle="-")

    return ax
