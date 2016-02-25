"""
Functions for plotting qq-plot.

Copytright (c) Shujia Huang
Date: 2016-01-28

"""
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

from ..util import is_numeric

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
        msg = "`a` can just be any float value in 0 < a < 1."
        raise ValueError(msg) 

    try:
        n = np.float(len(n))
    except TypeError:
        n = np.float(n)

    # ``n`` is a float value
    return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


def qqplot(data, other=None, ax=None, xlabel=None, ylabel=None, color=None, 
           ablinecolor='r', alpha=0.8, mlog10=True, **kwargs):
    """Creat Q-Q plot.
    **CAUSION: The x-axis(expected) is created from uniform distribution.**

    Parameters
    ----------
    data : array-like, ``Series`` of ``pandas`` or 1d array-like
        Data to be plotted

    other : array-like, ``Series`` of ``pandas`` or 1d array, or None, optional
        If provided, the sample quantiles of the `data` array-like object are 
        plotted against the sample quantiles of the `other` array-like object. 
        If not provided (default), the theoretical quantiles are used.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    xlabel: string, or None, optional
        Set the x axis label of the current axis.
        CAUSION: The x axis will always be the expected value.

    ylabel: string, or None, optional
        Set the y axis label of the current axis.
        CAUSION: The y axis will always be the observed value.

    color : matplotlib color, optional
        The dots color in the plot

    ablinecolor: matplotlib color, default is 'r' (red), optional
        Color for the abline in plot. if set ``ablinecolor=None`` 
        means do not plot the abline.

    alpha : float scalar, default is 0.8, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    mlog10 : bool, default is 'True', optional
        If true, -log10 of the y_value(always be the p-value) is plotted.
        It isn't very useful to plot raw p-values in GWAS QQ plot.

    kwargs : key, value pairings
        Other keyword arguments are passed to ``plt.scatter()``
        (in matplotlib.pyplot).


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the manhattanplot.


    Notes
    -----
    1. The X axis will always be the expected values and Y axis always be
       observed values in the plot.
    2. This plot function is not just suit for GWAS QQ plot, it could
       also be used for creating QQ plot for other data, which format 
       are list-like ::
        [value1, value2, ...] (all the values should between 0 and 1)


    Examples
    --------

    Plot a basic QQ plot:

    .. plot::
        :context: close-figs

        >>> import geneview as gv
        >>> df = gv.util.load_dataset('GOYA_preview')
        >>> gv.gwas.qqplot(df['pvalue'], 
        ...                xlabel="Expected p-value(-log10)",
        ...                ylabel="Observed p-value(-log10)")

    We could even create a QQ plot base on two different dataset:

    .. plot::
        :context: close-figs

        >>> import numpy as np
        >>> data1 = np.random.normal(size = 100)
        >>> data2 = np.random.normal(5.0, 1.0, size = 100)
        >>> gv.gwas.qqplot(data1, other=data2, mlog10=False, 
        ...                xlabel="Expected", ylabel="Observe")
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
        xlabel = 'Expected(-log10)' if other is None else '-Log10(value) of 2nd Sample'
    if ylabel is None:
        ylabel = 'Observed(-log10)' if other is None else '-Log10(value) of 1st Sample'

    data = np.array(data, dtype=float)
    # create observed and expected
    e = ppoints(len(data)) if other is None else sorted(other)
    if mlog10:
        o = -np.log10(sorted(data))
        e = -np.log10(e)
    else:
        o = np.array(sorted(data))
        e = np.array(e)

    ax = _do_plot(e, o, ax=ax, color=color, ablinecolor=ablinecolor, 
                  alpha=alpha, **kwargs)

    ax.set_xlabel(xlabel) 
    ax.set_ylabel(ylabel) 

    return ax


def qqnorm(data, ax=None, xlabel='Expected', ylabel='Observed', 
           color=None, ablinecolor='r', alpha=0.8, **kwargs):
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
        Axes object with the manhattanplot.


    Notes
    -----
    1. The X axis will always be the expected values and Y axis always be
       observed values in the plot.


    Examples
    --------

    Plot a basic QQ norm plot:

    .. plot::
        :context: close-figs

        >>> import numpy as np
        >>> import geneview as gv
        >>> data = np.random.normal(size = 100)
        >>> gv.gwas.qqnorm(data) 

    Plot a QQ norm plot with GOYA_preview data:

    .. plot::
        :context: close-figs

        >>> import geneview as gv
        >>> df = gv.util.load_dataset('GOYA_preview')
        >>> gv.gwas.qqnorm(df['pvalue'], 
        ...                xlabel="Expected value",
        ...                ylabel="Observed value")
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


def _do_plot(x, y, ax=None, color=None, ablinecolor='r', alpha=0.8, **kwargs):
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
        Axes object with the manhattanplot.
    """
    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()

    # Get the color from the current color cycle
    if color is None:
        line, = ax.plot(0, x.mean())
        color = line.get_color()
        line.remove()

    # x is for expected; y is for observed value
    ax.scatter(x, y, c=color, alpha=alpha, edgecolors='none', **kwargs)
    if ablinecolor:
        # plot the y=x line by expected: uniform distribution data
        ax.plot([x.min(), x.max()], [x.min(), x.max()], color=ablinecolor,
                linestyle='-')

    ax.set_xlim(xmin=x.min(), xmax=1.05 * x.max())
    ax.set_ylim(ymin=y.min())

    return ax

