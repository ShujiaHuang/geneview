"""
Functions for plotting qq-plot.

Copytright (c) Shujia Huang
Date: 2016-01-28

"""
import numpy as np
import matplotlib.pyplot as plt

from geneview.util import is_numeric

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
    
    a : float, optional, default: 0.5
        The offset fraction to be used; typically in (0, 1)
        And more often: a = 3.0/8.0 if n <= 10 else 0.5 


    Returns
    -------
        A list of value calculates by this formular:
        (1:n - a)/(n + (1-2a)


    Notes
    -----
    * Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) The New S Language. Wadsworth & Brooks/Cole.
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

def qq(pvalues, ax=None, mlog10=True, **kwargs):
    """
    Creat a Q-Q plot.

    Creates a quantile-quantile plot from p-values from a GWAS study.

    Parameters
    ----------
    pvalues : float. list-like and just 1-D
        A numeric list or array of p-values.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    mlog10 : bool, optional, default: True 
        If true, -log10 of the y_value(always be the p-value) is plotted. It
        isn't very useful to plot raw p-values in QQ plot.

    kwargs : key, value pairings
        Other keyword arguments are passed to ``plt.scatter()``
        (in matplotlib.pyplot) depending on being drawn.


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the manhattanplot.


    Notes
    1. This plot function is not just suit for GWAS QQ plot, it could
       also be used for creating QQ plot for other data, which format 
       are list-like ::
        [value1, value2, ...] (all the value should between 0 and 1)

    -----
    """
    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()
    
    if not all(map(is_numeric, pvalues)):
        msg = 'Input must all be numeric.'
        raise ValueError(msg)

    # Observed and expected
    if mlog10:
        o = -np.log10(sorted(pvalues, reverse=True))
        e = -np.log10(ppoints(len(pvalues))) # Why not reverse?! 
    else:
        o = sorted(pvalues)
        e = ppoints(len(pvalues))

    return ax





