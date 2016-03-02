"""Functions for visulization distributions from seaborn"""
from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from ..util import freedman_diaconis_bins
from ..palette import circos, color_palette, blend_palette 

def hist2d(x, y, ax=None, bins=None, norm_hist=False,
           xlabel=None, ylabel=None, cmap=None,
           **kwargs):
    """
    Make a 2D histogram plot.

    Parameters
    ----------
    x, y: 1d array-like
        The two input datasets
    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.
    bins : argument for matplotlib hist(), or None, optional
        Specification of hist bins(integer or array-like), or None to use
        Freedman-Diaconis rule.
    norm_hist : bool, optional
        If True, the histogram height shows a density rather than a count.
    xlabel: string, or None, optional
        Set the x axis label of the current axis.
    ylabel: string, or None, optional
        Set the y axis label of the current axis.
    kwargs : key, value pairings, optional
        Other keyword arguments are passed to ``plt.hist2d()``

    Returns
    -------
    ax : matplotlib Axes
        Axes object with the plot.

    Examples
    --------

    """
    if ax is None:
        ax = plt.gca()

    kwargs.setdefault("normed", norm_hist)

    x = x.astype(np.float64)
    y = y.astype(np.float64)

    if bins is None:
        x_bins = freedman_diaconis_bins(x)
        y_bins = freedman_diaconis_bins(y)
        bins = min(50, int(np.mean([x_bins, y_bins])))
    
    if cmap is None:
        colors = [circos['meth'+str(100-i)] for i in range(101)]
        cmap = blend_palette(color_palette(colors), as_cmap=True)        

    ax.hist2d(x, y, bins=bins, cmap=cmap, **kwargs)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)

    return ax
