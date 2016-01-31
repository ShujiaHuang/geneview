"""
    Plotting functions for manhattan plot.

    Copytright (c) Shujia Huang
    Date: 2016-01-23

    This model is based on brentp's script on github:
    https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py

    Thanks for Brentp's contributions

"""
from __future__ import print_function, division

from itertools import groupby, cycle
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np

##
from ..util import chr_id_cmp, despine
from ..palette import color_palette


def manhattanplot(data, ax=None, color=None, kind='scatter', 
                  xtick_label_set=None, CHR=None, alpha=0.8, 
                  mlog10=True, **kwargs):
    """
    Plot a manhattan plot.

    Parameters
    ----------
    data : float. 2D list or 2D numpy array. format [[id, x_val, y_val], ...]
        Input data for plot manhattan.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    color : matplotlib color, optional, default: color_palette('colorful', 4) 
        Color used for the plot elements. Could hex-code or rgb, 
        e.g: '#000000,#969696' or 'rb'

    kind : {'scatter' | 'line'}, optional
        Kind of plot to draw

    xtick_label_set : a set. 
        The x-labels for x-axis to draw in the figure

    CHR : str, optional, defualt: None
        Choice the specific chromosome to plot. And the x-axis will be the
        position of this chromosome instead of the chromosome id.

        CAUSION: this parameter could not be used with ``xtick_label_set``
                 together.

    alpha : scalar, optional, default: 0.8
        The alpha blending value, between 0(transparent) and 1(opaque)

    mlog10 : bool, optional, default: True
        If true, -log10 of the y_value(always be the p-value) is plotted. It
        isn't very useful to plot raw p-values, but plotting the raw value 
        could be useful for other genome-wide plots, for example peak heights,
        bayes factors, test statistics, other "scores", etc.

    kwargs : key, value pairings
        Other keyword arguments are passed to ``plt.scatter()`` or
        ``plt.vlines()`` (in matplotlib.pyplot) depending on whether 
        a scatter or line plot is being drawn.


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the manhattanplot.


    Notes
    -----
    1. This plot function is not just suit for GWAS manhattan plot,
       it could also be used for all the input data which format is ::

        [ [id1, x-value1, y-value1],
          [id2, x-value2, y-value2],
          ...
        ]

    2. The right and top spines of the plot have been setted to be 
       invisible by default.

    3. I'm going to add a parameter calls ``highlight`` to highlight a
       set of interesting positions (SNPs). And this parameter takes a 
       list-like value.

    """
    if CHR is not None and xtick_label_set is not None:
        msg = "``CHR`` and ``xtick_label_set`` can't be setted simultaneously."
        raise ValueError(msg)

    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()

    if color is None:
        color = color_palette("colorful", 4) 
        #color = color_palette("husl", 4) 

    if ',' in color: color = color.split(',')
    colors = cycle(color)
    
    last_x = 0
    xs_by_id = {} # use for collecting chromosome's position on x-axis
    x, y, c = [], [], []
    for seqid, rlist in groupby(data, key=itemgetter(0)):

        if CHR is not None and seqid != CHR: continue
            
        color = colors.next()
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        x.extend(region_xs)
        y.extend([r[2] for r in rlist])
        c.extend([color] * len(rlist))

        # ``xs_by_id`` is for setting up positions and ticks. Ticks should
        # be placed in the middle of a chromosome. The a new pos column is 
        # added that keeps a running sum of the positions of each successive 
        # chromsome.
        xs_by_id[seqid] = (region_xs[0] + region_xs[-1]) / 2
        last_x = x[-1]  # keep track so that chrs don't overlap in the plot.

    if not x:
        msg = ("zero-size array to reduction operation minimum which has no "
               "identity. This could be caused by zero-size array of ``x`` "
               "in the ``manhattanplot(...)`` function.")
        raise ValueError(msg)

    c = np.array(c)
    x = np.array(x)
    y = -np.log10(y) if mlog10 else np.array(y)

    if kind == 'scatter':
        ax.scatter(x, y, c=c, alpha=alpha, edgecolors='none', **kwargs)

    elif kind == 'line':
        ax.vlines(x, 0, y, colors=c, alpha=alpha, **kwargs)

    else:
        msg = "``kind`` must be either 'scatter' or 'line'"
        raise ValueError(msg)

    # Remove the 'top' and 'right' plot spines by default
    despine(ax=ax)

    ax.set_xlim(0, x[-1])
    ax.set_ylim(ymin=y.min())

    if xtick_label_set is None: 
        xtick_label_set = set(xs_by_id.keys())

    if CHR is None:
        xs_by_id = [(k, xs_by_id[k])
                     for k in sorted(xs_by_id.keys(), cmp=chr_id_cmp)
                     if k in xtick_label_set]

        ax.set_xticks([1] + [c[1] for c in xs_by_id])
        ax.set_xticklabels([''] + [c[0] for c in xs_by_id])

    else:
        # show the whole chromsome's position without scientific notation
        # if you are just interesting in this chromosome.
        ax.get_xaxis().get_major_formatter().set_scientific(False)

    return ax
