"""
    Copytright (c) Shujia Huang
    Date: 2016-01-23

    This model is based on brentp's script on github:
    https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py

    Plot a manhattan plot of the input file(s).
    python %prog [options] files
"""

import sys
import optparse
from itertools import groupby, cycle
from operator import itemgetter

from matplotlib import pyplot as plt
import numpy as np


def _chr_id_cmp(a, b):
    a = a.lower().replace('_', '')
    b = b.lower().replace('_', '')
    achr = a[3:] if a.startswith('chr') else a
    bchr = b[3:] if b.startswith('chr') else b

    try:
        return cmp(int(achr), int(bchr))
    except ValueError:
        if achr.isdigit() and not bchr.isdigit(): return -1
        if bchr.isdigit() and not achr.isdigit(): return 1
        # X Y
        return cmp(achr, bchr)


def manhattanplot(data, ax=None, color=None, mlog10=True, kind='scatter', 
                  xtick_label_set=None, alpha=0.8, **kwargs):
    """
    Plot a manhattan plot.

    Parameters
    ----------
    data : float. 2D list or 2D numpy array. format [[id, x_val, y_val], ...]
        Input data for plot manhattan.

    mlog10 : bool
        Set the y_value to be -log10 scale, optional, default: True 

    kind : {'scatter' | 'line'}, optional
        Kind of plot to draw

    color : matplotlib color, optional
        Color used for the plot elements. Could hex-code or rgb, 
        e.g: '#000000,#969696' or 'rb'

    xtick_label_set : a set. 
        The x-labels for x-axis to draw in the figure

    Returns
    -------
        ax : matplotlib Axes
            Axes object with the manhattanplot.
    """
    # Draw the plot and return the Axes
    if ax is None:
        ax = plt.gca()

    if color is None:
        color = '#6DC066,#FD482F,#8A2BE2,#3399FF' # It's colorful

    if ',' in color: color = color.split(',')
    colors = cycle(color)
    
    last_x = 0
    xs_by_id = {}
    x, y, c = [], [], []
    for seqid, rlist in groupby(data, key=itemgetter(0)):

        color = colors.next()
        rlist = list(rlist)
        region_xs = [last_x + r[1] for r in rlist]
        x.extend(region_xs)
        y.extend([r[2] for r in rlist])
        c.extend([color] * len(rlist))

        xs_by_id[seqid] = (region_xs[0] + region_xs[-1]) / 2

        # keep track so that chrs don't overlap.
        last_x = x[-1]

    c = np.array(c)
    x = np.array(x)
    y = -np.log10(y) if mlog10 else np.array(y)

    if kind == 'scatter':
        ax.scatter(x, y, s=20, c=c, alpha=alpha, edgecolors='none')

    elif kind == 'line':
        ax.vlines(x, 0, y, colors=c, alpha=alpha)

    else:
        msg = "``kind`` must be either 'scatter' or 'line'"
        raise ValueError(msg)

    ax.tick_params(labelsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_xlim(1, x[-1])
    ax.set_ylim(ymin=0)

    if xtick_label_set is None: 
        xtick_label_set = set(xs_by_id.keys())

    xs_by_id = [(k, xs_by_id[k])
                 for k in sorted(xs_by_id.keys(), cmp=_chr_id_cmp)
                 if k in xtick_label_set]

    ax.set_xticks([1] + [c[1] for c in xs_by_id])
    ax.set_xticklabels([''] + [c[0] for c in xs_by_id])

    return ax
