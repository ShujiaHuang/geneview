"""
Functions for visulization fastq data
"""
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt

from ..util import get_color_cycle
from ..io import Fastq, FastqReader

def fqqualplot(fqdata, phred=64, ax=None, title=None, 
               xlabel=None, ylabel=None, **kwargs):
    """
    Plotting fastq data quality distribution.

    Parameters
    ----------
    fqdata : Array like. Fastq type data in a list, array or Series
        The input fastq data for plot.

    phred : int, or 64, optional
        The phred value to convert the base quality from ASCII to int
        For Illumina fastq: 64(default)
        For Sanger data: 33

    ax : matplotlib axis, optional
        Axes to plot on, otherwise use current axis.

    title : string, or None, optional
        If not None, set title on the plot

    xlabel, ylabel : string,or None, optional
        Set the x/y axis label of the current axis.

    kwargs : key, values pairings, or None, optional
        Other keyword arguments are passed to boxplot in
        maplotlib.axis.Axes.boxplot.


    Returns
    -------
    ax : matplotlib Axes
        Axes object with the plot
    """
    if ax is None:
        ax = plt.gca()

    if len(fqdata) == 0:
        return ax

    print (kwargs)
    if 'showfliers' not in kwargs:
        kwargs.setdefault("showfliers", False)
    print (kwargs)

    data = []
    for r in fqdata:
        # Convert base quality from ASCII to be integer
        data.append([ord(b) - phred for b in r.qual])
    data = np.array(data)

    positions = [str(i) for i in range(1, len(data[0])+1)]
    ax.boxplot(data, labels=positions, **kwargs)

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    return ax


def fastqreport(fqfile, fqfilelist=None):
    """
    Create a report for fastq data by input one fastq file or 
    a fq file list.

    Parameters
    ----------

    Returns
    -------

    Examples
    --------

    """
    pass
