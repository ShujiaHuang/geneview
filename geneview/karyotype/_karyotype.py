"""
Plotting functions for create karyotype plots.

Copyright (c) Shujia Huang
Date: 2016-02-19
"""
from __future__ import division, print_function

from pandas import DataFrame
import pandas as pd
import matplotlib.pyplot as plt

from ..palette import circos  # ``circos`` is a color dict

def karyoplot(data, ax=None, CHR=None, alpha=0.8, color4none='#34728B', 
              **kwargs):
    """ Create karyotype plot.

    Parameters
    ----------
    data : string or array
        A karyotype information list or input file path, even more could be 
        any kind of URL link. e.g. AWS S3 link.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    CHR : string, optional, defualt: None
        Choice the specific chromosome to plot.

    alpha : scalar, optional, default: 0.8   
        The alpha blending value, between 0(transparent) and 1(opaque)

    color4none : matplotlib color, optional, default: '#34728B'(deep gray blue)
        The color for undefine band color of karyotype in the plot.

    Examples
    --------

    A basic karyotype plot get the input karyotype information from URL:

    .. plot::
        :context: close-figs

        >>> import geneview as gv
        >>> gv.karyoplot('https://github.com/ShujiaHuang/geneview-data/raw/'
        ...              'master/karyotype/karyotype_human_hg19.txt')

    """
    # Draw the plot and return the Axes 
    if ax is None:
        ax = plt.gca()

    if isinstance(data, string):
        # suppose to be a path to the input file or a url to the file
        data = pd.read_table(
            data, header=0, 
            names=['chrom', 'start', 'stop', 'name', 'gie_stain'])
    elif isinstance(data, DataFrame):
        # reset the columns
        data = DataFrame(
            data.values, 
            columns=['chrom', 'start', 'stop', 'name', 'gie_stain'])
    else:
        # convert to DataFrame of pandas
        data = DataFrame(
            data.values, 
            columns=['chrom', 'start', 'stop', 'name', 'gie_stain'])

    """
    # remove 'chr'
    if CHR: 
        CHR = CHR[3:] if CHR.startswith('chr') else CHR
    data['chrom'] = [c[3:] if c.startswith('chr') else c 
                     for c in data['chrom']]
    """

    return ax
