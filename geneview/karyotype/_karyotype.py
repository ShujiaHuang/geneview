"""
Plotting functions for create karyotype plots.

Copyright (c) Shujia Huang
Date: 2016-02-19
"""
import numpy as np
from pandas import DataFrame
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ..palette import circos  # ``circos`` is a color dict


def karyoplot(data, ax=None, width=0.5, CHR=None, alpha=0.8, color4none="#34728B", **kwargs):
    """ Create karyotype plot.

    Parameters
    ----------
    data : string or array
        A karyotype information list or input file path, even more could be 
        any kind of URL link. e.g. AWS S3 link.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    width : float, optional, default: 0.5
        Chromosom"s width in the plot

    CHR : string, optional, defualt: None
        Choice the specific chromosome to plot.

    alpha : scalar, optional, default: 0.8   
        The alpha blending value, between 0(transparent) and 1(opaque)

    color4none : matplotlib color, optional, default: "#34728B"(deep gray blue)
        The color for undefine band color of karyotype in the plot.

    kwargs : key, value pairings
        Other keyword arguments are passed to ``Rectangle`` in matplotlib.patches

    Examples
    --------

    A basic karyotype plot get the input karyotype information from URL:

    .. plot::
        :context: close-figs

        >>> import matplotlib.pyplot as plt
        >>> from geneview.utils import load_dataset
        >>> from geneview import karyoplot
        >>> fig, ax = plt.subplots(figsize=(20, 5))
        >>> k_fn = load_dataset("karyotype_human_hg19.txt")
        >>> _ = karyoplot(k_fn, ax=ax)

    """
    # Draw the plot and return the Axes 
    if ax is None:
        ax = plt.gca()

    if isinstance(data, str):
        # suppose to be a path to the input file or a url to the file
        data = pd.read_table(data, header=0,
                             names=["chrom", "start", "end", "name", "gie_stain"])
    elif isinstance(data, DataFrame):
        # reset the columns
        data = DataFrame(data.values, columns=["chrom", "start", "end", "name", "gie_stain"])
    else:
        # convert to DataFrame of pandas
        data = DataFrame(data, columns=["chrom", "start", "end", "name", "gie_stain"])

    yaxis = []
    for i, (chrom, kc_df) in enumerate(sorted(data.groupby("chrom"), key=lambda x: x[0])):

        if CHR is not None and chrom != CHR:
            continue

        yaxis.append(chrom)
        for _, r in kc_df.iterrows():
            band_color = circos[r.gie_stain] if r.gie_stain in circos else color4none
            band_rec = Rectangle((r.start, i), r.end - r.start, width,
                                 color=band_color, **kwargs)
            ax.add_patch(band_rec)

    xmax = data["end"].max() * 1.1
    xticks = np.arange(0, xmax, xmax / 10.)
    ax.set_xticks(xticks)
    ax.set_xticklabels(["{0}M".format(int(i / 10 ** 6)) for i in xticks])
    ax.set_xlim(0, xmax)

    ax.set_yticks([i + width / 2 for i in range(len(yaxis))])
    ax.set_yticklabels(yaxis)
    ax.set_ylim(0, len(yaxis))

    return ax
