"""
Plotting functions for create karyotype plots.

Copyright (c) Shujia Huang
Date: 2016-02-19
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ..palette import get_cytoband_color
from ..utils import read_cytoband


def _chrom_sort_key(chrom):
    """Sort chromosomes in natural biological order (chr1, chr2, ..., chr10, chrX, etc.)."""
    name = str(chrom)
    # Strip common prefixes for numeric extraction
    stripped = name.replace("chr", "").replace("Chr", "").replace("CHR", "")
    try:
        return (0, int(stripped), name)
    except ValueError:
        # Non-numeric chromosomes (X, Y, MT, etc.) sort after numeric ones
        return (1, 0, stripped)


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

    # Normalize any supported input (file/URL, DataFrame, or array of rows)
    # into the canonical chrom/chromStart/chromEnd/name/gieStain schema.
    data = read_cytoband(data)

    yaxis = []
    row = 0
    max_end = 0
    for chrom, kc_df in sorted(data.groupby("chrom"), key=lambda x: _chrom_sort_key(x[0])):

        if CHR is not None and chrom != CHR:
            continue

        yaxis.append(chrom)
        for _, r in kc_df.iterrows():
            band_color = get_cytoband_color(r.gieStain, default=color4none)
            # Use facecolor (not ``color``) so callers can override the band
            # border via ``edgecolor``/``linewidth`` kwargs; default the edge to
            # the fill color to keep the historical borderless appearance.
            band_style = dict(kwargs)
            band_style.setdefault("edgecolor", band_color)
            band_rec = Rectangle((r.chromStart, row), r.chromEnd - r.chromStart, width,
                                 facecolor=band_color, alpha=alpha, **band_style)
            ax.add_patch(band_rec)

        max_end = max(max_end, int(kc_df["chromEnd"].max()))
        row += 1

    # Scale the x-axis to the plotted chromosome(s); fall back to the full
    # dataset if ``CHR`` matched nothing so we never divide by zero below.
    xmax = (max_end if max_end else int(data["chromEnd"].max())) * 1.1
    xticks = np.arange(0, xmax, xmax / 10.)
    ax.set_xticks(xticks)
    ax.set_xticklabels(["{0}M".format(int(i / 10 ** 6)) for i in xticks])
    ax.set_xlim(0, xmax)

    ax.set_yticks([i + width / 2 for i in range(len(yaxis))])
    ax.set_yticklabels(yaxis)
    ax.set_ylim(0, len(yaxis))

    return ax
