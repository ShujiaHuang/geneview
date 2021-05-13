"""Plotting functions for manhattan plot.

Copyright (c) Shujia Huang
Date: 2021-02-21

This model is based on brentp's script on github:
https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py

Thanks for Brentp's contributions

"""
from itertools import cycle
from pandas import DataFrame
import numpy as np

from matplotlib.pyplot import subplots
from ..utils import adjust_text


# learn something from "https://github.com/reneshbedre/bioinfokit/blob/38fb4966827337f00421119a69259b92bb67a7d0/bioinfokit/visuz.py"
def manhattanplot(data, chrom="#CHROM", pos="POS", pv="P", snp="ID", logp=True, ax=None,
                  marker=".", color="#3B5488,#53BBD5", alpha=0.8,
                  title=None, xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$",
                  xtick_label_set=None, CHR=None, xticklabel_kws=None,
                  suggestiveline=1e-5, genomewideline=5e-8, sign_line_cols="#D62728,#2CA02C", hline_kws=None,
                  sign_marker_p=None, sign_marker_color="r",
                  is_annotate_topsnp=False, text_kws=None, ld_block_size=50000, **kwargs):
    """Creates a manhattan plot from PLINK assoc output (or any data frame with chromosome, position, and p-value).

    Parameters
    ----------
    data : DataFrame.
        A DataFrame with columns "#CHROM," "POS," "P," and optionally, "SNP."

    chrom : string, default is "#CHROM", optional
        A string denoting the column name for chromosome. Defaults to be PLINK2.x's "#CHROM".
        Said column must be a character.

    pos : string, default is "POS", optional.
        A string denoting the column name for chromosomal position. Default to PLINK2.x's "POS".
        Said column must be numeric.

    pv : string, default is "P", optional.
        A string denoting the column name for chromosomal p-value. Default to PLINK2.x's "P".
        Said column must be float type.

    snp : string, default is "ID", optional.
        A string denoting the column name for the SNP name (rs number) or the column which you want to
        represent the variants. Default to PLINK2.x's "P". Said column should be a character.

    logp : bool, optional
        If TRUE, the -log10 of the p-value is plotted. It isn't very useful
        to plot raw p-values, but plotting the raw value could be useful for
        other genome-wide plots, for example, peak heights, bayes factors, test
        statistics, other "scores," etc. default: True

    ax : matplotlib axis, optional
        Axis to plot on, otherwise uses current axis.

    marker : matplotlib markers for scatter plot, default is "o", optional

    color : matplotlib color, optional, defauft: "#3B5488,#53BBD5"
        Color used for the plot elements. Could hex-code or rgb,
        e.g: '#3B5488,#53BBD5' or 'rb'

    alpha : float scalar, default is 0.8, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    title : string, or None, optional
        Set the title of the current plot.

    xlabel: string, optional
        Set the x axis label of the current axis.

    ylabel: string, optional
        Set the y axis label of the current axis.

    xtick_label_set : a set. optional
        Set the current x axis ticks of the current axis.

    CHR : string, or None, optional
        Select a specific chromosome to plot. And the x-axis will be the
        position of this chromosome instead of the chromosome id.

        CAUTION: this parameter could not be used with ``xtick_label_set``
                 together.

    xticklabel_kws : key, value pairings, or None, optional
        Other keyword arguments are passed to set xtick labels in
        maplotlib.axis.Axes.set_xticklabels.

    suggestiveline : float or None, default is 1e-5, optional
        Where to draw a suggestive ax.axhline. Set None to be disable.

    genomewideline : float or None, default is 5e-8
        Where to draw a genome-wide significant ax.axhline. Set None to be disable.

    sign_line_cols : matplotlib color, default: "#D62728,#2CA02C", optional.
        Color used for ``suggestiveline`` and ``genomewideline``.
        Could be hex-code or rgb, e.g: "#D62728,#2CA02C" or 'rb'

    hline_kws : key, value pairings, or None, optional
        keyword arguments for plotting ax.axhline(``suggestiveline`` and ``genomewideline``)
        except the "color" key-pair.

    sign_marker_p : float or None, default None, optional.
        A P-value threshold (suggestive to be 1e-6) for marking the significant SNP sites.

    sign_marker_color : matplotlib color, default: "r", optional.
        Define a color code for significant SNP sites.

    is_annotate_topsnp : boolean, default is False, optional.
        Annotate the top SNP or not for the significant locus.

    text_kws: key, value pairings, or None, optional
        keyword arguments for plotting in`` matplotlib.axes.Axes.text(x, y, s, fontdict=None, **kwargs)``

    ld_block_size : integer, default is 50000, optional
        Set the size of LD block which for finding top SNP. And the top SNP's annotation represent the block.

    kwargs : key, value pairings, optional
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
       it could also be used for any input data which have [chromo-
       some, position and p-value] dataframe.

    2. The right and top spines of the plot have been set to be
       invisible by hand.

    Examples
    --------

    Plot a basic manhattan plot from PLINK2.x association output and reture the figure:

    .. plot::
        :context: close-figs

        >>> import pandas as pd
        >>> from geneview import manhattanplot
        >>> from geneview.utils import load_dataset
        >>> df = load_dataset("gwas")
        >>> ax = manhattanplot(data=df)

    Plot a basic manhattan plot with horizontal xtick labels:

    .. plot::
        :context: close-figs

        >>> xtick = set(['chr' + i for i in list(map(str, range(1, 10))) + ['11', '13', '15', '18', '21', 'X']])
        >>> ax = manhattanplot(data=df, xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$", xtick_label_set=xtick)

    Add a horizontal at y position=3 line with linestyle="--" and lingwidth=1.3
    across the axis:

    .. plot::
        :context: close-figs
    
        >>> ax = manhattanplot(data=df, hline_kws={"linestyle": "--", "lw": 1.3}, xlabel="Chromosome",
        ...                    ylabel=r"$-log_{10}{(P)}$",
        ...                    xtick_label_set = xtick)

    Rotate the x-axis ticklabel by setting ``xticklabel_kws``:

    .. plot::
        :context: close-figs

        >>> ax = manhattanplot(data=df,
        ...                    hline_kws={"linestyle": "--", "lw": 1.3},
        ...                    xlabel="Chromosome",
        ...                    ylabel=r"$-log_{10}{(P)}$",
        ...                    xticklabel_kws={"rotation": "vertical"})

    Plot a better one with genome-wide significant mark and annotate the Top SNP and save
    the figure to "output_manhattan_plot.png":

    .. plot::
        :context: close-figs

        >>> from matplotlib.pyplot import subplots
        >>> fig, ax = subplots(figsize=(12, 4), facecolor="w", edgecolor="k")  # define a plot
        >>> ax = manhattanplot(data=df,
        ...                    marker=".",
        ...                    sign_marker_p=1e-6,  # Genome wide significant p-value
        ...                    sign_marker_color="r",
        ...                    snp="ID",
        ...                    title="Test",
        ...                    xtick_label_set=xtick,
        ...                    xlabel="Chromosome",
        ...                    ylabel=r"$-log_{10}{(P)}$",
        ...                    sign_line_cols=["#D62728", "#2CA02C"],
        ...                    hline_kws={"linestyle": "--", "lw": 1.3},
        ...                    is_annotate_topsnp=True,
        ...                    ld_block_size=50000,  # 50000 bp
        ...                    text_kws={"fontsize": 12,  # The fontsize of text
        ...                              "arrowprops": dict(arrowstyle="-", color="k", alpha=0.6)},
        ...                    ax=ax)
    """
    if not isinstance(data, DataFrame):
        raise ValueError("[ERROR] Input data must be a pandas.DataFrame.")
    if chrom not in data:
        raise ValueError("[ERROR] Column \"%s\" not found!" % chrom)
    if pos not in data:
        raise ValueError("[ERROR] Column \"%s\" not found!" % pos)
    if pv not in data:
        raise ValueError("[ERROR] Column \"%s\" not found!" % pv)
    if is_annotate_topsnp and (snp not in data):
        raise ValueError("[ERROR] You're trying to annotate a set of SNPs but "
                         "NO SNP \"%s\" column found!" % snp)
    if CHR is not None and xtick_label_set is not None:
        raise ValueError("[ERROR] ``CHR`` and ``xtick_label_set`` can't be set simultaneously.")

    data[[chrom]] = data[[chrom]].astype(str)  # make sure all the chromosome id are character.

    # Draw the plot and return the Axes
    if ax is None:
        # ax = plt.gca()
        _, ax = subplots(figsize=(9, 3), facecolor="w", edgecolor="k")  # default

    if xticklabel_kws is None:
        xticklabel_kws = {}
    if hline_kws is None:
        hline_kws = {}
    if text_kws is None:
        text_kws = {}

    if "," in color:
        color = color.split(",")
    colors = cycle(color)

    last_xpos = 0
    xs_by_id = []  # use for collecting chromosome's position on x-axis
    x, y, c = [], [], []
    sign_snp_sites = []
    for seqid, group_data in data.groupby(by=chrom, sort=False):  # keep the raw order of chromosome

        if (CHR is not None) and (seqid != CHR):
            continue

        color = next(colors)
        for i, (site, p_value) in enumerate(zip(group_data[pos], group_data[pv])):
            y_value = -np.log10(p_value) if logp else p_value

            x.append(last_xpos + site)
            y.append(y_value)

            c.append(sign_marker_color if ((sign_marker_p is not None) and (p_value <= sign_marker_p)) else color)
            if (sign_marker_p is not None) and (p_value <= sign_marker_p):
                snp_id = group_data[snp].iloc[i]
                sign_snp_sites.append([last_xpos + site, y_value, snp_id])  # x_pos, y_value, text

        # ``xs_by_id`` is for setting up positions and ticks. Ticks should
        # be placed in the middle of a chromosome. The a new pos column is 
        # added that keeps a running sum of the positions of each successive 
        # chromsome.
        xs_by_id.append([seqid, last_xpos + (group_data[pos].iloc[0] + group_data[pos].iloc[-1]) / 2])
        last_xpos = x[-1]  # keep track so that chromosome will not overlap in the plot.

    if not x:
        raise ValueError("zero-size array to reduction operation minimum which has no "
                         "identity. This could be caused by zero-size array of ``x`` "
                         "in the ``manhattanplot(...)`` function.")

    if "marker" not in kwargs:
        kwargs["marker"] = marker

    # plot the main manhattan dot plot
    ax.scatter(x, y, c=c, alpha=alpha, edgecolors="none", **kwargs)

    if is_annotate_topsnp is not None:
        index = _find_SNPs_which_overlap_sign_neighbour_region(
            sign_snp_neighbour_region=_sign_snp_regions(sign_snp_sites, ld_block_size),
            x=x)

        # reset color for all SNPs which nearby the top SNPs.
        for i in index:
            ax.scatter(x[i], y[i], c=sign_marker_color, alpha=alpha, edgecolors="none", **kwargs)

    # Add GWAS significant lines
    if "color" in hline_kws:
        hline_kws.pop("color")

    sign_line_cols = sign_line_cols.split(",") if "," in sign_line_cols else sign_line_cols
    if suggestiveline is not None:
        ax.axhline(y=-np.log10(suggestiveline) if logp else suggestiveline, color=sign_line_cols[0], **hline_kws)
    if genomewideline is not None:
        ax.axhline(y=-np.log10(genomewideline) if logp else genomewideline, color=sign_line_cols[1], **hline_kws)

    # Plotting the Top SNP for each significant block
    if is_annotate_topsnp:
        sign_top_snp = _find_top_snp(sign_snp_sites, ld_block_size=ld_block_size, is_get_biggest=logp)
        if sign_top_snp:  # not empty
            texts = [ax.text(_x, _y, _text) for _x, _y, _text in sign_top_snp]
            adjust_text(texts, ax=ax, **text_kws)

    if CHR is None:

        if xtick_label_set is not None:
            ax.set_xticks([v for c, v in xs_by_id if c in xtick_label_set])
            ax.set_xticklabels([c for c, v in xs_by_id if c in xtick_label_set], **xticklabel_kws)
        else:
            ax.set_xticks([v for c, v in xs_by_id])
            ax.set_xticklabels([c for c, v in xs_by_id], **xticklabel_kws)

    else:
        # show the whole chromosomal position without scientific notation
        # if you are just interesting in this chromosome.
        ax.get_xaxis().get_major_formatter().set_scientific(False)

    ax.set_xlim(0, x[-1])
    ax.set_ylim(ymin=min(y), ymax=1.2 * max(y))

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return ax


def _find_top_snp(sign_snp_data, ld_block_size, is_get_biggest=True):
    """
    :param sign_snp_data:  A 2D array: [[xpos1, yvalue1, text1], [xpos2, yvalue2, text2], ...]
    """
    top_snp = []
    tmp_cube = []
    for i, (_x, _y, text) in enumerate(sign_snp_data):
        if i == 0:
            tmp_cube.append([_x, _y, text])
            continue

        if _x > tmp_cube[-1][0] + ld_block_size:
            # Sorted by y_value in increase/decrease order and only get the first value [0], which is the TopSNP.
            top_snp.append(sorted(tmp_cube, key=(lambda x: x[1]), reverse=is_get_biggest)[0])
            tmp_cube = []

        tmp_cube.append([_x, _y, text])

    if tmp_cube:  # deal the last one
        top_snp.append(sorted(tmp_cube, key=(lambda x: x[1]), reverse=True)[0])

    return top_snp


def _sign_snp_regions(sign_snp_data, ld_block_size):
    """Create region according to the coordinate of sign_snp_data."""
    regions = []
    for i, (_x, _y, _t) in enumerate(sign_snp_data):
        if i == 0:
            regions.append([_x - ld_block_size, _x])
            continue

        if _x > regions[-1][1] + ld_block_size:
            regions[-1][1] += ld_block_size
            regions.append([_x - ld_block_size, _x])
        else:
            regions[-1][1] = _x

    # The last
    if regions:
        regions[-1][1] += ld_block_size

    return regions


def _find_SNPs_which_overlap_sign_neighbour_region(sign_snp_neighbour_region, x):
    """
    """
    x_size = len(x)
    reg_size = len(sign_snp_neighbour_region)
    index = []
    tmp_index = 0
    for i in range(x_size):
        _x = x[i]

        is_overlap = False
        iter_index = range(tmp_index, reg_size)
        for j in iter_index:
            if _x > sign_snp_neighbour_region[j][1]: continue
            if _x < sign_snp_neighbour_region[j][0]: break

            tmp_index = j
            is_overlap = True
            break

        if is_overlap:
            index.append(i)

    # return the index
    return index
