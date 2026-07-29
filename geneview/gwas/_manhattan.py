"""Plotting functions for manhattan plot.

Copyright (c) Shujia Huang
Date: 2021-02-21

This model is based on brentp's script on github:
https://github.com/brentp/bio-playground/blob/master/plots/manhattan-plot.py

Thanks for Brentp's contributions

"""
import string
from pandas import DataFrame
import numpy as np

from ..utils import adjust_text
from ..utils._adjust_text import get_renderer
from .._core import styled_plot, color_cycle


# learn something from "https://github.com/reneshbedre/bioinfokit/blob/38fb4966827337f00421119a69259b92bb67a7d0/bioinfokit/visuz.py"
@styled_plot(figsize=(9, 3), subplot_kws={"facecolor": "w", "edgecolor": "k"})
def manhattanplot(data, chrom="#CHROM", pos="POS", pv="P", snp="ID", logp=True, ax=None,
                  marker=".", color="#3B5488,#53BBD5", alpha=0.8,
                  title=None, xlabel="Chromosome", ylabel=r"$-log_{10}{(P)}$",
                  xtick_label_set=None, CHR=None, xticklabel_kws=None,
                  suggestiveline=1e-5, genomewideline=5e-8, sign_line_cols="#D62728,#2CA02C", hline_kws=None,
                  sign_marker_p=None, sign_marker_color="r",
                  is_annotate_topsnp=False, text_kws=None, ld_block_size=50000,
                  annotate_fmt=None, annotate_layout="repel", adjust_text_kws=None,
                  style=None, **kwargs):
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
        keyword arguments forwarded to ``matplotlib.axes.Axes.text`` for styling
        the annotation labels (e.g. ``fontsize``, ``color``, ``rotation``). An
        ``arrowprops`` entry, if present, is used to draw the connecting arrows
        rather than styling the text.

    ld_block_size : integer, default is 50000, optional
        Set the size of LD block which for finding top SNP. And the top SNP's annotation represent the block.

    annotate_fmt : str, callable, or None, default None, optional
        Controls the content of each top-SNP label. ``None`` shows the SNP id
        only. A format string may reference the fields ``snp``, ``chrom``,
        ``pos``, ``p`` and ``log10p`` (e.g. ``"{snp}\n{p:.1e}"``). A callable
        receives those fields as keyword arguments and returns the label text.

    annotate_layout : str, default "repel", optional
        Strategy used to place the top-SNP labels:

        - ``"repel"``: iteratively push labels apart with ``adjust_text`` (good
          for a handful of labels).
        - ``"lane"``: lay the labels out in a single non-overlapping row near
          the top of the axes with leader lines back to each point. This is
          O(N log N), stays tidy and is much faster when there are many
          significant loci.

    adjust_text_kws : dict or None, optional
        Extra keyword arguments forwarded to ``adjust_text`` (only used by the
        ``"repel"`` layout), e.g. ``force_text``, ``expand_text``, ``only_move``
        or ``lim``. These override the sensible defaults chosen by geneview.

    style : str, PlotStyle, or None, optional
        Plot style to apply. Can be a registered style name (e.g. "nature",
        "science", "cell"), a PlotStyle object, or None (the default) to use
        the currently active style.

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
    _validate_annotate_fmt(annotate_fmt)  # fail fast on bad label format strings

    # ``ax`` is guaranteed non-None and the style context is already active:
    # both are handled by the ``@styled_plot`` decorator on this function.
    return _manhattanplot_impl(
        data, chrom, pos, pv, snp, logp, ax, marker, color, alpha,
        title, xlabel, ylabel, xtick_label_set, CHR, xticklabel_kws,
        suggestiveline, genomewideline, sign_line_cols, hline_kws,
        sign_marker_p, sign_marker_color, is_annotate_topsnp, text_kws,
        ld_block_size, annotate_fmt, annotate_layout, adjust_text_kws,
        **kwargs
    )


def _manhattanplot_impl(
    data, chrom, pos, pv, snp, logp, ax, marker, color, alpha,
    title, xlabel, ylabel, xtick_label_set, CHR, xticklabel_kws,
    suggestiveline, genomewideline, sign_line_cols, hline_kws,
    sign_marker_p, sign_marker_color, is_annotate_topsnp, text_kws,
    ld_block_size, annotate_fmt=None, annotate_layout="repel",
    adjust_text_kws=None, **kwargs
):
    """Internal implementation of manhattanplot, called within a style context."""
    data[[chrom]] = data[[chrom]].astype(str)  # make sure all the chromosome id are character.

    if xticklabel_kws is None:
        xticklabel_kws = {}
    if hline_kws is None:
        hline_kws = {}
    if text_kws is None:
        text_kws = {}

    colors = color_cycle(color)

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

            if (snp is not None) and (sign_marker_p is not None) and (p_value <= sign_marker_p):
                snp_id = group_data[snp].iloc[i]
                # x_pos, y_value, text, chrom, p_value, position
                sign_snp_sites.append([last_xpos + site, y_value, snp_id, seqid, p_value, site])

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

    if is_annotate_topsnp:
        index = _find_SNPs_which_overlap_sign_neighbour_region(
            sign_snp_neighbour_region=_sign_snp_regions(sign_snp_sites, ld_block_size),
            x=x)

        # reset color for all SNPs which nearby the top SNPs.
        xs, ys = [], []
        for i in index:
            xs.append(x[i])
            ys.append(y[i])

        ax.scatter(xs, ys, c=sign_marker_color, alpha=alpha, edgecolors="none", **kwargs)

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
    else:
        sign_top_snp = None

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

    # Annotate the top SNPs only after the axes limits are final; the label
    # placement works in display coordinates and needs the final dimensions.
    if is_annotate_topsnp and sign_top_snp:
        _annotate_top_snps(ax, sign_top_snp, logp=logp, annotate_fmt=annotate_fmt,
                           layout=annotate_layout, text_kws=text_kws,
                           adjust_text_kws=adjust_text_kws)

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    # Top/right spine visibility is enforced by the active style via the
    # ``@styled_plot`` decorator (the geneview default hides them).
    return ax


def _find_top_snp(sign_snp_data, ld_block_size, is_get_biggest=True):
    """Pick the representative (top) SNP for each LD block.

    :param sign_snp_data:  A 2D array where each row is at least
        ``[xpos, yvalue, text]`` and may carry extra trailing fields such as
        ``chrom, p_value, position``. The *whole* row is preserved in the
        returned records so downstream annotation can format richer labels.
    """
    top_snp = []
    tmp_cube = []
    current_chrom = None
    for i, item in enumerate(sign_snp_data):
        _x = item[0]
        _chrom = item[3] if len(item) > 3 else None

        if i == 0:
            tmp_cube.append(item)
            current_chrom = _chrom
            continue

        if _chrom != current_chrom or _x > tmp_cube[-1][0] + ld_block_size:
            # Sorted by y_value in increase/decrease order and only get the first value [0], which is the TopSNP.
            top_snp.append(sorted(tmp_cube, key=(lambda v: v[1]), reverse=is_get_biggest)[0])
            tmp_cube = []
            current_chrom = _chrom

        tmp_cube.append(item)

    if tmp_cube:  # deal the last one
        top_snp.append(sorted(tmp_cube, key=(lambda v: v[1]), reverse=is_get_biggest)[0])

    return top_snp


def _sign_snp_regions(sign_snp_data, ld_block_size):
    """Create region according to the coordinate of sign_snp_data."""
    regions = []
    for i, item in enumerate(sign_snp_data):
        _x = item[0]
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


_MISSING = object()


#: Field names that may be referenced in an ``annotate_fmt`` format string.
_ANNOTATE_FMT_FIELDS = ("snp", "chrom", "pos", "p", "log10p")


def _validate_annotate_fmt(annotate_fmt):
    """Fail fast on an unusable ``annotate_fmt`` format string.

    Unknown field names would otherwise surface as a cryptic ``KeyError``
    deep inside the annotation pass, after the whole plot has been drawn;
    invalid format specs similarly blow up late as ``ValueError``.  This
    check turns both into a single clear ``ValueError`` at call time.
    ``None`` and callables need no validation.
    """
    if annotate_fmt is None or callable(annotate_fmt):
        return
    if not isinstance(annotate_fmt, str):
        raise ValueError(
            "[ERROR] ``annotate_fmt`` must be None, a format string, or a "
            "callable; got %s." % type(annotate_fmt).__name__)

    for _, field, _, _ in string.Formatter().parse(annotate_fmt):
        if field is None:
            continue
        root = field.split(".")[0].split("[")[0]
        if root not in _ANNOTATE_FMT_FIELDS:
            raise ValueError(
                "[ERROR] Unknown field %r in ``annotate_fmt`` %r. Available "
                "fields: %s."
                % (root, annotate_fmt, ", ".join("{%s}" % f for f in _ANNOTATE_FMT_FIELDS)))

    # Render once against dummy values so an invalid format spec
    # (e.g. ``{p:xx}``) is rejected here instead of mid-plot.
    dummy = dict(snp="rs1", chrom="1", pos=1000, p=0.001, log10p=3.0)
    try:
        annotate_fmt.format(**dummy)
    except (ValueError, IndexError, KeyError) as exc:
        raise ValueError(
            "[ERROR] Invalid ``annotate_fmt`` %r: %s. Available fields: %s."
            % (annotate_fmt, exc, ", ".join("{%s}" % f for f in _ANNOTATE_FMT_FIELDS)))


def _build_label(item, logp, annotate_fmt):
    """Turn a top-SNP record into its annotation string.

    ``item`` is ``[x, y, snp, chrom, p, pos]`` with trailing fields optional.
    ``annotate_fmt`` may be None (snp id only), a format string referencing the
    fields ``snp/chrom/pos/p/log10p``, or a callable receiving them as kwargs.
    """
    snp = item[2] if len(item) > 2 else ""
    chrom = item[3] if len(item) > 3 else None
    p = item[4] if len(item) > 4 else None
    position = item[5] if len(item) > 5 else None

    if annotate_fmt is None:
        return str(snp)

    if p is not None and p > 0:
        log10p = -np.log10(p)
    else:
        log10p = item[1] if logp else None

    fields = dict(snp=snp, chrom=chrom, pos=position, p=p, log10p=log10p)
    if callable(annotate_fmt):
        return str(annotate_fmt(**fields))
    return annotate_fmt.format(**fields)


def _annotate_top_snps(ax, sign_top_snp, logp, annotate_fmt, layout, text_kws, adjust_text_kws):
    """Place the top-SNP labels using the requested layout strategy."""
    text_kws = dict(text_kws or {})
    adjust_text_kws = dict(adjust_text_kws or {})

    # ``arrowprops`` styles the connecting arrows, not the text; accept it from
    # either dict for convenience. An explicit ``None`` disables the arrows.
    arrowprops = text_kws.pop("arrowprops", _MISSING)
    if arrowprops is _MISSING:
        arrowprops = adjust_text_kws.pop("arrowprops", _MISSING)
    if arrowprops is _MISSING:
        arrowprops = dict(arrowstyle="-", color="0.5", lw=0.6, alpha=0.7)

    x_pos = [item[0] for item in sign_top_snp]
    y_pos = [item[1] for item in sign_top_snp]
    labels = [_build_label(item, logp, annotate_fmt) for item in sign_top_snp]

    if layout == "lane":
        return _layout_top_lane(ax, x_pos, y_pos, labels, text_kws, arrowprops)

    if layout != "repel":
        raise ValueError("[ERROR] ``annotate_layout`` must be 'repel' or 'lane', "
                         "got %r." % layout)

    texts = [ax.text(xx, yy, s, **text_kws) for xx, yy, s in zip(x_pos, y_pos, labels)]
    kws = dict(
        only_move={"points": "y", "text": "xy", "objects": "xy"},
        force_text=(0.3, 0.5),
        expand_text=(1.05, 1.4),
    )
    kws.update(adjust_text_kws)
    if arrowprops is not None:
        kws["arrowprops"] = arrowprops
    adjust_text(texts, ax=ax, **kws)
    return texts


def _layout_top_lane(ax, x_pos, y_pos, labels, text_kws, arrowprops):
    """Lay labels out in a single non-overlapping row near the top of the axes.

    Labels are sorted by x, then a one-dimensional sweep pushes them apart just
    enough to remove horizontal overlaps (O(N log N)). A leader line links each
    label back to its SNP. This stays tidy and fast even with many loci.
    """
    if not x_pos:
        return []

    order = sorted(range(len(x_pos)), key=lambda k: x_pos[k])
    rotation = text_kws.pop("rotation", 90)
    kws = dict(text_kws)
    kws.setdefault("va", "top")
    kws.setdefault("ha", "center")
    kws.setdefault("clip_on", False)

    ymin, ymax = ax.get_ylim()
    y_label = ymin + (ymax - ymin) * 0.98

    texts = [ax.text(x_pos[k], y_label, labels[k], rotation=rotation, **kws) for k in order]

    fig = ax.get_figure()
    fig.canvas.draw()
    r = get_renderer(fig)
    inv = ax.transData.inverted()
    widths = []
    for t in texts:
        bb = t.get_window_extent(r)
        (dx0, _), (dx1, _) = inv.transform([(bb.x0, 0.0), (bb.x1, 0.0)])
        widths.append(abs(dx1 - dx0))

    # Forward sweep: guarantee a minimum centre-to-centre gap between neighbours.
    new_xs = [x_pos[k] for k in order]
    for i in range(1, len(new_xs)):
        min_x = new_xs[i - 1] + (widths[i - 1] + widths[i]) / 2.0 * 1.15
        if new_xs[i] < min_x:
            new_xs[i] = min_x

    # If we ran past the right edge, clamp and sweep back to keep gaps intact.
    xlo, xhi = ax.get_xlim()
    if new_xs[-1] > xhi:
        new_xs[-1] = xhi
        for i in range(len(new_xs) - 2, -1, -1):
            max_x = new_xs[i + 1] - (widths[i] + widths[i + 1]) / 2.0 * 1.15
            if new_xs[i] > max_x:
                new_xs[i] = max_x

    for idx, k in enumerate(order):
        texts[idx].set_position((new_xs[idx], y_label))
        if arrowprops is not None:
            ax.annotate("", xy=(x_pos[k], y_pos[k]), xytext=(new_xs[idx], y_label),
                        arrowprops=arrowprops, annotation_clip=False)

    return texts
