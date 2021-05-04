"""Plot admixture

Author: Shujia Huang
Date: 2021-05-01
"""
import warnings
from itertools import cycle
import numpy as np
import pandas as pd

from matplotlib.pyplot import subplots

from ..algorithm import hierarchical_cluster
from ..palette import generate_colors_palette


def _draw_admixtureplot(
        data=None,
        group_order=None,
        palette="tab10",
        xticklabels=None,
        xticklabel_kws=None,
        ylabel=None,
        ylabel_kws=None,
        hierarchical_kws=None,
        ax=None
):
    """Plot admixture.

    Parameters
    ----------
    data : dict.
        Input data structure. The key of `data` is the name of sub population and the
        value of `data` is a pandas.DataFrame which contain the admixture output (.Q)
        of the specify population. e.g:
        data = {"sub-pop1": pandas.DataFrame("admixture result(.Q) for 'sub-pop1'"),
                "sub-pop2": pandas.DataFrame("admixture result(.Q) for 'sub-pop2'")}

    group_order : vector of strings, optional
        Specify the order of processing and plotting for the estimating sub populations.
    """
    if ax is None:
        _, ax = subplots(1, 1, figsize=(14, 2), facecolor="w", constrained_layout=True, dpi=300)

    if hierarchical_kws is None:
        hierarchical_kws = {"method": "average", "metric": "euclidean"}
    if "axis" not in hierarchical_kws:
        hierarchical_kws["axis"] = 0  # row axis (by sample) to use to calculate cluster.

    if group_order is None:
        group_order = list(set(data.keys()))

    n = 0
    for g in group_order:
        n += len(data[g])
        hc = hierarchical_cluster(data=data[g], **hierarchical_kws)
        data[g] = hc.data.iloc[hc.reordered_index]  # re-order the data.

    x = np.arange(n)
    base_y = np.zeros(len(x))

    column_names = data[group_order[0]].columns
    if isinstance(palette, list) and (len(palette) < len(column_names)):
        msg = ("The categories of colors setting by `palette` is less than "
               "the number of estimating sub populations (K) in admixture, "
               "which could cause a confuse plot. Please reset the palette.")
        warnings.warn(msg)

    colors = generate_colors_palette(cmap=palette, n_colors=len(column_names))
    palette = cycle(colors)

    for k in column_names:
        c = next(palette)  # one color for one 'k'

        start_g_pos = 0
        add_y = []
        for g in group_order:  # keep group order
            try:
                y = data[g][k]
                end_g_pos = start_g_pos + len(y)
                ax.bar(x[start_g_pos:end_g_pos], y, bottom=base_y[start_g_pos:end_g_pos],
                       color=c, width=1.0, linewidth=0)

                start_g_pos = end_g_pos
                add_y.append(y)
            except KeyError:
                raise KeyError("KeyError: '%s'. Missing in 'group_order'." % g)

        base_y += np.array(pd.concat(add_y, axis=0))

    g_pos = 0
    xticks_pos = []
    for g in group_order:  # make group order
        g_size = len(data[g])
        xticks_pos.append(g_pos + 0.5 * g_size)

        g_pos += g_size
        ax.axvline(x=g_pos, color="k", lw=1)

    ax.spines["left"].set_linewidth(1)
    ax.spines["top"].set_linewidth(1)
    ax.spines["right"].set_linewidth(1)
    ax.spines["bottom"].set_linewidth(1)

    ax.set_xlim([x[0], x[-1]])
    ax.set_ylim([0, 1.0])
    ax.tick_params(bottom=False, top=False, left=False, right=False)

    if xticklabel_kws is None:
        xticklabel_kws = {}

    if xticklabels and (len(group_order) != len(xticklabels)):
        raise ValueError("The size of 'xticklabels'(%d) is not the same with "
                         "'group_order'(%d)." % (len(xticklabels), len(group_order)))

    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(group_order if xticklabels is None else xticklabels,
                       **xticklabel_kws)

    if ylabel_kws is None:
        ylabel_kws = {}

    ax.set_ylabel(ylabel, **ylabel_kws)
    ax.set_yticks([])

    return ax


def _load_admixture_from_file(in_admixture_fname, in_sample_info_fname, shuffle_popsample_kws=None):

    if ("axis" in shuffle_popsample_kws) and (shuffle_popsample_kws["axis"] == 1):
        warnings.warn("axis=1 means sampling the data by columns, Which is not "
                      "allow and may not be right in admixture data.")

    df = pd.read_table(in_admixture_fname, sep=" ", header=None)
    sample_info = pd.read_table(in_sample_info_fname, header=None, names=["Group"])
    popset = set(sample_info["Group"])

    if len(sample_info) != len(df):
        raise ValueError("The size of sample_info(%d) and input admixture result(%d) "
                         "must be the same." % (len(sample_info), len(df)))

    data = {}
    for g in popset:
        g_data = df[sample_info["Group"] == g].copy()
        g_size = len(g_data)
        if shuffle_popsample_kws:
            shuffle_raw_n = shuffle_popsample_kws["n"] if "n" in shuffle_popsample_kws else None
            if shuffle_raw_n and shuffle_raw_n > g_size:
                shuffle_popsample_kws["n"] = g_size
                data[g] = g_data.sample(**shuffle_popsample_kws)
                shuffle_popsample_kws["n"] = shuffle_raw_n  # reset to the raw N
            else:
                data[g] = g_data.sample(**shuffle_popsample_kws)
        else:
            data[g] = g_data

    return data


def admixtureplot(
        data,
        population_info=None,
        shuffle_popsample_kws=None,
        group_order=None,
        palette="tab10",
        xticklabels=None,
        xticklabel_kws=None,
        ylabel=None,
        ylabel_kws=None,
        hierarchical_kws=None,
        ax=None
):
    """Plot admixture.

        Parameters
        ----------
        data : a file path to the admixture output(.Q) or a dict.
            Input data structure.
            The key of `data` is the name of sub population and the value of `data`
            is a pandas.DataFrame which contain the admixture output (.Q) of the
            specify population if `data` is a dict. e.g:
            data = {"sub-pop1": pandas.DataFrame("admixture result(.Q) for 'sub-pop1'"),
                    "sub-pop2": pandas.DataFrame("admixture result(.Q) for 'sub-pop2'")}

        population_info : String or None, optional
            File path to a sample information list, each line only contain a group information
            for admixture row.

        shuffle_popsample_kws : key, value pairings, or None, optional
            The keyword argument of `pandas.sample` method to sample rows (in random order) for
            each group of population. If None, keep the raw data.

        group_order : vector of strings, optional
            Specify the order of processing and plotting for the estimating sub populations.

        palette : string, list, or :class:`matplotlib.colors.Colormap`, optional, default: tab10.
            String values are passed to :func:`generate_colors`. List values imply categorical
            mapping, while a colormap object implies numeric mapping.

        xticklabels : list, optional
            The label texts of x-axis.

        xticklabel_kws : key, value pairings, or None, optional
            Other keyword arguments are passed to set xtick labels in
            maplotlib.axis.Axes.set_xticklabels.

        ylabel : string, optional
            Set the y axis label of the current axis. The label will set to be the
            column number of admixture output if ylabel is None.

        ylabel_kws : key, value pairings, or None, optional
            Other keyword arguments are passed to set y label in
            maplotlib.axis.Axes.set_ylabel.

        Returns
        -------
        ax : matplotlib Axes
            Axes object with the admixtureplot.

        Examples
        --------

        A simple admixture plot example.

        .. plot::
            :context: close-figs

            >>> from geneview import admixtureplot
            >>> ax = admixtureplot(data="../../examples/data/admixture_1KG.output.Q",
            ...                    population_info="../../examples/data/admixture_1KG_population.info")

        Plot the admixture by a specify population order

        .. plot::
            :context: close-figs

            >>> pop_group_1kg = ["KHV", "CDX", "CHS", "CHB", "JPT", "BEB", "STU", "ITU",
            ...                  "GIH", "PJL", "FIN", "CEU", "GBR", "IBS", "TSI", "PEL",
            ...                  "PUR", "MXL", "CLM", "ASW", "ACB", "GWD", "MSL", "YRI",
            ...                  "ESN", "LWK"]
            >>> ax = admixtureplot(data="../../examples/data/admixture.output.Q",
            ...                    population_info="../../examples/data/admixture_population.info")


        Setting the ``shuffle_popsample_kws`` argument if you wish to sample some samples for each
        population group in the admixture plot.

        only sampling 80 samples for each population group.
        .. plot::
            :context: close-figs
            >>> ax = admixtureplot(data="../../examples/data/admixture.output.Q",
            ...                    population_info="../../examples/data/admixture_1KG_population.info",
            ...                    shuffle_popsample_kws={"n": 80},
            ...                    group_order=pop_group_1kg)

        Or specifies the fraction of each group population.

        .. plot::
            :context: close-figs
            >>> ax = admixtureplot(data="../../examples/data/admixture.output.Q",
            ...                    population_info="../../examples/data/admixture_1KG_population.info",
            ...                    shuffle_popsample_kws={"frac": 0.2},   # 0.2 means 20% of the data (in random order).
            ...                    group_order=pop_group_1kg)


        Admixtureplot also allow you to define your own data (in dict type) by yourself.


        .. plot::
            :context: close-figs

            >>> df = pd.read_table("../../examples/data/admixture.output.Q", sep=" ", header=None)
            >>> sample_info = pd.read_table("../../examples/data/admixture_population.info", sep="\t", header=None, names=["Group"])
            >>> popset = set(sample_info["Group"])
            >>> data = {}
            >>> for g in popset:
            ...     g_data = df[sample_info["Group"]==g].copy()
            ...     # Sub sampling: keep less than 140 samples for each group
            ...     data[g] = g_data.sample(n=140, random_state=100) if len(g_data)>140 else g_data

            # Plot the figure
            >>> f, ax = plt.subplots(1, 1, figsize=(14, 3), facecolor="w", constrained_layout=True, dpi=300)
            >>> ax = admixtureplot(data=data,
            ...                    group_order=pop_group_1kg,
            ...                    palette="Set1",
            ...                    xticklabel_kws={"rotation": "vertical"},
            ...                    ylabel="K=11",
            ...                    ylabel_kws={"rotation": 0, "ha": "right"},
            ...                    ax=ax)
        """
    if shuffle_popsample_kws is None:
        shuffle_popsample_kws = {}

    if (population_info is not None) and (not isinstance(population_info, str)):
        raise ValueError("`population_info` must be a file path which contain "
                         "the population information for each row of admixture "
                         "result. One element per line.")

    if (population_info is None) and (not isinstance(data, dict)):
        raise ValueError("``data`` must be a dictionary and the key should be the group "
                         "information, values should be a dataframe which contain the result "
                         "of admixture.")

    if (not isinstance(data, dict)) and (not isinstance(data, str)):
        raise ValueError("`data` should be a dict or a file path to the admixture output(.Q).")

    if isinstance(data, str):
        data = _load_admixture_from_file(
            data, population_info, shuffle_popsample_kws=shuffle_popsample_kws
        )

    # infer ylabel if ylabel is None
    if ylabel is None:
        g = list(data.keys())[0]
        ylabel = "K=%d" % len(data[g].columns)

    return _draw_admixtureplot(data=data,
                               group_order=group_order,
                               palette=palette,
                               xticklabels=xticklabels,
                               xticklabel_kws=xticklabel_kws,
                               ylabel=ylabel,
                               ylabel_kws=ylabel_kws,
                               hierarchical_kws=hierarchical_kws,
                               ax=ax)
