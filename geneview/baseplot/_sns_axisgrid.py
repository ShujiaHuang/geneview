from __future__ import division

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from ..palette import color_palette
from .. import utils


class JointGrid(object):
    """Grid for drawing a bivariate plot with marginal univariate plots."""

    def __init__(self, x, y, data=None, size=6, ratio=5, space=.2,
                 dropna=True, xlim=None, ylim=None):
        """Set up the grid of subplots.

        Parameters
        ----------
        x, y : strings or vectors
            Data or names of variables in ``data``.
        data : DataFrame, optional
            DataFrame when ``x`` and ``y`` are variable names.
        size : numeric
            Size of each side of the figure in inches (it will be square).
        ratio : numeric
            Ratio of joint axes size to marginal axes height.
        space : numeric, optional
            Space between the joint and marginal axes
        dropna : bool, optional
            If True, remove observations that are missing from `x` and `y`.
        {x, y}lim : two-tuples, optional
            Axis limits to set before plotting.

        See Also
        --------
        jointplot : High-level interface for drawing bivariate plots with
                    several different default plot kinds.

        Examples
        --------

        Initialize the figure but don't draw any plots onto it:

        .. plot::
            :context: close-figs

            >>> import geneview as gv; gv.setup(style="ticks", color_codes=True)
            >>> tips = gv.load_dataset("tips")
            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)

        Add plots using default parameters:

        .. plot::
            :context: close-figs

            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot(gv.regplot, gv.distplot)

        Draw the join and marginal plots separately, which allows finer-level
        control other parameters:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter, color=".5", edgecolor="white")
            >>> g = g.plot_marginals(gv.distplot, kde=False, color=".5")

        Draw the two marginal plots separately:

        .. plot::
            :context: close-figs

            >>> import numpy as np
            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter, color="m", edgecolor="white")
            >>> _ = g.ax_marg_x.hist(tips["total_bill"], color="b", alpha=.6,
            ...                      bins=np.arange(0, 60, 5))
            >>> _ = g.ax_marg_y.hist(tips["tip"], color="r", alpha=.6,
            ...                      orientation="horizontal",
            ...                      bins=np.arange(0, 12, 1))

        Add an annotation with a statistic summarizing the bivariate
        relationship:

        .. plot::
            :context: close-figs

            >>> from scipy import stats
            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter,
            ...                  color="g", s=40, edgecolor="white")
            >>> g = g.plot_marginals(gv.distplot, kde=False, color="g")
            >>> g = g.annotate(stats.pearsonr)

        Use a custom function and formatting for the annotation

        .. plot::
            :context: close-figs

            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips)
            >>> g = g.plot_joint(plt.scatter,
            ...                  color="g", s=40, edgecolor="white")
            >>> g = g.plot_marginals(gv.distplot, kde=False, color="g")
            >>> rsquare = lambda a, b: stats.pearsonr(a, b)[0] ** 2
            >>> g = g.annotate(rsquare, template="{stat}: {val:.2f}",
            ...                stat="$R^2$", loc="upper left", fontsize=12)

        Remove the space between the joint and marginal axes:

        .. plot::
            :context: close-figs

            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips, space=0)
            >>> g = g.plot_joint(gv.kdeplot, cmap="Blues_d")
            >>> g = g.plot_marginals(gv.kdeplot, shade=True)

        Draw a smaller plot with relatively larger marginal axes:

        .. plot::
            :context: close-figs

            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips,
            ...                   size=5, ratio=2)
            >>> g = g.plot_joint(gv.kdeplot, cmap="Reds_d")
            >>> g = g.plot_marginals(gv.kdeplot, color="r", shade=True)

        Set limits on the axes:

        .. plot::
            :context: close-figs

            >>> g = gv.JointGrid(x="total_bill", y="tip", data=tips,
            ...                   xlim=(0, 50), ylim=(0, 8))
            >>> g = g.plot_joint(gv.kdeplot, cmap="Purples_d")
            >>> g = g.plot_marginals(gv.kdeplot, color="m", shade=True)

        """
        # Set up the subplot grid
        f = plt.figure(figsize=(size, size))
        gs = plt.GridSpec(ratio + 1, ratio + 1)

        ax_joint = f.add_subplot(gs[1:, :-1])
        ax_marg_x = f.add_subplot(gs[0, :-1], sharex=ax_joint)
        ax_marg_y = f.add_subplot(gs[1:, -1], sharey=ax_joint)

        self.fig = f
        self.ax_joint = ax_joint
        self.ax_marg_x = ax_marg_x
        self.ax_marg_y = ax_marg_y

        # Turn off tick visibility for the measure axis on the marginal plots
        plt.setp(ax_marg_x.get_xticklabels(), visible=False)
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)

        # Turn off the ticks on the density axis for the marginal plots
        plt.setp(ax_marg_x.yaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_x.yaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_majorticklines(), visible=False)
        plt.setp(ax_marg_y.xaxis.get_minorticklines(), visible=False)
        plt.setp(ax_marg_x.get_yticklabels(), visible=False)
        plt.setp(ax_marg_y.get_xticklabels(), visible=False)
        ax_marg_x.yaxis.grid(False)
        ax_marg_y.xaxis.grid(False)

        # Possibly extract the variables from a DataFrame
        if data is not None:
            if x in data:
                x = data[x]
            if y in data:
                y = data[y]

        # Possibly drop NA
        if dropna:
            not_na = pd.notnull(x) & pd.notnull(y)
            x = x[not_na]
            y = y[not_na]

        # Find the names of the variables
        if hasattr(x, "name"):
            xlabel = x.name
            ax_joint.set_xlabel(xlabel)
        if hasattr(y, "name"):
            ylabel = y.name
            ax_joint.set_ylabel(ylabel)

        # Convert the x and y data to arrays for plotting
        self.x = np.asarray(x)
        self.y = np.asarray(y)

        if xlim is not None:
            ax_joint.set_xlim(xlim)
        if ylim is not None:
            ax_joint.set_ylim(ylim)

        # Make the grid look nice
        utils.despine(f)
        utils.despine(ax=ax_marg_x, left=True)
        utils.despine(ax=ax_marg_y, bottom=True)
        f.tight_layout()
        f.subplots_adjust(hspace=space, wspace=space)

    def plot(self, joint_func, marginal_func, annot_func=None):
        """Shortcut to draw the full plot.

        Use `plot_joint` and `plot_marginals` directly for more control.

        Parameters
        ----------
        joint_func, marginal_func: callables
            Functions to draw the bivariate and univariate plots.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        self.plot_marginals(marginal_func)
        self.plot_joint(joint_func)
        if annot_func is not None:
            self.annotate(annot_func)
        return self

    def plot_joint(self, func, **kwargs):
        """Draw a bivariate plot of `x` and `y`.

        Parameters
        ----------
        func : plotting callable
            This must take two 1d arrays of data as the first two
            positional arguments, and it must plot on the "current" axes.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        plt.sca(self.ax_joint)
        func(self.x, self.y, **kwargs)

        return self

    def plot_marginals(self, func, **kwargs):
        """Draw univariate plots for `x` and `y` separately.

        Parameters
        ----------
        func : plotting callable
            This must take a 1d array of data as the first positional
            argument, it must plot on the "current" axes, and it must
            accept a "vertical" keyword argument to orient the measure
            dimension of the plot vertically.
        kwargs : key, value mappings
            Keyword argument are passed to the plotting function.

        Returns
        -------
        self : JointGrid instance
            Returns `self`.

        """
        kwargs["vertical"] = False
        plt.sca(self.ax_marg_x)
        func(self.x, **kwargs)

        kwargs["vertical"] = True
        plt.sca(self.ax_marg_y)
        func(self.y, **kwargs)

        return self

    def annotate(self, func, template=None, stat=None, loc="best", **kwargs):
        """Annotate the plot with a statistic about the relationship.

        Parameters
        ----------
        func : callable
            Statistical function that maps the x, y vectors either to (val, p)
            or to val.
        template : string format template, optional
            The template must have the format keys "stat" and "val";
            if `func` returns a p value, it should also have the key "p".
        stat : string, optional
            Name to use for the statistic in the annotation, by default it
            uses the name of `func`.
        loc : string or int, optional
            Matplotlib legend location code; used to place the annotation.
        kwargs : key, value mappings
            Other keyword arguments are passed to `ax.legend`, which formats
            the annotation.

        Returns
        -------
        self : JointGrid instance.
            Returns `self`.

        """
        default_template = "{stat} = {val:.2g}; p = {p:.2g}"

        # Call the function and determine the form of the return value(s)
        out = func(self.x, self.y)
        try:
            val, p = out
        except TypeError:
            val, p = out, None
            default_template, _ = default_template.split(";")

        # Set the default template
        if template is None:
            template = default_template

        # Default to name of the function
        if stat is None:
            stat = func.__name__

        # Format the annotation
        if p is None:
            annotation = template.format(stat=stat, val=val)
        else:
            annotation = template.format(stat=stat, val=val, p=p)

        # Draw an invisible plot and use the legend to draw the annotation
        # This is a bit of a hack, but `loc=best` works nicely and is not
        # easily abstracted.
        phantom, = self.ax_joint.plot(self.x, self.y, linestyle="", alpha=0)
        self.ax_joint.legend([phantom], [annotation], loc=loc, **kwargs)
        phantom.remove()

        return self

    def set_axis_labels(self, xlabel="", ylabel="", **kwargs):
        """Set the axis labels on the bivariate axes.

        Parameters
        ----------
        xlabel, ylabel : strings
            Label names for the x and y variables.
        kwargs : key, value mappings
            Other keyword arguments are passed to the set_xlabel or
            set_ylabel.

        Returns
        -------
        self : JointGrid instance
            returns `self`

        """
        self.ax_joint.set_xlabel(xlabel, **kwargs)
        self.ax_joint.set_ylabel(ylabel, **kwargs)
        return self

    def savefig(self, *args, **kwargs):
        """Wrap figure.savefig defaulting to tight bounding box."""
        kwargs.setdefault("bbox_inches", "tight")
        self.fig.savefig(*args, **kwargs)


class Grid(object):
    """Base class for grids of subplots."""
    _margin_titles = False
    _legend_out = True

    def set(self, **kwargs):
        """Set attributes on each subplot Axes."""
        for ax in self.axes.flat:
            ax.set(**kwargs)

        return self

    def savefig(self, *args, **kwargs):
        """Save the figure."""
        kwargs = kwargs.copy()
        kwargs.setdefault("bbox_inches", "tight")
        self.fig.savefig(*args, **kwargs)

    def add_legend(self, legend_data=None, title=None, label_order=None,
                   **kwargs):
        """Draw a legend, maybe placing it outside axes and resizing the figure.

        Parameters
        ----------
        legend_data : dict, optional
            Dictionary mapping label names to matplotlib artist handles. The
            default reads from ``self._legend_data``.
        title : string, optional
            Title for the legend. The default reads from ``self._hue_var``.
        label_order : list of labels, optional
            The order that the legend entries should appear in. The default
            reads from ``self.hue_names`` or sorts the keys in ``legend_data``.
        kwargs : key, value pairings
            Other keyword arguments are passed to the underlying legend methods
            on the Figure or Axes object.

        Returns
        -------
        self : Grid instance
            Returns self for easy chaining.

        """
        # Find the data for the legend
        legend_data = self._legend_data if legend_data is None else legend_data
        if label_order is None:
            if self.hue_names is None:
                label_order = np.sort(list(legend_data.keys()))
            else:
                label_order = list(map(str, self.hue_names))

        blank_handle = mpl.patches.Patch(alpha=0, linewidth=0)
        handles = [legend_data.get(l, blank_handle) for l in label_order]
        title = self._hue_var if title is None else title
        try:
            title_size = mpl.rcParams["axes.labelsize"] * .85
        except TypeError:  # labelsize is something like "large"
            title_size = mpl.rcParams["axes.labelsize"]

        # Set default legend kwargs
        kwargs.setdefault("scatterpoints", 1)

        if self._legend_out:
            # Draw a full-figure legend outside the grid
            figlegend = self.fig.legend(handles, label_order, "center right",
                                        **kwargs)
            self._legend = figlegend
            figlegend.set_title(title)

            # Set the title size a roundabout way to maintain
            # compatability with matplotlib 1.1
            prop = mpl.font_manager.FontProperties(size=title_size)
            figlegend._legend_title_box._text.set_font_properties(prop)

            # Draw the plot to set the bounding boxes correctly
            plt.draw()

            # Calculate and set the new width of the figure so the legend fits
            legend_width = figlegend.get_window_extent().width / self.fig.dpi
            figure_width = self.fig.get_figwidth()
            self.fig.set_figwidth(figure_width + legend_width)

            # Draw the plot again to get the new transformations
            plt.draw()

            # Now calculate how much space we need on the right side
            legend_width = figlegend.get_window_extent().width / self.fig.dpi
            space_needed = legend_width / (figure_width + legend_width)
            margin = .04 if self._margin_titles else .01
            self._space_needed = margin + space_needed
            right = 1 - self._space_needed

            # Place the subplot axes to give space for the legend
            self.fig.subplots_adjust(right=right)

        else:
            # Draw a legend in the first axis
            ax = self.axes.flat[0]
            leg = ax.legend(handles, label_order, loc="best", **kwargs)
            leg.set_title(title)

            # Set the title size a roundabout way to maintain
            # compatability with matplotlib 1.1
            prop = mpl.font_manager.FontProperties(size=title_size)
            leg._legend_title_box._text.set_font_properties(prop)

        return self

    def _clean_axis(self, ax):
        """Turn off axis labels and legend."""
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.legend_ = None
        return self

    def _update_legend_data(self, ax):
        """Extract the legend data from an axes object and save it."""
        handles, labels = ax.get_legend_handles_labels()
        data = {l: h for h, l in zip(handles, labels)}
        self._legend_data.update(data)

    def _get_palette(self, data, hue, hue_order, palette):
        """Get a list of colors for the hue variable."""
        if hue is None:
            palette = color_palette(n_colors=1)

        else:
            hue_names = utils.categorical_order(data[hue], hue_order)
            n_colors = len(hue_names)

            # By default use either the current color palette or HUSL
            if palette is None:
                current_palette = utils.get_color_cycle()
                if n_colors > len(current_palette):
                    colors = color_palette("husl", n_colors)
                else:
                    colors = color_palette(n_colors=n_colors)

            # Allow for palette to map from hue variable names
            elif isinstance(palette, dict):
                color_names = [palette[h] for h in hue_names]
                colors = color_palette(color_names, n_colors)

            # Otherwise act as if we just got a list of colors
            else:
                colors = color_palette(palette, n_colors)

            palette = color_palette(colors, n_colors)

        return palette


class PairGrid(Grid):
    """Subplot grid for plotting pairwise relationships in a dataset."""

    def __init__(self, data, hue=None, hue_order=None, palette=None,
                 hue_kws=None, vars=None, x_vars=None, y_vars=None,
                 diag_sharey=True, size=2.5, aspect=1,
                 despine=True, dropna=True):
        """Initialize the plot figure and PairGrid object.

        Parameters
        ----------
        data : DataFrame
            Tidy (long-form) dataframe where each column is a variable and
            each row is an observation.
        hue : string (variable name), optional
            Variable in ``data`` to map plot aspects to different colors.
        hue_order : list of strings
            Order for the levels of the hue variable in the palette
        palette : dict or geneview color palette
            Set of colors for mapping the ``hue`` variable. If a dict, keys
            should be values  in the ``hue`` variable.
        hue_kws : dictionary of param -> list of values mapping
            Other keyword arguments to insert into the plotting call to let
            other plot attributes vary across levels of the hue variable (e.g.
            the markers in a scatterplot).
        vars : list of variable names, optional
            Variables within ``data`` to use, otherwise use every column with
            a numeric datatype.
        {x, y}_vars : lists of variable names, optional
            Variables within ``data`` to use separately for the rows and
            columns of the figure; i.e. to make a non-square plot.
        size : scalar, optional
            Height (in inches) of each facet.
        aspect : scalar, optional
            Aspect * size gives the width (in inches) of each facet.
        despine : boolean, optional
            Remove the top and right spines from the plots.
        dropna : boolean, optional
            Drop missing values from the data before plotting.

        See Also
        --------
        pairplot : Easily drawing common uses of :class:`PairGrid`.
        FacetGrid : Subplot grid for plotting conditional relationships.

        Examples
        --------

        Draw a scatterplot for each pairwise relationship:

        .. plot::
            :context: close-figs

            >>> import matplotlib.pyplot as plt
            >>> import geneview as gv; gv.setup()
            >>> iris = gv.load_dataset("iris")
            >>> g = gv.PairGrid(iris)
            >>> g = g.map(plt.scatter)

        Show a univariate distribution on the diagonal:

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris)
            >>> g = g.map_diag(plt.hist)
            >>> g = g.map_offdiag(plt.scatter)

        (It's not actually necessary to catch the return value every time,
        as it is the same object, but it makes it easier to deal with the
        doctests).

        Color the points using a categorical variable:

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris, hue="species")
            >>> g = g.map(plt.scatter)
            >>> g = g.add_legend()

        Plot a subset of variables

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris, vars=["sepal_length", "sepal_width"])
            >>> g = g.map(plt.scatter)

        Pass additional keyword arguments to the functions

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris)
            >>> g = g.map_diag(plt.hist, edgecolor="w")
            >>> g = g.map_offdiag(plt.scatter, edgecolor="w", s=40)

        Use different variables for the rows and columns:

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris,
            ...                  x_vars=["sepal_length", "sepal_width"],
            ...                  y_vars=["petal_length", "petal_width"])
            >>> g = g.map(plt.scatter)

        Use different functions on the upper and lower triangles:

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris)
            >>> g = g.map_upper(plt.scatter)
            >>> g = g.map_lower(gv.kdeplot, cmap="Blues_d")
            >>> g = g.map_diag(gv.kdeplot, lw=3, legend=False)

        Use different colors and markers for each categorical level:

        .. plot::
            :context: close-figs

            >>> g = gv.PairGrid(iris, hue="species", palette="Set2",
            ...                 hue_kws={"marker": ["o", "s", "D"]})
            >>> g = g.map(plt.scatter, linewidths=1, edgecolor="w", s=40)
            >>> g = g.add_legend()

        """

        # Sort out the variables that define the grid
        if vars is not None:
            x_vars = list(vars)
            y_vars = list(vars)
        elif (x_vars is not None) or (y_vars is not None):
            if (x_vars is None) or (y_vars is None):
                raise ValueError("Must specify `x_vars` and `y_vars`")
        else:
            numeric_cols = self._find_numeric_cols(data)
            x_vars = numeric_cols
            y_vars = numeric_cols

        if np.isscalar(x_vars):
            x_vars = [x_vars]
        if np.isscalar(y_vars):
            y_vars = [y_vars]

        self.x_vars = list(x_vars)
        self.y_vars = list(y_vars)
        self.square_grid = self.x_vars == self.y_vars

        # Create the figure and the array of subplots
        figsize = len(x_vars) * size * aspect, len(y_vars) * size

        fig, axes = plt.subplots(len(y_vars), len(x_vars),
                                 figsize=figsize,
                                 sharex="col", sharey="row",
                                 squeeze=False)

        self.fig = fig
        self.axes = axes
        self.data = data

        # Save what we are going to do with the diagonal
        self.diag_sharey = diag_sharey
        self.diag_axes = None

        # Label the axes
        self._add_axis_labels()

        # Sort out the hue variable
        self._hue_var = hue
        if hue is None:
            self.hue_names = ["_nolegend_"]
            self.hue_vals = pd.Series(["_nolegend_"] * len(data),
                                      index=data.index)
        else:
            hue_names = utils.categorical_order(data[hue], hue_order)
            if dropna:
                # Filter NA from the list of unique hue names
                hue_names = list(filter(pd.notnull, hue_names))
            self.hue_names = hue_names
            self.hue_vals = data[hue]

        # Additional dict of kwarg -> list of values for mapping the hue var
        self.hue_kws = hue_kws if hue_kws is not None else {}

        self.palette = self._get_palette(data, hue, hue_order, palette)
        self._legend_data = {}

        # Make the plot look nice
        if despine:
            utils.despine(fig=fig)
        fig.tight_layout()

    def map(self, func, **kwargs):
        """Plot with the same function in every subplot.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, y_var in enumerate(self.y_vars):
            for j, x_var in enumerate(self.x_vars):
                hue_grouped = self.data.groupby(self.hue_vals)
                for k, label_k in enumerate(self.hue_names):

                    # Attempt to get data for this level, allowing for empty
                    try:
                        data_k = hue_grouped.get_group(label_k)
                    except KeyError:
                        data_k = pd.DataFrame(columns=self.data.columns,
                                              dtype=np.float)

                    ax = self.axes[i, j]
                    plt.sca(ax)

                    # Insert the other hue aesthetics if appropriate
                    for kw, val_list in self.hue_kws.items():
                        kwargs[kw] = val_list[k]

                    color = self.palette[k] if kw_color is None else kw_color
                    func(data_k[x_var], data_k[y_var],
                         label=label_k, color=color, **kwargs)

                self._clean_axis(ax)
                self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color
        self._add_axis_labels()

        return self

    def map_diag(self, func, **kwargs):
        """Plot with a univariate function on each diagonal subplot.

        Parameters
        ----------
        func : callable plotting function
            Must take an x array as a positional arguments and draw onto the
            "currently active" matplotlib Axes. There is a special case when
            using a ``hue`` variable and ``plt.hist``; the histogram will be
            plotted with stacked bars.

        """
        # Add special diagonal axes for the univariate plot
        if self.square_grid and self.diag_axes is None:
            diag_axes = []
            for i, (var, ax) in enumerate(zip(self.x_vars,
                                              np.diag(self.axes))):
                if i and self.diag_sharey:
                    diag_ax = ax._make_twin_axes(sharex=ax,
                                                 sharey=diag_axes[0],
                                                 frameon=False)
                else:
                    diag_ax = ax._make_twin_axes(sharex=ax, frameon=False)
                diag_ax.set_axis_off()
                diag_axes.append(diag_ax)
            self.diag_axes = np.array(diag_axes, np.object)

        # Plot on each of the diagonal axes
        for i, var in enumerate(self.x_vars):
            ax = self.diag_axes[i]
            hue_grouped = self.data[var].groupby(self.hue_vals)

            # Special-case plt.hist with stacked bars
            if func is plt.hist:
                plt.sca(ax)
                vals = []
                for label in self.hue_names:
                    # Attempt to get data for this level, allowing for empty
                    try:
                        vals.append(np.asarray(hue_grouped.get_group(label)))
                    except KeyError:
                        vals.append(np.array([]))
                func(vals, color=self.palette, histtype="barstacked",
                     **kwargs)
            else:
                for k, label_k in enumerate(self.hue_names):
                    # Attempt to get data for this level, allowing for empty
                    try:
                        data_k = hue_grouped.get_group(label_k)
                    except KeyError:
                        data_k = np.array([])
                    plt.sca(ax)
                    func(data_k, label=label_k,
                         color=self.palette[k], **kwargs)

            self._clean_axis(ax)

        self._add_axis_labels()

        return self

    def map_lower(self, func, **kwargs):
        """Plot with a bivariate function on the lower diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, j in zip(*np.tril_indices_from(self.axes, -1)):
            hue_grouped = self.data.groupby(self.hue_vals)
            for k, label_k in enumerate(self.hue_names):

                # Attempt to get data for this level, allowing for empty
                try:
                    data_k = hue_grouped.get_group(label_k)
                except KeyError:
                    data_k = pd.DataFrame(columns=self.data.columns,
                                          dtype=np.float)

                ax = self.axes[i, j]
                plt.sca(ax)

                x_var = self.x_vars[j]
                y_var = self.y_vars[i]

                # Insert the other hue aesthetics if appropriate
                for kw, val_list in self.hue_kws.items():
                    kwargs[kw] = val_list[k]

                color = self.palette[k] if kw_color is None else kw_color
                func(data_k[x_var], data_k[y_var], label=label_k,
                     color=color, **kwargs)

            self._clean_axis(ax)
            self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color
        self._add_axis_labels()

        return self

    def map_upper(self, func, **kwargs):
        """Plot with a bivariate function on the upper diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """
        kw_color = kwargs.pop("color", None)
        for i, j in zip(*np.triu_indices_from(self.axes, 1)):

            hue_grouped = self.data.groupby(self.hue_vals)

            for k, label_k in enumerate(self.hue_names):

                # Attempt to get data for this level, allowing for empty
                try:
                    data_k = hue_grouped.get_group(label_k)
                except KeyError:
                    data_k = pd.DataFrame(columns=self.data.columns,
                                          dtype=np.float)

                ax = self.axes[i, j]
                plt.sca(ax)

                x_var = self.x_vars[j]
                y_var = self.y_vars[i]

                # Insert the other hue aesthetics if appropriate
                for kw, val_list in self.hue_kws.items():
                    kwargs[kw] = val_list[k]

                color = self.palette[k] if kw_color is None else kw_color
                func(data_k[x_var], data_k[y_var], label=label_k,
                     color=color, **kwargs)

            self._clean_axis(ax)
            self._update_legend_data(ax)

        if kw_color is not None:
            kwargs["color"] = kw_color

        return self

    def map_offdiag(self, func, **kwargs):
        """Plot with a bivariate function on the off-diagonal subplots.

        Parameters
        ----------
        func : callable plotting function
            Must take x, y arrays as positional arguments and draw onto the
            "currently active" matplotlib Axes.

        """

        self.map_lower(func, **kwargs)
        self.map_upper(func, **kwargs)
        return self

    def _add_axis_labels(self):
        """Add labels to the left and bottom Axes."""
        for ax, label in zip(self.axes[-1, :], self.x_vars):
            ax.set_xlabel(label)
        for ax, label in zip(self.axes[:, 0], self.y_vars):
            ax.set_ylabel(label)

    def _find_numeric_cols(self, data):
        """Find which variables in a DataFrame are numeric."""
        # This can't be the best way to do this, but  I do not
        # know what the best way might be, so this seems ok
        numeric_cols = []
        for col in data:
            try:
                data[col].astype(np.float)
                numeric_cols.append(col)
            except (ValueError, TypeError):
                pass
        return numeric_cols
