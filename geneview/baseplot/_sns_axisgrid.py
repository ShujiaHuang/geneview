from __future__ import division

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from ..util import despine


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
        despine(f)
        despine(ax=ax_marg_x, left=True)
        despine(ax=ax_marg_y, bottom=True)
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
