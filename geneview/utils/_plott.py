"""
Small plotting-related utility functions.
=========================================

Copy from: the seaborn repo's utils.py on github, Then we make some modifies 
for geneview.

"""
import colorsys
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplcol

from distutils.version import LooseVersion
mpl_ge_150 = LooseVersion(mpl.__version__) >= "1.5.0"


class General(object):
    def __init__(self):
        pass

    @staticmethod
    def get_figure(is_show, fig_name=None, dpi="figure"):
        if fig_name:
            plt.savefig(fig_name, bbox_inches="tight", dpi=dpi)

        if is_show:
            plt.show()

        plt.clf()
        plt.close()

        return


def set_hls_values(color, h=None, l=None, s=None):
    """Independently manipulate the h, l, or s channels of a color.

    Parameters
    ----------
    color : matplotlib color
        hex, rgb-tuple, or html color name
    h, l, s : floats between 0 and 1, or None
        new values for each channel in hls space

    Returns
    -------
    new_color : rgb tuple
        new color code in RGB tuple representation

    """
    # Get rgb tuple representation
    rgb = mplcol.colorConverter.to_rgb(color)
    vals = list(colorsys.rgb_to_hls(*rgb))
    for i, val in enumerate([h, l, s]):
        if val is not None:
            vals[i] = val

    rgb = colorsys.hls_to_rgb(*vals)
    return rgb


def axlabel(xlabel, ylabel, **kwargs):
    """Grab current axis and label it."""
    ax = plt.gca()
    ax.set_xlabel(xlabel, **kwargs)
    ax.set_ylabel(ylabel, **kwargs)


def despine(fig=None, ax=None, top=True, right=True, left=False,
            bottom=False, offset=None, trim=False):
    """
    Remove the top and right spines from plot(s).

    Parametes
    ---------
    fig : matplotlib figure, optional
        Figure to despine all axes of, default uses current figure.

    ax : matplotlib axes, optional
        Specific axes object to despine.

    top, right, left, bottom : boolean, optional
        If True, remove that spine.

    offset : int or None  (default), optional
        Absolute distance, in points, spines should be moved away
        from the axes (negative values move spines inward).

    trim : bool, optional
        If true, limit spines to the smallest and largest major tick
        on each non-despined axis.

    Returns
    -------
    None
    """
    # Get references to the axes we want
    if fig is None and ax is None:
        axes = plt.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        for side in ["top", "right", "left", "bottom"]:
            # Toggle the spine objects
            is_visible = not locals()[side]
            ax_i.spines[side].set_visible(is_visible)
            if offset is not None and is_visible:
                _set_spine_position(ax_i.spines[side], ('outward', offset))

        # Set the ticks appropriately
        if bottom:
            ax_i.xaxis.tick_top()
        if top:
            ax_i.xaxis.tick_bottom()
        if left:
            ax_i.yaxis.tick_right()
        if right:
            ax_i.yaxis.tick_left()

        if trim:
            # clip off the parts of the spines that extend past major ticks
            xticks = ax_i.get_xticks()
            if xticks.size:
                firsttick = np.compress(xticks >= min(ax_i.get_xlim()), 
                                        xticks)[0]
                lasttick = np.compress(xticks <= max(ax_i.get_xlim()),
                                       xticks)[-1]
                ax_i.spines['bottom'].set_bounds(firsttick, lasttick)
                ax_i.spines['top'].set_bounds(firsttick, lasttick)
                newticks = xticks.compress(xticks <= lasttick)
                newticks = newticks.compress(newticks >= firsttick)
                ax_i.set_xticks(newticks)

            yticks = ax_i.get_yticks()
            if yticks.size:
                firsttick = np.compress(yticks >= min(ax_i.get_ylim()),
                                        yticks)[0]
                lasttick = np.compress(yticks <= max(ax_i.get_ylim()),
                                       yticks)[-1]
                ax_i.spines['left'].set_bounds(firsttick, lasttick)
                ax_i.spines['right'].set_bounds(firsttick, lasttick)
                newticks = yticks.compress(yticks <= lasttick)
                newticks = newticks.compress(newticks >= firsttick)
                ax_i.set_yticks(newticks)


def offset_spines(offset=10, fig=None, ax=None):
    """
    Simple function to offset spines away from axes.
    Use this immediately after creating figure and axes objects.
    Offsetting spines after plotting or manipulating the axes
    objects may result in loss of labels, ticks, and formatting.
    Parameters
    ----------
    offset : int, optional
        Absolute distance, in points, spines should be moved away
        from the axes (negative values move spines inward).
    fig : matplotlib figure, optional
        Figure to despine all axes of, default uses current figure.
    ax : matplotlib axes, optional
        Specific axes object to despine
    Returns
    -------
    None
    """
    warn_msg = "`offset_spines` is deprecated and will be removed in v0.5"
    warnings.warn(warn_msg, UserWarning)

    # Get references to the axes we want
    if fig is None and ax is None:
        axes = plt.gcf().axes
    elif fig is not None:
        axes = fig.axes
    elif ax is not None:
        axes = [ax]

    for ax_i in axes:
        for spine in ax_i.spines.values():
            _set_spine_position(spine, ('outward', offset))


def _set_spine_position(spine, position):
    """
    Set the spine's position without resetting an associated axis.
    As of matplotlib v. 1.0.0, if a spine has an associated axis, then
    spine.set_position() calls axis.cla(), which resets locators, formatters,
    etc.  We temporarily replace that call with axis.reset_ticks(), which is
    sufficient for our purposes.
    """
    axis = spine.axis
    if axis is not None:
        cla = axis.cla
        axis.cla = axis.reset_ticks
    spine.set_position(position)
    if axis is not None:
        axis.cla = cla


def get_color_cycle():
    if mpl_ge_150:
        cyl = mpl.rcParams['axes.prop_cycle']
        # matplotlib 1.5 verifies that axes.prop_cycle *is* a cycler
        # but no garuantee that there's a `color` key.
        # so users could have a custom rcParmas w/ no color...
        try:
            return [x['color'] for x in cyl]
        except KeyError:
            pass  # just return axes.color style below
    return mpl.rcParams['axes.color_cycle']
