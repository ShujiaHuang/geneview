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

