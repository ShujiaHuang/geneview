""" plot 2-6 sets venn diagram

Author: Shujia Huang
Date: 2021-04-26 12:03:00

Thanks to the code from tctianchi and LankyCyril: https://github.com/tctianchi/pyvenn
"""
from matplotlib.pyplot import subplots
from matplotlib.colors import to_rgba
from matplotlib.patches import Ellipse, Polygon

from ..palette import generate_colors_palette

__all__ = ["venn", "generate_petal_labels"]

SHAPE_COORDS = {
    2: [(.375, .500), (.625, .500)],
    3: [(.333, .633), (.666, .633), (.500, .310)],
    4: [(.350, .400), (.450, .500), (.544, .500), (.644, .400)],
    5: [(.428, .449), (.469, .543), (.558, .523), (.578, .432), (.489, .383)],
    6: [(.637, .921, .649, .274, .188, .667),
        (.981, .769, .335, .191, .393, .671),
        (.941, .397, .292, .475, .456, .747),
        (.662, .119, .316, .548, .662, .700),
        (.309, .081, .374, .718, .681, .488),
        (.016, .626, .726, .687, .522, .327)]
}

SHAPE_DIMS = {
    2: [(.50, .50), (.50, .50)],
    3: [(.50, .50), (.50, .50), (.50, .50)],
    4: [(.72, .45), (.72, .45), (.72, .45), (.72, .45)],
    5: [(.87, .50), (.87, .50), (.87, .50), (.87, .50), (.87, .50)],
    6: [(None,)] * 6
}

SHAPE_ANGLES = {
    2: [0, 0],
    3: [0, 0, 0],
    4: [140, 140, 40, 40],
    5: [155, 82, 10, 118, 46],
    6: [None] * 6
}

PETAL_LABEL_COORDS = {
    2: {"01": (.74, .50), "10": (.26, .50), "11": (.50, .50)},
    3: {"001": (.500, .270), "010": (.730, .650), "011": (.610, .460),
        "100": (.270, .650), "101": (.390, .460), "110": (.500, .650),
        "111": (.500, .508)},
    4: {"0001": (.85, .42), "0010": (.68, .72), "0011": (.77, .59),
        "0100": (.32, .72), "0101": (.71, .30), "0110": (.50, .66),
        "0111": (.65, .50), "1000": (.14, .42), "1001": (.50, .17),
        "1010": (.29, .30), "1011": (.39, .24), "1100": (.23, .59),
        "1101": (.61, .24), "1110": (.35, .50), "1111": (.50, .38)},
    5: {"00001": (.27, .11), "00010": (.72, .11), "00011": (.55, .13),
        "00100": (.91, .58), "00101": (.78, .64), "00110": (.84, .41),
        "00111": (.76, .55), "01000": (.51, .90), "01001": (.39, .15),
        "01010": (.42, .78), "01011": (.50, .15), "01100": (.67, .76),
        "01101": (.70, .71), "01110": (.51, .74), "01111": (.64, .67),
        "10000": (.10, .61), "10001": (.20, .31), "10010": (.76, .25),
        "10011": (.65, .23), "10100": (.18, .50), "10101": (.21, .37),
        "10110": (.81, .37), "10111": (.74, .40), "11000": (.27, .70),
        "11001": (.34, .25), "11010": (.33, .72), "11011": (.51, .22),
        "11100": (.25, .58), "11101": (.28, .39), "11110": (.36, .66),
        "11111": (.51, .47)},
    6: {"000001": (.212, .562), "000010": (.430, .249), "000011": (.356, .444),
        "000100": (.609, .255), "000101": (.323, .546), "000110": (.513, .316),
        "000111": (.523, .348), "001000": (.747, .458), "001001": (.325, .492),
        "001010": (.670, .481), "001011": (.359, .478), "001100": (.653, .444),
        "001101": (.344, .526), "001110": (.653, .466), "001111": (.363, .503),
        "010000": (.750, .616), "010001": (.682, .654), "010010": (.402, .310),
        "010011": (.392, .421), "010100": (.653, .691), "010101": (.651, .644),
        "010110": (.490, .340), "010111": (.468, .399), "011000": (.692, .545),
        "011001": (.666, .592), "011010": (.665, .496), "011011": (.374, .470),
        "011100": (.653, .537), "011101": (.652, .579), "011110": (.653, .488),
        "011111": (.389, .486), "100000": (.553, .806), "100001": (.313, .604),
        "100010": (.388, .694), "100011": (.375, .633), "100100": (.605, .359),
        "100101": (.334, .555), "100110": (.582, .397), "100111": (.542, .372),
        "101000": (.468, .708), "101001": (.355, .572), "101010": (.420, .679),
        "101011": (.375, .597), "101100": (.641, .436), "101101": (.348, .538),
        "101110": (.635, .453), "101111": (.370, .548), "110000": (.594, .689),
        "110001": (.579, .670), "110010": (.398, .670), "110011": (.395, .653),
        "110100": (.633, .682), "110101": (.616, .656), "110110": (.587, .427),
        "110111": (.526, .415), "111000": (.495, .677), "111001": (.505, .648),
        "111010": (.428, .663), "111011": (.430, .631), "111100": (.639, .524),
        "111101": (.591, .604), "111110": (.622, .477), "111111": (.501, .523)}
}

DATASET_LEGEND_COORDS = {
    # n_set => [x, y, horizontalalignment, verticalalignment]
    2: {
        0: (0.20, 0.76, "right", "bottom"),
        1: (0.80, 0.76, "left", "bottom"),
    },
    3: {
        0: (0.15, 0.87, "right", "bottom"),
        1: (0.85, 0.87, "left", "bottom"),
        2: (0.50, 0.02, "center", "top"),
    },
    4: {
        0: (0.13, 0.18, "right", "center"),
        1: (0.18, 0.83, "right", "bottom"),
        2: (0.82, 0.83, "left", "bottom"),
        3: (0.87, 0.18, "left", "top"),
    },
    5: {
        0: (0.02, 0.72, "right", "center"),
        1: (0.72, 0.94, "center", "bottom"),
        2: (0.97, 0.74, "left", "center"),
        3: (0.88, 0.05, "left", "center"),
        4: (0.12, 0.05, "right", "center"),
    },
    6: {
        0: (0.674, 0.824, "center", "center"),
        1: (0.747, 0.751, "center", "center"),
        2: (0.739, 0.396, "center", "center"),
        3: (0.700, 0.247, "center", "center"),
        4: (0.291, 0.255, "center", "center"),
        5: (0.203, 0.484, "center", "center"),
    }
}


def generate_colors(cmap="viridis", n_colors=6, alpha=.4):
    """Generate colors from matplotlib colormap; pass list to use exact colors"""
    if not isinstance(n_colors, int) or (n_colors < 2) or (n_colors > 6):
        raise ValueError("n_colors must be an integer between 2 and 6")

    colors = generate_colors_palette(cmap=cmap, n_colors=n_colors, alpha=alpha)
    return colors[:n_colors]


def less_transparent_color(color, alpha_factor=2):
    """Bump up color's alpha"""
    new_alpha = (1 + to_rgba(color)[3]) / alpha_factor
    return to_rgba(color, alpha=new_alpha)


def draw_ellipse(ax, x, y, w, h, a, color):
    """Wrapper for drawing ellipse; called like `draw_ellipse(ax, *coords, *dims, angle, color)`"""
    ax.add_patch(
        Ellipse(
            xy=(x, y),
            width=w,
            height=h,
            angle=a,
            facecolor=color,
            edgecolor=less_transparent_color(color),
            lw=1,
        )
    )


def draw_triangle(ax, x1, y1, x2, y2, x3, y3, _dim, _angle, color):
    """Wrapper for drawing triangle; called like `draw_triangle(ax, *coords, None, None, color)`"""
    ax.add_patch(
        Polygon(
            xy=[(x1, y1), (x2, y2), (x3, y3)],
            closed=True,
            facecolor=color,
            edgecolor=less_transparent_color(color),
            lw=1,
        )
    )
    return


def draw_text(ax, x, y, text, color="black", fontsize=14, ha="center", va="center"):
    """Wrapper for drawing text"""
    ax.text(
        x, y, text,
        fontsize=fontsize,
        color=color,
        horizontalalignment=ha,
        verticalalignment=va,
    )
    return


def _generate_logics(n_sets):
    """Generate intersection identifiers in binary (0010 etc)"""
    for i in range(1, 2 ** n_sets):
        yield bin(i)[2:].zfill(n_sets)


def generate_petal_labels(datasets, fmt="{size}"):
    """Generate petal descriptions for venn diagram based on set sizes."""
    datasets = list(datasets)
    n_sets = len(datasets)
    if n_sets < 2 or n_sets > 6:
        raise ValueError("Number of sets must be between 2 and 6.")

    dataset_union = set.union(*datasets)
    universe_size = len(dataset_union)
    petal_labels = {}
    for logic in _generate_logics(n_sets):
        included_sets = [datasets[i] for i in range(n_sets) if logic[i] == "1"]
        excluded_sets = [datasets[i] for i in range(n_sets) if logic[i] == "0"]
        petal_set = (
                (dataset_union & set.intersection(*included_sets)) -
                set.union(set(), *excluded_sets)
        )
        petal_labels[logic] = fmt.format(
            logic=logic,
            size=len(petal_set),
            percentage=(100 * len(petal_set) / max(universe_size, 1))
        )
    return petal_labels


def _init_axes(ax):
    """Create axes if do not exist, set axes parameters"""
    if ax is None:
        _, ax = subplots(nrows=1, ncols=1, figsize=(7, 7))

    ax.set_axis_off()
    return ax


def _get_n_sets(petal_labels, dataset_labels):
    """Infer number of sets, check consistency"""
    n_sets = len(dataset_labels)
    petal_labels_set = set(_generate_logics(n_sets))
    for logic in petal_labels.keys():
        if len(logic) != n_sets:
            raise ValueError("Inconsistent petal and dataset labels: %d, %d" % (len(logic), n_sets))
        if not (set(logic) <= {"0", "1"}):
            raise KeyError("Key not understood: " + logic)
        if logic not in petal_labels_set:
            raise KeyError("'%s' is not a legal key." % logic)
    return n_sets


def _draw_venn(data, names=None, palette=None, alpha=0.4, fontsize=14,
               legend_use_petal_color=False, legend_loc=None, ax=None):
    """Draw Venn diagram, annotate petals and dataset labels.
    """
    DEFAULT_COLORS = [
        # r, g, b, a
        [0.361, 0.753, 0.384, 0.5],
        [0.353, 0.608, 0.831, 0.5],
        [0.965, 0.925, 0.337, 0.6],
        [0.945, 0.353, 0.376, 0.4],
        [1.000, 0.459, 0.000, 0.3],
        [0.322, 0.322, 0.745, 0.2]
    ]

    if not isinstance(names, list):
        raise ValueError("Names of sets should be a list and must not be empty.")

    n_sets = _get_n_sets(data, names)
    if 2 <= n_sets < 6:
        draw_shape = draw_ellipse
    elif n_sets == 6:
        draw_shape = draw_triangle
    else:
        raise ValueError("Number of sets must be between 2 and 6")

    ax = _init_axes(ax)
    if palette is None:
        palette = [DEFAULT_COLORS[i] for i in range(n_sets)]
    else:
        palette = generate_colors(n_colors=n_sets, cmap=palette, alpha=alpha)

    shape_params = zip(
        SHAPE_COORDS[n_sets],
        SHAPE_DIMS[n_sets],
        SHAPE_ANGLES[n_sets],
        palette
    )

    # Draw the shape for venn diagram
    for coords, dims, angle, color in shape_params:
        draw_shape(ax, *coords, *dims, angle, color)

    # annotate the value for each petal of venn plot
    for k, value in data.items():
        if k in PETAL_LABEL_COORDS[n_sets]:
            x, y = PETAL_LABEL_COORDS[n_sets][k]
            draw_text(ax, x, y, value, fontsize=fontsize)

    if legend_loc is not None:
        ax.legend(names, loc=legend_loc, prop={"size": fontsize})
    else:
        # plot the legend name for each dataset
        for i in range(n_sets):
            x, y, ha, va = DATASET_LEGEND_COORDS[n_sets][i]
            if legend_use_petal_color:
                c = palette[i]
                c[-1] = 1.0  # Set alpha blending value to opaque for legend text.
            else:
                c = "k"
            draw_text(ax, x, y, names[i], fontsize=fontsize + 2, color=c, ha=ha, va=va)

    return ax


def is_valid_dataset_dict(data):
    """Validate passed data (must be dictionary of sets)"""
    if not (hasattr(data, "keys") and hasattr(data, "values")):
        return False
    for dataset in data.values():
        if not isinstance(dataset, set):
            return False
    else:
        return True


def is_already_venn_dataset(petal_labels, dataset_labels):
    if not isinstance(dataset_labels, list):
        return False

    if not (hasattr(petal_labels, "keys") and hasattr(petal_labels, "values")):
        return False

    n_sets = len(dataset_labels)
    petal_labels_set = set(_generate_logics(n_sets))
    valid = True
    for logic in petal_labels.keys():
        if not isinstance(petal_labels[logic], str):
            valid = False
        if len(logic) != n_sets:
            valid = False
        if not (set(logic) <= {"0", "1"}):
            valid = False
        if logic not in petal_labels_set:
            valid = False
        if not valid:
            break

    return valid


def vennx(data, names=None, palette=None, alpha=0.4, fontsize=14,
          legend_use_petal_color=False, legend_loc=None, ax=None):
    """Generate venn diagram by input petal labels data.

    Parameters
    ----------
    data : dict, require
        A dictionaries of labels for dataset that will produce the venn diagram.
        The ``data`` variable should looks like:
        {'001': '0',
         '010': '5',
         '011': '0',
         '100': '3',
         '101': '2',
         '110': '2',
         '111': '3'}

    names : list, optional, default ``list(data.keys())``
        The label names for each petal in venn plot.

    palette : string, list, or :class:`matplotlib.colors.Colormap`, optional, default: DEFAULT_COLOR.
        String values are passed to :func:`generate_colors`. List values imply categorical
        mapping, while a colormap object implies numeric mapping.

    alpha : float scalar, default is 0.4, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    fontsize : integer, optional, default: 14
        Set the fontsize for plot.

    legend_use_petal_color : bool, optional, default: False
        Determine to use the same color as petals for legend or not.

    legend_loc : string. optional
        Place a legend on the axes. String value is passed to matplotlib :rc:`legend.loc`.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise create a default axis by plt.subplots() with figsize=(7, 7)
        in ``_draw_venn()``.

    Returns
    -------
    ax : matplotlib Axes
        Axes object with the venn plot.

    Examples
    --------

    Plot a three sets venn diagram.

    .. plot::
        :context: close-figs

        >>> from numpy.random import choice
        >>> from geneview.baseplot._venn import vennx, generate_petal_labels
        >>> dataset_dict = {name:set(choice(1000, 250, replace=False)) for name in list("ABCD")}
        >>> petal_labels = generate_petal_labels(dataset_dict.values(), fmt="{size}\\n({percentage:.1f}%)")
        >>> ax = vennx(data=petal_labels, names=list(dataset_dict.keys()))

    Set the labels for each petal by the keys of dataset.

    .. plot::
        :context: close-figs
        >>> ax = vennx(data=petal_labels, names=list(dataset_dict.keys()))

    Or set the labels by customer.

    .. plot::
        :context: close-figs
        >>> ax = vennx(data=petal_labels, names=["set 1", "set 2", "set 3", "set 4"],
        ...            legend_use_petal_color=True)
    """
    if not isinstance(names, list):
        raise ValueError("Names of sets should be a list and must not be empty.")

    if not is_already_venn_dataset(data, names):
        raise TypeError("``data`` is not a dict or the value is not a string. ")

    return _draw_venn(data=data,
                      names=names,
                      palette=palette,
                      alpha=alpha,
                      fontsize=fontsize,
                      legend_use_petal_color=legend_use_petal_color,
                      legend_loc=legend_loc,
                      ax=ax)


def venn(data, names=None, fmt="{size}", palette="viridis", alpha=0.4, fontsize=14,
         legend_use_petal_color=False, legend_loc=None, ax=None):
    """Check input, generate petal labels, draw venn diagram.

    Parameters
    ----------
    data : dict, require
        A dictionaries of dataset that will produce the venn diagram.

    names : list, optional, default ``list(data.keys())``
        The label names for each petal in venn plot.

    fmt : string, optional, default: "{size}"
        A Python 3 style format string that understands {size}, {percentage}, and {logic}.
        Here set the data format for petal datasets' lebals.

    palette : string, list, or :class:`matplotlib.colors.Colormap`, optional, default: "viridis".
        String values are passed to :func:`generate_colors`. List values imply categorical
        mapping, while a colormap object implies numeric mapping.

    alpha : float scalar, default is 0.4, optional
        The alpha blending value, between 0(transparent) and 1(opaque)

    fontsize : integer, optional, default: 14
        Set the fontsize for plot.

    legend_use_petal_color : bool, optional, default: False
        Determine to use the same color as petals for legend or not.

    legend_loc : string. optional
        Place a legend on the axes. String value is passed to matplotlib :rc:`legend.loc`.

    ax : matplotlib axis, optional
        Axis to plot on, otherwise create a default axis by plt.subplots() with figsize=(7, 7)
        in ``_draw_venn()``.

    Returns
    -------
    ax : matplotlib Axes
        Axes object with the venn plot.

    Examples
    --------

    Plot a minimal venn plot example.

    .. plot::
        :context: close-figs

        >>> from geneview import venn
        >>> data = {
        ...     "Dataset 1": {"A", "B", "D", "E"},
        ...     "Dataset 2": {"C", "F", "B", "G"},
        ...     "Dataset 3": {"J", "C", "K"}
        ... }
        >>> ax = venn(data)

   Rename the labels for each petal by manual and set the color of legend text to be the
   same as petal by setting bool variable as ``legend_use_petal_color=True``.

    .. plot::
        :context: close-figs

        >>> # Keep in mind that the order of ["A", "B", "C"] represent the patel: 001, 010 and 100, respectively.
        >>> ax = venn(data, names=["A", "B", "C"], legend_use_petal_color=True)

    `venn()` function could also use the prepared Venn data, directly.

    .. plot::
        :context: close-figs

        >>> data = {'011': 'ns', '101': '48', '110': '50', '111': 'ns'}  # A three sets of venn data
        >>> data_name = ["α", "β", "γ"]  # The names of the three data sets.
        >>> ax = venn(data=data, names=data_name, palette="plasma")

    Examples of Venn diagrams for various numbers of sets.

    Venn diagrams can be plotted for 2, 3, 4, or 5 sets using ellipses, and for 6 sets using triangles.
    The venn() function accepts optional arguments data, names, fmt, palette, alpha, fontsize,
    legend_loc and ax.

    .. plot::
        :context: close-figs

        >>> from itertools import chain, islice
        >>> from string import ascii_uppercase
        >>> from numpy.random import choice
        >>> import matplotlib.pyplot as plt
        >>> from geneview import venn
        >>> _, top_axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 5))
        >>> _, bot_axs = plt.subplots(ncols=2, nrows=1, figsize=(18, 8))
        >>> cmaps = ["cool", list("rgb"), "plasma", "viridis", "Set1"]
        >>> letters = iter(ascii_uppercase)
        >>> for n_sets, cmap, ax in zip(range(2, 7), cmaps, chain(top_axs, bot_axs)):
        ...    dataset_dict = {name: set(choice(1000, 700, replace=False)) for name in islice(letters, n_sets)}
        ...    _ = venn(dataset_dict,
        ...             fmt="{percentage:.1f}%",  # "{size}", "{logic}"
        ...             palette=cmap,
        ...             fontsize=12,
        ...             legend_use_petal_color=True,
        ...             legend_loc="upper left",
        ...             ax=ax)

    If necessary, the labels on the petals (i.e., various intersections in the Venn diagram) can be adjusted manually.
    For this, generate_petal_labels() can be called first to get the petal_labels dictionary, or you can deal the data
    by yourself, just make sure that the data is store in a dictionary and the format should be like:

        # The keys of dict is consist by 0 and 1, and the value should be str type.
        {'001': '0',
         '010': '5',
         '011': '0',
         '100': '3',
         '101': '2',
         '110': '2',
         '111': '3'}

    After modification, pass the dictionary to function venn().

    >>> from geneview import generate_petal_labels
    >>> dataset_dict = {name: set(choice(1000, 250, replace=False)) for name in list("ABCD")}
    >>> petal_labels = generate_petal_labels(dataset_dict.values(), fmt="{logic}\\n({percentage:.1f}%)")
    >>> ax = venn(data=petal_labels, names=list(dataset_dict.keys()), legend_use_petal_color=True)
    """
    if is_already_venn_dataset(data, names):
        return vennx(data=data,
                     names=names,
                     palette=palette,
                     alpha=alpha,
                     fontsize=fontsize,
                     legend_use_petal_color=legend_use_petal_color,
                     legend_loc=legend_loc, ax=ax)

    elif not is_valid_dataset_dict(data):
        raise TypeError("Only dictionaries of sets are understood")

    return vennx(data=generate_petal_labels(data.values(), fmt=fmt),
                 names=list(data.keys()) if names is None else names,
                 palette=palette,
                 alpha=alpha,
                 fontsize=fontsize,
                 legend_use_petal_color=legend_use_petal_color,
                 legend_loc=legend_loc,
                 ax=ax)
