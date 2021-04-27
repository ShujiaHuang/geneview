""" plot 2-6 sets venn diagram

Author: Shujia Huang
Date: 2021-04-26 12:03:00

Thanks to the code from tctianchi and LankyCyril: https://github.com/tctianchi/pyvenn
"""
from matplotlib.pyplot import gca
from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Ellipse, Polygon

__all__ = ["venn"]

SHAPE_COORDS = {
    2: [(.375, .500), (.625, .500)],
    3: [(.333, .633), (.666, .633), (.500, .310)],
    4: [(.350, .400), (.450, .500), (.544, .500), (.644, .400)],
    5: [(.428, .449), (.469, .543), (.558, .523), (.578, .432), (.489, .383)],
    6: [
        (.637, .921, .649, .274, .188, .667),
        (.981, .769, .335, .191, .393, .671),
        (.941, .397, .292, .475, .456, .747),
        (.662, .119, .316, .548, .662, .700),
        (.309, .081, .374, .718, .681, .488),
        (.016, .626, .726, .687, .522, .327)
    ]
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
    2: {
        "01": (.74, .50), "10": (.26, .50), "11": (.50, .50)
    },
    3: {
        "001": (.500, .270), "010": (.730, .650), "011": (.610, .460),
        "100": (.270, .650), "101": (.390, .460), "110": (.500, .650),
        "111": (.500, .508)
    },
    4: {
        "0001": (.85, .42), "0010": (.68, .72), "0011": (.77, .59),
        "0100": (.32, .72), "0101": (.71, .30), "0110": (.50, .66),
        "0111": (.65, .50), "1000": (.14, .42), "1001": (.50, .17),
        "1010": (.29, .30), "1011": (.39, .24), "1100": (.23, .59),
        "1101": (.61, .24), "1110": (.35, .50), "1111": (.50, .38)
    },
    5: {
        "00001": (.27, .11), "00010": (.72, .11), "00011": (.55, .13),
        "00100": (.91, .58), "00101": (.78, .64), "00110": (.84, .41),
        "00111": (.76, .55), "01000": (.51, .90), "01001": (.39, .15),
        "01010": (.42, .78), "01011": (.50, .15), "01100": (.67, .76),
        "01101": (.70, .71), "01110": (.51, .74), "01111": (.64, .67),
        "10000": (.10, .61), "10001": (.20, .31), "10010": (.76, .25),
        "10011": (.65, .23), "10100": (.18, .50), "10101": (.21, .37),
        "10110": (.81, .37), "10111": (.74, .40), "11000": (.27, .70),
        "11001": (.34, .25), "11010": (.33, .72), "11011": (.51, .22),
        "11100": (.25, .58), "11101": (.28, .39), "11110": (.36, .66),
        "11111": (.51, .47)
    },
    6: {
        "000001": (.212, .562), "000010": (.430, .249), "000011": (.356, .444),
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
        "111101": (.591, .604), "111110": (.622, .477), "111111": (.501, .523)
    }
}

LABEL_LEGEND_COORDS = {
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

COLORS = [
    # r, g, b, a
    [0.361, 0.753, 0.384, 0.5],
    [0.353, 0.608, 0.831, 0.5],
    [0.965, 0.925, 0.337, 0.6],
    [0.945, 0.353, 0.376, 0.4],
    [1.000, 0.459, 0.000, 0.3],
    [0.322, 0.322, 0.745, 0.2]
]


def generate_colors(cmap="viridis", n_colors=6, alpha=.4):
    """Generate colors from matplotlib colormap; pass list to use exact colors"""
    if not isinstance(n_colors, int) or (n_colors < 2) or (n_colors > 6):
        raise ValueError("n_colors must be an integer between 2 and 6")
    if isinstance(cmap, list):
        colors = [to_rgba(color, alpha=alpha) for color in cmap]
    else:
        scalar_mappable = ScalarMappable(cmap=cmap)
        colors = scalar_mappable.to_rgba(range(n_colors), alpha=alpha).tolist()
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


def _init_axes(ax):
    """Create axes if do not exist, set axes parameters"""
    if ax is None:
        ax = gca()

    ax.set_axis_off()
    return ax


def get_n_sets(petal_labels, dataset_labels):
    """Infer number of sets, check consistency"""
    n_sets = len(dataset_labels)
    petal_labels_set = _generate_logics(n_sets)
    for logic in petal_labels.keys():
        if len(logic) != n_sets:
            raise ValueError("Inconsistent petal and dataset labels: %d, %d" % (len(logic), n_sets))
        if not (set(logic) <= {"0", "1"}):
            raise KeyError("Key not understood: " + logic)
        if logic not in petal_labels_set:
            raise KeyError("Key is not legal: " + logic)
    return n_sets


def venn(dataset, dataset_labels, palette=None, alpha=0.4, fontsize=14, legend_loc=None, ax=None):
    """Draw Venn diagram, annotate petals and dataset labels.
    """
    n_sets = get_n_sets(dataset, dataset_labels)
    if 2 <= n_sets < 6:
        draw_shape = draw_ellipse
    elif n_sets == 6:
        draw_shape = draw_triangle
    else:
        raise ValueError("Number of sets must be between 2 and 6")

    ax = _init_axes(ax)
    if palette is None:
        palette = [COLORS[i] for i in range(n_sets)]
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
    for k, value in dataset.items():
        if k in PETAL_LABEL_COORDS[n_sets]:
            x, y = PETAL_LABEL_COORDS[n_sets][k]
            draw_text(ax, x, y, value, fontsize=fontsize)

    # plot the name for each dataset
    for i in range(n_sets):
        x, y, ha, va = LABEL_LEGEND_COORDS[n_sets][i]
        draw_text(ax, x, y, dataset_labels[i], fontsize=fontsize + 2, ha=ha, va=va)

    if legend_loc is not None:
        ax.legend(dataset_labels, loc=legend_loc, bbox_to_anchor=(1.0, 0.5), fontsize=fontsize, fancybox=True)

    return ax
