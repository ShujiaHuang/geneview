"""
Shape-drawing primitives for mutation plots.

Provides functions to draw circles, pie charts, fans, flags and pins at
specific *genomic* x-coordinates with y in axes-fraction [0, 1].

Shapes use **Ellipse patches in data coordinates** whose width/height are
computed so they appear perfectly circular on screen regardless of the
data-axis aspect ratio.
"""

from typing import Optional, List, Tuple

import numpy as np
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.transforms as mtransforms


def _blended(ax):
    """Return a transform: x in *data* coords, y in *axes* fraction [0,1]."""
    return mtransforms.blended_transform_factory(ax.transData, ax.transAxes)


def _data_radius(ax, radius_axes: float) -> Tuple[float, float]:
    """Return (rx_data, ry_data) for a circle of *radius_axes* axes-fraction.

    rx_data : half-width in data x-units
    ry_data : half-height in data y-units (y range is [0,1] so == axes frac)
    """
    fig = ax.figure
    pos = ax.get_position()
    w_inches = pos.width * fig.get_figwidth()
    h_inches = pos.height * fig.get_figheight()
    xlim = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]  # normally 1

    # Desired pixel radius
    r_px = radius_axes * h_inches * fig.dpi
    # Convert to data units
    rx = r_px / (w_inches * fig.dpi / x_range)
    ry = r_px / (h_inches * fig.dpi / y_range)
    return rx, ry


def draw_circle(ax, x, y, radius_axes, facecolor, edgecolor="black",
                linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a circle centred at (x-data, y-axes) with *radius_axes* size.

    Uses an Ellipse patch in data coordinates that appears circular on screen.
    """
    rx, ry = _data_radius(ax, radius_axes)
    ellipse = mpatches.Ellipse(
        (x, y), width=2 * rx, height=2 * ry,
        facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(ellipse)


def draw_square(ax, x, y, radius_axes, facecolor, edgecolor="black",
                linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a square centred at (x-data, y-axes) with *radius_axes* size."""
    rx, ry = _data_radius(ax, radius_axes)
    rect = mpatches.Rectangle(
        (x - rx, y - ry), 2 * rx, 2 * ry,
        facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(rect)


def draw_diamond(ax, x, y, radius_axes, facecolor, edgecolor="black",
                 linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a diamond (rotated square) centred at (x-data, y-axes)."""
    rx, ry = _data_radius(ax, radius_axes)
    verts = [(x, y - ry), (x + rx, y), (x, y + ry), (x - rx, y), (x, y - ry)]
    codes = [mpath.Path.MOVETO, mpath.Path.LINETO, mpath.Path.LINETO,
             mpath.Path.LINETO, mpath.Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(
        path, facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(patch)


def draw_triangle_up(ax, x, y, radius_axes, facecolor, edgecolor="black",
                     linewidth=1.0, alpha=1.0, zorder=5):
    """Draw an upward-pointing triangle centred at (x-data, y-axes)."""
    rx, ry = _data_radius(ax, radius_axes)
    verts = [(x, y + ry), (x - rx, y - ry), (x + rx, y - ry), (x, y + ry)]
    codes = [mpath.Path.MOVETO, mpath.Path.LINETO, mpath.Path.LINETO,
             mpath.Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(
        path, facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(patch)


def draw_triangle_down(ax, x, y, radius_axes, facecolor, edgecolor="black",
                       linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a downward-pointing triangle centred at (x-data, y-axes)."""
    rx, ry = _data_radius(ax, radius_axes)
    verts = [(x, y - ry), (x - rx, y + ry), (x + rx, y + ry), (x, y - ry)]
    codes = [mpath.Path.MOVETO, mpath.Path.LINETO, mpath.Path.LINETO,
             mpath.Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(
        path, facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(patch)


def draw_node_label(ax, x, y, radius_axes, text, color="white",
                    fontsize=6, zorder=7):
    """Render *text* centred inside a shape at (x-data, y-axes)."""
    ax.text(
        x, y, str(text),
        ha="center", va="center",
        fontsize=fontsize, color=color,
        zorder=zorder, clip_on=False,
    )


# Registry mapping shape names to draw functions
_SHAPE_DISPATCH = {
    "circle": draw_circle,
    "square": draw_square,
    "diamond": draw_diamond,
    "triangle_point_up": draw_triangle_up,
    "triangle_point_down": draw_triangle_down,
}


def draw_shape(ax, shape_name, x, y, radius_axes, **kwargs):
    """Generic shape dispatcher.

    Parameters
    ----------
    shape_name : str
        One of ``'circle'``, ``'square'``, ``'diamond'``,
        ``'triangle_point_up'``, ``'triangle_point_down'``.
    **kwargs
        Passed to the underlying draw function (``facecolor``, ``edgecolor``,
        ``linewidth``, ``alpha``, ``zorder``).
    """
    fn = _SHAPE_DISPATCH.get(shape_name, draw_circle)
    fn(ax, x, y, radius_axes, **kwargs)


def draw_pie(ax, x, y, radius_axes, values, colors, edgecolor="black",
             linewidth=0.8, alpha=1.0, zorder=5):
    """Draw a small pie chart at (x-data, y-axes).

    Parameters
    ----------
    values : list of float
        Relative sizes of pie slices.
    colors : list of str
        Colours for each slice.
    """
    rx, ry = _data_radius(ax, radius_axes)
    total = sum(values)
    if total <= 0:
        values = [1.0]
        colors = ["#CCCCCC"]
        total = 1.0

    theta1 = 0.0
    for val, col in zip(values, colors):
        theta2 = theta1 + 360.0 * val / total
        # Build pie slice as a path in data coords (ellipse-aware)
        n_pts = max(3, int(abs(theta2 - theta1) / 5))
        angles = np.linspace(np.radians(theta1), np.radians(theta2), n_pts)
        verts = [(x, y)]
        for a in angles:
            verts.append((x + rx * np.cos(a), y + ry * np.sin(a)))
        verts.append((x, y))
        codes = [mpath.Path.MOVETO] + \
                [mpath.Path.LINETO] * (len(verts) - 2) + \
                [mpath.Path.CLOSEPOLY]
        path = mpath.Path(verts, codes)
        wedge = mpatches.PathPatch(
            path, facecolor=col, edgecolor=edgecolor,
            linewidth=linewidth, alpha=alpha, zorder=zorder,
        )
        ax.add_patch(wedge)
        theta1 = theta2


def draw_fan(ax, x, y, radius_axes, score, color, edgecolor="black",
             linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a fan (sector) at (x-data, y-axes).

    Parameters
    ----------
    score : float
        Fill fraction [0, 1].  Controls the angular span of the fan.
    """
    score = max(0.0, min(1.0, score))
    rx, ry = _data_radius(ax, radius_axes)
    half_angle = 180.0 * score

    # Use a wedge shape
    n_pts = max(3, int(half_angle / 3))
    angles = np.linspace(np.radians(90 - half_angle),
                         np.radians(90 + half_angle), n_pts)
    verts = [(x, y)]
    for a in angles:
        verts.append((x + rx * np.cos(a), y + ry * np.sin(a)))
    verts.append((x, y))

    codes = [mpath.Path.MOVETO] + \
            [mpath.Path.LINETO] * (len(verts) - 2) + \
            [mpath.Path.CLOSEPOLY]
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(
        path, facecolor=color if color else "white",
        edgecolor=edgecolor if edgecolor else "black",
        linewidth=linewidth, alpha=alpha, zorder=zorder,
    )
    ax.add_patch(patch)


def draw_pin(ax, x, y, radius_axes, facecolor, edgecolor="black",
             linewidth=1.0, alpha=1.0, zorder=5):
    """Draw a map-pin / teardrop at (x-data, y-axes).

    The pin has a circular body above and a pointed tip at (x, y).
    Mirrors trackViewer's ``map-pin-red.xml`` shape.
    """
    rx, ry = _data_radius(ax, radius_axes)

    # Pin proportions (matching R: width=2r, height=3r)
    body_ry = ry * 1.0       # vertical radius of the circular body
    body_rx = rx * 1.0       # horizontal radius of the circular body
    body_cy = y + body_ry * 2.0  # body centre sits 2*ry above the tip

    # Build the teardrop path: point at bottom, circle at top
    # Tangent lines from tip (x, y) to the ellipse body
    n_arc = 40
    # Arc from ~210° around through 90° (top) to ~330° (avoiding bottom)
    angles = np.linspace(np.radians(210), np.radians(330) + 2 * np.pi, n_arc)
    # Normalise to [0, 2pi]
    angles = np.mod(angles, 2 * np.pi)

    verts = [(x, y)]  # start at the tip
    codes = [mpath.Path.MOVETO]
    for a in angles:
        verts.append((x + body_rx * np.cos(a), body_cy + body_ry * np.sin(a)))
        codes.append(mpath.Path.LINETO)
    verts.append((x, y))  # back to tip
    codes.append(mpath.Path.CLOSEPOLY)

    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(
        path, facecolor=facecolor, edgecolor=edgecolor,
        linewidth=linewidth, alpha=alpha, zorder=zorder,
        clip_on=False,
    )
    ax.add_patch(patch)

    # Inner highlight dot (white circle inside the body)
    dot_rx = body_rx * 0.35
    dot_ry = body_ry * 0.35
    dot = mpatches.Ellipse(
        (x, body_cy), width=2 * dot_rx, height=2 * dot_ry,
        facecolor="white", edgecolor="white",
        linewidth=0, alpha=alpha * 0.7, zorder=zorder + 1,
        clip_on=False,
    )
    ax.add_patch(dot)


def draw_flag(ax, x, y, radius_axes, color, edgecolor="black",
              linewidth=1.0, cex=1.0, zorder=5, label=""):
    """Draw a flag banner at (x-data, y-axes).

    A rectangular banner with optional text label inside, matching
    trackViewer's flag type (a horizontal rectangle, not a pennant).
    """
    rx, ry = _data_radius(ax, radius_axes)
    fig = ax.figure
    pos = ax.get_position()
    xlim = ax.get_xlim()
    x_range = xlim[1] - xlim[0]
    w_inches = pos.width * fig.get_figwidth()
    h_inches = pos.height * fig.get_figheight()

    # Banner dimensions in data-y units
    banner_h = ry * 2.5          # banner height
    # Banner width: estimate from label text (≈ character width)
    char_px = 5.0 * cex          # approximate character width in pixels
    label_len = max(len(str(label)), 3)
    banner_w_px = label_len * char_px + 6  # padding
    # Convert pixel width to data-x units
    banner_dx = banner_w_px / (w_inches * fig.dpi) * x_range

    # Banner bottom-left at (x, y), extending right and up
    bx, by = x, y
    rect = mpatches.FancyBboxPatch(
        (bx, by), banner_dx, banner_h,
        boxstyle="round,pad=0",
        facecolor=color if color else "#0080FF",
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
        clip_on=False,
    )
    ax.add_patch(rect)

    # Draw label text inside the banner
    if label:
        ax.text(
            bx + banner_dx * 0.5, by + banner_h * 0.5,
            str(label),
            ha="center", va="center",
            fontsize=5 * cex,
            color="white" if _is_dark(color) else "black",
            zorder=zorder + 1,
            clip_on=False,
        )


def _is_dark(color):
    """Rough check whether *color* is dark (for choosing text colour)."""
    try:
        import matplotlib.colors as mcolors
        rgb = mcolors.to_rgb(color)
        luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
        return luminance < 0.5
    except Exception:
        return True


def draw_stem(ax, x, y_start, y_end, color="#CCCCCC", linewidth=0.8,
              linestyle="--", zorder=2):
    """Draw a vertical stem line at x (data) from y_start to y_end (axes)."""
    # y is in axes fraction [0,1], use blended transform
    ax.plot(
        [x, x], [y_start, y_end],
        color=color, linewidth=linewidth,
        linestyle=linestyle, zorder=zorder,
        transform=_blended(ax),
        clip_on=False,
    )


def draw_stem_bent(ax, x0, y0, x1, y1, x2, y2, color="black",
                   linewidth=1.0, zorder=2):
    """Draw a bent stem: (x0,y0) -> (x1,y1) -> (x2,y2).

    x values are in data coords; y values in axes fraction.
    """
    ax.plot(
        [x0, x1, x2], [y0, y1, y2],
        color=color, linewidth=linewidth,
        zorder=zorder, clip_on=False,
        transform=_blended(ax),
    )
