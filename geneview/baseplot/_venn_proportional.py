"""Area-proportional 2- and 3-set Venn diagrams (circle-based).

Pure matplotlib, no extra dependencies.  Circle areas are proportional to the
set sizes and pairwise overlap areas are solved numerically (bisection on the
circle--circle "lens" area) so that the drawn overlaps match the real
intersection sizes.

For 3 sets the layout is the standard approximation used by area-proportional
Venn tools: circle radii encode set sizes, pairwise centre distances encode
pairwise overlaps, and the third circle is placed by triangulation.  Region
labels are positioned at numerically-estimated region centroids, which keeps
each label inside the region it describes regardless of the geometry.

Author: Shujia Huang
"""
from math import acos, sqrt, pi

import numpy as np
from matplotlib.patches import Circle

from ._venn import less_transparent_color, draw_text


def _clip(v):
    """Clamp a value to ``[-1, 1]`` for safe :func:`math.acos` input."""
    return max(-1.0, min(1.0, v))


def _lens_area(r1, r2, d):
    """Intersection area of two circles (radii r1, r2, centre distance d)."""
    if d <= abs(r1 - r2):
        return pi * min(r1, r2) ** 2
    if d >= r1 + r2:
        return 0.0
    d1 = (d * d - r2 * r2 + r1 * r1) / (2 * d)
    d2 = d - d1
    a1 = r1 * r1 * acos(_clip(d1 / r1)) - d1 * sqrt(max(0.0, r1 * r1 - d1 * d1))
    a2 = r2 * r2 * acos(_clip(d2 / r2)) - d2 * sqrt(max(0.0, r2 * r2 - d2 * d2))
    return a1 + a2


def _solve_distance(r1, r2, target):
    """Centre distance ``d`` such that the lens area equals *target*.

    ``_lens_area`` decreases monotonically with ``d`` on ``[|r1-r2|, r1+r2]``,
    so a simple bisection converges quickly.
    """
    if r1 <= 0 or r2 <= 0 or target <= 0:
        return r1 + r2  # disjoint (or degenerate) -> tangent / separated
    if target >= pi * min(r1, r2) ** 2:
        return abs(r1 - r2)  # full containment
    lo, hi = abs(r1 - r2), r1 + r2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if _lens_area(r1, r2, mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _radius(area):
    """Radius of a circle whose area (in count units) equals *area*."""
    return sqrt(max(area, 0.0) / pi)


def _region_label_positions(centers, radii, logics, bbox, n=400):
    """Estimate a representative interior point for each region via grid sampling.

    A region is the set of grid points inside exactly the circles whose bit is
    ``"1"`` in its *logic* string; the centroid of those points is a robust
    label anchor that stays inside the region for any geometry.
    """
    xmin, xmax, ymin, ymax = bbox
    xs = np.linspace(xmin, xmax, n)
    ys = np.linspace(ymin, ymax, n)
    X, Y = np.meshgrid(xs, ys)
    inside = [((X - cx) ** 2 + (Y - cy) ** 2) <= r * r
              for (cx, cy), r in zip(centers, radii)]

    positions = {}
    for logic in logics:
        mask = np.ones(X.shape, dtype=bool)
        for i, bit in enumerate(logic):
            mask &= inside[i] if bit == "1" else ~inside[i]
        if mask.any():
            positions[logic] = (float(X[mask].mean()), float(Y[mask].mean()))
    return positions


def _layout_venn2(sizes):
    """Return circle ``centers`` and ``radii`` for a 2-set proportional Venn."""
    t0 = sizes["10"] + sizes["11"]
    t1 = sizes["01"] + sizes["11"]
    r0, r1 = _radius(t0), _radius(t1)
    d = _solve_distance(r0, r1, sizes["11"])
    return [(0.0, 0.0), (d, 0.0)], [r0, r1]


def _layout_venn3(sizes):
    """Return circle ``centers`` and ``radii`` for a 3-set proportional Venn."""
    t0 = sizes["100"] + sizes["110"] + sizes["101"] + sizes["111"]
    t1 = sizes["010"] + sizes["110"] + sizes["011"] + sizes["111"]
    t2 = sizes["001"] + sizes["101"] + sizes["011"] + sizes["111"]
    r0, r1, r2 = _radius(t0), _radius(t1), _radius(t2)

    d01 = _solve_distance(r0, r1, sizes["110"] + sizes["111"])
    d02 = _solve_distance(r0, r2, sizes["101"] + sizes["111"])
    d12 = _solve_distance(r1, r2, sizes["011"] + sizes["111"])

    if d01 < 1e-9:
        d01 = max(r0, r1, 1e-3)
    c0 = (0.0, 0.0)
    c1 = (d01, 0.0)
    # Triangulate the third centre from its distances to the first two.
    x2 = (d02 * d02 - d12 * d12 + d01 * d01) / (2 * d01)
    y2 = sqrt(max(0.0, d02 * d02 - x2 * x2))  # clamp if triangle inequality fails
    return [c0, c1, (x2, y2)], [r0, r1, r2]


def _draw_circles(ax, centers, radii, names, colors):
    """Add the filled circle patches (labelled for an optional legend)."""
    for (cx, cy), r, color, name in zip(centers, radii, colors, names):
        ax.add_patch(Circle(
            (cx, cy), r,
            facecolor=color,
            edgecolor=less_transparent_color(color),
            lw=1,
            label=name,
        ))


def _place_names(ax, centers, radii, names, colors, fontsize,
                 legend_use_petal_color):
    """Place each set name just outside its circle, radiating from the centroid."""
    cxm = sum(c[0] for c in centers) / len(centers)
    cym = sum(c[1] for c in centers) / len(centers)
    offset = 0.15 * max(radii + [1e-6])
    for (cx, cy), r, color, name in zip(centers, radii, colors, names):
        dx, dy = cx - cxm, cy - cym
        norm = sqrt(dx * dx + dy * dy)
        if norm < 1e-9:  # single circle / coincident centres -> push upward
            dx, dy, norm = 0.0, 1.0, 1.0
        ux, uy = dx / norm, dy / norm
        tx, ty = cx + ux * (r + offset), cy + uy * (r + offset)
        if legend_use_petal_color:
            c = list(color)
            c[-1] = 1.0  # opaque for legible text
        else:
            c = "k"
        draw_text(ax, tx, ty, name, color=c, fontsize=fontsize + 2)


def draw_proportional_venn(ax, sizes, labels, names, colors, fontsize=14,
                           legend_use_petal_color=False, legend_loc=None):
    """Draw an area-proportional Venn diagram for 2 or 3 sets.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes (assumed to already have its frame turned off).
    sizes : dict
        ``{logic: int}`` element counts per petal (see ``_generate_petal_sizes``).
    labels : dict
        ``{logic: str}`` formatted labels drawn inside each region.
    names : list
        Set names, in set order.
    colors : list
        RGBA fill colours, one per set.
    """
    n_sets = len(names)
    if n_sets == 2:
        centers, radii = _layout_venn2(sizes)
        logics = ["10", "01", "11"]
    elif n_sets == 3:
        centers, radii = _layout_venn3(sizes)
        logics = ["100", "010", "001", "110", "101", "011", "111"]
    else:
        raise ValueError("Proportional Venn supports only 2 or 3 sets.")

    _draw_circles(ax, centers, radii, names, colors)

    xmin = min(cx - r for (cx, cy), r in zip(centers, radii))
    xmax = max(cx + r for (cx, cy), r in zip(centers, radii))
    ymin = min(cy - r for (cx, cy), r in zip(centers, radii))
    ymax = max(cy + r for (cx, cy), r in zip(centers, radii))
    # Guard against a degenerate (all-empty) bounding box.
    if xmax - xmin < 1e-6:
        xmin, xmax = -1.0, 1.0
    if ymax - ymin < 1e-6:
        ymin, ymax = -1.0, 1.0
    pad = 0.25 * max(xmax - xmin, ymax - ymin)

    positions = _region_label_positions(centers, radii, logics,
                                         (xmin, xmax, ymin, ymax))
    for logic in logics:
        if logic in positions and logic in labels:
            x, y = positions[logic]
            draw_text(ax, x, y, labels[logic], fontsize=fontsize)

    if legend_loc is not None:
        ax.legend(names, loc=legend_loc, prop={"size": fontsize})
    else:
        _place_names(ax, centers, radii, names, colors, fontsize,
                     legend_use_petal_color)

    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_aspect("equal")
    return ax
