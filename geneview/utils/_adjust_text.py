"""Automatic adjustment of text-label positions in matplotlib to minimize overlaps.

This is a vectorized re-implementation inspired by Phlya's ``adjustText``
(https://github.com/Phlya/adjustText) and
<https://stackoverflow.com/questions/19073683/matplotlib-overlapping-annotations-text>.

The original iterative core queried the renderer for the window extent of every
label on *every* iteration (``text.get_window_extent(renderer)``), which made it
prohibitively slow when many labels were present. A label's width and height in
display coordinates do **not** change when the label is merely translated, so
this version measures each label extent only *once* and then tracks label
positions analytically with NumPy arrays. The per-iteration cost therefore
becomes pure vectorized array math (O(N^2) broadcast operations) rather than
O(N) expensive renderer calls, which is typically 10-100x faster.

New capabilities compared with the previous vendored copy:

* ``only_move`` is honoured to restrict movement to a single axis per overlap
  type (e.g. keep a label locked to its point's x-position and only move it up).
* Repulsion uses a minimum-translation heuristic (push along the axis of least
  overlap) which produces much tidier layouts.
* Connecting arrows are drawn from the final label position back to the point.

Modified by Shujia Huang.
"""
from collections import deque

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.path import get_path_collection_extents


# --------------------------------------------------------------------------- #
# Renderer / bounding-box helpers (kept close to the original public helpers)
# --------------------------------------------------------------------------- #
def get_renderer(fig):
    try:
        return fig.canvas.get_renderer()
    except AttributeError:
        return fig.canvas.renderer


def get_bboxes_pathcollection(sc, ax):
    """Return a list of display-coordinate bounding boxes for a scatter plot.

    Thanks to ImportanceOfBeingErnest: https://stackoverflow.com/a/55007838
    """
    transform = sc.get_transform()
    transOffset = sc.get_offset_transform()
    offsets = sc._offsets
    paths = sc.get_paths()
    transforms = sc.get_transforms()

    if not transform.is_affine:
        paths = [transform.transform_path_non_affine(p) for p in paths]
        transform = transform.get_affine()
    if not transOffset.is_affine:
        offsets = transOffset.transform_non_affine(offsets)
        transOffset = transOffset.get_affine()

    if isinstance(offsets, np.ma.MaskedArray):
        offsets = offsets.filled(np.nan)

    bboxes = []
    if len(paths) and len(offsets):
        if len(paths) < len(offsets):
            paths = [paths[0]] * len(offsets)
        if len(transforms) < len(offsets):
            transforms = [transforms[0]] * len(offsets)
        for p, o, t in zip(paths, offsets, transforms):
            result = get_path_collection_extents(
                transform.frozen(), [p], [t], [o], transOffset.frozen()
            )
            bboxes.append(result)
    return bboxes


def get_bboxes(objs, r=None, expand=(1, 1), ax=None, transform=None):
    """Return a list of display-coordinate bounding boxes for ``objs``.

    Works with a list of matplotlib artists, a list of ``Bbox`` objects or a
    ``PathCollection`` (e.g. the return value of ``ax.scatter``).
    """
    ax = ax or plt.gca()
    r = r or get_renderer(ax.get_figure())
    try:
        return [i.get_window_extent(r).expanded(*expand) for i in objs]
    except (AttributeError, TypeError):
        try:
            if all(isinstance(obj, matplotlib.transforms.BboxBase) for obj in objs):
                return objs
            raise ValueError("Something is wrong")
        except TypeError:
            return get_bboxes_pathcollection(objs, ax)


def get_text_position(text, ax):
    """Return the display-coordinate position of a text object."""
    x, y = text.get_position()
    x = ax.convert_xunits(x)
    y = ax.convert_yunits(y)
    return text.get_transform().transform((x, y))


def set_text_position(text, t_x, t_y):
    """Set a text object's position from display coordinates."""
    x, y = text.get_transform().inverted().transform((t_x, t_y))
    text.set_position((x, y))


def get_orig_coords(transform, t_x, t_y):
    return transform.inverted().transform((t_x, t_y))


def get_midpoint(bbox):
    return (bbox.x0 + bbox.x1) / 2, (bbox.y0 + bbox.y1) / 2


def float_to_tuple(a):
    """Coerce a scalar or 2-sequence into a ``(x, y)`` float tuple."""
    try:
        a = float(a)
        return (a, a)
    except TypeError:
        assert len(a) == 2
        try:
            return float(a[0]), float(a[1])
        except TypeError:
            raise TypeError("Force values must be castable to floats")


# --------------------------------------------------------------------------- #
# Vectorized geometry core
#
# All coordinate arrays below use the layout ``coords[i] = [x0, x1, y0, y1]``
# in *display* coordinates, i.e. the left/right/bottom/top edges of box ``i``.
# --------------------------------------------------------------------------- #
def _coords_from_bboxes(bboxes):
    if not bboxes:
        return np.zeros((0, 4), dtype=float)
    return np.array([[b.xmin, b.xmax, b.ymin, b.ymax] for b in bboxes], dtype=float)


def _shift(base, disp):
    """Translate every box in ``base`` (N, 4) by ``disp`` (N, 2)."""
    out = base.copy()
    out[:, 0:2] += disp[:, 0:1]
    out[:, 2:4] += disp[:, 1:2]
    return out


def _expand_coords(coords, ex, ey):
    """Grow every box about its centre by factors ``(ex, ey)``."""
    if len(coords) == 0:
        return coords
    dx = (coords[:, 1] - coords[:, 0]) * (ex - 1.0) / 2.0
    dy = (coords[:, 3] - coords[:, 2]) * (ey - 1.0) / 2.0
    out = coords.copy()
    out[:, 0] -= dx
    out[:, 1] += dx
    out[:, 2] -= dy
    out[:, 3] += dy
    return out


def _centers(coords):
    return (coords[:, 0] + coords[:, 1]) / 2.0, (coords[:, 2] + coords[:, 3]) / 2.0


def _repel_text(coords):
    """Pairwise text-text repulsion using minimum-translation displacement.

    Returns per-text ``(dx, dy)`` displacement arrays plus the total ``x`` and
    ``y`` overlap magnitudes used as a convergence metric.
    """
    n = coords.shape[0]
    if n < 2:
        return np.zeros(n), np.zeros(n), 0.0, 0.0

    x0, x1, y0, y1 = coords[:, 0], coords[:, 1], coords[:, 2], coords[:, 3]
    ox = np.minimum(x1[:, None], x1[None, :]) - np.maximum(x0[:, None], x0[None, :])
    oy = np.minimum(y1[:, None], y1[None, :]) - np.maximum(y0[:, None], y0[None, :])
    overlap = (ox > 0) & (oy > 0)
    np.fill_diagonal(overlap, False)

    cx, cy = _centers(coords)
    sx = np.sign(cx[:, None] - cx[None, :])
    sy = np.sign(cy[:, None] - cy[None, :])
    # Break direction ties (labels sharing an exact centre) deterministically by
    # index order so perfectly-stacked labels still separate instead of drifting
    # together in the same direction.
    idx = np.arange(n)
    tie = np.sign(idx[:, None] - idx[None, :])
    sx = np.where(sx == 0, tie, sx)
    sy = np.where(sy == 0, tie, sy)
    # Push apart along the axis with the *smaller* overlap (minimum translation)
    # so labels separate cleanly instead of drifting diagonally.
    push_x = overlap & (ox <= oy)
    push_y = overlap & (oy < ox)
    mx = np.where(push_x, ox * sx, 0.0)
    my = np.where(push_y, oy * sy, 0.0)

    qx = np.where(overlap, ox, 0.0).sum() / 2.0
    qy = np.where(overlap, oy, 0.0).sum() / 2.0
    return mx.sum(axis=1), my.sum(axis=1), qx, qy


def _repel_from_boxes(coords, others):
    """Repel each box in ``coords`` away from the fixed boxes in ``others``."""
    n = coords.shape[0]
    if n == 0 or len(others) == 0:
        return np.zeros(n), np.zeros(n), 0.0, 0.0

    ox = (np.minimum(coords[:, 1][:, None], others[:, 1][None, :])
          - np.maximum(coords[:, 0][:, None], others[:, 0][None, :]))
    oy = (np.minimum(coords[:, 3][:, None], others[:, 3][None, :])
          - np.maximum(coords[:, 2][:, None], others[:, 2][None, :]))
    overlap = (ox > 0) & (oy > 0)

    cx, cy = _centers(coords)
    ocx, ocy = _centers(others)
    sx = np.sign(cx[:, None] - ocx[None, :])
    sy = np.sign(cy[:, None] - ocy[None, :])
    push_x = overlap & (ox <= oy)
    push_y = overlap & (oy < ox)
    mx = np.where(push_x, ox * np.where(sx == 0, 1.0, sx), 0.0)
    my = np.where(push_y, oy * np.where(sy == 0, 1.0, sy), 0.0)
    return mx.sum(axis=1), my.sum(axis=1), np.abs(mx).sum(), np.abs(my).sum()


def _repel_from_points(coords, px, py):
    """Repel each box in ``coords`` away from any point it contains."""
    n = coords.shape[0]
    if n == 0 or len(px) == 0:
        return np.zeros(n), np.zeros(n), 0.0, 0.0

    x0 = coords[:, 0][:, None]
    x1 = coords[:, 1][:, None]
    y0 = coords[:, 2][:, None]
    y1 = coords[:, 3][:, None]
    px = np.asarray(px)[None, :]
    py = np.asarray(py)[None, :]

    inside = (px > x0) & (px < x1) & (py > y0) & (py < y1)
    cx = (x0 + x1) / 2.0
    cy = (y0 + y1) / 2.0
    # Displacement needed to slide the box off the point along each axis.
    dx = np.where(cx >= px, px - x0, px - x1)
    dy = np.where(cy >= py, py - y0, py - y1)
    use_x = np.abs(dx) <= np.abs(dy)
    mx = np.where(inside & use_x, dx, 0.0)
    my = np.where(inside & ~use_x, dy, 0.0)
    return mx.sum(axis=1), my.sum(axis=1), np.abs(mx).sum(), np.abs(my).sum()


def _into_axes(coords, ax_bbox):
    """Displacement that brings each out-of-axes box back inside."""
    n = coords.shape[0]
    dx = np.zeros(n)
    dy = np.zeros(n)
    dx = np.where(coords[:, 0] < ax_bbox.xmin, ax_bbox.xmin - coords[:, 0], dx)
    dx = np.where(coords[:, 1] > ax_bbox.xmax, ax_bbox.xmax - coords[:, 1], dx)
    dy = np.where(coords[:, 2] < ax_bbox.ymin, ax_bbox.ymin - coords[:, 2], dy)
    dy = np.where(coords[:, 3] > ax_bbox.ymax, ax_bbox.ymax - coords[:, 3], dy)
    return dx, dy


def _mask_axes(dx, dy, spec):
    """Zero-out a displacement component when ``spec`` forbids that axis."""
    if "x" not in spec:
        dx = np.zeros_like(dx)
    if "y" not in spec:
        dy = np.zeros_like(dy)
    return dx, dy


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def adjust_text(
        texts,
        x=None,
        y=None,
        add_objects=None,
        ax=None,
        expand_text=(1.05, 1.2),
        expand_points=(1.05, 1.2),
        expand_objects=(1.05, 1.2),
        expand_align=(1.05, 1.2),
        autoalign="xy",
        va="center",
        ha="center",
        force_text=(0.2, 0.4),
        force_points=(0.2, 0.5),
        force_objects=(0.1, 0.25),
        lim=200,
        precision=0.01,
        only_move=None,
        avoid_text=True,
        avoid_points=True,
        avoid_self=True,
        min_arrow_len=5.0,
        *args,
        **kwargs):
    """Iteratively adjust the positions of ``texts`` to minimize overlaps.

    Call this *after* all plotting (especially anything that changes the axes
    limits) is complete, because the algorithm works in display coordinates and
    needs the final axes dimensions.

    Parameters
    ----------
    texts : list of matplotlib.text.Text
        The text objects to adjust.
    x, y : array_like, optional
        Coordinates of points to repel from (data coordinates). If not given,
        the original text positions are used when ``avoid_self`` is True.
    add_objects : list or PathCollection, optional
        Additional matplotlib objects (or Bboxes / a PathCollection) to avoid.
    ax : matplotlib.axes.Axes, optional
        Defaults to the current axes.
    expand_text, expand_points, expand_objects : (float, float)
        Multipliers used to pad the text bounding boxes when repelling from
        texts / points / objects respectively.
    force_text, force_points, force_objects : (float, float)
        Per-axis scaling of the repulsion displacement from each source.
    lim : int, default 200
        Maximum number of iterations.
    precision : float, default 0.01
        Convergence threshold, as a fraction of the summed label widths/heights.
    only_move : dict, optional
        Restrict movement per overlap type, e.g.
        ``{"points": "y", "text": "xy", "objects": "xy"}``. Valid values are
        ``""``, ``"x"``, ``"y"`` and ``"xy"``.
    avoid_text, avoid_points, avoid_self : bool
        Toggle repulsion between texts, from points, and from a label's own
        original position.
    min_arrow_len : float, default 5.0
        Minimum label displacement (in display units) below which no connecting
        arrow is drawn (avoids tiny stub arrows under the label).
    args, kwargs :
        ``arrowprops`` (and any extra kwargs) are forwarded to ``ax.annotate``
        to draw the connecting arrows once optimization is complete.

    Returns
    -------
    int
        The number of iterations performed.
    """
    if not len(texts):
        return 0

    if only_move is None:
        only_move = {"points": "xy", "text": "xy", "objects": "xy"}

    plt.draw()
    ax = ax or plt.gca()
    r = get_renderer(ax.get_figure())
    transform = texts[0].get_transform()

    force_text = float_to_tuple(force_text)
    force_points = float_to_tuple(force_points)
    force_objects = float_to_tuple(force_objects)

    for text in texts:
        text.set_va(va)
        text.set_ha(ha)

    # Measure every label once; widths/heights are invariant under translation.
    anchors = np.array([get_text_position(t, ax) for t in texts], dtype=float)
    base = _coords_from_bboxes(get_bboxes(texts, r, (1.0, 1.0), ax))
    disp = np.zeros((len(texts), 2), dtype=float)

    # Points to repel from.
    if x is not None and y is not None:
        pts = np.array([transform.transform((xi, yi)) for xi, yi in zip(x, y)],
                       dtype=float)
        px, py = pts[:, 0], pts[:, 1]
    elif avoid_self:
        px, py = anchors[:, 0].copy(), anchors[:, 1].copy()
    else:
        px, py = np.array([]), np.array([])

    # Objects to repel from.
    if add_objects is None:
        obj_coords = np.zeros((0, 4))
    else:
        try:
            obj_coords = _coords_from_bboxes(get_bboxes(add_objects, r, (1, 1), ax))
        except Exception:
            raise ValueError("Can't get bounding boxes from add_objects - is it a "
                             "flat list of matplotlib objects?")

    ax_bbox = ax.patch.get_extents()
    sum_width = float(np.sum(base[:, 1] - base[:, 0])) if len(base) else 1.0
    sum_height = float(np.sum(base[:, 3] - base[:, 2])) if len(base) else 1.0
    move_any_x = any("x" in v for v in only_move.values())
    move_any_y = any("y" in v for v in only_move.values())
    precision_x = precision * sum_width if move_any_x else np.inf
    precision_y = precision * sum_height if move_any_y else np.inf

    history = deque([(np.inf, np.inf)] * 5, maxlen=5)
    n_iter = 0
    for n_iter in range(1, lim + 1):
        coords = _shift(base, disp)
        step = np.zeros((len(texts), 2), dtype=float)
        qx = qy = 0.0

        if avoid_text:
            dx, dy, ox, oy = _repel_text(_expand_coords(coords, *expand_text))
            dx, dy = _mask_axes(dx, dy, only_move.get("text", "xy"))
            step[:, 0] += dx * force_text[0]
            step[:, 1] += dy * force_text[1]
            qx += ox
            qy += oy

        if avoid_points and len(px):
            dx, dy, ox, oy = _repel_from_points(_expand_coords(coords, *expand_points), px, py)
            dx, dy = _mask_axes(dx, dy, only_move.get("points", "xy"))
            step[:, 0] += dx * force_points[0]
            step[:, 1] += dy * force_points[1]
            qx += ox
            qy += oy

        if len(obj_coords):
            dx, dy, ox, oy = _repel_from_boxes(_expand_coords(coords, *expand_objects), obj_coords)
            dx, dy = _mask_axes(dx, dy, only_move.get("objects", "xy"))
            step[:, 0] += dx * force_objects[0]
            step[:, 1] += dy * force_objects[1]
            qx += ox
            qy += oy

        disp += step
        # Keep every label inside the axes.
        adx, ady = _into_axes(_shift(base, disp), ax_bbox)
        disp[:, 0] += adx
        disp[:, 1] += ady

        history.append((qx, qy))
        if qx < precision_x and qy < precision_y:
            break
        # Bail out if we are no longer making progress (oscillating / stuck).
        hist_max = np.max(np.array(history), axis=0)
        if qx >= hist_max[0] and qy >= hist_max[1] and n_iter > len(history):
            break

    # Commit the final positions back onto the text objects.
    final = anchors + disp
    for t, (fx, fy) in zip(texts, final):
        set_text_position(t, fx, fy)

    # Draw connecting arrows from each label back to its original point.
    if "arrowprops" in kwargs:
        kwap = kwargs.pop("arrowprops")
        coords = _shift(base, disp)
        for j, t in enumerate(texts):
            move = np.hypot(disp[j, 0], disp[j, 1])
            if move < min_arrow_len:
                continue
            ap = {"patchA": t}
            ap.update(kwap)
            mid = ((coords[j, 0] + coords[j, 1]) / 2.0,
                   (coords[j, 2] + coords[j, 3]) / 2.0)
            ax.annotate(
                "",
                xy=get_orig_coords(transform, anchors[j, 0], anchors[j, 1]),
                xytext=transform.inverted().transform(mid),
                arrowprops=ap,
                xycoords=transform,
                textcoords=transform,
                *args,
                **kwargs
            )

    return n_iter
