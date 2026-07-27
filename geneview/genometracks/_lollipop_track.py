"""
LolliplotTrack — lollipop-style variant track for genometracks.

Ported from trackViewer's ``lolliplot()`` (R/Bioconductor).  Variants are
drawn as stems rising from a gene-feature baseline to shapes (circles, pies,
pins, flags) whose height encodes a score.  Works both as a composable track
inside :func:`plot_tracks` and as a standalone plot via :func:`lolliplot`.

Per-SNP columns supported (all optional):
    ``cex``, ``node_label``, ``node_label_color``, ``node_label_size``,
    ``label_rotation``, ``label_color``, ``dashline_col``, ``side``,
    ``shape``, ``stack_factor``.
"""

from typing import Optional, Dict, Any, List, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.patches as mpatches

from ._base import Track, GenomicInterval
from ._mutation_features import draw_features, rescale_position
from ._mutation_shapes import (
    draw_circle, draw_pie, draw_fan, draw_flag, draw_pin, draw_stem,
    draw_stem_zbend, draw_shape, draw_node_label, draw_square, draw_diamond,
    draw_triangle_up, draw_triangle_down,
    line_coupled_radius, _blended, _data_radius,
    resolve_style_visuals,
)


# ---------------------------------------------------------------------------
# Default colour palette for variant shapes
# ---------------------------------------------------------------------------
_DEFAULT_SNP_COLORS = [
    "#0080FF", "#E69F00", "#009E73", "#D55E00",
    "#CC79A7", "#56B4E9", "#F0E442",
]

SHAPE_NAMES = (
    "circle", "square", "diamond",
    "triangle_point_up", "triangle_point_down",
)


def _validate_snp_data(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise and validate a variant DataFrame."""
    df = df.copy()
    df.columns = [c.lower() for c in df.columns]
    required = {"start"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"snp_data missing required columns: {missing}")
    if "chrom" not in df.columns:
        df["chrom"] = "chr1"
    df["start"] = df["start"].astype(int)
    if "end" not in df.columns:
        df["end"] = df["start"] + 1
    else:
        df["end"] = df["end"].astype(int)
    if "score" not in df.columns:
        df["score"] = 1
    df["score"] = pd.to_numeric(df["score"], errors="coerce").fillna(1)
    return df


def _validate_features(df: Optional[pd.DataFrame]) -> Optional[pd.DataFrame]:
    if df is None:
        return None
    df = df.copy()
    df.columns = [c.lower() for c in df.columns]
    required = {"start", "end"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"features missing required columns: {missing}")
    if "chrom" not in df.columns:
        df["chrom"] = "chr1"
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def _snp_val(snp, col: str, default):
    """Safely read an optional per-row value with NaN fallback."""
    if col not in snp.index:
        return default
    val = snp[col]
    if isinstance(val, float) and np.isnan(val):
        return default
    if isinstance(val, str) and val.strip() == "":
        return default
    return val


def _is_list_val(val) -> bool:
    """Check if a cell value is a list (for multi-shape Tanghulu)."""
    return isinstance(val, (list, np.ndarray))


def _spread_x_positions(positions, min_gap: float, lo: float, hi: float):
    """Spread overlapping node positions apart horizontally.

    Mirrors trackViewer's ``lolliplot`` behaviour: nodes (lollipops) that
    sit closer than one node diameter are pushed apart so their shapes and
    labels remain legible, while distant nodes keep their true genomic
    coordinate.  Each group of mutually-overlapping nodes is re-centred on
    the group's original centroid and clamped to ``[lo, hi]``.

    Parameters
    ----------
    positions : array-like of float
        Original (draw) x-positions in data units.
    min_gap : float
        Minimum centre-to-centre spacing in data units (≈ node diameter).
    lo, hi : float
        Lower/upper bounds of the visible region (data units).

    Returns
    -------
    numpy.ndarray
        Adjusted x-positions, in the same order as the input.
    """
    positions = np.asarray(positions, dtype=float)
    n = len(positions)
    if n <= 1 or min_gap <= 0:
        return positions.copy()

    order = np.argsort(positions, kind="stable")
    xs = positions[order]
    out = xs.copy()

    # Group positions into clusters of mutually-overlapping nodes.
    clusters = [[0]]
    for i in range(1, n):
        if xs[i] - xs[i - 1] < min_gap:
            clusters[-1].append(i)
        else:
            clusters.append([i])

    for cl in clusters:
        k = len(cl)
        if k == 1:
            continue
        center = xs[cl].mean()
        total = (k - 1) * min_gap
        start = center - total / 2.0
        # Clamp the spread block within [lo, hi].
        if start < lo:
            start = lo
        if start + total > hi:
            start = hi - total
        for j, idx in enumerate(cl):
            out[idx] = start + j * min_gap

    result = np.empty_like(out)
    result[order] = out
    return result


# Physical size of one text line (inches) for a given ``cex`` — mirrors
# trackViewer's LINE-based sizing and matches ``line_coupled_radius``.
_LINE_PT = 11.0


def _label_inches(snp_df: pd.DataFrame, cex: float) -> float:
    """Vertical room (inches) needed by the longest 90°-rotated label."""
    labels: List[str] = []
    for col in ("label", "name"):
        if col in snp_df.columns:
            vals = [str(v) for v in snp_df[col].dropna() if str(v).strip()]
            if vals:
                labels = vals
                break
    max_chars = max((len(s) for s in labels), default=0)
    # Rotated text: its length becomes vertical extent (~0.6 em per char).
    return max_chars * (6.0 * cex) * 0.6 / 72.0


def _estimate_lolliplot_inches(snp_df: pd.DataFrame, plot_type: str,
                               cex: float, switch_limit: int) -> float:
    """Estimate the vertical content height (inches) of a lolliplot panel.

    Ported in spirit from trackViewer's ``getHeight()``: the panel needs to
    fit the gene feature, the stem/stack (score-dependent), and the labels.
    Because Stage-1 nodes are a fixed *physical* size (≈ one text line), a
    panel sized to this estimate lets the stack fill it instead of leaving a
    large empty band above the variants.
    """
    line_in = _LINE_PT * cex / 72.0  # node diameter / stacking pitch
    if len(snp_df) == 0:
        return 1.6
    scores = pd.to_numeric(snp_df["score"], errors="coerce").fillna(1.0)
    score_max = max(float(scores.max()), 1.0)
    score_max_int = int(np.ceil(score_max))
    score_type_int = all(float(s) == int(float(s)) for s in scores)
    use_tanghulu = (plot_type == "circle" and score_type_int
                    and score_max_int <= switch_limit)
    # Tanghulu stacks ``score_max`` nodes; other styles use a nominal stem.
    stack_in = score_max_int * line_in if use_tanghulu else 8.0 * line_in
    label_in = _label_inches(snp_df, cex)
    feature_in = 0.30
    gap_in = 0.18
    if "side" in snp_df.columns:  # caterpillar: content on both sides
        return feature_in + 2.0 * (stack_in + label_in) + gap_in
    return feature_in + stack_in + label_in + gap_in


# ---------------------------------------------------------------------------
# LolliplotTrack
# ---------------------------------------------------------------------------

class LolliplotTrack(Track):
    """Lollipop-style variant / mutation track.

    Draws variants as stems rising from a gene-feature baseline to shapes
    (circles, pies, pins, flags) at the top.  The stem height encodes a
    numeric *score* (e.g. mutation count, pathogenicity score).

    Can be used inside :func:`plot_tracks` alongside other tracks, or
    standalone via the :func:`lolliplot` convenience function.

    Parameters
    ----------
    snp_data : pd.DataFrame
        Variant data with columns: ``chrom``, ``start``.
        Optional columns: ``score``, ``color``, ``fill``, ``border``,
        ``label``, ``alpha``, ``shape``, ``side`` ('top' or 'bottom'),
        ``pie_values`` (list), ``pie_colors`` (list),
        ``cex``, ``node_label``, ``node_label_color``, ``node_label_size``,
        ``label_rotation``, ``label_color``, ``dashline_col``,
        ``stack_factor``.
    features : pd.DataFrame or None
        Gene-feature / domain data with columns: ``chrom``, ``start``,
        ``end``.  Optional: ``name``, ``color``, ``fill``, ``height``,
        ``feature``, ``feature_layer_id``.
    type : str
        Shape type: ``'circle'``, ``'pie'``, ``'pin'``, ``'flag'``,
        ``'pie.stack'``.
    cex : float
        Shape-size scaling factor (global default; overridden by per-SNP
        ``cex`` column).
    node_size : str or float
        Node radius control.  ``"auto"`` (default) couples the radius to the
        text line height (mirroring trackViewer's ``LINEW/2``) so nodes stay a
        sensible, panel-independent size.  A float is an explicit radius given
        as a fraction of the axes height (the historical default was ``0.04``)
        and is further scaled by ``cex``.
    ratio_yx : str or float
        Vertical-to-horizontal scaling of node stacking.  ``"auto"`` (default)
        keeps stacked nodes touching (round-circle pitch).  A float rescales
        the vertical pitch/offsets: ``>1`` adds vertical breathing room,
        ``<1`` packs stacked nodes tighter.
    adaptive_height : bool
        When ``True`` (default) and ``height`` is left unset, the track's
        relative height is derived from its content (feature + stem/stack +
        labels), mirroring trackViewer's ``getHeight()`` so the panel stays
        compact instead of leaving a large empty band above the variants.
        The standalone :func:`lolliplot` also auto-sizes the figure height
        from this estimate when ``figsize`` is not given.
    dashline_col : str
        Dotted guide-line colour (global default; overridden by per-SNP
        ``dashline_col`` column).
    show_yaxis : bool
        Whether to show a small y-axis indicating score scale.
    label_on_feature : bool
        Place feature labels on the rectangles.
    lollipop_style_switch_limit : int
        When max score exceeds this limit, switch from stacked-circle
        ("Tanghulu") style to single-circle style.
    legend : dict or None
        Legend specification: ``{"labels": [...], "fill": [...]}``.
        Optional keys: ``ncol``, ``loc``, ``fontsize``.
    jitter : str or None
        Label placement mode.  ``"auto"`` (default): automatically uses
        aligned placement when labels would overlap; ``"label"``: always
        use aligned labels with connectors; ``None``: inline labels
        directly above nodes (may overlap).
    yaxis : list or None
        Custom y-axis tick positions (e.g. ``[0, 5, 10]``).
    ylab : str or None
        Y-axis label text (e.g. ``"# evidences"``).
    rescale : list or None
        Coordinate remapping: list of ``(from_start, from_end, to_start,
        to_end)`` tuples.
    name : str
        Track name shown in the title panel.
    height : float
        Relative track height.
    """

    TYPES = ("circle", "pie", "pin", "flag", "pie.stack")

    def __init__(
        self,
        snp_data: pd.DataFrame,
        features: Optional[pd.DataFrame] = None,
        type: str = "circle",
        cex: float = 1.0,
        node_size: Union[str, float] = "auto",
        ratio_yx: Union[str, float] = "auto",
        dashline_col: str = "#CCCCCC",
        show_yaxis: bool = True,
        label_on_feature: bool = False,
        lollipop_style_switch_limit: int = 10,
        legend: Optional[Dict] = None,
        jitter: Optional[str] = "auto",
        yaxis: Optional[list] = None,
        ylab: Optional[str] = None,
        rescale: Optional[list] = None,
        name: str = "Lolliplot",
        adaptive_height: bool = True,
        height: Optional[float] = None,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if type not in self.TYPES:
            raise ValueError(f"type must be one of {self.TYPES}, got '{type}'")
        snp_df = _validate_snp_data(snp_data)
        # Resolve a content-aware relative height before initialising the base
        # Track (which stores ``height`` for the plot_tracks layout).
        if height is None:
            if adaptive_height:
                content_in = _estimate_lolliplot_inches(
                    snp_df, type, cex, lollipop_style_switch_limit)
                resolved_height = max(1.0, min(content_in / 1.2, 6.0))
            else:
                resolved_height = 3.0
        else:
            resolved_height = height
        super().__init__(name=name, height=resolved_height,
                         display_params=display_params)
        self.snp_data = snp_df
        self.features = _validate_features(features)
        self.type = type
        self.cex = cex
        self.node_size = node_size
        self.ratio_yx = ratio_yx
        self.adaptive_height = adaptive_height
        self.dashline_col = dashline_col
        self.show_yaxis = show_yaxis
        self.label_on_feature = label_on_feature
        self.lollipop_style_switch_limit = lollipop_style_switch_limit
        self.legend = legend
        self.jitter = jitter
        self.yaxis = yaxis
        self.ylab = ylab
        self.rescale = rescale

    # ------------------------------------------------------------------
    # Node-size helpers
    # ------------------------------------------------------------------

    def _base_radius_axes(self, ax, cex: Optional[float] = None) -> float:
        """Node radius (axes-height fraction) honouring ``node_size``/``cex``."""
        if cex is None:
            cex = self.cex
        return line_coupled_radius(ax, cex=cex, node_size=self.node_size)

    def _layout_ry(self, ax, radius_axes: float) -> float:
        """Vertical radius for stacking pitch/offsets, honouring ``ratio_yx``."""
        _, ry = _data_radius(ax, radius_axes)
        if self.ratio_yx != "auto":
            ry = ry * float(self.ratio_yx)
        return ry

    def _figure_height_inches(self) -> float:
        """Content-aware figure height (inches) for standalone rendering."""
        content_in = _estimate_lolliplot_inches(
            self.snp_data, self.type, self.cex,
            self.lollipop_style_switch_limit)
        # Add room for x-axis ticks and figure margins.
        return max(1.8, min(content_in + 0.6, 9.0))

    def _label_fontsize(self, ax, cex: Optional[float] = None) -> float:
        """Panel-aware label fontsize (pt) — avoids overwhelming small panels."""
        if cex is None:
            cex = self.cex
        base_pt = 6.0 * cex
        # Cap at ~60% of physical panel height to prevent giant labels
        fig = ax.figure
        pos = ax.get_position()
        h_inches = pos.height * fig.get_figheight()
        max_pt = h_inches * 72.0 * 0.08  # labels shouldn't exceed 8% of panel
        return min(base_pt, max(4.0, max_pt))

    # ------------------------------------------------------------------
    # Track interface
    # ------------------------------------------------------------------

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the full genomic extent of the data + features."""
        all_starts = []
        all_ends = []
        if len(self.snp_data) > 0:
            all_starts.append(int(self.snp_data["start"].min()))
            all_ends.append(int(self.snp_data["end"].max()))
        if self.features is not None and len(self.features) > 0:
            all_starts.append(int(self.features["start"].min()))
            all_ends.append(int(self.features["end"].max()))
        if not all_starts:
            return None
        chrom = str(self.snp_data["chrom"].iloc[0])
        return GenomicInterval(
            chrom=chrom,
            start=min(all_starts),
            end=max(all_ends),
        )

    def draw(self, ax, region: GenomicInterval) -> None:
        """Render the lollipop track onto the given Axes."""
        r_start = region.start
        r_end = region.end

        # Resolve style-driven visuals (line widths, colours, font sizes).
        # Falls back to the legacy constants when no plotstyle is active, so
        # un-styled output is unchanged.
        vis = resolve_style_visuals(getattr(self, "_gv_style", None))
        self._vis = vis

        # Fix the x-range up front so geometry-dependent helpers (node
        # de-overlap, circle radii) see the true region width instead of the
        # default (0, 1) autoscale limits that matplotlib would otherwise
        # apply only at render time.
        ax.set_xlim(r_start, r_end)
        ax.set_ylim(0, 1)
        ax.get_yaxis().set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        # ---- Detect caterpillar mode ----
        has_side = "side" in self.snp_data.columns
        if has_side:
            baseline_y = 0.45
        else:
            baseline_y = 0.08

        # ---- Draw features ----
        feature_height = draw_features(
            ax, self.features, r_start, r_end,
            baseline_y=baseline_y,
            label_on_feature=self.label_on_feature,
            rescale_map=self.rescale,
            baseline_lw=vis["stem_lw"],
            baseline_color=vis["edge_color"],
            label_fontsize=vis["label_fontsize"],
            label_color=vis["text_color"],
        )

        # ---- Filter SNPs to region ----
        snps = self.snp_data[
            (self.snp_data["start"] >= r_start)
            & (self.snp_data["start"] <= r_end)
        ].copy()

        if len(snps) == 0:
            if self.legend:
                self._draw_legend(ax)
            return

        snps = snps.sort_values("start").reset_index(drop=True)

        # ---- Apply rescale to positions ----
        if self.rescale is not None:
            snps["_draw_x"] = snps["start"].apply(
                lambda p: rescale_position(float(p), self.rescale))
        else:
            snps["_draw_x"] = snps["start"].astype(float)

        # ---- Compute layout parameters ----
        score_max = max(float(snps["score"].max()), 1.0)
        score_max_int = int(np.ceil(score_max))
        score_type_int = all(
            float(s) == int(float(s)) for s in snps["score"]
        )

        use_tanghulu = (
            self.type == "circle"
            and score_type_int
            and score_max_int <= self.lollipop_style_switch_limit
        )
        if use_tanghulu:
            score_max_display = float(score_max_int)
        else:
            score_max_display = score_max

        feature_top = baseline_y + feature_height
        stem_base_top = feature_top + 0.005
        stem_top_max = 0.82
        radius_axes_base = self._base_radius_axes(ax)
        rx_base, _ = _data_radius(ax, radius_axes_base)
        ry_base = self._layout_ry(ax, radius_axes_base)
        v_range_top = stem_top_max - stem_base_top

        # ---- Spread overlapping node x-positions (trackViewer-style) ----
        # Nodes closer than one diameter are pushed apart so shapes and
        # labels stay legible; distant nodes keep their true coordinate.
        # ``pie.stack`` intentionally stacks multiple pies at one position,
        # so spreading is skipped there.
        if self.type != "pie.stack":
            min_gap = 2.0 * rx_base * 1.15
            snps["_node_x"] = _spread_x_positions(
                snps["_draw_x"].values, min_gap, r_start, r_end,
            )
        else:
            snps["_node_x"] = snps["_draw_x"].astype(float)

        # Bottom layout for caterpillar
        if has_side:
            bot_feature_bottom = baseline_y - feature_height
            stem_base_bot = bot_feature_bottom - 0.005
            stem_bot_max = 0.02
            v_range_bot = stem_base_bot - stem_bot_max

        # ---- pie.stack mode ----
        if self.type == "pie.stack":
            self._draw_pie_stack(ax, snps, stem_base_top, stem_top_max,
                                 radius_axes_base, feature_top)
            if self.legend:
                self._draw_legend(ax)
            return

        # ---- Decide label mode (auto-jitter) ----
        use_aligned = (self.jitter == "label")
        if self.jitter == "auto" and self.type != "flag":
            # Check if labeled nodes are densely packed enough to warrant
            # aligned labels.  Criterion: any two adjacent labelled nodes
            # are within 3% of the x-range — dense enough that rotated
            # labels form a cluttered wall of text.
            has_labels = any(
                ("label" in snps.columns and pd.notna(snps["label"].iloc[i])
                 and str(snps["label"].iloc[i]).strip())
                or ("name" in snps.columns and pd.notna(snps["name"].iloc[i])
                    and str(snps["name"].iloc[i]).strip())
                for i in range(len(snps))
            )
            if has_labels and len(snps) > 1:
                node_xs = snps["_node_x"].values
                sorted_xs = np.sort(node_xs)
                min_dist = float(np.min(np.diff(sorted_xs)))
                x_range = float(r_end - r_start) or 1.0
                # Threshold: 3% of region — tight enough for labels to clash
                if min_dist < x_range * 0.03:
                    use_aligned = True

        # ---- Draw lollipops (first pass: shapes + stems) ----
        label_items = []  # collect for aligned-label post-processing
        for snp_idx, (_, snp) in enumerate(snps.iterrows()):
            x = float(snp["_node_x"])
            true_x = float(snp["_draw_x"])
            raw_score = int(float(snp["score"]))

            # Per-SNP properties
            snp_cex = float(_snp_val(snp, "cex", self.cex))
            fill = _snp_val(snp, "fill", "#0080FF")
            border = _snp_val(snp, "border", "black")
            alpha_val = float(_snp_val(snp, "alpha", 1.0))
            snp_dashline = _snp_val(snp, "dashline_col", self.dashline_col)
            side = _snp_val(snp, "side", "top")
            shape_name = _snp_val(snp, "shape", "circle")

            # Label properties
            label_rot = float(_snp_val(snp, "label_rotation", 90))
            label_col = _snp_val(snp, "label_color", "black")

            # Node label
            node_lbl = _snp_val(snp, "node_label", None)
            node_lbl_col = _snp_val(snp, "node_label_color", "white")
            node_lbl_size = float(_snp_val(snp, "node_label_size",
                                           6 * snp_cex))

            radius_axes = self._base_radius_axes(ax, snp_cex)
            ry = self._layout_ry(ax, radius_axes)

            # Determine layout direction
            if side == "bottom" and has_side:
                base_y = stem_base_bot
                top_y = stem_bot_max
                v_range = v_range_bot
                direction = -1
                feat_top_draw = bot_feature_bottom
            else:
                base_y = stem_base_top
                top_y = stem_top_max
                v_range = v_range_top
                direction = 1
                feat_top_draw = feature_top

            # Jitter (unused — deferred to _draw_aligned_labels)

            if use_tanghulu:
                label_y = self._draw_tanghulu(
                    ax, snp, x, raw_score, base_y, ry, radius_axes,
                    fill, border, alpha_val, snp_dashline,
                    score_max_display, top_y, direction,
                    feat_top_draw, shape_name,
                    node_lbl, node_lbl_col, node_lbl_size, true_x,
                )
            else:
                label_y = self._draw_non_tanghulu(
                    ax, snp, x, base_y, top_y, v_range,
                    radius_axes, ry, fill, border, alpha_val,
                    snp_dashline, direction, feat_top_draw,
                    shape_name, score_max_display,
                    node_lbl, node_lbl_col, node_lbl_size, true_x,
                )

            # External label
            label = _snp_val(snp, "label", _snp_val(snp, "name", ""))
            if self.type != "flag" and label and str(label).strip():
                lbl_fs = self._label_fontsize(ax, snp_cex)
                if use_aligned and direction == 1:
                    # Defer label drawing for aligned placement
                    label_items.append({
                        "x": x,
                        "shape_top_y": label_y,
                        "text": str(label),
                        "rotation": label_rot,
                        "color": label_col,
                        "fontsize": lbl_fs,
                        "cex": snp_cex,
                        "dashline_col": snp_dashline,
                    })
                else:
                    # Top-side labels grow upward from their anchor; bottom-side
                    # (caterpillar) labels must grow downward so they sit clear
                    # of the nodes instead of overlapping them.
                    label_va = "bottom" if direction == 1 else "top"
                    ax.text(
                        x, label_y, str(label),
                        ha="center", va=label_va, rotation=label_rot,
                        fontsize=lbl_fs, color=label_col,
                        zorder=6, clip_on=False,
                    )

        # ---- Aligned labels (jitter="label") ----
        if label_items:
            self._draw_aligned_labels(
                ax, label_items, stem_top_max, ry_base,
            )

        # ---- Y-axis ----
        if self.show_yaxis and score_max > 1:
            if use_tanghulu:
                y_top = min(stem_base_top + score_max_display * 2 * ry_base,
                            stem_top_max)
            else:
                y_top = stem_top_max
            self._draw_yaxis(ax, score_max, stem_base_top, y_top)

        # ---- Legend ----
        if self.legend:
            self._draw_legend(ax)

    # ------------------------------------------------------------------
    # Aligned labels (trackViewer-style jitter="label")
    # ------------------------------------------------------------------

    def _draw_aligned_labels(self, ax, label_items, stem_top_max, ry_base):
        """Draw labels aligned at a uniform y-level with anti-overlap.

        Mirrors trackViewer's behaviour: all labels are placed near the
        top of the plot area; thin dashed connector lines link each shape
        to its label.  When two adjacent labels are too close
        horizontally, the later one is shifted upward to avoid overlap.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        label_items : list[dict]
            Each dict has keys: x, shape_top_y, text, rotation, color,
            fontsize, cex, dashline_col.
        stem_top_max : float
            Upper y-limit of the stem region (axes fraction).
        ry_base : float
            Base circle radius in data-y units (for connector spacing).
        """
        if not label_items:
            return

        fig = ax.figure
        pos = ax.get_position()
        xlim = ax.get_xlim()
        x_range = xlim[1] - xlim[0]
        w_inches = pos.width * fig.get_figwidth()
        h_inches = pos.height * fig.get_figheight()

        # Base y for the first row of labels (above all shapes)
        label_base_y = stem_top_max + 0.04
        row_step = 0.035  # vertical step between stagger rows

        def _label_width_data(text, fontsize, rotation):
            """Estimate the horizontal footprint of a label in data units."""
            char_w_px = fontsize * 0.55
            text_w_px = max(len(str(text)), 1) * char_w_px
            text_h_px = fontsize * 1.1
            rot_rad = np.radians(abs(rotation))
            effective_px = (
                abs(text_w_px * np.sin(rot_rad))
                + abs(text_h_px * np.cos(rot_rad))
            )
            effective_px = max(effective_px, char_w_px * 2)
            return effective_px / (w_inches * fig.dpi) * x_range

        # Sort by x (they should already be sorted)
        label_items.sort(key=lambda d: d["x"])

        # Greedy row assignment: each label starts at row 0
        n = len(label_items)
        rows = [0] * n
        placed = []  # list of (x, width_data, row) for overlap checks

        for i, item in enumerate(label_items):
            w = _label_width_data(item["text"], item["fontsize"],
                                  item["rotation"])
            half_w = w * 0.55  # slight padding
            best_row = 0
            for px, pw, pr in placed:
                if abs(item["x"] - px) < (half_w + pw * 0.55) and pr == best_row:
                    best_row = pr + 1
                    # Re-check with new row
                    for px2, pw2, pr2 in placed:
                        if (abs(item["x"] - px2) < (half_w + pw2 * 0.55)
                                and pr2 == best_row):
                            best_row = pr2 + 1
            rows[i] = best_row
            placed.append((item["x"], w, best_row))

        # Draw each label with connector line
        for i, item in enumerate(label_items):
            x = item["x"]
            shape_top_y = item["shape_top_y"]
            label_y = label_base_y + rows[i] * row_step

            # Clamp to plot area
            label_y = min(label_y, 0.97)

            # Connector line from shape top to label
            if label_y > shape_top_y + ry_base * 0.5:
                draw_stem(
                    ax, x, shape_top_y, label_y,
                    color=item["dashline_col"],
                    linewidth=self._vis["connector_lw"], linestyle=":",
                )

            # Draw label text
            ax.text(
                x, label_y, item["text"],
                ha="center", va="bottom",
                rotation=item["rotation"],
                fontsize=item["fontsize"],
                color=item["color"],
                zorder=6, clip_on=False,
            )

    # ------------------------------------------------------------------
    # Tanghulu drawing
    # ------------------------------------------------------------------

    def _draw_tanghulu(
        self, ax, snp, x, raw_score, base_y, ry, radius_axes,
        fill, border, alpha_val, dashline_col,
        score_max_display, top_y, direction, feat_top,
        shape_name, node_lbl, node_lbl_col, node_lbl_size, true_x=None,
    ) -> float:
        """Draw a Tanghulu (stacked circles) lollipop. Returns label_y."""
        if true_x is None:
            true_x = x
        # Handle list-valued shape/fill
        shapes_list = _snp_val(snp, "shape", None)
        fills_list = _snp_val(snp, "fill", None)
        use_multi_shape = _is_list_val(shapes_list)
        use_multi_fill = _is_list_val(fills_list)

        # Draw solid stem from feature top to first circle, rooted at the
        # variant's true coordinate and deflected to the (possibly displaced)
        # node column so the stem stays connected to the feature below.
        stem_col = border if not isinstance(border, list) else border[0]
        if direction > 0:
            first_cy = base_y + ry
        else:
            first_cy = base_y - ry
        draw_stem_zbend(ax, true_x, x, feat_top, first_cy,
                        color=stem_col, linewidth=self._vis["stem_lw"],
                        linestyle="-")

        for i in range(1, raw_score + 1):
            if direction > 0:
                cy = base_y + (i - 0.5) * 2 * ry
            else:
                cy = base_y - (i - 0.5) * 2 * ry

            # Per-circle shape and fill
            circle_shape = shape_name
            circle_fill = fill
            if use_multi_shape:
                circle_shape = shapes_list[(i - 1) % len(shapes_list)]
            if use_multi_fill:
                circle_fill = fills_list[(i - 1) % len(fills_list)]

            draw_shape(
                ax, circle_shape, x, cy, radius_axes,
                facecolor=circle_fill, edgecolor=border,
                alpha=alpha_val, linewidth=self._vis["border_lw"],
            )

            # Node label on each circle
            if node_lbl is not None:
                draw_node_label(ax, x, cy, radius_axes,
                                str(node_lbl), node_lbl_col,
                                node_lbl_size)

        if direction > 0:
            top_circle_y = base_y + (raw_score - 0.5) * 2 * ry
            guide_top = base_y + score_max_display * 2 * ry
            guide_top = min(guide_top, top_y)
            if guide_top > top_circle_y + ry:
                draw_stem(ax, x, top_circle_y + ry, guide_top,
                          color=dashline_col, linewidth=self._vis["guide_lw"],
                          linestyle=":")
            return top_circle_y + ry + 0.005
        else:
            top_circle_y = base_y - (raw_score - 0.5) * 2 * ry
            return top_circle_y - ry - 0.005

    # ------------------------------------------------------------------
    # Non-Tanghulu drawing
    # ------------------------------------------------------------------

    def _draw_non_tanghulu(
        self, ax, snp, x, base_y, top_y, v_range,
        radius_axes, ry, fill, border, alpha_val,
        dashline_col, direction, feat_top,
        shape_name, score_max_display,
        node_lbl, node_lbl_col, node_lbl_size, true_x=None,
    ) -> float:
        """Draw a non-Tanghulu lollipop. Returns label_y."""
        if true_x is None:
            true_x = x
        score = float(snp["score"])
        if score_max_display > 0:
            frac = score / score_max_display
        else:
            frac = 1.0

        if direction > 0:
            shape_y = base_y + frac * v_range
        else:
            shape_y = base_y - frac * v_range

        draw_stem_zbend(ax, true_x, x, feat_top, shape_y,
                        color=border, linewidth=self._vis["stem_lw"],
                        linestyle="-")

        if self.type == "circle" or self.type == "pie.stack":
            draw_shape(
                ax, shape_name, x, shape_y, radius_axes,
                facecolor=fill, edgecolor=border, alpha=alpha_val,
                linewidth=self._vis["border_lw"],
            )
        elif self.type == "pie":
            pie_vals = _snp_val(snp, "pie_values", [1, 1])
            pie_cols = _snp_val(snp, "pie_colors", ["#0080FF", "#E69F00"])
            if isinstance(pie_vals, str):
                pie_vals = [float(v) for v in pie_vals.split(",")]
            if isinstance(pie_cols, str):
                pie_cols = pie_cols.split(",")
            draw_pie(ax, x, shape_y, radius_axes, pie_vals, pie_cols,
                     edgecolor=border, linewidth=self._vis["border_lw"])
        elif self.type == "pin":
            draw_pin(ax, x, shape_y, radius_axes,
                     facecolor=fill if fill else "#D55E00",
                     edgecolor=border, alpha=alpha_val,
                     linewidth=self._vis["border_lw"])
        elif self.type == "flag":
            label_text = _snp_val(snp, "label", _snp_val(snp, "name", ""))
            if pd.isna(label_text):
                label_text = ""
            draw_flag(ax, x, shape_y, radius_axes,
                      color=fill, edgecolor=border,
                      cex=float(_snp_val(snp, "cex", self.cex)),
                      label=str(label_text),
                      linewidth=self._vis["border_lw"])

        # Node label inside shape
        if node_lbl is not None and self.type != "flag":
            draw_node_label(ax, x, shape_y, radius_axes,
                            str(node_lbl), node_lbl_col, node_lbl_size)

        # Guide line for circle type
        if self.type == "circle" and score_max_display > 1:
            if direction > 0:
                guide_top = top_y
                if shape_y + ry < guide_top:
                    draw_stem(ax, x, shape_y + ry, guide_top,
                              color=dashline_col,
                              linewidth=self._vis["guide_lw"], linestyle=":")
            # (no guide for bottom direction)

        if direction > 0:
            label_y = shape_y + ry + 0.005
            if self.type == "pin":
                label_y = shape_y + ry * 3.2
            return label_y
        else:
            label_y = shape_y - ry - 0.005
            if self.type == "pin":
                label_y = shape_y - ry * 3.2
            return label_y

    # ------------------------------------------------------------------
    # pie.stack drawing
    # ------------------------------------------------------------------

    def _draw_pie_stack(self, ax, snps, stem_base, stem_top_max,
                        radius_axes_base, feature_top):
        """Draw pie.stack layout: stacked pies at each position."""
        has_stack_factor = "stack_factor" in snps.columns

        # Group SNPs by position
        grouped = snps.groupby("_draw_x")
        max_stacks = 1
        if has_stack_factor:
            for _, group in grouped:
                max_stacks = max(max_stacks, len(group))

        score_max = max(float(snps["score"].max()), 1.0)
        v_range = stem_top_max - stem_base

        for x_pos, group in grouped:
            x = float(x_pos)
            n_stacks = len(group)

            # Stem height = number of stacks
            norm_h = n_stacks / max(max_stacks, 1)
            stem_top_y = stem_base + norm_h * v_range * 0.8

            draw_stem(ax, x, feature_top, stem_top_y,
                      color=self._vis["edge_color"],
                      linewidth=self._vis["stem_lw"], linestyle="-")

            # Draw stacked pies
            pie_r = radius_axes_base * 0.8
            _, ry_pie = _data_radius(ax, pie_r)

            for si, (_, snp) in enumerate(group.iterrows()):
                cy = stem_base + (si + 0.5) * 2 * ry_pie
                pie_vals = _snp_val(snp, "pie_values", [1, 1])
                pie_cols = _snp_val(snp, "pie_colors", ["#0080FF", "#E69F00"])
                if isinstance(pie_vals, str):
                    pie_vals = [float(v) for v in pie_vals.split(",")]
                if isinstance(pie_cols, str):
                    pie_cols = pie_cols.split(",")
                draw_pie(ax, x, cy, pie_r, pie_vals, pie_cols,
                         linewidth=self._vis["border_lw"])

                # Label
                label = _snp_val(snp, "label", _snp_val(snp, "name", ""))
                if label and str(label).strip():
                    ax.text(x, cy + ry_pie + 0.003, str(label),
                            ha="center", va="bottom", rotation=90,
                            fontsize=5, zorder=6, clip_on=False)

            # Guide line
            top_pie_y = stem_base + (n_stacks - 0.5) * 2 * ry_pie
            if top_pie_y + ry_pie < stem_top_max:
                draw_stem(ax, x, top_pie_y + ry_pie, stem_top_max,
                          color=self._vis["guide_color"],
                          linewidth=self._vis["guide_lw"], linestyle=":")

    # ------------------------------------------------------------------
    # Y-axis
    # ------------------------------------------------------------------

    def _draw_yaxis(self, ax, score_max, y_bottom, y_top):
        """Draw a small y-axis on the left indicating score scale."""
        vis = getattr(self, "_vis", None) or resolve_style_visuals(None)
        ax.annotate(
            "", xy=(0, y_bottom), xytext=(0, y_top),
            xycoords=ax.get_yaxis_transform(),
            arrowprops=dict(arrowstyle="-", color=vis["axis_color"],
                            lw=vis["axis_lw"]),
        )

        if self.yaxis is not None:
            # Custom tick positions
            for tick_val in self.yaxis:
                frac = tick_val / max(score_max, 1)
                tick_y = y_bottom + frac * (y_top - y_bottom)
                label = str(tick_val)
                if isinstance(self.yaxis, dict):
                    label = self.yaxis.get(tick_val, str(tick_val))
                ax.text(-0.02, tick_y, label,
                        transform=ax.transAxes,
                        ha="right", va="center",
                        fontsize=vis["tick_fontsize"], color=vis["axis_color"])
        else:
            ax.text(-0.02, y_bottom, "0",
                    transform=ax.transAxes,
                    ha="right", va="center",
                    fontsize=vis["tick_fontsize"], color=vis["axis_color"])
            ax.text(-0.02, y_top, str(int(score_max)),
                    transform=ax.transAxes,
                    ha="right", va="center",
                    fontsize=vis["tick_fontsize"], color=vis["axis_color"])

        # ylab
        if self.ylab:
            mid_y = (y_bottom + y_top) / 2
            ax.text(-0.06, mid_y, self.ylab,
                    transform=ax.transAxes,
                    ha="center", va="center", rotation=90,
                    fontsize=vis["label_fontsize"], color=vis["text_color"])

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------

    def _draw_legend(self, ax):
        """Draw a legend based on the ``self.legend`` dict."""
        if not self.legend:
            return
        labels = self.legend.get("labels", [])
        fills = self.legend.get("fill", [])
        if not labels or not fills:
            return

        handles = []
        for fill_col in fills:
            handles.append(mpatches.Patch(
                facecolor=fill_col,
                edgecolor=self._vis["edge_color"] if getattr(self, "_vis", None)
                else "gray",
                linewidth=self._vis["connector_lw"] if getattr(self, "_vis", None)
                else 0.5))

        ncol = self.legend.get("ncol", max(1, len(labels) // 4 + 1))
        loc = self.legend.get("loc", "upper right")
        fontsize = self.legend.get("fontsize", 7)

        ax.legend(handles, labels, loc=loc, ncol=ncol,
                  fontsize=fontsize, framealpha=0.8,
                  bbox_to_anchor=(1.0, 1.0))

    # ------------------------------------------------------------------

    def __repr__(self):
        return (
            f"LolliplotTrack(name='{self.name}', height={self.height}, "
            f"type='{self.type}', n_variants={len(self.snp_data)})"
        )


# ---------------------------------------------------------------------------
# Standalone convenience function
# ---------------------------------------------------------------------------

def lolliplot(
    snp_data: pd.DataFrame,
    features: Optional[pd.DataFrame] = None,
    type: str = "circle",
    region: Optional[GenomicInterval] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    show_title: bool = False,
    legend: Optional[Dict] = None,
    style: Optional[str] = None,
    ax=None,
    **kwargs,
):
    """Convenience function: create a :class:`LolliplotTrack` and render it.

    When called without ``ax``, creates a new figure via
    :func:`plot_tracks`.  When ``ax`` is provided, draws directly
    onto it (standalone mode).

    Parameters
    ----------
    snp_data : pd.DataFrame
        Variant data (see :class:`LolliplotTrack`).
    features : pd.DataFrame or None
        Gene-feature data.
    type : str
        Shape type: ``'circle'``, ``'pie'``, ``'pin'``, ``'flag'``,
        ``'pie.stack'``.
    region : GenomicInterval, optional
        Genomic region.  Auto-determined from data if *None*.
    figsize : tuple, optional
        Figure size in inches.  When *None* (default), width is 12"
        and height is auto-computed from the content.
    title : str, optional
        Plot title.
    show_title : bool
        Whether to show the track title panel.  Default ``False``.
    legend : dict, optional
        Legend specification passed to :class:`LolliplotTrack`.
    style : str, optional
        Name of a registered plot style (e.g. ``"nature"``, ``"science"``,
        ``"cell"``).  Routes node borders, stems, connector/guide lines and
        axis rules through the journal style.
    ax : matplotlib Axes, optional
        Draw into this axes instead of creating a new figure.
    **kwargs
        Additional keyword arguments passed to :class:`LolliplotTrack`.

    Returns
    -------
    matplotlib.axes.Axes or list of Axes
    """
    track = LolliplotTrack(
        snp_data, features=features, type=type, legend=legend, **kwargs,
    )

    # Resolve the requested style so standalone rendering can route its
    # visuals through the journal style too.
    resolved_style = None
    if style is not None:
        from ..plotstyle import get_style as _get_style
        resolved_style = _get_style(style)
        track._gv_style = resolved_style

    # Draw into an existing axes (fully standalone)
    if ax is not None:
        if region is None:
            region = track.get_region()
            if region is None:
                raise ValueError("Cannot determine region from data.")
            span = region.end - region.start
            region = GenomicInterval(
                region.chrom,
                max(0, region.start - int(span * 0.02)),
                region.end + int(span * 0.02),
            )
        ax.set_xlim(region.start, region.end)
        track.draw(ax, region)
        return ax

    # Use plot_tracks for a polished layout
    from ._track_plot import plot_tracks

    if region is None:
        region = track.get_region()
        if region is not None:
            span = region.end - region.start
            region = GenomicInterval(
                region.chrom,
                max(0, region.start - int(span * 0.02)),
                region.end + int(span * 0.02),
            )

    axes = plot_tracks(
        [track],
        region=region,
        figsize=figsize or (12, track._figure_height_inches()),
        title=title,
        show_title=show_title,
        style=style,
    )
    return axes[0] if isinstance(axes, list) else axes
