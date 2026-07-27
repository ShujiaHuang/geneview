"""
Shared gene feature / protein-domain drawing utilities.

Provides :func:`draw_features` which renders a horizontal baseline with
coloured rectangles representing gene features (exons, protein domains, etc.)
on a matplotlib Axes.  Used by both :class:`Lolliplot` and
:class:`DandelionPlot`.
"""

from typing import Optional, Dict, Any, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.patches as mpatches

from ._mutation_shapes import _is_dark, text_data_width


# Default colours for common feature types
_DEFAULT_FEATURE_FILLS = {
    "exon": "#51C6E6",
    "CDS": "#0080FF",
    "UTR": "#ADD8E6",
    "domain": "#FF8833",
    "motif": "#DFA32D",
    "default": "#B0C4DE",
}


def rescale_position(pos: float, rescale_map: Optional[List[Tuple]]) -> float:
    """Remap a genomic position through a coordinate mapping.

    Parameters
    ----------
    pos : float
        Original genomic position.
    rescale_map : list of (from_start, from_end, to_start, to_end) or None
        Coordinate remapping rules.  If *pos* falls within a ``from`` range
        it is linearly mapped to the corresponding ``to`` range.  Positions
        outside all ranges are returned unchanged.

    Returns
    -------
    float
        Remapped position.
    """
    if rescale_map is None:
        return pos
    for f_start, f_end, t_start, t_end in rescale_map:
        if f_start <= pos <= f_end:
            frac = (pos - f_start) / max(f_end - f_start, 1)
            return t_start + frac * (t_end - t_start)
    return pos


def draw_features(
    ax,
    feature_data: Optional[pd.DataFrame],
    region_start: int,
    region_end: int,
    baseline_y: float = 0.08,
    default_height: float = 0.05,
    show_labels: bool = True,
    label_fontsize: float = 7.0,
    label_on_feature: bool = False,
    rescale_map: Optional[List[Tuple]] = None,
    baseline_lw: float = 0.8,
    baseline_color: str = "black",
    label_color: str = "black",
) -> float:
    """Draw gene-feature rectangles on a horizontal baseline.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes (y-limits should already be set to [0, 1]).
    feature_data : pd.DataFrame or None
        DataFrame with columns ``start``, ``end``, and optionally ``name``,
        ``color``, ``fill``, ``height``, ``feature``, ``feature_layer_id``.
    region_start, region_end : int
        Genomic coordinate range (used for baseline extent).
    baseline_y : float
        Y position of the baseline in axes fraction [0, 1].
    default_height : float
        Default rectangle height in axes fraction.
    show_labels : bool
        Whether to draw feature name labels.
    label_fontsize : float
        Font size for feature labels.
    label_on_feature : bool
        If True, prefer drawing labels centred on the rectangle; otherwise
        labels are centred inside the rectangle when they fit and moved above
        it (externalised) when the rectangle is too narrow.
    rescale_map : list of (from_start, from_end, to_start, to_end), optional
        Coordinate remapping rules applied to feature start/end before drawing.
    baseline_lw : float
        Line width of the horizontal backbone / baseline.
    baseline_color : str
        Colour of the horizontal backbone / baseline.
    label_color : str
        Colour for externalised (above-rectangle) feature labels.  On-feature
        labels pick black/white automatically for contrast with the fill.

    Returns
    -------
    float
        Total vertical space consumed (in axes fraction) above the baseline,
        i.e. the top of the tallest feature rectangle.
    """
    if feature_data is None or len(feature_data) == 0:
        # Always draw the baseline even with no features
        ax.axhline(
            y=baseline_y, color=baseline_color, linewidth=baseline_lw,
            zorder=2, xmin=0.0, xmax=1.0,
        )
        return default_height + 0.01

    feats = feature_data.copy()
    feats.columns = [c.lower() for c in feats.columns]

    # Check for multi-layer support
    has_layers = "feature_layer_id" in feats.columns

    if has_layers:
        return _draw_multi_layer_features(
            ax, feats, region_start, region_end,
            baseline_y=baseline_y,
            default_height=default_height,
            show_labels=show_labels,
            label_fontsize=label_fontsize,
            label_on_feature=label_on_feature,
            rescale_map=rescale_map,
            baseline_lw=baseline_lw,
            baseline_color=baseline_color,
            label_color=label_color,
        )
    else:
        return _draw_single_layer_features(
            ax, feats, region_start, region_end,
            baseline_y=baseline_y,
            default_height=default_height,
            show_labels=show_labels,
            label_fontsize=label_fontsize,
            label_on_feature=label_on_feature,
            rescale_map=rescale_map,
            baseline_lw=baseline_lw,
            baseline_color=baseline_color,
            label_color=label_color,
        )


def _draw_single_layer_features(
    ax, feats, region_start, region_end,
    baseline_y, default_height, show_labels, label_fontsize,
    label_on_feature, rescale_map,
    baseline_lw=0.8, baseline_color="black", label_color="black",
) -> float:
    """Draw features on a single baseline."""
    # Always draw the baseline
    ax.axhline(
        y=baseline_y, color=baseline_color, linewidth=baseline_lw,
        zorder=2, xmin=0.0, xmax=1.0,
    )

    # Apply rescale
    if rescale_map is not None:
        feats = feats.copy()
        feats["start"] = feats["start"].apply(
            lambda p: rescale_position(p, rescale_map))
        feats["end"] = feats["end"].apply(
            lambda p: rescale_position(p, rescale_map))

    # Clip to visible region
    feats = feats[
        (feats["start"] < region_end) & (feats["end"] > region_start)
    ]
    if len(feats) == 0:
        return default_height + 0.01

    max_top = 0.0
    for _, row in feats.iterrows():
        h = float(row.get("height", default_height)) if "height" in feats.columns else default_height
        fill = (
            row.get("fill", _DEFAULT_FEATURE_FILLS.get(
                row.get("feature", "default"), _DEFAULT_FEATURE_FILLS["default"]))
            if "fill" in feats.columns
            else _DEFAULT_FEATURE_FILLS.get(
                row.get("feature", "default"), _DEFAULT_FEATURE_FILLS["default"])
        )
        edge_color = (
            row.get("color", "black") if "color" in feats.columns else "black"
        )
        lwd = float(row.get("lwd", 1)) if "lwd" in feats.columns else 1.0

        rect = mpatches.FancyBboxPatch(
            (row["start"], baseline_y - h / 2),
            row["end"] - row["start"],
            h,
            boxstyle="round,pad=0",
            facecolor=fill,
            edgecolor=edge_color,
            linewidth=lwd,
            zorder=3,
            clip_on=True,
        )
        ax.add_patch(rect)

        # Label
        name = row.get("name", "") if "name" in feats.columns else ""
        if show_labels and pd.notna(name) and str(name).strip():
            name_s = str(name)
            # A label is drawn inside the rectangle (with automatic black/white
            # contrast against the fill) whenever it fits; otherwise it is
            # externalised above the rectangle so it is never overplotted by
            # the block edge, the stems or the nodes above it.
            block_w = abs(row["end"] - row["start"])
            fits_inside = (
                text_data_width(ax, name_s, label_fontsize) <= block_w * 0.9
            )
            if fits_inside:
                ax.text(
                    (row["start"] + row["end"]) / 2,
                    baseline_y,
                    name_s,
                    ha="center",
                    va="center",
                    fontsize=label_fontsize,
                    fontweight="bold" if label_on_feature else "normal",
                    color="white" if _is_dark(fill) else "black",
                    zorder=4,
                    clip_on=True,
                )
            else:
                ax.text(
                    (row["start"] + row["end"]) / 2,
                    baseline_y + h / 2 + 0.005,
                    name_s,
                    ha="center",
                    va="bottom",
                    fontsize=label_fontsize,
                    color=label_color,
                    zorder=4,
                    clip_on=False,
                )

        top = h / 2
        if top > max_top:
            max_top = top

    return max_top + 0.01


def _draw_multi_layer_features(
    ax, feats, region_start, region_end,
    baseline_y, default_height, show_labels, label_fontsize,
    label_on_feature, rescale_map,
    baseline_lw=0.8, baseline_color="black", label_color="black",
) -> float:
    """Draw features grouped by ``feature_layer_id`` on separate baselines."""
    layers = feats.groupby("feature_layer_id", sort=False)
    layer_gap = 0.02  # vertical gap between layers
    current_y = baseline_y
    overall_max_top = 0.0

    for layer_id, layer_feats in layers:
        layer_top = _draw_single_layer_features(
            ax, layer_feats.drop(columns=["feature_layer_id"]),
            region_start, region_end,
            baseline_y=current_y,
            default_height=default_height,
            show_labels=show_labels,
            label_fontsize=label_fontsize,
            label_on_feature=label_on_feature,
            rescale_map=rescale_map,
            baseline_lw=baseline_lw,
            baseline_color=baseline_color,
            label_color=label_color,
        )
        overall_max_top = max(overall_max_top, current_y + layer_top - baseline_y)
        current_y += layer_top + layer_gap

    return overall_max_top
