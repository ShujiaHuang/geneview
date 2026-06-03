"""
AnnotationTrack - A track for displaying genomic range annotations.

Displays genomic intervals as rectangles, arrows, or ellipses, with support
for grouped features, strand direction, and multiple stacking modes.

Ported from Gviz's AnnotationTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union, Callable
from itertools import cycle

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from ._base import StackedTrack, GenomicInterval
from ._stacking import compute_stacking, get_stack_heights, get_num_stacks


# Default feature colors by feature type
_DEFAULT_FEATURE_COLORS = {
    "exon": "#3C5488",
    "CDS": "#3C5488",
    "UTR": "#8DB4E2",
    "five_prime_UTR": "#8DB4E2",
    "three_prime_UTR": "#8DB4E2",
    "intron": "#AAAAAA",
    "gene": "#5B8DB8",
    "transcript": "#5B8DB8",
    "default": "#5B8DB8",
}


class AnnotationTrack(StackedTrack):
    """A track for displaying genomic range annotations.

    Features are drawn as rectangles, arrows, or ellipses along the genomic
    coordinate axis. Overlapping features can be stacked using various modes
    (squish, pack, dense, full, hide).

    Parameters
    ----------
    data : pd.DataFrame or str
        DataFrame with columns: chrom, start, end, and optionally:
        strand, name, group, feature, id.
        If a string, interpreted as a file path (BED/GFF auto-detected).
    stacking : str
        Stacking mode: 'squish', 'pack', 'dense', 'full', or 'hide'.
        Default is 'squish'.
    feature_colors : dict, optional
        Mapping of feature types to colors. E.g. {"exon": "blue", "CDS": "red"}.
    shape : str
        Feature shape: 'box', 'ellipse', or 'arrow'. Default is 'box'.
    show_label : bool
        Whether to show feature name labels. Default is False.
    label_pos : str
        Label position: 'above', 'below', or 'inside'. Default is 'above'.
    arrow_direction : bool
        If True, show strand direction with arrows. Default is True.
    name : str
        Track name for the title panel. Default is "Annotation".
    height : float
        Relative track height. Default is 1.0.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> import pandas as pd
    >>> from geneview.genometracks import AnnotationTrack, plot_tracks, GenomicInterval
    >>> data = pd.DataFrame({
    ...     "chrom": ["chr7"] * 4,
    ...     "start": [2000000, 2070000, 2100000, 2160000],
    ...     "end":   [2050000, 2130000, 2150000, 2170000],
    ...     "strand": ["-", "+", "-", "-"],
    ...     "name": ["feat1", "feat2", "feat3", "feat4"],
    ...     "group": ["grp1", "grp2", "grp1", "grp3"],
    ... })
    >>> track = AnnotationTrack(data, stacking="squish")
    >>> plot_tracks([track], region=GenomicInterval("chr7", 1900000, 2200000))
    """

    SHAPES = ("box", "ellipse", "arrow")

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        stacking: str = "squish",
        feature_colors: Optional[Dict[str, str]] = None,
        shape: str = "box",
        show_label: bool = False,
        label_pos: str = "above",
        arrow_direction: bool = True,
        name: str = "Annotation",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if shape not in self.SHAPES:
            raise ValueError(f"Shape must be one of {self.SHAPES}, got '{shape}'.")

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=display_params,
        )

        self._feature_colors = dict(_DEFAULT_FEATURE_COLORS)
        if feature_colors:
            self._feature_colors.update(feature_colors)

        self.shape = shape
        self.show_label = show_label
        self.label_pos = label_pos
        self.arrow_direction = arrow_direction

    def _get_feature_color(self, feature: Optional[str], idx: int, color_cycle) -> str:
        """Get color for a feature, using feature_colors dict or cycling."""
        if feature is not None and feature in self._feature_colors:
            return self._feature_colors[feature]
        if feature is not None and pd.notna(feature):
            # Assign color from cycle for unknown features
            return next(color_cycle)
        return self.get_param("fill", "#5B8DB8")

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the annotation track.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        sub = self.subset(region)
        if sub is None or len(sub) == 0:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        # Ensure strand column
        if "strand" not in sub.columns:
            sub["strand"] = "*"

        # Compute stacking
        stacks = self.compute_stacks(sub)
        stack_info = get_stack_heights(stacks, mode=self.stacking, stack_height_frac=self.stack_height)

        total_rows = stack_info["total_rows"]
        row_height = stack_info["row_height"]
        y_positions = stack_info["y_positions"]

        # Set y limits based on stacking
        if total_rows == 0:
            ax.set_ylim(0, 1)
        else:
            ax.set_ylim(-0.05, 1.05)

        ax.set_xlim(region.start, region.end)

        # Get feature column for coloring
        has_feature = "feature" in sub.columns
        has_name = "name" in sub.columns or "id" in sub.columns
        label_col = "name" if "name" in sub.columns else ("id" if "id" in sub.columns else None)

        # Color cycle for unknown features
        default_colors = ["#3C5488", "#5B8DB8", "#DC0000", "#F8766D",
                          "#00BA38", "#619CFF", "#FDB462", "#B79F00"]
        color_cycle = cycle(default_colors)
        # Pre-assign colors to unique features
        feature_color_map = {}
        if has_feature:
            for feat in sub["feature"].dropna().unique():
                if feat not in self._feature_colors:
                    feature_color_map[feat] = next(color_cycle)

        span = region.end - region.start

        for i, (idx, row) in enumerate(sub.iterrows()):
            stack_row = stacks[i]
            y_center = y_positions[stack_row] if total_rows > 0 else 0.5
            feat_height = row_height * 0.85

            # Clip to visible region
            x_start = max(row["start"], region.start)
            x_end = min(row["end"], region.end)
            width = x_end - x_start

            if width <= 0:
                continue

            # Determine color
            if has_feature and pd.notna(row.get("feature")):
                feat = row["feature"]
                if feat in self._feature_colors:
                    color = self._feature_colors[feat]
                elif feat in feature_color_map:
                    color = feature_color_map[feat]
                else:
                    color = self.get_param("fill", "#5B8DB8")
            else:
                color = self.get_param("fill", "#5B8DB8")

            edge_color = self.get_param("col_border", "#333333")
            lwd = self.get_param("lwd", 0.5)
            alpha = self.get_param("alpha", 1.0)

            # Draw the feature shape
            if self.shape == "box":
                self._draw_box(ax, x_start, y_center, width, feat_height,
                               color, edge_color, lwd, alpha, row.get("strand", "*"))
            elif self.shape == "ellipse":
                self._draw_ellipse(ax, x_start, y_center, width, feat_height,
                                   color, edge_color, lwd, alpha)
            elif self.shape == "arrow":
                self._draw_arrow(ax, x_start, y_center, width, feat_height,
                                 color, edge_color, lwd, alpha, row.get("strand", "*"))

            # Draw strand direction indicator for boxes
            if self.shape == "box" and self.arrow_direction and row.get("strand", "*") in ("+", "-"):
                self._draw_strand_arrow(ax, x_start, x_end, y_center, feat_height,
                                        row["strand"], span)

            # Draw label if requested
            if self.show_label and label_col and pd.notna(row.get(label_col)):
                label = str(row[label_col])
                fontsize = self.get_param("fontsize", 8)
                if self.label_pos == "inside":
                    ax.text(x_start + width / 2, y_center, label,
                            ha="center", va="center", fontsize=fontsize,
                            color="white", fontweight="bold",
                            clip_on=True, zorder=5)
                elif self.label_pos == "below":
                    ax.text(x_start + width / 2, y_center - feat_height / 2 - 0.02,
                            label, ha="center", va="top", fontsize=fontsize,
                            color="#333333", clip_on=True, zorder=5)
                else:  # above
                    ax.text(x_start + width / 2, y_center + feat_height / 2 + 0.02,
                            label, ha="center", va="bottom", fontsize=fontsize,
                            color="#333333", clip_on=True, zorder=5)

        ax.axis("off")

    def _draw_box(self, ax, x, y, w, h, color, edge_color, lwd, alpha, strand="*"):
        """Draw a rectangular feature."""
        rect = FancyBboxPatch(
            (x, y - h / 2), w, h,
            boxstyle="round,pad=0",
            facecolor=color, edgecolor=edge_color,
            linewidth=lwd, alpha=alpha, zorder=3,
        )
        ax.add_patch(rect)

    def _draw_ellipse(self, ax, x, y, w, h, color, edge_color, lwd, alpha):
        """Draw an elliptical feature."""
        ellipse = mpatches.Ellipse(
            (x + w / 2, y), w, h,
            facecolor=color, edgecolor=edge_color,
            linewidth=lwd, alpha=alpha, zorder=3,
        )
        ax.add_patch(ellipse)

    def _draw_arrow(self, ax, x, y, w, h, color, edge_color, lwd, alpha, strand="*"):
        """Draw an arrow-shaped feature indicating strand direction."""
        arrow_width = min(w * 0.15, h * 2)

        if strand == "-":
            # Arrow pointing left
            verts = [
                (x + arrow_width, y - h / 2),
                (x + w, y - h / 2),
                (x + w, y + h / 2),
                (x + arrow_width, y + h / 2),
                (x, y),
                (x + arrow_width, y - h / 2),
            ]
        else:
            # Arrow pointing right (default, also for unstranded)
            verts = [
                (x, y - h / 2),
                (x + w - arrow_width, y - h / 2),
                (x + w, y),
                (x + w - arrow_width, y + h / 2),
                (x, y + h / 2),
                (x, y - h / 2),
            ]

        polygon = mpatches.Polygon(
            verts, closed=True,
            facecolor=color, edgecolor=edge_color,
            linewidth=lwd, alpha=alpha, zorder=3,
        )
        ax.add_patch(polygon)

    def _draw_strand_arrow(self, ax, x_start, x_end, y_center, feat_height, strand, span):
        """Draw a small strand direction indicator inside a box feature."""
        mid_x = (x_start + x_end) / 2
        arrow_len = min((x_end - x_start) * 0.3, span * 0.01)
        if arrow_len < span * 0.001:
            return  # Too small to draw

        if strand == "+":
            dx = arrow_len
        else:
            dx = -arrow_len

        ax.annotate(
            "",
            xy=(mid_x + dx / 2, y_center),
            xytext=(mid_x - dx / 2, y_center),
            arrowprops=dict(
                arrowstyle="->", color="white",
                lw=1.0, alpha=0.8,
            ),
            zorder=4,
        )
