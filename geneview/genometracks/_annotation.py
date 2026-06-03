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


# Default feature colors by feature type (Gviz palette)
_DEFAULT_FEATURE_COLORS = {
    "exon": "lightblue",
    "CDS": "#0080FF",
    "UTR": "#ADD8E6",
    "five_prime_UTR": "#ADD8E6",
    "three_prime_UTR": "#ADD8E6",
    "intron": "#808080",
    "gene": "lightblue",
    "transcript": "lightblue",
    "default": "lightblue",
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
        Feature shape: 'box', 'ellipse', 'arrow', 'fixedArrow', or 'smallArrow'.
        Default is 'box'.
    show_label : bool
        Whether to show feature name labels. Default is False.
    label_pos : str
        Label position: 'above', 'below', or 'inside'. Default is 'above'.
    arrow_direction : bool
        If True, show strand direction with arrows. Default is True.
    group_annotation : str or None
        Label groups by: 'id', 'group', 'feature', or None.
        When set, draw connecting lines and group labels. Default is None.
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

    SHAPES = ("box", "ellipse", "arrow", "fixedArrow", "smallArrow")

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        stacking: str = "squish",
        feature_colors: Optional[Dict[str, str]] = None,
        shape: str = "arrow",
        show_label: bool = False,
        label_pos: str = "above",
        arrow_direction: bool = True,
        group_annotation: Optional[str] = None,
        name: str = "Annotation",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if shape not in self.SHAPES:
            raise ValueError(f"Shape must be one of {self.SHAPES}, got '{shape}'.")

        dp = {
            "fill": "lightblue",  # Gviz AnnotationTrack default
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=dp,
        )

        self._feature_colors = dict(_DEFAULT_FEATURE_COLORS)
        if feature_colors:
            self._feature_colors.update(feature_colors)

        self.shape = shape
        self.show_label = show_label
        self.label_pos = label_pos
        self.arrow_direction = arrow_direction
        self.group_annotation = group_annotation

    def _get_feature_color(self, feature: Optional[str], idx: int, color_cycle) -> str:
        """Get color for a feature, using feature_colors dict or cycling."""
        if feature is not None and feature in self._feature_colors:
            return self._feature_colors[feature]
        if feature is not None and pd.notna(feature):
            # Assign color from cycle for unknown features
            return next(color_cycle)
        return self.get_param("fill", "lightblue")

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

        # Color cycle for unknown features (Gviz-inspired palette)
        default_colors = ["lightblue", "#0080FF", "#DC0000", "#F8766D",
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
                    color = self.get_param("fill", "lightblue")
            else:
                color = self.get_param("fill", "lightblue")

            edge_color = self.get_param("col", "#333333")  # Bug 4 fix: use 'col' not 'col_border'
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
            elif self.shape == "fixedArrow":
                self._draw_fixed_arrow(ax, x_start, y_center, width, feat_height,
                                       color, edge_color, lwd, alpha, row.get("strand", "*"))
            elif self.shape == "smallArrow":
                self._draw_small_arrow(ax, x_start, y_center, width, feat_height,
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

        # Draw group annotations with connecting lines (Task 13)
        if self.group_annotation and "group" in sub.columns:
            self._draw_group_annotations(ax, sub, region, stacks, stack_info)

    def _draw_group_annotations(self, ax, data, region, stacks, stack_info):
        """Draw group labels and connecting lines for grouped features."""
        if "group" not in data.columns:
            return

        y_positions = stack_info["y_positions"]
        total_rows = stack_info["total_rows"]
        line_color = self.get_param("col_line", "#888888")
        fontsize = self.get_param("fontsize", 7)

        # Determine label column
        label_col_map = {"id": "name", "group": "group", "feature": "feature"}
        label_col = label_col_map.get(self.group_annotation, "group")
        if label_col not in data.columns:
            label_col = "group"
        if label_col not in data.columns:
            return

        groups = data.groupby("group")
        for grp_name, grp_data in groups:
            if pd.isna(grp_name) or grp_name == ".":
                continue

            # Get group extent
            grp_start = grp_data["start"].min()
            grp_end = grp_data["end"].max()

            # Draw connecting line between features in the group
            if len(grp_data) > 1:
                y_line = 0.02  # Near bottom
                ax.plot([max(grp_start, region.start), min(grp_end, region.end)],
                        [y_line, y_line],
                        color=line_color, linewidth=0.5, alpha=0.5, zorder=1)

            # Draw group label
            actual_col = label_col if label_col in grp_data.columns else "group"
            label = str(grp_data[actual_col].iloc[0]) if actual_col in grp_data.columns else str(grp_name)
            x_label = (max(grp_start, region.start) + min(grp_end, region.end)) / 2
            ax.text(x_label, -0.03, label,
                    ha="center", va="top", fontsize=fontsize,
                    color="#555555", fontstyle="italic", clip_on=True, zorder=5)

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

    def _draw_fixed_arrow(self, ax, x, y, w, h, color, edge_color, lwd, alpha, strand="*"):
        """Draw a fixed-width arrow feature (arrowHeadWidth param controls head size)."""
        head_width = self.get_param("arrow_head_width", None)
        if head_width is None:
            # Default: 20% of feature width or h*2, whichever is smaller
            head_width = min(w * 0.2, h * 2)
        else:
            # Convert pixel-like param to data coords (approximate)
            span = ax.get_xlim()
            pixel_scale = (span[1] - span[0]) / 800  # approximate pixels to data
            head_width = head_width * pixel_scale

        if strand == "-":
            verts = [
                (x + head_width, y - h / 2),
                (x + w, y - h / 2),
                (x + w, y + h / 2),
                (x + head_width, y + h / 2),
                (x, y),
                (x + head_width, y - h / 2),
            ]
        else:
            verts = [
                (x, y - h / 2),
                (x + w - head_width, y - h / 2),
                (x + w, y),
                (x + w - head_width, y + h / 2),
                (x, y + h / 2),
                (x, y - h / 2),
            ]

        polygon = mpatches.Polygon(
            verts, closed=True,
            facecolor=color, edgecolor=edge_color,
            linewidth=lwd, alpha=alpha, zorder=3,
        )
        ax.add_patch(polygon)

    def _draw_small_arrow(self, ax, x, y, w, h, color, edge_color, lwd, alpha, strand="*"):
        """Draw a small arrow (50% of regular arrow head)."""
        arrow_width = min(w * 0.08, h * 1.0)  # Half the regular arrow

        if strand == "-":
            verts = [
                (x + arrow_width, y - h / 2),
                (x + w, y - h / 2),
                (x + w, y + h / 2),
                (x + arrow_width, y + h / 2),
                (x, y),
                (x + arrow_width, y - h / 2),
            ]
        else:
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
