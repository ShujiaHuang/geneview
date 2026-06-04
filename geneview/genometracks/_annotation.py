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
    >>> _ = plot_tracks([track], region=GenomicInterval("chr7", 1900000, 2200000))
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
        just_group: str = "left",
        show_overplotting: bool = False,
        merge_groups: bool = False,
        col_line: Optional[str] = None,
        name: str = "Annotation",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        if shape not in self.SHAPES:
            raise ValueError(f"Shape must be one of {self.SHAPES}, got '{shape}'.")

        dp = {
            "fill": "lightblue",  # Gviz AnnotationTrack default
            "col": "transparent",  # Gviz: AnnotationTrack col = "transparent"
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=dp, **kwargs,
        )

        self._feature_colors = dict(_DEFAULT_FEATURE_COLORS)
        if feature_colors:
            self._feature_colors.update(feature_colors)

        self.shape = shape
        self.show_label = show_label
        self.label_pos = label_pos
        self.arrow_direction = arrow_direction
        self.group_annotation = group_annotation
        self.just_group = just_group
        self.show_overplotting = show_overplotting
        self.merge_groups = merge_groups
        if col_line is not None:
            self.set_param("col_line", col_line)

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
        stack_info = get_stack_heights(
            stacks, mode=self.stacking, stack_height_frac=self.stack_height,
            reverse_stacking=self.reverse_stacking,
        )

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

            edge_color = self.get_param("col", "transparent")  # Gviz: col = "transparent"
            # Matplotlib doesn't understand "transparent", convert to "none"
            if edge_color == "transparent":
                edge_color = "none"
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

    @classmethod
    def from_bam(
        cls,
        filepath: str,
        region: Optional[GenomicInterval] = None,
        name: str = "BAM Annotation",
        **kwargs,
    ) -> "AnnotationTrack":
        """Create an AnnotationTrack from BAM file read coordinates.

        Requires the optional ``pysam`` package.

        Parameters
        ----------
        filepath : str
            Path to a BAM or CRAM file.
        region : GenomicInterval, optional
            Genomic region to extract reads from.
        name : str
            Track name.
        **kwargs
            Additional keyword arguments passed to ``AnnotationTrack``.

        Returns
        -------
        AnnotationTrack
        """
        try:
            import pysam
        except ImportError:
            raise ImportError(
                "The 'pysam' package is required for from_bam(). "
                "Install it with: pip install pysam"
            )

        aln = pysam.AlignmentFile(filepath, "rb")
        try:
            if region is not None:
                reads = aln.fetch(region.chrom, region.start, region.end)
            else:
                reads = aln.fetch()

            rows = []
            for read in reads:
                if read.is_unmapped:
                    continue
                rows.append({
                    "chrom": read.reference_name,
                    "start": read.reference_start,
                    "end": read.reference_end,
                    "strand": "-" if read.is_reverse else "+",
                    "name": read.query_name,
                })
        finally:
            aln.close()

        df = pd.DataFrame(rows)
        if len(df) == 0:
            df = pd.DataFrame(columns=["chrom", "start", "end", "strand", "name"])

        return cls(data=df, name=name, **kwargs)


class DetailsAnnotationTrack(AnnotationTrack):
    """AnnotationTrack with detail panels drawn below each feature.

    For each selected feature, a custom function is called to draw a
    detail panel in the lower portion of the track, connected to the
    feature by a thin line.

    Parameters
    ----------
    data : pd.DataFrame or str
        Same as AnnotationTrack.
    fun : callable, optional
        Function ``fun(ax, identifier, region, **kwargs)`` called for
        each detail panel.  If None, a default density-like plot is drawn.
    select_fun : callable, optional
        Function ``select_fun(row) -> bool`` to decide which features
        get detail panels.  Default: all features.
    details_size : float
        Fraction of track height reserved for detail panels.  Default 0.4.
    details_connector_col : str
        Color for connector lines.  Default ``"gray"``.
    details_connector_lty : str
        Line style for connectors.  Default ``"--"``.
    details_connector_lwd : float
        Line width for connectors.  Default 0.5.
    details_border_col : str
        Border color for detail panels.  Default ``"gray"``.
    details_border_fill : str
        Fill color for detail panel background.  Default ``"white"``.
    details_min_width : int
        Minimum feature width (bp) to qualify for a detail panel.
        Default 0.
    details_ratio : float
        Width ratio of detail panel relative to feature width.
        Default 1.0.
    group_details : bool
        If True, draw one detail panel per group instead of per feature.
        Default False.
    details_fun_args : dict, optional
        Extra keyword arguments passed to *fun*.
    **kwargs
        Additional keyword arguments passed to AnnotationTrack.
    """

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        fun: Optional[Callable] = None,
        select_fun: Optional[Callable] = None,
        details_size: float = 0.4,
        details_connector_col: str = "gray",
        details_connector_lty: str = "--",
        details_connector_lwd: float = 0.5,
        details_border_col: str = "gray",
        details_border_fill: str = "white",
        details_min_width: int = 0,
        details_ratio: float = 1.0,
        group_details: bool = False,
        details_fun_args: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        super().__init__(data=data, **kwargs)
        self.fun = fun
        self.select_fun = select_fun
        self.details_size = details_size
        self.details_connector_col = details_connector_col
        self.details_connector_lty = details_connector_lty
        self.details_connector_lwd = details_connector_lwd
        self.details_border_col = details_border_col
        self.details_border_fill = details_border_fill
        self.details_min_width = details_min_width
        self.details_ratio = details_ratio
        self.group_details = group_details
        self.details_fun_args = details_fun_args or {}

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw annotation features and detail panels.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        # Draw the main annotation track in the upper portion
        sub = self.subset(region)
        if sub is None or len(sub) == 0:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        # Adjust y-axis to make room for details
        ax.set_xlim(region.start, region.end)
        ax.set_ylim(-0.05, 1.05)

        # Delegate to parent draw for features
        super().draw(ax, region)

        # Collect features that qualify for detail panels
        detail_features = []
        if "strand" not in sub.columns:
            sub["strand"] = "*"

        for _, row in sub.iterrows():
            feat_width = row["end"] - row["start"]
            if feat_width < self.details_min_width:
                continue
            if self.select_fun is not None and not self.select_fun(row):
                continue
            detail_features.append(row)

        if not detail_features:
            return

        # Draw detail panels
        n_details = len(detail_features)
        if n_details == 0:
            return

        detail_y_bottom = -0.3
        detail_y_top = 0.0
        detail_height = (detail_y_top - detail_y_bottom) * self.details_size

        span = region.end - region.start
        detail_width_bp = span / max(n_details, 1) * self.details_ratio

        for i, row in enumerate(detail_features):
            feat_mid = (row["start"] + row["end"]) / 2
            detail_start = feat_mid - detail_width_bp / 2
            detail_end = feat_mid + detail_width_bp / 2

            # Draw connector line
            ax.plot(
                [feat_mid, feat_mid],
                [0.0, detail_y_top],
                color=self.details_connector_col,
                linestyle=self.details_connector_lty,
                linewidth=self.details_connector_lwd,
                alpha=0.5, zorder=1,
            )

            # Draw detail panel background
            rect = mpatches.Rectangle(
                (detail_start, detail_y_bottom),
                detail_width_bp,
                detail_height + abs(detail_y_bottom),
                facecolor=self.details_border_fill,
                edgecolor=self.details_border_col,
                linewidth=0.5, alpha=0.9, zorder=2,
            )
            ax.add_patch(rect)

            # Call custom function or draw default
            identifier = row.get("name", row.get("id", f"feat_{i}"))
            if self.fun is not None:
                try:
                    self.fun(ax, identifier, region, **self.details_fun_args)
                except Exception:
                    pass
            else:
                self._draw_default_detail(ax, detail_start, detail_end,
                                          detail_y_bottom, detail_height)

    def _draw_default_detail(self, ax, x_start, x_end, y_bottom, height):
        """Draw a default detail panel (simple density-like bar)."""
        mid_x = (x_start + x_end) / 2
        w = (x_end - x_start) * 0.6
        h = abs(height) * 0.6
        y_center = y_bottom + abs(height) * 0.5

        rect = mpatches.FancyBboxPatch(
            (mid_x - w / 2, y_center - h / 2), w, h,
            boxstyle="round,pad=0",
            facecolor=self.get_param("fill", "lightblue"),
            edgecolor=self.details_border_col,
            linewidth=0.3, alpha=0.7, zorder=3,
        )
        ax.add_patch(rect)


# Register class-specific display parameter overrides
from ._base import _CLASS_DISPLAY_PARAM_OVERRIDES
_CLASS_DISPLAY_PARAM_OVERRIDES["AnnotationTrack"] = {
    "fill": "lightblue",
    "col": "transparent",
}
