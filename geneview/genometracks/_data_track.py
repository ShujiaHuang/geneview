"""
DataTrack - A track for visualizing numeric genomic data.

Supports multiple plot types: line, histogram, heatmap, polygon (area),
gradient, points, and mountain. Multi-sample data with grouping is supported.

Ported from Gviz's DataTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PolyCollection

from ._base import NumericTrack, GenomicInterval


# Supported plot types
_PLOT_TYPES = (
    "line", "histogram", "heatmap", "polygon", "gradient",
    "points", "mountain", "boxplot",
    "b", "s", "S",  # Gviz-compatible: combined, stair-post, stair-pre
)

# Aggregation functions
_AGG_FUNCS = {
    "mean": np.nanmean,
    "median": np.nanmedian,
    "sum": np.nansum,
    "min": np.nanmin,
    "max": np.nanmax,
}

# Default color palettes
_HEATMAP_CMAP = "RdYlBu_r"
_GRADIENT_COLORS = ("white", "#3C5488")


class DataTrack(NumericTrack):
    """A track for visualizing numeric genomic data along coordinates.

    Supports multiple visualization modes including line plots, histograms,
    heatmaps, area plots (polygon), gradients, and scatter points.

    Parameters
    ----------
    data : pd.DataFrame or str
        DataFrame with 'chrom', 'start', 'end' plus numeric value column(s).
        If a string, interpreted as a file path (BigWig/bedGraph auto-detected).
    type : str
        Plot type: 'line', 'histogram', 'heatmap', 'polygon', 'gradient',
        'points', 'mountain', or 'boxplot'. Default is 'line'.
    value_columns : list of str, optional
        Column names to use as values. Auto-detected if None.
    groups : list or None, optional
        Group assignments for multi-sample data.
    ylim : tuple, optional
        Y-axis limits. Auto-computed if None.
    baseline : float, optional
        Baseline value for polygon/mountain plots. Default is 0.
    col : str or list, optional
        Color(s) for the data. Default is '#3C5488'.
    fill : str, optional
        Fill color for histograms/polygons. Default is '#5B8DB8'.
    gradient : tuple, optional
        Two colors for gradient/heatmap. Default is ('white', '#3C5488').
    ncolor : int, optional
        Number of gradient levels. Default is 100.
    show_sample_names : bool, optional
        For heatmap type, show sample names on the y-axis. Default is True.
    transformation : callable, optional
        Function applied to value arrays before plotting (e.g. np.log2).
    window : int or str, optional
        Binning/windowing control. Positive int: bin into N windows.
        Negative int: running window of that size. "auto": auto-detect.
    window_size : int, optional
        Fixed running window size in bins.
    aggregation : str
        Aggregation method for windowed/grouped data: "mean", "median",
        "sum", "min", "max". Default is "mean".
    legend : bool
        If True and groups are set, show a legend. Default is False.
    name : str
        Track name. Default is "Data".
    height : float
        Relative track height. Default is 1.5.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from geneview.genometracks import DataTrack, plot_tracks, GenomicInterval
    >>> np.random.seed(42)
    >>> data = pd.DataFrame({
    ...     "chrom": ["chr7"] * 50,
    ...     "start": np.arange(2000000, 2050000, 1000),
    ...     "end":   np.arange(2001000, 2051000, 1000),
    ...     "value": np.random.randn(50).cumsum(),
    ... })
    >>> track = DataTrack(data, type="histogram")
    >>> plot_tracks([track], region=GenomicInterval("chr7", 2000000, 2050000))
    """

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        type: str = "line",
        value_columns: Optional[List[str]] = None,
        groups: Optional[List] = None,
        ylim: Optional[Tuple[float, float]] = None,
        baseline: float = 0.0,
        col: Union[str, List[str]] = "#3C5488",
        fill: str = "#5B8DB8",
        gradient: Tuple[str, str] = ("white", "#3C5488"),
        ncolor: int = 100,
        show_sample_names: bool = True,
        transformation: Optional[Any] = None,
        window: Optional[Union[int, str]] = None,
        window_size: Optional[int] = None,
        aggregation: str = "mean",
        legend: bool = False,
        name: str = "Data",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if type not in _PLOT_TYPES:
            raise ValueError(
                f"Plot type must be one of {_PLOT_TYPES}, got '{type}'."
            )

        dp = {
            "col": col,
            "fill": fill,
            "fontsize": 8,
            "fontcolor": "#555555",
            "col_axis": "#666666",
            "lwd": 1.5,
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, value_columns=value_columns, name=name,
            height=height, display_params=dp,
        )

        self.plot_type = type
        self.groups = groups
        self._ylim = ylim
        self.baseline = baseline
        self.gradient_colors = gradient
        self.ncolor = ncolor
        self.show_sample_names = show_sample_names
        self.transformation = transformation
        self.window = window
        self.window_size = window_size
        self.aggregation = aggregation
        self.legend = legend

    def get_ylim(self) -> Tuple[float, float]:
        """Get y-axis limits, using user-provided or auto-computed."""
        if self._ylim is not None:
            return self._ylim
        return super().get_ylim()

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the data track.

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
            ax.set_ylim(*self.get_ylim())
            self._format_y_axis(ax)
            return

        ax.set_xlim(region.start, region.end)

        # Dispatch to the appropriate plot method
        draw_methods = {
            "line": self._draw_line,
            "histogram": self._draw_histogram,
            "heatmap": self._draw_heatmap,
            "polygon": self._draw_polygon,
            "gradient": self._draw_gradient,
            "points": self._draw_points,
            "mountain": self._draw_mountain,
            "boxplot": self._draw_boxplot,
            "b": self._draw_combined,
            "s": self._draw_stairs_post,
            "S": self._draw_stairs_pre,
        }

        # Apply windowing if configured
        plot_data, windowed_values = self._apply_window(sub)

        draw_method = draw_methods.get(self.plot_type, self._draw_line)
        if self.window is not None and self.plot_type not in ("heatmap", "boxplot"):
            # Use windowed data for applicable types
            draw_method(ax, plot_data, region, _precomputed=windowed_values)
        else:
            draw_method(ax, sub, region)

        # Draw grid if requested
        if self.get_param("grid", False):
            ax.grid(axis="y", color=self.get_param("col_grid", "#DDDDDD"),
                    linewidth=0.3, alpha=0.7, zorder=1)

        self._format_y_axis(ax)

        # Draw legend if requested
        if self.legend and self.groups is not None:
            ax.legend(loc="upper right", fontsize=6, framealpha=0.8)

    def _format_y_axis(self, ax):
        """Format the y-axis with appropriate labels."""
        ylim = self.get_ylim()
        ax.set_ylim(ylim)

        # Show y-axis with a few tick marks
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color(self.get_param("col_axis", "#666666"))
        ax.spines["left"].set_linewidth(0.5)
        ax.spines["bottom"].set_visible(False)

        ax.yaxis.set_major_locator(plt.MaxNLocator(3))
        ax.tick_params(axis="y", labelsize=7, colors=self.get_param("col_axis", "#666666"))
        ax.set_xticklabels([])

    def _get_value_arrays(self, data: pd.DataFrame) -> List[Tuple[str, np.ndarray]]:
        """Get value columns and their arrays, applying transformation if set."""
        result = []
        for col in self.value_columns:
            if col in data.columns:
                vals = data[col].values.astype(float)
                if self.transformation is not None:
                    vals = self.transformation(vals)
                result.append((col, vals))
        return result

    def _apply_window(self, data: pd.DataFrame) -> Tuple[pd.DataFrame, List[Tuple[str, np.ndarray]]]:
        """Apply windowing/smoothing to data if configured.

        Returns modified data and value arrays.
        """
        if self.window is None:
            return data, self._get_value_arrays(data)

        value_arrays = self._get_value_arrays(data)
        if not value_arrays:
            return data, value_arrays

        agg_func = _AGG_FUNCS.get(self.aggregation, np.nanmean)
        n = len(data)

        if self.window == "auto":
            win = max(1, n // 50)
        elif isinstance(self.window, str):
            win = self.window_size or max(1, n // 50)
        elif isinstance(self.window, int):
            if self.window > 0:
                # Bin into N equal windows
                win = max(1, n // self.window)
            else:
                win = abs(self.window)
        else:
            return data, value_arrays

        if win <= 1:
            return data, value_arrays

        # Build windowed data
        new_starts = []
        new_ends = []
        new_vals_dict = {name: [] for name, _ in value_arrays}
        starts_arr = data["start"].values
        ends_arr = data["end"].values

        for i in range(0, n, win):
            chunk_end = min(i + win, n)
            new_starts.append(int(starts_arr[i]))
            new_ends.append(int(ends_arr[chunk_end - 1]))
            for name, vals in value_arrays:
                chunk = vals[i:chunk_end]
                new_vals_dict[name].append(agg_func(chunk))

        new_data = pd.DataFrame({
            "chrom": [data["chrom"].iloc[0]] * len(new_starts),
            "start": new_starts,
            "end": new_ends,
        })
        new_value_arrays = [(name, np.array(new_vals_dict[name])) for name, _ in value_arrays]
        return new_data, new_value_arrays

    def _get_midpoints(self, data: pd.DataFrame) -> np.ndarray:
        """Get midpoint positions for data bins."""
        return (data["start"].values + data["end"].values) / 2.0

    def _draw_line(self, ax, data, region, _precomputed=None):
        """Draw data as line plot."""
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)

        lwd = self.get_param("lwd", 1.5)
        alpha = self.get_param("alpha", 1.0)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            ax.plot(midpoints, values, color=color, linewidth=lwd,
                    alpha=alpha, zorder=3, label=col_name)

    def _draw_histogram(self, ax, data, region, _precomputed=None):
        """Draw data as histogram (vertical bars)."""
        fill = self.get_param("fill", "#5B8DB8")
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        alpha = self.get_param("alpha", 0.8)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)

        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            # Draw vertical bars centered on each bin
            for j, (s, e, v) in enumerate(zip(data["start"], data["end"], values)):
                if np.isnan(v):
                    continue
                width = e - s
                # Fix Bug 2: draw bars from baseline correctly
                bar_bottom = min(self.baseline, v)
                bar_height = abs(v - self.baseline)
                bar_alpha = alpha if v >= self.baseline else alpha * 0.7
                ax.bar(s + width / 2, bar_height,
                       width=width * 0.9, bottom=bar_bottom,
                       color=color, edgecolor="none", alpha=bar_alpha, zorder=3)

            # Draw baseline
            ax.axhline(y=self.baseline, color="#888888", linewidth=0.5,
                       linestyle="--", zorder=2)

    def _draw_polygon(self, ax, data, region, _precomputed=None):
        """Draw data as filled polygon (area plot)."""
        fill = self.get_param("fill", "#5B8DB8")
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        alpha = self.get_param("alpha", 0.5)
        lwd = self.get_param("lwd", 1.0)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            # Create polygon: baseline -> data -> baseline
            x = np.concatenate([[midpoints[0]], midpoints, [midpoints[-1]]])
            y = np.concatenate([[self.baseline], values, [self.baseline]])

            ax.fill(x, y, color=color, alpha=alpha, zorder=3)
            ax.plot(midpoints, values, color=color, linewidth=lwd,
                    alpha=min(1.0, alpha + 0.3), zorder=4)

    def _draw_heatmap(self, ax, data, region):
        """Draw data as a heatmap."""
        cmap_name = self.get_param("heatmap_cmap", _HEATMAP_CMAP)
        cmap = plt.get_cmap(cmap_name)

        value_arrays = self._get_value_arrays(data)
        if not value_arrays:
            return

        # Build 2D array (samples x positions)
        matrix = np.array([vals for _, vals in value_arrays])

        # Normalize to 0-1 for colormap
        vmin = np.nanmin(matrix)
        vmax = np.nanmax(matrix)
        if vmin == vmax:
            vmax = vmin + 1

        n_samples = matrix.shape[0]
        n_positions = matrix.shape[1]

        # Draw as image
        starts = data["start"].values
        ends = data["end"].values
        extent = [starts[0], ends[-1], 0, n_samples]

        ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                  extent=extent, origin="lower", interpolation="nearest", zorder=3)

        # Show sample names on y-axis if requested
        if self.show_sample_names and n_samples > 1:
            sample_names = [name for name, _ in value_arrays]
            ax.set_yticks(np.arange(n_samples) + 0.5)
            ax.set_yticklabels(sample_names, fontsize=6)
        else:
            ax.set_yticks([])

    def _draw_gradient(self, ax, data, region):
        """Draw collapsed average as a color gradient."""
        colors = self.gradient_colors
        cmap = LinearSegmentedColormap.from_list("gradient", colors, N=self.ncolor)

        value_arrays = self._get_value_arrays(data)
        if not value_arrays:
            return

        # Average across samples
        matrix = np.array([vals for _, vals in value_arrays])
        avg = np.nanmean(matrix, axis=0)

        vmin = np.nanmin(avg)
        vmax = np.nanmax(avg)
        if vmin == vmax:
            vmax = vmin + 1

        # Draw as a single-row image
        starts = data["start"].values
        ends = data["end"].values
        extent = [starts[0], ends[-1], 0, 1]

        ax.imshow(avg.reshape(1, -1), aspect="auto", cmap=cmap,
                  vmin=vmin, vmax=vmax, extent=extent,
                  origin="lower", interpolation="nearest", zorder=3)
        ax.set_yticks([])

    def _draw_points(self, ax, data, region, _precomputed=None):
        """Draw data as scatter points."""
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        alpha = self.get_param("alpha", 0.8)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            mask = np.isfinite(values)
            ax.scatter(midpoints[mask], values[mask], c=color,
                       s=10, alpha=alpha, zorder=3, label=col_name)

    def _draw_mountain(self, ax, data, region, _precomputed=None):
        """Draw smoothed mountain plot (area above baseline)."""
        fill = self.get_param("fill", "#5B8DB8")
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        alpha = self.get_param("alpha", 0.5)
        lwd = self.get_param("lwd", 1.0)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            # Simple smoothing via running average
            kernel_size = max(3, len(values) // 20)
            if kernel_size % 2 == 0:
                kernel_size += 1
            kernel = np.ones(kernel_size) / kernel_size
            smoothed = np.convolve(values, kernel, mode="same")

            # Area above baseline
            baseline_vals = np.maximum(smoothed, self.baseline)
            x = np.concatenate([[midpoints[0]], midpoints, [midpoints[-1]]])
            y = np.concatenate([[self.baseline], baseline_vals, [self.baseline]])

            ax.fill(x, y, color=color, alpha=alpha, zorder=3)
            ax.plot(midpoints, baseline_vals, color=color, linewidth=lwd,
                    alpha=min(1.0, alpha + 0.3), zorder=4)

            ax.axhline(y=self.baseline, color="#888888", linewidth=0.5,
                       linestyle="--", zorder=2)

    def _draw_boxplot(self, ax, data, region):
        """Draw data as boxplots at each genomic position."""
        col = self.get_param("col", "#3C5488")
        alpha = self.get_param("alpha", 0.8)
        midpoints = self._get_midpoints(data)

        value_arrays = self._get_value_arrays(data)
        if not value_arrays:
            return

        # Build matrix (samples x positions) and draw boxplot per position
        matrix = np.array([vals for _, vals in value_arrays])

        positions_to_plot = []
        data_to_plot = []
        for j in range(matrix.shape[1]):
            col_data = matrix[:, j]
            col_data = col_data[np.isfinite(col_data)]
            if len(col_data) >= 1:
                positions_to_plot.append(midpoints[j])
                data_to_plot.append(col_data)

        if data_to_plot:
            bp = ax.boxplot(data_to_plot, positions=positions_to_plot,
                            widths=(midpoints[1] - midpoints[0]) * 0.6 if len(midpoints) > 1 else 100,
                            patch_artist=True, showfliers=False, zorder=3)
            for patch in bp["boxes"]:
                patch.set_facecolor(col if isinstance(col, str) else col[0])
                patch.set_alpha(alpha)

    def _draw_combined(self, ax, data, region, _precomputed=None):
        """Draw data as combined line + points (Gviz type 'b')."""
        self._draw_points(ax, data, region, _precomputed=_precomputed)
        self._draw_line(ax, data, region, _precomputed=_precomputed)

    def _draw_stairs_post(self, ax, data, region, _precomputed=None):
        """Draw stair-step plot, horizontal first then vertical (Gviz type 's')."""
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        lwd = self.get_param("lwd", 1.5)
        alpha = self.get_param("alpha", 1.0)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            ax.step(midpoints, values, where="post", color=color,
                    linewidth=lwd, alpha=alpha, zorder=3, label=col_name)

    def _draw_stairs_pre(self, ax, data, region, _precomputed=None):
        """Draw stair-step plot, vertical first then horizontal (Gviz type 'S')."""
        col = self.get_param("col", "#3C5488")
        if isinstance(col, list):
            colors = col
        else:
            colors = [col] * len(self.value_columns)
        lwd = self.get_param("lwd", 1.5)
        alpha = self.get_param("alpha", 1.0)
        midpoints = self._get_midpoints(data)

        value_arrays = _precomputed if _precomputed else self._get_value_arrays(data)
        for i, (col_name, values) in enumerate(value_arrays):
            color = colors[i % len(colors)]
            ax.step(midpoints, values, where="pre", color=color,
                    linewidth=lwd, alpha=alpha, zorder=3, label=col_name)
