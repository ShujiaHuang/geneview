"""
GenomeAxisTrack - A genomic coordinate axis track.

Displays a customizable genomic coordinate ruler with tick marks and labels.
Automatically formats labels (bp/kb/Mb/Gb) based on the zoom level.
Optionally shows a scale bar instead of a full axis.

Ported from Gviz's GenomeAxisTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from ._base import Track, GenomicInterval, _format_genomic_position


# "Nice" tick intervals for genomic coordinates (in base pairs)
_NICE_INTERVALS = [
    1, 2, 5, 10, 20, 50, 100, 200, 500,
    1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000,
    1_000_000, 2_000_000, 5_000_000, 10_000_000, 20_000_000, 50_000_000,
    100_000_000, 200_000_000, 500_000_000, 1_000_000_000,
]


def _choose_tick_interval(span: int, target_ticks: int = 8) -> int:
    """Choose a human-readable tick interval for the given span."""
    ideal = span / target_ticks
    for interval in _NICE_INTERVALS:
        if interval >= ideal:
            return interval
    return _NICE_INTERVALS[-1]


class GenomeAxisTrack(Track):
    """A track displaying a genomic coordinate axis.

    Creates a ruler-like axis showing genomic positions with automatically
    formatted labels (bp/kb/Mb/Gb). Can optionally show a simple scale
    bar instead of a full axis.

    Parameters
    ----------
    ranges : pd.DataFrame, optional
        DataFrame with 'chrom', 'start', 'end' columns to highlight
        specific regions on the axis. Default is None (no highlights).
    name : str
        Track name. Default is "Axis".
    height : float
        Relative track height. Default is 0.3.
    scale : float or None, optional
        If not None, show a scale bar instead of a full axis.
        If 0 < scale <= 1, interpreted as a fraction of the plotted region.
        If scale > 1, interpreted as an absolute length in base pairs.
    label_pos : str
        Position of tick labels: 'above', 'below', or 'alternating'.
        Default is 'alternating'.
    little_ticks : bool
        If True, add smaller sub-ticks between major ticks. Default is False.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> from geneview.genometracks import GenomeAxisTrack, plot_tracks, GenomicInterval
    >>> ax_track = GenomeAxisTrack()
    >>> region = GenomicInterval("chr7", 2000000, 2200000)
    >>> plot_tracks([ax_track], region=region)
    """

    def __init__(
        self,
        ranges: Optional[pd.DataFrame] = None,
        name: str = "Axis",
        height: float = 0.3,
        scale: Optional[float] = None,
        label_pos: str = "alternating",
        little_ticks: bool = False,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        dp = {
            "background_title": "transparent",
            "col_border_title": "transparent",
            "show_title": False,
            "col": "#555555",
            "col_range": "#8B8378",
            "fill_range": "#CDC8B1",
            "lwd": 2,
            "fontsize": 9,
            "fontcolor": "#555555",
        }
        if display_params:
            dp.update(display_params)

        super().__init__(name=name, height=height, display_params=dp)

        self._ranges = ranges
        self.scale = scale
        self.label_pos = label_pos
        self.little_ticks = little_ticks

    def get_region(self) -> Optional[GenomicInterval]:
        """GenomeAxisTrack doesn't define its own region."""
        return None

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the genome axis on the given axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        if self.scale is not None:
            self._draw_scale_bar(ax, region)
        else:
            self._draw_full_axis(ax, region)

    def _draw_full_axis(self, ax, region: GenomicInterval) -> None:
        """Draw a full genomic coordinate axis with ticks and labels."""
        span = region.end - region.start
        color = self.get_param("col", "#555555")
        lwd = self.get_param("lwd", 2)
        fontsize = self.get_param("fontsize", 9)
        fontcolor = self.get_param("fontcolor", "#555555")

        # Draw the main axis line (in data coordinates; ylim is set to [0,1] below)
        y_center = 0.5
        ax.plot([region.start, region.end], [y_center, y_center],
                color=color, linewidth=lwd, zorder=1)

        # Draw highlighted ranges if provided
        if self._ranges is not None and len(self._ranges) > 0:
            range_data = self._ranges[
                (self._ranges["chrom"] == region.chrom) &
                (self._ranges["start"] < region.end) &
                (self._ranges["end"] > region.start)
            ]
            rfill = self.get_param("fill_range", "#CDC8B1")
            rcolor = self.get_param("col_range", "#8B8378")
            for _, row in range_data.iterrows():
                s = max(row["start"], region.start)
                e = min(row["end"], region.end)
                ax.axvspan(s, e, alpha=0.4, color=rfill, zorder=1)

        # Calculate tick positions
        tick_interval = _choose_tick_interval(span)
        first_tick = np.ceil(region.start / tick_interval) * tick_interval
        ticks = np.arange(first_tick, region.end, tick_interval)

        # Draw ticks and labels
        tick_height = 0.15
        label_positions = []

        for i, tick_pos in enumerate(ticks):
            # Major tick
            ax.plot([tick_pos, tick_pos],
                    [y_center - tick_height, y_center + tick_height],
                    color=color, linewidth=lwd * 0.7, zorder=2)

            # Determine label position based on label_pos setting
            if self.label_pos == "alternating":
                above = (i % 2 == 0)
            elif self.label_pos == "above":
                above = True
            else:  # below
                above = False

            label = _format_genomic_position(tick_pos, span)
            y_label = y_center + tick_height + 0.05 if above else y_center - tick_height - 0.05
            va = "bottom" if above else "top"

            ax.text(tick_pos, y_label, label,
                    ha="center", va=va, fontsize=fontsize,
                    color=fontcolor, zorder=3)
            label_positions.append((above, y_label))

        # Draw minor ticks if requested
        if self.little_ticks and tick_interval > 1:
            minor_interval = tick_interval / 5
            first_minor = np.ceil(region.start / minor_interval) * minor_interval
            minor_ticks = np.arange(first_minor, region.end, minor_interval)
            # Remove major tick positions
            major_set = set(ticks)
            minor_ticks = [t for t in minor_ticks if t not in major_set]

            minor_height = tick_height * 0.5
            for tick_pos in minor_ticks:
                ax.plot([tick_pos, tick_pos],
                        [y_center - minor_height, y_center + minor_height],
                        color=color, linewidth=lwd * 0.4, alpha=0.6, zorder=2)

        # Configure axes
        ax.set_xlim(region.start, region.end)
        ax.set_ylim(0, 1)
        ax.axis("off")

    def _draw_scale_bar(self, ax, region: GenomicInterval) -> None:
        """Draw a simple scale bar instead of a full axis."""
        span = region.end - region.start
        color = self.get_param("col", "#555555")
        lwd = self.get_param("lwd", 2)
        fontsize = self.get_param("fontsize", 9)
        fontcolor = self.get_param("fontcolor", "#555555")

        # Calculate scale bar length
        scale_len = self.scale
        if 0 < scale_len <= 1:
            scale_len = span * scale_len
            # Round to a nice number
            exp = np.floor(np.log10(scale_len))
            scale_len = round(scale_len, -int(exp))
            if scale_len == 0:
                scale_len = span * self.scale

        if scale_len > span:
            scale_len = span * 0.05

        # Position scale bar at left side with some margin
        x_offset = region.start + span * 0.03
        y_center = 0.5
        tick_height = 0.2

        # Draw scale bar line
        ax.plot([x_offset, x_offset + scale_len],
                [y_center, y_center],
                color=color, linewidth=lwd, zorder=2)

        # Draw end caps
        for x in [x_offset, x_offset + scale_len]:
            ax.plot([x, x],
                    [y_center - tick_height, y_center + tick_height],
                    color=color, linewidth=lwd, zorder=2)

        # Draw label
        label = _format_genomic_position(scale_len, scale_len)
        label_x = x_offset + scale_len / 2

        if self.label_pos == "above":
            y_label = y_center + tick_height + 0.05
            va = "bottom"
        elif self.label_pos == "below":
            y_label = y_center - tick_height - 0.05
            va = "top"
        else:  # beside
            label_x = x_offset + scale_len + span * 0.01
            y_label = y_center
            va = "center"

        ax.text(label_x, y_label, label,
                ha="center", va=va, fontsize=fontsize,
                color=fontcolor, fontweight="bold", zorder=3)

        # Configure axes
        ax.set_xlim(region.start, region.end)
        ax.set_ylim(0, 1)
        ax.axis("off")
