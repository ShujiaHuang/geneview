"""
HighlightTrack - A container track for adding highlighted regions across tracks.

Allows overlaying colored highlight regions on multiple tracks simultaneously.
The highlight regions are defined by genomic intervals and drawn as
semi-transparent rectangles spanning the full height of each contained track.

Ported from Gviz's HighlightTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ._base import Track, RangeTrack, GenomicInterval


class HighlightTrack(RangeTrack):
    """A track that highlights genomic regions across multiple tracks.

    HighlightTrack wraps other track objects and adds colored rectangular
    overlays at specified genomic positions. The highlight is drawn across
    all contained tracks during rendering.

    This track is used in conjunction with ``plot_tracks()``, which handles
    expanding the contained tracks and applying highlights.

    Parameters
    ----------
    track_list : list of Track
        List of Track objects to wrap with highlighting.
    regions : pd.DataFrame
        DataFrame with 'chrom', 'start', 'end' columns defining
        highlight regions.
    fill : str or list of str
        Fill color(s) for highlight rectangles. If a list, one color per region.
        Default is '#FFE3E6' (light red).
    col : str or list of str
        Edge color(s) for highlight rectangles. If a list, one color per region.
        Default is 'red'.
    alpha : float
        Transparency of highlight (0=transparent, 1=opaque). Default is 0.3.
    in_background : bool
        If True, draw highlights behind track content. Default is True.
    name : str
        Track name. Default is "HighlightTrack".
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> from geneview.genometracks import (
    ...     HighlightTrack, GenomeAxisTrack, AnnotationTrack,
    ...     plot_tracks, GenomicInterval,
    ... )
    >>> ax_track = GenomeAxisTrack()
    >>> hl_regions = pd.DataFrame({
    ...     "chrom": ["chr7"], "start": [2050000], "end": [2100000]
    ... })
    >>> ht = HighlightTrack(
    ...     track_list=[ax_track],
    ...     regions=hl_regions,
    ...     fill="#FFE3E6",
    ... )
    """

    def __init__(
        self,
        track_list: Optional[List[Track]] = None,
        regions: Optional[Union[pd.DataFrame, List[GenomicInterval]]] = None,
        fill: Union[str, List[str]] = "#FFE3E6",
        col: Union[str, List[str]] = "red",
        alpha: float = 0.3,
        in_background: bool = True,
        name: str = "HighlightTrack",
        display_params: Optional[Dict[str, Any]] = None,
    ):
        dp = {
            "fill": fill,
            "col": col,
            "alpha": alpha,
        }
        if display_params:
            dp.update(display_params)

        # Convert regions to DataFrame if needed
        if isinstance(regions, list) and all(isinstance(r, GenomicInterval) for r in regions):
            regions = pd.DataFrame([
                {"chrom": r.chrom, "start": r.start, "end": r.end}
                for r in regions
            ])

        super().__init__(
            data=regions, name=name, height=1.0,
            display_params=dp,
        )

        self.track_list = track_list or []
        self.in_background = in_background

    @property
    def fill(self) -> Union[str, List[str]]:
        return self.get_param("fill", "#FFE3E6")

    @property
    def col(self) -> Union[str, List[str]]:
        return self.get_param("col", "red")

    @property
    def alpha(self) -> float:
        return self.get_param("alpha", 0.3)

    def get_region(self) -> Optional[GenomicInterval]:
        """Return region from contained tracks or highlight regions."""
        # First check contained regions
        region = super().get_region()
        if region is not None:
            return region
        # Then check contained tracks
        for track in self.track_list:
            r = track.get_region()
            if r is not None:
                return r
        return None

    def draw(self, ax, region: GenomicInterval) -> None:
        """HighlightTrack does not draw itself - it's handled by plot_tracks."""
        pass

    def draw_highlights(self, ax, region: GenomicInterval) -> None:
        """Draw the highlight rectangles on the given axes.

        This method is called by the track plot orchestrator.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw highlights on.
        region : GenomicInterval
            The current genomic region being displayed.
        """
        sub = self.subset(region)
        if sub is None or len(sub) == 0:
            return

        fill = self.fill
        col = self.col
        alpha = self.alpha

        # Get current y-limits
        ylim = ax.get_ylim()

        for row_idx, (_, row) in enumerate(sub.iterrows()):
            x_start = max(row["start"], region.start)
            x_end = min(row["end"], region.end)
            width = x_end - x_start

            if width <= 0:
                continue

            # Support per-region colors (Task 9)
            row_fill = fill[row_idx % len(fill)] if isinstance(fill, list) else fill
            row_col = col[row_idx % len(col)] if isinstance(col, list) else col

            rect = mpatches.Rectangle(
                (x_start, ylim[0]), width, ylim[1] - ylim[0],
                facecolor=row_fill, edgecolor=row_col,
                linewidth=0.5, alpha=alpha, zorder=1,
            )
            ax.add_patch(rect)

    def __repr__(self):
        return (
            f"HighlightTrack(tracks={len(self.track_list)}, "
            f"regions={len(self._data) if self._data is not None else 0})"
        )
