"""
OverlayTrack - A container track for overlaying multiple tracks on the same panel.

Allows compositing multiple tracks (e.g. DataTracks) on a single plot area.
The first track's title and axis configuration are used. Alpha blending
is useful for transparency when overlaying.

Ported from Gviz's OverlayTrack-class.R.
"""

from typing import Optional, List, Dict, Any

from ._base import Track, GenomicInterval


class OverlayTrack(Track):
    """A container track that overlays multiple tracks on the same plot panel.

    All contained tracks are drawn on the same matplotlib Axes, allowing
    visual comparison (e.g. overlapping DataTracks with different colors).
    The first track's title and y-axis configuration are used.

    Parameters
    ----------
    track_list : list of Track
        List of Track objects to overlay on the same panel.
    name : str
        Track name for the title panel. Defaults to the first track's name.
    height : float
        Relative track height. Default is 1.5.
    display_params : dict, optional
        Additional display parameters applied to all contained tracks.

    Examples
    --------
    >>> from geneview.genometracks import (
    ...     OverlayTrack, DataTrack, GenomeAxisTrack, plot_tracks, GenomicInterval,
    ... )
    >>> import pandas as pd, numpy as np
    >>> rng = np.random.RandomState(42)
    >>> n = 50
    >>> starts = np.linspace(2_000_000, 2_050_000, n, dtype=int)
    >>> data1 = pd.DataFrame({"chrom": ["chr7"]*n, "start": starts,
    ...     "end": starts + 1000, "value": rng.randn(n).cumsum()})
    >>> data2 = pd.DataFrame({"chrom": ["chr7"]*n, "start": starts,
    ...     "end": starts + 1000, "value": rng.randn(n).cumsum()})
    >>> dt1 = DataTrack(data1, type="line", col="blue", name="Sample A")
    >>> dt2 = DataTrack(data2, type="line", col="red", name="Sample B")
    >>> ot = OverlayTrack(track_list=[dt1, dt2], name="Overlay")
    >>> _ = plot_tracks([GenomeAxisTrack(), ot],
    ...             region=GenomicInterval("chr7", 2_000_000, 2_050_000))
    """

    def __init__(
        self,
        track_list: Optional[List[Track]] = None,
        name: Optional[str] = None,
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        self.track_list = track_list or []
        # Use first track's name if not specified
        if name is None and self.track_list:
            name = self.track_list[0].name
        elif name is None:
            name = "Overlay"

        super().__init__(name=name, height=height, display_params=display_params)

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw all contained tracks on the same axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        if not self.track_list:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            return

        # Draw each track on the same axes
        for track in self.track_list:
            track.draw(ax, region)

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the region from the first contained track that has one."""
        for track in self.track_list:
            r = track.get_region()
            if r is not None:
                return r
        return None

    def __repr__(self):
        return (
            f"OverlayTrack(name='{self.name}', "
            f"tracks={len(self.track_list)})"
        )
