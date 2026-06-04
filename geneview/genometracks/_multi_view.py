"""
Multi-view layouts for genome tracks.

Provides :func:`plot_tracks_grid` for rendering multiple independent
genomic views side-by-side (e.g. comparing two loci) and
:func:`plot_tracks_multi` for stacking tracks from different genomic
regions in a single figure.

Ported from ``genomeview``'s ``ViewRow`` (side-by-side) and ``Document``
(stacked multi-region) concepts.
"""

from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from ._base import Track, GenomicInterval
from ._track_plot import plot_tracks, _expand_tracks, _determine_chromosome, _determine_region


def plot_tracks_grid(
    views: List[List[Track]],
    regions: Optional[List[Optional[GenomicInterval]]] = None,
    columns: int = 2,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    fontsize_main: float = 14,
    show_title: bool = True,
    title_width: float = 0.08,
    style: Optional[str] = None,
    **kwargs,
) -> List:
    """Render multiple genomic views side-by-side in a grid layout.

    Each *view* is an independent list of tracks rendered in its own
    column, allowing comparison of different genomic regions or datasets.

    Parameters
    ----------
    views : list of list of Track
        Each element is a list of :class:`Track` objects representing
        one genomic view (column).
    regions : list of GenomicInterval, optional
        One region per view.  If ``None``, regions are auto-derived
        from each view's track data.
    columns : int
        Number of columns in the grid.  Default is 2.
    figsize : tuple of float, optional
        Figure size ``(width, height)`` in inches.  Auto-computed if
        ``None``.
    title : str, optional
        Main title for the entire figure.
    fontsize_main : float
        Font size for the main title.
    show_title : bool
        Whether to show per-track title panels.  Default is True.
    title_width : float
        Fraction of each column used for title panels.  Default is 0.08.
    style : str, optional
        Name of a registered plot style.
    **kwargs
        Additional display parameters applied to all tracks.

    Returns
    -------
    list of list of Axes
        Axes objects for each view, in the same order as *views*.

    Examples
    --------
    >>> from geneview.genometracks import (                # doctest: +SKIP
    ...     plot_tracks_grid, GenomeAxisTrack, AnnotationTrack,
    ...     GenomicInterval,
    ... )
    >>> view1 = [GenomeAxisTrack(), AnnotationTrack(data1)]
    >>> view2 = [GenomeAxisTrack(), AnnotationTrack(data2)]
    >>> axes = plot_tracks_grid(                           # doctest: +SKIP
    ...     [view1, view2],
    ...     regions=[GenomicInterval("chr1", 0, 1e6),
    ...              GenomicInterval("chr2", 0, 1e6)],
    ... )
    """
    if regions is None:
        regions = [None] * len(views)

    n_views = len(views)
    n_rows = (n_views + columns - 1) // columns

    # Compute figure size
    if figsize is None:
        width = 12
        # Estimate height from the tallest view
        max_tracks = max(len(v) for v in views) if views else 1
        height = 1.0 + max_tracks * 1.2 * n_rows
        figsize = (width, height)

    fig = plt.figure(figsize=figsize, facecolor="white")

    has_main_title = title is not None and title != ""
    n_grid_rows = n_rows + (1 if has_main_title else 0)

    outer_gs = gridspec.GridSpec(
        n_grid_rows, columns,
        hspace=0.15, wspace=0.15,
    )

    # Main title
    row_offset = 0
    if has_main_title:
        title_ax = fig.add_subplot(outer_gs[0, :])
        title_ax.text(0.5, 0.5, title, ha="center", va="center",
                      fontsize=fontsize_main, fontweight="bold",
                      transform=title_ax.transAxes)
        title_ax.axis("off")
        row_offset = 1

    all_axes = []

    for idx, (view_tracks, view_region) in enumerate(zip(views, regions)):
        col = idx % columns
        row = row_offset + idx // columns

        sub_gs = outer_gs[row, col].subgridspec(
            len(view_tracks), 1, hspace=0.05,
        )

        view_axes = []
        for ti, track in enumerate(view_tracks):
            ax = fig.add_subplot(sub_gs[ti])
            if view_region is not None:
                track.draw(ax, view_region)
            else:
                # Auto-derive region
                expanded, highlights = _expand_tracks([track])
                chrom = _determine_chromosome(expanded)
                region = _determine_region(expanded, chrom)
                if region is not None:
                    track.draw(ax, region)
            view_axes.append(ax)

        all_axes.append(view_axes)

    plt.subplots_adjust(
        left=0.03, right=0.97, top=0.95 if not has_main_title else 0.92,
        bottom=0.03,
    )

    return all_axes


def plot_tracks_multi(
    sections: List[Tuple[List[Track], Optional[GenomicInterval]]],
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    fontsize_main: float = 14,
    show_title: bool = True,
    title_width: float = 0.08,
    style: Optional[str] = None,
    **kwargs,
) -> List:
    """Stack multiple genomic sections (potentially different regions) vertically.

    This is analogous to ``genomeview.Document`` which can contain
    multiple ``GenomeView`` objects, each covering a different region.

    Parameters
    ----------
    sections : list of (track_list, region) tuples
        Each element is a ``(tracks, region)`` pair.  ``region`` may be
        ``None`` to auto-derive from track data.
    figsize : tuple of float, optional
        Figure size ``(width, height)`` in inches.
    title : str, optional
        Main title for the entire figure.
    fontsize_main : float
        Font size for the main title.
    show_title : bool
        Whether to show per-track title panels.
    title_width : float
        Fraction of figure width for title panels.
    style : str, optional
        Name of a registered plot style.
    **kwargs
        Additional display parameters applied to all tracks.

    Returns
    -------
    list of Axes
        All Axes objects from every section.

    Examples
    --------
    >>> from geneview.genometracks import (                # doctest: +SKIP
    ...     plot_tracks_multi, GenomeAxisTrack, AnnotationTrack,
    ...     GenomicInterval,
    ... )
    >>> section1 = ([GenomeAxisTrack(), AnnotationTrack(data1)],
    ...             GenomicInterval("chr1", 0, 1e6))
    >>> section2 = ([GenomeAxisTrack(), AnnotationTrack(data2)],
    ...             GenomicInterval("chr2", 0, 1e6))
    >>> axes = plot_tracks_multi([section1, section2])     # doctest: +SKIP
    """
    all_axes = []

    for tracks, region in sections:
        axes = plot_tracks(
            tracks,
            region=region,
            title=title if not all_axes else None,  # Title only on first section
            fontsize_main=fontsize_main,
            show_title=show_title,
            title_width=title_width,
            style=style,
            figsize=figsize,
            **kwargs,
        )
        all_axes.extend(axes)

    return all_axes
