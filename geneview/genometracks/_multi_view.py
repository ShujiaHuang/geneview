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
from ._genome_axis import GenomeAxisTrack
from ._ideogram import IdeogramTrack
from ._track_plot import (
    plot_tracks,
    add_panel_labels,
    _expand_tracks,
    _determine_chromosome,
    _determine_region,
    _apply_style_to_tracks,
    _draw_title_panel,
)


def _draw_section(fig, cell, tracks, region, show_title, title_width):
    """Draw a list of tracks stacked vertically inside one GridSpec *cell*.

    Shared by :func:`plot_tracks_grid` (one section per grid cell) and
    :func:`plot_tracks_multi` (one section per stacked region).  Handles
    ``HighlightTrack`` expansion, per-section region auto-derivation,
    optional left-hand title panels, and x-axis labelling on the bottom
    track — mirroring :func:`plot_tracks`'s single-view layout.

    Returns the list of *data* axes (title panels excluded), so callers can
    treat the first entry as the section's anchor (e.g. for panel labels).
    """
    from ._base import _genomic_position_formatter
    from ..plotstyle import get_active_style

    expanded, highlights = _expand_tracks(list(tracks))
    for hl in highlights:
        hl._target_ids = set(id(t) for t in hl.track_list)

    if region is None:
        chrom = _determine_chromosome(expanded)
        region = _determine_region(expanded, chrom)

    n = len(expanded)
    if n == 0:
        return []

    ncol = 2 if show_title else 1
    sub = cell.subgridspec(
        n, ncol,
        width_ratios=[title_width, 1.0 - title_width] if show_title else None,
        hspace=0.05, wspace=0.02,
    )

    active = get_active_style()
    tick_fs = active.tracks_tick_fontsize if active is not None else 7
    axis_color = active.tracks_axis_color if active is not None else "darkgray"
    axis_lw = active.tracks_axis_linewidth if active is not None else 0.8

    data_axes = []
    for ti, track in enumerate(expanded):
        if show_title:
            ax_title = fig.add_subplot(sub[ti, 0])
            _draw_title_panel(ax_title, track, region)
            ax_data = fig.add_subplot(sub[ti, 1])
        else:
            ax_data = fig.add_subplot(sub[ti])

        if region is not None:
            track.draw(ax_data, region)

        is_last = ti == n - 1
        is_axis_like = isinstance(track, (GenomeAxisTrack, IdeogramTrack))
        if is_axis_like:
            pass  # these tracks render their own coordinate labels
        elif is_last and region is not None:
            span = region.end - region.start
            ax_data.xaxis.set_major_formatter(_genomic_position_formatter(span))
            ax_data.tick_params(axis="x", labelsize=tick_fs)
            ax_data.spines["bottom"].set_visible(True)
            ax_data.spines["bottom"].set_color(axis_color)
            ax_data.spines["bottom"].set_linewidth(axis_lw)
        else:
            ax_data.set_xticklabels([])

        for hl in highlights:
            if id(track) in getattr(hl, "_target_ids", set()):
                hl.draw_highlights(ax_data, region)

        data_axes.append(ax_data)

    return data_axes


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
    panel_labels: Union[bool, List[str], None] = None,
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
        Whether to show per-track title panels on the left of each view.
        Default is True.
    title_width : float
        Fraction of each column used for title panels.  Default is 0.08.
    style : str, optional
        Name of a registered plot style.
    panel_labels : bool or list of str, optional
        When truthy, draw sequential ``a/b/c`` labels on each view (one
        label per grid cell, anchored at the view's first track).  Pass a
        list of strings to use custom labels.  Label case / size follow the
        active plot style (see :func:`add_panel_labels`).  Default is None
        (no labels).
    **kwargs
        Additional display parameters applied to all tracks in every view.

    Returns
    -------
    list of list of Axes
        Data axes for each view, in the same order as *views*.

    Examples
    --------
    >>> from geneview.genometracks import (                # doctest: +SKIP
    ...     plot_tracks_grid, GenomeAxisTrack, AnnotationTrack,
    ...     GenomicInterval,
    ... )
    >>> view1 = [GenomeAxisTrack(), AnnotationTrack(data1)]  # doctest: +SKIP
    >>> view2 = [GenomeAxisTrack(), AnnotationTrack(data2)]  # doctest: +SKIP
    >>> axes = plot_tracks_grid(                             # doctest: +SKIP
    ...     [view1, view2],
    ...     regions=[GenomicInterval("chr1", 0, 1e6),
    ...              GenomicInterval("chr2", 0, 1e6)],
    ... )
    """
    if regions is None:
        regions = [None] * len(views)

    # Resolve the requested style and cascade its track-parameter overrides
    # into every view's tracks.  The figure is then created and drawn under
    # ``use_style`` so the style's rcParams (fonts, axes / tick styling)
    # apply to this grid too (previously ``style`` was silently ignored).
    from ..plotstyle import use_style as _use_style, get_style as _get_style
    resolved_style = _get_style(style) if style is not None else None
    for view_tracks in views:
        if resolved_style is not None:
            _apply_style_to_tracks(view_tracks, resolved_style)
        # Apply user display parameters to every track (was silently dropped).
        if kwargs:
            for track in view_tracks:
                track.set_params(kwargs)

    with _use_style(resolved_style):
        all_axes = _plot_tracks_grid_impl(
            views, regions, columns, figsize, title, fontsize_main,
            show_title, title_width,
        )
        if panel_labels:
            custom = panel_labels if isinstance(panel_labels, (list, tuple)) else None
            add_panel_labels(all_axes, labels=custom)
        return all_axes


def _plot_tracks_grid_impl(
    views, regions, columns, figsize, title, fontsize_main,
    show_title=True, title_width=0.08,
):
    """Build the grid figure (called inside the active-style context)."""
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
        cell = outer_gs[row, col]
        view_axes = _draw_section(
            fig, cell, view_tracks, view_region, show_title, title_width,
        )
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
    panel_labels: Union[bool, List[str], None] = None,
    **kwargs,
) -> List:
    """Stack multiple genomic sections (potentially different regions) vertically.

    All sections are drawn into a **single** figure — analogous to
    ``genomeview.Document`` which can contain multiple ``GenomeView``
    objects, each covering a different region.

    Parameters
    ----------
    sections : list of (track_list, region) tuples
        Each element is a ``(tracks, region)`` pair.  ``region`` may be
        ``None`` to auto-derive from track data.
    figsize : tuple of float, optional
        Figure size ``(width, height)`` in inches.  Auto-computed if
        ``None``.
    title : str, optional
        Main title for the entire figure (drawn once, at the top).
    fontsize_main : float
        Font size for the main title.
    show_title : bool
        Whether to show per-track title panels.  Default is True.
    title_width : float
        Fraction of figure width for title panels.  Default is 0.08.
    style : str, optional
        Name of a registered plot style.
    panel_labels : bool or list of str, optional
        When truthy, draw sequential ``a/b/c`` labels on each section
        (one label per section, anchored at its first track).  Pass a list
        of strings for custom labels.  Default is None (no labels).
    **kwargs
        Additional display parameters applied to all tracks in every
        section.

    Returns
    -------
    list of Axes
        All data axes from every section, top to bottom.

    Examples
    --------
    >>> from geneview.genometracks import (                # doctest: +SKIP
    ...     plot_tracks_multi, GenomeAxisTrack, AnnotationTrack,
    ...     GenomicInterval,
    ... )
    >>> section1 = ([GenomeAxisTrack(), AnnotationTrack(data1)],   # doctest: +SKIP
    ...             GenomicInterval("chr1", 0, 1e6))
    >>> section2 = ([GenomeAxisTrack(), AnnotationTrack(data2)],   # doctest: +SKIP
    ...             GenomicInterval("chr2", 0, 1e6))
    >>> axes = plot_tracks_multi([section1, section2])     # doctest: +SKIP
    """
    from ..plotstyle import use_style as _use_style, get_style as _get_style
    resolved_style = _get_style(style) if style is not None else None

    # Cascade style + user display parameters into every section's tracks.
    for tracks, _ in sections:
        if resolved_style is not None:
            _apply_style_to_tracks(list(tracks), resolved_style)
        if kwargs:
            for track in tracks:
                track.set_params(kwargs)

    # Per-section relative height = sum of that section's track heights.
    section_heights = []
    for tracks, _ in sections:
        expanded, _hl = _expand_tracks(list(tracks))
        section_heights.append(sum(t.height for t in expanded) or 1.0)

    if figsize is None:
        if resolved_style is not None:
            width = resolved_style.tracks_figsize_width
            height_per = resolved_style.tracks_height_per_track
        else:
            width = 12
            height_per = 1.2
        height = 1.0 + sum(section_heights) * height_per
        figsize = (width, height)

    has_main_title = title is not None and title != ""

    with _use_style(resolved_style):
        fig = plt.figure(figsize=figsize, facecolor="white")

        height_ratios = ([0.3] if has_main_title else []) + section_heights
        outer_gs = gridspec.GridSpec(
            len(height_ratios), 1, height_ratios=height_ratios, hspace=0.3,
        )

        row_offset = 0
        if has_main_title:
            title_ax = fig.add_subplot(outer_gs[0])
            title_ax.text(0.5, 0.5, title, ha="center", va="center",
                          fontsize=fontsize_main, fontweight="bold",
                          transform=title_ax.transAxes)
            title_ax.axis("off")
            row_offset = 1

        all_axes = []
        section_anchors = []
        for si, (tracks, region) in enumerate(sections):
            cell = outer_gs[row_offset + si]
            data_axes = _draw_section(
                fig, cell, tracks, region, show_title, title_width,
            )
            if data_axes:
                section_anchors.append(data_axes[0])
            all_axes.extend(data_axes)

        plt.subplots_adjust(
            left=0.03, right=0.97,
            top=0.95 if not has_main_title else 0.92, bottom=0.03,
        )

        if panel_labels:
            custom = panel_labels if isinstance(panel_labels, (list, tuple)) else None
            add_panel_labels(section_anchors, labels=custom)

    return all_axes
