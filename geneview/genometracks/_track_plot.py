"""
plot_tracks - Main orchestrator for rendering genome tracks.

Provides the ``plot_tracks()`` function which takes a list of Track objects
and renders them as a vertically stacked genome browser view, similar to
Gviz's plotTracks() function and UCSC Genome Browser output.

Uses matplotlib's GridSpec for layout management.

Ported from Gviz's plotTracks.R.
"""

from typing import Optional, List, Dict, Any, Union, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch

from ._base import Track, GenomicInterval
from ._highlight import HighlightTrack
from ._genome_axis import GenomeAxisTrack
from ._ideogram import IdeogramTrack


def _apply_style_to_tracks(tracks, style):
    """Apply plotstyle track-parameter overrides to a list of tracks.

    Each track's display parameters are updated with the style's
    ``to_track_params()`` overrides, but only for keys the track has
    not already customised (i.e. the style is a *floor*, not a ceiling).
    """
    if style is None:
        return
    overrides = style.to_track_params()
    for track in tracks:
        for key, value in overrides.items():
            # Only override if the track is still using its class default
            # (i.e. the user did not explicitly set this param).
            track.set_param(key, value)


def plot_tracks(
    track_list: Union[Track, List[Track]],
    region: Optional[GenomicInterval] = None,
    sizes: Optional[List[float]] = None,
    title: Optional[str] = None,
    title_width: float = 0.08,
    extend_left: Union[int, float] = 0,
    extend_right: Union[int, float] = 0,
    figsize: Optional[Tuple[float, float]] = None,
    ax: Optional[Any] = None,
    main: Optional[str] = None,
    fontsize_main: float = 14,
    chromosome: Optional[str] = None,
    show_title: bool = True,
    reverse_strand: bool = False,
    style: Optional[str] = None,
    cex: Optional[float] = None,
    add: bool = False,
    ylim: Optional[Tuple[float, float]] = None,
    scheme: Optional[str] = None,
    panel_only: bool = False,
    margin: float = 6,
    inner_margin: float = 3,
    **kwargs,
) -> List:
    """Plot one or more genome tracks in a vertically stacked layout.

    This is the main entry point for rendering genome browser-style plots.
    Each track is drawn in its own panel, stacked vertically and sharing
    the same genomic x-axis coordinates. The output resembles the UCSC
    Genome Browser display.

    Parameters
    ----------
    track_list : Track or list of Track
        A single Track or list of Track objects to plot.
    region : GenomicInterval, optional
        The genomic region to display. If None, auto-derived from tracks.
    sizes : list of float, optional
        Relative vertical sizes for each track. If None, uses each track's
        ``height`` attribute.
    title : str, optional
        Main title for the entire plot.
    title_width : float
        Fraction of figure width used for track title panels. Default is 0.08.
    extend_left : int or float
        Extend the plotting range to the left. If a float between -1 and 1
        (exclusive), interpreted as a fraction of the region span.
    extend_right : int or float
        Extend the plotting range to the right. Same fractional semantics.
    show_title : bool
        Whether to show track title panels on the left. Default is True.
    reverse_strand : bool
        If True, reverse the x-axis so 3' is on the left. Default is False.
    figsize : tuple of float, optional
        Figure size (width, height) in inches. Auto-computed if None.
    ax : matplotlib Axes, optional
        If provided, plot into this single axes (no title panels or stacking).
    main : str, optional
        Alias for ``title`` parameter.
    fontsize_main : float
        Font size for the main title. Default is 14.
    chromosome : str, optional
        Force a specific chromosome for all tracks.
    style : str, optional
        Name of a registered plot style (e.g. ``"nature"``, ``"science"``,
        ``"cell"``).  When provided, the style's font sizes, colours, figure
        dimensions, and track-panel settings are applied to all tracks and
        the overall layout.  ``None`` (default) uses the currently active
        matplotlib rcParams.
    cex : float, optional
        Global font expansion factor applied to all track text elements.
        Multiplies the track's existing fontsize.  Default is None (no change).
    add : bool
        If True, plot into an existing axes without creating a new figure.
        Requires ``ax`` to be provided.  Default is False.
    ylim : tuple of float, optional
        Global y-axis limits applied to all DataTrack panels.
    scheme : str, optional
        Name of a color scheme to apply to GeneRegionTrack and
        AnnotationTrack objects (e.g. ``"genes"``, ``"transcripts"``).
    panel_only : bool
        If True, draw only the data panels without title panels and
        without creating a new figure. Useful for embedding tracks in
        larger matplotlib layouts. Implies ``add=True``. Default is False.
    margin : float
        Pixel margin around the entire plot. Translated to
        ``subplots_adjust`` values. Default is 6 (Gviz-compatible).
    inner_margin : float
        Inner margin between title and data panels in pixels.
        Controls ``wspace`` in the GridSpec. Default is 3.
    **kwargs
        Additional display parameters applied to all tracks.

    Returns
    -------
    list
        List of matplotlib Axes objects, one per track.

    Examples
    --------
    >>> import pandas as pd
    >>> from geneview.genometracks import (
    ...     plot_tracks, GenomeAxisTrack, AnnotationTrack,
    ...     GeneRegionTrack, DataTrack, GenomicInterval,
    ... )
    >>> data = pd.DataFrame({
    ...     "chrom": ["chr7"] * 4,
    ...     "start": [2000000, 2070000, 2100000, 2160000],
    ...     "end":   [2050000, 2130000, 2150000, 2170000],
    ...     "strand": ["-", "+", "-", "-"],
    ...     "name": ["feat1", "feat2", "feat3", "feat4"],
    ... })
    >>> ax_track = GenomeAxisTrack()
    >>> ann_track = AnnotationTrack(data)
    >>> region = GenomicInterval("chr7", 2000000, 2200000)
    >>> axes = plot_tracks([ax_track, ann_track], region=region)

    Plot with custom figure size and title:

    >>> axes = plot_tracks(
    ...     [ax_track, ann_track],
    ...     region=region,
    ...     figsize=(12, 6),
    ...     title="My Genomic Region",
    ... )
    """
    # Normalize input
    if isinstance(track_list, Track):
        track_list = [track_list]

    if main is not None and title is None:
        title = main

    # Resolve the style (if requested)
    from ..plotstyle import use_style as _use_style, get_style as _get_style
    resolved_style = None
    if style is not None:
        resolved_style = _get_style(style)

    # Expand HighlightTrack objects
    expanded_tracks, highlights = _expand_tracks(track_list)

    if not expanded_tracks:
        raise ValueError("No tracks to plot.")

    # Build per-highlight target-panel mapping (Bug 1 fix)
    for hl in highlights:
        hl._target_ids = set(id(t) for t in hl.track_list)

    # Apply style-level track parameter overrides (before user kwargs)
    _apply_style_to_tracks(expanded_tracks, resolved_style)

    # Apply global cex factor to all tracks
    if cex is not None:
        for track in expanded_tracks:
            current_fs = track.get_param("fontsize", 10)
            track.set_param("fontsize", current_fs * cex)

    # Apply color scheme to applicable tracks
    if scheme is not None:
        from ._schemes import apply_scheme
        for track in expanded_tracks:
            from ._gene_region import GeneRegionTrack
            from ._annotation import AnnotationTrack
            if isinstance(track, (GeneRegionTrack, AnnotationTrack)):
                apply_scheme(track, scheme)

    # Apply global ylim to DataTrack objects
    if ylim is not None:
        from ._data_track import DataTrack
        for track in expanded_tracks:
            if isinstance(track, DataTrack):
                track._ylim = ylim

    # Apply global display parameters to all tracks
    if kwargs:
        for track in expanded_tracks:
            track.set_params(kwargs)

    # Determine chromosome
    if chromosome is not None:
        target_chrom = chromosome
    else:
        target_chrom = _determine_chromosome(expanded_tracks)

    # Determine plotting region
    if region is None:
        region = _determine_region(expanded_tracks, target_chrom)
        if region is None:
            raise ValueError(
                "Cannot determine plotting region. Please provide a 'region' "
                "parameter or ensure tracks contain data."
            )

    # Apply extensions (support fractional values like Gviz)
    span = region.end - region.start
    ext_l = int(extend_left * span) if isinstance(extend_left, float) and -1 < extend_left < 1 else int(extend_left)
    ext_r = int(extend_right * span) if isinstance(extend_right, float) and -1 < extend_right < 1 else int(extend_right)
    region = region.extend(left=ext_l, right=ext_r)

    # Determine sizes
    if sizes is None:
        sizes = [t.height for t in expanded_tracks]
    total_size = sum(sizes)
    norm_sizes = [s / total_size for s in sizes]

    # Simple mode: plotting into a provided axes or add mode
    if ax is not None or add or panel_only:
        if ax is None:
            ax = plt.gca()
        axes_list = _plot_single_ax(expanded_tracks, highlights, ax, region)
        if reverse_strand:
            for a in axes_list:
                a.invert_xaxis()
        return axes_list

    # Full layout mode with title panels
    axes_list = _plot_full_layout(
        expanded_tracks, highlights, region, norm_sizes,
        title_width, title, figsize, fontsize_main,
        show_title=show_title,
        style=resolved_style,
        margin=margin,
        inner_margin=inner_margin,
    )
    if reverse_strand:
        for a in axes_list:
            a.invert_xaxis()
    return axes_list


def _expand_tracks(track_list: List[Track]) -> Tuple[List[Track], List[HighlightTrack]]:
    """Expand HighlightTrack objects, returning flat track list and highlights."""
    expanded = []
    highlights = []

    for track in track_list:
        if isinstance(track, HighlightTrack):
            # Extract contained tracks and record the highlight
            highlights.append(track)
            if track.track_list:
                expanded.extend(track.track_list)
            # If no contained tracks, the highlight still needs a region reference
        else:
            expanded.append(track)

    return expanded, highlights


def _determine_chromosome(tracks: List[Track]) -> Optional[str]:
    """Determine the chromosome from the track list."""
    for track in tracks:
        region = track.get_region()
        if region is not None:
            return region.chrom
    return None


def _determine_region(tracks: List[Track], chrom: Optional[str]) -> Optional[GenomicInterval]:
    """Determine the plotting region from track data."""
    if chrom is None:
        return None

    min_start = None
    max_end = None

    for track in tracks:
        region = track.get_region()
        if region is not None and region.chrom == chrom:
            if min_start is None or region.start < min_start:
                min_start = region.start
            if max_end is None or region.end > max_end:
                max_end = region.end

    if min_start is not None and max_end is not None:
        return GenomicInterval(chrom=chrom, start=min_start, end=max_end)
    return None


def _plot_single_ax(
    tracks: List[Track],
    highlights: List[HighlightTrack],
    ax,
    region: GenomicInterval,
) -> List:
    """Plot all tracks into a single axes (simple mode)."""
    track_id_to_ax = {id(t): ax for t in tracks}

    for track in tracks:
        track.draw(ax, region)

    # Apply highlights only to targeted track panels (Bug 1 fix)
    for hl in highlights:
        target_ids = getattr(hl, '_target_ids', set())
        if any(tid in track_id_to_ax for tid in target_ids):
            hl.draw_highlights(track_id_to_ax[next(iter(target_ids))], region)

    return [ax]


def _plot_full_layout(
    tracks: List[Track],
    highlights: List[HighlightTrack],
    region: GenomicInterval,
    sizes: List[float],
    title_width: float,
    title: Optional[str],
    figsize: Optional[Tuple[float, float]],
    fontsize_main: float,
    show_title: bool = True,
    style=None,
    margin: float = 6,
    inner_margin: float = 3,
) -> List:
    """Plot tracks with full layout: title panels + data panels."""
    n_tracks = len(tracks)

    # Build track-id → panel-index mapping (Bug 1 fix)
    track_id_to_idx = {id(t): i for i, t in enumerate(tracks)}

    # Compute figure size
    if figsize is None:
        # Use style-specific dimensions when available
        if style is not None:
            width = style.tracks_figsize_width
            height_per_track = style.tracks_height_per_track
        else:
            width = 12
            height_per_track = 1.2
        base_height = 1.0  # For margins/title
        height = base_height + sum(s * height_per_track for s in sizes)
        figsize = (width, height)

    fig = plt.figure(figsize=figsize, facecolor="white")

    # Create GridSpec: 2 columns (title, data), n_tracks rows
    has_title = show_title and title is not None and title != ""
    n_rows = n_tracks + (1 if has_title else 0)

    # When show_title is False, collapse the title column
    if show_title:
        width_ratios = [title_width, 1.0 - title_width]
    else:
        width_ratios = [0, 1.0]

    gs = gridspec.GridSpec(
        n_rows, 2,
        width_ratios=width_ratios,
        height_ratios=([0.3] if has_title else []) + list(sizes),
        hspace=0.05,
        wspace=max(0, inner_margin / 100.0) if show_title else 0,
    )

    axes_list = []
    title_axes = []

    row_offset = 1 if has_title else 0

    # Add main title if specified
    if has_title:
        title_ax = fig.add_subplot(gs[0, :])
        title_ax.text(0.5, 0.5, title, ha="center", va="center",
                      fontsize=fontsize_main, fontweight="bold",
                      transform=title_ax.transAxes)
        title_ax.axis("off")

    # Create track panels
    for i, track in enumerate(tracks):
        row = row_offset + i

        # Title panel (left)
        ax_title = fig.add_subplot(gs[row, 0])
        title_axes.append(ax_title)

        # Data panel (right)
        ax_data = fig.add_subplot(gs[row, 1])
        axes_list.append(ax_data)

        # Draw track content
        track.draw(ax_data, region)

        # Draw title panel (skip if show_title is False)
        if show_title:
            _draw_title_panel(ax_title, track, region)
        else:
            ax_title.axis("off")

        # Apply highlights only to targeted panels (Bug 1 fix)
        for hl in highlights:
            target_ids = getattr(hl, '_target_ids', set())
            if id(track) in target_ids:
                hl.draw_highlights(ax_data, region)

        # Share x-axis between tracks (except the last one shows labels)
        is_last = (i == n_tracks - 1)
        is_genome_axis = isinstance(track, GenomeAxisTrack)
        is_ideogram = isinstance(track, IdeogramTrack)

        # Use style-specific tick font size when available
        tick_fs = style.tracks_tick_fontsize if style is not None else 7

        if not is_last and not is_genome_axis and not is_ideogram:
            ax_data.set_xticklabels([])
        elif is_last and not is_genome_axis and not is_ideogram:
            # Show x-axis with genomic position labels
            span = region.end - region.start
            from ._base import _genomic_position_formatter
            ax_data.xaxis.set_major_formatter(_genomic_position_formatter(span))
            ax_data.tick_params(axis="x", labelsize=tick_fs)
            ax_data.spines["bottom"].set_visible(True)
            ax_data.spines["bottom"].set_color(
                style.tracks_axis_color if style is not None else "darkgray"
            )
            ax_data.spines["bottom"].set_linewidth(
                style.tracks_axis_linewidth if style is not None else 0.8
            )
        else:
            ax_data.set_xticklabels([])

    # Convert pixel margin to figure fraction
    dpi = fig.dpi
    margin_frac_h = margin / (figsize[1] * dpi)
    margin_frac_w = margin / (figsize[0] * dpi)

    plt.subplots_adjust(
        left=max(0.01, margin_frac_w),
        right=min(0.99, 1.0 - margin_frac_w),
        top=min(0.97 if not has_title else 0.94, 1.0 - margin_frac_h),
        bottom=max(0.03, margin_frac_h),
    )

    return axes_list


def _draw_title_panel(ax, track: Track, region: GenomicInterval) -> None:
    """Draw the title panel for a track (left side label).

    Mimics Gviz's title panel: lightgray background, white bold text
    rotated 90°, transparent border by default.
    """
    show_title = track.get_param("show_title", True)
    if not show_title:
        ax.axis("off")
        return

    bg_color = track.get_param("background_title", "#D3D3D3")
    text_color = track.get_param("col_title", "white")
    fontsize = track.get_param("fontsize_title", 10)
    fontface = track.get_param("fontface_title", "bold")
    rotation = track.get_param("rotation_title", 90)

    # Background
    ax.set_facecolor(bg_color)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Track name text
    weight = "bold" if fontface == "bold" else "normal"
    style = "italic" if fontface == "italic" else "normal"

    if rotation == 90:
        ax.text(0.5, 0.5, track.name,
                ha="center", va="center", rotation=90,
                fontsize=fontsize, color=text_color,
                fontweight=weight, fontstyle=style,
                transform=ax.transAxes)
    elif rotation != 0:
        ax.text(0.5, 0.5, track.name,
                ha="center", va="center", rotation=rotation,
                fontsize=fontsize, color=text_color,
                fontweight=weight, fontstyle=style,
                transform=ax.transAxes)
    else:
        ax.text(0.1, 0.5, track.name,
                ha="left", va="center",
                fontsize=fontsize, color=text_color,
                fontweight=weight, fontstyle=style,
                transform=ax.transAxes)

    # Border — Gviz default is transparent (no visible border)
    border_color = track.get_param("col_border_title", "transparent")
    # Matplotlib doesn't understand "transparent"; convert to "none"
    if border_color == "transparent":
        border_color = "none"
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(border_color)
        spine.set_linewidth(0.5)

    ax.set_xticks([])
    ax.set_yticks([])
