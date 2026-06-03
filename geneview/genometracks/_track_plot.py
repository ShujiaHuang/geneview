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
    **kwargs
        Additional display parameters applied to all tracks.

    Returns
    -------
    list
        List of matplotlib Axes objects, one per track.

    Examples
    --------
    >>> from geneview.genometracks import (
    ...     plot_tracks, GenomeAxisTrack, AnnotationTrack,
    ...     GeneRegionTrack, DataTrack, GenomicInterval,
    ... )
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

    # Expand HighlightTrack objects
    expanded_tracks, highlights = _expand_tracks(track_list)

    if not expanded_tracks:
        raise ValueError("No tracks to plot.")

    # Build per-highlight target-panel mapping (Bug 1 fix)
    for hl in highlights:
        hl._target_ids = set(id(t) for t in hl.track_list)

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

    # Simple mode: plotting into a provided axes
    if ax is not None:
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
) -> List:
    """Plot tracks with full layout: title panels + data panels."""
    n_tracks = len(tracks)

    # Build track-id → panel-index mapping (Bug 1 fix)
    track_id_to_idx = {id(t): i for i, t in enumerate(tracks)}

    # Compute figure size
    if figsize is None:
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
        wspace=0.02 if show_title else 0,
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

        if not is_last and not is_genome_axis:
            ax_data.set_xticklabels([])
        elif is_last and not is_genome_axis:
            # Show x-axis with genomic position labels
            span = region.end - region.start
            from ._base import _genomic_position_formatter
            ax_data.xaxis.set_major_formatter(_genomic_position_formatter(span))
            ax_data.tick_params(axis="x", labelsize=7)
            ax_data.spines["bottom"].set_visible(True)
            ax_data.spines["bottom"].set_color("#666666")
            ax_data.spines["bottom"].set_linewidth(0.5)
        else:
            ax_data.set_xticklabels([])

    plt.subplots_adjust(
        left=0.01, right=0.99,
        top=0.97 if not has_title else 0.94,
        bottom=0.03,
    )

    return axes_list


def _draw_title_panel(ax, track: Track, region: GenomicInterval) -> None:
    """Draw the title panel for a track (left side label)."""
    show_title = track.get_param("show_title", True)
    if not show_title:
        ax.axis("off")
        return

    bg_color = track.get_param("background_title", "#E8E8E8")
    text_color = track.get_param("col_title", "#333333")
    fontsize = track.get_param("fontsize_title", 9)
    fontface = track.get_param("fontface_title", "bold")
    rotation = track.get_param("rotation_title", 0)

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

    # Border
    border_color = track.get_param("col_border_title", "#CCCCCC")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_color(border_color)
    ax.spines["left"].set_linewidth(0.5)

    ax.set_xticks([])
    ax.set_yticks([])
