"""CLI subcommand for creating genome track plots.

Usage::

    $ geneview tracks --region <chr:start-end> [options]

This subcommand creates genome browser-style track plots from the command line.
It supports IdeogramTrack (chromosome ideogram), AnnotationTrack, GeneRegionTrack,
DataTrack, HighlightTrack, and GenomeAxisTrack.

Examples
--------

Basic ideogram view::

    $ geneview tracks --region chr7:20000000-60000000 --ideogram -o ideogram.png

Add annotation and data tracks::

    $ geneview tracks --region chr7:26490000-26720000 \\
        --ideogram \\
        -a cpg_islands.bed \\
        -g gene_models.gtf \\
        -d coverage.bedgraph \\
        -o tracks.png

Customize track appearance::

    $ geneview tracks --region chr7:26490000-26720000 \\
        -d signal.bedgraph --data-type line --data-color blue \\
        -a features.bed --annotation-shape box \\
        --highlight regions.bed --highlight-fill yellow \\
        -o custom_tracks.png

Author: Shujia Huang
"""
import sys
import re
import argparse

from .utils import add_common_figure_args, add_style_arg, save_figure, resolve_output_path


def register(subparsers):
    """Register the ``tracks`` subcommand.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        The subparsers action object from the parent parser.
    """
    p = subparsers.add_parser(
        "tracks",
        help="Create a genome track plot from genomic data files.",
        description="Create a genome browser-style track plot with ideogram, "
                    "annotation, gene region, and data tracks. Supports BED, GFF/GTF, "
                    "bedGraph, and BigWig input formats.",
        epilog="Example: geneview tracks --region chr7:20M-60M --ideogram -a features.bed -o tracks.png",
    )

    # --- Required arguments ---
    p.add_argument(
        "--region",
        required=True,
        help="Genomic region to display in format chr:start-end (e.g., chr7:20000000-60000000). "
             "Supports M/m and K/k suffixes for megabases/kilobases.")

    # --- Ideogram and axis ---
    p.add_argument(
        "--ideogram",
        action="store_true", default=False,
        help="Add a chromosome ideogram track (auto-loads karyotype data).")
    p.add_argument(
        "--genome-build",
        choices=["hg38", "hg19"], default="hg38",
        help="Genome build for ideogram karyotype. (default: hg38)")
    p.add_argument(
        "--no-axis",
        action="store_true", default=False,
        help="Disable the default genome axis track.")

    # --- Track inputs (repeatable, order-preserving) ---
    p.add_argument(
        "-a", "--annotation",
        action="append", default=[], dest="annotations", metavar="FILE",
        help="BED file for AnnotationTrack. Can be specified multiple times. "
             "Tracks are stacked in the order flags are given.")
    p.add_argument(
        "-g", "--gene-region",
        action="append", default=[], dest="gene_regions", metavar="FILE",
        help="GFF/GTF file for GeneRegionTrack. Can be specified multiple times.")
    p.add_argument(
        "-d", "--data",
        action="append", default=[], dest="data_files", metavar="FILE",
        help="bedGraph or BigWig file for DataTrack. Can be specified multiple times.")
    p.add_argument(
        "--highlight",
        default=None, metavar="FILE",
        help="BED file defining regions to highlight across all tracks.")

    # --- Track appearance options ---
    p.add_argument(
        "--data-type",
        choices=["line", "histogram", "heatmap", "polygon", "gradient", "points"],
        default="histogram",
        help="Plot type for DataTrack. (default: histogram)")
    p.add_argument(
        "--data-color",
        default="#5B8DB8",
        help="Color for DataTrack. (default: #5B8DB8)")
    p.add_argument(
        "--annotation-shape",
        choices=["box", "arrow", "ellipse", "fixedArrow", "smallArrow"],
        default="arrow",
        help="Feature shape for AnnotationTrack. (default: arrow)")
    p.add_argument(
        "--collapse-transcripts",
        choices=["gene", "longest", "shortest", "meta"],
        default="longest",
        help="How to collapse transcripts in GeneRegionTrack. (default: longest)")
    p.add_argument(
        "--highlight-fill",
        default="#FFF3BF",
        help="Fill color for highlight regions. (default: #FFF3BF)")
    p.add_argument(
        "--highlight-alpha",
        type=float, default=0.3,
        help="Transparency for highlight regions, 0-1. (default: 0.3)")

    # --- Track name options ---
    p.add_argument(
        "--annotation-name",
        default="Annotation",
        help="Track name for annotation tracks. (default: Annotation)")
    p.add_argument(
        "--gene-region-name",
        default="GeneRegion",
        help="Track name for gene region tracks. (default: GeneRegion)")
    p.add_argument(
        "--data-name",
        default="Data",
        help="Track name for data tracks. (default: Data)")

    # --- Style ---
    add_style_arg(p)

    # --- Common figure arguments ---
    add_common_figure_args(p)

    p.set_defaults(func=run)


def _parse_region(region_str):
    """Parse a genomic region string into (chrom, start, end).

    Parameters
    ----------
    region_str : str
        Region in format chr:start-end. Supports M/m and K/k suffixes.
        Examples: chr7:20000000-60000000, chr7:20M-60M, chr7:20000K-60000K

    Returns
    -------
    tuple of (str, int, int)
        (chromosome, start, end)

    Raises
    ------
    ValueError
        If the region string cannot be parsed.
    """
    # Match pattern: chrom:start-end
    match = re.match(
        r"^(chr\w+|[\w]+):(\d+[\d.]*[MmKk]?)-(\d+[\d.]*[MmKk]?)$",
        region_str.strip()
    )
    if not match:
        raise ValueError(
            f"Cannot parse region '{region_str}'. "
            f"Expected format: chr:start-end (e.g., chr7:20000000-60000000 or chr7:20M-60M)"
        )

    chrom = match.group(1)
    start_str = match.group(2)
    end_str = match.group(3)

    def _parse_position(s):
        """Parse a position string with optional M/K suffix."""
        s = s.strip()
        if s.endswith(("M", "m")):
            return int(float(s[:-1]) * 1_000_000)
        elif s.endswith(("K", "k")):
            return int(float(s[:-1]) * 1_000)
        else:
            return int(s)

    start = _parse_position(start_str)
    end = _parse_position(end_str)

    if start >= end:
        raise ValueError(f"Start position ({start}) must be less than end position ({end}).")

    return chrom, start, end


def run(args):
    """Execute the ``tracks`` subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    from geneview.genometracks import (
        plot_tracks, GenomicInterval, GenomeAxisTrack,
        IdeogramTrack, AnnotationTrack, GeneRegionTrack,
        DataTrack, HighlightTrack,
        read_bed, read_auto,
    )
    import matplotlib.pyplot as plt

    # --- Parse region ---
    chrom, start, end = _parse_region(args.region)
    region = GenomicInterval(chrom, start, end)

    # --- Build track list ---
    tracks = []

    # Add ideogram if requested
    if args.ideogram:
        try:
            itrack = IdeogramTrack(chromosome=chrom, genome_build=args.genome_build)
            tracks.append(itrack)
        except Exception as e:
            sys.stderr.write(f"[WARNING] Could not load ideogram: {e}\n")

    # Add genome axis (unless disabled)
    if not args.no_axis:
        tracks.append(GenomeAxisTrack())

    # Load highlight regions if provided
    hl_regions = None
    if args.highlight:
        try:
            hl_regions = read_bed(args.highlight)
        except Exception as e:
            sys.stderr.write(f"[WARNING] Could not load highlight file: {e}\n")

    # Collect all data/annotation tracks to apply highlights
    data_tracks = []

    # Add annotation tracks
    for bed_file in args.annotations:
        try:
            data = read_auto(bed_file)
            track = AnnotationTrack(
                data,
                shape=args.annotation_shape,
                name=args.annotation_name,
            )
            tracks.append(track)
            data_tracks.append(track)
        except Exception as e:
            raise RuntimeError(
                "Could not load annotation file '%s': %s" % (bed_file, e)
            )

    # Add gene region tracks
    for gff_file in args.gene_regions:
        try:
            data = read_auto(gff_file)
            track = GeneRegionTrack(
                data,
                collapse_transcripts=args.collapse_transcripts,
                name=args.gene_region_name,
            )
            tracks.append(track)
            data_tracks.append(track)
        except Exception as e:
            raise RuntimeError(
                "Could not load gene region file '%s': %s" % (gff_file, e)
            )

    # Add data tracks
    for data_file in args.data_files:
        try:
            data = read_auto(data_file)
            track = DataTrack(
                data,
                type=args.data_type,
                col=args.data_color,
                name=args.data_name,
            )
            tracks.append(track)
            data_tracks.append(track)
        except Exception as e:
            raise RuntimeError(
                "Could not load data file '%s': %s" % (data_file, e)
            )

    # Wrap tracks with HighlightTrack if highlights are provided
    if hl_regions is not None and data_tracks:
        ht = HighlightTrack(
            track_list=data_tracks,
            regions=hl_regions,
            fill=args.highlight_fill,
            alpha=args.highlight_alpha,
        )
        # Replace the data_tracks in the track list with the HighlightTrack wrapper
        # The HighlightTrack will expand during plot_tracks
        for dt in data_tracks:
            if dt in tracks:
                tracks.remove(dt)
        tracks.append(ht)

    # --- Check that we have at least one track ---
    if not tracks:
        raise RuntimeError(
            "No tracks to plot. Use --ideogram, -a, -g, or -d to add tracks."
        )

    # --- Create figure and plot ---
    figsize = args.figsize if args.figsize else (14, max(3, len(tracks) * 1.5))

    try:
        axes = plot_tracks(tracks, region=region, figsize=figsize, style=args.style)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to plot tracks: {e}\n")
        return 1

    # --- Save output ---
    output_path = resolve_output_path(args, "genome_tracks.png")
    fig = plt.gcf()
    save_figure(fig, output_path, dpi=args.dpi)
    plt.close(fig)

    return 0
