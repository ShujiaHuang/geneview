"""Shared utility functions for geneview CLI subcommands.

This module provides common helper functions used across CLI subcommands,
including file I/O, argument validation, and figure output handling.

Author: Shujia Huang
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for CLI
import matplotlib.pyplot as plt


def get_figure_output_format(output_path):
    """Infer the output figure format from the file extension.

    Parameters
    ----------
    output_path : str
        Path to the output figure file.

    Returns
    -------
    fmt : str
        The file extension (e.g., 'png', 'pdf', 'svg') without the leading dot.
        Defaults to 'png' if the extension is not recognized.
    """
    _, ext = os.path.splitext(output_path)
    ext = ext.lstrip(".").lower()
    supported = {"png", "pdf", "svg", "eps", "ps", "jpg", "jpeg", "tiff", "tif"}
    if ext in supported:
        return ext
    return "png"


def save_figure(fig, output_path, dpi=300):
    """Save a matplotlib figure to the specified path.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.

    output_path : str
        Path to the output file. The format is inferred from the extension.

    dpi : int, optional
        Resolution in dots per inch. Default is 300.
    """
    fmt = get_figure_output_format(output_path)
    fig.savefig(output_path, format=fmt, dpi=dpi, bbox_inches="tight")
    sys.stderr.write("[INFO] Figure saved to %s\n" % output_path)


def add_common_figure_args(parser):
    """Add common figure-related arguments to a subcommand parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subcommand parser to augment.
    """
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output figure path. Supported formats: png, pdf, svg, eps. "
             "(default: <subcommand>.png)")
    parser.add_argument(
        "--figsize",
        nargs=2, type=float, default=None, metavar=("WIDTH", "HEIGHT"),
        help="Figure size in inches: WIDTH HEIGHT (e.g., 12 4).")
    parser.add_argument(
        "--dpi",
        type=int, default=300,
        help="Figure resolution in dots per inch. (default: 300)")
    parser.add_argument(
        "--facecolor",
        default="w",
        help="Figure face color. (default: w)")


def create_figure(args, default_figsize=(12, 4)):
    """Create a matplotlib figure and axes from CLI args.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments (must contain ``figsize`` and ``facecolor``).

    default_figsize : tuple, optional
        Default figure size if ``args.figsize`` is None. Default is (12, 4).

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    figsize = args.figsize if args.figsize else default_figsize
    fig, ax = plt.subplots(figsize=figsize, facecolor=args.facecolor, edgecolor="k")
    return fig, ax


def resolve_output_path(args, default_name):
    """Resolve the output path, falling back to a default name.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments (must contain ``output``).

    default_name : str
        Default file name to use if ``args.output`` is None.

    Returns
    -------
    output_path : str
    """
    if args.output is None:
        return default_name
    return args.output
