"""Shared utility functions for geneview CLI subcommands.

This module provides common helper functions used across CLI subcommands,
including file I/O, argument validation, and figure output handling.

Author: Shujia Huang
"""
import os
import sys
import json

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for CLI
import matplotlib.pyplot as plt


# Valid plot-style names that can be passed via --style.
VALID_STYLES = ("geneview", "nature", "science", "cell")


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


def add_style_arg(parser):
    """Add a ``--style`` argument to a subcommand parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subcommand parser to augment.
    """
    parser.add_argument(
        "--style",
        choices=VALID_STYLES,
        default=None,
        help="Apply a built-in plot style to the figure. "
             "Choices: geneview (default), nature, science, cell. "
             "Each style configures fonts, sizes, colours, and export settings "
             "to comply with the corresponding journal's guidelines. "
             "(default: None, uses the currently active style)")


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


def add_table_input_args(parser, chrom=False, pos=False, snp=False):
    """Add tabular-input column-name arguments shared by table readers.

    Always adds ``--sep`` and ``--pv``; adds ``--chrom``, ``--pos`` and
    ``--snp`` only when requested, so ``manhattan`` (all five) and ``qq``
    (separator + p-value only) share one definition.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subcommand parser to augment.
    chrom : bool, optional
        Add a ``--chrom`` argument. Default False.
    pos : bool, optional
        Add a ``--pos`` argument. Default False.
    snp : bool, optional
        Add a ``--snp`` argument. Default False.
    """
    parser.add_argument("--sep", default="\t",
                        help="Column separator in the input file. (default: tab)")
    if chrom:
        parser.add_argument("--chrom", default="#CHROM",
                            help="Column name for chromosome. (default: #CHROM)")
    if pos:
        parser.add_argument("--pos", default="POS",
                            help="Column name for position. (default: POS)")
    parser.add_argument("--pv", default="P",
                        help="Column name for p-value. (default: P)")
    if snp:
        parser.add_argument("--snp", default="ID",
                            help="Column name for SNP identifier. (default: ID)")


def add_set_arg(parser):
    """Add a repeatable ``--set KEY=VALUE`` argument for nested plot kwargs.

    This exposes the underlying ``**kwargs``-style dictionaries (e.g.
    ``text_kws``, ``adjust_text_kws``, ``hline_kws``) that would otherwise be
    unreachable from the command line, avoiding one-flag-per-knob explosion.

    Parameters
    ----------
    parser : argparse.ArgumentParser
        The subcommand parser to augment.
    """
    parser.add_argument(
        "--set", dest="extra_kws", action="append", default=[],
        metavar="KEY=VALUE",
        help="Override/extend a nested plot kwarg as a dotted KEY=VALUE pair "
             "(repeatable). VALUE is auto-typed: int, float, true/false, "
             "none/null, a comma list (0.5,0.8), or JSON ([..]/{..}). "
             "Examples: --set adjust_text_kws.force_text=0.5,0.8 "
             "--set adjust_text_kws.lim=300 --set hline_kws.lw=2. "
             "Later --set values override the dedicated flags.")


def _coerce_set_value(raw):
    """Coerce a ``--set`` VALUE string into a native Python object.

    Rules (first match wins): ``none``/``null``->None, ``true``/``false``->bool,
    comma list->tuple of coerced items, ``[``/``{``/``(`` container->JSON, then
    int, then float, else the raw string.
    """
    s = raw.strip()
    low = s.lower()
    if low in ("none", "null"):
        return None
    if low == "true":
        return True
    if low == "false":
        return False
    if "," in s:  # comma-separated -> tuple (good for (x, y) pairs)
        return tuple(_coerce_set_value(p) for p in s.split(","))
    if s[:1] in ("[", "{", "("):
        try:
            return json.loads(s.replace("(", "[").replace(")", "]"))
        except ValueError:
            return s
    for cast in (int, float):
        try:
            return cast(s)
        except ValueError:
            continue
    return s


def parse_set_overrides(items):
    """Parse a list of ``KEY=VALUE`` strings into a nested dict.

    Parameters
    ----------
    items : list of str
        Raw ``--set`` values, e.g. ``["adjust_text_kws.force_text=0.5,0.8"]``.

    Returns
    -------
    dict
        Nested mapping, e.g. ``{"adjust_text_kws": {"force_text": (0.5, 0.8)}}``.

    Raises
    ------
    ValueError
        If an item is not in ``KEY=VALUE`` form or has an empty key.
    """
    overrides = {}
    for item in items:
        if "=" not in item:
            raise ValueError("Invalid --set value %r; expected KEY=VALUE." % item)
        key, raw = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError("Invalid --set value %r; empty KEY." % item)
        value = _coerce_set_value(raw)
        parts = key.split(".")
        node = overrides
        for part in parts[:-1]:
            nxt = node.get(part)
            if not isinstance(nxt, dict):
                nxt = {}
                node[part] = nxt
            node = nxt
        node[parts[-1]] = value
    return overrides


def deep_merge(base, override):
    """Recursively merge ``override`` into ``base`` (in place) and return it.

    Nested dicts are merged key-by-key; any non-dict value (or a dict whose
    counterpart is not a dict) overwrites the base value.
    """
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_merge(base[key], value)
        else:
            base[key] = value
    return base

