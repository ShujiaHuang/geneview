"""CLI subcommand for creating Q-Q plots.

Usage::

    $ geneview qq -i <input_file> [options]

The input file should be a tab-delimited text file (or CSV) containing a
column of p-values. The default column name is ``P``.

Examples
--------

Basic usage::

    $ geneview qq -i gwas_results.assoc -o qq.png

Specify column name and customize appearance::

    $ geneview qq -i gwas.csv --sep "," --pv PVAL \\
        --title "My GWAS QQ" --marker "o" \\
        -o qq_custom.png

Author: Shujia Huang
"""
import sys
import pandas as pd

from .utils import add_common_figure_args, create_figure, save_figure, resolve_output_path


def register(subparsers):
    """Register the ``qq`` subcommand.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        The subparsers action object from the parent parser.
    """
    p = subparsers.add_parser(
        "qq",
        help="Create a Q-Q plot from GWAS association results.",
        description="Create a Q-Q (quantile-quantile) plot from a file "
                    "containing p-values. The x-axis shows expected quantiles "
                    "from a uniform distribution, and the y-axis shows observed "
                    "quantiles.",
        epilog="Example: geneview qq -i gwas.assoc -o qq.png",
    )

    # --- Required arguments ---
    p.add_argument(
        "-i", "--input",
        required=True, dest="input_file",
        help="Input file path. Tab-delimited by default (use --sep to change). "
             "Must contain a column of p-values.")

    # --- Column name arguments ---
    p.add_argument("--sep", default="\t",
                   help="Column separator in the input file. (default: tab)")
    p.add_argument("--pv", default="P",
                   help="Column name for p-value. (default: P)")

    # --- Plot content arguments ---
    p.add_argument("--title", default=None,
                   help="Plot title. The genomic inflation factor (lambda) will "
                        "be appended automatically. (default: None)")
    p.add_argument("--xlabel", default=None,
                   help="X-axis label. (default: auto-generated)")
    p.add_argument("--ylabel", default=None,
                   help="Y-axis label. (default: auto-generated)")
    p.add_argument("--marker", default="o",
                   help="Matplotlib marker style. (default: o)")
    p.add_argument("--color", default=None,
                   help="Color for the data points. (default: auto)")
    p.add_argument("--alpha", type=float, default=0.8,
                   help="Marker transparency, 0 (transparent) to 1 (opaque). "
                        "(default: 0.8)")
    p.add_argument("--no-logp", dest="logp", action="store_false", default=True,
                   help="Plot raw p-values instead of -log10(p).")
    p.add_argument("--ablinecolor", default="r",
                   help="Color for the y=x reference line. "
                        "Set 'none' to disable. (default: r)")

    # --- Common figure arguments ---
    add_common_figure_args(p)

    p.set_defaults(func=run)


def run(args):
    """Execute the ``qq`` subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    from geneview import qqplot

    # --- Read input ---
    df = pd.read_csv(args.input_file, sep=args.sep)

    if args.pv not in df.columns:
        raise ValueError("Column '%s' not found in input file. "
                         "Available columns: %s" % (args.pv, list(df.columns)))

    data = df[args.pv].dropna()

    # --- Handle abline color ---
    ablinecolor = args.ablinecolor
    if ablinecolor.lower() == "none":
        ablinecolor = None

    # --- Create figure ---
    fig, ax = create_figure(args, default_figsize=(6, 6))

    ax = qqplot(
        data=data,
        logp=args.logp,
        marker=args.marker,
        color=args.color,
        alpha=args.alpha,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        ablinecolor=ablinecolor,
        ax=ax,
    )

    # --- Save output ---
    output_path = resolve_output_path(args, "qq.png")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=args.dpi)
