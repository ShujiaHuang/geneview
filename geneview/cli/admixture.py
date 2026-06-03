"""CLI subcommand for creating Admixture plots.

Usage::

    $ geneview admixture -i <Q_file> -p <population_info_file> [options]

The input Q file is the standard output of ADMIXTURE (``.Q`` file), where each
row represents a sample and each column represents an ancestral component (K).
Values are space-delimited ancestry proportions.

The population info file is a plain text file with one population label per line,
corresponding to the rows of the Q file.

Examples
--------

Basic usage::

    $ geneview admixture -i output.3.Q -p population.txt -o admixture.png

Customize appearance::

    $ geneview admixture -i output.5.Q -p population.txt \\
        --palette Set1 --edgewidth 2.0 \\
        --xtick-rotation 45 --set-xticklabel-top \\
        -o admixture_K5.png

Specify group order and sample a fraction of each population::

    $ geneview admixture -i output.3.Q -p population.txt \\
        --group-order POP1 POP2 POP3 \\
        --shuffle-frac 0.5 \\
        -o admixture_sampled.png

Author: Shujia Huang
"""
import sys

from .utils import add_common_figure_args, create_figure, save_figure, resolve_output_path


def register(subparsers):
    """Register the ``admixture`` subcommand.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        The subparsers action object from the parent parser.
    """
    p = subparsers.add_parser(
        "admixture",
        help="Create an Admixture plot from ADMIXTURE .Q output.",
        description="Create a stacked bar plot (Admixture plot) from the "
                    "output of ADMIXTURE software. Requires a .Q file "
                    "(ancestry proportions) and a population info file.",
        epilog="Example: geneview admixture -i output.3.Q -p pop.txt -o admixture.png",
    )

    # --- Required arguments ---
    p.add_argument(
        "-i", "--input",
        required=True, dest="input_file",
        help="Path to the ADMIXTURE output .Q file. Space-delimited, no header. "
             "Each row is a sample; each column is an ancestral component.")
    p.add_argument(
        "-p", "--population-info",
        required=True, dest="population_info",
        help="Path to a text file with one population label per line, "
             "corresponding to the rows of the .Q file.")

    # --- Plot content arguments ---
    p.add_argument("--group-order", nargs="+", default=None,
                   help="Specify the order of populations in the plot. "
                        "If not set, uses the order found in the data.")
    p.add_argument("--palette", default="tab10",
                   help="Color palette name or comma-separated hex colors. "
                        "(default: tab10)")
    p.add_argument("--linewidth", type=float, default=1.0,
                   help="Width of vertical lines between populations. "
                        "(default: 1.0)")
    p.add_argument("--edgewidth", type=float, default=1.0,
                   help="Width of the figure frame edges. (default: 1.0)")
    p.add_argument("--ylabel", default=None,
                   help="Y-axis label. Auto-generated as 'K=<n>' if not set.")
    p.add_argument("--set-xticklabel-top", action="store_true", default=False,
                   help="Place population labels at the top of the figure.")
    p.add_argument("--xtick-rotation", type=float, default=None,
                   help="Rotation angle for x-tick labels (e.g., 45).")
    p.add_argument("--xtick-labels", nargs="+", default=None,
                   help="Custom x-tick labels (must match number of groups).")

    # --- Sampling arguments ---
    p.add_argument("--shuffle-n", type=int, default=None,
                   help="Randomly sample N individuals per population.")
    p.add_argument("--shuffle-frac", type=float, default=None,
                   help="Randomly sample a fraction (0-1) of each population.")

    # --- Common figure arguments ---
    add_common_figure_args(p)

    p.set_defaults(func=run)


def run(args):
    """Execute the ``admixture`` subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    from geneview import admixtureplot

    # --- Build shuffle kwargs ---
    shuffle_kws = {}
    if args.shuffle_n is not None:
        shuffle_kws["n"] = args.shuffle_n
    if args.shuffle_frac is not None:
        shuffle_kws["frac"] = args.shuffle_frac
    shuffle_kws = shuffle_kws if shuffle_kws else None

    # --- Build xticklabel_kws ---
    xticklabel_kws = {}
    if args.xtick_rotation is not None:
        xticklabel_kws["rotation"] = args.xtick_rotation

    # --- Parse palette ---
    palette = args.palette
    if "," in palette:
        palette = palette.split(",")

    # --- Create figure ---
    fig, ax = create_figure(args, default_figsize=(14, 2))

    ax = admixtureplot(
        data=args.input_file,
        population_info=args.population_info,
        shuffle_popsample_kws=shuffle_kws,
        group_order=args.group_order,
        linewidth=args.linewidth,
        edgewidth=args.edgewidth,
        palette=palette,
        xticklabels=args.xtick_labels,
        xticklabel_kws=xticklabel_kws if xticklabel_kws else None,
        ylabel=args.ylabel,
        ylabel_kws={"rotation": 0, "ha": "right"} if args.ylabel else None,
        set_xticklabel_top=args.set_xticklabel_top,
        ax=ax,
    )

    # --- Save output ---
    output_path = resolve_output_path(args, "admixture.png")
    save_figure(fig, output_path, dpi=args.dpi)
