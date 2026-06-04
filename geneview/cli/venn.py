"""CLI subcommand for creating Venn diagrams.

Usage::

    $ geneview venn -i <file1> <file2> [<file3> ...] [options]

Each input file should contain one identifier per line (e.g., gene names,
variant IDs). The command supports 2 to 6 input files.

Examples
--------

Compare two gene lists::

    $ geneview venn -i genes_A.txt genes_B.txt -o venn2.png

Compare three datasets with custom names and colors::

    $ geneview venn -i DEG_list1.txt DEG_list2.txt DEG_list3.txt \\
        --names "Study A" "Study B" "Study C" \\
        --palette plasma \\
        --fmt "{size}\\n({percentage:.0f}%)" \\
        -o venn3.png

Author: Shujia Huang
"""
import os
import sys

from .utils import add_common_figure_args, add_style_arg, create_figure, save_figure, resolve_output_path


def _read_set_from_file(filepath):
    """Read a set of identifiers from a text file (one per line).

    Parameters
    ----------
    filepath : str
        Path to the input file.

    Returns
    -------
    data : set
        A set of stripped, non-empty line strings.
    """
    with open(filepath, "r") as f:
        return set(line.strip() for line in f if line.strip())


def register(subparsers):
    """Register the ``venn`` subcommand.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        The subparsers action object from the parent parser.
    """
    p = subparsers.add_parser(
        "venn",
        help="Create a Venn diagram from 2-6 input files.",
        description="Create a Venn diagram by comparing sets of identifiers "
                    "from 2 to 6 input files. Each file should contain one "
                    "identifier per line.",
        epilog="Example: geneview venn -i genes_A.txt genes_B.txt genes_C.txt -o venn.png",
    )

    # --- Required arguments ---
    p.add_argument(
        "-i", "--input",
        required=True, nargs="+", dest="input_files",
        help="Input file paths (2-6 files). Each file should contain one "
             "identifier per line.")

    # --- Plot content arguments ---
    p.add_argument("--names", nargs="+", default=None,
                   help="Labels for each dataset. If not provided, file names "
                        "(without extension) will be used.")
    p.add_argument("--fmt", default="{size}",
                   help="Format string for petal labels. Supports {size}, "
                        "{percentage}, and {logic}. (default: {size})")
    p.add_argument("--palette", default="viridis",
                   help="Color palette name (e.g., viridis, plasma, Set1) or "
                        "comma-separated hex colors. (default: viridis)")
    p.add_argument("--alpha", type=float, default=0.4,
                   help="Alpha blending for petal colors, 0-1. (default: 0.4)")
    p.add_argument("--fontsize", type=int, default=14,
                   help="Font size for petal labels. (default: 14)")
    p.add_argument("--legend-use-petal-color", action="store_true", default=False,
                   help="Color legend text with the same color as petals.")
    p.add_argument("--legend-loc", default=None,
                   help="Legend location (e.g., 'upper left', 'lower right'). "
                        "If not set, labels are placed around the diagram.")

    # --- Style ---
    add_style_arg(p)

    # --- Common figure arguments ---
    add_common_figure_args(p)

    p.set_defaults(func=run)


def run(args):
    """Execute the ``venn`` subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    from geneview import venn

    # --- Validate input ---
    n_files = len(args.input_files)
    if n_files < 2 or n_files > 6:
        raise ValueError("Venn diagram requires 2 to 6 input files, got %d." % n_files)

    for f in args.input_files:
        if not os.path.isfile(f):
            raise FileNotFoundError("Input file not found: %s" % f)

    # --- Read data ---
    datasets = {}
    names = args.names if args.names else []

    for i, filepath in enumerate(args.input_files):
        dataset = _read_set_from_file(filepath)
        if args.names and i < len(args.names):
            name = args.names[i]
        else:
            name = os.path.splitext(os.path.basename(filepath))[0]
        datasets[name] = dataset

    if len(datasets) < 2:
        raise ValueError("At least 2 distinct datasets are required.")

    # --- Parse palette ---
    palette = args.palette
    if "," in palette:
        palette = palette.split(",")

    # --- Create figure ---
    fig, ax = create_figure(args, default_figsize=(7, 7))

    ax = venn(
        data=datasets,
        fmt=args.fmt,
        palette=palette,
        alpha=args.alpha,
        fontsize=args.fontsize,
        legend_use_petal_color=args.legend_use_petal_color,
        legend_loc=args.legend_loc,
        style=args.style,
        ax=ax,
    )

    # --- Save output ---
    output_path = resolve_output_path(args, "venn.png")
    save_figure(fig, output_path, dpi=args.dpi)
