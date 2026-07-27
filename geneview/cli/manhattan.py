"""CLI subcommand for creating Manhattan plots.

Usage::

    $ geneview manhattan -i <input_file> [options]

The input file should be a tab-delimited text file (or CSV) with at least
three columns: chromosome, position, and p-value. The default column names
follow the PLINK2.x convention (``#CHROM``, ``POS``, ``P``).

Examples
--------

Basic usage::

    $ geneview manhattan -i gwas_results.assoc -o manhattan.png

Specify column names and add significance annotations::

    $ geneview manhattan -i gwas.csv --sep "," \\
        --chrom CHROM --pos BP --pv PVAL \\
        --sign-marker-p 1e-6 --annotate-topsnp \\
        -o manhattan_annotated.png

Customise annotation labels and layout::

    $ geneview manhattan -i gwas_results.assoc \\
        --sign-marker-p 1e-6 --annotate-topsnp \\
        --annotate-fmt '{snp}\nP={p:.1e}' \\
        --annotate-layout lane --annotate-rotation 90 \\
        -o manhattan_lane.png

Plot only chromosome 8 with a custom title::

    $ geneview manhattan -i gwas_results.assoc \\
        --chr chr8 --title "Chromosome 8 GWAS" \\
        -o manhattan_chr8.png

Author: Shujia Huang
"""
import sys
import pandas as pd

from .utils import (add_common_figure_args, add_style_arg, add_table_input_args,
                    add_set_arg, create_figure, save_figure, resolve_output_path,
                    parse_set_overrides, deep_merge)


def register(subparsers):
    """Register the ``manhattan`` subcommand.

    Parameters
    ----------
    subparsers : argparse._SubParsersAction
        The subparsers action object from the parent parser.
    """
    p = subparsers.add_parser(
        "manhattan",
        help="Create a Manhattan plot from GWAS association results.",
        description="Create a Manhattan plot from PLINK assoc output or any "
                    "tabular file with chromosome, position, and p-value columns.",
        epilog="Example: geneview manhattan -i gwas.assoc -o manhattan.png",
    )

    # --- Required arguments ---
    p.add_argument(
        "-i", "--input",
        required=True, dest="input_file",
        help="Input file path. Tab-delimited by default (use --sep to change). "
             "Must contain columns for chromosome, position, and p-value.")

    # --- Column name arguments ---
    add_table_input_args(p, chrom=True, pos=True, snp=True)

    # --- Plot content arguments ---
    p.add_argument("--title", default=None,
                   help="Plot title. (default: None)")
    p.add_argument("--xlabel", default="Chromosome",
                   help="X-axis label. (default: Chromosome)")
    p.add_argument("--ylabel", default=r"$-log_{10}{(P)}$",
                   help="Y-axis label. (default: $-log_{10}{(P)})$)")
    p.add_argument("--marker", default=".",
                   help="Matplotlib marker style. (default: .)")
    p.add_argument("--color", default="#3B5488,#53BBD5",
                   help="Comma-separated colors for chromosomes. "
                        "(default: #3B5488,#53BBD5)")
    p.add_argument("--alpha", type=float, default=0.8,
                   help="Marker transparency, 0 (transparent) to 1 (opaque). "
                        "(default: 0.8)")
    p.add_argument("--no-logp", dest="logp", action="store_false", default=True,
                   help="Plot raw p-values instead of -log10(p).")

    # --- Chromosome selection ---
    p.add_argument("--chr", dest="chromosome", default=None,
                   help="Plot only a specific chromosome (e.g., chr8). "
                        "Cannot be used together with --xtick-labels.")
    p.add_argument("--xtick-labels", nargs="+", default=None,
                   help="Only show these chromosome labels on the x-axis. "
                        "Cannot be used together with --chr.")
    p.add_argument("--xtick-rotation", type=float, default=None,
                   help="Rotation angle for x-tick labels (e.g., 45).")

    # --- Significance lines and markers ---
    p.add_argument("--suggestiveline", type=float, default=1e-5,
                   help="Suggestive significance threshold. Set 0 to disable. "
                        "(default: 1e-5)")
    p.add_argument("--genomewideline", type=float, default=5e-8,
                   help="Genome-wide significance threshold. Set 0 to disable. "
                        "(default: 5e-8)")
    p.add_argument("--sign-line-colors", default="#D62728,#2CA02C",
                   help="Comma-separated colors for suggestive and genome-wide lines. "
                        "(default: #D62728,#2CA02C)")
    p.add_argument("--hline-linestyle", default="--",
                   help="Line style for significance lines. (default: --)")
    p.add_argument("--hline-lw", type=float, default=1.3,
                   help="Line width for significance lines. (default: 1.3)")

    # --- Annotation ---
    p.add_argument("--sign-marker-p", type=float, default=None,
                   help="P-value threshold for highlighting significant SNPs.")
    p.add_argument("--sign-marker-color", default="r",
                   help="Color for significant SNP markers. (default: r)")
    p.add_argument("--annotate-topsnp", action="store_true", default=False,
                   help="Annotate the top SNP in each significant locus.")
    p.add_argument("--ld-block-size", type=int, default=50000,
                   help="LD block size (bp) for top SNP annotation. "
                        "(default: 50000)")
    p.add_argument("--annotate-fmt", default=None,
                   help="Label content for annotations. A Python format string "
                        "over the fields {snp}, {chrom}, {pos}, {p}, {log10p} "
                        r"(e.g. '{snp}\nP={p:.1e}'). If omitted, the SNP id is used.")
    p.add_argument("--annotate-layout", choices=["repel", "lane"], default="repel",
                   help="Annotation placement strategy. 'repel' iteratively pushes "
                        "overlapping labels apart; 'lane' lays them out tidily along "
                        "the top with leader lines. (default: repel)")
    p.add_argument("--annotate-rotation", type=float, default=None,
                   help="Rotation angle (degrees) for annotation label text "
                        "(e.g. 90 for vertical). (default: horizontal / layout default)")
    p.add_argument("--annotate-color", default=None,
                   help="Color for annotation label text. (default: matplotlib default)")
    p.add_argument("--no-annotate-arrow", action="store_true", default=False,
                   help="Do not draw connecting arrows/leader lines to annotations.")
    p.add_argument("--text-fontsize", type=int, default=12,
                   help="Font size for SNP text annotations. (default: 12)")

    # --- Nested-kwarg overrides (e.g. adjust_text_kws) ---
    add_set_arg(p)

    # --- Style ---
    add_style_arg(p)

    # --- Common figure arguments ---
    add_common_figure_args(p)

    p.set_defaults(func=run)


def run(args):
    """Execute the ``manhattan`` subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    from geneview import manhattanplot

    # --- Read input ---
    df = pd.read_csv(args.input_file, sep=args.sep)

    # --- Validate the requested single chromosome exists (clear error) ---
    if args.chromosome is not None:
        if args.chrom not in df.columns:
            raise ValueError("Chromosome column '%s' not found in input file. "
                             "Available columns: %s" % (args.chrom, list(df.columns)))
        chrom_values = set(df[args.chrom].astype(str))
        if str(args.chromosome) not in chrom_values:
            raise ValueError(
                "Chromosome '%s' (from --chr) is not present in the input file. "
                "Available chromosomes: %s"
                % (args.chromosome, ", ".join(sorted(chrom_values))))

    # --- Build kwargs ---
    xtick_label_set = set(args.xtick_labels) if args.xtick_labels else None
    suggestiveline = args.suggestiveline if args.suggestiveline > 0 else None
    genomewideline = args.genomewideline if args.genomewideline > 0 else None

    xticklabel_kws = {}
    if args.xtick_rotation is not None:
        xticklabel_kws["rotation"] = args.xtick_rotation

    hline_kws = {"linestyle": args.hline_linestyle, "lw": args.hline_lw}

    # --- Decode escape sequences in --annotate-fmt ---
    # The shell passes a literal "\n"/"\t" (backslash + letter); translate
    # these into real newlines/tabs so multi-line labels render correctly.
    annotate_fmt = args.annotate_fmt
    if annotate_fmt is not None:
        annotate_fmt = annotate_fmt.replace("\\n", "\n").replace("\\t", "\t")

    # --- Build annotation label styling (text_kws) ---
    text_kws = {"fontsize": args.text_fontsize}
    if args.annotate_rotation is not None:
        text_kws["rotation"] = args.annotate_rotation
    if args.annotate_color is not None:
        text_kws["color"] = args.annotate_color
    if args.no_annotate_arrow:
        text_kws["arrowprops"] = None
    else:
        text_kws["arrowprops"] = dict(arrowstyle="-", color="k", alpha=0.6)

    # --- Apply generic --set overrides (deep-merged; later --set wins) ---
    overrides = parse_set_overrides(args.extra_kws)
    adjust_text_kws = overrides.pop("adjust_text_kws", None)
    deep_merge(text_kws, overrides.pop("text_kws", {}))
    deep_merge(hline_kws, overrides.pop("hline_kws", {}))
    deep_merge(xticklabel_kws, overrides.pop("xticklabel_kws", {}))
    extra_kwargs = overrides  # any remaining top-level overrides

    # --- Create figure ---
    fig, ax = create_figure(args, default_figsize=(12, 4))

    ax = manhattanplot(
        data=df,
        chrom=args.chrom,
        pos=args.pos,
        pv=args.pv,
        snp=args.snp,
        logp=args.logp,
        marker=args.marker,
        color=args.color,
        alpha=args.alpha,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        xtick_label_set=xtick_label_set,
        CHR=args.chromosome,
        xticklabel_kws=xticklabel_kws if xticklabel_kws else None,
        suggestiveline=suggestiveline,
        genomewideline=genomewideline,
        sign_line_cols=args.sign_line_colors,
        hline_kws=hline_kws,
        sign_marker_p=args.sign_marker_p,
        sign_marker_color=args.sign_marker_color,
        is_annotate_topsnp=args.annotate_topsnp,
        text_kws=text_kws if args.annotate_topsnp else None,
        ld_block_size=args.ld_block_size,
        annotate_fmt=annotate_fmt,
        annotate_layout=args.annotate_layout,
        adjust_text_kws=adjust_text_kws,
        style=args.style,
        ax=ax,
        **extra_kwargs,
    )

    # --- Save output ---
    output_path = resolve_output_path(args, "manhattan.png")
    fig.tight_layout()
    save_figure(fig, output_path, dpi=args.dpi)
