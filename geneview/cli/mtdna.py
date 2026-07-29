"""CLI subcommands for mitochondrial (mtDNA) visualization.

Registers a family of ``mito-*`` subcommands driven by MitoQuest-style inputs
(multi/single-sample VCF, ``copynum`` TSV, BAM/CRAM):

    $ geneview mito-map          -i cohort.mt.vcf.gz -o mito_map.png
    $ geneview mito-heteroplasmy -i cohort.mt.vcf.gz -o het_scatter.png
    $ geneview mito-heatmap      -i cohort.mt.vcf.gz -o het_heatmap.png
    $ geneview mito-coverage     -i s1.bam s2.bam    -o mito_cov.png
    $ geneview mito-copynumber   -i *.cn.tsv         -o mito_cn.png

Author: Shujia Huang
"""
import os

import matplotlib.pyplot as plt

from .utils import (
    add_common_figure_args, add_style_arg, create_figure,
    save_figure, resolve_output_path,
)


def register(subparsers):
    """Register all ``mito-*`` subcommands on *subparsers*."""
    _register_map(subparsers)
    _register_heteroplasmy(subparsers)
    _register_heatmap(subparsers)
    _register_coverage(subparsers)
    _register_copynumber(subparsers)


# ---------------------------------------------------------------------------
# mito-map
# ---------------------------------------------------------------------------
def _register_map(subparsers):
    p = subparsers.add_parser(
        "mito-map",
        help="Circular rCRS map of the mitochondrial genome with variant lollipops.",
        description="Draw a circular map of the 16,569 bp human mtDNA (rCRS): "
                    "genes coloured by feature type, with variants from a VCF "
                    "overlaid as lollipops whose length encodes heteroplasmy (VAF).",
        epilog="Example: geneview mito-map -i cohort.mt.vcf.gz -o mito_map.png",
    )
    p.add_argument("-i", "--vcf", dest="vcf", default=None,
                   help="Input VCF/BCF (e.g. from 'mitoquest caller'). Omit to "
                        "draw the bare gene ring.")
    p.add_argument("--region", default=None,
                   help="Restrict to a region (chrom or chrom:start-end).")
    p.add_argument("--samples", nargs="+", default=None,
                   help="Restrict to these sample names.")
    p.add_argument("--min-vaf", type=float, default=0.0,
                   help="Drop variants with VAF below this value. (default: 0.0)")
    p.add_argument("--gene-label", choices=["large", "all", "none"], default="large",
                   help="Which genes to label. (default: large)")
    p.add_argument("--color-by", choices=["feature", "status", "var_type"],
                   default="feature",
                   help="Colour lollipops by feature/status/var_type. (default: feature)")
    p.add_argument("--title", default=None, help="Text drawn in the ring centre.")
    add_style_arg(p)
    add_common_figure_args(p)
    p.set_defaults(func=_run_map)


def _run_map(args):
    from geneview import mtdna

    variants = None
    if args.vcf:
        _require_file(args.vcf)
        variants = mtdna.read_mito_vcf(args.vcf, samples=args.samples,
                                       region=args.region, min_vaf=args.min_vaf)
    figsize = args.figsize if args.figsize else (8, 8)
    fig, ax = plt.subplots(figsize=figsize, facecolor=args.facecolor,
                           subplot_kw={"projection": "polar"})
    mtdna.mito_genome_map(variants, gene_label=args.gene_label,
                          color_by=args.color_by, title=args.title,
                          style=args.style, ax=ax)
    save_figure(fig, resolve_output_path(args, "mito_map.png"), dpi=args.dpi)


# ---------------------------------------------------------------------------
# mito-heteroplasmy (scatter)
# ---------------------------------------------------------------------------
def _register_heteroplasmy(subparsers):
    p = subparsers.add_parser(
        "mito-heteroplasmy",
        help="Linear position-vs-VAF heteroplasmy landscape from a VCF.",
        description="Scatter the heteroplasmy fraction (VAF) of every variant "
                    "against its mtDNA position, with a gene strip and a "
                    "heteroplasmy threshold line.",
        epilog="Example: geneview mito-heteroplasmy -i cohort.mt.vcf.gz --hue sample -o het.png",
    )
    _add_vcf_args(p)
    p.add_argument("--hue", choices=["feature", "status", "var_type", "sample"],
                   default="feature", help="Point colour encoding. (default: feature)")
    p.add_argument("--het-threshold", type=float, default=0.03,
                   help="VAF threshold line. Use a negative value to hide. (default: 0.03)")
    p.add_argument("--no-gene-band", action="store_true", default=False,
                   help="Hide the gene strip below the axis.")
    add_style_arg(p)
    add_common_figure_args(p)
    p.set_defaults(func=_run_heteroplasmy)


def _run_heteroplasmy(args):
    from geneview import mtdna

    _require_file(args.vcf)
    df = mtdna.read_mito_vcf(args.vcf, samples=args.samples,
                             region=args.region, min_vaf=args.min_vaf)
    fig, ax = create_figure(args, default_figsize=(12, 3.8))
    threshold = None if args.het_threshold is not None and args.het_threshold < 0 \
        else args.het_threshold
    mtdna.heteroplasmy_scatter(df, hue=args.hue, het_threshold=threshold,
                               show_gene_band=not args.no_gene_band,
                               style=args.style, ax=ax)
    save_figure(fig, resolve_output_path(args, "mito_heteroplasmy.png"), dpi=args.dpi)


# ---------------------------------------------------------------------------
# mito-heatmap
# ---------------------------------------------------------------------------
def _register_heatmap(subparsers):
    p = subparsers.add_parser(
        "mito-heatmap",
        help="Samples x variant-sites heteroplasmy (VAF) heatmap from a VCF.",
        description="Draw a cohort heatmap of heteroplasmy fractions "
                    "(rows = samples, columns = variant sites).",
        epilog="Example: geneview mito-heatmap -i cohort.mt.vcf.gz -o heatmap.png",
    )
    _add_vcf_args(p)
    p.add_argument("--cmap", default="rocket_r",
                   help="Colormap for the VAF scale. (default: rocket_r)")
    p.add_argument("--site-label", choices=["pos", "variant"], default="pos",
                   help="Column label style. (default: pos)")
    add_style_arg(p)
    add_common_figure_args(p)
    p.set_defaults(func=_run_heatmap)


def _run_heatmap(args):
    from geneview import mtdna

    _require_file(args.vcf)
    df = mtdna.read_mito_vcf(args.vcf, samples=args.samples,
                             region=args.region, min_vaf=args.min_vaf)
    fig, ax = create_figure(args, default_figsize=(12, 5))
    mtdna.heteroplasmy_heatmap(df, cmap=args.cmap, site_label=args.site_label,
                               style=args.style, ax=ax)
    save_figure(fig, resolve_output_path(args, "mito_heatmap.png"), dpi=args.dpi)


# ---------------------------------------------------------------------------
# mito-coverage
# ---------------------------------------------------------------------------
def _register_coverage(subparsers):
    p = subparsers.add_parser(
        "mito-coverage",
        help="Sequencing depth across the mitochondrial genome from BAM/CRAM.",
        description="Compute and plot per-bin mtDNA coverage for one or more "
                    "indexed BAM/CRAM files.",
        epilog="Example: geneview mito-coverage -i s1.bam s2.bam -o mito_cov.png",
    )
    p.add_argument("-i", "--input", nargs="+", required=True, dest="inputs",
                   help="Input BAM/CRAM file(s), one per sample (indexed).")
    p.add_argument("--reference", default=None,
                   help="Reference FASTA (required to decode CRAM).")
    p.add_argument("--sample-names", nargs="+", default=None,
                   help="Sample labels aligned with -i (default: from file names).")
    p.add_argument("--bins", type=int, default=1000,
                   help="Number of bins across the contig. (default: 1000)")
    p.add_argument("--contig", default=None,
                   help="Force the mitochondrial contig name (default: auto-detect).")
    p.add_argument("--log", action="store_true", default=False,
                   help="Use a logarithmic depth axis.")
    add_style_arg(p)
    add_common_figure_args(p)
    p.set_defaults(func=_run_coverage)


def _run_coverage(args):
    from geneview import mtdna

    for f in args.inputs:
        _require_file(f)
    cov = mtdna.read_mito_coverage(args.inputs, sample_names=args.sample_names,
                                   bins=args.bins, reference=args.reference,
                                   contig=args.contig)
    fig, ax = create_figure(args, default_figsize=(12, 3.8))
    mtdna.mito_coverage_plot(cov, log=args.log, style=args.style, ax=ax)
    save_figure(fig, resolve_output_path(args, "mito_coverage.png"), dpi=args.dpi)


# ---------------------------------------------------------------------------
# mito-copynumber
# ---------------------------------------------------------------------------
def _register_copynumber(subparsers):
    p = subparsers.add_parser(
        "mito-copynumber",
        help="Per-sample mtDNA copy number with 95%% CI from copynum TSVs.",
        description="Read one or more 'mitoquest copynum' TSV files and plot "
                    "the per-sample mtDNA copy number with its confidence interval.",
        epilog="Example: geneview mito-copynumber -i *.cn.tsv -o mito_cn.png",
    )
    p.add_argument("-i", "--input", nargs="+", required=True, dest="inputs",
                   help="Input 'mitoquest copynum' TSV file(s), one per sample.")
    p.add_argument("--sample-names", nargs="+", default=None,
                   help="Sample labels aligned with -i (default: from file names).")
    p.add_argument("--orient", choices=["v", "h"], default="v",
                   help="Bar orientation: v (vertical) or h (horizontal). (default: v)")
    p.add_argument("--no-ci", action="store_true", default=False,
                   help="Hide the 95%% confidence-interval error bars.")
    p.add_argument("--baseline", type=float, default=None,
                   help="Draw a reference line at this copy number.")
    add_style_arg(p)
    add_common_figure_args(p)
    p.set_defaults(func=_run_copynumber)


def _run_copynumber(args):
    from geneview import mtdna

    for f in args.inputs:
        _require_file(f)
    df = mtdna.read_mito_copynumber(args.inputs, sample_names=args.sample_names)
    if len(df) == 0:
        raise ValueError("No mitochondrial rows found in the given copynum TSV(s).")
    fig, ax = create_figure(args, default_figsize=(7, 4.5))
    mtdna.mito_copynumber_plot(df, orient=args.orient, show_ci=not args.no_ci,
                               baseline=args.baseline, style=args.style, ax=ax)
    save_figure(fig, resolve_output_path(args, "mito_copynumber.png"), dpi=args.dpi)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _add_vcf_args(p):
    """Add the VCF-input arguments shared by the VCF-driven subcommands."""
    p.add_argument("-i", "--vcf", dest="vcf", required=True,
                   help="Input VCF/BCF file (e.g. from 'mitoquest caller').")
    p.add_argument("--region", default=None,
                   help="Restrict to a region (chrom or chrom:start-end).")
    p.add_argument("--samples", nargs="+", default=None,
                   help="Restrict to these sample names.")
    p.add_argument("--min-vaf", type=float, default=0.0,
                   help="Drop variants with VAF below this value. (default: 0.0)")


def _require_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input file not found: {path}")
