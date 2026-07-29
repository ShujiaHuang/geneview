"""Mitochondrial DNA (mtDNA) visualization module for geneview.

A purpose-built toolkit for the figures that matter most in human mtDNA
analysis, driven by the outputs of pipelines such as MitoQuest
(https://github.com/ShujiaHuang/mitoquest): multi/single-sample VCFs with
per-sample heteroplasmy fractions, ``copynum`` TSV files, and BAM/CRAM
alignments.

Reference backbone
------------------
get_mt_genes, gene_at, genes_in_range, is_mt_contig, MT_LENGTH
    The rCRS (NC_012920.1) 16,569 bp gene map (13 protein-coding, 22 tRNA,
    2 rRNA genes + the D-loop control region) and lookup helpers.

Readers
-------
read_mito_vcf
    Single/multi-sample VCF -> tidy long table (sample x site x allele) with
    heteroplasmy fraction (VAF), depth, genotype and HET/HOM status.
read_mito_copynumber
    One or more ``mitoquest copynum`` TSVs -> per-sample copy number + 95% CI.
read_mito_coverage
    BAM/CRAM files -> binned depth across the mitochondrial contig.

Plots
-----
mito_genome_map
    Circular rCRS map: genes coloured by type + variant lollipops (length =
    heteroplasmy).
heteroplasmy_scatter
    Linear position-vs-VAF landscape with a gene strip.
heteroplasmy_heatmap
    Samples x variant-sites VAF heatmap for cohorts.
mito_coverage_plot
    Sequencing depth across the mitochondrial genome.
mito_copynumber_plot
    Per-sample mtDNA copy number with confidence intervals.

Author: Shujia Huang
"""
from ._reference import (
    MT_LENGTH,
    MT_CONTIG_ALIASES,
    MT_FEATURE_COLORS,
    MT_FEATURE_LABELS,
    MT_HYPERVARIABLE_REGIONS,
    get_mt_genes,
    gene_at,
    genes_in_range,
    is_mt_contig,
)
from ._io import (
    read_mito_vcf,
    read_mito_copynumber,
    read_mito_coverage,
    MITO_VCF_COLUMNS,
)
from ._genome_map import mito_genome_map
from ._heteroplasmy import heteroplasmy_scatter, heteroplasmy_heatmap
from ._coverage import mito_coverage_plot, mito_copynumber_plot

__all__ = [
    # Reference
    "MT_LENGTH",
    "MT_CONTIG_ALIASES",
    "MT_FEATURE_COLORS",
    "MT_FEATURE_LABELS",
    "MT_HYPERVARIABLE_REGIONS",
    "get_mt_genes",
    "gene_at",
    "genes_in_range",
    "is_mt_contig",
    # Readers
    "read_mito_vcf",
    "read_mito_copynumber",
    "read_mito_coverage",
    "MITO_VCF_COLUMNS",
    # Plots
    "mito_genome_map",
    "heteroplasmy_scatter",
    "heteroplasmy_heatmap",
    "mito_coverage_plot",
    "mito_copynumber_plot",
]
