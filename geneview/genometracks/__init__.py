"""
Genome track visualization module for geneview.

This module provides a genome browser-style track plotting system, inspired
by Gviz (R/Bioconductor) and adapted for Python/matplotlib.

Track Types
-----------
GenomeAxisTrack
    Genomic coordinate axis/ruler.
AnnotationTrack
    Generic genomic range annotations (BED/GFF features).
GeneRegionTrack
    Gene models with exons, UTRs, and introns.
DataTrack
    Numeric data visualization (BigWig, bedGraph).
IdeogramTrack
    Chromosome ideogram with cytoband coloring.
HighlightTrack
    Cross-track region highlighting.
OverlayTrack
    Overlay multiple tracks on the same panel.
AlignmentsTrack
    BAM/CRAM read alignments with coverage, pileup, and sashimi modes.
BAMCoverageTrack
    Standalone per-base coverage line/fill track from BAM/CRAM.
GroupedAlignmentsTrack
    BAM/CRAM reads split into groups (e.g. by haplotype tag).
VCFTrack
    VCF/BCF variant display with custom coloring.
LolliplotTrack
    Lollipop-style variant plot: stems + shapes on gene features.
DandelionTrack
    Dandelion-style clustered variant plot: fanning stems from groups.

Core Functions
--------------
plot_tracks
    Main function for rendering stacked genome tracks.
find_tracks
    Retrieve tracks by name or type from a track list.
plot_tracks_grid
    Side-by-side comparison of multiple genomic views.
plot_tracks_multi
    Stacked sections from different genomic regions.
GenomicInterval
    Data class representing a genomic region.
visualize_files
    Quick convenience function: file paths -> track list.
lolliplot
    Convenience function for :class:`LolliplotTrack`.
dandelion_plot
    Convenience function for :class:`DandelionTrack`.

File I/O
--------
read_bed, read_gff, read_bedgraph, read_bigwig, read_bigbed, read_bam_coverage, read_cram_coverage
    Readers for common genomic file formats.

Utilities
---------
match_chrom_format
    Normalise chromosome names (e.g. ``chr1`` <-> ``1``).
get_ticks
    Standalone utility for computing nicely-spaced axis ticks.
is_paired_end, is_long_frag_dataset
    BAM file introspection helpers.
MismatchCounts
    Quick-consensus mismatch tallying for long-read filtering.

Examples
--------
>>> from geneview.genometracks import (
...     plot_tracks, GenomeAxisTrack, AnnotationTrack,
...     GeneRegionTrack, DataTrack, GenomicInterval,
... )
>>> import pandas as pd
>>> region = GenomicInterval("chr7", 2000000, 2200000)
>>> ax_track = GenomeAxisTrack()
>>> _ = plot_tracks([ax_track], region=region)
"""

from ._base import (
    GenomicInterval,
    Track,
    RangeTrack,
    StackedTrack,
    NumericTrack,
    available_display_params,
)
from ._genome_axis import GenomeAxisTrack, get_ticks
from ._annotation import AnnotationTrack, DetailsAnnotationTrack
from ._gene_region import GeneRegionTrack
from ._data_track import DataTrack
from ._highlight import HighlightTrack
from ._overlay import OverlayTrack
from ._ideogram import IdeogramTrack
from ._sequence_track import SequenceTrack
from ._alignments_track import AlignmentsTrack, BAMCoverageTrack
from ._grouped_alignments import GroupedAlignmentsTrack, get_group_by_tag_fn
from ._vcf_track import VCFTrack
from ._lollipop_track import LolliplotTrack, lolliplot
from ._dandelion_track import DandelionTrack, dandelion_plot
from ._mismatch_counts import MismatchCounts
from ._biomart import BiomartGeneRegionTrack
from ._ucsc import UcscTrack
from ._track_plot import plot_tracks, find_tracks
from ._io import (
    read_bed,
    read_gff,
    read_bedgraph,
    read_bigwig,
    read_bigbed,
    read_bam_coverage,
    read_cram_coverage,
    read_wig,
    read_fasta,
    read_2bit,
    read_auto,
    GenomeSource,
    FastaGenomeSource,
)
from ._export import export_tracks, save_figure
from ._schemes import apply_scheme
from ._utils import match_chrom_format, is_paired_end, is_long_frag_dataset, reverse_comp
from ._convenience import visualize_files
from ._multi_view import plot_tracks_grid, plot_tracks_multi

__all__ = [
    # Track types
    "GenomeAxisTrack",
    "AnnotationTrack",
    "DetailsAnnotationTrack",
    "GeneRegionTrack",
    "DataTrack",
    "IdeogramTrack",
    "HighlightTrack",
    "OverlayTrack",
    "SequenceTrack",
    "AlignmentsTrack",
    "BAMCoverageTrack",
    "GroupedAlignmentsTrack",
    "VCFTrack",
    "LolliplotTrack",
    "DandelionTrack",
    "BiomartGeneRegionTrack",
    "UcscTrack",
    # Core
    "plot_tracks",
    "find_tracks",
    "plot_tracks_grid",
    "plot_tracks_multi",
    "GenomicInterval",
    "available_display_params",
    # Convenience
    "visualize_files",
    "lolliplot",
    "dandelion_plot",
    "get_group_by_tag_fn",
    # Base classes (for advanced usage)
    "Track",
    "RangeTrack",
    "StackedTrack",
    "NumericTrack",
    # I/O utilities
    "read_bed",
    "read_gff",
    "read_bedgraph",
    "read_bigwig",
    "read_bigbed",
    "read_bam_coverage",
    "read_cram_coverage",
    "read_wig",
    "read_fasta",
    "read_2bit",
    "read_auto",
    "GenomeSource",
    "FastaGenomeSource",
    # Export
    "export_tracks",
    "save_figure",
    # Schemes
    "apply_scheme",
    # Utilities
    "match_chrom_format",
    "is_paired_end",
    "is_long_frag_dataset",
    "MismatchCounts",
    "reverse_comp",
    "get_ticks",
]
