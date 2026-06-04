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

Core Functions
--------------
plot_tracks
    Main function for rendering stacked genome tracks.
GenomicInterval
    Data class representing a genomic region.

File I/O
--------
read_bed, read_gff, read_bedgraph, read_bigwig, read_bam_coverage, read_cram_coverage
    Readers for common genomic file formats.

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
)
from ._genome_axis import GenomeAxisTrack
from ._annotation import AnnotationTrack
from ._gene_region import GeneRegionTrack
from ._data_track import DataTrack
from ._highlight import HighlightTrack
from ._overlay import OverlayTrack
from ._ideogram import IdeogramTrack
from ._track_plot import plot_tracks
from ._io import (
    read_bed,
    read_gff,
    read_bedgraph,
    read_bigwig,
    read_bam_coverage,
    read_cram_coverage,
    read_auto,
)

__all__ = [
    # Track types
    "GenomeAxisTrack",
    "AnnotationTrack",
    "GeneRegionTrack",
    "DataTrack",
    "IdeogramTrack",
    "HighlightTrack",
    "OverlayTrack",
    # Core
    "plot_tracks",
    "GenomicInterval",
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
    "read_bam_coverage",
    "read_cram_coverage",
    "read_auto",
]
