"""
Convenience helpers for quickly visualising genomic data files.

Provides :func:`visualize_files` which accepts a list of file paths
(BAM, BED, BigWig, etc.) and a genomic region, auto-detects file types,
constructs the appropriate tracks, and returns a ready-to-render track
list that can be passed directly to :func:`plot_tracks`.

Ported from ``genomeview.convenience.visualize_data``.
"""

import os
from typing import List, Optional, Union, Dict

from ._base import Track, GenomicInterval
from ._genome_axis import GenomeAxisTrack


def visualize_files(
    file_paths: Union[List[str], Dict[str, str]],
    region: GenomicInterval,
    reference: Optional[str] = None,
    axis_on_top: bool = False,
    show_mismatches: bool = True,
    quick_consensus: bool = True,
    consensus_threshold: float = 0.2,
) -> List[Track]:
    """Build a track list for the given files and genomic region.

    Auto-detects file types from extensions and creates the appropriate
    track objects.  Supported extensions:

    * ``.bam``, ``.cram`` -- :class:`AlignmentsTrack`
    * ``.bed``, ``.bed.gz`` -- :class:`AnnotationTrack`
    * ``.gff``, ``.gff3``, ``.gtf`` -- :class:`AnnotationTrack`
    * ``.bigwig``, ``.bw`` -- :class:`DataTrack`
    * ``.bedgraph``, ``.bdg`` -- :class:`DataTrack`

    Parameters
    ----------
    file_paths : list of str or dict
        File paths to visualise.  If a ``dict``, keys are used as track
        names and values as file paths.
    region : GenomicInterval
        The genomic region to display.
    reference : str, optional
        Path to a reference FASTA for BAM mismatch display.
    axis_on_top : bool
        If True, place the genome axis at the top of the track list.
        Default is False (axis at bottom).
    show_mismatches : bool
        Passed to :class:`AlignmentsTrack`.  Default is True.
    quick_consensus : bool
        Passed to :class:`AlignmentsTrack`.  Default is True.
    consensus_threshold : float
        Passed to :class:`AlignmentsTrack`.  Default is 0.2.

    Returns
    -------
    list of Track
        A list of track objects ready to be passed to :func:`plot_tracks`.

    Examples
    --------
    >>> from geneview.genometracks import visualize_files, plot_tracks, GenomicInterval
    >>> tracks = visualize_files(                          # doctest: +SKIP
    ...     ["sample.bam", "genes.bed.gz"],
    ...     region=GenomicInterval("chr1", 100000, 200000),
    ... )
    >>> axes = plot_tracks(tracks)                         # doctest: +SKIP
    """
    from ._alignments_track import AlignmentsTrack
    from ._annotation import AnnotationTrack
    from ._data_track import DataTrack
    from ._io import read_auto
    from ._utils import is_paired_end, is_long_frag_dataset

    if isinstance(file_paths, dict):
        names = list(file_paths.keys())
        paths = list(file_paths.values())
    else:
        names = [None] * len(file_paths)
        paths = list(file_paths)

    tracks: List[Track] = []

    def _add_axis():
        tracks.append(GenomeAxisTrack())

    if axis_on_top:
        _add_axis()

    for name, path in zip(names, paths):
        ext = _normalise_ext(path)
        track_name = name or os.path.basename(path).split(".")[0]

        if ext in (".bam", ".cram"):
            is_paired = False
            min_indel_size = 0
            try:
                is_paired = is_paired_end(path)
                if not is_paired:
                    try:
                        if is_long_frag_dataset(path):
                            min_indel_size = 5
                    except Exception:
                        pass
            except Exception:
                pass

            sub = AlignmentsTrack(
                filepath=path,
                is_paired=is_paired,
                reference=reference,
                show_mismatches=show_mismatches and reference is not None,
                quick_consensus=quick_consensus,
                consensus_threshold=consensus_threshold,
                name=track_name,
            )
            tracks.append(sub)

        elif ext in (".bed", ".bed.gz", ".bigbed", ".bb",
                      ".gff", ".gff3", ".gtf"):
            data = read_auto(path)
            sub = AnnotationTrack(data, name=track_name)
            tracks.append(sub)

        elif ext in (".bigwig", ".bw", ".bedgraph", ".bdg"):
            data = read_auto(path)
            sub = DataTrack(data, name=track_name)
            tracks.append(sub)

        else:
            raise ValueError(
                f"Unsupported file extension for '{path}'. "
                f"Supported: .bam, .cram, .bed, .bed.gz, .gff, .gff3, "
                f".gtf, .bigwig, .bw, .bedgraph, .bdg"
            )

    if not axis_on_top:
        _add_axis()

    return tracks


def _normalise_ext(path: str) -> str:
    """Return the normalised lower-case extension of *path*."""
    base = os.path.basename(path).lower()
    # Handle double extensions like .bed.gz
    if base.endswith(".gz"):
        base = base[:-3]
    _, ext = os.path.splitext(base)
    return ext
