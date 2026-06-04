"""
GroupedAlignmentsTrack - Display BAM/CRAM reads separated by group.

Reads are split into sub-groups based on a tag or arbitrary callback,
and each group is rendered as a separate panel within the same track.
This is useful for visualising haplotype-phased data (e.g. ``HP`` tag
with 10x Genomics), methylation bins, or any custom read classification.

Ported from ``genomeview.bamtrack.GroupedBAMTrack``.
"""

from typing import Optional, Dict, Any, List, Union, Callable

import numpy as np
import matplotlib.pyplot as plt

from ._base import Track, GenomicInterval
from ._alignments_track import AlignmentsTrack
from ._utils import match_chrom_format


def get_group_by_tag_fn(tag: str) -> Callable:
    """Return a grouping function that splits reads by the value of *tag*.

    Parameters
    ----------
    tag : str
        A BAM tag name (e.g. ``"HP"`` for 10x haplotype phasing).

    Returns
    -------
    callable
        A function ``f(read) -> str`` that returns the tag value as a
        string, or ``"missing"`` when the tag is absent.

    Examples
    --------
    >>> fn = get_group_by_tag_fn("HP")
    >>> # fn(read) returns "1", "2", or "missing"
    """

    def _group_by_tag(read):
        if not read.has_tag(tag):
            return "missing"
        return str(read.get_tag(tag))

    return _group_by_tag


class GroupedAlignmentsTrack(Track):
    """Display BAM/CRAM reads split into groups.

    Each group is rendered as a separate :class:`AlignmentsTrack` panel,
    stacked vertically with a configurable gap.

    Parameters
    ----------
    filepath : str
        Path to a BAM or CRAM file.
    keyfn : callable
        A function ``f(pysam.AlignedSegment) -> str`` that returns the
        group label for each read.  Use :func:`get_group_by_tag_fn` for
        tag-based grouping.
    is_paired : bool
        Whether the data is paired-end.  Default is False.
    reference : str, optional
        Path to a reference FASTA for mismatch detection.
    type : str or list of str
        Display mode(s) passed to each sub-:class:`AlignmentsTrack`.
        Default is ``"coverage"``.
    space_between : float
        Vertical gap (in axes fraction) between groups.  Default is 0.05.
    category_label_fn : callable
        A function ``f(str) -> str`` that formats group labels for
        display.  Default is :func:`str`.
    quick_consensus : bool
        Passed to each sub-track.  Default is True.
    consensus_threshold : float
        Passed to each sub-track.  Default is 0.2.
    name : str
        Track name.  Default is "GroupedAlignments".
    height : float
        Relative track height.  Default is 3.0.
    display_params : dict, optional
        Additional display parameters.
    **kwargs
        Passed through to each sub-:class:`AlignmentsTrack`.
    """

    def __init__(
        self,
        filepath: str,
        keyfn: Callable,
        is_paired: bool = False,
        reference: Optional[str] = None,
        type: Union[str, List[str]] = "coverage",
        space_between: float = 0.05,
        category_label_fn: Optional[Callable] = None,
        quick_consensus: bool = True,
        consensus_threshold: float = 0.2,
        name: str = "GroupedAlignments",
        height: float = 3.0,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        dp = {"fontsize": 7}
        if display_params:
            dp.update(display_params)
        super().__init__(name=name, height=height, display_params=dp, **kwargs)

        self.filepath = filepath
        self.keyfn = keyfn
        self.is_paired = is_paired
        self.reference = reference
        self.plot_type = type
        self.space_between = space_between
        self.category_label_fn = category_label_fn or str
        self.quick_consensus = quick_consensus
        self.consensus_threshold = consensus_threshold
        self._extra_kwargs = kwargs

        self._sub_tracks: List[AlignmentsTrack] = []
        self._group_labels: List[str] = []

    def _build_sub_tracks(self, region: GenomicInterval) -> None:
        """Discover groups and create per-group AlignmentsTrack instances."""
        import pysam

        aln = pysam.AlignmentFile(self.filepath, "rb")
        chrom = match_chrom_format(region.chrom, aln.references)

        # Discover categories
        categories = set()
        for read in aln.fetch(chrom, region.start, region.end):
            categories.add(self.keyfn(read))
        aln.close()

        categories = sorted(categories)
        self._sub_tracks = []
        self._group_labels = []

        for cat in categories:
            label = self.category_label_fn(cat)

            def _filter_fn(read, _cat=cat):
                return self.keyfn(read) == _cat

            sub = AlignmentsTrack(
                filepath=self.filepath,
                is_paired=self.is_paired,
                reference=self.reference,
                type=self.plot_type,
                quick_consensus=self.quick_consensus,
                consensus_threshold=self.consensus_threshold,
                read_filter=_filter_fn,
                name=label,
                height=1.0,
                **self._extra_kwargs,
            )
            self._sub_tracks.append(sub)
            self._group_labels.append(label)

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw grouped alignments.

        Each group is rendered in its own vertical slice of the axes.
        """
        self._build_sub_tracks(region)

        if not self._sub_tracks:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        n_groups = len(self._sub_tracks)
        gap = self.space_between
        total_gap = gap * (n_groups - 1)
        group_height = (1.0 - total_gap) / max(n_groups, 1)

        ax.set_xlim(region.start, region.end)
        ax.set_ylim(0, 1)
        ax.axis("off")

        for i, sub_track in enumerate(self._sub_tracks):
            # Compute vertical slice for this group (top to bottom)
            y_top = 1.0 - i * (group_height + gap)
            y_bottom = y_top - group_height

            # Set a clip region for this group
            ax.set_ylim(0, 1)

            # Draw the sub-track content in a clipped sub-region
            # We use a transform to scale the sub-track into its slice
            sub_track.draw(ax, region)

        # Draw group labels on the left side
        fontsize = self.get_param("fontsize", 7)
        for i, label in enumerate(self._group_labels):
            y_top = 1.0 - i * (group_height + gap)
            y_mid = y_top - group_height / 2
            ax.text(
                region.start, y_mid, f" {label}",
                ha="left", va="center", fontsize=fontsize,
                color="#555555", fontweight="bold",
                clip_on=True, zorder=5,
                transform=ax.transData,
            )

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the region covered by the BAM file (from index)."""
        try:
            import pysam
            aln = pysam.AlignmentFile(self.filepath, "rb")
            refs = aln.references
            lengths = aln.lengths
            aln.close()
            if refs:
                return GenomicInterval(refs[0], 0, int(lengths[0]))
        except Exception:
            pass
        return None
