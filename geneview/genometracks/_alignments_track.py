"""
AlignmentsTrack - A track for displaying read alignments from BAM/CRAM files.

Displays sequencing reads as coverage histograms, pileup diagrams, or sashimi
plots, with optional mismatch and indel highlighting.

Ported from Gviz's AlignmentsTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union, Callable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, Arc
from matplotlib.collections import PolyCollection

from ._base import Track, StackedTrack, GenomicInterval
from ._utils import match_chrom_format


# CIGAR operation codes (from pysam / SAM spec)
_CIGAR_M = 0   # alignment match (can be seq match or mismatch)
_CIGAR_I = 1   # insertion to reference
_CIGAR_D = 2   # deletion from reference
_CIGAR_N = 3   # skipped region (intron)
_CIGAR_S = 4   # soft clipping
_CIGAR_H = 5   # hard clipping
_CIGAR_P = 6   # padding
_CIGAR_EQ = 7  # sequence match
_CIGAR_X = 8   # sequence mismatch

# Nucleotide colors (matching genomeview: A=blue, C=orange, G=green, T=black)
_NUC_COLORS = {
    "A": "blue", "C": "orange", "G": "green", "T": "black",
    "N": "gray",
}


class AlignmentsTrack(StackedTrack):
    """Display BAM/CRAM read alignments.

    Requires the optional ``pysam`` package.

    Parameters
    ----------
    filepath : str
        Path to a BAM or CRAM file.
    is_paired : bool
        Whether the data is paired-end.  Default is False.
    show_mismatches : bool
        If True, color-code mismatched bases against the reference.
        Requires *reference*.  Default is True.
    show_indels : bool
        If True, draw insertions as vertical bars and deletions as
        bridging lines.  Default is True.
    reference : str, optional
        Path to a reference FASTA (indexed) for mismatch detection.
    type : str or list of str
        Display mode(s): ``"coverage"``, ``"pileup"``, or ``"sashimi"``.
        Can be a list for combined display (e.g. ``["coverage", "pileup"]``).
        Default is ``"coverage"``.
    coverage_height : float
        Relative height of the coverage panel.  Default is 0.4.
    sashimi_height : float
        Relative height of the sashimi arc panel.  Default is 0.3.
    sashimi_score : int
        Minimum number of reads supporting a junction to draw it.
        Default is 1.
    sashimi_filter : pd.DataFrame, optional
        DataFrame with columns ``chrom``, ``start``, ``end`` defining
        junction regions to include in the sashimi plot.
    sashimi_filter_tolerance : int
        Tolerance in bp for matching sashimi filter regions.  Default is 0.
    reverse_stacking : bool
        If True, reverse the stacking order.  Default is False.
    col_mates : str
        Color for mate-pair connectors.  Default is ``"lightblue"``.
    col_gap : str
        Color for gap (intron) lines in reads.  Default is ``"lightgray"``.
    col_deletion : str
        Color for deletion indicators.  Default is ``"red"``.
    col_insertion : str
        Color for insertion indicators.  Default is ``"blue"``.
    fill_coverage : str
        Fill color for coverage histogram.  Default is ``"#5B8DB8"``.
    fill_reads : str
        Fill color for pileup reads.  Default is ``"#BDBDBD"``.
    alpha_reads : float
        Alpha for pileup reads.  Default is 0.8.
    alpha_mismatch : float
        Alpha for mismatch highlights.  Default is 0.6.
    transformation : callable, optional
        Function applied to coverage values before plotting.
    quick_consensus : bool
        If True (default), use per-position allele-frequency tallies to
        filter mismatches: a mismatch is only drawn when the alt allele is
        supported by at least *consensus_threshold* fraction of reads at
        that position.  Essential for noisy long-read data (PacBio / ONT).
    consensus_threshold : float
        Minimum alt-allele fraction required to draw a mismatch when
        *quick_consensus* is enabled.  Default is 0.2.
    read_filter : callable, optional
        A function ``f(pysam.AlignedSegment) -> bool``; only reads for
        which it returns ``True`` are included in pileup and coverage
        calculations.  Useful for filtering by MAPQ, tag, etc.
    show_clipping : bool
        If True, draw soft/hard clipped bases as coloured blocks at read
        edges.  Default is False.
    col_clipping : str
        Color for clipping indicators.  Default is ``"cyan"``.
    min_indel_size : int
        Minimum indel length (in bp) to draw.  Indels smaller than this
        are silently skipped.  Useful for long-read data where small
        indels are common noise.  Default is 0 (show all).
    show_insertion_labels : bool
        If True, draw insertion lengths as text labels when there is
        sufficient space.  Default is False.
    color_by_strand : bool
        If True, colour forward-strand reads differently from reverse-
        strand reads (``fill_reads_fwd`` / ``fill_reads_rev``).
        Default is False.
    fill_reads_fwd : str
        Fill color for forward-strand reads when ``color_by_strand`` is
        True.  Default is ``"#E89E9D"``.
    fill_reads_rev : str
        Fill color for reverse-strand reads when ``color_by_strand`` is
        True.  Default is ``"#8C8FCE"``.
    include_secondary : bool
        If True (default), include secondary and supplementary alignments
        in the pileup view.  Set to False to show only primary alignments.
    overlap_color : str, optional
        Color used to highlight overlapping portions of paired-end read
        mates.  When set, regions where both mates overlap are drawn as a
        thin overlay bar on top of the reads.  Default is None (disabled).
    draw_read_labels : bool
        If True, draw read query names as text labels to the right of each
        read in pileup mode.  Default is False.
    min_insertion_label_size : int
        Minimum insertion length (in bp) for which an insertion-size label
        is drawn automatically, even when ``show_insertion_labels`` is
        False.  Insertions shorter than this still show the vertical bar
        but omit the label.  Default is 5.
    color_fn : callable, optional
        A function ``f(pysam.AlignedSegment) -> str`` that returns a
        matplotlib color string for each read.  When set, this overrides
        ``fill_reads``, ``color_by_strand``, and related color parameters.
        Useful for coloring reads by insert size, mapping quality, or
        any other per-read property (equivalent to genomeview's
        ``track.color_fn``).
    name : str
        Track name.  Default is "Alignments".
    height : float
        Relative track height.  Default is 2.0.
    display_params : dict, optional
        Additional display parameters.
    """

    def __init__(
        self,
        filepath: str,
        is_paired: bool = False,
        show_mismatches: bool = True,
        show_indels: bool = True,
        reference: Optional[str] = None,
        type: Union[str, List[str]] = "coverage",
        coverage_height: float = 0.4,
        sashimi_height: float = 0.3,
        sashimi_score: int = 1,
        sashimi_filter: Optional[pd.DataFrame] = None,
        sashimi_filter_tolerance: int = 0,
        reverse_stacking: bool = False,
        col_mates: str = "gray",
        col_gap: str = "lightgray",
        col_deletion: str = "black",
        col_insertion: str = "purple",
        fill_coverage: str = "#5B8DB8",
        fill_reads: str = "#BDBDBD",
        alpha_reads: float = 0.8,
        alpha_mismatch: float = 0.6,
        transformation: Optional[Callable] = None,
        quick_consensus: bool = True,
        consensus_threshold: float = 0.2,
        read_filter: Optional[Callable] = None,
        show_clipping: bool = False,
        col_clipping: str = "cyan",
        min_indel_size: int = 0,
        show_insertion_labels: bool = False,
        color_by_strand: bool = True,
        fill_reads_fwd: str = "#E89E9D",
        fill_reads_rev: str = "#8C8FCE",
        include_secondary: bool = True,
        overlap_color: Optional[str] = None,
        draw_read_labels: bool = False,
        min_insertion_label_size: int = 5,
        color_fn: Optional[Callable] = None,
        name: str = "Alignments",
        height: float = 2.0,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        dp = {"fontsize": 7}
        if display_params:
            dp.update(display_params)
        super().__init__(data=None, name=name, height=height, display_params=dp,
                         **kwargs)

        self.filepath = filepath
        self.is_paired = is_paired
        self.show_mismatches = show_mismatches
        self.show_indels = show_indels
        self.reference = reference
        self.coverage_height = coverage_height
        self.sashimi_height = sashimi_height
        self.sashimi_score = sashimi_score
        self.sashimi_filter = sashimi_filter
        self.sashimi_filter_tolerance = sashimi_filter_tolerance
        self.reverse_stacking = reverse_stacking
        self.col_mates = col_mates
        self.col_gap = col_gap
        self.col_deletion = col_deletion
        self.col_insertion = col_insertion
        self.fill_coverage = fill_coverage
        self.fill_reads = fill_reads
        self.alpha_reads = alpha_reads
        self.alpha_mismatch = alpha_mismatch
        self.transformation = transformation
        self.quick_consensus = quick_consensus
        self.consensus_threshold = consensus_threshold
        self.read_filter = read_filter
        self.show_clipping = show_clipping
        self.col_clipping = col_clipping
        self.min_indel_size = min_indel_size
        self.show_insertion_labels = show_insertion_labels
        self.color_by_strand = color_by_strand
        self.fill_reads_fwd = fill_reads_fwd
        self.fill_reads_rev = fill_reads_rev
        self.include_secondary = include_secondary
        self.overlap_color = overlap_color
        self.draw_read_labels = draw_read_labels
        self.min_insertion_label_size = min_insertion_label_size
        self.color_fn = color_fn

        if isinstance(type, str):
            self.plot_types = [type]
        else:
            self.plot_types = list(type)

    # ------------------------------------------------------------------
    # pysam helpers
    # ------------------------------------------------------------------

    def _import_pysam(self):
        """Import pysam or raise a clear error."""
        try:
            import pysam
            return pysam
        except ImportError:
            raise ImportError(
                "The 'pysam' package is required for AlignmentsTrack. "
                "Install it with: pip install pysam"
            )

    def _fetch_reads(self, region: GenomicInterval):
        """Fetch reads overlapping the region. Returns pysam AlignmentFile and list.

        When ``self.read_filter`` is set, only reads passing the filter
        callback are included.
        """
        pysam = self._import_pysam()
        aln = pysam.AlignmentFile(self.filepath, "rb")
        chrom = match_chrom_format(region.chrom, aln.references)
        raw_reads = aln.fetch(chrom, region.start, region.end)
        if self.read_filter is not None:
            reads = [r for r in raw_reads if self.read_filter(r)]
        else:
            reads = list(raw_reads)
        # Filter secondary/supplementary if requested
        if not self.include_secondary:
            reads = [r for r in reads
                     if not r.is_secondary and not r.is_supplementary]
        return aln, reads

    def _build_mismatch_counts(self, region: GenomicInterval):
        """Build a MismatchCounts tally for the region (quick-consensus mode).

        Returns ``None`` when quick_consensus is disabled or no reference
        is available.
        """
        if not self.quick_consensus:
            return None
        pysam = self._import_pysam()
        from ._mismatch_counts import MismatchCounts
        aln = pysam.AlignmentFile(self.filepath, "rb")
        try:
            mc = MismatchCounts(region.chrom, region.start, region.end)
            mc.tally_reads(aln)
            return mc
        finally:
            aln.close()

    def _compute_coverage(self, region: GenomicInterval) -> np.ndarray:
        """Compute per-base coverage across the region."""
        pysam = self._import_pysam()
        aln = pysam.AlignmentFile(self.filepath, "rb")
        chrom = match_chrom_format(region.chrom, aln.references)
        cov = aln.count_coverage(
            chrom, region.start, region.end,
            quality_threshold=0,
        )
        # cov is (A, C, G, T) arrays; sum for total coverage
        total = np.array(cov[0]) + np.array(cov[1]) + np.array(cov[2]) + np.array(cov[3])
        aln.close()
        if self.transformation is not None:
            total = self.transformation(total.astype(float))
        return total

    def _get_reference_seq(self, region: GenomicInterval) -> Optional[str]:
        """Load reference sequence for the region (if reference is set)."""
        if self.reference is None:
            return None
        pysam = self._import_pysam()
        fa = pysam.FastaFile(self.reference)
        try:
            chrom = match_chrom_format(region.chrom, fa.references)
            return fa.fetch(chrom, region.start, region.end).upper()
        finally:
            fa.close()

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the alignments track.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        ax.set_xlim(region.start, region.end)
        ax.set_ylim(0, 1)
        ax.axis("off")

        if "coverage" in self.plot_types:
            self._draw_coverage(ax, region)

        if "pileup" in self.plot_types:
            self._draw_pileup(ax, region)

        if "sashimi" in self.plot_types:
            self._draw_sashimi(ax, region)

    def _draw_coverage(self, ax, region: GenomicInterval) -> None:
        """Draw a coverage depth histogram."""
        cov = self._compute_coverage(region)
        if len(cov) == 0:
            return

        positions = np.arange(region.start, region.start + len(cov))
        max_cov = np.max(cov) if len(cov) > 0 else 1
        if max_cov == 0:
            max_cov = 1

        # Normalize to fill upper portion of axes
        y_scale = self.coverage_height
        norm_cov = (cov / max_cov) * y_scale

        # Draw as filled area
        ax.fill_between(
            positions, 1.0 - norm_cov, 1.0,
            color=self.fill_coverage, alpha=0.7, step="mid", zorder=3,
        )
        ax.plot(
            positions, 1.0 - norm_cov,
            color="#3C5488", linewidth=0.5, zorder=4,
        )

        # Coverage axis label
        ax.text(
            region.start, 1.02, f"max={int(max_cov)}",
            ha="left", va="bottom", fontsize=6, color="#555555",
        )

    def _draw_pileup(self, ax, region: GenomicInterval) -> None:
        """Draw individual reads as boxes with gaps and mismatches."""
        try:
            aln, reads = self._fetch_reads(region)
        except Exception:
            return

        if not reads:
            return

        # Reference sequence for mismatch detection
        ref_seq = self._get_reference_seq(region) if self.show_mismatches else None

        # Quick-consensus mismatch tallies
        mismatch_counts = self._build_mismatch_counts(region) if self.show_mismatches else None

        # Compute stacking for reads
        read_starts = []
        read_ends = []
        for read in reads:
            if read.is_unmapped:
                continue
            rs = read.reference_start
            re = read.reference_end
            if rs is not None and re is not None:
                read_starts.append(rs)
                read_ends.append(re)

        if not read_starts:
            aln.close()
            return

        from ._stacking import compute_stacking
        starts_arr = np.array(read_starts)
        ends_arr = np.array(read_ends)

        if self.reverse_stacking:
            starts_arr = starts_arr[::-1]
            ends_arr = ends_arr[::-1]
            reads = [r for r in reads if not r.is_unmapped][::-1]
        else:
            reads = [r for r in reads if not r.is_unmapped]

        stacks = compute_stacking(starts_arr, ends_arr, mode="squish", min_distance=0)
        n_rows = int(stacks.max()) + 1 if len(stacks) > 0 else 1

        # Pileup occupies lower portion of axes
        pileup_bottom = 0.0
        pileup_top = 1.0 - self.coverage_height if "coverage" in self.plot_types else 1.0
        row_height = (pileup_top - pileup_bottom) / max(n_rows, 1) * 0.85

        span = region.end - region.start
        h = row_height * 0.8

        for i, read in enumerate(reads):
            if i >= len(stacks):
                break
            stack_row = stacks[i]
            y_center = pileup_bottom + (pileup_top - pileup_bottom) - (stack_row + 0.5) * row_height

            self._draw_single_read(
                ax, read, region, y_center, row_height, ref_seq, span,
                mismatch_counts=mismatch_counts,
            )

            # Draw read label to the right of the read
            if self.draw_read_labels and read.query_name:
                read_end_x = min(read.reference_end or region.end, region.end)
                ax.text(
                    read_end_x + span * 0.005, y_center,
                    read.query_name,
                    ha="left", va="center",
                    fontsize=3.5, color="#555555", alpha=0.8,
                    zorder=5,
                )

            # Draw overlap highlight for paired-end reads
            if (self.overlap_color and self.is_paired
                    and read.is_proper_pair
                    and read.next_reference_start is not None):
                mate_start = read.next_reference_start
                read_end = read.reference_end or 0
                if mate_start < read_end:
                    # Overlap region exists
                    ov_start = max(mate_start, region.start)
                    ov_end = min(read_end, region.end)
                    if ov_end > ov_start:
                        ov_rect = mpatches.Rectangle(
                            (ov_start, y_center - h * 0.4),
                            ov_end - ov_start, h * 0.8,
                            facecolor=self.overlap_color,
                            edgecolor="none",
                            alpha=0.35, zorder=5,
                        )
                        ax.add_patch(ov_rect)

            # Draw mate connector for paired-end reads
            if self.is_paired and read.is_proper_pair and read.next_reference_start is not None:
                mate_start = read.next_reference_start
                if region.start <= mate_start <= region.end:
                    ax.plot(
                        [read.reference_end, mate_start],
                        [y_center, y_center],
                        color=self.col_mates, linewidth=0.3,
                        alpha=0.5, zorder=1,
                    )

        aln.close()

    def _draw_read_body(self, ax, read, region, x_start, x_end, y_center, h,
                         read_color):
        """Draw read body as a block arrow (genomeview style).

        The arrow tip sits exactly at the read's alignment boundary
        (``x_end`` for forward, ``x_start`` for reverse) and the body
        is cut inward so the visual extent never exceeds the true
        alignment coordinates.

        The arrow depth is computed in *genomic* coordinates using the
        axes aspect ratio so it is visually ~half the read height.
        """
        x_start = max(x_start, region.start)
        x_end = min(x_end, region.end)
        if x_end <= x_start:
            return

        width = x_end - x_start
        half_h = h / 2

        # Convert arrow depth from visual (y) to genomic (x) coordinates
        ax_pos = ax.get_position()
        fig = ax.get_figure()
        fig_w, fig_h = fig.get_size_inches()
        axes_w_inches = ax_pos.width * fig_w
        axes_h_inches = ax_pos.height * fig_h
        span = region.end - region.start
        if axes_h_inches > 0 and axes_w_inches > 0:
            genomic_per_visual = (span / axes_w_inches) / (1.0 / axes_h_inches)
            arrow_w = h * 0.5 * genomic_per_visual
        else:
            arrow_w = span * 0.003

        # Cap arrow depth to at most 30% of read width
        arrow_w = min(arrow_w, width * 0.3)

        if arrow_w > 0 and width > span * 0.001:
            if not read.is_reverse:
                # Forward: arrow tip at x_end, body cut inward
                verts = [
                    (x_start, y_center - half_h),
                    (x_end - arrow_w, y_center - half_h),
                    (x_end, y_center),
                    (x_end - arrow_w, y_center + half_h),
                    (x_start, y_center + half_h),
                ]
            else:
                # Reverse: arrow tip at x_start, body cut inward
                verts = [
                    (x_start + arrow_w, y_center - half_h),
                    (x_end, y_center - half_h),
                    (x_end, y_center + half_h),
                    (x_start + arrow_w, y_center + half_h),
                    (x_start, y_center),
                ]
            poly = mpatches.Polygon(
                verts, closed=True,
                facecolor=read_color, edgecolor="none",
                alpha=self.alpha_reads, zorder=3,
            )
            ax.add_patch(poly)
        else:
            rect = mpatches.Rectangle(
                (x_start, y_center - half_h), width, h,
                facecolor=read_color, edgecolor="none",
                alpha=self.alpha_reads, zorder=3,
            )
            ax.add_patch(rect)

    def _draw_single_read(self, ax, read, region, y_center, row_height,
                          ref_seq, span, mismatch_counts=None):
        """Draw a single read with CIGAR-aware blocks, clipping, and indels."""
        h = row_height * 0.8

        # Custom color function takes priority
        if self.color_fn is not None:
            try:
                read_color = self.color_fn(read)
            except Exception:
                read_color = self.fill_reads
        elif self.color_by_strand:
            read_color = self.fill_reads_rev if read.is_reverse else self.fill_reads_fwd
        else:
            read_color = self.fill_reads

        if read.cigartuples is None:
            return

        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples:
            if op in (_CIGAR_M, _CIGAR_EQ, _CIGAR_X):
                self._draw_read_body(ax, read, region, ref_pos, ref_pos + length,
                                     y_center, h, read_color)
                if self.show_mismatches and ref_seq is not None:
                    self._draw_mismatches(
                        ax, read, ref_seq, ref_pos, query_pos, length,
                        region, y_center, h,
                        mismatch_counts=mismatch_counts,
                    )
                ref_pos += length
                query_pos += length

            elif op == _CIGAR_I:
                if (self.show_indels
                        and length > self.min_indel_size
                        and region.start <= ref_pos <= region.end):
                    self._draw_insertion_gv(ax, ref_pos, y_center, h, length, region)
                query_pos += length

            elif op == _CIGAR_D:
                x_start = max(ref_pos, region.start)
                x_end = min(ref_pos + length, region.end)
                if (self.show_indels
                        and length > self.min_indel_size
                        and x_end > x_start):
                    self._draw_deletion_gv(ax, x_start, x_end, y_center, h)
                ref_pos += length

            elif op == _CIGAR_N:
                x_start = max(ref_pos, region.start)
                x_end = min(ref_pos + length, region.end)
                if x_end > x_start:
                    ax.plot(
                        [x_start, x_end], [y_center, y_center],
                        color=self.col_gap, linewidth=0.5, zorder=2,
                    )
                ref_pos += length

            elif op == _CIGAR_S:
                if self.show_clipping and length >= 5:
                    clip_x = ref_pos
                    if region.start <= clip_x <= region.end:
                        clip_w = max(1, span / 500)
                        rect = mpatches.Rectangle(
                            (clip_x - clip_w / 2, y_center - h / 2),
                            clip_w, h,
                            facecolor=self.col_clipping,
                            edgecolor="none",
                            alpha=0.7, zorder=4,
                        )
                        ax.add_patch(rect)
                query_pos += length

    def _draw_insertion_gv(self, ax, ref_pos, y_center, h, length, region):
        """Draw insertion in genomeview style: purple I-beam with size label."""
        span = region.end - region.start
        ax.plot([ref_pos, ref_pos],
                [y_center - h * 0.6, y_center + h * 0.6],
                color=self.col_insertion, linewidth=1.5, zorder=4)
        cap_w = max(span * 0.0005, 1)
        ax.plot([ref_pos - cap_w, ref_pos + cap_w],
                [y_center - h * 0.6, y_center - h * 0.6],
                color=self.col_insertion, linewidth=1.0, zorder=4)
        ax.plot([ref_pos - cap_w, ref_pos + cap_w],
                [y_center + h * 0.6, y_center + h * 0.6],
                color=self.col_insertion, linewidth=1.0, zorder=4)
        draw_label = self.show_insertion_labels or length >= self.min_insertion_label_size
        if draw_label:
            ax.text(ref_pos, y_center, str(length),
                    ha="center", va="center", fontsize=3.5,
                    color="white", fontweight="bold", zorder=5)

    def _draw_deletion_gv(self, ax, x_start, x_end, y_center, h):
        """Draw deletion in genomeview style: white rect + black line."""
        width = x_end - x_start
        rect = mpatches.Rectangle(
            (x_start, y_center - h / 2), width, h,
            facecolor="white", edgecolor="none", zorder=3.5)
        ax.add_patch(rect)
        ax.plot([x_start, x_end], [y_center, y_center],
                color="black", linewidth=1.0, zorder=4)

    def _draw_mismatches(self, ax, read, ref_seq, ref_pos, query_pos,
                         length, region, y_center, h,
                         mismatch_counts=None):
        """Color-code mismatched bases in a read.

        When *mismatch_counts* is provided (quick-consensus mode), a
        mismatch is only drawn when the alt allele is supported by at
        least ``self.consensus_threshold`` fraction of reads at that
        position.
        """
        try:
            query_seq = read.query_sequence
            if query_seq is None:
                return
        except Exception:
            return

        for offset in range(length):
            rp = ref_pos + offset
            qp = query_pos + offset
            if rp < region.start or rp >= region.end:
                continue
            if rp - region.start >= len(ref_seq):
                continue
            if qp >= len(query_seq):
                continue

            ref_base = ref_seq[rp - region.start]
            query_base = query_seq[qp].upper()

            if query_base != ref_base and query_base in _NUC_COLORS:
                # Skip if quick-consensus is active and the alt allele
                # does not have sufficient support at this position.
                if mismatch_counts is not None:
                    if not mismatch_counts.query(
                        query_base, rp, threshold=self.consensus_threshold
                    ):
                        continue
                color = _NUC_COLORS[query_base]
                bp_width = max(1, (region.end - region.start) / 500)
                rect = mpatches.Rectangle(
                    (rp, y_center - h / 2), bp_width, h,
                    facecolor=color, edgecolor="none",
                    alpha=self.alpha_mismatch, zorder=4,
                )
                ax.add_patch(rect)

    def _draw_sashimi(self, ax, region: GenomicInterval) -> None:
        """Draw sashimi plot (splice junction arcs with read counts)."""
        try:
            aln, reads = self._fetch_reads(region)
        except Exception:
            return

        if not reads:
            return

        # Collect junction information from reads
        junctions = {}  # (start, end) -> count

        for read in reads:
            if read.is_unmapped or read.cigartuples is None:
                continue

            ref_pos = read.reference_start
            for op, length in read.cigartuples:
                if op == _CIGAR_N and length > 0:
                    j_start = ref_pos
                    j_end = ref_pos + length
                    key = (j_start, j_end)
                    junctions[key] = junctions.get(key, 0) + 1
                    ref_pos += length
                elif op in (_CIGAR_M, _CIGAR_EQ, _CIGAR_X, _CIGAR_D):
                    ref_pos += length

        aln.close()

        # Apply sashimi filter if set
        if self.sashimi_filter is not None:
            filtered = {}
            for (js, je), count in junctions.items():
                for _, frow in self.sashimi_filter.iterrows():
                    tol = self.sashimi_filter_tolerance
                    if (abs(js - frow["start"]) <= tol and
                            abs(je - frow["end"]) <= tol):
                        filtered[(js, je)] = count
                        break
            junctions = filtered

        # Filter by minimum score
        junctions = {k: v for k, v in junctions.items() if v >= self.sashimi_score}

        if not junctions:
            return

        # Draw arcs
        max_count = max(junctions.values())
        max_arc_height = self.sashimi_height

        for (js, je), count in junctions.items():
            mid = (js + je) / 2
            width = je - js
            arc_height = max_arc_height * (count / max_count)

            # Draw arc
            arc = Arc(
                (mid, 0.05), width, arc_height * 2,
                angle=0, theta1=0, theta2=180,
                color="#E63946", linewidth=1.0, alpha=0.7, zorder=3,
            )
            ax.add_patch(arc)

            # Label with count
            ax.text(
                mid, 0.05 + arc_height + 0.02, str(count),
                ha="center", va="bottom", fontsize=5,
                color="#E63946", zorder=5,
            )

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the region covered by the BAM file (from index)."""
        try:
            pysam = self._import_pysam()
            aln = pysam.AlignmentFile(self.filepath, "rb")
            refs = aln.references
            lengths = aln.lengths
            aln.close()
            if refs:
                return GenomicInterval(refs[0], 0, int(lengths[0]))
        except Exception:
            pass
        return None


class BAMCoverageTrack(Track):
    """Display per-base coverage from a BAM/CRAM file as a line or filled area.

    Unlike ``AlignmentsTrack(type="coverage")`` which shows coverage inside
    a pileup panel, ``BAMCoverageTrack`` is a standalone track that draws
    coverage as a continuous line or filled-area plot — similar to
    genomeview's ``BAMCoverageTrack``.

    Requires the optional ``pysam`` package.

    Parameters
    ----------
    filepath : str
        Path to a BAM or CRAM file.
    type : str
        Display style: ``"line"`` (default) or ``"fill"``.
    col : str
        Line / fill color.  Default is ``"#5B8DB8"``.
    alpha : float
        Transparency.  Default is 0.7.
    linewidth : float
        Line width for ``"line"`` mode.  Default is 1.0.
    transformation : callable, optional
        Function applied to coverage values before plotting
        (e.g. ``np.log1p``).
    name : str
        Track name.  Default is derived from the file name.
    height : float
        Relative track height.  Default is 1.5.
    display_params : dict, optional
        Additional display parameters.
    """

    def __init__(
        self,
        filepath: str,
        type: str = "line",
        col: str = "#5B8DB8",
        alpha: float = 0.7,
        linewidth: float = 1.0,
        transformation: Optional[Callable] = None,
        name: Optional[str] = None,
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        if name is None:
            import os
            name = os.path.basename(filepath).split(".")[0]
        dp = {"fontsize": 8}
        if display_params:
            dp.update(display_params)
        super().__init__(name=name, height=height, display_params=dp, **kwargs)

        self.filepath = filepath
        self.plot_type = type
        self.col = col
        self.alpha = alpha
        self.linewidth = linewidth
        self.transformation = transformation

    def _import_pysam(self):
        try:
            import pysam
            return pysam
        except ImportError:
            raise ImportError(
                "The 'pysam' package is required for BAMCoverageTrack. "
                "Install it with: pip install pysam"
            )

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw coverage line or filled area."""
        pysam = self._import_pysam()
        aln = pysam.AlignmentFile(self.filepath, "rb")
        chrom = match_chrom_format(region.chrom, aln.references)

        # Compute per-base coverage
        cov = aln.count_coverage(
            chrom, region.start, region.end, quality_threshold=0,
        )
        total = np.array(cov[0]) + np.array(cov[1]) + np.array(cov[2]) + np.array(cov[3])
        aln.close()

        if len(total) == 0:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        if self.transformation is not None:
            total = self.transformation(total.astype(float))

        positions = np.arange(region.start, region.start + len(total))
        max_cov = float(np.max(total)) if len(total) > 0 else 1.0
        if max_cov == 0:
            max_cov = 1.0

        ax.set_xlim(region.start, region.end)
        ax.set_ylim(0, max_cov * 1.1)

        if self.plot_type == "fill":
            ax.fill_between(
                positions, 0, total,
                color=self.col, alpha=self.alpha, step="mid", zorder=3,
            )
            ax.plot(
                positions, total,
                color=self.col, linewidth=self.linewidth * 0.5, zorder=4,
            )
        else:
            ax.plot(
                positions, total,
                color=self.col, linewidth=self.linewidth, alpha=self.alpha,
                zorder=3,
            )

        # Annotate max coverage
        ax.text(
            region.start, max_cov * 1.05,
            f"max={int(max_cov)}",
            ha="left", va="bottom", fontsize=6, color="#555555",
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the region covered by the BAM file (from index)."""
        try:
            pysam = self._import_pysam()
            aln = pysam.AlignmentFile(self.filepath, "rb")
            refs = aln.references
            lengths = aln.lengths
            aln.close()
            if refs:
                return GenomicInterval(refs[0], 0, int(lengths[0]))
        except Exception:
            pass
        return None
