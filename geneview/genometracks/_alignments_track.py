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

from ._base import StackedTrack, GenomicInterval


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

# Nucleotide colors (matching SequenceTrack)
_NUC_COLORS = {
    "A": "#009E73", "C": "#0072B2", "G": "#E69F00", "T": "#D55E00",
    "N": "#999999",
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
        col_mates: str = "lightblue",
        col_gap: str = "lightgray",
        col_deletion: str = "red",
        col_insertion: str = "blue",
        fill_coverage: str = "#5B8DB8",
        fill_reads: str = "#BDBDBD",
        alpha_reads: float = 0.8,
        alpha_mismatch: float = 0.6,
        transformation: Optional[Callable] = None,
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
        """Fetch reads overlapping the region. Returns pysam AlignmentFile and iterator."""
        pysam = self._import_pysam()
        aln = pysam.AlignmentFile(self.filepath, "rb")
        reads = list(aln.fetch(region.chrom, region.start, region.end))
        return aln, reads

    def _compute_coverage(self, region: GenomicInterval) -> np.ndarray:
        """Compute per-base coverage across the region."""
        pysam = self._import_pysam()
        aln = pysam.AlignmentFile(self.filepath, "rb")
        cov = aln.count_coverage(
            region.chrom, region.start, region.end,
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
            return fa.fetch(region.chrom, region.start, region.end).upper()
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

        for i, read in enumerate(reads):
            if i >= len(stacks):
                break
            stack_row = stacks[i]
            y_center = pileup_bottom + (pileup_top - pileup_bottom) - (stack_row + 0.5) * row_height

            # Draw read body from aligned pairs
            self._draw_single_read(
                ax, read, region, y_center, row_height, ref_seq, span,
            )

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

    def _draw_single_read(self, ax, read, region, y_center, row_height,
                          ref_seq, span):
        """Draw a single read with CIGAR-aware blocks."""
        h = row_height * 0.8
        read_color = self.fill_reads

        # Iterate over CIGAR blocks
        if read.cigartuples is None:
            return

        ref_pos = read.reference_start
        query_pos = 0

        for op, length in read.cigartuples:
            if op in (_CIGAR_M, _CIGAR_EQ, _CIGAR_X):
                # Alignment match
                x_start = max(ref_pos, region.start)
                x_end = min(ref_pos + length, region.end)
                if x_end > x_start:
                    rect = mpatches.Rectangle(
                        (x_start, y_center - h / 2), x_end - x_start, h,
                        facecolor=read_color, edgecolor="none",
                        alpha=self.alpha_reads, zorder=3,
                    )
                    ax.add_patch(rect)

                    # Mismatch coloring
                    if self.show_mismatches and ref_seq is not None:
                        self._draw_mismatches(
                            ax, read, ref_seq, ref_pos, query_pos, length,
                            region, y_center, h,
                        )

                ref_pos += length
                query_pos += length

            elif op == _CIGAR_I:
                # Insertion
                if self.show_indels and region.start <= ref_pos <= region.end:
                    ax.plot(
                        [ref_pos, ref_pos],
                        [y_center - h * 0.6, y_center + h * 0.6],
                        color=self.col_insertion, linewidth=1.5, zorder=4,
                    )
                query_pos += length

            elif op == _CIGAR_D:
                # Deletion
                x_start = max(ref_pos, region.start)
                x_end = min(ref_pos + length, region.end)
                if self.show_indels and x_end > x_start:
                    ax.plot(
                        [x_start, x_end],
                        [y_center, y_center],
                        color=self.col_deletion, linewidth=2.0, zorder=4,
                    )
                ref_pos += length

            elif op == _CIGAR_N:
                # Skipped region (intron) — draw thin line
                x_start = max(ref_pos, region.start)
                x_end = min(ref_pos + length, region.end)
                if x_end > x_start:
                    ax.plot(
                        [x_start, x_end],
                        [y_center, y_center],
                        color=self.col_gap, linewidth=0.5, zorder=2,
                    )
                ref_pos += length

            elif op == _CIGAR_S:
                query_pos += length
            # H, P: no consumption

    def _draw_mismatches(self, ax, read, ref_seq, ref_pos, query_pos,
                         length, region, y_center, h):
        """Color-code mismatched bases in a read."""
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
