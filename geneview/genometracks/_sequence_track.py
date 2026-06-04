"""
SequenceTrack - A track for displaying nucleotide sequences.

Displays genomic DNA sequence as colored letters, boxes, or a single line
depending on the zoom level.  Supports complement display, 5'->3' direction
arrows, and custom nucleotide colors.

Ported from Gviz's SequenceTrack-class.R.
"""

from typing import Optional, Dict, Any, Union

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ._base import Track, GenomicInterval


# Default nucleotide colors (Gviz/biovizBase-inspired)
_DEFAULT_NUC_COLORS = {
    "A": "#009E73",   # green
    "C": "#0072B2",   # blue
    "G": "#E69F00",   # yellow/orange
    "T": "#D55E00",   # red
    "U": "#D55E00",   # red (RNA)
    "N": "#999999",   # gray
}

# Complement mapping
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT)[::-1]


class SequenceTrack(Track):
    """Display nucleotide sequence along the genomic axis.

    At high zoom (small window), individual colored letters are drawn.
    At medium zoom, colored boxes per nucleotide are shown.
    At low zoom (large window), a single colored line is drawn.

    Parameters
    ----------
    sequence : str or pd.DataFrame, optional
        Nucleotide sequence string or DataFrame with columns
        ``chrom``, ``start``, ``end``, ``seq``.
    fasta_path : str, optional
        Path to an indexed FASTA file (requires pysam).
    twobit_path : str, optional
        Path to a 2bit file (requires py2bit).
    chromosome : str, optional
        Chromosome name (used with *fasta_path* or *twobit_path*).
    complement : bool
        If True, show the reverse complement.  Default is False.
    add53 : bool
        If True, draw a 5'→3' direction arrow.  Default is False.
    fontcolor : dict, optional
        Mapping of nucleotide → color.  Default uses Gviz-like palette.
    noLetters : bool
        If True, draw colored boxes instead of letters even at high zoom.
        Default is False.
    cex : float
        Font size expansion factor.  Default is 1.0.
    min_width : int
        Minimum window width (bp) below which letters are drawn.
        Default is 200.
    name : str
        Track name.  Default is "Sequence".
    height : float
        Relative track height.  Default is 0.5.
    display_params : dict, optional
        Additional display parameters.
    """

    def __init__(
        self,
        sequence: Optional[Union[str, Any]] = None,
        fasta_path: Optional[str] = None,
        twobit_path: Optional[str] = None,
        chromosome: Optional[str] = None,
        complement: bool = False,
        add53: bool = False,
        fontcolor: Optional[Dict[str, str]] = None,
        noLetters: bool = False,
        cex: float = 1.0,
        min_width: int = 200,
        name: str = "Sequence",
        height: float = 0.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        dp = {"fontsize": 8}
        if display_params:
            dp.update(display_params)
        super().__init__(name=name, height=height, display_params=dp)

        self._sequence = sequence
        self._fasta_path = fasta_path
        self._twobit_path = twobit_path
        self._chromosome = chromosome
        self.complement = complement
        self.add53 = add53
        self.noLetters = noLetters
        self.cex = cex
        self.min_width = min_width

        self._nuc_colors = dict(_DEFAULT_NUC_COLORS)
        if fontcolor:
            self._nuc_colors.update(fontcolor)

    # ------------------------------------------------------------------
    # Sequence loading
    # ------------------------------------------------------------------

    def _load_sequence(self, region: GenomicInterval) -> str:
        """Load the sequence for the given region."""
        if isinstance(self._sequence, str):
            seq = self._sequence
        elif self._fasta_path is not None:
            from ._io import read_fasta
            seq = read_fasta(self._fasta_path, region)
        elif self._twobit_path is not None:
            from ._io import read_2bit
            seq = read_2bit(self._twobit_path, region)
        else:
            # Try DataFrame
            try:
                import pandas as pd
                if isinstance(self._sequence, pd.DataFrame):
                    df = self._sequence
                    mask = (
                        (df["chrom"] == region.chrom)
                        & (df["start"] < region.end)
                        & (df["end"] > region.start)
                    )
                    sub = df[mask]
                    if len(sub) > 0:
                        seq = "".join(sub["seq"].astype(str).tolist())
                    else:
                        seq = ""
                else:
                    seq = ""
            except Exception:
                seq = ""

        if self.complement:
            seq = _reverse_complement(seq)

        return seq.upper()

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the track's genomic extent if known."""
        if self._chromosome:
            return GenomicInterval(self._chromosome, 0, 1_000_000)
        return None

    # ------------------------------------------------------------------
    # Drawing
    # ------------------------------------------------------------------

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the sequence track.

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

        seq = self._load_sequence(region)
        if not seq:
            return

        span = region.end - region.start
        bp_per_pixel = span / max(1, ax.get_window_extent().width if hasattr(ax, "get_window_extent") else 800)

        if span <= self.min_width and not self.noLetters:
            self._draw_letters(ax, seq, region)
        elif span <= self.min_width * 10:
            self._draw_boxes(ax, seq, region)
        else:
            self._draw_line(ax, seq, region)

        # 5'->3' direction arrow
        if self.add53:
            self._draw_53_arrow(ax, region)

    def _draw_letters(self, ax, seq: str, region: GenomicInterval) -> None:
        """Draw individual colored nucleotide letters."""
        span = region.end - region.start
        bp_width = span / max(len(seq), 1)
        fontsize = self.get_param("fontsize", 8) * self.cex

        for i, nuc in enumerate(seq):
            x = region.start + (i + 0.5) * bp_width
            color = self._nuc_colors.get(nuc, "#999999")
            ax.text(
                x, 0.5, nuc,
                ha="center", va="center",
                fontsize=fontsize, color=color,
                fontfamily="monospace", fontweight="bold",
                clip_on=True, zorder=3,
            )

    def _draw_boxes(self, ax, seq: str, region: GenomicInterval) -> None:
        """Draw colored boxes per nucleotide."""
        span = region.end - region.start
        bp_width = span / max(len(seq), 1)

        for i, nuc in enumerate(seq):
            x = region.start + i * bp_width
            color = self._nuc_colors.get(nuc, "#999999")
            rect = mpatches.FancyBboxPatch(
                (x, 0.2), bp_width * 0.95, 0.6,
                boxstyle="round,pad=0",
                facecolor=color, edgecolor="none",
                alpha=0.8, zorder=3,
            )
            ax.add_patch(rect)

    def _draw_line(self, ax, seq: str, region: GenomicInterval) -> None:
        """Draw a single colored line representing the sequence."""
        # Compute overall GC content to color the line
        if len(seq) > 0:
            gc = sum(1 for c in seq if c in ("G", "C")) / len(seq)
        else:
            gc = 0.5

        # Interpolate color based on GC content
        color = "#0072B2" if gc > 0.5 else "#E69F00"
        ax.axhline(y=0.5, xmin=0, xmax=1, color=color,
                   linewidth=3, alpha=0.8, zorder=3)

    def _draw_53_arrow(self, ax, region: GenomicInterval) -> None:
        """Draw a 5'->3' direction arrow."""
        span = region.end - region.start
        arrow_len = span * 0.05
        x_start = region.end - arrow_len * 2
        x_end = region.end - arrow_len

        ax.annotate(
            "",
            xy=(x_end, 0.9),
            xytext=(x_start, 0.9),
            arrowprops=dict(
                arrowstyle="->", color="#333333", lw=1.2,
            ),
            zorder=5,
        )
        ax.text(
            x_start - arrow_len * 0.2, 0.9, "5'",
            ha="right", va="center", fontsize=6, color="#333333", zorder=5,
        )
        ax.text(
            x_end + arrow_len * 0.2, 0.9, "3'",
            ha="left", va="center", fontsize=6, color="#333333", zorder=5,
        )
