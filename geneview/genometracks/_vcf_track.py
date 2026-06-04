"""
VCFTrack - A track for displaying genomic variants from VCF/BCF files.

Displays variants as colored rectangles, stacked when overlapping.
Supports custom coloring via a ``color_fn`` callback, similar to
genomeview's VCFTrack extension example.

Ported from the genomeview examples notebook (VCFTrack extension).
"""

from typing import Optional, Callable, Dict, Any, List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from ._base import Track, GenomicInterval
from ._utils import match_chrom_format


# Default nucleotide colors for variant alt alleles
_NUC_COLORS: Dict[str, str] = {
    "A": "#0072B2",   # blue
    "C": "#E69F00",   # orange
    "G": "#009E73",   # green
    "T": "#D55E00",   # vermillion
    "N": "#999999",   # gray
}


def _default_variant_color(variant) -> str:
    """Color a variant by its first alt allele nucleotide."""
    try:
        alt = str(variant.alts[0]).upper()
        return _NUC_COLORS.get(alt, "#999999")
    except (IndexError, TypeError, AttributeError):
        return "#999999"


class VCFTrack(Track):
    """Display genomic variants from a VCF or BCF file.

    Variants are drawn as colored rectangles at their genomic positions.
    Overlapping variants are stacked vertically.  A custom ``color_fn``
    can be provided to color variants by alt allele, genotype, quality,
    or any other property accessible on a ``pysam.VariantRecord``.

    Requires the optional ``pysam`` package.

    Parameters
    ----------
    filepath : str
        Path to a VCF or BCF file (must be indexed with ``.tbi`` or
        ``.csi``).
    color_fn : callable, optional
        A function ``f(pysam.VariantRecord) -> str`` that returns a
        matplotlib color string for each variant.  Default colors by
        first alt allele nucleotide (A=blue, C=orange, G=green,
        T=vermillion).
    min_variant_width : float
        Minimum width (in data coordinates) for a variant rectangle,
        expressed as a fraction of the visible region width.  Ensures
        variants remain visible when zoomed out.  Default is 0.002
        (0.2% of the region).
    row_height : float
        Height of each variant row as a fraction of the axes height.
        Default is 0.12.
    margin_y : float
        Vertical gap between stacked rows as a fraction of axes height.
        Default is 0.02.
    alpha : float
        Opacity of variant rectangles.  Default is 0.85.
    name : str
        Track name.  Default is "Variants".
    height : float
        Relative track height.  Default is 1.0.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> from geneview.genometracks import VCFTrack              # doctest: +SKIP
    >>> vtrack = VCFTrack("sample.vcf.gz", name="SNPs")         # doctest: +SKIP
    >>> # Custom coloring by variant quality
    >>> def color_by_qual(v):                                    # doctest: +SKIP
    ...     return "red" if v.qual and v.qual < 30 else "blue"
    >>> vtrack = VCFTrack("sample.vcf.gz", color_fn=color_by_qual)  # doctest: +SKIP
    """

    def __init__(
        self,
        filepath: str,
        color_fn: Optional[Callable] = None,
        min_variant_width: float = 0.002,
        row_height: float = 0.12,
        margin_y: float = 0.02,
        alpha: float = 0.85,
        name: str = "Variants",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        dp = {"fontsize": 7}
        if display_params:
            dp.update(display_params)
        super().__init__(name=name, height=height, display_params=dp, **kwargs)

        self.filepath = filepath
        self.color_fn = color_fn or _default_variant_color
        self.min_variant_width = min_variant_width
        self.row_height = row_height
        self.margin_y = margin_y
        self.alpha = alpha

    def _import_pysam(self):
        """Import pysam or raise a clear error."""
        try:
            import pysam
            return pysam
        except ImportError:
            raise ImportError(
                "The 'pysam' package is required for VCFTrack. "
                "Install it with: pip install pysam"
            )

    def _fetch_variants(self, region: GenomicInterval):
        """Fetch variants overlapping the region.

        Returns a list of pysam.VariantRecord objects.
        """
        pysam = self._import_pysam()
        vcf = pysam.VariantFile(self.filepath)
        try:
            chrom = match_chrom_format(region.chrom, list(vcf.header.contigs))
            records = list(vcf.fetch(chrom, region.start, region.end))
            return records
        except Exception:
            return []
        finally:
            vcf.close()

    def _stack_variants(self, variants, region: GenomicInterval) -> List[int]:
        """Assign each variant to a row to avoid overlaps.

        Uses a simple greedy interval scheduling algorithm.

        Returns
        -------
        list of int
            Row index for each variant.
        """
        if not variants:
            return []

        # Sort by start position
        indexed = sorted(enumerate(variants), key=lambda x: x[1].start)
        row_ends: List[int] = []  # End position of last variant in each row
        assignments = [0] * len(variants)

        for orig_idx, var in indexed:
            placed = False
            for row_idx, row_end in enumerate(row_ends):
                if var.start >= row_end:
                    row_ends[row_idx] = var.stop + 1
                    assignments[orig_idx] = row_idx
                    placed = True
                    break
            if not placed:
                assignments[orig_idx] = len(row_ends)
                row_ends.append(var.stop + 1)

        return assignments

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the VCF track.

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

        variants = self._fetch_variants(region)
        if not variants:
            ax.text(
                0.5, 0.5, "No variants",
                ha="center", va="center",
                fontsize=8, color="#999999",
                transform=ax.transAxes,
            )
            return

        rows = self._stack_variants(variants, region)
        n_rows = max(rows) + 1 if rows else 1

        span = region.end - region.start
        min_w = span * self.min_variant_width

        for var, row_idx in zip(variants, rows):
            var_start = var.start
            var_end = var.stop + 1  # pysam VCF stop is inclusive
            var_width = max(var_end - var_start, min_w)

            # Clamp to region
            x_start = max(var_start, region.start)
            x_end = min(var_start + var_width, region.end)
            if x_end <= x_start:
                continue

            # Compute y position
            y_top = 1.0 - row_idx * (self.row_height + self.margin_y)
            y_bottom = y_top - self.row_height
            if y_bottom < 0:
                continue

            color = self.color_fn(var)

            rect = mpatches.Rectangle(
                (x_start, y_bottom), x_end - x_start, self.row_height,
                facecolor=color, edgecolor="none",
                alpha=self.alpha, zorder=3,
            )
            ax.add_patch(rect)

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the region span of the VCF file's first contig."""
        pysam = self._import_pysam()
        try:
            vcf = pysam.VariantFile(self.filepath)
            contigs = list(vcf.header.contigs)
            if contigs:
                chrom = contigs[0]
                length = vcf.header.contigs[chrom].length
                if length and length > 0:
                    return GenomicInterval(chrom, 0, length)
            vcf.close()
        except Exception:
            pass
        return None

    def __repr__(self):
        return f"VCFTrack(filepath='{self.filepath}', name='{self.name}')"
