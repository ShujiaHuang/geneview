"""
GeneRegionTrack - A track for displaying gene models.

Displays gene structures with exons as thick boxes, UTRs as thin boxes,
and introns as connecting lines. Supports transcript grouping, strand-aware
rendering, and gene/transcript/exon labels.

Ported from Gviz's GeneRegionTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union
from collections import defaultdict
from packaging.version import Version

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch

from ._base import StackedTrack, GenomicInterval
from ._stacking import compute_stacking, get_stack_heights, get_num_stacks


# Features that are drawn as thin boxes (UTRs and non-coding features)
_THIN_BOX_FEATURES = {
    "utr", "utr5", "utr3", "five_prime_utr", "three_prime_utr",
    "start_codon", "stop_codon", "ncrna", "nc_rna",
}

# Feature type normalization
_FEATURE_TYPE_MAP = {
    "cds": "CDS",
    "exon": "exon",
    "five_prime_utr": "utr5",
    "three_prime_utr": "utr3",
    "5utr": "utr5",
    "3utr": "utr3",
    "utr": "utr",
}


def _normalize_feature(feat: str) -> str:
    """Normalize feature type string."""
    if pd.isna(feat):
        return "exon"
    return _FEATURE_TYPE_MAP.get(str(feat).lower().strip(), str(feat).strip())


def _is_thin_feature(feat: str) -> bool:
    """Check if a feature should be drawn as a thin box (UTR)."""
    return _normalize_feature(feat).lower() in _THIN_BOX_FEATURES


class GeneRegionTrack(StackedTrack):
    """A track for displaying gene model annotations.

    Displays gene structures with exons as thick boxes, UTRs as thinner boxes,
    and introns as connecting lines between exons. Multiple transcripts per
    gene can be shown with stacking.

    Parameters
    ----------
    data : pd.DataFrame or str
        DataFrame with columns: chrom, start, end, and optionally:
        strand, transcript_id, gene_id, gene_name, exon_id, feature.
        The 'feature' column should contain 'CDS', 'exon', 'UTR', etc.
        If a string, interpreted as a GFF/GTF file path.
    stacking : str
        Stacking mode. Default is 'squish'.
    collapse_transcripts : str or bool
        How to handle multiple transcripts:
        - False: show all transcripts separately
        - True or 'gene': collapse to gene level
        - 'longest': show only longest transcript
        - 'shortest': show only shortest transcript
        - 'meta': union of all exons (meta-transcript)
        Default is False.
    show_id : str or None
        What identifier to show: 'gene', 'transcript', 'exon', or None.
        Default is 'gene'.
    thin_box_features : set, optional
        Feature types to draw as thin boxes. Default includes UTRs.
    name : str
        Track name. Default is "GeneRegion".
    height : float
        Relative track height. Default is 1.5.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    >>> import pandas as pd
    >>> from geneview.genometracks import GeneRegionTrack, plot_tracks, GenomicInterval
    >>> data = pd.DataFrame({
    ...     "chrom": ["chr7"] * 6,
    ...     "start": [1000, 1500, 2000, 3000, 3200, 4000],
    ...     "end":   [1200, 1700, 2200, 3100, 3400, 4200],
    ...     "strand": ["+"] * 6,
    ...     "feature": ["UTR", "CDS", "CDS", "CDS", "CDS", "UTR"],
    ...     "transcript_id": ["tx1"] * 3 + ["tx2"] * 3,
    ...     "gene_name": ["GeneA"] * 6,
    ... })
    >>> track = GeneRegionTrack(data)
    >>> plot_tracks([track], region=GenomicInterval("chr7", 500, 5000))
    """

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        stacking: str = "squish",
        collapse_transcripts: Union[str, bool] = False,
        show_id: Optional[str] = "gene",
        thin_box_features: Optional[set] = None,
        name: str = "GeneRegion",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        dp = {
            "col": "#3C5488",
            "fill": "#3C5488",
            "fill_utr": "#8DB4E2",
            "col_intron": "#888888",
            "lwd": 0.8,
            "fontsize": 8,
            "fontcolor": "#333333",
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=dp,
        )

        self.collapse_transcripts = collapse_transcripts
        self.show_id = show_id
        self._thin_features = thin_box_features or _THIN_BOX_FEATURES

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the gene region track.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        sub = self.subset(region)
        if sub is None or len(sub) == 0:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        # Ensure required columns
        if "strand" not in sub.columns:
            sub["strand"] = "*"
        if "feature" not in sub.columns:
            sub["feature"] = "exon"

        # Handle transcript collapsing
        if self.collapse_transcripts:
            sub = self._collapse(sub)

        # Determine grouping level for stacking
        group_col = self._get_group_column(sub)

        # Compute stacking at the group (transcript) level
        if group_col and group_col in sub.columns:
            groups = sub.groupby(group_col)
            group_starts = []
            group_ends = []
            group_keys = []
            for key, grp in groups:
                group_starts.append(grp["start"].min())
                group_ends.append(grp["end"].max())
                group_keys.append(key)

            stacks = compute_stacking(
                np.array(group_starts), np.array(group_ends),
                mode=self.stacking,
                min_distance=int(self.get_param("min_distance", 1)),
            )
            # Map group keys to stack rows
            group_to_stack = dict(zip(group_keys, stacks))
        else:
            # No grouping - stack individual features
            stacks = self.compute_stacks(sub)
            group_to_stack = None

        # Calculate layout
        if group_to_stack is not None:
            total_rows = get_num_stacks(stacks)
        else:
            total_rows = get_num_stacks(stacks)

        if total_rows == 0:
            total_rows = 1

        row_height = 0.85 / total_rows
        y_positions = np.linspace(1.0 - row_height / 2, row_height / 2, total_rows)

        ax.set_xlim(region.start, region.end)
        ax.set_ylim(-0.05, 1.05)

        fill_color = self.get_param("fill", "#3C5488")
        fill_utr = self.get_param("fill_utr", "#8DB4E2")
        intron_color = self.get_param("col_intron", "#888888")
        lwd = self.get_param("lwd", 0.8)
        fontsize = self.get_param("fontsize", 8)
        fontcolor = self.get_param("fontcolor", "#333333")
        alpha = self.get_param("alpha", 1.0)

        # Draw by group (transcript)
        if group_col and group_col in sub.columns:
            for group_key, grp in sub.groupby(group_col):
                stack_row = group_to_stack[group_key]
                y_center = y_positions[stack_row]
                self._draw_transcript(
                    ax, grp, region, y_center, row_height,
                    fill_color, fill_utr, intron_color, lwd, alpha,
                )
                # Draw gene/transcript label
                if self.show_id:
                    self._draw_label(
                        ax, grp, region, y_center, row_height,
                        fontsize, fontcolor,
                    )
        else:
            # No grouping - draw individual features
            for i, (idx, row) in enumerate(sub.iterrows()):
                stack_row = stacks[i] if len(stacks) > i else 0
                y_center = y_positions[stack_row]
                self._draw_single_feature(
                    ax, row, region, y_center, row_height,
                    fill_color, fill_utr, lwd, alpha,
                )

        ax.axis("off")

    def _get_group_column(self, data: pd.DataFrame) -> Optional[str]:
        """Determine which column to use for grouping."""
        for col in ["transcript_id", "transcript", "gene_id", "gene", "group"]:
            if col in data.columns:
                return col
        return None

    def _collapse(self, data: pd.DataFrame) -> pd.DataFrame:
        """Collapse transcripts based on settings."""
        if self.collapse_transcripts in (True, "gene"):
            # Collapse to gene level - merge all exons per gene
            gene_col = None
            for col in ["gene_id", "gene", "gene_name"]:
                if col in data.columns:
                    gene_col = col
                    break

            if gene_col is None:
                return data

            # Get unique gene extents
            result_rows = []
            for gene, grp in data.groupby(gene_col):
                grp = grp.copy()  # Bug 3 fix: avoid view modifications
                # Take the first strand
                strand = grp["strand"].iloc[0] if "strand" in grp.columns else "*"
                # Merge overlapping exons
                exon_ranges = grp[["start", "end"]].sort_values("start").values
                merged = [list(exon_ranges[0])]
                for s, e in exon_ranges[1:]:
                    if s <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], e)
                    else:
                        merged.append([s, e])

                for s, e in merged:
                    result_rows.append({
                        "chrom": grp["chrom"].iloc[0],
                        "start": s, "end": e,
                        "strand": strand,
                        "feature": "CDS",
                        gene_col: gene,
                    })
            return pd.DataFrame(result_rows)

        elif self.collapse_transcripts == "longest":
            tx_col = self._get_group_column(data)
            if tx_col is None:
                return data
            gene_col = None
            for col in ["gene_id", "gene", "gene_name"]:
                if col in data.columns:
                    gene_col = col
                    break
            if gene_col is None:
                return data

            result_rows = []
            for gene, grp in data.groupby(gene_col):
                grp = grp.copy()  # Bug 3 fix
                if tx_col not in grp.columns:
                    result_rows.append(grp)
                    continue
                # Find longest transcript
                _apply_kwargs = {} if Version(pd.__version__) < Version("2.2.0") else {"include_groups": False}
                tx_lengths = grp.groupby(tx_col).apply(
                    lambda g: g["end"].max() - g["start"].min(),
                    **_apply_kwargs,
                )
                longest_tx = tx_lengths.idxmax()
                result_rows.append(grp[grp[tx_col] == longest_tx])
            return pd.concat(result_rows) if result_rows else data

        elif self.collapse_transcripts == "meta":
            # Meta-transcript: union of all exon positions per gene
            gene_col = None
            for col in ["gene_id", "gene", "gene_name"]:
                if col in data.columns:
                    gene_col = col
                    break
            if gene_col is None:
                return data

            result_rows = []
            for gene, grp in data.groupby(gene_col):
                grp = grp.copy()  # Bug 3 fix
                strand = grp["strand"].iloc[0] if "strand" in grp.columns else "*"
                chrom = grp["chrom"].iloc[0]

                # Merge all exon ranges into a meta-transcript
                exon_ranges = grp[["start", "end"]].sort_values("start").values
                merged = [list(exon_ranges[0])]
                for s, e in exon_ranges[1:]:
                    if s <= merged[-1][1]:
                        merged[-1][1] = max(merged[-1][1], e)
                    else:
                        merged.append([s, e])

                # Create a single transcript_id for the meta-transcript
                tx_col = self._get_group_column(grp)
                tx_name = f"{gene}_meta" if tx_col else gene

                for s, e in merged:
                    row = {
                        "chrom": chrom,
                        "start": s, "end": e,
                        "strand": strand,
                        "feature": "CDS",
                        gene_col: gene,
                    }
                    if tx_col:
                        row[tx_col] = tx_name
                    result_rows.append(row)
            return pd.DataFrame(result_rows) if result_rows else data

        return data

    def _draw_transcript(self, ax, grp, region, y_center, row_height,
                         fill_color, fill_utr, intron_color, lwd, alpha):
        """Draw a complete transcript with intron lines and exon boxes."""
        # Sort exons by position
        grp = grp.sort_values("start")

        strand = grp["strand"].iloc[0] if "strand" in grp.columns else "*"
        thick_frac = 0.7  # Fraction of row height for CDS
        thin_frac = 0.4   # Fraction of row height for UTR

        # Draw intron connecting line first
        if len(grp) > 1:
            x_intron_start = grp["start"].iloc[0] + (grp["end"].iloc[0] - grp["start"].iloc[0]) / 2
            x_intron_end = grp["start"].iloc[-1] + (grp["end"].iloc[-1] - grp["start"].iloc[-1]) / 2
            ax.plot([x_intron_start, x_intron_end],
                    [y_center, y_center],
                    color=intron_color, linewidth=lwd * 0.8, alpha=alpha, zorder=2)

        # Draw each exon/feature
        for _, row in grp.iterrows():
            x_start = max(row["start"], region.start)
            x_end = min(row["end"], region.end)
            width = x_end - x_start
            if width <= 0:
                continue

            feat = _normalize_feature(row.get("feature", "exon"))
            is_thin = _is_thin_feature(str(row.get("feature", "exon")))

            if is_thin:
                h = row_height * thin_frac
                color = fill_utr
            else:
                h = row_height * thick_frac
                color = fill_color

            # Draw exon box
            rect = FancyBboxPatch(
                (x_start, y_center - h / 2), width, h,
                boxstyle="round,pad=0",
                facecolor=color, edgecolor="#333333",
                linewidth=lwd * 0.5, alpha=alpha, zorder=3,
            )
            ax.add_patch(rect)

            # Draw strand direction arrows inside exon
            if strand in ("+", "-"):
                self._draw_strand_chevrons(ax, x_start, x_end, y_center, h, strand, region)

    def _draw_single_feature(self, ax, row, region, y_center, row_height,
                             fill_color, fill_utr, lwd, alpha):
        """Draw a single feature without transcript context."""
        x_start = max(row["start"], region.start)
        x_end = min(row["end"], region.end)
        width = x_end - x_start
        if width <= 0:
            return

        is_thin = _is_thin_feature(str(row.get("feature", "exon")))
        h = row_height * (0.4 if is_thin else 0.7)
        color = fill_utr if is_thin else fill_color

        rect = FancyBboxPatch(
            (x_start, y_center - h / 2), width, h,
            boxstyle="round,pad=0",
            facecolor=color, edgecolor="#333333",
            linewidth=lwd * 0.5, alpha=alpha, zorder=3,
        )
        ax.add_patch(rect)

    def _draw_strand_chevrons(self, ax, x_start, x_end, y_center, height, strand, region):
        """Draw small chevron arrows inside exon boxes to indicate strand."""
        span = region.end - region.start
        box_width = x_end - x_start
        if box_width < span * 0.005:
            return  # Box too small for chevrons

        chevron_size = min(height * 0.3, box_width * 0.1)
        if chevron_size < 1:
            return

        # Draw 1-3 chevrons depending on box width
        n_chevrons = min(3, max(1, int(box_width / (chevron_size * 4))))
        mid_x = (x_start + x_end) / 2
        spacing = chevron_size * 2

        for i in range(n_chevrons):
            cx = mid_x + (i - (n_chevrons - 1) / 2) * spacing

            if strand == "+":
                verts = [
                    (cx - chevron_size / 2, y_center - chevron_size / 2),
                    (cx + chevron_size / 2, y_center),
                    (cx - chevron_size / 2, y_center + chevron_size / 2),
                ]
            else:
                verts = [
                    (cx + chevron_size / 2, y_center - chevron_size / 2),
                    (cx - chevron_size / 2, y_center),
                    (cx + chevron_size / 2, y_center + chevron_size / 2),
                ]

            line = mpatches.Polygon(
                verts, closed=False,
                facecolor="none", edgecolor="white",
                linewidth=1.0, alpha=0.7, zorder=4,
            )
            ax.add_patch(line)

    def _draw_label(self, ax, grp, region, y_center, row_height,
                    fontsize, fontcolor):
        """Draw a gene or transcript label."""
        if self.show_id == "gene":
            for col in ["gene_name", "gene", "gene_id", "symbol"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label) and label != "unknown":
                        break
            else:
                return
        elif self.show_id == "transcript":
            for col in ["transcript_id", "transcript"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label):
                        break
            else:
                return
        else:
            return

        # Position label below the transcript
        x_pos = (grp["start"].min() + grp["end"].max()) / 2
        x_pos = np.clip(x_pos, region.start, region.end)
        y_pos = y_center - row_height * 0.55

        ax.text(x_pos, y_pos, str(label),
                ha="center", va="top", fontsize=fontsize,
                color=fontcolor, fontstyle="italic",
                clip_on=True, zorder=5)
