"""
GeneRegionTrack - A track for displaying gene models.

Displays gene structures in UCSC Genome Browser style:
- Coding exons (CDS) as thick solid blocks
- UTRs as thinner blocks (~half CDS height)
- Introns as thin horizontal connecting lines with directional chevron arrows
- Gene labels positioned to the LEFT of the transcript

Supports transcript grouping, strand-aware coloring, and gene/transcript/exon labels.

Ported from Gviz's GeneRegionTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union
from collections import defaultdict
from packaging.version import Version

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    >>> _ = plot_tracks([track], region=GenomicInterval("chr7", 500, 5000))
    """
    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        stacking: str = "squish",
        collapse_transcripts: Union[str, bool] = False,
        show_id: Optional[str] = "gene",
        thin_box_features: Optional[set] = None,
        exon_annotation: Optional[str] = None,
        gene_symbols: bool = False,
        transcript_annotation: Optional[str] = None,
        name: str = "GeneRegion",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        dp = {
            "col": "none",                   # genomeview: no edge stroke
            "fill": "#E89E9D",               # genomeview: forward strand color
            "fill_rev": "#8C8FCE",           # genomeview: reverse strand color
            "fill_utr": None,                # genomeview: same color, just thinner
            "col_intron": None,              # genomeview: same color as exon
            "lwd": 0.8,
            "fontsize": 8,
            "fontcolor": "#333333",
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=dp, **kwargs,
        )

        self.collapse_transcripts = collapse_transcripts
        self.show_id = show_id
        self._thin_features = thin_box_features or _THIN_BOX_FEATURES
        self.exon_annotation = exon_annotation
        self.gene_symbols = gene_symbols
        self.transcript_annotation = transcript_annotation

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

        # Remove redundant exon features when CDS/UTR features exist
        # (GFF/GTF includes both exon + CDS/UTR rows at the same positions)
        group_col = self._get_group_column(sub)
        if group_col and group_col in sub.columns:
            dedup_parts = []
            for _, grp in sub.groupby(group_col):
                dedup_parts.append(self._deduplicate_features(grp))
            sub = pd.concat(dedup_parts)
        else:
            sub = self._deduplicate_features(sub)

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

        # Support reverse stacking order
        if self.reverse_stacking:
            y_positions = y_positions[::-1]

        ax.set_xlim(region.start, region.end)
        ax.set_ylim(-0.05, 1.05)

        fill_color = self.get_param("fill", "#E89E9D")
        fill_rev = self.get_param("fill_rev", "#8C8FCE")
        fill_utr = self.get_param("fill_utr", None)  # None = same as CDS
        edge_color = self.get_param("col", "none")
        intron_color = self.get_param("col_intron", None)  # None = same as exon
        self._fill_rev = fill_rev
        self._fill_utr_override = fill_utr
        self._intron_color_override = intron_color
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
                    fill_color, fill_utr, edge_color, intron_color, lwd, alpha,
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
                    fill_color, fill_utr, edge_color, lwd, alpha,
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

        elif self.collapse_transcripts == "shortest":
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
                grp = grp.copy()
                if tx_col not in grp.columns:
                    result_rows.append(grp)
                    continue
                _apply_kwargs = {} if Version(pd.__version__) < Version("2.2.0") else {"include_groups": False}
                tx_lengths = grp.groupby(tx_col).apply(
                    lambda g: g["end"].max() - g["start"].min(),
                    **_apply_kwargs,
                )
                shortest_tx = tx_lengths.idxmin()
                result_rows.append(grp[grp[tx_col] == shortest_tx])
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

    @staticmethod
    def _deduplicate_features(grp: pd.DataFrame) -> pd.DataFrame:
        """Remove redundant exon rows when CDS/UTR features exist.

        GFF/GTF files typically include both 'exon' rows and their
        constituent 'CDS'/'UTR' rows at the same coordinates.  genomeview's
        BEDTrack model only uses exon boundaries with coding_start/end to
        determine thick vs thin, so drawing both would produce double boxes.

        This method filters out 'exon' features that overlap with CDS/UTR
        features in the same transcript, keeping only the sub-features.
        When no CDS/UTR features exist, exon features are kept as-is.
        """
        features = grp["feature"].astype(str).str.lower().str.strip() if "feature" in grp.columns else pd.Series(["exon"] * len(grp))
        has_sub = features.isin({"cds", "utr", "utr5", "utr3",
                                  "five_prime_utr", "three_prime_utr",
                                  "5utr", "3utr", "ncrna", "nc_rna"}).any()
        if not has_sub:
            return grp

        # Drop 'exon' rows whose (start, end) overlaps a sub-feature
        sub_mask = features.isin({"cds", "utr", "utr5", "utr3",
                                   "five_prime_utr", "three_prime_utr",
                                   "5utr", "3utr", "ncrna", "nc_rna",
                                   "start_codon", "stop_codon"})
        sub_ranges = set(zip(grp.loc[sub_mask, "start"], grp.loc[sub_mask, "end"]))

        keep = []
        for idx, row in grp.iterrows():
            feat = str(row.get("feature", "exon")).lower().strip()
            if feat == "exon":
                # Keep exon only if no sub-feature covers this exact range
                if (row["start"], row["end"]) not in sub_ranges:
                    keep.append(idx)
            else:
                keep.append(idx)
        return grp.loc[keep]

    def _draw_transcript(self, ax, grp, region, y_center, row_height,
                         fill_color, fill_utr, edge_color, intron_color, lwd, alpha):
        """Draw a complete transcript (UCSC Genome Browser style).

        UCSC style:
        - CDS: thick solid blocks (no edge stroke by default)
        - UTR: thinner blocks (~half CDS height)
        - Introns: thin horizontal connecting lines with directional chevron arrows
        - Single strand-based color per transcript
        - Gene label to the LEFT of the transcript
        Drawing order: thin first, then thick on top.
        """
        grp = self._deduplicate_features(grp).sort_values("start")
        strand = grp["strand"].iloc[0] if "strand" in grp.columns else "*"

        # UCSC: strand-based coloring
        if strand == "-" and hasattr(self, '_fill_rev'):
            tx_color = self._fill_rev
        else:
            tx_color = fill_color

        # Override colors if explicitly set
        effective_intron = self._intron_color_override if hasattr(self, '_intron_color_override') and self._intron_color_override else tx_color
        effective_utr = self._fill_utr_override if hasattr(self, '_fill_utr_override') and self._fill_utr_override else tx_color
        effective_edge = edge_color if edge_color and edge_color != "none" else "none"

        # UCSC proportions: CDS ~full, UTR ~half of CDS
        thick_frac = 0.7   # CDS
        thin_frac = 0.35   # UTR (slightly taller than genomeview)

        # --- Step 1: Draw intron connecting lines with directional chevron arrows ---
        if len(grp) > 1:
            exon_starts = grp["start"].values
            exon_ends = grp["end"].values
            for i in range(len(grp) - 1):
                intron_start = exon_ends[i]
                intron_end = exon_starts[i + 1]
                x_intron_start = max(intron_start, region.start)
                x_intron_end = min(intron_end, region.end)
                if x_intron_end <= x_intron_start:
                    continue
                # Thin intron line
                ax.plot([x_intron_start, x_intron_end],
                        [y_center, y_center],
                        color=effective_intron, linewidth=min(lwd, 0.8),
                        alpha=alpha, zorder=2, solid_capstyle='butt')
                # Directional chevron arrows along the intron
                if strand in ("+", "-"):
                    self._draw_intron_arrows(
                        ax, x_intron_start, x_intron_end, y_center,
                        row_height, strand, effective_intron, lwd, alpha, region,
                    )

        # --- Step 2: Draw thin features first (UTR/non-coding) ---
        for _, row in grp.iterrows():
            feat = str(row.get("feature", "exon"))
            if not _is_thin_feature(feat):
                continue
            x_start = max(row["start"], region.start)
            x_end = min(row["end"], region.end)
            width = x_end - x_start
            if width <= 0:
                continue
            h = row_height * thin_frac
            rect = mpatches.Rectangle(
                (x_start, y_center - h / 2), width, h,
                facecolor=effective_utr, edgecolor=effective_edge,
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(rect)
            if self.exon_annotation:
                self._draw_exon_label(ax, row, x_start, x_end, y_center, h, region)

        # --- Step 3: Draw thick features on top (CDS/exon) ---
        for _, row in grp.iterrows():
            feat = str(row.get("feature", "exon"))
            if _is_thin_feature(feat):
                continue
            x_start = max(row["start"], region.start)
            x_end = min(row["end"], region.end)
            width = x_end - x_start
            if width <= 0:
                continue
            h = row_height * thick_frac
            rect = mpatches.Rectangle(
                (x_start, y_center - h / 2), width, h,
                facecolor=tx_color, edgecolor=effective_edge,
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(rect)
            if self.exon_annotation:
                self._draw_exon_label(ax, row, x_start, x_end, y_center, h, region)

    def _draw_single_feature(self, ax, row, region, y_center, row_height,
                             fill_color, fill_utr, edge_color, lwd, alpha):
        """Draw a single feature without transcript context.

        Note: called only when there is no transcript grouping.  Redundant
        exon rows should already have been filtered by the caller.
        """
        x_start = max(row["start"], region.start)
        x_end = min(row["end"], region.end)
        width = x_end - x_start
        if width <= 0:
            return

        is_thin = _is_thin_feature(str(row.get("feature", "exon")))
        h = row_height * (0.35 if is_thin else 0.7)
        effective_utr = self._fill_utr_override if hasattr(self, '_fill_utr_override') and self._fill_utr_override else fill_color
        color = effective_utr if is_thin else fill_color
        effective_edge = edge_color if edge_color and edge_color != "none" else "none"

        rect = mpatches.Rectangle(
            (x_start, y_center - h / 2), width, h,
            facecolor=color, edgecolor=effective_edge,
            linewidth=0, alpha=alpha, zorder=3,
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

    def _draw_intron_arrows(self, ax, x_start, x_end, y_center, row_height,
                            strand, color, lwd, alpha, region):
        """Draw directional chevron arrows along an intron line (genomeview style).

        Matches genomeview's arrow() SVG path geometry:
          right: (x-2.5s, y-5s) -> (x+2.5s, y) -> (x-2.5s, y+5s)
          left:  (x+2.5s, y-5s) -> (x-2.5s, y) -> (x+2.5s, y+5s)
        where s = arrow_scale (= thinnest_width * 0.4 in genomeview).

        Arrow vertical extent is scaled to row_height * thin_frac in visual
        coordinates and converted to genomic coordinates via axes aspect
        ratio.  Arrow horizontal extent follows genomeview's 2.5:5 ratio.
        """
        span = region.end - region.start
        intron_width = x_end - x_start
        if intron_width < span * 0.002:
            return

        thin_frac = 0.3  # UTR height fraction — matches intron arrow extent

        # Convert arrow vertical extent from visual to genomic coordinates
        ax_pos = ax.get_position()
        fig = ax.get_figure()
        fig_w, fig_h = fig.get_size_inches()
        axes_w_inches = ax_pos.width * fig_w
        axes_h_inches = ax_pos.height * fig_h
        if axes_h_inches > 0 and axes_w_inches > 0:
            genomic_per_visual = (span / axes_w_inches) / (1.0 / axes_h_inches)
            # Half-height in genomic coords (genomeview: 5*scale)
            half_h_g = row_height * thin_frac * genomic_per_visual
            # Half-width follows genomeview's 2.5:5 ratio
            half_w_g = half_h_g * 0.5
        else:
            half_h_g = row_height * thin_frac * span * 0.001
            half_w_g = half_h_g * 0.5

        # Number of arrows: one per row_height*0.75 of intron length
        n_arrows = max(1, int(round(intron_width / (row_height * 0.75))))
        n_arrows = min(n_arrows, 20)

        for i in range(1, n_arrows + 1):
            frac = i / (n_arrows + 1)
            cx = x_start + (x_end - x_start) * frac

            if strand == "+":
                # genomeview right: (x-2.5s, y-5s) -> (x+2.5s, y) -> (x-2.5s, y+5s)
                verts = [
                    (cx - half_w_g, y_center - half_h_g),
                    (cx + half_w_g, y_center),
                    (cx - half_w_g, y_center + half_h_g),
                ]
            else:
                # genomeview left: (x+2.5s, y-5s) -> (x-2.5s, y) -> (x+2.5s, y+5s)
                verts = [
                    (cx + half_w_g, y_center - half_h_g),
                    (cx - half_w_g, y_center),
                    (cx + half_w_g, y_center + half_h_g),
                ]

            arrow = mpatches.Polygon(
                verts, closed=False,
                facecolor="none", edgecolor=color,
                linewidth=lwd * 0.75, alpha=alpha, zorder=2.5,
            )
            ax.add_patch(arrow)

    def _draw_exon_label(self, ax, row, x_start, x_end, y_center, h, region):
        """Draw a label on an exon based on the exon_annotation setting."""
        span = region.end - region.start
        box_width = x_end - x_start
        if box_width < span * 0.005:
            return  # Too small for labels

        fontsize = self.get_param("fontsize", 8) * 0.75
        fontcolor = self.get_param("fontcolor", "#333333")
        mid_x = (x_start + x_end) / 2

        mode = self.exon_annotation
        label = None

        if mode == "exon":
            label = str(row.get("feature", "exon"))
        elif mode == "symbol":
            for col in ["gene_name", "symbol", "gene"]:
                if col in row.index and pd.notna(row[col]):
                    label = str(row[col])
                    break
        elif mode == "gene":
            for col in ["gene_id", "gene", "gene_name"]:
                if col in row.index and pd.notna(row[col]):
                    label = str(row[col])
                    break
        elif mode == "transcript":
            for col in ["transcript_id", "transcript"]:
                if col in row.index and pd.notna(row[col]):
                    label = str(row[col])
                    break
        elif mode == "feature":
            label = str(row.get("feature", ""))

        if label:
            ax.text(mid_x, y_center + h / 2 + 0.01, label,
                    ha="center", va="bottom", fontsize=fontsize,
                    color=fontcolor, clip_on=True, zorder=5)

    def _draw_label(self, ax, grp, region, y_center, row_height,
                    fontsize, fontcolor):
        """Draw a gene or transcript label (UCSC style: label to the LEFT)."""
        # Handle transcript_annotation alias
        effective_show_id = self.show_id
        if self.transcript_annotation is not None:
            if self.transcript_annotation == "symbol":
                effective_show_id = "symbol"
            else:
                effective_show_id = self.transcript_annotation

        if effective_show_id == "gene":
            for col in ["gene_name", "gene", "gene_id", "symbol"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label) and label != "unknown":
                        break
            else:
                return
        elif effective_show_id == "transcript":
            for col in ["transcript_id", "transcript"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label):
                        break
            else:
                return
        elif effective_show_id == "symbol":
            # Use gene_name instead of gene_id
            for col in ["gene_name", "symbol", "gene"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label):
                        break
            else:
                return
        elif effective_show_id == "exon":
            return  # Exon labels handled by exon_annotation
        else:
            return

        # Override with gene_symbols if requested
        if self.gene_symbols and effective_show_id == "gene":
            for col in ["gene_name", "symbol"]:
                if col in grp.columns:
                    label = grp[col].iloc[0]
                    if pd.notna(label):
                        break

        # UCSC style: position label to the LEFT of the transcript
        # Gene name is drawn to the left with right-alignment (ha="right")
        span = region.end - region.start
        label_offset = span * 0.005

        x_pos = max(grp["start"].min() - label_offset, region.start)
        ha = "right"

        x_pos = np.clip(x_pos, region.start, region.end)
        y_pos = y_center

        ax.text(x_pos, y_pos, str(label),
                ha=ha, va="center", fontsize=fontsize,
                color=fontcolor,
                clip_on=True, zorder=5)


# Register class-specific display parameter overrides
from ._base import _CLASS_DISPLAY_PARAM_OVERRIDES
_CLASS_DISPLAY_PARAM_OVERRIDES["GeneRegionTrack"] = {
    "col": "none",
    "fill": "#E89E9D",
    "fill_rev": "#8C8FCE",
    "fill_utr": None,
    "col_intron": None,
    "lwd": 0.8,
    "fontsize": 8,
    "fontcolor": "#333333",
}
