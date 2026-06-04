"""
GeneRegionTrack - A track for displaying gene models.

Displays gene structures with multiple drawing styles:

**UCSC style** (default):
- Backbone line from transcript start to end
- Coding exons (CDS) as thick solid blocks
- UTRs as thinner blocks (~half CDS height)
- Stepped polygons for CDS/UTR transitions within an exon block
- Introns as thin horizontal connecting lines with directional chevron arrows
- Gene labels positioned to the LEFT of the transcript

**flybase style**:
- Backbone line from transcript start to end
- Last exon drawn as a filled directional arrow polygon
- UTR with configurable height (``height_utr``)

**tssarrow style**:
- Vertical line + arrow at the transcription start site (TSS)
- Half-height exon boxes with L-shaped intron connections
- Configurable arrow length via ``arrow_length`` display param

**exonarrows style**:
- Full-height exon boxes with directional arrows drawn inside
- Filled rectangle intron connectors between exons

Supports transcript grouping, strand-aware coloring, and gene/transcript/exon labels.

Ported from Gviz's GeneRegionTrack-class.R and pyGenomeTracks' BedTrack.py.
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
from matplotlib.lines import Line2D


# Valid gene drawing styles
_VALID_STYLES = {"UCSC", "flybase", "tssarrow", "exonarrows"}


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
    style : str
        Drawing style for gene models. One of:
        - 'UCSC' (default): thick CDS, thin UTR, intron chevron arrows
        - 'flybase': backbone line, last exon as filled arrow
        - 'tssarrow': TSS arrow + vertical line, half-height exons
        - 'exonarrows': full-height exons with arrows inside
    name : str
        Track name. Default is "GeneRegion".
    height : float
        Relative track height. Default is 1.5.
    display_params : dict, optional
        Additional display parameters. Supports pyGenomeTracks-style params:
        ``color_utr``, ``color_backbone``, ``color_arrow``, ``height_utr``,
        ``height_intron``, ``arrow_interval``, ``arrowhead_fraction``,
        ``arrowhead_included``, ``arrow_length``.

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
        style: str = "UCSC",
        name: str = "GeneRegion",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
        **kwargs,
    ):
        if style not in _VALID_STYLES:
            raise ValueError(
                f"Invalid style '{style}'. Must be one of {sorted(_VALID_STYLES)}."
            )

        dp = {
            "col": "none",                   # genomeview: no edge stroke
            "fill": "#E89E9D",               # genomeview: forward strand color
            "fill_rev": "#8C8FCE",           # genomeview: reverse strand color
            "fill_utr": None,                # genomeview: same color, just thinner
            "col_intron": None,              # genomeview: same color as exon
            "lwd": 0.8,
            "fontsize": 8,
            "fontcolor": "#333333",
            # --- pyGenomeTracks-style display params ---
            "color_utr": None,               # None = same as exon; or explicit color
            "color_backbone": None,          # None = same as intron/exon color
            "color_arrow": None,             # None = same as exon color
            "height_utr": 1.0,               # UTR height relative to CDS (1=same)
            "height_intron": 0.5,            # intron connector height relative to CDS
            "arrow_interval": 2,             # distance between arrows (in units of arrow size)
            "arrowhead_fraction": 0.004,     # arrowhead size as fraction of region span
            "arrowhead_included": False,     # whether arrow tip is at interval boundary
            "arrow_length": None,            # tssarrow: explicit arrow length in bp (None=auto)
        }
        if display_params:
            dp.update(display_params)

        super().__init__(
            data=data, stacking=stacking, name=name, height=height,
            display_params=dp, **kwargs,
        )

        self.collapse_transcripts = collapse_transcripts
        self.show_id = show_id
        self.style = style
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

    @staticmethod
    def _split_to_blocks(grp: pd.DataFrame):
        """Split a transcript group into (start, end, type) blocks.

        Mirrors pyGenomeTracks' ``_split_bed_to_blocks`` but works with
        geneview's GFF-style DataFrame.  Returns a list of tuples
        ``(start, end, feature_type)`` where *feature_type* is ``"coding"``
        or ``"UTR"``.

        Supports two data layouts:
        1. **GFF-style** – each row is an exon/CDS/UTR feature.
        2. **BED12-style** – a single row with ``thick_start``, ``thick_end``,
           ``block_count``, ``block_sizes``, ``block_starts`` columns.
        """
        # --- BED12 path ---
        if "block_count" in grp.columns and "block_sizes" in grp.columns:
            first = grp.iloc[0]
            try:
                block_count = int(first["block_count"])
                block_sizes = [
                    int(x) for x in str(first["block_sizes"]).rstrip(",").split(",") if x
                ]
                block_starts = [
                    int(x) for x in str(first["block_starts"]).rstrip(",").split(",") if x
                ]
            except (ValueError, TypeError):
                return []

            tx_start = int(first["start"])
            thick_start = (
                int(first["thick_start"])
                if "thick_start" in first.index and pd.notna(first.get("thick_start"))
                else None
            )
            thick_end = (
                int(first["thick_end"])
                if "thick_end" in first.index and pd.notna(first.get("thick_end"))
                else None
            )

            n_blocks = min(block_count, len(block_sizes), len(block_starts))
            positions = []
            for i in range(n_blocks):
                x0 = tx_start + block_starts[i]
                x1 = x0 + block_sizes[i]

                # No coding region at all
                if thick_start is None or thick_end is None or thick_start == thick_end:
                    positions.append((x0, x1, "UTR"))
                    continue

                # CDS start within this block
                if x0 < thick_start < x1:
                    positions.append((x0, thick_start, "UTR"))
                    x0 = thick_start

                # CDS end within this block
                if x0 < thick_end < x1:
                    positions.append((x0, thick_end, "coding"))
                    x0 = thick_end

                if x1 <= thick_start or x0 >= thick_end:
                    ftype = "UTR"
                else:
                    ftype = "coding"
                positions.append((x0, x1, ftype))

            return positions

        # --- GFF-style path ---
        grp_sorted = grp.sort_values("start")
        features = grp_sorted["feature"].astype(str).str.lower().str.strip() if "feature" in grp_sorted.columns else pd.Series(["exon"] * len(grp_sorted))
        has_cds = features.isin({"cds"}).any()
        has_utr = features.isin({"utr", "utr5", "utr3", "five_prime_utr",
                                  "three_prime_utr", "5utr", "3utr",
                                  "ncrna", "nc_rna"}).any()

        if not has_cds and not has_utr:
            # All features are generic 'exon' – treat as coding
            return [
                (int(r["start"]), int(r["end"]), "coding")
                for _, r in grp_sorted.iterrows()
            ]

        positions = []
        for (_, r), feat in zip(grp_sorted.iterrows(), features):
            s, e = int(r["start"]), int(r["end"])
            if feat in {"cds"}:
                positions.append((s, e, "coding"))
            elif feat in {"utr", "utr5", "utr3", "five_prime_utr",
                          "three_prime_utr", "5utr", "3utr",
                          "ncrna", "nc_rna", "start_codon", "stop_codon"}:
                positions.append((s, e, "UTR"))
            else:
                # generic exon – check overlap with CDS/UTR already added
                positions.append((s, e, "coding"))

        return positions

    def _draw_transcript(self, ax, grp, region, y_center, row_height,
                         fill_color, fill_utr, edge_color, intron_color, lwd, alpha):
        """Draw a complete transcript, dispatching to style-specific methods."""
        grp = self._deduplicate_features(grp).sort_values("start")
        strand = grp["strand"].iloc[0] if "strand" in grp.columns else "*"

        # Strand-based coloring
        if strand == "-" and hasattr(self, '_fill_rev'):
            tx_color = self._fill_rev
        else:
            tx_color = fill_color

        # Resolve effective colours
        effective_intron = self._intron_color_override if hasattr(self, '_intron_color_override') and self._intron_color_override else tx_color
        effective_utr = self._fill_utr_override if hasattr(self, '_fill_utr_override') and self._fill_utr_override else tx_color
        effective_edge = edge_color if edge_color and edge_color != "none" else "none"

        # pyGenomeTracks-style colour overrides
        color_utr = self.get_param("color_utr", None) or effective_utr
        color_backbone = self.get_param("color_backbone", None) or effective_intron
        color_arrow = self.get_param("color_arrow", None) or tx_color
        height_utr = self.get_param("height_utr", 1.0)
        height_intron = self.get_param("height_intron", 0.5)
        arrow_interval = self.get_param("arrow_interval", 2)
        arrowhead_fraction = self.get_param("arrowhead_fraction", 0.004)
        arrowhead_included = self.get_param("arrowhead_included", False)

        span = region.end - region.start
        small_relative = arrowhead_fraction * span

        if self.style == "flybase":
            self._draw_transcript_flybase(
                ax, grp, region, y_center, row_height,
                tx_color, color_utr, effective_edge, color_backbone,
                lwd, alpha, strand, height_utr, small_relative, arrowhead_included,
            )
        elif self.style == "tssarrow":
            self._draw_transcript_tssarrow(
                ax, grp, region, y_center, row_height,
                tx_color, color_utr, effective_edge,
                lwd, alpha, strand, height_utr, small_relative,
            )
        elif self.style == "exonarrows":
            self._draw_transcript_exonarrows(
                ax, grp, region, y_center, row_height,
                tx_color, color_utr, effective_edge, color_arrow,
                lwd, alpha, strand, height_intron, arrow_interval,
                small_relative, region,
            )
        else:
            # UCSC (default)
            self._draw_transcript_ucsc(
                ax, grp, region, y_center, row_height,
                tx_color, effective_utr, effective_edge, effective_intron,
                lwd, alpha, strand,
            )

    def _draw_transcript_ucsc(self, ax, grp, region, y_center, row_height,
                               tx_color, effective_utr, edge_color, intron_color,
                               lwd, alpha, strand):
        """Draw transcript in UCSC Genome Browser style.

        - Backbone line from transcript start to end (pyGenomeTracks)
        - CDS: thick solid blocks (no edge stroke by default)
        - UTR: thinner blocks (~half CDS height)
        - Stepped polygons for CDS/UTR transitions within an exon block
        - Introns: thin horizontal connecting lines with directional chevron arrows
        """
        thick_frac = 0.7   # CDS
        thin_frac = 0.35   # UTR

        positions = self._split_to_blocks(grp)
        if not positions:
            return

        # --- Step 1: Draw backbone line from transcript start to end ---
        tx_start = int(grp["start"].min())
        tx_end = int(grp["end"].max())
        color_backbone = self.get_param("color_backbone", None) or intron_color
        bx_start = max(tx_start, region.start)
        bx_end = min(tx_end, region.end)
        if bx_end > bx_start:
            ax.plot([bx_start, bx_end], [y_center, y_center],
                    color=color_backbone, linewidth=max(0.5, lwd),
                    alpha=alpha, zorder=1)

        # --- Step 2: Draw blocks with stepped polygons for CDS/UTR transitions ---
        last_block_end = None
        for s, e, ft in positions:
            x_s = max(s, region.start)
            x_e = min(e, region.end)
            if x_e <= x_s:
                if last_block_end is None or e > last_block_end:
                    last_block_end = e
                continue

            # Draw backbone line in gap between previous block and this one
            if last_block_end is not None and last_block_end < s:
                gap_s = max(last_block_end, region.start)
                gap_e = min(s, region.end)
                if gap_e > gap_s:
                    ax.plot([gap_s, gap_e], [y_center, y_center],
                            color=color_backbone, linewidth=max(0.5, lwd),
                            alpha=alpha, zorder=1)
                    if strand in ("+", "-"):
                        self._draw_intron_arrows(
                            ax, gap_s, gap_e, y_center,
                            row_height, strand, intron_color, lwd, alpha, region,
                        )

            is_utr = (ft == "UTR")

            if is_utr:
                h = row_height * thin_frac
                rect = mpatches.Rectangle(
                    (x_s, y_center - h / 2), x_e - x_s, h,
                    facecolor=effective_utr, edgecolor=edge_color,
                    linewidth=0, alpha=alpha, zorder=3,
                )
                ax.add_patch(rect)
            else:
                # CDS: check for stepped polygon (CDS/UTR transition within this block)
                next_cds_start = None
                prev_cds_end = None
                for ns, ne, nft in positions:
                    if nft == "UTR" and ns == e:
                        next_cds_start = ns
                    if nft == "UTR" and ne == s:
                        prev_cds_end = ne

                y_thick_top = y_center + row_height * thick_frac / 2
                y_thick_bot = y_center - row_height * thick_frac / 2
                y_thin_top = y_center + row_height * thin_frac / 2
                y_thin_bot = y_center - row_height * thin_frac / 2

                if next_cds_start is not None and prev_cds_end is not None:
                    # CDS in middle with UTR on both sides: stepped polygon
                    vertices = [
                        (x_s, y_thin_top), (x_s, y_thin_bot),
                        (x_s, y_thick_bot), (x_s, y_thick_top),
                        (x_e, y_thick_top), (x_e, y_thick_bot),
                        (x_e, y_thin_bot), (x_e, y_thin_top),
                    ]
                    poly = mpatches.Polygon(
                        vertices, closed=True, fill=True,
                        facecolor=tx_color, edgecolor="none",
                        linewidth=0, alpha=alpha, zorder=3,
                    )
                    ax.add_patch(poly)
                elif next_cds_start is not None:
                    # CDS transitions to UTR at right end
                    vertices = [
                        (x_s, y_thick_top), (x_s, y_thick_bot),
                        (x_e, y_thick_bot), (x_e, y_thin_bot),
                        (x_e, y_thin_top), (x_e, y_thick_top),
                    ]
                    poly = mpatches.Polygon(
                        vertices, closed=True, fill=True,
                        facecolor=tx_color, edgecolor="none",
                        linewidth=0, alpha=alpha, zorder=3,
                    )
                    ax.add_patch(poly)
                elif prev_cds_end is not None:
                    # UTR transitions to CDS at left end
                    vertices = [
                        (x_s, y_thin_top), (x_s, y_thick_top),
                        (x_e, y_thick_top), (x_e, y_thick_bot),
                        (x_e, y_thin_bot), (x_s, y_thin_bot),
                    ]
                    poly = mpatches.Polygon(
                        vertices, closed=True, fill=True,
                        facecolor=tx_color, edgecolor="none",
                        linewidth=0, alpha=alpha, zorder=3,
                    )
                    ax.add_patch(poly)
                else:
                    # Pure CDS: simple thick rectangle
                    rect = mpatches.Rectangle(
                        (x_s, y_center - row_height * thick_frac / 2),
                        x_e - x_s, row_height * thick_frac,
                        facecolor=tx_color, edgecolor=edge_color,
                        linewidth=0, alpha=alpha, zorder=3,
                    )
                    ax.add_patch(rect)

            last_block_end = e

        # --- Step 3: Draw exon annotation labels (if configured) ---
        if self.exon_annotation:
            for _, row in grp.iterrows():
                x_start = max(row["start"], region.start)
                x_end = min(row["end"], region.end)
                if x_end <= x_start:
                    continue
                feat = str(row.get("feature", "exon"))
                h = row_height * (thin_frac if _is_thin_feature(feat) else thick_frac)
                self._draw_exon_label(ax, row, x_start, x_end, y_center, h, region)

    # ------------------------------------------------------------------
    # pyGenomeTracks-style drawing helpers
    # ------------------------------------------------------------------

    def _make_bed_arrow(self, start, end, strand, y_center, half_height,
                        small_relative, arrowhead_included):
        """Return polygon vertices for a filled directional arrow (BED convention).

        Ported from pyGenomeTracks ``BedTrack._draw_arrow``.
        """
        y0 = y_center - half_height
        y1 = y_center + half_height

        if strand == "+":
            x0 = start
            if arrowhead_included:
                x1 = max(start, end - small_relative)
                x2 = end
            else:
                x1 = end
                x2 = end + small_relative
            vertices = [(x0, y0), (x0, y1), (x1, y1),
                        (x2, y_center), (x1, y0)]
        else:
            if arrowhead_included:
                x0 = min(end, start + small_relative)
                xb = start
            else:
                x0 = start
                xb = start - small_relative
            x1 = end
            vertices = [(x0, y0), (xb, y_center), (x0, y1),
                        (x1, y1), (x1, y0)]
        return vertices

    def _plot_small_arrow_on_interval(self, ax, start_pos, end_pos, y_top, y_bottom,
                                       strand, color, lwd, alpha, small_relative,
                                       arrow_interval):
        """Plot evenly distributed small chevron arrows on an interval.

        Ported from pyGenomeTracks ``BedTrack.plot_arrows_on_interval``.
        """
        interval_length = end_pos - start_pos + 1
        if interval_length <= small_relative:
            return
        interval_center = start_pos + interval_length / 2
        step = max(1, int(arrow_interval * small_relative))
        pos = np.arange(start_pos + small_relative,
                        end_pos + small_relative, step)
        if len(pos) == 0:
            return
        pos = pos + interval_center - pos.mean()

        for xpos in pos:
            if strand == ".":
                continue
            if strand == "+":
                xdata = [xpos - small_relative / 4,
                         xpos + small_relative / 4,
                         xpos - small_relative / 4]
            else:
                xdata = [xpos + small_relative / 4,
                         xpos - small_relative / 4,
                         xpos + small_relative / 4]
            ydata = [y_top, (y_bottom + y_top) / 2, y_bottom]
            ax.add_line(Line2D(xdata, ydata, color=color,
                               linewidth=max(0.3, lwd * 0.6), alpha=alpha, zorder=4))

    # ------------------------------------------------------------------
    # flybase style
    # ------------------------------------------------------------------

    def _draw_transcript_flybase(self, ax, grp, region, y_center, row_height,
                                  tx_color, utr_color, edge_color, backbone_color,
                                  lwd, alpha, strand, height_utr,
                                  small_relative, arrowhead_included):
        """Draw transcript in Flybase Gbrowse style.

        - Backbone line from transcript start to end
        - Last exon drawn as filled directional arrow polygon
        - UTR uses configurable height (``height_utr``) and colour
        """
        positions = self._split_to_blocks(grp)
        if not positions:
            return

        tx_start = int(grp["start"].min())
        tx_end = int(grp["end"].max())
        thick_frac = 0.7

        # Simple case: single exon, no introns, fully coding
        if len(positions) == 1:
            s, e, ft = positions[0]
            x_s = max(s, region.start)
            x_e = min(e, region.end)
            if x_e <= x_s:
                return
            h = row_height * thick_frac
            if strand in ("+", "-"):
                verts = self._make_bed_arrow(
                    x_s, x_e, strand, y_center, h / 2,
                    small_relative, arrowhead_included,
                )
                poly = mpatches.Polygon(
                    verts, closed=True, fill=True,
                    facecolor=tx_color, edgecolor=edge_color,
                    linewidth=0, alpha=alpha, zorder=3,
                )
                ax.add_patch(poly)
            else:
                rect = mpatches.Rectangle(
                    (x_s, y_center - h / 2), x_e - x_s, h,
                    facecolor=tx_color, edgecolor=edge_color,
                    linewidth=0, alpha=alpha, zorder=3,
                )
                ax.add_patch(rect)
            return

        # Draw backbone line
        bx_start = max(tx_start, region.start)
        bx_end = min(tx_end, region.end)
        if bx_end > bx_start:
            ax.plot([bx_start, bx_end], [y_center, y_center],
                    color=backbone_color, linewidth=max(0.5, lwd),
                    alpha=alpha, zorder=1)

        # Order positions for arrow drawing: the *last* drawn exon gets the arrow
        draw_positions = list(positions)
        if strand == "-":
            draw_positions = draw_positions[::-1]

        # The first element in draw_positions (after reversal for -) is the
        # "leading" exon that gets the arrow head.
        first_pos = draw_positions.pop()

        # Draw the leading exon as a filled arrow
        s, e, ft = first_pos
        if ft == "UTR":
            _color = utr_color
            hh = row_height * thick_frac * height_utr / 2
        else:
            _color = tx_color
            hh = row_height * thick_frac / 2

        x_s = max(s, region.start)
        x_e = min(e, region.end)
        if x_e > x_s:
            verts = self._make_bed_arrow(
                x_s, x_e, strand, y_center, hh,
                small_relative, arrowhead_included,
            )
            poly = mpatches.Polygon(
                verts, closed=True, fill=True,
                facecolor=_color, edgecolor=edge_color,
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(poly)

        # Draw remaining exons as rectangles
        for s, e, ft in draw_positions:
            if ft == "UTR":
                _color = utr_color
                h = row_height * thick_frac * height_utr
            else:
                _color = tx_color
                h = row_height * thick_frac

            x_s = max(s, region.start)
            x_e = min(e, region.end)
            if x_e <= x_s:
                continue
            rect = mpatches.Rectangle(
                (x_s, y_center - h / 2), x_e - x_s, h,
                facecolor=_color, edgecolor=edge_color,
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(rect)

    # ------------------------------------------------------------------
    # tssarrow style
    # ------------------------------------------------------------------

    def _draw_transcript_tssarrow(self, ax, grp, region, y_center, row_height,
                                   tx_color, utr_color, edge_color,
                                   lwd, alpha, strand, height_utr,
                                   small_relative):
        """Draw transcript in TSS-arrow style.

        - Vertical line + arrow at TSS (start for +, end for -)
        - Half-height exon boxes with L-shaped intron connections
        """
        positions = self._split_to_blocks(grp)
        if not positions:
            return

        thick_frac = 0.7
        half_box_h = row_height * thick_frac * 0.5  # half-height exons

        y_bottom = y_center + half_box_h
        y_top_intron = y_center - row_height * thick_frac * 0.25

        # Draw TSS arrow + vertical line
        if strand in ("+", "-"):
            arrow_length_dp = self.get_param("arrow_length", None)
            if arrow_length_dp is not None:
                arrow_length = float(arrow_length_dp)
            else:
                arrow_length = 10 * small_relative
            y_arrow = y_center - row_height * thick_frac * 0.4
            head_width = row_height * thick_frac * 0.25
            head_length = small_relative * 3

            if strand == "+":
                x_tss = grp["start"].min()
                dx = arrow_length
            else:
                x_tss = grp["end"].max()
                dx = -arrow_length

            y_bottom_line = y_center + half_box_h
            ax.add_line(Line2D(
                (x_tss, x_tss), (y_bottom_line, y_arrow),
                color=tx_color, linewidth=max(0.5, lwd), alpha=alpha, zorder=3,
            ))
            ax.arrow(x_tss, y_arrow, dx, 0, overhang=1, width=0,
                     head_width=head_width, head_length=head_length,
                     length_includes_head=True,
                     color=tx_color, linewidth=max(0.3, lwd * 0.5),
                     alpha=alpha, zorder=3)

        # Draw exon boxes and L-shaped intron connections
        last_corner = None
        for s, e, ft in positions:
            if ft == "UTR":
                _color = utr_color
                h = half_box_h * height_utr
            else:
                _color = tx_color
                h = half_box_h

            x_s = max(s, region.start)
            x_e = min(e, region.end)
            if x_e <= x_s:
                continue

            y0 = y_center + h / 2
            vertices = [(x_s, y0), (x_s, y0 - h),
                        (x_e, y0 - h), (x_e, y0)]
            poly = mpatches.Polygon(
                vertices, closed=True, fill=True,
                facecolor=_color, edgecolor="none",
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(poly)

            # L-shaped intron connection
            if last_corner is not None and last_corner[0] < x_s:
                xdata = (last_corner[0],
                         (last_corner[0] + x_s) / 2,
                         x_s)
                ydata = (last_corner[1], y_top_intron, y_center - h / 2 + h)
                ax.add_line(Line2D(xdata, ydata, color=last_corner[2],
                                   linewidth=max(0.3, lwd * 0.5),
                                   alpha=alpha, zorder=2))

            last_corner = (x_e, y_center - h / 2, _color)

    # ------------------------------------------------------------------
    # exonarrows style
    # ------------------------------------------------------------------

    def _draw_transcript_exonarrows(self, ax, grp, region, y_center, row_height,
                                     tx_color, utr_color, edge_color, arrow_color,
                                     lwd, alpha, strand, height_intron,
                                     arrow_interval, small_relative, plot_region):
        """Draw transcript with arrows inside exon boxes.

        ::
            ___________          ____________
            |          |________|            |
            | >  >  >  |________|  >  >  >   |
            |__________|        |____________|
        """
        positions = self._split_to_blocks(grp)
        if not positions:
            return

        thick_frac = 0.7
        full_h = row_height * thick_frac
        intron_h = full_h * height_intron

        y_top = y_center + full_h / 2
        y_bottom = y_center - full_h / 2
        y_intron_top = y_center + intron_h / 2
        y_intron_bottom = y_center - intron_h / 2

        last_exon_end = None
        for s, e, ft in positions:
            if ft == "UTR":
                _color = utr_color
            else:
                _color = tx_color

            x_s = max(s, region.start)
            x_e = min(e, region.end)
            if x_e <= x_s:
                continue

            # Draw exon rectangle
            vertices = [(x_s, y_top), (x_s, y_bottom),
                        (x_e, y_bottom), (x_e, y_top)]
            poly = mpatches.Polygon(
                vertices, closed=True, fill=True,
                facecolor=_color, edgecolor=edge_color,
                linewidth=0, alpha=alpha, zorder=3,
            )
            ax.add_patch(poly)

            # Draw directional arrows inside the exon
            self._plot_small_arrow_on_interval(
                ax, x_s, x_e,
                y_center + full_h * 0.25, y_center - full_h * 0.25,
                strand, arrow_color, lwd, alpha,
                small_relative, arrow_interval,
            )

            # Draw intron connector (filled rectangle) between exons
            if last_exon_end is not None and last_exon_end < x_s:
                ix_s = max(last_exon_end, region.start)
                ix_e = min(x_s, region.end)
                if ix_e > ix_s:
                    verts = [
                        (ix_s, y_intron_top), (ix_s, y_intron_bottom),
                        (ix_e, y_intron_bottom), (ix_e, y_intron_top),
                    ]
                    intron_poly = mpatches.Polygon(
                        verts, closed=True, fill=True,
                        facecolor=_color, edgecolor=edge_color,
                        linewidth=0, alpha=alpha, zorder=3,
                    )
                    ax.add_patch(intron_poly)

            last_exon_end = x_e

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
    # pyGenomeTracks-style params
    "color_utr": None,
    "color_backbone": None,
    "color_arrow": None,
    "height_utr": 1.0,
    "height_intron": 0.5,
    "arrow_interval": 2,
    "arrowhead_fraction": 0.004,
    "arrowhead_included": False,
    "arrow_length": None,
}
