"""
IdeogramTrack - A track for displaying chromosome ideograms.

Displays a schematic representation of a chromosome with cytoband coloring,
centromere position, and a highlight box indicating the current genomic region.

Ported from Gviz's IdeogramTrack-class.R.
"""

from typing import Optional, List, Dict, Any, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.path import Path

from ._base import Track, GenomicInterval
from ..utils._dataset import load_dataset


# Cytoband colors (Gviz: biovizBase CYTOBAND palette / circos)
_CYTOBAND_COLORS = {
    "gneg": "#FFFFFF",
    "gpos25": "#C8C8C8",
    "gpos33": "#D2D2D2",
    "gpos50": "#C8C8C8",
    "gpos66": "#A0A0A0",
    "gpos75": "#828282",
    "gpos100": "#000000",
    "gpos": "#000000",
    "acen": "#D92F27",
    "stalk": "#647FA4",
    "gvar": "#DCDCDC",
}

# Default karyotype dataset names in geneview-data repo
_KARYOTYPE_DATASETS = {
    "hg38": "karyotype_human_hg38.txt",
    "hg19": "karyotype_human_hg19.txt",
}


def _get_band_color(gie_stain: str) -> str:
    """Get the color for a cytoband gieStain type."""
    return _CYTOBAND_COLORS.get(str(gie_stain).strip(), "#C8C8C8")


class IdeogramTrack(Track):
    """A track displaying a chromosome ideogram.

    Shows a schematic representation of a chromosome with cytoband coloring,
    centromere position (as a triangle or circle), rounded caps at the ends,
    and a highlight box indicating the current genomic region.

    Mimics Gviz's IdeogramTrack.

    Parameters
    ----------
    bands : pd.DataFrame or str, optional
        Cytoband data with columns: chrom, chromStart, chromEnd, name, gieStain.
        If a string, interpreted as a file path (tab-separated).
        If None (default), automatically loads the human karyotype from the
        geneview-data repository based on ``genome_build``.
    chromosome : str, optional
        Chromosome to display. If None, uses the first chromosome in the data.
    genome_build : str
        Genome build for default cytoband data: ``"hg38"`` or ``"hg19"``.
        Only used when ``bands=None``. Default is ``"hg38"``.
    show_id : bool
        Whether to show the chromosome name label. Default is True.
    show_band_id : bool
        Whether to show band name labels inside bands. Default is False.
    centromere_shape : str
        Shape of the centromere: 'triangle' or 'circle'. Default is 'triangle'.
    outline : bool
        Whether to draw black outline around bands. Default is False.
    name : str
        Track name for the title panel. Default is chromosome name.
    height : float
        Relative track height. Default is 0.5.
    display_params : dict, optional
        Additional display parameters.

    Examples
    --------
    Auto-load human karyotype (hg38) — no data needed:

    >>> from geneview.genometracks import IdeogramTrack, plot_tracks, GenomicInterval
    >>> track = IdeogramTrack(chromosome="chr7")
    >>> plot_tracks([track], region=GenomicInterval("chr7", 20000000, 60000000))

    Use hg19 build:

    >>> track = IdeogramTrack(chromosome="chr7", genome_build="hg19")

    Or provide custom cytoband data:

    >>> import pandas as pd
    >>> bands = pd.DataFrame({
    ...     "chrom": ["chr7"] * 5,
    ...     "chromStart": [0, 20000000, 40000000, 58000000, 60000000],
    ...     "chromEnd": [20000000, 40000000, 58000000, 60000000, 159000000],
    ...     "name": ["p15", "p14", "p13", "p12", "p11"],
    ...     "gieStain": ["gneg", "gpos25", "gneg", "gpos50", "acen"],
    ... })
    >>> track = IdeogramTrack(bands, chromosome="chr7")
    """

    def __init__(
        self,
        bands: Union[pd.DataFrame, str, None] = None,
        chromosome: Optional[str] = None,
        genome_build: str = "hg38",
        show_id: bool = True,
        show_band_id: bool = False,
        centromere_shape: str = "triangle",
        outline: bool = False,
        name: Optional[str] = None,
        height: float = 0.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        dp = {
            "background_title": "transparent",
            "col_border_title": "transparent",
            "show_title": False,
            "col": "red",            # Gviz: red highlight box border
            "fill": "#FFE3E6",       # Gviz: pink highlight box fill
            "lwd": 1.0,
            "fontsize": 10,
            "fontcolor": "#808080",   # Gviz: DEFAULT_SHADED_COL
        }
        if display_params:
            dp.update(display_params)

        super().__init__(name=name or "Ideogram", height=height, display_params=dp)

        # Auto-load from geneview-data when bands is not provided
        if bands is None:
            dataset_name = _KARYOTYPE_DATASETS.get(genome_build)
            if dataset_name is None:
                raise ValueError(
                    f"Unknown genome_build '{genome_build}'. "
                    f"Supported: {list(_KARYOTYPE_DATASETS.keys())}"
                )
            bands = load_dataset(dataset_name)

        self._band_table = self._load_bands(bands)
        self.chromosome = chromosome
        self.genome_build = genome_build
        self.show_id = show_id
        self.show_band_id = show_band_id
        self.centromere_shape = centromere_shape
        self.outline = outline

        # Filter bands to the selected chromosome
        if self._band_table is not None and chromosome is not None:
            chrom_clean = str(chromosome).replace("chr", "")
            self._chr_bands = self._band_table[
                self._band_table["chrom"].astype(str).str.replace("chr", "") == chrom_clean
            ].copy()
            if len(self._chr_bands) == 0:
                # Try exact match
                self._chr_bands = self._band_table[
                    self._band_table["chrom"] == chromosome
                ].copy()
        elif self._band_table is not None:
            # Use first chromosome
            first_chrom = self._band_table["chrom"].iloc[0]
            self.chromosome = first_chrom
            self._chr_bands = self._band_table[
                self._band_table["chrom"] == first_chrom
            ].copy()
        else:
            self._chr_bands = pd.DataFrame()

        # Update name if not explicitly set
        if name is None and self.chromosome is not None:
            self.name = str(self.chromosome)

    def _load_bands(self, bands) -> Optional[pd.DataFrame]:
        """Load and validate cytoband data."""
        if bands is None:
            return None

        if isinstance(bands, str):
            # Read tab-separated file.
            # UCSC/Gviz karyotype files use '#chrom' as the first header line;
            # plain TSV files use 'chrom'.  Peek at the first character to decide.
            with open(bands) as fh:
                first_char = fh.read(1)
            if first_char == "#":
                bands = pd.read_table(
                    bands, comment="#", header=None,
                    names=["chrom", "chromStart", "chromEnd", "name", "gieStain"],
                )
            else:
                bands = pd.read_table(bands)

        if not isinstance(bands, pd.DataFrame):
            raise TypeError(f"Expected DataFrame or file path, got {type(bands)}")

        # Normalize column names
        col_map = {}
        for c in bands.columns:
            cl = c.lower().strip()
            if cl in ("chrom", "chromosome", "chr"):
                col_map[c] = "chrom"
            elif cl in ("chromstart", "start"):
                col_map[c] = "chromStart"
            elif cl in ("chromend", "end"):
                col_map[c] = "chromEnd"
            elif cl in ("name", "band", "band_name"):
                col_map[c] = "name"
            elif cl in ("giestain", "gie_stain", "stain", "type"):
                col_map[c] = "gieStain"
        bands = bands.rename(columns=col_map)

        required = {"chrom", "chromStart", "chromEnd", "name", "gieStain"}
        missing = required - set(bands.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        bands["chromStart"] = bands["chromStart"].astype(int)
        bands["chromEnd"] = bands["chromEnd"].astype(int)
        return bands

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the full extent of the chromosome."""
        if self._chr_bands is None or len(self._chr_bands) == 0:
            return None
        return GenomicInterval(
            chrom=str(self.chromosome),
            start=int(self._chr_bands["chromStart"].min()),
            end=int(self._chr_bands["chromEnd"].max()),
        )

    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the ideogram track.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region to display.
        """
        if self._chr_bands is None or len(self._chr_bands) == 0:
            ax.set_xlim(region.start, region.end)
            ax.set_ylim(0, 1)
            ax.axis("off")
            return

        chr_len = int(self._chr_bands["chromEnd"].max())
        margin_y = 0.25  # Vertical margin

        ax.set_xlim(0, chr_len)
        ax.set_ylim(0, 1)

        # Draw highlight box for current region (behind the chromosome)
        hl_fill = self.get_param("fill", "#FFE3E6")
        hl_col = self.get_param("col", "red")
        hl_lwd = self.get_param("lwd", 1.0)
        region_start = max(0, region.start)
        region_end = min(chr_len, region.end)
        hl_rect = mpatches.FancyBboxPatch(
            (region_start, 0.1), region_end - region_start, 0.8,
            boxstyle="round,pad=0",
            facecolor=hl_fill, edgecolor="none",
            zorder=1,
        )
        ax.add_patch(hl_rect)

        # Draw chromosome bands
        bands = self._chr_bands.sort_values("chromStart")
        lcol = "black"
        lwd = 0.5

        # Find centromere bands
        cent_bands = bands[bands["gieStain"] == "acen"]
        has_centromere = len(cent_bands) >= 2

        # Separate bands before and after centromere
        if has_centromere:
            cent_start = cent_bands["chromStart"].min()
            cent_end = cent_bands["chromEnd"].max()
            bef_bands = bands[bands["chromEnd"] <= cent_start]
            aft_bands = bands[bands["chromStart"] >= cent_end]
        else:
            bef_bands = bands
            aft_bands = pd.DataFrame()

        # Draw normal bands (before centromere)
        for _, row in bef_bands.iterrows():
            self._draw_band(ax, row, margin_y, lcol if self.outline else None, lwd)

        # Draw normal bands (after centromere)
        for _, row in aft_bands.iterrows():
            self._draw_band(ax, row, margin_y, lcol if self.outline else None, lwd)

        # Draw centromere
        if has_centromere and self.centromere_shape == "triangle":
            self._draw_centromere_triangle(ax, cent_bands, chr_len, margin_y)
        elif has_centromere and self.centromere_shape == "circle":
            self._draw_centromere_circle(ax, cent_bands, chr_len, margin_y)
        elif has_centromere:
            # Fallback: draw as regular bands
            for _, row in cent_bands.iterrows():
                self._draw_band(ax, row, margin_y, lcol, lwd)

        # Draw rounded caps
        self._draw_cap(ax, 0, margin_y, "left", bef_bands.iloc[0] if len(bef_bands) > 0 else None)
        last_band = aft_bands.iloc[-1] if len(aft_bands) > 0 else (
            bef_bands.iloc[-1] if len(bef_bands) > 0 else None
        )
        self._draw_cap(ax, chr_len, margin_y, "right", last_band)

        # Draw chromosome outline
        self._draw_outline(ax, chr_len, margin_y, has_centromere, cent_bands, lcol, lwd)

        # Draw highlight box border (on top)
        hl_border_rect = mpatches.FancyBboxPatch(
            (region_start, 0.1), region_end - region_start, 0.8,
            boxstyle="round,pad=0",
            facecolor="none", edgecolor=hl_col,
            linewidth=hl_lwd, zorder=5,
        )
        ax.add_patch(hl_border_rect)

        # Show chromosome name
        if self.show_id:
            chr_label = str(self.chromosome).replace("chr", "")
            chr_label = f"Chromosome {chr_label}"
            fontsize = self.get_param("fontsize", 10)
            fontcolor = self.get_param("fontcolor", "#808080")
            ax.text(chr_len * 0.5, 1.05, chr_label,
                    ha="center", va="bottom", fontsize=fontsize,
                    color=fontcolor, zorder=5)

        # Show band names if requested
        if self.show_band_id:
            self._draw_band_labels(ax, bands, margin_y)

        ax.set_xlim(0, chr_len)
        ax.set_ylim(-0.05, 1.15)
        ax.axis("off")

    def _draw_band(self, ax, row, margin_y, edge_col, lwd):
        """Draw a single cytoband rectangle."""
        color = _get_band_color(row["gieStain"])
        s = row["chromStart"]
        e = row["chromEnd"]
        rect = mpatches.Rectangle(
            (s, margin_y), e - s, 1 - 2 * margin_y,
            facecolor=color,
            edgecolor=edge_col if edge_col else color,
            linewidth=lwd if edge_col else 0,
            zorder=3,
        )
        ax.add_patch(rect)

    def _draw_centromere_triangle(self, ax, cent_bands, chr_len, margin_y):
        """Draw centromere as two triangles forming a pinched waist."""
        acen_color = _CYTOBAND_COLORS["acen"]

        cent_start = cent_bands["chromStart"].min()
        cent_end = cent_bands["chromEnd"].max()
        mid = (cent_start + cent_end) / 2

        # Left triangle (p-arm side)
        verts_left = [
            (cent_start, margin_y),
            (cent_start + (mid - cent_start), 0.5),
            (cent_start, 1 - margin_y),
            (cent_start, margin_y),
        ]
        tri_left = Polygon(verts_left, closed=True,
                           facecolor=acen_color, edgecolor=acen_color,
                           linewidth=0.5, zorder=4)
        ax.add_patch(tri_left)

        # Right triangle (q-arm side)
        verts_right = [
            (mid, 0.5),
            (cent_end, margin_y),
            (cent_end, 1 - margin_y),
            (mid, 0.5),
        ]
        tri_right = Polygon(verts_right, closed=True,
                            facecolor=acen_color, edgecolor=acen_color,
                            linewidth=0.5, zorder=4)
        ax.add_patch(tri_right)

    def _draw_centromere_circle(self, ax, cent_bands, chr_len, margin_y):
        """Draw centromere as a circle."""
        acen_color = _CYTOBAND_COLORS["acen"]
        cent_start = cent_bands["chromStart"].min()
        cent_end = cent_bands["chromEnd"].max()
        mid = (cent_start + cent_end) / 2
        radius = min((cent_end - cent_start) / 2, (1 - 2 * margin_y) / 2)

        circle = mpatches.Circle(
            (mid, 0.5), radius,
            facecolor=acen_color, edgecolor="black",
            linewidth=0.5, zorder=4,
        )
        ax.add_patch(circle)

    def _draw_cap(self, ax, x, margin_y, side, nearest_band):
        """Draw a rounded cap at the chromosome end."""
        color = _get_band_color(nearest_band["gieStain"]) if nearest_band is not None else "#FFFFFF"
        cap_width = 0.01  # Small cap width as fraction of chromosome

        if side == "left":
            # Left rounded cap
            theta = np.linspace(np.pi / 2, 3 * np.pi / 2, 50)
            r = (1 - 2 * margin_y) / 2
            cx = x
            cy = 0.5
            xs = cx + r * 0.02 * np.cos(theta)
            ys = cy + r * np.sin(theta)
        else:
            # Right rounded cap
            theta = np.linspace(-np.pi / 2, np.pi / 2, 50)
            r = (1 - 2 * margin_y) / 2
            cx = x
            cy = 0.5
            xs = cx + r * 0.02 * np.cos(theta)
            ys = cy + r * np.sin(theta)

        verts = list(zip(xs, ys))
        verts.append(verts[0])
        cap = Polygon(verts, closed=True,
                      facecolor=color, edgecolor="black",
                      linewidth=0.5, zorder=4)
        ax.add_patch(cap)

    def _draw_outline(self, ax, chr_len, margin_y, has_centromere, cent_bands, color, lwd):
        """Draw the chromosome outline."""
        # Top and bottom lines
        if has_centromere:
            cent_start = cent_bands["chromStart"].min()
            cent_end = cent_bands["chromEnd"].max()
            mid = (cent_start + cent_end) / 2

            # Top outline
            ax.plot([0, cent_start], [1 - margin_y, 1 - margin_y],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([cent_start, mid], [1 - margin_y, 0.5],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([mid, cent_end], [0.5, 1 - margin_y],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([cent_end, chr_len], [1 - margin_y, 1 - margin_y],
                    color=color, linewidth=lwd, zorder=5)

            # Bottom outline
            ax.plot([0, cent_start], [margin_y, margin_y],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([cent_start, mid], [margin_y, 0.5],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([mid, cent_end], [0.5, margin_y],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([cent_end, chr_len], [margin_y, margin_y],
                    color=color, linewidth=lwd, zorder=5)
        else:
            ax.plot([0, chr_len], [1 - margin_y, 1 - margin_y],
                    color=color, linewidth=lwd, zorder=5)
            ax.plot([0, chr_len], [margin_y, margin_y],
                    color=color, linewidth=lwd, zorder=5)

        # End caps (vertical lines)
        ax.plot([0, 0], [margin_y, 1 - margin_y],
                color=color, linewidth=lwd, zorder=5)
        ax.plot([chr_len, chr_len], [margin_y, 1 - margin_y],
                color=color, linewidth=lwd, zorder=5)

    def _draw_band_labels(self, ax, bands, margin_y):
        """Draw band name labels inside bands."""
        chr_len = int(self._chr_bands["chromEnd"].max())
        fontsize = self.get_param("fontsize", 7) * 0.7

        for _, row in bands.iterrows():
            s = row["chromStart"]
            e = row["chromEnd"]
            mid_x = (s + e) / 2
            band_width = e - s
            label = str(row["name"])

            # Skip if band is too narrow
            if band_width < chr_len * 0.01:
                continue

            # Choose text color based on band darkness
            gie = str(row["gieStain"]).strip()
            text_color = "white" if gie in ("gpos100", "gpos75", "gpos") else "black"

            ax.text(mid_x, 0.5, label,
                    ha="center", va="center",
                    fontsize=fontsize, color=text_color,
                    zorder=6)
