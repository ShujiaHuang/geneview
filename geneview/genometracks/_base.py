"""
Core base classes and data structures for genome track visualization.

This module provides the foundational classes for the genome track system,
inspired by Gviz (R/Bioconductor) but adapted for Python/matplotlib.

Classes:
    GenomicInterval: A single genomic region (chrom, start, end, strand).
    Track: Abstract base class for all track types.
    RangeTrack: Base class for tracks containing genomic ranges.
    StackedTrack: Base class for tracks with stacking of overlapping features.
    NumericTrack: Base class for tracks with numeric data.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Union

import numpy as np
import pandas as pd


@dataclass
class GenomicInterval:
    """Represents a genomic region with chromosome, start, end, and optional strand.

    This is the Python equivalent of R's GRanges for a single interval.

    Parameters
    ----------
    chrom : str
        Chromosome identifier (e.g. "chr1", "1").
    start : int
        Start position (0-based, inclusive).
    end : int
        End position (exclusive).
    strand : str, optional
        Strand information: "+", "-", or "*" (unstranded). Default is "*".
    """
    chrom: str
    start: int
    end: int
    strand: str = "*"

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)
        if self.end < self.start:
            raise ValueError(
                f"End position ({self.end}) cannot be less than start ({self.start})."
            )
        if self.strand not in ("+", "-", "*"):
            raise ValueError(f"Strand must be '+', '-', or '*', got '{self.strand}'.")

    @property
    def width(self) -> int:
        """Width of the interval in base pairs."""
        return self.end - self.start

    def overlaps(self, other: "GenomicInterval") -> bool:
        """Check if this interval overlaps with another on the same chromosome."""
        return (self.chrom == other.chrom and
                self.start < other.end and
                self.end > other.start)

    def contains(self, position: int, chrom: Optional[str] = None) -> bool:
        """Check if a position is within this interval."""
        if chrom is not None and self.chrom != chrom:
            return False
        return self.start <= position < self.end

    def extend(self, left: int = 0, right: int = 0) -> "GenomicInterval":
        """Return a new interval extended by the given amounts."""
        return GenomicInterval(
            chrom=self.chrom,
            start=max(0, self.start - left),
            end=self.end + right,
            strand=self.strand,
        )

    def __repr__(self):
        strand_str = f", strand={self.strand}" if self.strand != "*" else ""
        return f"GenomicInterval({self.chrom}:{self.start}-{self.end}{strand_str})"


# Default display parameter values mirroring Gviz's GdObject prototype
_DEFAULT_DISPLAY_PARAMS: Dict[str, Any] = {
    "alpha": 1.0,
    "background_panel": "white",
    "background_title": "#D3D3D3",      # Gviz: lightgray
    "col": "#0080FF",                    # Gviz: DEFAULT_SYMBOL_COL
    "fill": "lightgray",                 # Gviz: DEFAULT_FILL_COL
    "col_axis": "#A9A9A9",               # Gviz: darkgray
    "col_border": "transparent",          # Gviz: transparent by default
    "col_grid": "#808080",               # Gviz: medium gray
    "col_title": "white",                # Gviz: white text on gray
    "col_border_title": "white",          # Gviz: white border on title
    "fontface": "normal",
    "fontface_title": "bold",
    "fontsize": 10,
    "fontsize_title": 10,
    "frame": False,
    "grid": False,
    "lwd": 1.0,
    "lty": "-",
    "show_title": True,
    "reverse_strand": False,
    "rotation_title": 90,                # Gviz: rotated 90 degrees
    "cex": 1.0,
    "min_width": 1,     # minimum feature width in pixels
    "min_height": 3,    # minimum feature height in pixels
    "min_distance": 1,  # minimum distance between features for collapsing
    "collapse": True,
    "stack_height": 0.75,  # Gviz: use 75% of available row height
}


class Track(ABC):
    """Abstract base class for all genome track types.

    Every track has a name, a relative height, display parameters controlling
    appearance, and a ``draw`` method that renders onto a matplotlib Axes.

    Parameters
    ----------
    name : str
        Track name displayed in the title panel.
    height : float
        Relative height of the track (default 1.0).
    display_params : dict, optional
        Dictionary of display parameters overriding defaults.
    """

    def __init__(
        self,
        name: str = "Track",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        self.name = name
        self.height = height
        self._dp = dict(_DEFAULT_DISPLAY_PARAMS)
        if display_params:
            self._dp.update(display_params)

    def get_param(self, key: str, default: Any = None) -> Any:
        """Get a display parameter value."""
        return self._dp.get(key, default)

    def set_param(self, key: str, value: Any) -> None:
        """Set a display parameter value."""
        self._dp[key] = value

    def set_params(self, params: Dict[str, Any]) -> None:
        """Set multiple display parameters at once."""
        self._dp.update(params)

    @abstractmethod
    def draw(self, ax, region: GenomicInterval) -> None:
        """Draw the track content on the given matplotlib Axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        region : GenomicInterval
            The genomic region currently displayed.
        """
        pass

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the full genomic extent of this track, or None if not applicable."""
        return None

    def _format_axis(self, ax, region: GenomicInterval, show_xlabel: bool = True):
        """Common axis formatting for genomic tracks."""
        ax.set_xlim(region.start, region.end)
        ax.get_yaxis().set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        if not show_xlabel:
            ax.set_xticklabels([])
        else:
            ax.xaxis.set_major_formatter(
                _genomic_position_formatter(region.end - region.start)
            )

    def __repr__(self):
        return f"{self.__class__.__name__}(name='{self.name}', height={self.height})"


class RangeTrack(Track):
    """Base class for tracks that contain genomic ranges.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with at minimum 'chrom', 'start', 'end' columns.
    name : str
        Track name.
    height : float
        Relative track height.
    display_params : dict, optional
        Display parameters.
    """

    REQUIRED_COLUMNS = {"chrom", "start", "end"}

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        name: str = "RangeTrack",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(name=name, height=height, display_params=display_params)
        self._data = self._validate_data(data)

    def _validate_data(self, data) -> Optional[pd.DataFrame]:
        """Validate and normalize input data."""
        if data is None:
            return None
        if isinstance(data, str):
            # File path - delegate to I/O utilities
            from ._io import read_auto
            data = read_auto(data)
        if not isinstance(data, pd.DataFrame):
            raise TypeError(f"Expected DataFrame or file path, got {type(data)}")
        # Copy to avoid mutating the caller's DataFrame
        data = data.copy()
        # Normalize column names to lowercase
        data.columns = [c.lower().strip() for c in data.columns]
        missing = self.REQUIRED_COLUMNS - set(data.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        data["start"] = data["start"].astype(int)
        data["end"] = data["end"].astype(int)
        data["chrom"] = data["chrom"].astype(str)
        return data

    @property
    def data(self) -> Optional[pd.DataFrame]:
        """The underlying DataFrame of genomic ranges."""
        return self._data

    def subset(self, region: GenomicInterval) -> Optional[pd.DataFrame]:
        """Return data rows overlapping the given region."""
        if self._data is None or len(self._data) == 0:
            return None
        mask = (
            (self._data["chrom"] == region.chrom) &
            (self._data["start"] < region.end) &
            (self._data["end"] > region.start)
        )
        result = self._data[mask].copy()
        return result if len(result) > 0 else None

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the full genomic extent of the data."""
        if self._data is None or len(self._data) == 0:
            return None
        chrom = self._data["chrom"].iloc[0]
        return GenomicInterval(
            chrom=chrom,
            start=int(self._data["start"].min()),
            end=int(self._data["end"].max()),
        )


class StackedTrack(RangeTrack):
    """Base class for tracks with stacking of overlapping features.

    Stacking determines how overlapping features are arranged vertically.
    Supported modes: 'squish', 'pack', 'dense', 'full', 'hide'.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with genomic range data.
    stacking : str
        Stacking mode. Default is 'squish'.
    stack_height : float
        Fraction of available height per stack row. Default is 0.75.
    name : str
        Track name.
    height : float
        Relative track height.
    display_params : dict, optional
        Display parameters.
    """

    STACKING_MODES = ("hide", "dense", "squish", "pack", "full")

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        stacking: str = "squish",
        stack_height: float = 0.75,
        name: str = "StackedTrack",
        height: float = 1.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if stacking not in self.STACKING_MODES:
            raise ValueError(
                f"Invalid stacking mode '{stacking}'. "
                f"Must be one of {self.STACKING_MODES}."
            )
        super().__init__(data=data, name=name, height=height, display_params=display_params)
        self.stacking = stacking
        self.stack_height = stack_height
        self._stacks: Optional[np.ndarray] = None

    def compute_stacks(self, data: pd.DataFrame) -> np.ndarray:
        """Compute stack assignments for overlapping features.

        Returns an array of stack row indices (0-based) for each feature.
        """
        from ._stacking import compute_stacking
        return compute_stacking(
            data["start"].values,
            data["end"].values,
            mode=self.stacking,
            min_distance=int(self.get_param("min_distance", 1)),
        )


class NumericTrack(RangeTrack):
    """Base class for tracks containing numeric data values.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame with 'chrom', 'start', 'end' plus numeric value column(s).
    value_columns : list of str, optional
        Column names containing numeric values. If None, auto-detected.
    name : str
        Track name.
    height : float
        Relative track height.
    display_params : dict, optional
        Display parameters.
    """

    def __init__(
        self,
        data: Union[pd.DataFrame, str, None] = None,
        value_columns: Optional[List[str]] = None,
        name: str = "NumericTrack",
        height: float = 1.5,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(data=data, name=name, height=height, display_params=display_params)
        self._value_columns = value_columns
        if self._data is not None and self._value_columns is None:
            self._value_columns = self._detect_value_columns()

    def _detect_value_columns(self) -> List[str]:
        """Auto-detect numeric value columns."""
        if self._data is None:
            return []
        non_range_cols = [
            c for c in self._data.columns
            if c not in ("chrom", "start", "end", "strand", "name",
                         "id", "group", "feature", "score", "thick_start",
                         "thick_end", "item_rgb", "block_count", "block_sizes",
                         "block_starts", "gene_id", "gene_name", "transcript_id",
                         "exon_id")
        ]
        return [c for c in non_range_cols if pd.api.types.is_numeric_dtype(self._data[c])]

    @property
    def value_columns(self) -> List[str]:
        """Column names containing numeric values."""
        return self._value_columns or []

    def get_ylim(self) -> tuple:
        """Compute y-axis limits from data."""
        if self._data is None or not self.value_columns:
            return (0, 1)
        vals = self._data[self.value_columns].values.ravel()
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            return (0, 1)
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        margin = (vmax - vmin) * 0.05
        return (vmin - margin, vmax + margin)


def _genomic_position_formatter(span: int):
    """Return a tick formatter appropriate for the genomic span."""
    from matplotlib.ticker import FuncFormatter

    def _fmt(x, _pos):
        if span >= 1e9:
            return f"{x / 1e9:.1f}Gb"
        elif span >= 1e6:
            return f"{x / 1e6:.1f}Mb"
        elif span >= 1e3:
            return f"{x / 1e3:.1f}kb"
        else:
            return f"{int(x)}bp"
    return FuncFormatter(_fmt)


def _format_genomic_position(value: float, span: int) -> str:
    """Format a single genomic position value."""
    if span >= 1e9:
        return f"{value / 1e9:.1f}Gb"
    elif span >= 1e6:
        return f"{value / 1e6:.1f}Mb"
    elif span >= 1e3:
        return f"{value / 1e3:.1f}kb"
    else:
        return f"{int(value)}bp"
