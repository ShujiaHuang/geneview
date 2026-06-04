"""
Export track data to standard genomic file formats.

Supports BED, GFF, bedGraph, and WIG export from RangeTrack and
NumericTrack subclasses.

Ported from Gviz's exportTracks.R.
"""

import os
from typing import Optional, Union

import pandas as pd

from ._base import RangeTrack, NumericTrack, Track


def export_tracks(
    track: Union[Track, "pd.DataFrame"],
    filepath: str,
    fmt: str = "bed",
    source: str = "geneview",
) -> str:
    """Export track data to a standard genomic file format.

    Parameters
    ----------
    track : Track or pd.DataFrame
        A track object or DataFrame to export.  Must contain the appropriate
        columns for the target format (e.g. chrom, start, end for BED).
    filepath : str
        Output file path.
    fmt : str
        Output format: ``"bed"``, ``"gff"``, ``"bedgraph"``, or ``"wig"``.
    source : str
        Source name used in the GFF ``source`` column.  Default is
        ``"geneview"``.

    Returns
    -------
    str
        The path to the written file.

    Examples
    --------
    >>> import pandas as pd
    >>> from geneview.genometracks import AnnotationTrack, export_tracks
    >>> data = pd.DataFrame({
    ...     "chrom": ["chr1", "chr1"],
    ...     "start": [100, 500],
    ...     "end":   [200, 600],
    ...     "name":  ["feat1", "feat2"],
    ...     "strand": ["+", "-"],
    ... })
    >>> track = AnnotationTrack(data)
    >>> export_tracks(track, "/tmp/test.bed", fmt="bed")  # doctest: +SKIP
    '/tmp/test.bed'
    """
    fmt = fmt.lower().strip()
    valid = ("bed", "gff", "bedgraph", "wig")
    if fmt not in valid:
        raise ValueError(f"Format must be one of {valid}, got '{fmt}'.")

    # Extract the DataFrame from the track
    if isinstance(track, pd.DataFrame):
        df = track.copy()
    elif isinstance(track, (RangeTrack, NumericTrack)):
        if track.data is None:
            raise ValueError("Track has no data to export.")
        df = track.data.copy()
    elif hasattr(track, "data") and isinstance(track.data, pd.DataFrame):
        df = track.data.copy()
    else:
        raise TypeError(
            f"Expected a Track with a .data DataFrame or a pd.DataFrame, "
            f"got {type(track).__name__}."
        )

    # Normalize column names
    df.columns = [c.lower().strip() for c in df.columns]

    if fmt == "bed":
        _write_bed(df, filepath)
    elif fmt == "gff":
        _write_gff(df, filepath, source=source)
    elif fmt == "bedgraph":
        _write_bedgraph(df, filepath)
    elif fmt == "wig":
        _write_wig(df, filepath)

    return filepath


def _write_bed(df: pd.DataFrame, filepath: str) -> None:
    """Write a BED file (BED3 through BED6 depending on available columns)."""
    cols = ["chrom", "start", "end"]
    if "name" in df.columns:
        cols.append("name")
    if "score" in df.columns:
        cols.append("score")
    if "strand" in df.columns:
        cols.append("strand")

    out = df[cols].copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)

    # Fill missing optional columns with BED defaults
    if "name" not in out.columns:
        out["name"] = "."
    if "score" not in out.columns:
        out["score"] = "0"
    if "strand" not in out.columns:
        out["strand"] = "."

    out.to_csv(filepath, sep="\t", header=False, index=False)


def _write_gff(df: pd.DataFrame, filepath: str, source: str = "geneview") -> None:
    """Write a GFF3 file."""
    # GFF is 1-based
    rows = []
    for _, row in df.iterrows():
        seqname = str(row.get("chrom", "."))
        src = source
        feature = str(row.get("feature", "."))
        start = int(row["start"]) + 1  # Convert 0-based to 1-based
        end = int(row["end"])
        score = str(row.get("score", "."))
        strand = str(row.get("strand", "."))
        frame = str(row.get("frame", "."))

        # Build attributes
        attrs = []
        for key in ("name", "id", "gene_id", "transcript_id", "gene_name"):
            if key in row.index and pd.notna(row[key]):
                attrs.append(f"{key}={row[key]}")
        attr_str = ";".join(attrs) if attrs else "."

        rows.append([seqname, src, feature, start, end, score, strand, frame, attr_str])

    out = pd.DataFrame(rows, columns=[
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes",
    ])

    with open(filepath, "w") as fh:
        fh.write("##gff-version 3\n")
    out.to_csv(filepath, sep="\t", header=False, index=False, mode="a")


def _write_bedgraph(df: pd.DataFrame, filepath: str) -> None:
    """Write a bedGraph file."""
    if "value" not in df.columns:
        # Try to find a numeric value column
        for col in df.columns:
            if col not in ("chrom", "start", "end", "strand", "name", "id"):
                if pd.api.types.is_numeric_dtype(df[col]):
                    df = df.rename(columns={col: "value"})
                    break
        if "value" not in df.columns:
            raise ValueError("No numeric value column found for bedGraph export.")

    out = df[["chrom", "start", "end", "value"]].copy()
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out.to_csv(filepath, sep="\t", header=False, index=False)


def _write_wig(df: pd.DataFrame, filepath: str) -> None:
    """Write a WIG file (variableStep format)."""
    if "value" not in df.columns:
        for col in df.columns:
            if col not in ("chrom", "start", "end", "strand", "name", "id"):
                if pd.api.types.is_numeric_dtype(df[col]):
                    df = df.rename(columns={col: "value"})
                    break
        if "value" not in df.columns:
            raise ValueError("No numeric value column found for WIG export.")

    with open(filepath, "w") as fh:
        for chrom, grp in df.groupby("chrom"):
            # Determine span (use the most common span)
            spans = grp["end"].values - grp["start"].values
            span = int(pd.Series(spans).mode().iloc[0]) if len(spans) > 0 else 1

            fh.write(f"variableStep chrom={chrom} span={span}\n")
            for _, row in grp.iterrows():
                pos = int(row["start"]) + 1  # WIG is 1-based
                fh.write(f"{pos}\t{row['value']}\n")


# ---------------------------------------------------------------------------
# Figure export helpers
# ---------------------------------------------------------------------------

def save_figure(
    axes,
    filepath: str,
    dpi: int = 150,
    fmt: Optional[str] = None,
    bbox_inches: str = "tight",
    **kwargs,
) -> str:
    """Save a track figure to disk with auto-detected format.

    Parameters
    ----------
    axes : list of matplotlib.axes.Axes
        The axes returned by :func:`plot_tracks` (or any track plot).
    filepath : str
        Output file path.  The format is inferred from the extension
        (``.png``, ``.pdf``, ``.svg``, ``.eps``).  Override with *fmt*.
    dpi : int
        Resolution for raster formats (PNG, JPEG).  Default is 150.
    fmt : str, optional
        Explicit output format (``"png"``, ``"pdf"``, ``"svg"``, ``"eps"``).
        If ``None``, inferred from *filepath* extension.
    bbox_inches : str
        Passed to ``Figure.savefig``.  Default is ``"tight"``.
    **kwargs
        Additional keyword arguments passed to ``Figure.savefig``.

    Returns
    -------
    str
        The path to the saved file.

    Examples
    --------
    >>> axes = plot_tracks(tracks, region=region)          # doctest: +SKIP
    >>> save_figure(axes, "output.pdf")                    # doctest: +SKIP
    'output.pdf'
    """
    if not axes:
        raise ValueError("No axes provided.")
    fig = axes[0].figure if isinstance(axes, list) else axes.figure
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches, format=fmt, **kwargs)
    return filepath
