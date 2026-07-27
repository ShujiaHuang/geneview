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

# Formats used by the ``vector`` one-line toggle in :func:`save_figure`.
_VECTOR_FORMATS = {"pdf", "svg", "eps", "ps"}
_RASTER_FORMATS = {"png", "jpg", "jpeg", "tif", "tiff"}


def save_figure(
    axes,
    filepath: str,
    dpi: Optional[int] = None,
    fmt: Optional[str] = None,
    bbox_inches: Optional[str] = None,
    vector: Optional[bool] = None,
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
    dpi : int, optional
        Resolution for raster formats (PNG, JPEG).  When ``None`` (default),
        the active plot style's ``savefig_dpi`` is used (e.g. 300 for the
        journal styles), falling back to 150 if no style is active.
    fmt : str, optional
        Explicit output format (``"png"``, ``"pdf"``, ``"svg"``, ``"eps"``).
        If ``None``, inferred from *filepath* extension.  When *fmt* differs
        from the path's extension, the on-disk extension is rewritten to
        match (so ``save_figure(axes, "fig.png", fmt="pdf")`` writes
        ``fig.pdf``).  This makes *fmt* a one-line raster/vector switch.
    bbox_inches : str, optional
        Passed to ``Figure.savefig``.  When ``None`` (default), the active
        style's ``savefig_bbox`` is used, falling back to ``"tight"``.
    vector : bool, optional
        One-line vector toggle.  When ``True``, the figure is exported as a
        vector format (keeping an already-vector extension, else defaulting
        to PDF).  When ``False``, a raster format is forced (keeping a
        raster extension, else defaulting to PNG).  When ``None`` (default),
        *fmt* / the path extension decide.  An explicit *fmt* wins over
        *vector*.
    **kwargs
        Additional keyword arguments passed to ``Figure.savefig``.

    Returns
    -------
    str
        The path to the saved file (with the extension rewritten to match
        the resolved format when necessary).

    Examples
    --------
    >>> axes = plot_tracks(tracks, region=region)          # doctest: +SKIP
    >>> save_figure(axes, "output.pdf")                    # doctest: +SKIP
    'output.pdf'
    >>> save_figure(axes, "output.png", fmt="svg")         # doctest: +SKIP
    'output.svg'
    >>> save_figure(axes, "output.png", vector=True)       # doctest: +SKIP
    'output.pdf'
    """
    if not axes:
        raise ValueError("No axes provided.")
    # Honor the active journal style for export defaults so a single
    # ``set_style(...)`` also raises the export DPI / bbox strategy.
    from ..plotstyle import get_active_style
    active = get_active_style()
    if dpi is None:
        dpi = active.savefig_dpi if active is not None else 150
    if bbox_inches is None:
        bbox_inches = active.savefig_bbox if active is not None else "tight"

    # Resolve the target format.  Precedence: explicit ``fmt`` > ``vector``
    # toggle > the path's own extension > PNG.
    root, ext = os.path.splitext(filepath)
    ext_fmt = ext[1:].lower() if ext else None
    if fmt is None and vector is not None:
        if vector:
            fmt = ext_fmt if ext_fmt in _VECTOR_FORMATS else "pdf"
        else:
            fmt = ext_fmt if ext_fmt in _RASTER_FORMATS else "png"
    target_fmt = (fmt or ext_fmt or "png").lower()

    # Keep the on-disk extension consistent with the actual format, but only
    # when the path already carries a (differing) extension -- an explicit
    # extensionless path is left untouched (the caller chose it deliberately).
    if ext_fmt is not None and target_fmt != ext_fmt:
        filepath = root + "." + target_fmt

    fig = axes[0].figure if isinstance(axes, list) else axes.figure
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches,
                format=target_fmt, **kwargs)
    return filepath
