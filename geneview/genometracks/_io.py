"""
File I/O utilities for genome track data.

Provides readers for common genomic file formats:
    - BED (Browser Extensible Data)
    - GFF/GTF (General Feature Format)
    - BigWig (binary wiggle format, optional pyBigWig dependency)
    - BAM coverage (optional pysam dependency)

All readers return pandas DataFrames with standardized column names.
"""

import os
from typing import Optional, Union

import pandas as pd
import numpy as np

from ._base import GenomicInterval


# Standard BED column names
_BED_COLUMNS = [
    "chrom", "start", "end", "name", "score", "strand",
    "thick_start", "thick_end", "item_rgb", "block_count",
    "block_sizes", "block_starts",
]

# Standard GFF columns
_GFF_COLUMNS = [
    "chrom", "source", "feature", "start", "end",
    "score", "strand", "frame", "attributes",
]


def read_bed(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Read a BED file into a DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the BED file.
    nrows : int, optional
        Maximum number of rows to read.

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized BED column names.
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
    )
    # Assign standard column names based on number of columns
    ncols = df.shape[1]
    df.columns = _BED_COLUMNS[:ncols]

    # Ensure proper types
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["chrom"] = df["chrom"].astype(str)

    if "strand" in df.columns:
        df["strand"] = df["strand"].fillna("*")

    return df


def read_gff(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Read a GFF/GTF file into a DataFrame.

    Supports both GFF3 and GTF formats. Parses the attributes column
    to extract common fields (gene_id, transcript_id, gene_name, etc.).

    Parameters
    ----------
    filepath : str
        Path to the GFF/GTF file.
    nrows : int, optional
        Maximum number of rows to read.

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized GFF column names plus parsed attributes.
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
    )
    df.columns = _GFF_COLUMNS

    # GFF is 1-based, convert to 0-based half-open
    df["start"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)
    df["chrom"] = df["chrom"].astype(str)
    df["strand"] = df["strand"].fillna("*")

    # Parse attributes column
    parsed = _parse_gff_attributes(df["attributes"])
    df = pd.concat([df.drop(columns=["attributes"]), parsed], axis=1)

    return df


def _parse_gff_attributes(attr_series: pd.Series) -> pd.DataFrame:
    """Parse GFF/GTF attributes column into separate DataFrame columns.

    Handles both GFF3 (key=value;) and GTF (key "value";) formats.
    """
    records = []
    for attrs in attr_series:
        record = {}
        if pd.isna(attrs) or attrs == ".":
            records.append(record)
            continue

        attrs_str = str(attrs)
        # Detect format: GFF3 uses '=', GTF uses quotes
        if "=" in attrs_str and '"' not in attrs_str:
            # GFF3 format: key=value;key=value
            for item in attrs_str.split(";"):
                item = item.strip()
                if "=" in item:
                    key, value = item.split("=", 1)
                    record[key.strip()] = value.strip()
        else:
            # GTF format: key "value"; key "value"
            for item in attrs_str.split(";"):
                item = item.strip()
                if " " in item:
                    parts = item.split(None, 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        value = parts[1].strip().strip('"').strip("'")
                        record[key] = value
        records.append(record)

    return pd.DataFrame(records)


def read_bedgraph(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Read a bedGraph file into a DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the bedGraph file.
    nrows : int, optional
        Maximum number of rows to read.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value.
    """
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        comment="#",
        nrows=nrows,
    )
    if df.shape[1] >= 4:
        df = df.iloc[:, :4]
        df.columns = ["chrom", "start", "end", "value"]
    else:
        raise ValueError(
            f"bedGraph file must have at least 4 columns, got {df.shape[1]}."
        )

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["chrom"] = df["chrom"].astype(str)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")

    return df


def read_bigwig(
    filepath: str,
    region: Optional[GenomicInterval] = None,
    bins: int = 1000,
) -> pd.DataFrame:
    """Read a BigWig file into a DataFrame.

    Requires the optional ``pyBigWig`` package.

    Parameters
    ----------
    filepath : str
        Path to the BigWig file.
    region : GenomicInterval, optional
        Genomic region to read. If None, reads the entire first chromosome.
    bins : int
        Number of bins to use when summarizing the signal.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value.
    """
    try:
        import pyBigWig
    except ImportError:
        raise ImportError(
            "The 'pyBigWig' package is required to read BigWig files. "
            "Install it with: pip install pyBigWig"
        )

    bw = pyBigWig.open(filepath)
    try:
        chroms = bw.chroms()
        if region is not None:
            chrom = region.chrom
            start = region.start
            end = region.end
        else:
            chrom = list(chroms.keys())[0]
            start = 0
            end = chroms[chrom]

        values = bw.stats(chrom, start, end, nBins=bins, type="mean")
        bin_size = (end - start) / bins

        records = []
        for i, val in enumerate(values):
            bin_start = int(start + i * bin_size)
            bin_end = int(start + (i + 1) * bin_size)
            records.append({
                "chrom": chrom,
                "start": bin_start,
                "end": bin_end,
                "value": val if val is not None else np.nan,
            })
        return pd.DataFrame(records)
    finally:
        bw.close()


def read_bam_coverage(
    filepath: str,
    region: Optional[GenomicInterval] = None,
    bins: int = 1000,
) -> pd.DataFrame:
    """Compute coverage from a BAM file.

    Requires the optional ``pysam`` package.

    Parameters
    ----------
    filepath : str
        Path to the BAM file (must be indexed).
    region : GenomicInterval, optional
        Genomic region to compute coverage for.
    bins : int
        Number of bins for summarizing coverage.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value (coverage depth).
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "The 'pysam' package is required to read BAM files. "
            "Install it with: pip install pysam"
        )

    bam = pysam.AlignmentFile(filepath, "rb")
    try:
        if region is not None:
            chrom = region.chrom
            start = region.start
            end = region.end
        else:
            # Use first reference
            chrom = bam.references[0]
            start = 0
            end = bam.get_reference_length(chrom)

        # Use pysam count to get coverage
        bin_size = max(1, (end - start) // bins)
        records = []
        for pos in range(start, end, bin_size):
            bin_end = min(pos + bin_size, end)
            count = bam.count(chrom, pos, bin_end)
            coverage = count / max(1, bin_end - pos)
            records.append({
                "chrom": chrom,
                "start": pos,
                "end": bin_end,
                "value": coverage,
            })
        return pd.DataFrame(records)
    finally:
        bam.close()


def read_auto(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Auto-detect file format and read genomic data.

    Detects format based on file extension:
        - .bed -> read_bed
        - .gff, .gff3, .gtf -> read_gff
        - .bedgraph, .bdg -> read_bedgraph
        - .bw, .bigwig -> read_bigwig
        - .bam -> read_bam_coverage

    Parameters
    ----------
    filepath : str
        Path to the genomic data file.
    nrows : int, optional
        Maximum number of rows (not applicable for binary formats).

    Returns
    -------
    pd.DataFrame
        DataFrame with standardized column names.

    Raises
    ------
    ValueError
        If the file format cannot be determined.
    """
    ext = os.path.splitext(filepath)[1].lower()

    readers = {
        ".bed": read_bed,
        ".gff": read_gff,
        ".gff3": read_gff,
        ".gtf": read_gff,
        ".bedgraph": read_bedgraph,
        ".bdg": read_bedgraph,
    }

    if ext in readers:
        return readers[ext](filepath, nrows=nrows)
    elif ext in (".bw", ".bigwig"):
        return read_bigwig(filepath)
    elif ext == ".bam":
        return read_bam_coverage(filepath)
    else:
        # Try as tab-separated with header
        try:
            df = pd.read_csv(filepath, sep="\t", nrows=nrows)
            df.columns = [c.lower().strip() for c in df.columns]
            if {"chrom", "start", "end"}.issubset(set(df.columns)):
                return df
        except Exception:
            pass

        raise ValueError(
            f"Cannot determine file format for '{filepath}'. "
            f"Supported extensions: .bed, .gff, .gff3, .gtf, .bedgraph, .bdg, .bw, .bigwig, .bam"
        )
