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
    return _read_alignment_coverage(filepath, region=region, bins=bins, mode="rb")


def read_cram_coverage(
    filepath: str,
    region: Optional[GenomicInterval] = None,
    bins: int = 1000,
    reference: Optional[str] = None,
) -> pd.DataFrame:
    """Compute coverage from a CRAM file.

    Requires the optional ``pysam`` package. CRAM files are reference-based
    compressed alignments; a reference FASTA is typically needed to decode
    them. If the reference path is embedded in the CRAM header and
    accessible on the local filesystem, ``reference`` can be omitted.

    Parameters
    ----------
    filepath : str
        Path to the CRAM file (must be indexed with ``samtools index``).
    region : GenomicInterval, optional
        Genomic region to compute coverage for.
    bins : int
        Number of bins for summarizing coverage.
    reference : str, optional
        Path to the reference FASTA file (must be indexed with ``samtools faidx``).
        Required unless the reference is already accessible via the CRAM header.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value (coverage depth).

    Examples
    --------
    >>> from geneview.genometracks import read_cram_coverage, GenomicInterval
    >>> region = GenomicInterval("chr7", 2000000, 2050000)
    >>> df = read_cram_coverage("sample.cram", region=region,
    ...                         reference="hg38.fa")  # doctest: +SKIP
    """
    return _read_alignment_coverage(
        filepath, region=region, bins=bins, mode="rc", reference=reference,
    )


def _read_alignment_coverage(
    filepath: str,
    region: Optional[GenomicInterval] = None,
    bins: int = 1000,
    mode: str = "rb",
    reference: Optional[str] = None,
) -> pd.DataFrame:
    """Shared implementation for BAM and CRAM coverage computation.

    Parameters
    ----------
    filepath : str
        Path to the BAM/CRAM file.
    region : GenomicInterval, optional
        Genomic region to compute coverage for.
    bins : int
        Number of bins for summarizing coverage.
    mode : str
        pysam open mode: ``"rb"`` for BAM, ``"rc"`` for CRAM.
    reference : str, optional
        Path to the reference FASTA (required for CRAM).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value (coverage depth).
    """
    try:
        import pysam
    except ImportError:
        fmt = "CRAM" if mode == "rc" else "BAM"
        raise ImportError(
            f"The 'pysam' package is required to read {fmt} files. "
            "Install it with: pip install pysam"
        )

    open_kwargs = {}
    if reference is not None:
        open_kwargs["reference_filename"] = reference

    alignment = pysam.AlignmentFile(filepath, mode, **open_kwargs)
    try:
        if region is not None:
            chrom = region.chrom
            start = region.start
            end = region.end
        else:
            # Use first reference
            chrom = alignment.references[0]
            start = 0
            end = alignment.get_reference_length(chrom)

        # Use pysam count to get coverage
        bin_size = max(1, (end - start) // bins)
        records = []
        for pos in range(start, end, bin_size):
            bin_end = min(pos + bin_size, end)
            count = alignment.count(chrom, pos, bin_end)
            coverage = count / max(1, bin_end - pos)
            records.append({
                "chrom": chrom,
                "start": pos,
                "end": bin_end,
                "value": coverage,
            })
        return pd.DataFrame(records)
    finally:
        alignment.close()


def read_wig(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Read a WIG (Wiggle) file into a DataFrame.

    Supports both ``fixedStep`` and ``variableStep`` WIG formats.

    Parameters
    ----------
    filepath : str
        Path to the WIG file.
    nrows : int, optional
        Maximum number of data rows to read.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, start, end, value.
    """
    records = []
    chrom = None
    start = 0
    step = 1
    span = 1
    mode = None  # "fixedStep" or "variableStep"
    count = 0

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            if line.startswith("fixedStep"):
                mode = "fixedStep"
                params = _parse_wig_header(line)
                chrom = params.get("chrom", chrom)
                start = int(params.get("start", 1)) - 1  # WIG is 1-based
                step = int(params.get("step", 1))
                span = int(params.get("span", 1))
                continue

            if line.startswith("variableStep"):
                mode = "variableStep"
                params = _parse_wig_header(line)
                chrom = params.get("chrom", chrom)
                span = int(params.get("span", 1))
                continue

            if mode == "fixedStep":
                try:
                    value = float(line)
                except ValueError:
                    continue
                records.append({
                    "chrom": chrom,
                    "start": start,
                    "end": start + span,
                    "value": value,
                })
                start += step
                count += 1
                if nrows and count >= nrows:
                    break

            elif mode == "variableStep":
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        pos = int(parts[0]) - 1  # WIG is 1-based
                        value = float(parts[1])
                    except ValueError:
                        continue
                    records.append({
                        "chrom": chrom,
                        "start": pos,
                        "end": pos + span,
                        "value": value,
                    })
                    count += 1
                    if nrows and count >= nrows:
                        break

    if not records:
        return pd.DataFrame(columns=["chrom", "start", "end", "value"])

    df = pd.DataFrame(records)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def _parse_wig_header(line: str) -> dict:
    """Parse a WIG header line (fixedStep or variableStep) into a dict."""
    params = {}
    # Remove the keyword prefix
    parts = line.split(None, 1)
    if len(parts) < 2:
        return params
    for item in parts[1].split():
        if "=" in item:
            key, value = item.split("=", 1)
            params[key.strip()] = value.strip()
    return params


def read_fasta(
    filepath: str,
    region: Optional[GenomicInterval] = None,
) -> str:
    """Read a genomic sequence from a FASTA file.

    Requires the optional ``pysam`` package for indexed FASTA access.

    Parameters
    ----------
    filepath : str
        Path to the FASTA file (must be indexed with ``samtools faidx``).
    region : GenomicInterval, optional
        Genomic region to extract.  If *None*, the entire first sequence
        is returned.

    Returns
    -------
    str
        Nucleotide sequence string (uppercase ACGTN).
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "The 'pysam' package is required to read FASTA files. "
            "Install it with: pip install pysam"
        )

    fa = pysam.FastaFile(filepath)
    try:
        if region is not None:
            seq = fa.fetch(region.chrom, region.start, region.end)
        else:
            seq = fa.fetch(fa.references[0])
        return seq.upper()
    finally:
        fa.close()


def read_2bit(
    filepath: str,
    region: Optional[GenomicInterval] = None,
) -> str:
    """Read a genomic sequence from a 2bit file.

    Requires the optional ``py2bit`` package.

    Parameters
    ----------
    filepath : str
        Path to the 2bit file.
    region : GenomicInterval, optional
        Genomic region to extract.  If *None*, the entire first sequence
        is returned.

    Returns
    -------
    str
        Nucleotide sequence string (uppercase ACGTN).
    """
    try:
        import py2bit
    except ImportError:
        raise ImportError(
            "The 'py2bit' package is required to read 2bit files. "
            "Install it with: pip install py2bit"
        )

    tb = py2bit.open(filepath)
    try:
        if region is not None:
            seq = tb.sequence(region.chrom, region.start, region.end)
        else:
            chroms = tb.chroms()
            chrom = next(iter(chroms))
            seq = tb.sequence(chrom, 0, chroms[chrom])
        return seq.upper()
    finally:
        tb.close()


def read_auto(filepath: str, nrows: Optional[int] = None) -> pd.DataFrame:
    """Auto-detect file format and read genomic data.

    Detects format based on file extension:
        - .bed -> read_bed
        - .gff, .gff3, .gtf -> read_gff
        - .bedgraph, .bdg -> read_bedgraph
        - .wig -> read_wig
        - .bw, .bigwig -> read_bigwig
        - .bam -> read_bam_coverage
        - .cram -> read_cram_coverage

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
        ".wig": read_wig,
    }

    if ext in readers:
        return readers[ext](filepath, nrows=nrows)
    elif ext in (".bw", ".bigwig"):
        return read_bigwig(filepath)
    elif ext == ".bam":
        return read_bam_coverage(filepath)
    elif ext == ".cram":
        return read_cram_coverage(filepath)
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
            f"Supported extensions: .bed, .gff, .gff3, .gtf, .bedgraph, .bdg, "
            f".wig, .bw, .bigwig, .bam, .cram"
        )
