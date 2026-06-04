"""
Shared utilities for genometracks modules.

Provides chromosome-name normalisation and BAM-file introspection helpers
ported from the ``genomeview`` library.
"""

from typing import Iterable, Optional


def match_chrom_format(chrom: str, keys: Iterable[str]) -> str:
    """Return *chrom* or its ``chr``-prefixed/stripped variant that exists in *keys*.

    Many genomic tools use different chromosome naming conventions (``chr1`` vs
    ``1``).  This helper tries the original name first, then the alternate
    form, and falls back to the original if neither matches.

    Parameters
    ----------
    chrom : str
        Chromosome name to look up.
    keys : iterable of str
        Available chromosome names (e.g. from a BAM header or FASTA index).

    Returns
    -------
    str
        The variant of *chrom* present in *keys*, or *chrom* unchanged if
        neither form is found.

    Examples
    --------
    >>> match_chrom_format("chr1", ["1", "2", "3"])
    '1'
    >>> match_chrom_format("1", ["chr1", "chr2"])
    'chr1'
    >>> match_chrom_format("chrM", ["MT"])
    'chrM'
    """
    keys_set = set(keys)
    if chrom in keys_set:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
    else:
        alt = f"chr{chrom}"
    if alt in keys_set:
        return alt
    # Also handle chrM <-> MT
    if chrom == "chrM" and "MT" in keys_set:
        return "MT"
    if chrom == "MT" and "chrM" in keys_set:
        return "chrM"
    return chrom


def is_paired_end(bam_path: str, n: int = 100) -> bool:
    """Detect whether a BAM/CRAM file contains paired-end reads.

    Scans up to *n* reads and returns ``True`` as soon as a paired read is
    encountered.

    Parameters
    ----------
    bam_path : str
        Path to the BAM or CRAM file.
    n : int
        Maximum number of reads to inspect.  Default is 100.

    Returns
    -------
    bool

    Raises
    ------
    ImportError
        If ``pysam`` is not installed.
    """
    import pysam

    aln = pysam.AlignmentFile(bam_path, "rb")
    try:
        for i, read in enumerate(aln.fetch()):
            if read.is_paired:
                return True
            if i >= n:
                break
        return False
    finally:
        aln.close()


def is_long_frag_dataset(bam_path: str, n: int = 1000) -> bool:
    """Detect whether a BAM file comes from a long-fragment sequencing protocol.

    Returns ``True`` when the first *n* reads appear to be single-end with a
    query length exceeding 1000 bp (characteristic of PacBio or ONT data).

    Parameters
    ----------
    bam_path : str
        Path to the BAM or CRAM file.
    n : int
        Maximum number of reads to inspect.  Default is 1000.

    Returns
    -------
    bool

    Raises
    ------
    ImportError
        If ``pysam`` is not installed.
    """
    import pysam

    aln = pysam.AlignmentFile(bam_path, "rb")
    try:
        for i, read in enumerate(aln.fetch()):
            if read.is_paired:
                return False
            if read.query_length > 1000:
                return True
            if i >= n:
                break
        return False
    finally:
        aln.close()


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

_COMPLEMENT_TABLE = str.maketrans("ATCGNatcgn", "TAGCNtagcn")


def reverse_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Non-ACGT characters are left unchanged.

    Parameters
    ----------
    seq : str
        A DNA sequence string.

    Returns
    -------
    str
        The reverse complement.

    Examples
    --------
    >>> reverse_comp("ATCG")
    'CGAT'
    >>> reverse_comp("AATTCC")
    'GGAATT'
    """
    return seq[::-1].translate(_COMPLEMENT_TABLE)
