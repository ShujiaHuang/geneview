"""
Quick-consensus mismatch filtering for BAM/CRAM alignments.

Tallies per-base nucleotide and indel frequencies across all reads covering a
genomic region, then exposes a simple ``query()`` method that answers:
"is there sufficient evidence (above a configurable threshold) for this
alternative allele at this position?"

This is essential when visualising long-read data (PacBio / ONT) where
individual reads have high per-base error rates; without consensus filtering
the pileup view is overwhelmed by noise.

Ported from ``genomeview.quickconsensus.MismatchCounts``.
"""

from typing import Optional

import numpy as np

from ._utils import match_chrom_format


# Nucleotide / event type to row index in the counts matrix.
_TYPES_TO_ID = {"A": 0, "C": 1, "G": 2, "T": 3, "DEL": 5}
_N_TYPES = 6  # rows in the counts matrix (A, C, G, T, unused, DEL)


class MismatchCounts:
    """Accumulate per-position base and indel tallies from a BAM file.

    Parameters
    ----------
    chrom : str
        Chromosome name.
    start : int
        Region start (0-based, inclusive).
    end : int
        Region end (exclusive).

    Attributes
    ----------
    counts : np.ndarray
        Shape ``(6, end - start)`` matrix; rows correspond to A, C, G, T,
        (unused), DEL.
    insertions : np.ndarray
        Shape ``(end - start)`` array of insertion event counts.
    """

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)

        length = self.end - self.start
        self.counts: np.ndarray = np.zeros((_N_TYPES, length), dtype=np.float64)
        self.insertions: np.ndarray = np.zeros(length, dtype=np.float64)

    # ------------------------------------------------------------------
    # Tallying
    # ------------------------------------------------------------------

    def tally_reads(self, bam) -> None:
        """Walk every read overlapping the region and accumulate tallies.

        Parameters
        ----------
        bam : pysam.AlignmentFile
            An *open* alignment file.  The chromosome name is normalised
            automatically via :func:`match_chrom_format`.
        """
        chrom = match_chrom_format(self.chrom, bam.references)
        for pileup_col in bam.pileup(chrom, self.start, self.end, truncate=True):
            for pileup_read in pileup_col.pileups:
                if pileup_read.is_refskip:
                    continue
                if pileup_read.is_del:
                    self._add_count(pileup_col.pos, "DEL")
                else:
                    nuc = pileup_read.alignment.query_sequence[pileup_read.query_position]
                    if nuc != "N":
                        self._add_count(pileup_col.pos, nuc)
                if pileup_read.indel > 0:
                    self._add_count(pileup_col.pos, "INS")

    def _add_count(self, position: int, type_: str) -> None:
        idx = position - self.start
        if idx < 0 or idx >= self.end - self.start:
            return
        if type_ == "INS":
            self.insertions[idx] += 1
        else:
            row = _TYPES_TO_ID.get(type_)
            if row is not None:
                self.counts[row, idx] += 1

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def query(
        self,
        type_: str,
        start: int,
        end: Optional[int] = None,
        threshold: float = 0.2,
        del_threshold: float = 0.3,
    ) -> bool:
        """Return True if *type_* has sufficient support in [start, end].

        Parameters
        ----------
        type_ : str
            One of ``"A"``, ``"C"``, ``"G"``, ``"T"``, ``"DEL"``, ``"INS"``.
        start : int
            Genomic start of the query window.
        end : int, optional
            Genomic end of the query window (defaults to *start*).
        threshold : float
            Minimum allele fraction for SNVs and insertions.  Default 0.2.
        del_threshold : float
            Minimum allele fraction for deletions.  Default 0.3.

        Returns
        -------
        bool
        """
        if start < self.start or start >= self.end:
            return False

        s = start - self.start
        if end is None:
            e = s
        else:
            e = min(end - self.start, self.end - self.start - 1)

        total = self.counts[:, s : e + 1].sum(axis=0)
        total = total.astype(float)

        if type_ == "INS":
            ins = self.insertions[s : e + 1]
            with np.errstate(divide="ignore", invalid="ignore"):
                frac = np.where(total > 0, ins / total, 0.0)
            return bool((frac > threshold).any())

        row = _TYPES_TO_ID.get(type_)
        if row is None:
            return False

        this_type = self.counts[row, s : e + 1]
        thr = del_threshold if type_ == "DEL" else threshold
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total > 0, this_type / total, 0.0)
        return bool((frac > thr).any())
