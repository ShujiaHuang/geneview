"""Readers for mitochondrial (mtDNA) analysis inputs.

The functions here normalise the three common MitoQuest / mtDNA-pipeline
outputs into tidy :class:`pandas.DataFrame` objects that the plotting layer
consumes directly:

* :func:`read_mito_vcf`        — single- or multi-sample VCF (e.g. from
  ``mitoquest caller``) -> one row per *sample x site x ALT allele*, carrying
  the heteroplasmy fraction (VAF), depth, genotype and HET/HOM status.
* :func:`read_mito_copynumber` — one or more ``mitoquest copynum`` TSV files
  -> one row per sample with the mtDNA copy number and its 95% CI.
* :func:`read_mito_coverage`   — per-base / binned depth across the
  mitochondrial contig from BAM/CRAM files (reuses geneview's alignment
  reader, so BAM and CRAM, single or multi sample, are all handled).

MitoQuest VCF field notes
-------------------------
``mitoquest caller`` writes a per-sample ``FORMAT`` block whose central
quantity is the *heteroplasmy fraction* (the fraction of mtDNA molecules
carrying the ALT allele).  Across MitoQuest versions this has been exposed as
``AF`` (current) or ``HF`` (older builds), always one value per ALT allele.
:func:`read_mito_vcf` therefore probes a list of candidate field names and
uses whichever is present, so both layouts parse without configuration.

Author: Shujia Huang
"""
import os
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._reference import is_mt_contig

# Candidate per-sample FORMAT keys that carry the heteroplasmy fraction / VAF,
# in priority order.  MitoQuest current builds use ``AF``; older ones ``HF``.
_VAF_FORMAT_KEYS = ("AF", "HF", "VAF", "FREQ")
# Candidate per-sample depth keys, in priority order.
_DP_FORMAT_KEYS = ("DP", "AD_SUM", "NR")

# Standard column order for the tidy long-format variant table.
MITO_VCF_COLUMNS = [
    "sample", "chrom", "pos", "ref", "alt",
    "vaf", "depth", "gt", "status", "var_type", "variant_id",
]


def _classify_var_type(ref: str, alt: str) -> str:
    """Classify an ALT allele as SNV / INS / DEL relative to *ref*."""
    if ref is None or alt is None:
        return "UNK"
    if len(ref) == len(alt):
        return "SNV" if len(ref) == 1 else "MNV"
    return "DEL" if len(ref) > len(alt) else "INS"


def _first_present(mapping, keys):
    """Return the value of the first key in *keys* present (and not None)."""
    for key in keys:
        try:
            val = mapping.get(key)
        except (KeyError, TypeError):
            val = None
        if val is not None:
            return key, val
    return None, None


def _as_float(value, index=None):
    """Coerce a (possibly per-allele tuple) FORMAT value to a float."""
    if value is None:
        return np.nan
    if isinstance(value, (tuple, list)):
        if not value:
            return np.nan
        picked = value[index] if (index is not None and index < len(value)) else value[0]
        return _as_float(picked)
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def read_mito_vcf(
    filepath: str,
    samples: Optional[Sequence[str]] = None,
    region: Optional[str] = None,
    min_vaf: float = 0.0,
    include_ref: bool = False,
    mt_only: bool = True,
) -> pd.DataFrame:
    """Read an mtDNA VCF into a tidy long-format variant table.

    Parses a single- or multi-sample VCF (such as the output of
    ``mitoquest caller``) with :mod:`pysam`, expanding every non-reference
    genotype into one row per *sample x site x ALT allele* and attaching the
    heteroplasmy fraction (VAF), depth, genotype string and HET/HOM status.

    Parameters
    ----------
    filepath : str
        Path to the VCF/BCF file.  Plain ``.vcf``, bgzipped ``.vcf.gz`` and
        ``.bcf`` are all accepted.
    samples : sequence of str, optional
        Restrict the output to these sample names.  ``None`` keeps all
        samples present in the VCF header.
    region : str, optional
        Restrict parsing to a ``chrom`` / ``chrom:start-end`` region.  Requires
        the file to be indexed (``.tbi``/``.csi``).  ``None`` scans the whole
        file.
    min_vaf : float, optional
        Drop rows whose VAF is below this threshold (default 0.0 keeps all).
    include_ref : bool, optional
        When True, also emit rows for reference (``0``) alleles.  Default
        False (only variant alleles are returned).
    mt_only : bool, optional
        When True (default), keep only records on a mitochondrial contig
        (``chrM``/``MT``/...); set False to keep every contig.

    Returns
    -------
    pandas.DataFrame
        Columns: ``sample, chrom, pos, ref, alt, vaf, depth, gt, status,
        var_type, variant_id``.  ``pos`` is 1-based; ``status`` is ``"HET"``
        (heteroplasmic: more than one distinct non-missing allele) or
        ``"HOM"`` (homoplasmic).  Empty (0-row) frame when nothing passes the
        filters.

    Raises
    ------
    ImportError
        If :mod:`pysam` is not installed.
    FileNotFoundError
        If *filepath* does not exist.
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "The 'pysam' package is required to read VCF files. "
            "Install it with: pip install pysam"
        )
    if not os.path.exists(filepath):
        raise FileNotFoundError("VCF file not found: %s" % filepath)

    keep = set(samples) if samples is not None else None
    records = []

    with pysam.VariantFile(str(filepath)) as vcf:
        vcf_samples = list(vcf.header.samples)
        if keep is not None:
            vcf_samples = [s for s in vcf_samples if s in keep]

        iterator = vcf.fetch(region=region) if region else vcf
        for rec in iterator:
            chrom = rec.chrom
            if mt_only and not is_mt_contig(chrom):
                continue
            pos = rec.pos
            ref = rec.ref if rec.ref is not None else "."
            alts = list(rec.alts) if rec.alts else []
            rsid = rec.id if rec.id is not None else "."

            for sample_id in vcf_samples:
                call = rec.samples[sample_id]
                gt = call.get("GT")
                gt_alleles = [a for a in gt if a is not None] if gt else []
                gt_str = "/".join("." if a is None else str(a) for a in gt) if gt else "."
                # Heteroplasmy is encoded as more than one distinct allele.
                status = "HET" if len(set(gt_alleles)) > 1 else "HOM"

                vaf_key, vaf_raw = _first_present(call, _VAF_FORMAT_KEYS)
                dp_key, dp_raw = _first_present(call, _DP_FORMAT_KEYS)
                depth = _as_float(dp_raw)

                # Which ALT allele indices are actually carried by this sample.
                alt_indices = sorted({a for a in gt_alleles if a and a > 0})
                if not alt_indices and include_ref:
                    alt_indices = [0]
                for ai in alt_indices:
                    if ai == 0:
                        alt_seq, vtype = ref, "REF"
                        vaf = _as_float(vaf_raw, 0)
                    else:
                        alt_seq = alts[ai - 1] if ai - 1 < len(alts) else "."
                        vtype = _classify_var_type(ref, alt_seq)
                        # VAF field is per-ALT (index ai-1); AD-style per-allele
                        # (index ai) is also tolerated by falling back.
                        vaf = _as_float(vaf_raw, ai - 1)
                    if not include_ref and ai == 0:
                        continue
                    if np.isfinite(vaf) and vaf < min_vaf:
                        continue
                    records.append({
                        "sample": sample_id,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt_seq,
                        "vaf": vaf,
                        "depth": depth,
                        "gt": gt_str,
                        "status": status,
                        "var_type": vtype,
                        "variant_id": rsid,
                    })

    if not records:
        return pd.DataFrame(columns=MITO_VCF_COLUMNS)
    return pd.DataFrame(records, columns=MITO_VCF_COLUMNS)


def read_mito_copynumber(
    filepaths: Union[str, Sequence[str]],
    sample_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """Read one or more ``mitoquest copynum`` TSV files into a tidy table.

    The ``mitoquest copynum`` output has one row per contig; this reader keeps
    the mitochondrial row from each file and returns a per-sample table of the
    mtDNA copy number with its 95% confidence interval.

    Parameters
    ----------
    filepaths : str or sequence of str
        Path(s) to ``copynum`` TSV files.  A single string is treated as one
        file; a list/tuple is treated as one file per sample.
    sample_names : sequence of str, optional
        Sample labels aligned with *filepaths*.  When omitted, each label is
        derived from the file name (leading component before the first dot).

    Returns
    -------
    pandas.DataFrame
        Columns: ``sample, chrom, copy_number, ci_low, ci_high``.  Sorted by
        ``copy_number`` descending.  Empty frame when no mitochondrial row is
        found.
    """
    if isinstance(filepaths, str):
        filepaths = [filepaths]
    filepaths = list(filepaths)

    if sample_names is not None:
        if len(sample_names) != len(filepaths):
            raise ValueError(
                "sample_names length (%d) must match filepaths length (%d)."
                % (len(sample_names), len(filepaths))
            )
    records = []
    for i, path in enumerate(filepaths):
        if not os.path.exists(path):
            raise FileNotFoundError("copynum TSV not found: %s" % path)
        if sample_names is not None:
            sample = sample_names[i]
        else:
            sample = os.path.basename(path).split(".")[0]

        df = pd.read_csv(path, sep="\t", comment=None)
        # The header line starts with '#Chromosome'; normalise it.
        df = df.rename(columns=lambda c: c.lstrip("#").strip())
        chrom_col = _match_column(df, ["Chromosome", "chrom", "contig", "CHROM"])
        cn_col = _match_column(df, ["CopyNum", "copy_number", "CN"])
        lo_col = _match_column(df, ["CopyNum-CI95-Lower", "ci_low", "CI95_Lower"])
        hi_col = _match_column(df, ["CopyNum-CI95-Upper", "ci_high", "CI95_Upper"])
        if chrom_col is None or cn_col is None:
            raise ValueError(
                "copynum TSV %s missing required columns (need a chromosome "
                "column and a CopyNum column)." % path
            )
        mt_rows = df[df[chrom_col].apply(is_mt_contig)]
        for _, row in mt_rows.iterrows():
            records.append({
                "sample": sample,
                "chrom": row[chrom_col],
                "copy_number": float(row[cn_col]),
                "ci_low": float(row[lo_col]) if lo_col else np.nan,
                "ci_high": float(row[hi_col]) if hi_col else np.nan,
            })

    columns = ["sample", "chrom", "copy_number", "ci_low", "ci_high"]
    if not records:
        return pd.DataFrame(columns=columns)
    out = pd.DataFrame(records, columns=columns)
    return out.sort_values("copy_number", ascending=False).reset_index(drop=True)


def _match_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    """Return the first column in *df* matching any candidate (case-insensitive)."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def read_mito_coverage(
    filepaths: Union[str, Sequence[str]],
    sample_names: Optional[Sequence[str]] = None,
    bins: int = 1000,
    reference: Optional[str] = None,
    contig: Optional[str] = None,
) -> pd.DataFrame:
    """Compute binned mtDNA coverage from one or more BAM/CRAM files.

    Wraps geneview's alignment coverage reader, auto-detecting the
    mitochondrial contig name in each file's header, so a mixed set of BAM
    and CRAM files (single or multi sample) all reduce to one long-format
    depth table suitable for :func:`geneview.mtdna.mito_coverage_plot`.

    Parameters
    ----------
    filepaths : str or sequence of str
        Path(s) to indexed BAM/CRAM files (one per sample).
    sample_names : sequence of str, optional
        Sample labels aligned with *filepaths*.  Derived from file names when
        omitted.
    bins : int, optional
        Number of bins spanning the mitochondrial contig (default 1000).
    reference : str, optional
        Reference FASTA (required for CRAM decoding; ignored for BAM).
    contig : str, optional
        Force a specific mitochondrial contig name.  When ``None`` (default),
        the reader picks the first header contig that looks mitochondrial.

    Returns
    -------
    pandas.DataFrame
        Columns: ``sample, chrom, start, end, pos, depth`` where ``pos`` is
        the bin midpoint (1-based) and ``depth`` the mean coverage in the bin.

    Raises
    ------
    ImportError
        If :mod:`pysam` is not installed.
    ValueError
        If a file has no recognisable mitochondrial contig.
    """
    from ..genometracks._base import GenomicInterval
    from ..genometracks._io import open_alignment_file, _read_alignment_coverage

    if isinstance(filepaths, str):
        filepaths = [filepaths]
    filepaths = list(filepaths)
    if sample_names is not None and len(sample_names) != len(filepaths):
        raise ValueError(
            "sample_names length (%d) must match filepaths length (%d)."
            % (len(sample_names), len(filepaths))
        )

    frames = []
    for i, path in enumerate(filepaths):
        if not os.path.exists(path):
            raise FileNotFoundError("alignment file not found: %s" % path)
        sample = sample_names[i] if sample_names is not None \
            else os.path.basename(path).split(".")[0]
        mode = "rc" if str(path).lower().endswith(".cram") else "rb"

        # Resolve the mitochondrial contig name from the header.
        mt_contig = contig
        aln = open_alignment_file(path, reference=reference, mode=mode)
        try:
            refs = list(aln.references)
            length = None
            if mt_contig is None:
                for name in refs:
                    if is_mt_contig(name):
                        mt_contig = name
                        break
            if mt_contig is None or mt_contig not in refs:
                raise ValueError(
                    "No mitochondrial contig found in %s (contigs: %s). "
                    "Pass contig=... explicitly." % (path, ", ".join(refs[:5]))
                )
            length = aln.get_reference_length(mt_contig)
        finally:
            aln.close()

        region = GenomicInterval(mt_contig, 0, length)
        cov = _read_alignment_coverage(
            path, region=region, bins=bins, mode=mode, reference=reference,
        )
        cov = cov.copy()
        cov["sample"] = sample
        # Bin midpoint as 1-based position for plotting against the gene map.
        cov["pos"] = ((cov["start"] + cov["end"]) / 2.0) + 1
        cov = cov.rename(columns={"value": "depth"})
        frames.append(cov[["sample", "chrom", "start", "end", "pos", "depth"]])

    if not frames:
        return pd.DataFrame(columns=["sample", "chrom", "start", "end", "pos", "depth"])
    return pd.concat(frames, ignore_index=True)
