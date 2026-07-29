"""Human mitochondrial genome (rCRS / NC_012920.1) reference backbone.

This module hard-codes the canonical human mitochondrial DNA (mtDNA) gene
map used by every plot in :mod:`geneview.mtdna`.  Coordinates follow the
revised Cambridge Reference Sequence (rCRS, GenBank ``NC_012920.1``), the
16,569 bp reference that MITOMAP, HelixMTdb, gnomAD (chrM) and MitoQuest are
all coordinated to.  Positions are **1-based, inclusive** — the convention
used in VCF and in the mtDNA literature.

The 37 genes comprise 13 protein-coding genes, 22 transfer RNAs (tRNA) and
2 ribosomal RNAs (rRNA); the non-coding control region (D-loop, which spans
the origin) is included as a distinct feature because most disease- and
population-relevant hypervariable variation lives there.

Author: Shujia Huang
"""
from typing import Dict, List, Optional

# Length of the rCRS mitochondrial genome (bp).
MT_LENGTH = 16569

# Common aliases used for the mitochondrial contig across references.
MT_CONTIG_ALIASES = ("chrM", "MT", "chrMT", "M", "NC_012920.1", "rCRS")

# Feature-type -> default colour (colour-blind-friendly, journal-neutral).
# Consumed as a fall-back when the active plot style does not override it.
MT_FEATURE_COLORS = {
    "protein_coding": "#4C72B0",   # blue
    "tRNA": "#DD8452",             # orange
    "rRNA": "#55A868",             # green
    "control_region": "#C44E52",   # red (D-loop / HVR)
}

# Human-readable feature-type labels for legends.
MT_FEATURE_LABELS = {
    "protein_coding": "Protein-coding",
    "tRNA": "tRNA",
    "rRNA": "rRNA",
    "control_region": "Control region (D-loop)",
}

# ---------------------------------------------------------------------------
# The 37 mtDNA genes + control region (rCRS, 1-based inclusive).
#
# Each entry: (name, start, end, strand, feature_type)
#   * strand ``"+"`` = heavy strand, ``"-"`` = light strand.
#   * The control region spans the origin (16024..16569 and 1..576); it is
#     represented here as two arcs so downstream code never has to special-
#     case the wrap-around.
# ---------------------------------------------------------------------------
_MT_GENE_TABLE = [
    # name,            start,   end,  strand, feature_type
    ("D-loop",             1,     576, "+", "control_region"),
    ("MT-TF",            577,     647, "+", "tRNA"),          # tRNA-Phe
    ("MT-RNR1",          648,    1601, "+", "rRNA"),          # 12S rRNA
    ("MT-TV",           1602,    1670, "+", "tRNA"),          # tRNA-Val
    ("MT-RNR2",         1671,    3229, "+", "rRNA"),          # 16S rRNA
    ("MT-TL1",          3230,    3304, "+", "tRNA"),          # tRNA-Leu(UUR)
    ("MT-ND1",          3307,    4262, "+", "protein_coding"),
    ("MT-TI",           4263,    4331, "+", "tRNA"),          # tRNA-Ile
    ("MT-TQ",           4329,    4400, "-", "tRNA"),          # tRNA-Gln
    ("MT-TM",           4402,    4469, "+", "tRNA"),          # tRNA-Met
    ("MT-ND2",          4470,    5511, "+", "protein_coding"),
    ("MT-TW",           5512,    5579, "+", "tRNA"),          # tRNA-Trp
    ("MT-TA",           5587,    5655, "-", "tRNA"),          # tRNA-Ala
    ("MT-TN",           5657,    5729, "-", "tRNA"),          # tRNA-Asn
    ("MT-TC",           5761,    5826, "-", "tRNA"),          # tRNA-Cys
    ("MT-TY",           5826,    5891, "-", "tRNA"),          # tRNA-Tyr
    ("MT-CO1",          5904,    7445, "+", "protein_coding"),
    ("MT-TS1",          7446,    7514, "-", "tRNA"),          # tRNA-Ser(UCN)
    ("MT-TD",           7518,    7585, "+", "tRNA"),          # tRNA-Asp
    ("MT-CO2",          7586,    8269, "+", "protein_coding"),
    ("MT-TK",           8295,    8364, "+", "tRNA"),          # tRNA-Lys
    ("MT-ATP8",         8366,    8572, "+", "protein_coding"),
    ("MT-ATP6",         8527,    9207, "+", "protein_coding"),
    ("MT-CO3",          9207,    9990, "+", "protein_coding"),
    ("MT-TG",           9991,   10058, "+", "tRNA"),          # tRNA-Gly
    ("MT-ND3",         10059,   10404, "+", "protein_coding"),
    ("MT-TR",          10405,   10469, "+", "tRNA"),          # tRNA-Arg
    ("MT-ND4L",        10470,   10766, "+", "protein_coding"),
    ("MT-ND4",         10760,   12137, "+", "protein_coding"),
    ("MT-TH",          12138,   12206, "+", "tRNA"),          # tRNA-His
    ("MT-TS2",         12207,   12265, "+", "tRNA"),          # tRNA-Ser(AGY)
    ("MT-TL2",         12266,   12336, "+", "tRNA"),          # tRNA-Leu(CUN)
    ("MT-ND5",         12337,   14148, "+", "protein_coding"),
    ("MT-ND6",         14149,   14673, "-", "protein_coding"),
    ("MT-TE",          14674,   14742, "-", "tRNA"),          # tRNA-Glu
    ("MT-CYB",         14747,   15887, "+", "protein_coding"),
    ("MT-TT",          15888,   15953, "+", "tRNA"),          # tRNA-Thr
    ("MT-TP",          15956,   16023, "-", "tRNA"),          # tRNA-Pro
    ("D-loop",         16024,   16569, "+", "control_region"),
]

# The three hypervariable segments of the control region (rCRS, 1-based) in MITOMAP (https://www.mitomap.org/MITOMAP/GenomeLoci).
MT_HYPERVARIABLE_REGIONS = {
    "HVR1": (16024, 16383),
    "HVR2": (57, 372),
    "HVR3": (438, 574),
}


def get_mt_genes(as_dataframe: bool = False):
    """Return the canonical rCRS mtDNA gene map.

    Parameters
    ----------
    as_dataframe : bool, optional
        When True, return a :class:`pandas.DataFrame` with columns
        ``name, start, end, strand, feature_type``.  When False (default),
        return a list of plain dicts (keeps the module import-light for code
        that only needs coordinates).

    Returns
    -------
    list of dict or pandas.DataFrame
        One record per gene / control-region arc, ordered by ``start``.
        The control region appears as two arcs (1..576 and 16024..16569).
    """
    records = [
        {"name": name, "start": start, "end": end,
         "strand": strand, "feature_type": ftype}
        for (name, start, end, strand, ftype) in _MT_GENE_TABLE
    ]
    if as_dataframe:
        import pandas as pd
        return pd.DataFrame(records,
                            columns=["name", "start", "end", "strand", "feature_type"])
    return records


def gene_at(position: int) -> Optional[Dict]:
    """Return the gene/feature overlapping a 1-based mtDNA *position*.

    Parameters
    ----------
    position : int
        A 1-based coordinate on the rCRS (1..16569).

    Returns
    -------
    dict or None
        The first matching feature record (``name, start, end, strand,
        feature_type``), or ``None`` when the position falls in an
        intergenic gap.  When features overlap (e.g. MT-ATP8/MT-ATP6), the
        earliest-starting feature is returned.
    """
    for (name, start, end, strand, ftype) in _MT_GENE_TABLE:
        if start <= position <= end:
            return {"name": name, "start": start, "end": end,
                    "strand": strand, "feature_type": ftype}
    return None


def genes_in_range(start: int, end: int) -> List[Dict]:
    """Return all features overlapping the closed interval ``[start, end]``.

    Parameters
    ----------
    start, end : int
        1-based inclusive bounds on the rCRS.

    Returns
    -------
    list of dict
        Feature records overlapping the interval, ordered by ``start``.
    """
    lo, hi = (start, end) if start <= end else (end, start)
    hits = []
    for (name, gs, ge, strand, ftype) in _MT_GENE_TABLE:
        if gs <= hi and ge >= lo:
            hits.append({"name": name, "start": gs, "end": ge,
                         "strand": strand, "feature_type": ftype})
    return hits


def is_mt_contig(name: str) -> bool:
    """Return True when *name* looks like a mitochondrial contig name.

    Recognises the common aliases (``chrM``, ``MT``, ``chrMT``, ``M``,
    ``NC_012920.1``, ...) case-insensitively.
    """
    if name is None:
        return False
    return str(name).strip().lower() in {a.lower() for a in MT_CONTIG_ALIASES}
