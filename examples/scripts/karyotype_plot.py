"""Whole-genome karyotype plots with ``geneview.karyoplot``.

:func:`geneview.karyoplot` draws a whole-genome ideogram overview: every
chromosome is stacked as a horizontal bar and cytobands are colored by their
``gieStain`` code (shared with :class:`~geneview.genometracks.IdeogramTrack`
through the canonical ``geneview.palette`` color map).

It accepts the same cytoband inputs as the rest of geneview via the shared
``read_cytoband`` loader: a UCSC/Gviz karyotype file, a ``pandas.DataFrame``,
or an array of ``chrom, chromStart, chromEnd, name, gieStain`` rows.

Run:  python examples/scripts/karyotype_plot.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geneview as gv
from geneview.utils import read_cytoband

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Main assembly chromosomes (drop the alt/random/unplaced scaffolds so the
# whole-genome overview stays readable).
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def _save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, name), dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {name}")
    plt.close(fig)


# Load the human hg38 cytoband table from the geneview-data repository. The
# file uses the UCSC '#chrom' header, which ``read_cytoband`` handles directly.
karyotype_fn = gv.load_dataset("karyotype_human_hg38.txt")
bands = read_cytoband(karyotype_fn)
main_bands = bands[bands["chrom"].isin(MAIN_CHROMS)]

# ---------------------------------------------------------------------------
# 1. Whole-genome karyotype (all main chromosomes)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(16, 8))
gv.karyoplot(main_bands, ax=ax)
ax.set_title("Human karyotype (hg38) — whole genome", fontsize=12, loc="left")
_save(fig, "karyotype_whole_genome.png")

# ---------------------------------------------------------------------------
# 2. A single chromosome via the ``CHR`` parameter
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 2.2))
gv.karyoplot(main_bands, ax=ax, CHR="chr7", width=0.6)
ax.set_title("Single chromosome (chr7)", fontsize=12, loc="left")
_save(fig, "karyotype_single_chr.png")

# ---------------------------------------------------------------------------
# 3. Custom styling: bar width, band outlines, and a highlight color for
#    stains that are absent from the standard palette.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(16, 8))
gv.karyoplot(
    main_bands, ax=ax,
    width=0.7,
    color4none="#F4A582",      # color for any undefined gieStain
    edgecolor="#444444",       # forwarded to matplotlib Rectangle via **kwargs
    linewidth=0.3,
)
ax.set_title("Custom styling (outlined bands, wider bars)",
             fontsize=12, loc="left")
_save(fig, "karyotype_custom_style.png")

# ---------------------------------------------------------------------------
# 4. Plotting from an inline DataFrame / array of rows (no file needed)
# ---------------------------------------------------------------------------
toy = [
    ["chr1", 0, 50_000_000, "p", "gneg"],
    ["chr1", 50_000_000, 55_000_000, "cen", "acen"],
    ["chr1", 55_000_000, 120_000_000, "q", "gpos75"],
    ["chr2", 0, 40_000_000, "p", "gpos50"],
    ["chr2", 40_000_000, 44_000_000, "cen", "acen"],
    ["chr2", 44_000_000, 100_000_000, "q", "gneg"],
]
fig, ax = plt.subplots(figsize=(12, 2.4))
gv.karyoplot(toy, ax=ax, width=0.6)
ax.set_title("Karyotype from an inline array of bands", fontsize=12, loc="left")
_save(fig, "karyotype_from_array.png")

print("[INFO] All karyotype figures generated successfully.")
