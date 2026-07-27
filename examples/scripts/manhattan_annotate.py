"""Manhattan plot top-SNP annotation -- exhaustive examples.

This script demonstrates *every* adjustable knob of the significant-site
annotation system used by ``geneview.manhattanplot``:

  * ``is_annotate_topsnp``  -- turn annotation on
  * ``annotate_fmt``        -- label content (None / format string / callable)
  * ``annotate_layout``     -- "repel" (iterative de-overlap) or
                               "lane" (tidy top row + leader lines)
  * ``text_kws``            -- label *styling* (fontsize, color, rotation, ...)
                               plus an optional ``arrowprops`` for the leaders
  * ``adjust_text_kws``     -- fine-tuning of the repel engine
                               (force_text, expand_text, only_move, lim, ...)
  * ``sign_marker_p`` / ``ld_block_size`` -- which sites become labels

Every figure is written to ``examples/figures/`` as ``manhattan_annotate_*.png``.

Run:  python examples/scripts/manhattan_annotate.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import geneview as gv

# ---------------------------------------------------------------------------
# Journal style switch -- change this ONE line to render every figure below
# in a journal-compliant theme.
#   "nature"  -> Nature Research guidelines (thin rules, compact fonts)
#   "science" -> Science / AAAS guidelines
#   "cell"    -> Cell Press guidelines
#   None      -> default geneview look
# ---------------------------------------------------------------------------
STYLE = "nature"

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

XTICK = set(["chr" + c for c in list(map(str, range(1, 15))) + ["16", "18", "20", "22"]])


def _save(ax, name):
    ax.figure.savefig(os.path.join(OUT_DIR, name), dpi=150, bbox_inches="tight")
    print("[INFO] Saved %s" % name)
    plt.close("all")


def _make_dense_gwas(n_chroms=22, per_chrom=600, n_loci=45, seed=7):
    """Synthesize a GWAS table with many independent significant loci.

    The bundled ``gwas`` dataset only carries a couple of genome-wide hits,
    which is not enough to stress-test the "many labels" layouts, so we build a
    denser one here (background noise + ``n_loci`` scattered strong signals).
    """
    rng = np.random.RandomState(seed)
    rows = []
    for c in range(1, n_chroms + 1):
        pos = np.sort(rng.randint(1, 2_000_000, size=per_chrom))
        pval = rng.uniform(1e-3, 1.0, size=per_chrom)
        for i, (p, v) in enumerate(zip(pos, pval)):
            rows.append({"#CHROM": "chr%d" % c, "POS": int(p),
                         "P": float(v), "ID": "rs%d_%d" % (c, i)})
    df = pd.DataFrame(rows)
    hit_idx = rng.choice(len(df), size=n_loci, replace=False)
    df.loc[hit_idx, "P"] = rng.uniform(1e-15, 1e-8, size=n_loci)
    return df


# Real bundled dataset (few significant sites -- good for the "typical" cases).
GWAS = gv.utils.load_dataset("gwas").loc[:, ["#CHROM", "POS", "P", "ID"]]
# Dense synthetic dataset (many loci -- good for the scaling / lane demos).
DENSE = _make_dense_gwas()

COMMON = dict(
    xlabel="Chromosome",
    ylabel=r"$-\log_{10}{(P)}$",
    hline_kws={"linestyle": "--", "lw": 1.0},
    xticklabel_kws={"rotation": 45},
    style=STYLE,
)


# ===========================================================================
# 1. Default annotation -- repel layout, SNP id only, default light arrows
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-6, is_annotate_topsnp=True,
                 xtick_label_set=XTICK, title="1. Default annotation (repel + SNP id)",
                 **COMMON)
_save(ax, "manhattan_annotate_01_default.png")


# ===========================================================================
# 2. Custom label via a FORMAT STRING (fields: snp/chrom/pos/p/log10p)
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-5, is_annotate_topsnp=True,
                 annotate_fmt="{snp}\nP={p:.1e}",
                 text_kws={"fontsize": 8},
                 xtick_label_set=XTICK, title="2. Format-string label (SNP + P-value)",
                 **COMMON)
_save(ax, "manhattan_annotate_02_fmt_string.png")


# ===========================================================================
# 3. Custom label via a CALLABLE formatter
# ===========================================================================
def _locus_label(snp, chrom, pos, p, log10p):
    return "%s:%d\n(-logP=%.1f)" % (chrom, int(pos), log10p)


_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-5, is_annotate_topsnp=True,
                 annotate_fmt=_locus_label,
                 text_kws={"fontsize": 8, "color": "#B22222"},
                 xtick_label_set=XTICK, title="3. Callable label formatter (chrom:pos)",
                 **COMMON)
_save(ax, "manhattan_annotate_03_callable.png")


# ===========================================================================
# 4. Text STYLING via text_kws -- fontsize, color, weight, rotation
#    (rotation now actually applies to the label text; previously ignored)
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-5, is_annotate_topsnp=True,
                 text_kws={"fontsize": 10, "color": "navy",
                           "fontweight": "bold", "rotation": 30},
                 xtick_label_set=XTICK, title="4. Styled text (bold, navy, rotated 30 deg)",
                 **COMMON)
_save(ax, "manhattan_annotate_04_text_style.png")


# ===========================================================================
# 5. Custom ARROWS via text_kws["arrowprops"]
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-5, is_annotate_topsnp=True,
                 annotate_fmt="{snp}",
                 text_kws={"fontsize": 8,
                           "arrowprops": dict(arrowstyle="->", color="#2CA02C",
                                              lw=1.0, connectionstyle="arc3,rad=0.2")},
                 xtick_label_set=XTICK, title="5. Custom curved green arrows",
                 **COMMON)
_save(ax, "manhattan_annotate_05_arrows.png")


# ===========================================================================
# 6. NO arrows -- pass arrowprops=None explicitly
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, sign_marker_p=1e-5, is_annotate_topsnp=True,
                 text_kws={"fontsize": 9, "arrowprops": None},
                 xtick_label_set=XTICK, title="6. No connecting arrows",
                 **COMMON)
_save(ax, "manhattan_annotate_06_no_arrows.png")


# ===========================================================================
# 7. Tune the REPEL ENGINE via adjust_text_kws
#    only_move locks labels to their point's x and moves them only up;
#    stronger force_text + wider expand_text spreads them further apart.
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(DENSE, ax=ax, sign_marker_p=1e-8, is_annotate_topsnp=True,
                 annotate_fmt="{snp}", annotate_layout="repel",
                 text_kws={"fontsize": 7},
                 adjust_text_kws={"only_move": {"points": "y", "text": "xy", "objects": "xy"},
                                  "force_text": (0.5, 0.8),
                                  "expand_text": (1.1, 1.6),
                                  "lim": 300},
                 title="7. Repel engine tuned via adjust_text_kws (dense)",
                 **COMMON)
_save(ax, "manhattan_annotate_07_adjust_kws.png")


# ===========================================================================
# 8. LANE layout -- tidy vertical labels along the top with leader lines
#    (default rotation=90); scales cleanly to many loci.
# ===========================================================================
_, ax = plt.subplots(figsize=(12, 4))
gv.manhattanplot(DENSE, ax=ax, sign_marker_p=1e-8, is_annotate_topsnp=True,
                 annotate_fmt="{snp}", annotate_layout="lane",
                 text_kws={"fontsize": 7},
                 title="8. Lane layout, vertical labels (dense)",
                 **COMMON)
_save(ax, "manhattan_annotate_08_lane_vertical.png")


# ===========================================================================
# 9. LANE layout, HORIZONTAL labels (rotation=0) with P-value + custom arrow
# ===========================================================================
_, ax = plt.subplots(figsize=(14, 4))
gv.manhattanplot(DENSE, ax=ax, sign_marker_p=1e-8, is_annotate_topsnp=True,
                 annotate_fmt="{snp} ({p:.0e})", annotate_layout="lane",
                 text_kws={"fontsize": 6, "rotation": 0,
                           "arrowprops": dict(arrowstyle="-", color="0.4", lw=0.4)},
                 title="9. Lane layout, horizontal labels + P-value (dense)",
                 **COMMON)
_save(ax, "manhattan_annotate_09_lane_horizontal.png")


# ===========================================================================
# 10. Side-by-side: repel vs lane on the SAME dense data
# ===========================================================================
fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(13, 8))
gv.manhattanplot(DENSE, ax=ax_top, sign_marker_p=1e-8, is_annotate_topsnp=True,
                 annotate_fmt="{snp}", annotate_layout="repel",
                 text_kws={"fontsize": 6}, title="repel layout", **COMMON)
gv.manhattanplot(DENSE, ax=ax_bot, sign_marker_p=1e-8, is_annotate_topsnp=True,
                 annotate_fmt="{snp}", annotate_layout="lane",
                 text_kws={"fontsize": 6}, title="lane layout", **COMMON)
fig.suptitle("10. repel vs lane (identical dense data)", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "manhattan_annotate_10_repel_vs_lane.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved manhattan_annotate_10_repel_vs_lane.png")
plt.close("all")


# ===========================================================================
# 11. Single-chromosome zoom (CHR=) with annotation + LD-block control
#     ld_block_size groups nearby significant SNPs so each block gets one label.
# ===========================================================================
_, ax = plt.subplots(figsize=(10, 4))
gv.manhattanplot(GWAS, ax=ax, CHR="chr8", sign_marker_p=1e-4,
                 is_annotate_topsnp=True, annotate_fmt="{snp}\nP={p:.1e}",
                 ld_block_size=100000, text_kws={"fontsize": 8},
                 xlabel="Position on chr8", ylabel=r"$-\log_{10}{(P)}$",
                 hline_kws={"linestyle": "--", "lw": 1.0},
                 title="11. Single chromosome (chr8) with block labels",
                 style=STYLE)
_save(ax, "manhattan_annotate_11_single_chrom.png")


print("[INFO] All manhattan annotation figures generated successfully.")
