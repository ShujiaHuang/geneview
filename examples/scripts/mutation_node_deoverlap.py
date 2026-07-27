"""trackViewer-style automatic node de-overlap for mutation tracks.

When several variants sit closer together than one node diameter, their
shapes and labels would collide and become unreadable.  Following
trackViewer's ``lolliplot()``, :class:`LolliplotTrack` automatically spreads
such nodes apart horizontally (re-centred on the cluster centroid, clamped to
the visible region) while distant nodes keep their true genomic coordinate.
:class:`DandelionTrack` handles the same situation by fanning clustered nodes
out from a shared stem.

This script highlights that behaviour on deliberately dense clusters and
across every shape type / layout.

Run:  python examples/scripts/mutation_node_deoverlap.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from geneview.genometracks import (
    LolliplotTrack,
    DandelionTrack,
    lolliplot,
    dandelion_plot,
    GenomicInterval,
)
from geneview.plotstyle import get_style

# ---------------------------------------------------------------------------
# Journal style switch -- change this ONE line to render every figure below
# in a journal-compliant theme.  Node borders, stems, connector/guide lines,
# axis rules and font sizes are all routed through the chosen style.
#   "nature"  -> Nature Research guidelines (thin 0.4pt rules, compact fonts)
#   "science" -> Science / AAAS guidelines
#   "cell"    -> Cell Press guidelines
#   None      -> default geneview look (thicker rules, larger fonts)
# ---------------------------------------------------------------------------
STYLE = "nature"
_STYLE_OBJ = get_style(STYLE) if STYLE else None

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Shared data: three tight clusters (100/105/108, 400/410/420, 1400/1402)
# plus a few isolated variants that should NOT move.
# ---------------------------------------------------------------------------
SNP = [10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402]

snps = pd.DataFrame({
    "chrom": ["chr1"] * len(SNP),
    "start": SNP,
    "label": [f"snp{s}" for s in SNP],
})

features = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 501, 1001],
    "end":   [120, 900, 1405],
    "name":  ["block1", "block2", "block3"],
    "fill":  ["#FF8833", "#51C6E6", "#DFA32D"],
})

region = GenomicInterval("chr1", 0, 1500)


def _save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, name), dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {name}")
    plt.close("all")


# ---------------------------------------------------------------------------
# 1. Circle lolliplot: clustered nodes spread, labels stay legible
# ---------------------------------------------------------------------------
ax = lolliplot(snps, features=features,
               title="Node de-overlap — close variants spread apart",
               style=STYLE)
_save(ax.figure, "mutation_deoverlap_circle.png")

# ---------------------------------------------------------------------------
# 2. Tanghulu stacking: each variant keeps its own stack after spreading
# ---------------------------------------------------------------------------
np.random.seed(42)
snps_score = snps.copy()
snps_score["score"] = np.random.randint(1, 6, len(SNP))
ax = lolliplot(snps_score, features=features,
               title="Node de-overlap — Tanghulu stacking",
               style=STYLE)
_save(ax.figure, "mutation_deoverlap_tanghulu.png")

# ---------------------------------------------------------------------------
# 3. Every shape type keeps clusters readable
# ---------------------------------------------------------------------------
snps_pie = snps.copy()
vals = np.random.randint(20, 100, len(SNP))
snps_pie["pie_values"] = [[int(v), 100 - int(v)] for v in vals]
snps_pie["pie_colors"] = [["#87CEFA", "#98CE31"]] * len(SNP)

fig, axes = plt.subplots(3, 1, figsize=(12, 10))
for ax_i, shape_type in zip(axes, ["pin", "pie", "flag"]):
    data = snps_pie if shape_type == "pie" else snps
    ax_i.set_xlim(region.start, region.end)
    track = LolliplotTrack(data, features=features, type=shape_type,
                           name=shape_type)
    track._gv_style = _STYLE_OBJ
    track.draw(ax_i, region)
    ax_i.set_title(f"type='{shape_type}'", fontsize=10, fontweight="bold")
fig.suptitle("Node de-overlap across shape types",
             fontsize=13, fontweight="bold")
fig.tight_layout()
_save(fig, "mutation_deoverlap_shape_types.png")

# ---------------------------------------------------------------------------
# 4. Caterpillar (two-sided): top labels grow up, bottom labels grow down
# ---------------------------------------------------------------------------
snps_side = snps_score.copy()
snps_side["side"] = ["top", "bottom"] * (len(SNP) // 2) + ["top"] * (len(SNP) % 2)
ax = lolliplot(snps_side, features=features, figsize=(12, 5),
               title="Node de-overlap — caterpillar (two-sided) layout",
               style=STYLE)
_save(ax.figure, "mutation_deoverlap_caterpillar.png")

# ---------------------------------------------------------------------------
# 5. Dandelion: clusters fan out from a shared stem (complementary approach)
# ---------------------------------------------------------------------------
ax = dandelion_plot(snps, features=features, type="fan",
                    title="Dandelion — clusters fan out from a shared stem",
                    style=STYLE)
_save(ax.figure, "mutation_deoverlap_dandelion.png")

print("[INFO] All node de-overlap figures generated successfully.")
