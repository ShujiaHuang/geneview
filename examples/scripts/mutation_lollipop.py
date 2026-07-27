"""Lollipop plot examples -- circle, pie, pin, flag types.

Demonstrates LolliplotTrack from geneview.genometracks with various
shape types and customisation options.  Works both standalone and as
a composable track inside plot_tracks().

Run:  python examples/scripts/mutation_lollipop.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from geneview.genometracks import LolliplotTrack, lolliplot, GenomicInterval
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

# ---------------------------------------------------------------------------
# Locate example data
# ---------------------------------------------------------------------------
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "mutation")
SNP_FILE = os.path.join(DATA_DIR, "snps_example.tsv")
FEAT_FILE = os.path.join(DATA_DIR, "features_example.tsv")

snps = pd.read_csv(SNP_FILE, sep="\t")
feats = pd.read_csv(FEAT_FILE, sep="\t")

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Basic circle lolliplot (standalone via convenience function)
# ---------------------------------------------------------------------------
ax = lolliplot(snps, features=feats, type="circle",
               title="Lolliplot (circle)", style=STYLE)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_lollipop_circle.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_circle.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Pie chart lolliplot
# ---------------------------------------------------------------------------
ax = lolliplot(snps, features=feats, type="pie",
               title="Lolliplot (pie)", style=STYLE)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_lollipop_pie.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_pie.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Pin lolliplot
# ---------------------------------------------------------------------------
ax = lolliplot(snps, features=feats, type="pin",
               title="Lolliplot (pin)", style=STYLE)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_lollipop_pin.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_pin.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Flag lolliplot
# ---------------------------------------------------------------------------
ax = lolliplot(snps, features=feats, type="flag",
               title="Lolliplot (flag)", style=STYLE)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_lollipop_flag.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_flag.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. Customised lolliplot with larger shapes
# ---------------------------------------------------------------------------
ax = lolliplot(
    snps, features=feats, type="circle",
    cex=1.5, dashline_col="#AAAAAA",
    label_on_feature=True,
    figsize=(14, 5), title="Lolliplot (customised, cex=1.5)",
    style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_lollipop_custom.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_custom.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 6. Side-by-side comparison of all types (using draw() on existing axes)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(4, 1, figsize=(12, 14))
region = GenomicInterval("chr1", 0, 1500)
for ax_i, t in zip(axes, ["circle", "pie", "pin", "flag"]):
    track = LolliplotTrack(snps, features=feats, type=t, name=f"Type: {t}")
    track._gv_style = _STYLE_OBJ
    ax_i.set_xlim(region.start, region.end)
    track.draw(ax_i, region)
    ax_i.set_title(f"Type: {t}", fontsize=10, fontweight="bold")
fig.suptitle("Lolliplot — All Shape Types", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "mutation_lollipop_all_types.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_lollipop_all_types.png")
plt.close("all")

print("[INFO] All lollipop figures generated successfully.")
