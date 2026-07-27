"""Dandelion plot examples -- fan, circle, pie, pin types.

Demonstrates DandelionTrack from geneview.genometracks with various
shape types, clustering parameters, and customisation options.
Works both standalone and as a composable track inside plot_tracks().

Run:  python examples/scripts/mutation_dandelion.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from geneview.genometracks import DandelionTrack, dandelion_plot, GenomicInterval
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
# 1. Basic fan dandelion plot (standalone via convenience function)
# ---------------------------------------------------------------------------
ax = dandelion_plot(
    snps, features=feats, type="fan",
    title="Dandelion Plot (fan)", style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_dandelion_fan.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_fan.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 2. Circle dandelion plot
# ---------------------------------------------------------------------------
ax = dandelion_plot(
    snps, features=feats, type="circle",
    title="Dandelion Plot (circle)", style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_dandelion_circle.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_circle.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 3. Pin dandelion plot
# ---------------------------------------------------------------------------
ax = dandelion_plot(
    snps, features=feats, type="pin",
    title="Dandelion Plot (pin)", style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_dandelion_pin.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_pin.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 4. Custom maxgaps (tighter clustering)
# ---------------------------------------------------------------------------
ax = dandelion_plot(
    snps, features=feats, type="fan",
    maxgaps=0.05,
    title="Dandelion Plot (maxgaps=0.05)", style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_dandelion_maxgaps.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_maxgaps.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 5. Custom height method (sum of scores)
# ---------------------------------------------------------------------------
ax = dandelion_plot(
    snps, features=feats, type="fan",
    height_method=lambda scores: sum(scores),
    title="Dandelion Plot (height=sum of scores)", style=STYLE,
)
ax.figure.savefig(os.path.join(OUT_DIR, "mutation_dandelion_height_sum.png"),
                  dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_height_sum.png")
plt.close("all")

# ---------------------------------------------------------------------------
# 6. Side-by-side comparison of all types (using draw() on existing axes)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(4, 1, figsize=(12, 14))
region = GenomicInterval("chr1", 0, 1500)
for ax_i, t in zip(axes, ["fan", "circle", "pie", "pin"]):
    track = DandelionTrack(snps, features=feats, type=t, name=f"Type: {t}")
    track._gv_style = _STYLE_OBJ
    ax_i.set_xlim(region.start, region.end)
    track.draw(ax_i, region)
    ax_i.set_title(f"Type: {t}", fontsize=10, fontweight="bold")
fig.suptitle("Dandelion Plot — All Shape Types", fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "mutation_dandelion_all_types.png"),
            dpi=150, bbox_inches="tight")
print("[INFO] Saved mutation_dandelion_all_types.png")
plt.close("all")

print("[INFO] All dandelion figures generated successfully.")
