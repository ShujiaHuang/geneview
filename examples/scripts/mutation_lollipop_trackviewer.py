"""Lolliplot examples mirroring the trackViewer vignette.

Reproduces the examples from
https://jianhong.github.io/trackViewer/articles/lollipopPlot.html
using geneview's :class:`LolliplotTrack` and :func:`lolliplot`.

Run:  python examples/scripts/mutation_lollipop_trackviewer.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from geneview.genometracks import (
    LolliplotTrack,
    lolliplot,
    GenomeAxisTrack,
    AnnotationTrack,
    plot_tracks,
    GenomicInterval,
)

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# ===================================================================
# Shared data (mirrors trackViewer vignette data)
# ===================================================================
# SNP positions — same as the R vignette
SNP = [10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402]

sample_df = pd.DataFrame({
    "chrom": ["chr1"] * len(SNP),
    "start": SNP,
    "label": [f"snp{s}" for s in SNP],
})

# Features (protein domains) — same as the R vignette
features_df = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 501, 1001],
    "end":   [120, 900, 1405],
    "name":  ["block1", "block2", "block3"],
})

region_full = GenomicInterval("chr1", 0, 1500)


def _save(fig, name):
    path = os.path.join(OUT_DIR, name)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Saved {name}")
    plt.close("all")


# ===================================================================
# 1. Basic lolliplot  (vignette §Lolliplot)
# ===================================================================
ax = lolliplot(sample_df, features=features_df,
               figsize=(12, 4), title="Basic Lolliplot")
_save(ax.figure, "mutation_tv_01_basic.png")

# ===================================================================
# 2. Define the range  (vignette: ranges = GRanges("chr1", IRanges(104, 109)))
# ===================================================================
region_zoom = GenomicInterval("chr1", 104, 109)
ax = lolliplot(sample_df, features=features_df,
               region=region_zoom,
               figsize=(12, 4), title="Lolliplot — Zoomed to 104-109")
_save(ax.figure, "mutation_tv_02_range.png")

# ===================================================================
# 3. Change the color of the features  (vignette §Change the color of the features)
# ===================================================================
features_colored = features_df.copy()
features_colored["fill"] = ["#FF8833", "#51C6E6", "#DFA32D"]
ax = lolliplot(sample_df, features=features_colored,
               figsize=(12, 4), title="Custom Feature Colors")
_save(ax.figure, "mutation_tv_03_feature_colors.png")

# ===================================================================
# 4. Change the color and opacity of the lollipop
#    (vignette §Change the color and opacity of the lollipop)
# ===================================================================
np.random.seed(42)
palette = ["#0080FF", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9"]
sample_color = sample_df.copy()
sample_color["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
sample_color["border"] = np.random.choice(["#CCCCCC", "#4D4D4D"], len(SNP))
sample_color["alpha"] = np.random.uniform(0.4, 1.0, len(SNP))
ax = lolliplot(sample_color, features=features_colored,
               figsize=(12, 4), title="Custom SNP Colors & Opacity")
_save(ax.figure, "mutation_tv_04_snp_colors.png")

# ===================================================================
# 5. Change the height of the features  (vignette §Change the height)
# ===================================================================
features_height = features_colored.copy()
features_height["height"] = [0.03, 0.06, 0.09]
ax = lolliplot(sample_color, features=features_height,
               figsize=(12, 4), title="Custom Feature Heights")
_save(ax.figure, "mutation_tv_05_feature_height.png")

# ===================================================================
# 6. Change the height of a lollipop plot  (vignette §Change the height)
#    score < 10 → Tanghulu stacking
# ===================================================================
np.random.seed(42)
sample_score = sample_df.copy()
sample_score["score"] = np.random.randint(1, 6, len(SNP))
ax = lolliplot(sample_score, features=features_colored,
               figsize=(12, 4), title="Score-based Height (Tanghulu)")
_save(ax.figure, "mutation_tv_06_score_tanghulu.png")

# ===================================================================
# 7. Remove y-axis  (vignette: yaxis=FALSE)
# ===================================================================
track = LolliplotTrack(sample_score, features=features_colored,
                       show_yaxis=False)
axes = plot_tracks([track], region=region_full, figsize=(12, 4),
                   title="Lolliplot — No Y-axis")
_save(axes[0].figure, "mutation_tv_07_no_yaxis.png")

# ===================================================================
# 8. Score > style switch limit  (vignette: lollipop_style_switch_limit)
# ===================================================================
np.random.seed(42)
sample_big = sample_df.copy()
sample_big["score"] = np.random.randint(1, 16, len(SNP))
ax = lolliplot(sample_big, features=features_colored,
               figsize=(12, 4), title="Large Scores (Non-Tanghulu)")
_save(ax.figure, "mutation_tv_08_score_large.png")

# Increase the switch limit so Tanghulu is still used
ax = lolliplot(sample_big, features=features_colored,
               lollipop_style_switch_limit=20,
               figsize=(12, 4), title="Large Scores (Tanghulu, limit=20)")
_save(ax.figure, "mutation_tv_08b_score_limit20.png")

# ===================================================================
# 9. Customise the x-axis  (vignette §Customize the x-axis label position)
#    Use GenomeAxisTrack alongside LolliplotTrack
# ===================================================================
tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(sample_score, features=features_colored),
]
axes = plot_tracks(tracks, region=region_full, figsize=(12, 5),
                   title="Lolliplot + GenomeAxis")
_save(axes[0].figure, "mutation_tv_09_xaxis.png")

# ===================================================================
# 10. Label on feature  (vignette: label_on_feature=TRUE)
# ===================================================================
ax = lolliplot(sample_score, features=features_colored,
               label_on_feature=True,
               figsize=(12, 4), title="Labels on Features")
_save(ax.figure, "mutation_tv_10_label_on_feature.png")

# ===================================================================
# 11. Google pin type  (vignette: type="pin")
# ===================================================================
sample_pin = sample_color.copy()
sample_pin["score"] = np.random.randint(1, 6, len(SNP))
ax = lolliplot(sample_pin, features=features_colored,
               type="pin", figsize=(12, 4), title="Lolliplot (pin)")
_save(ax.figure, "mutation_tv_11_pin.png")

# ===================================================================
# 12. Flag type  (vignette: type="flag")
# ===================================================================
sample_flag = sample_color.copy()
sample_flag["label"] = [f"snp{s}" for s in SNP]
ax = lolliplot(sample_flag, features=features_colored,
               type="flag", figsize=(14, 4), title="Lolliplot (flag)")
_save(ax.figure, "mutation_tv_12_flag.png")

# ===================================================================
# 13. Pie plot  (vignette: type="pie")
# ===================================================================
np.random.seed(42)
sample_pie = sample_df.copy()
x = np.random.randint(20, 100, len(SNP))
sample_pie["pie_values"] = [[int(v), 100 - int(v)] for v in x]
sample_pie["pie_colors"] = [["#87CEFA", "#98CE31"]] * len(SNP)
sample_pie["border"] = "#4D4D4D"
ax = lolliplot(sample_pie, features=features_colored,
               type="pie", figsize=(12, 4), title="Lolliplot (pie)")
_save(ax.figure, "mutation_tv_13_pie.png")

# ===================================================================
# 14. Change the node size  (vignette: cex)
# ===================================================================
ax = lolliplot(sample_score, features=features_colored,
               cex=0.5, figsize=(12, 4), title="Lolliplot (cex=0.5)")
_save(ax.figure, "mutation_tv_14_cex_small.png")

ax = lolliplot(sample_score, features=features_colored,
               cex=1.8, figsize=(12, 5), title="Lolliplot (cex=1.8)")
_save(ax.figure, "mutation_tv_14b_cex_large.png")

# ===================================================================
# 15. Plot the lollipop plot with annotation and axis tracks
#     (vignette §Plot the lollipop plot with the coverage and annotation tracks)
# ===================================================================
tracks = [
    GenomeAxisTrack(),
    AnnotationTrack(features_colored, name="Domains", shape="box",
                    show_label=True),
    LolliplotTrack(sample_score, features=features_colored,
                   name="Mutations"),
]
axes = plot_tracks(tracks, region=region_full, figsize=(12, 7),
                   title="Lolliplot + Annotation + Axis")
_save(axes[0].figure, "mutation_tv_15_with_tracks.png")

# ===================================================================
# 16. Multiple LolliplotTracks stacked  (vignette §Multiple layers)
# ===================================================================
np.random.seed(123)
sample2 = sample_df.copy()
sample2["score"] = np.random.randint(1, 8, len(SNP))
sample2["fill"] = "#DB7575"
sample2["border"] = "#4D4D4D"

tracks = [
    GenomeAxisTrack(),
    LolliplotTrack(sample_score, features=features_colored,
                   name="Sample A"),
    LolliplotTrack(sample2, features=features_colored,
                   name="Sample B"),
]
axes = plot_tracks(tracks, region=region_full, figsize=(12, 8),
                   title="Two Lolliplot Tracks (Multi-layer)")
_save(axes[0].figure, "mutation_tv_16_multi_layer.png")

# ===================================================================
# 17. Side-by-side comparison of all types
# ===================================================================
fig, axes = plt.subplots(4, 1, figsize=(14, 16))
sample_all = sample_color.copy()
sample_all["score"] = np.random.randint(1, 6, len(SNP))
sample_all["label"] = [f"snp{s}" for s in SNP]

x_pie = np.random.randint(20, 100, len(SNP))
sample_all["pie_values"] = [[int(v), 100 - int(v)] for v in x_pie]
sample_all["pie_colors"] = [["#87CEFA", "#98CE31"]] * len(SNP)

for ax_i, t in zip(axes, ["circle", "pie", "pin", "flag"]):
    track = LolliplotTrack(sample_all, features=features_colored,
                           type=t, name=f"Type: {t}")
    ax_i.set_xlim(region_full.start, region_full.end)
    track.draw(ax_i, region_full)
    ax_i.set_title(f"Type: {t}", fontsize=10, fontweight="bold")

fig.suptitle("Lolliplot — All Shape Types (trackViewer style)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
_save(fig, "mutation_tv_17_all_types.png")

# ===================================================================
# 18. Comprehensive example: scores + custom colors + tracks
# ===================================================================
np.random.seed(7)
sample_full = pd.DataFrame({
    "chrom": ["chr1"] * len(SNP),
    "start": SNP,
    "label": [f"snp{s}" for s in SNP],
    "score": np.random.randint(1, 12, len(SNP)),
    "fill": [palette[i % len(palette)] for i in range(len(SNP))],
    "border": np.random.choice(["#CCCCCC", "#4D4D4D"], len(SNP)),
    "alpha": np.random.uniform(0.6, 1.0, len(SNP)),
})

features_full = pd.DataFrame({
    "chrom": ["chr1", "chr1", "chr1"],
    "start": [1, 501, 1001],
    "end":   [120, 900, 1405],
    "name":  ["block1", "block2", "block3"],
    "fill":  ["#FF8833", "#51C6E6", "#DFA32D"],
    "height": [0.04, 0.06, 0.08],
})

tracks = [
    GenomeAxisTrack(),
    AnnotationTrack(features_full, name="Domains", shape="box",
                    show_label=True),
    LolliplotTrack(sample_full, features=features_full,
                   lollipop_style_switch_limit=15,
                   name="Mutations"),
]
axes = plot_tracks(tracks, region=region_full, figsize=(14, 8),
                   title="Comprehensive Lolliplot (trackViewer style)")
_save(axes[0].figure, "mutation_tv_18_comprehensive.png")

print("[INFO] All trackViewer-style lolliplot figures generated successfully.")


# ===================================================================
# NEW FEATURES — geneview-specific enhancements
# ===================================================================

# ===================================================================
# 19. Node labels inside circles
# ===================================================================
np.random.seed(42)
sample_nl = sample_df.copy()
sample_nl["score"] = np.random.randint(1, 6, len(SNP))
sample_nl["node_label"] = [str(s) for s in SNP]
sample_nl["node_label_color"] = "white"
sample_nl["node_label_size"] = 5
sample_nl["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_nl, features=features_colored,
               figsize=(12, 4), title="Node Labels Inside Circles")
_save(ax.figure, "mutation_tv_19_node_labels.png")

# ===================================================================
# 20. Per-SNP cex variation
# ===================================================================
sample_pcex = sample_df.copy()
sample_pcex["score"] = np.random.randint(1, 6, len(SNP))
sample_pcex["cex"] = np.linspace(0.4, 2.0, len(SNP))
sample_pcex["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_pcex, features=features_colored,
               figsize=(12, 4), title="Per-SNP cex Variation")
_save(ax.figure, "mutation_tv_20_per_snp_cex.png")

# ===================================================================
# 21. Caterpillar layout (side=top/bottom)
# ===================================================================
sample_cat = sample_df.copy()
sample_cat["score"] = np.random.randint(1, 5, len(SNP))
sample_cat["side"] = ["top" if i % 2 == 0 else "bottom" for i in range(len(SNP))]
sample_cat["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_cat, features=features_colored,
               figsize=(12, 5), title="Caterpillar Layout (top/bottom)")
_save(ax.figure, "mutation_tv_21_caterpillar.png")

# ===================================================================
# 22. Multiple shapes in Tanghulu stack
# ===================================================================
sample_ms = pd.DataFrame({
    "chrom": ["chr1"] * 4,
    "start": [100, 300, 600, 900],
    "score": [3, 4, 2, 5],
    "fill": [["#FF0000", "#0000FF", "#00FF00"],
             ["#FF00FF", "#FFFF00"],
             ["#0080FF", "#E69F00", "#009E73"],
             ["#D55E00", "#CC79A7"]],
    "shape": [["circle", "square", "diamond"],
              ["triangle_point_up", "circle"],
              ["square", "diamond", "circle"],
              ["triangle_point_down", "triangle_point_up"]],
})
ax = lolliplot(sample_ms, features=features_colored,
               figsize=(12, 4), title="Multi-shape Tanghulu Stack")
_save(ax.figure, "mutation_tv_22_multi_shape_tanghulu.png")

# ===================================================================
# 23. Legend
# ===================================================================
sample_leg = sample_df.copy()
sample_leg["score"] = np.random.randint(1, 6, len(SNP))
sample_leg["fill"] = ["#87CEFA" if i % 2 == 0 else "#98CE31" for i in range(len(SNP))]
legend = {"labels": ["Wild Type", "Mutant"], "fill": ["#87CEFA", "#98CE31"]}
ax = lolliplot(sample_leg, features=features_colored,
               legend=legend, figsize=(12, 4), title="Lolliplot with Legend")
_save(ax.figure, "mutation_tv_23_legend.png")

# ===================================================================
# 24. Label jitter (aligned labels with anti-overlap)
# ===================================================================
sample_jit = sample_df.copy()
np.random.seed(42)
sample_jit["score"] = np.random.randint(1, 6, len(SNP))
sample_jit["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
# Per-SNP dashline colour (matches trackViewer: sample.gr$dashline.col <- sample.gr$color)
sample_jit["dashline_col"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_jit, features=features_colored,
               jitter="label", figsize=(12, 4),
               title="Aligned Labels with Anti-overlap (jitter='label')")
_save(ax.figure, "mutation_tv_24_jitter.png")

# ===================================================================
# 25. Custom y-axis ticks + ylab
# ===================================================================
sample_ya = sample_df.copy()
sample_ya["score"] = np.random.randint(1, 11, len(SNP))
ax = lolliplot(sample_ya, features=features_colored,
               yaxis=[0, 5, 10], ylab="# evidences",
               figsize=(12, 4), title="Custom Y-axis + ylab")
_save(ax.figure, "mutation_tv_25_yaxis_ylab.png")

# ===================================================================
# 26. pie.stack type
# ===================================================================
sample_ps = pd.DataFrame({
    "chrom": ["chr1"] * 6,
    "start": [100, 100, 100, 400, 400, 400],
    "score": [1, 1, 1, 1, 1, 1],
    "pie_values": [[3, 7], [5, 5], [2, 8], [6, 4], [1, 9], [4, 6]],
    "pie_colors": [["#0080FF", "#E69F00"]] * 6,
    "stack_factor": [1, 2, 3, 1, 2, 3],
    "label": ["a", "b", "c", "d", "e", "f"],
})
ax = lolliplot(sample_ps, type="pie.stack",
               figsize=(12, 4), title="pie.stack Type")
_save(ax.figure, "mutation_tv_26_pie_stack.png")

# ===================================================================
# 27. Rescale coordinate mapping
# ===================================================================
sample_rs = pd.DataFrame({
    "chrom": ["chr1"] * 5,
    "start": [100, 300, 500, 700, 900],
    "score": [3, 5, 2, 7, 4],
    "label": ["a", "b", "c", "d", "e"],
    "fill": [palette[i % len(palette)] for i in range(5)],
})
feats_rs = pd.DataFrame({
    "chrom": ["chr1", "chr1"],
    "start": [50, 600],
    "end":   [400, 950],
    "name":  ["region A", "region B"],
})
rescale_map = [(0, 500, 0, 250), (500, 1000, 300, 600)]
ax = lolliplot(sample_rs, features=feats_rs, rescale=rescale_map,
               figsize=(12, 4), title="Rescale Coordinate Mapping")
_save(ax.figure, "mutation_tv_27_rescale.png")

# ===================================================================
# 28. Multi-layer features
# ===================================================================
feats_ml = pd.DataFrame({
    "chrom": ["chr1"] * 5,
    "start": [1, 501, 1001, 200, 700],
    "end":   [120, 900, 1405, 400, 850],
    "name":  ["Domain A", "Kinase", "DNA Bind", "Motif X", "Motif Y"],
    "fill":  ["#FF8833", "#51C6E6", "#DFA32D", "#CC79A7", "#009E73"],
    "feature_layer_id": [1, 1, 1, 2, 2],
})
ax = lolliplot(sample_score, features=feats_ml,
               figsize=(12, 5), title="Multi-layer Features")
_save(ax.figure, "mutation_tv_28_multi_layer_features.png")

# ===================================================================
# 29. Different shapes per SNP
# ===================================================================
sample_sh = sample_df.copy()
sample_sh["score"] = np.random.randint(1, 6, len(SNP))
sample_sh["shape"] = ["circle", "square", "diamond",
                      "triangle_point_up", "triangle_point_down",
                      "circle", "square", "diamond",
                      "triangle_point_up", "triangle_point_down",
                      "circle", "square", "diamond"]
sample_sh["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_sh, features=features_colored,
               figsize=(12, 4), title="Different Shapes Per SNP")
_save(ax.figure, "mutation_tv_29_shapes.png")

# ===================================================================
# 30. Per-SNP label rotation + color
# ===================================================================
sample_lr = sample_df.copy()
sample_lr["score"] = np.random.randint(1, 6, len(SNP))
sample_lr["label_rotation"] = [0, 30, 45, 60, 90, 0, 30, 45, 60, 90, 0, 30, 45]
sample_lr["label_color"] = ["red", "blue", "green", "orange", "purple",
                            "red", "blue", "green", "orange", "purple",
                            "red", "blue", "green"]
sample_lr["fill"] = [palette[i % len(palette)] for i in range(len(SNP))]
ax = lolliplot(sample_lr, features=features_colored,
               figsize=(12, 4), title="Per-SNP Label Rotation & Color")
_save(ax.figure, "mutation_tv_30_label_rotation.png")

print("[INFO] All new-feature figures generated successfully.")
