"""Plotstyle examples — applying journal-specific styles to geneview plots.

Demonstrates how to use geneview's built-in plot styles (geneview, nature,
science, cell) with Manhattan plots, Q-Q plots, Venn diagrams, and
Admixture plots.

Usage::

    python examples/scripts/plotstyle.py

"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geneview as gv
from geneview.plotstyle import apply_style, use_style, list_styles


# ---------------------------------------------------------------------------
# 0. Show available styles
# ---------------------------------------------------------------------------
print("Available plot styles:", list_styles())
# Expected: ['cell', 'geneview', 'nature', 'science']


# ---------------------------------------------------------------------------
# 1. Manhattan plots in different styles
# ---------------------------------------------------------------------------
df = gv.utils.load_dataset("gwas")

for style_name in ["geneview", "nature", "science", "cell"]:
    fig, ax = plt.subplots(figsize=(10, 3.5), facecolor="w", edgecolor="k")
    ax = gv.manhattanplot(
        data=df,
        marker=".",
        sign_marker_p=1e-6,
        sign_marker_color="r",
        hline_kws={"linestyle": "--", "lw": 1.0},
        xlabel="Chromosome",
        ylabel=r"$-log_{10}{(P)}$",
        title=f"Manhattan plot — {style_name} style",
        style=style_name,
        ax=ax,
    )
    fig.tight_layout()
    fig.savefig(f"examples/figures/manhattan_style_{style_name}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: manhattan_style_{style_name}.png")


# ---------------------------------------------------------------------------
# 2. Q-Q plots in different styles
# ---------------------------------------------------------------------------
for style_name in ["geneview", "nature", "science", "cell"]:
    fig, ax = plt.subplots(figsize=(5, 5), facecolor="w", edgecolor="k")
    ax = gv.qqplot(
        data=df["P"],
        marker="o",
        title=f"QQ plot — {style_name}",
        xlabel=r"Expected $-log_{10}{(P)}$",
        ylabel=r"Observed $-log_{10}{(P)}$",
        style=style_name,
        ax=ax,
    )
    fig.tight_layout()
    fig.savefig(f"examples/figures/qq_style_{style_name}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: qq_style_{style_name}.png")


# ---------------------------------------------------------------------------
# 3. Venn diagrams in different styles
# ---------------------------------------------------------------------------
venn_data = {
    "DEG Study A": {"TP53", "BRCA1", "EGFR", "KRAS", "MYC", "PTEN", "RB1"},
    "DEG Study B": {"BRCA1", "EGFR", "ALK", "BRAF", "MYC", "PIK3CA"},
    "DEG Study C": {"TP53", "KRAS", "BRAF", "NRAS", "IDH1"},
}

for style_name in ["geneview", "nature", "science", "cell"]:
    fig, ax = plt.subplots(figsize=(6, 6))
    ax = gv.venn(
        data=venn_data,
        palette="Set2",
        fontsize=11,
        legend_use_petal_color=True,
        style=style_name,
        ax=ax,
    )
    fig.tight_layout()
    fig.savefig(f"examples/figures/venn_style_{style_name}.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: venn_style_{style_name}.png")


# ---------------------------------------------------------------------------
# 4. Global style application
# ---------------------------------------------------------------------------
apply_style("nature")
print("\nGlobal style set to 'nature'.")
fig, ax = plt.subplots(figsize=(10, 3.5), facecolor="w", edgecolor="k")
ax = gv.manhattanplot(data=df, marker=".", ax=ax)
fig.tight_layout()
fig.savefig("examples/figures/manhattan_nature_global.png",
            dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved: manhattan_nature_global.png")

# Restore default
apply_style("geneview")


# ---------------------------------------------------------------------------
# 5. Context-manager (temporary) style
# ---------------------------------------------------------------------------
print("\nUsing context manager for temporary style application...")
with use_style("science"):
    fig, ax = plt.subplots(figsize=(5, 5), facecolor="w", edgecolor="k")
    ax = gv.qqplot(data=df["P"], marker="o", ax=ax)
    fig.tight_layout()
    fig.savefig("examples/figures/qq_science_context.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved: qq_science_context.png")

print("\nDone — all plotstyle examples generated.")
