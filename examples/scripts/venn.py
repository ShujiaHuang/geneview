"""Venn diagram examples for geneview.

Demonstrates:

1. Schematic Venn diagrams for 2-6 sets (ellipses for 2-5, triangles for 6).
2. Area-proportional Venn diagrams for 2-3 sets (circle areas and overlaps
   match the real set / intersection sizes) via ``proportional=True``.
3. Journal style integration: with ``palette=None`` the petal colours follow
   the active plot style (``geneview``/``nature``/``science``/``cell``).
4. Manual petal-label adjustment through ``generate_petal_labels``.

Usage::

    python examples/scripts/venn.py

Change ``STYLE`` below to switch the journal look for every figure at once.
"""
from itertools import islice
from string import ascii_uppercase

from numpy.random import default_rng

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geneview as gv

# One-line style switch (None | "geneview" | "nature" | "science" | "cell").
STYLE = "nature"
gv.set_style(STYLE)

OUTDIR = "examples/figures"
rng = default_rng(42)


# ---------------------------------------------------------------------------
# 1. Schematic Venn diagrams for 2-6 sets  ->  venn.png
# ---------------------------------------------------------------------------
fig_all = plt.figure(figsize=(18, 13))
combined_axs = [fig_all.add_subplot(2, 3, i + 1) for i in range(5)]
letters = iter(ascii_uppercase)
for n_sets, ax in zip(range(2, 7), combined_axs):
    dataset_dict = {
        name: set(rng.choice(1000, 700, replace=False))
        for name in islice(letters, n_sets)
    }
    gv.venn(dataset_dict,
            fmt="{percentage:.1f}%",   # "{size}", "{logic}"
            fontsize=12,
            legend_use_petal_color=True,
            ax=ax)
    ax.set_title(f"{n_sets} sets", fontsize=13)

fig_all.tight_layout()
fig_all.savefig(f"{OUTDIR}/venn.png", dpi=150, bbox_inches="tight")
plt.close(fig_all)
print(f"  Saved: {OUTDIR}/venn.png")


# ---------------------------------------------------------------------------
# 2. Minimal 3-set example  ->  venn3.png
# ---------------------------------------------------------------------------
table = {
    "Dataset 1": {"A", "B", "D", "E"},
    "Dataset 2": {"C", "F", "B", "G"},
    "Dataset 3": {"J", "C", "K"},
}
fig, ax = plt.subplots(figsize=(6, 6))
gv.venn(table, legend_use_petal_color=True, ax=ax)
fig.tight_layout()
fig.savefig(f"{OUTDIR}/venn3.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/venn3.png")


# ---------------------------------------------------------------------------
# 3. Manual adjustment of petal labels  ->  venn4.png
# ---------------------------------------------------------------------------
dataset_dict = {
    name: set(rng.choice(1000, 250, replace=False))
    for name in list("ABCD")
}
petal_labels = gv.generate_petal_labels(
    dataset_dict.values(), fmt="{logic}\n({percentage:.1f}%)")
fig, ax = plt.subplots(figsize=(6, 6))
gv.venn(data=petal_labels, names=list(dataset_dict.keys()),
        legend_use_petal_color=True, ax=ax)
fig.tight_layout()
fig.savefig(f"{OUTDIR}/venn4.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/venn4.png")


# ---------------------------------------------------------------------------
# 4. Area-proportional Venn (2-3 sets)  ->  venn_proportional.png
# ---------------------------------------------------------------------------
prop2 = {"Large study": set(range(0, 200)), "Small study": set(range(180, 235))}
prop3 = {
    "Cohort A": set(range(0, 120)),
    "Cohort B": set(range(80, 160)),
    "Cohort C": set(range(60, 110)),
}

fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18, 6))
gv.venn(prop2, proportional=True, legend_use_petal_color=True, ax=axs[0])
axs[0].set_title("Proportional (2 sets)", fontsize=13)
gv.venn(prop3, proportional=True, legend_use_petal_color=True, ax=axs[1])
axs[1].set_title("Proportional (3 sets)", fontsize=13)
gv.venn(prop3, legend_use_petal_color=True, ax=axs[2])
axs[2].set_title("Schematic (3 sets, for comparison)", fontsize=13)
fig.tight_layout()
fig.savefig(f"{OUTDIR}/venn_proportional.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/venn_proportional.png")

print("\nDone — all venn examples generated.")
