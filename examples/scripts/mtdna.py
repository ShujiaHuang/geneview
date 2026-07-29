"""Mitochondrial (mtDNA) visualization examples for geneview.

Demonstrates the five purpose-built mtDNA figures, all driven by the kind of
data a MitoQuest run produces (https://github.com/ShujiaHuang/mitoquest):

1. ``mito_genome_map``      -> circular rCRS map with variant lollipops.
2. ``heteroplasmy_scatter`` -> linear position-vs-VAF landscape.
3. ``heteroplasmy_heatmap`` -> samples x variant-sites VAF matrix.
4. ``mito_coverage_plot``   -> sequencing depth across the genome.
5. ``mito_copynumber_plot`` -> per-sample mtDNA copy number with 95% CI.

The script fabricates a small synthetic cohort so it runs with no external
files.  In practice you would replace the synthetic frames with::

    variants = gv.read_mito_vcf("cohort.mt.vcf.gz")      # from mitoquest caller
    copynum  = gv.read_mito_copynumber(glob("*.cn.tsv"))  # from mitoquest copynum
    coverage = gv.read_mito_coverage(["s1.bam", "s2.cram"], reference="rCRS.fa")

Usage::

    python examples/scripts/mtdna.py

Change ``STYLE`` below to switch the journal look for every figure at once.

Author: Shujia Huang
"""
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import geneview as gv

# One-line style switch (None | "geneview" | "nature" | "science" | "cell").
# Note: the "nature" palette starts with black (#000000), so under it the
# protein-coding gene arcs render black.  "cell" gives protein-coding a clear
# blue instead — see docs/mtdna_guide.md for how feature colours are resolved.
STYLE = "cell"
gv.set_style(STYLE)

OUTDIR = "examples/figures"
rng = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# Synthetic MitoQuest-style inputs
# ---------------------------------------------------------------------------
def make_cohort_variants(n_samples=12):
    """Fabricate a tidy long-format variant table like read_mito_vcf()."""
    # A mix of common homoplasmic markers and low-level heteroplasmies,
    # spread across control region, tRNA, rRNA and protein-coding genes.
    sites = [
        # pos,  ref, alt, base_vaf, prevalence
        (263,   "A", "G", 0.99, 0.9),   # D-loop, common homoplasmy
        (750,   "A", "G", 0.99, 0.8),   # MT-RNR1 (12S)
        (1438,  "A", "G", 0.98, 0.7),   # MT-RNR1
        (3243,  "A", "G", 0.15, 0.3),   # MT-TL1, classic MELAS heteroplasmy
        (8344,  "A", "G", 0.10, 0.2),   # MT-TK, MERRF heteroplasmy
        (8860,  "A", "G", 0.99, 0.9),   # MT-ATP6
        (11719, "G", "A", 0.98, 0.6),   # MT-ND4
        (14747, "C", "T", 0.06, 0.25),  # MT-CYB, low-level
        (16189, "T", "C", 0.99, 0.5),   # HVR1
        (16519, "T", "C", 0.99, 0.7),   # HVR1
    ]
    samples = ["S%02d" % (i + 1) for i in range(n_samples)]
    rows = []
    for pos, ref, alt, base_vaf, prev in sites:
        for s in samples:
            if rng.random() > prev:
                continue
            # Homoplasmic sites cluster near the base VAF; low-level ones jitter.
            if base_vaf > 0.9:
                vaf = float(np.clip(rng.normal(0.99, 0.01), 0.85, 1.0))
                status = "HOM"
            else:
                vaf = float(np.clip(base_vaf + rng.normal(0, 0.04), 0.01, 0.6))
                status = "HET"
            rows.append({
                "sample": s, "chrom": "chrM", "pos": pos, "ref": ref,
                "alt": alt, "vaf": round(vaf, 3),
                "depth": int(rng.integers(800, 2500)),
                "gt": "1" if status == "HOM" else "0/1",
                "status": status, "var_type": "SNV", "variant_id": ".",
            })
    return pd.DataFrame(rows)


def make_copynumber(samples):
    """Fabricate per-sample mtDNA copy number with a 95% CI."""
    cn = rng.normal(450, 130, len(samples)).clip(80, None)
    ci = cn * 0.06
    return pd.DataFrame({
        "sample": samples,
        "chrom": "chrM",
        "copy_number": cn.round(1),
        "ci_low": (cn - ci).round(1),
        "ci_high": (cn + ci).round(1),
    })


def make_coverage(samples, bins=1000):
    """Fabricate binned depth curves with a realistic D-loop dip."""
    pos = np.linspace(1, gv.MT_LENGTH, bins)
    frames = []
    for i, s in enumerate(samples):
        base = 600 + i * 120
        depth = base * (1 + 0.15 * np.sin(pos / 900.0))
        # Control-region / origin coverage dip near the ends.
        depth *= np.where((pos < 600) | (pos > 16000), 0.55, 1.0)
        depth += rng.normal(0, base * 0.05, bins)
        frames.append(pd.DataFrame({
            "sample": s, "chrom": "chrM", "pos": pos,
            "depth": depth.clip(0),
        }))
    return pd.concat(frames, ignore_index=True)


variants = make_cohort_variants()
samples = sorted(variants["sample"].unique())
copynum = make_copynumber(samples)
coverage = make_coverage(samples[:4])


# ---------------------------------------------------------------------------
# 1. Circular rCRS genome map  ->  mtdna_genome_map.png
# ---------------------------------------------------------------------------
# The polar projection must be set on the axes, so create it explicitly.
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": "polar"})
gv.mito_genome_map(
    variants, gene_label="large", color_by="feature",
    title="Human mtDNA\n(rCRS, 16,569 bp)", ax=ax,
)
fig.savefig(f"{OUTDIR}/mtdna_genome_map.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/mtdna_genome_map.png")


# ---------------------------------------------------------------------------
# 2. Heteroplasmy landscape (position vs VAF)  ->  mtdna_heteroplasmy.png
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 3.8))
gv.heteroplasmy_scatter(variants, hue="status", het_threshold=0.03, ax=ax)
ax.set_title("Heteroplasmy landscape across the mitochondrial genome")
fig.savefig(f"{OUTDIR}/mtdna_heteroplasmy.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/mtdna_heteroplasmy.png")


# ---------------------------------------------------------------------------
# 3. Cohort heteroplasmy heatmap  ->  mtdna_heatmap.png
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 5))
gv.heteroplasmy_heatmap(variants, site_label="variant", ax=ax)
ax.set_title("Cohort heteroplasmy (VAF) by sample and site")
fig.savefig(f"{OUTDIR}/mtdna_heatmap.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/mtdna_heatmap.png")


# ---------------------------------------------------------------------------
# 4. Sequencing depth across the genome  ->  mtdna_coverage.png
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 3.8))
gv.mito_coverage_plot(coverage, hue="sample", ax=ax)
ax.set_title("mtDNA sequencing depth")
fig.savefig(f"{OUTDIR}/mtdna_coverage.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/mtdna_coverage.png")


# ---------------------------------------------------------------------------
# 5. Per-sample mtDNA copy number  ->  mtdna_copynumber.png
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
gv.mito_copynumber_plot(
    copynum, orient="v", show_ci=True,
    baseline=float(copynum["copy_number"].median()), ax=ax,
)
ax.set_title("mtDNA copy number per sample (95% CI; dashed = cohort median)")
fig.savefig(f"{OUTDIR}/mtdna_copynumber.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTDIR}/mtdna_copynumber.png")

print("\nDone — all mtDNA examples generated.")
