# Mitochondrial DNA (mtDNA) Guide

The `geneview.mtdna` module is a purpose-built toolkit for the figures that
matter most in human mitochondrial genome analysis. It turns the routine
outputs of an mtDNA pipeline — especially
[MitoQuest](https://github.com/ShujiaHuang/mitoquest) — into a small set of
publication-ready plots, with a matching command-line interface.

Everything is anchored to the **revised Cambridge Reference Sequence** (rCRS,
GenBank `NC_012920.1`), the 16,569 bp circular reference that MITOMAP,
HelixMTdb, gnomAD (chrM) and MitoQuest are all coordinated to. Positions are
**1-based, inclusive**, matching the VCF and mtDNA-literature convention.

---

## Table of Contents

1. [Why these figures?](#why-these-figures)
2. [Inputs and readers](#inputs-and-readers)
3. [The reference backbone](#the-reference-backbone)
4. [Plots](#plots)
   - [mito_genome_map](#mito_genome_map)
   - [heteroplasmy_scatter](#heteroplasmy_scatter)
   - [heteroplasmy_heatmap](#heteroplasmy_heatmap)
   - [mito_coverage_plot](#mito_coverage_plot)
   - [mito_copynumber_plot](#mito_copynumber_plot)
5. [Journal styles](#journal-styles)
6. [Command-Line Interface](#command-line-interface)
7. [End-to-end example](#end-to-end-example)
8. [Tips and gotchas](#tips-and-gotchas)

---

## Why these figures?

Across the human mtDNA literature (disease diagnostics, ageing, cancer,
forensics and population genetics) a handful of graphical questions recur.
The module maps each one to a single function:

| Question a paper needs to answer | Plot |
|----------------------------------|------|
| *Where* on the genome do the variants fall, and *how heteroplasmic* are they? | `mito_genome_map` (circular rCRS map) |
| What is the heteroplasmy (VAF) *landscape* along the genome? Which variants clear the calling threshold? | `heteroplasmy_scatter` |
| Which variants are *shared vs private* across a cohort, and at what heteroplasmy level? | `heteroplasmy_heatmap` |
| Is coverage even? Are there NUMT / reference artefacts or drop-outs? | `mito_coverage_plot` |
| How does mtDNA *copy number* vary across samples/tissues? | `mito_copynumber_plot` |

The single most important mtDNA measurement is the **heteroplasmy fraction**
(VAF): the fraction of mtDNA molecules in a sample that carry a given allele.
Because a cell holds hundreds–thousands of mtDNA copies, an allele can sit
anywhere between 0 and 1, and that fraction is often the biologically and
clinically decisive number. Three of the five plots put VAF front and centre.

---

## Inputs and readers

All readers live in `geneview.mtdna` (also re-exported at the top level as
`geneview.read_mito_*`). Each returns a tidy `pandas.DataFrame`.

### `read_mito_vcf` — single/multi-sample VCF

```python
df = gv.read_mito_vcf(
    "cohort.mt.vcf.gz",
    samples=None,      # restrict to a list of sample names
    region=None,       # "chrM" or "chrM:1-576" (needs a .tbi/.csi index)
    min_vaf=0.0,       # drop rows below this heteroplasmy fraction
    include_ref=False, # also emit rows for reference alleles
    mt_only=True,      # keep only mitochondrial contigs
)
```

The output has one row per **sample × site × ALT allele**:

| Column | Meaning |
|--------|---------|
| `sample` | sample name from the VCF header |
| `chrom`, `pos` | contig and 1-based position |
| `ref`, `alt` | reference and this ALT allele |
| `vaf` | heteroplasmy fraction for this ALT allele |
| `depth` | read depth at the site |
| `gt` | genotype string (mtDNA is haploid; heteroplasmy is encoded as multiple alleles, e.g. `0/1`, `1/2`) |
| `status` | `HET` (heteroplasmic — more than one distinct non-missing allele) or `HOM` (homoplasmic) |
| `var_type` | `SNV` / `MNV` / `INS` / `DEL` |
| `variant_id` | VCF ID column (e.g. an rsID) |

**MitoQuest field compatibility.** `mitoquest caller` writes the heteroplasmy
fraction into a per-sample `FORMAT` field that has been named `AF` (current
builds) or `HF` (older builds), always one value per ALT allele. The reader
probes a candidate list (`AF`, `HF`, `VAF`, `FREQ`) and uses whichever is
present, so both layouts parse with no configuration. Multi-allelic calls
(`1/2`) expand into one row per ALT with the corresponding per-ALT VAF.

### `read_mito_copynumber` — `mitoquest copynum` TSV

```python
cn = gv.read_mito_copynumber(["s1.cn.tsv", "s2.cn.tsv"],
                             sample_names=None)  # default: from file names
```

Keeps the mitochondrial row from each TSV and returns
`sample, chrom, copy_number, ci_low, ci_high`, sorted by copy number. The
reader normalises the `#Chromosome …` header and matches the `CopyNum` /
`CopyNum-CI95-Lower` / `CopyNum-CI95-Upper` columns case-insensitively.

### `read_mito_coverage` — BAM/CRAM

```python
cov = gv.read_mito_coverage(
    ["s1.bam", "s2.cram"],
    sample_names=None,   # default: from file names
    bins=1000,           # bins spanning the contig
    reference="rCRS.fa", # required to decode CRAM
    contig=None,         # force a contig name; default auto-detect
)
```

Auto-detects the mitochondrial contig in each file's header (`chrM`, `MT`,
`chrMT`, `M`, `NC_012920.1`, …), so mixed BAM/CRAM, single or multi sample,
all reduce to one long-format table `sample, chrom, start, end, pos, depth`
(`pos` = bin midpoint). This wraps geneview's shared alignment reader, so it
behaves identically to the genome-tracks coverage machinery.

---

## The reference backbone

The rCRS gene map (13 protein-coding, 22 tRNA, 2 rRNA genes + the D-loop
control region) is available directly:

```python
gv.get_mt_genes()                 # list of dicts (name, start, end, strand, feature_type)
gv.get_mt_genes(as_dataframe=True)
gv.gene_at(3308)                  # -> {"name": "MT-ND1", ...}
gv.mtdna.genes_in_range(577, 1601)
gv.mtdna.is_mt_contig("chrM")     # True
gv.MT_LENGTH                      # 16569
```

The D-loop spans the origin, so it is represented as two arcs (`1..576` and
`16024..16569`); downstream code never has to special-case the wrap-around.
The three hypervariable segments are exposed as
`gv.mtdna.MT_HYPERVARIABLE_REGIONS`.

Feature-type colours default to a colour-blind-friendly palette
(`gv.mtdna.MT_FEATURE_COLORS`) but automatically follow the active journal
style — see [Journal styles](#journal-styles).

---

## Plots

Every plotting function is decorated with geneview's `styled_plot`, so it:

- accepts an optional `ax=` (created automatically when omitted),
- accepts an optional `style=` (`"geneview"`/`"nature"`/`"science"`/`"cell"`
  or a `PlotStyle`), applied for the duration of the call,
- returns the `matplotlib.axes.Axes` it drew on.

### `mito_genome_map`

The flagship figure: a circular rCRS map with the gene ring coloured by
feature type and variants drawn as inward **lollipops** whose length encodes
VAF.

```python
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": "polar"})
gv.mito_genome_map(
    variants,                 # tidy frame from read_mito_vcf (or None for a bare ring)
    gene_label="large",       # "large" | "all" | "none"
    color_by="feature",       # "feature" | "status" | "var_type"
    feature_colors=None,      # dict override, e.g. {"tRNA": "#999"}
    tick_interval=2000,       # bp between position ticks (0 hides)
    title="Human mtDNA",      # text in the ring centre
    ax=ax,
)
```

> **Polar axes.** This plot needs a polar projection. When you pass your own
> `ax`, create it with `subplot_kw={"projection": "polar"}`. When `ax` is
> omitted, `mito_genome_map` builds the polar axes itself.

Key parameters: `gene_ring=(inner, outer)` and `lollipop_base` /
`lollipop_min` control the ring and stem geometry (radii in axes units 0..1;
stems point inward so a *smaller* tip radius = a *longer*, higher-VAF stem).

### `heteroplasmy_scatter`

The linear counterpart: heteroplasmy fraction (y) against genome position (x),
with a gene strip below the axis and a dashed heteroplasmy-threshold line.

```python
gv.heteroplasmy_scatter(
    variants,
    hue="status",         # "feature" | "status" | "var_type" | "sample" | None
    het_threshold=0.03,   # dashed line at the calling threshold; None hides it
    show_gene_band=True,  # coloured gene strip beneath the axis
    marker_size=26,
)
```

Colour by `status` to separate heteroplasmic from homoplasmic calls, by
`sample` for a small cohort overlay, or by `feature` to see which gene classes
carry variation.

### `heteroplasmy_heatmap`

A cohort matrix: rows = samples, columns = variant sites, cell colour = VAF.
Ideal for spotting shared vs private and low-level heteroplasmy.

```python
gv.heteroplasmy_heatmap(
    variants,
    cmap="rocket_r",       # falls back to "Reds" if unavailable
    site_label="variant",  # "pos" | "variant" (pos ref>alt)
    vmin=0.0, vmax=1.0,
    colorbar=True,
)
```

Missing (uncalled) sample×site cells are drawn in light grey, so absence is
visually distinct from a VAF of 0.

### `mito_coverage_plot`

Sequencing depth across the genome — the QC backbone of any call set, and the
view that exposes uneven coverage, NUMT artefacts and origin drop-outs.

```python
cov = gv.read_mito_coverage(["s1.bam", "s2.bam"])
gv.mito_coverage_plot(
    cov,
    hue="sample",          # per-sample lines; None for a single colour
    log=False,             # log depth axis (disables the gene band)
    fill=True,             # fill under a single-sample curve
    show_gene_band=True,
)
```

### `mito_copynumber_plot`

Per-sample mtDNA copy number with 95% confidence-interval whiskers, straight
from `mitoquest copynum`.

```python
cn = gv.read_mito_copynumber(glob.glob("*.cn.tsv"))
gv.mito_copynumber_plot(
    cn,
    orient="v",            # "v" (vertical) | "h" (better for many samples)
    show_ci=True,          # 95% CI error bars
    baseline=cn["copy_number"].median(),  # reference line; None hides it
    sort=True,
)
```

---

## Journal styles

Colours follow the active plot style. Set it globally, per call, or with a
context manager:

```python
gv.set_style("nature")                      # global
gv.mito_genome_map(variants, style="cell")  # per call
with gv.use_style("science"):
    gv.heteroplasmy_scatter(variants)
```

`resolve_feature_colors` maps the first four palette colours of the active
style onto protein-coding / rRNA / tRNA / control-region, so the gene ring and
gene strips stay consistent with the rest of a figure. An explicit
`feature_colors=` dict always wins.

### Feature colours: resolution order and how to change them

The four feature-type colours (protein-coding, rRNA, tRNA, control region) are
resolved with a three-tier priority — highest wins:

1. **An explicit `feature_colors=` dict** passed to the plotting function.
   Partial dicts are merged over the lower tiers, so you can recolour just one
   type.
2. **The active plot style's palette.** When a style is set (via
   `gv.set_style(...)`, `style=` or `gv.use_style(...)`), its first four
   colours are mapped in the fixed order
   `protein_coding → rRNA → tRNA → control_region`.
3. **The built-in defaults** in
   `geneview.mtdna._reference.MT_FEATURE_COLORS`
   (protein-coding `#4C72B0` blue, rRNA `#55A868` green, tRNA `#DD8452`
   orange, control region `#C44E52` red). These apply only when **no** style is active.

> **Why are protein-coding genes black under `nature`?**
> The `nature` style uses the colour-blind-safe **Okabe–Ito** palette, whose
> *first* colour is pure black (`#000000`). Because tier 2 maps that first
> colour onto protein-coding (the first feature type), the 13 protein-coding
> arcs render black under `nature`. This is expected, not a bug — it is simply
> the palette's ordering. `science`, `cell` and `geneview` start with a
> coloured hue, so protein-coding is coloured under them (the bundled example
> `examples/scripts/mtdna.py` uses `cell`, giving protein-coding a blue arc).

To change the colours, pick whichever tier fits your need:

```python
# (a) Recolour one or more feature types for a single call — highest priority.
gv.mito_genome_map(variants, feature_colors={"protein_coding": "#4C72B0"})

# (b) Keep the nature style but give protein-coding a colour instead of black.
with gv.use_style("nature"):
    gv.mito_genome_map(variants, feature_colors={"protein_coding": "#3B5488"})

# (c) Provide the full mapping (order does not matter here).
colors = {
    "protein_coding": "#0072B2",
    "rRNA": "#D55E00",
    "tRNA": "#009E73",
    "control_region": "#CC79A7",
}
gv.mito_genome_map(variants, feature_colors=colors)

# (d) Switch to a style whose palette does not start with black.
gv.set_style("cell")   # or "science" / "geneview"
```

`feature_colors=` is accepted by `mito_genome_map` and `heteroplasmy_scatter`
(the plots that draw the feature-coloured gene ring/strip). The
`heteroplasmy_heatmap` colours cells by VAF through its `cmap=` argument
(default `rocket_r`) instead. To change the no-style fallback globally, edit
`MT_FEATURE_COLORS` in `geneview/mtdna/_reference.py`.

---

## Command-Line Interface

Five `mito-*` subcommands mirror the Python API and share the
[common figure options](./cli.md#common-options-all-subcommands)
(`-o/--output`, `--figsize`, `--dpi`, `--facecolor`, `--style`).

```bash
# Circular rCRS map with variant lollipops
geneview mito-map -i cohort.mt.vcf.gz -o mito_map.png \
    --gene-label large --color-by feature --title "Cohort mtDNA"

# Position-vs-VAF heteroplasmy landscape
geneview mito-heteroplasmy -i cohort.mt.vcf.gz -o het.png \
    --hue status --het-threshold 0.03

# Cohort samples x sites heatmap
geneview mito-heatmap -i cohort.mt.vcf.gz -o heatmap.png --site-label variant

# Sequencing depth from BAM/CRAM (multi-sample)
geneview mito-coverage -i s1.bam s2.cram -o cov.png \
    --reference rCRS.fa --bins 1000 --log

# Per-sample copy number with 95% CI
geneview mito-copynumber -i *.cn.tsv -o cn.png --orient h --baseline 250
```

VCF-driven commands (`mito-map`, `mito-heteroplasmy`, `mito-heatmap`) share
`--region`, `--samples` and `--min-vaf`. Run `geneview <subcommand> --help`
for the full option list.

---

## End-to-end example

A runnable script that fabricates a synthetic MitoQuest-style cohort and
produces all five figures ships with the repository:

```bash
python examples/scripts/mtdna.py
# -> examples/figures/mtdna_genome_map.png, mtdna_heteroplasmy.png,
#    mtdna_heatmap.png, mtdna_coverage.png, mtdna_copynumber.png
```

Change the single `STYLE = "cell"` line at the top to re-theme every figure
at once, and swap the `make_*` helpers for the real `gv.read_mito_*` readers
to run it on your own data.

---

## Tips and gotchas

- **Polar axes for the genome map.** Only `mito_genome_map` uses a polar
  projection; the other four are ordinary Cartesian axes. Passing a non-polar
  `ax` to `mito_genome_map` will fail — let it create its own, or build the
  axes with `subplot_kw={"projection": "polar"}`.
- **CRAM needs a reference.** `read_mito_coverage` requires `reference=` (a
  FASTA) to decode CRAM; BAM does not.
- **Contig naming.** The readers recognise `chrM`, `MT`, `chrMT`, `M` and
  `NC_012920.1` interchangeably. If your BAM uses an unusual name, pass
  `contig=` to `read_mito_coverage` explicitly.
- **Heteroplasmy threshold.** `het_threshold` is cosmetic (a guide line); the
  actual low-VAF filtering happens in the reader via `min_vaf`, mirroring
  `mitoquest caller`'s calling threshold.
- **Large cohorts.** For the heatmap and copy-number bar chart, prefer
  `orient="h"` / a taller figure when the sample count is high.
