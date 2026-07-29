# Changelog

All notable changes to **geneview** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Mitochondrial (mtDNA) visualization module (`geneview.mtdna`):** a
  purpose-built toolkit for human mtDNA analysis, driven by
  [MitoQuest](https://github.com/ShujiaHuang/mitoquest)-style inputs
  (single/multi-sample VCF with per-sample heteroplasmy fractions,
  `copynum` TSV, and BAM/CRAM). Provides:
  - **Readers** `read_mito_vcf`, `read_mito_copynumber`, `read_mito_coverage`
    that normalise those inputs into tidy `DataFrame`s (the VCF reader auto-
    detects the heteroplasmy field across MitoQuest versions, `AF`/`HF`).
  - **Plots** `mito_genome_map` (circular rCRS map with variant lollipops),
    `heteroplasmy_scatter` (position-vs-VAF landscape), `heteroplasmy_heatmap`
    (samples x sites VAF matrix), `mito_coverage_plot` (sequencing depth) and
    `mito_copynumber_plot` (per-sample copy number with 95% CI).
  - **Reference backbone** `get_mt_genes`, `gene_at`, `genes_in_range`,
    `is_mt_contig`, `MT_LENGTH` — the rCRS (NC_012920.1) 16,569 bp, 37-gene
    map plus the D-loop control region.
  - **CLI** `geneview mito-map`, `mito-heteroplasmy`, `mito-heatmap`,
    `mito-coverage`, `mito-copynumber`. See `docs/mtdna_guide.md`.

- **Area-proportional Venn diagrams:** `venn(..., proportional=True)` draws
  2- or 3-set diagrams whose circle areas and overlaps scale with the real
  set / intersection sizes (radii and centre distances are solved so the drawn
  areas match the counts). For other set counts the flag is ignored with a
  warning and the schematic layout is used. Exposed on the CLI via
  `geneview venn --proportional`.

### Changed
- **Venn colour palette follows the active plot style:** `venn()`'s `palette`
  now defaults to `None`, in which case colours are taken from the active
  `PlotStyle` (e.g. `nature`/`science`/`cell`), falling back to a built-in
  curated palette. An explicit `palette` still takes precedence. The
  `geneview venn --palette` default changed from `viridis` to the active
  `--style` palette accordingly.

## [0.8.2] - 2026-07-27

### Fixed
- **CRAM pileup warning (`_mismatch_counts.tally_reads`):** pysam's pileup
  iterator defaults `multiple_iterators=True`, which is unsupported for CRAM
  and triggered `UserWarning: multiple_iterators not implemented for CRAM`.
  The tally now opens a dedicated handle and requests
  `multiple_iterators=False`, which is correct and warning-free for both BAM
  and CRAM.

[0.8.2]: https://github.com/ShujiaHuang/geneview/compare/v0.8.1...v0.8.2

## [0.8.1] - 2026-07-27

### Fixed
- **CRAM reference handling (`geneview tracks` and the track API):** BAM/CRAM
  files are now opened through a single `open_alignment_file()` helper that
  auto-detects the pysam mode (`.cram` -> `"rc"`, otherwise `"rb"`) and, when
  a `--reference` FASTA is supplied, passes it to pysam as `reference_filename`.
  CRAM input is decoded against the user-supplied reference instead of the
  (often stale, machine-specific) path embedded in the CRAM header. The
  reference is now threaded through `AlignmentsTrack`, `BAMCoverageTrack`,
  `GroupedAlignmentsTrack`, `AnnotationTrack.from_bam`, and the
  `is_paired_end` / `is_long_frag_dataset` helpers.

[0.8.1]: https://github.com/ShujiaHuang/geneview/compare/v0.8.0...v0.8.1

## [0.8.0] - 2026-07-27

Highlights: a full mutation-visualization suite (Lollipop/Dandelion tracks),
richly styled UCSC-grade gene models, a more robust and expert-styled alignment
track, and a hardened command-line interface.

### Added
- **Mutation visualization (trackViewer parity):** new `LolliplotTrack` and
  `DandelionTrack` with Tanghulu-style stacking, node de-overlap, and automatic
  contrast-aware labels; shared mutation feature/shape primitives.
- **GeneRegionTrack:** four drawing styles (`UCSC`, `flybase`, `tssarrow`,
  `exonarrows`), UCSC-style backbone line, stepped exon polygons, intron chevron
  arrows, read-direction arrows, and a configurable `arrow_length`.
- **Genome source abstraction:** new `GenomeSource` / `FastaGenomeSource`
  classes for lazy reference access with `close()` / context-manager support.
- **CLI — global `--debug`:** re-raises the full traceback on error; otherwise
  prints a concise `[ERROR]` plus a hint to re-run with `--debug`.
- **CLI — Manhattan annotation:** `--annotate-fmt` label templating over
  `{snp, chrom, pos, p, log10p}` and `--annotate-layout {repel,lane}` placement
  with rotation/color options; `annotate_fmt` is validated at call time so
  unknown fields / bad format specs fail fast with a clear message.
- **Styling:** journal palettes and style presets (`nature`, `science`, `cell`),
  panel labels, and multi-view layout support.

### Changed
- **AlignmentsTrack:** sashimi arcs are strand-colored (prefers the `XS` tag)
  with line width scaled by √(read support), matching the IGV/Gviz convention;
  pileup nucleotide colors unified with `SequenceTrack` (Okabe-Ito,
  colorblind-safe); added `color_fn` support.
- **Karyotype/Ideogram:** unified cytoband color loading; each chromosome is
  outlined so white `gneg` bands remain visible.
- Vector export improvements.
- `release.py` modernized: builds sdist + wheel via `python -m build` and
  validates with `twine check`; publishing is a separate explicit step.

### Fixed
- **AlignmentsTrack / grouped alignments:** BAM/CRAM file handles are wrapped in
  `try/finally` and never leak on failure; pileup/sashimi draw failures now warn
  and leave a blank panel instead of silently swallowing the exception.
- **DataTrack:** heatmap/gradient cells are drawn at true genomic coordinates.

### Docs
- New CLI reference (`docs/cli.md`), Manhattan annotation demo + figures,
  karyoplot / palette-comparison / node-deoverlap examples, refreshed tutorial
  notebooks and READMEs, and a new geneview logo.

[0.8.0]: https://github.com/ShujiaHuang/geneview/compare/v0.7.0...v0.8.0
