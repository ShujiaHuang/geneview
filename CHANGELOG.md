# Changelog

All notable changes to **geneview** are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
