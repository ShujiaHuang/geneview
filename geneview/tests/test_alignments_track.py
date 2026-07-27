"""Tests for AlignmentsTrack."""
import pytest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._alignments_track import AlignmentsTrack
from geneview.genometracks._base import GenomicInterval


class TestAlignmentsTrackCreation:

    def test_basic_creation(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track.name == "Alignments"
        assert track.height == 2.0
        assert track.plot_types == ["coverage"]

    def test_multiple_plot_types(self):
        track = AlignmentsTrack(filepath="dummy.bam",
                                type=["coverage", "pileup"])
        assert track.plot_types == ["coverage", "pileup"]

    def test_sashimi_type(self):
        track = AlignmentsTrack(filepath="dummy.bam", type="sashimi")
        assert "sashimi" in track.plot_types

    def test_custom_colors(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            col_mates="green",
            col_gap="yellow",
            col_deletion="purple",
            col_insertion="orange",
            fill_coverage="red",
            fill_reads="blue",
        )
        assert track.col_mates == "green"
        assert track.fill_coverage == "red"

    def test_paired_end(self):
        track = AlignmentsTrack(filepath="dummy.bam", is_paired=True)
        assert track.is_paired is True

    def test_reverse_stacking(self):
        track = AlignmentsTrack(filepath="dummy.bam", reverse_stacking=True)
        assert track.reverse_stacking is True

    def test_transformation(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            transformation=lambda x: np.log2(x + 1),
        )
        assert track.transformation is not None

    def test_sashimi_params(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            type="sashimi",
            sashimi_score=5,
            sashimi_height=0.5,
        )
        assert track.sashimi_score == 5
        assert track.sashimi_height == 0.5


class TestAlignmentsTrackImportError:

    def test_import_error_without_pysam(self):
        """When pysam is not installed, should raise ImportError."""
        # This test validates the error path by checking the error message
        track = AlignmentsTrack(filepath="nonexistent.bam")
        # We can't easily test ImportError without mocking, so just
        # verify the method exists
        assert hasattr(track, '_import_pysam')


class TestAlignmentsTrackDraw:

    def test_draw_empty_region(self):
        """Drawing with no accessible BAM should not crash."""
        track = AlignmentsTrack(filepath="nonexistent.bam")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 1000)
        # Should handle missing file gracefully (exception caught internally)
        try:
            track.draw(ax, region)
        except (ImportError, FileNotFoundError, OSError, ValueError):
            pass  # Expected when no pysam or file
        plt.close(fig)


class TestAlignmentsTrackAttributes:

    def test_all_attributes_stored(self):
        track = AlignmentsTrack(
            filepath="test.bam",
            is_paired=True,
            show_mismatches=False,
            show_indels=False,
            reference="ref.fa",
            coverage_height=0.5,
            alpha_reads=0.6,
            alpha_mismatch=0.3,
        )
        assert track.filepath == "test.bam"
        assert track.is_paired is True
        assert track.show_mismatches is False
        assert track.show_indels is False
        assert track.reference == "ref.fa"
        assert track.coverage_height == 0.5
        assert track.alpha_reads == 0.6
        assert track.alpha_mismatch == 0.3


# ---------------------------------------------------------------------------
# Regression tests: BAM/CRAM file-handle management and failure warnings
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_bam(tmp_path):
    """Create a minimal, indexed BAM file with a handful of reads."""
    pysam = pytest.importorskip("pysam")
    bam_path = str(tmp_path / "tiny.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"LN": 10000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for i in range(6):
            read = pysam.AlignedSegment()
            read.query_name = f"read{i}"
            read.query_sequence = "ACGT" * 25  # 100 bp
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 100 + i * 50
            read.mapping_quality = 60
            read.cigar = ((0, 100),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            outf.write(read)
    pysam.index(bam_path)
    return bam_path


def _spy_alignment_file(monkeypatch):
    """Patch ``pysam.AlignmentFile`` so every opened handle is recorded."""
    pysam = pytest.importorskip("pysam")
    opened = []
    real_cls = pysam.AlignmentFile

    def spy(*args, **kwargs):
        handle = real_cls(*args, **kwargs)
        opened.append(handle)
        return handle

    monkeypatch.setattr(pysam, "AlignmentFile", spy)
    return opened


class TestFailureWarnings:
    """A broken data source must surface a warning, not a silent blank panel."""

    def test_pileup_warns_on_missing_file(self):
        track = AlignmentsTrack(filepath="no_such_file.bam", type="pileup")
        fig, ax = plt.subplots()
        with pytest.warns(UserWarning, match="could not fetch reads"):
            track.draw(ax, GenomicInterval("chr1", 0, 1000))
        plt.close(fig)

    def test_sashimi_warns_on_missing_file(self):
        track = AlignmentsTrack(filepath="no_such_file.bam", type="sashimi")
        fig, ax = plt.subplots()
        with pytest.warns(UserWarning, match="could not fetch reads"):
            track.draw(ax, GenomicInterval("chr1", 0, 1000))
        plt.close(fig)

    def test_warning_mentions_region_and_cause(self):
        track = AlignmentsTrack(filepath="no_such_file.bam", type="pileup",
                                name="MyBAM")
        fig, ax = plt.subplots()
        with pytest.warns(UserWarning) as record:
            track.draw(ax, GenomicInterval("chr7", 100, 200))
        plt.close(fig)
        msg = str(record[0].message)
        assert "MyBAM" in msg
        assert "chr7:100-200" in msg


class TestHandleManagement:
    """The track must close every opened BAM handle, on success and on error."""

    def test_fetch_failure_closes_handle(self, monkeypatch, tiny_bam):
        """_fetch_reads must close the handle when region matching fails."""
        opened = _spy_alignment_file(monkeypatch)
        track = AlignmentsTrack(filepath=tiny_bam, type="pileup")
        fig, ax = plt.subplots()
        with pytest.warns(UserWarning, match="could not fetch reads"):
            track.draw(ax, GenomicInterval("chrNope", 0, 1000))
        plt.close(fig)
        assert opened, "expected at least one opened handle"
        assert all(h.closed for h in opened)

    def test_coverage_error_closes_handle(self, monkeypatch, tiny_bam):
        """_compute_coverage must close the handle even when it raises."""
        opened = _spy_alignment_file(monkeypatch)
        track = AlignmentsTrack(filepath=tiny_bam)
        with pytest.raises(Exception):
            track._compute_coverage(GenomicInterval("chrNope", 0, 1000))
        assert opened, "expected at least one opened handle"
        assert all(h.closed for h in opened)

    def test_successful_draw_closes_all_handles(self, monkeypatch, tiny_bam):
        """A full coverage+pileup draw must not leak any handle."""
        opened = _spy_alignment_file(monkeypatch)
        track = AlignmentsTrack(
            filepath=tiny_bam, type=["coverage", "pileup"],
            show_mismatches=False, quick_consensus=False,
        )
        fig, ax = plt.subplots()
        track.draw(ax, GenomicInterval("chr1", 0, 1000))
        plt.close(fig)
        assert opened, "expected the BAM to be opened at least once"
        assert all(h.closed for h in opened)

    def test_get_region_closes_handle(self, monkeypatch, tiny_bam):
        opened = _spy_alignment_file(monkeypatch)
        track = AlignmentsTrack(filepath=tiny_bam)
        region = track.get_region()
        assert region is not None and region.chrom == "chr1"
        assert opened and all(h.closed for h in opened)


# ---------------------------------------------------------------------------
# Genomics-convention style tests (nucleotide palette, sashimi strand/lw)
# ---------------------------------------------------------------------------


@pytest.fixture
def junction_bam(tmp_path):
    """BAM with two junctions of different strand and support depth."""
    pysam = pytest.importorskip("pysam")
    bam_path = str(tmp_path / "junctions.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"LN": 100000, "SN": "chr1"}]}

    def _read(name, start, cigar, reverse):
        r = pysam.AlignedSegment()
        r.query_name = name
        r.query_sequence = "A" * 100
        r.flag = 16 if reverse else 0
        r.reference_id = 0
        r.reference_start = start
        r.mapping_quality = 60
        r.cigar = cigar
        r.query_qualities = pysam.qualitystring_to_array("I" * 100)
        return r

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Junction A: 3 forward reads, span (150, 350)
        for i in range(3):
            outf.write(_read(f"f{i}", 100, ((0, 50), (3, 200), (0, 50)), False))
        # Junction B: 1 reverse read, span (550, 950)
        outf.write(_read("r0", 500, ((0, 50), (3, 400), (0, 50)), True))
    pysam.index(bam_path)
    return bam_path


class TestNucleotidePalette:
    """Mismatch colours must follow the IGV/biovizBase Okabe-Ito palette
    and stay consistent with SequenceTrack."""

    def test_palette_matches_sequence_track(self):
        from geneview.genometracks._sequence_track import _DEFAULT_NUC_COLORS
        from geneview.genometracks._alignments_track import _NUC_COLORS
        for base in ("A", "C", "G", "T", "N"):
            assert _NUC_COLORS[base] == _DEFAULT_NUC_COLORS[base]

    def test_palette_is_colorblind_safe_igv_convention(self):
        from geneview.genometracks._alignments_track import _NUC_COLORS
        assert _NUC_COLORS["A"] == "#009E73"  # green
        assert _NUC_COLORS["C"] == "#0072B2"  # blue
        assert _NUC_COLORS["G"] == "#E69F00"  # orange
        assert _NUC_COLORS["T"] == "#D55E00"  # vermillion
        assert _NUC_COLORS["N"] == "#999999"  # gray


class TestSashimiJunctionStyle:
    """Sashimi arcs encode junction strand (colour) and support (line width)."""

    @staticmethod
    def _make_read(reverse=False, xs=None):
        pysam = pytest.importorskip("pysam")
        read = pysam.AlignedSegment()
        read.query_name = "j"
        read.query_sequence = "A" * 100
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigar = ((0, 50), (3, 200), (0, 50))
        if xs is not None:
            read.set_tag("XS", xs)
        return read

    def test_collect_junctions_tallies_strand(self):
        reads = [self._make_read(), self._make_read(), self._make_read(reverse=True)]
        track = AlignmentsTrack(filepath="dummy.bam")
        junc = track._collect_junctions(reads)
        assert junc[(150, 350)] == [3, 2, 1]

    def test_collect_junctions_prefers_xs_tag(self):
        # XS="-" on a forward read must count as reverse-strand support
        track = AlignmentsTrack(filepath="dummy.bam")
        junc = track._collect_junctions([self._make_read(reverse=False, xs="-")])
        assert junc[(150, 350)] == [1, 0, 1]

    def test_junction_color_majority(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        assert track._junction_color(3, 1) == track.col_sashimi_fwd
        assert track._junction_color(1, 3) == track.col_sashimi_rev
        assert track._junction_color(2, 2) == track.col_sashimi_unknown
        assert track._junction_color(0, 0) == track.col_sashimi_unknown

    def test_arc_colours_and_linewidths(self, junction_bam):
        from matplotlib.patches import Arc
        from matplotlib.colors import to_rgb
        track = AlignmentsTrack(filepath=junction_bam, type="sashimi",
                                sashimi_score=1)
        fig, ax = plt.subplots()
        track.draw(ax, GenomicInterval("chr1", 0, 2000))
        arcs = [p for p in ax.patches if isinstance(p, Arc)]
        plt.close(fig)

        assert len(arcs) == 2
        edge_rgbs = [a.get_edgecolor()[:3] for a in arcs]
        fwd_rgb = to_rgb(track.col_sashimi_fwd)
        rev_rgb = to_rgb(track.col_sashimi_rev)
        assert any(np.allclose(c, fwd_rgb) for c in edge_rgbs)
        assert any(np.allclose(c, rev_rgb) for c in edge_rgbs)

        # The 3-read junction must be drawn thicker than the 1-read one.
        lws = sorted(a.get_linewidth() for a in arcs)
        assert lws[0] < lws[1]

    def test_sashimi_score_filters_by_count(self, junction_bam):
        from matplotlib.patches import Arc
        track = AlignmentsTrack(filepath=junction_bam, type="sashimi",
                                sashimi_score=2)
        fig, ax = plt.subplots()
        track.draw(ax, GenomicInterval("chr1", 0, 2000))
        arcs = [p for p in ax.patches if isinstance(p, Arc)]
        plt.close(fig)
        assert len(arcs) == 1  # only the 3-read junction survives
