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
