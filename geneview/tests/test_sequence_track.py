"""Tests for SequenceTrack."""
import pytest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._sequence_track import SequenceTrack, _reverse_complement
from geneview.genometracks._base import GenomicInterval
from geneview.genometracks._track_plot import plot_tracks


class TestSequenceTrackHelpers:

    def test_reverse_complement(self):
        assert _reverse_complement("ACGT") == "ACGT"  # palindrome
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("GCGC") == "GCGC"
        assert _reverse_complement("ACGTacgt") == "acgtACGT"

    def test_reverse_complement_empty(self):
        assert _reverse_complement("") == ""


class TestSequenceTrackCreation:

    def test_basic_string_sequence(self):
        track = SequenceTrack(sequence="ACGTACGT")
        assert track.name == "Sequence"
        assert track.height == 0.5

    def test_custom_name_height(self):
        track = SequenceTrack(sequence="ACGT", name="MySeq", height=1.0)
        assert track.name == "MySeq"
        assert track.height == 1.0

    def test_default_nuc_colors(self):
        track = SequenceTrack(sequence="ACGT")
        assert "A" in track._nuc_colors
        assert "C" in track._nuc_colors
        assert "G" in track._nuc_colors
        assert "T" in track._nuc_colors

    def test_custom_fontcolor(self):
        track = SequenceTrack(sequence="ACGT", fontcolor={"A": "red"})
        assert track._nuc_colors["A"] == "red"


class TestSequenceTrackDraw:

    def test_draw_letters_small_window(self):
        """At small window (< 200bp), letters should be drawn."""
        seq = "ACGTACGT" * 10
        track = SequenceTrack(sequence=seq)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, len(seq))
        track.draw(ax, region)
        assert ax.get_xlim() == (0, len(seq))
        plt.close(fig)

    def test_draw_boxes_medium_window(self):
        """At medium window (200-2000bp), colored boxes should be drawn."""
        seq = "ACGTACGT" * 100
        track = SequenceTrack(sequence=seq)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, len(seq))
        track.draw(ax, region)
        assert ax.get_xlim() == (0, len(seq))
        plt.close(fig)

    def test_draw_line_large_window(self):
        """At large window (> 2000bp), a single colored line."""
        seq = "ACGTACGT" * 500
        track = SequenceTrack(sequence=seq)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, len(seq))
        track.draw(ax, region)
        assert ax.get_xlim() == (0, len(seq))
        plt.close(fig)

    def test_draw_noLetters(self):
        """noLetters=True forces boxes even at small window."""
        seq = "ACGTACGT"
        track = SequenceTrack(sequence=seq, noLetters=True)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, len(seq))
        track.draw(ax, region)
        plt.close(fig)


class TestSequenceTrackComplement:

    def test_complement_mode(self):
        seq = "AAAA"
        track = SequenceTrack(sequence=seq, complement=True)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 4)
        track.draw(ax, region)
        plt.close(fig)


class TestSequenceTrackAdd53:

    def test_add53_arrow(self):
        seq = "ACGTACGT" * 10
        track = SequenceTrack(sequence=seq, add53=True)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, len(seq))
        track.draw(ax, region)
        plt.close(fig)


class TestSequenceTrackIntegration:

    def test_plot_tracks_integration(self):
        seq = "ACGTACGT" * 20
        track = SequenceTrack(sequence=seq)
        region = GenomicInterval("chr1", 0, len(seq))
        axes = plot_tracks([track], region=region, figsize=(8, 2))
        assert len(axes) == 1
        plt.close("all")

    def test_empty_sequence(self):
        track = SequenceTrack(sequence="")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 100)
        track.draw(ax, region)
        plt.close(fig)
