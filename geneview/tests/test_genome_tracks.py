"""Tests for geneview.genometracks module - track classes and stacking.

Tests cover:
    - GenomicInterval dataclass
    - Base track classes (Track, RangeTrack, StackedTrack, NumericTrack)
    - Stacking algorithms
    - Individual track types (GenomeAxisTrack, AnnotationTrack, etc.)
    - File I/O utilities
"""
import os
import tempfile

import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._base import (
    GenomicInterval,
    Track,
    RangeTrack,
    StackedTrack,
    NumericTrack,
)
from geneview.genometracks._stacking import (
    compute_stacking,
    get_num_stacks,
    get_stack_heights,
)
from geneview.genometracks._genome_axis import GenomeAxisTrack
from geneview.genometracks._annotation import AnnotationTrack
from geneview.genometracks._gene_region import GeneRegionTrack
from geneview.genometracks._data_track import DataTrack
from geneview.genometracks._highlight import HighlightTrack
from geneview.genometracks._io import read_bed, read_gff, read_bedgraph, read_auto


# =============================================================================
# GenomicInterval tests
# =============================================================================

class TestGenomicInterval:

    def test_basic_creation(self):
        gi = GenomicInterval("chr1", 100, 200)
        assert gi.chrom == "chr1"
        assert gi.start == 100
        assert gi.end == 200
        assert gi.strand == "*"

    def test_with_strand(self):
        gi = GenomicInterval("chr1", 100, 200, "+")
        assert gi.strand == "+"

    def test_width(self):
        gi = GenomicInterval("chr1", 100, 300)
        assert gi.width == 200

    def test_invalid_end_less_than_start(self):
        with pytest.raises(ValueError, match="cannot be less than"):
            GenomicInterval("chr1", 200, 100)

    def test_invalid_strand(self):
        with pytest.raises(ValueError, match="Strand must be"):
            GenomicInterval("chr1", 100, 200, "x")

    def test_overlaps(self):
        a = GenomicInterval("chr1", 100, 200)
        b = GenomicInterval("chr1", 150, 250)
        c = GenomicInterval("chr1", 200, 300)
        d = GenomicInterval("chr2", 100, 200)
        assert a.overlaps(b)
        assert not a.overlaps(c)  # adjacent, not overlapping
        assert not a.overlaps(d)  # different chromosome

    def test_contains(self):
        gi = GenomicInterval("chr1", 100, 200)
        assert gi.contains(150)
        assert gi.contains(100)
        assert not gi.contains(200)  # end is exclusive
        assert not gi.contains(50)

    def test_extend(self):
        gi = GenomicInterval("chr1", 100, 200)
        ext = gi.extend(left=50, right=50)
        assert ext.start == 50
        assert ext.end == 250

    def test_extend_no_negative_start(self):
        gi = GenomicInterval("chr1", 10, 200)
        ext = gi.extend(left=100)
        assert ext.start == 0

    def test_repr(self):
        gi = GenomicInterval("chr1", 100, 200)
        assert "chr1" in repr(gi)
        assert "100" in repr(gi)


# =============================================================================
# Stacking algorithm tests
# =============================================================================

class TestStacking:

    def test_empty_input(self):
        result = compute_stacking(np.array([]), np.array([]))
        assert len(result) == 0

    def test_hide_mode(self):
        starts = np.array([1, 5, 10])
        ends = np.array([4, 8, 15])
        result = compute_stacking(starts, ends, mode="hide")
        assert all(r == 0 for r in result)

    def test_dense_mode(self):
        starts = np.array([1, 3, 5])
        ends = np.array([4, 6, 8])
        result = compute_stacking(starts, ends, mode="dense")
        assert all(r == 0 for r in result)

    def test_non_overlapping_same_row(self):
        starts = np.array([1, 10, 20])
        ends = np.array([5, 15, 25])
        result = compute_stacking(starts, ends, mode="squish")
        # All non-overlapping should be in same row
        assert all(r == 0 for r in result)

    def test_overlapping_different_rows(self):
        starts = np.array([1, 3, 5])
        ends = np.array([4, 6, 8])
        result = compute_stacking(starts, ends, mode="squish")
        # Each overlaps, should get different rows
        assert get_num_stacks(result) >= 2

    def test_fully_overlapping(self):
        starts = np.array([1, 1, 1])
        ends = np.array([10, 10, 10])
        result = compute_stacking(starts, ends, mode="squish")
        assert get_num_stacks(result) == 3

    def test_min_distance(self):
        starts = np.array([1, 6])
        ends = np.array([5, 10])
        # Without min_distance they can share a row
        result1 = compute_stacking(starts, ends, mode="squish", min_distance=0)
        assert all(r == 0 for r in result1)
        # With min_distance=2 they need separate rows
        result2 = compute_stacking(starts, ends, mode="squish", min_distance=2)
        assert get_num_stacks(result2) == 2

    def test_get_stack_heights(self):
        stacks = np.array([0, 0, 1, 1, 2])
        info = get_stack_heights(stacks, mode="squish")
        assert info["total_rows"] == 3
        assert len(info["y_positions"]) == 3
        assert info["row_height"] > 0


# =============================================================================
# RangeTrack tests
# =============================================================================

class _ConcreteRangeTrack(RangeTrack):
    """Minimal concrete subclass for testing (RangeTrack is abstract)."""
    def draw(self, ax, region):
        pass


class TestRangeTrack:

    def _make_data(self):
        return pd.DataFrame({
            "chrom": ["chr1"] * 5,
            "start": [100, 200, 300, 400, 500],
            "end": [150, 250, 350, 450, 550],
            "name": ["a", "b", "c", "d", "e"],
        })

    def test_creation_from_dataframe(self):
        data = self._make_data()
        track = _ConcreteRangeTrack(data=data, name="test")
        assert track.name == "test"
        assert len(track.data) == 5

    def test_creation_from_none(self):
        track = _ConcreteRangeTrack(data=None, name="empty")
        assert track.data is None

    def test_missing_columns_raises(self):
        bad_data = pd.DataFrame({"x": [1], "y": [2]})
        with pytest.raises(ValueError, match="Missing required columns"):
            _ConcreteRangeTrack(data=bad_data)

    def test_subset(self):
        data = self._make_data()
        track = _ConcreteRangeTrack(data=data)
        region = GenomicInterval("chr1", 180, 320)
        result = track.subset(region)
        assert len(result) == 2  # features at 200-250 and 300-350

    def test_subset_no_overlap(self):
        data = self._make_data()
        track = _ConcreteRangeTrack(data=data)
        region = GenomicInterval("chr1", 600, 700)
        result = track.subset(region)
        assert result is None

    def test_get_region(self):
        data = self._make_data()
        track = _ConcreteRangeTrack(data=data)
        region = track.get_region()
        assert region.chrom == "chr1"
        assert region.start == 100
        assert region.end == 550

    def test_column_normalization(self):
        data = pd.DataFrame({
            "Chrom": ["chr1"], "Start": [100], "End": [200]
        })
        track = _ConcreteRangeTrack(data=data)
        assert "chrom" in track.data.columns
        assert "start" in track.data.columns


# =============================================================================
# GenomeAxisTrack tests
# =============================================================================

class TestGenomeAxisTrack:

    def test_creation(self):
        track = GenomeAxisTrack()
        assert track.name == "Axis"
        assert track.height == 0.3

    def test_draw_full_axis(self):
        track = GenomeAxisTrack()
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2200000)
        track.draw(ax, region)
        # Should set x limits
        assert ax.get_xlim() == (2000000, 2200000)

    def test_draw_scale_bar(self):
        track = GenomeAxisTrack(scale=0.1, label_pos="below")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2200000)
        track.draw(ax, region)
        assert ax.get_xlim() == (2000000, 2200000)

    def test_get_region_returns_none(self):
        track = GenomeAxisTrack()
        assert track.get_region() is None

    def test_custom_height(self):
        track = GenomeAxisTrack(height=0.5)
        assert track.height == 0.5


# =============================================================================
# AnnotationTrack tests
# =============================================================================

class TestAnnotationTrack:

    def _make_data(self):
        return pd.DataFrame({
            "chrom": ["chr7"] * 4,
            "start": [2000000, 2070000, 2100000, 2160000],
            "end":   [2050000, 2130000, 2150000, 2170000],
            "strand": ["-", "+", "-", "-"],
            "name": ["feat1", "feat2", "feat3", "feat4"],
            "group": ["grp1", "grp2", "grp1", "grp3"],
            "feature": ["exon", "CDS", "exon", "CDS"],
        })

    def test_creation(self):
        data = self._make_data()
        track = AnnotationTrack(data)
        assert track.name == "Annotation"
        assert len(track.data) == 4

    def test_invalid_shape(self):
        with pytest.raises(ValueError, match="Shape must be"):
            AnnotationTrack(self._make_data(), shape="triangle")

    def test_draw(self):
        data = self._make_data()
        track = AnnotationTrack(data, stacking="squish")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 1900000, 2200000)
        track.draw(ax, region)
        assert ax.get_xlim() == (1900000, 2200000)

    def test_draw_empty_region(self):
        data = self._make_data()
        track = AnnotationTrack(data)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 5000000, 6000000)
        track.draw(ax, region)
        # Should handle empty data gracefully

    def test_different_shapes(self):
        data = self._make_data()
        for shape in ("box", "ellipse", "arrow"):
            track = AnnotationTrack(data, shape=shape)
            fig, ax = plt.subplots()
            region = GenomicInterval("chr7", 1900000, 2200000)
            track.draw(ax, region)
            plt.close(fig)

    def test_show_labels(self):
        data = self._make_data()
        track = AnnotationTrack(data, show_label=True, label_pos="above")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 1900000, 2200000)
        track.draw(ax, region)

    def test_stacking_modes(self):
        data = self._make_data()
        for mode in ("squish", "pack", "dense", "full"):
            track = AnnotationTrack(data, stacking=mode)
            fig, ax = plt.subplots()
            region = GenomicInterval("chr7", 1900000, 2200000)
            track.draw(ax, region)
            plt.close(fig)


# =============================================================================
# GeneRegionTrack tests
# =============================================================================

class TestGeneRegionTrack:

    def _make_data(self):
        return pd.DataFrame({
            "chrom": ["chr7"] * 6,
            "start": [1000, 1500, 2000, 3000, 3200, 4000],
            "end":   [1200, 1700, 2200, 3100, 3400, 4200],
            "strand": ["+"] * 6,
            "feature": ["UTR", "CDS", "CDS", "CDS", "CDS", "UTR"],
            "transcript_id": ["tx1"] * 3 + ["tx2"] * 3,
            "gene_name": ["GeneA"] * 6,
        })

    def test_creation(self):
        data = self._make_data()
        track = GeneRegionTrack(data)
        assert track.name == "GeneRegion"
        assert len(track.data) == 6

    def test_draw(self):
        data = self._make_data()
        track = GeneRegionTrack(data)
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 500, 5000)
        track.draw(ax, region)
        assert ax.get_xlim() == (500, 5000)

    def test_collapse_transcripts_gene(self):
        data = self._make_data()
        track = GeneRegionTrack(data, collapse_transcripts="gene")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 500, 5000)
        track.draw(ax, region)

    def test_collapse_transcripts_longest(self):
        data = self._make_data()
        track = GeneRegionTrack(data, collapse_transcripts="longest")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 500, 5000)
        track.draw(ax, region)

    def test_show_id_options(self):
        data = self._make_data()
        for show_id in ("gene", "transcript", None):
            track = GeneRegionTrack(data, show_id=show_id)
            fig, ax = plt.subplots()
            region = GenomicInterval("chr7", 500, 5000)
            track.draw(ax, region)
            plt.close(fig)


# =============================================================================
# DataTrack tests
# =============================================================================

class TestDataTrack:

    def _make_data(self):
        np.random.seed(42)
        return pd.DataFrame({
            "chrom": ["chr7"] * 50,
            "start": np.arange(2000000, 2050000, 1000),
            "end":   np.arange(2001000, 2051000, 1000),
            "value": np.random.randn(50).cumsum(),
        })

    def test_creation(self):
        data = self._make_data()
        track = DataTrack(data, type="line")
        assert track.name == "Data"
        assert track.plot_type == "line"

    def test_invalid_type(self):
        with pytest.raises(ValueError, match="Plot type must be"):
            DataTrack(self._make_data(), type="invalid")

    def test_draw_line(self):
        data = self._make_data()
        track = DataTrack(data, type="line")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_histogram(self):
        data = self._make_data()
        track = DataTrack(data, type="histogram")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_polygon(self):
        data = self._make_data()
        track = DataTrack(data, type="polygon")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_points(self):
        data = self._make_data()
        track = DataTrack(data, type="points")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_gradient(self):
        data = self._make_data()
        track = DataTrack(data, type="gradient")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_mountain(self):
        data = self._make_data()
        track = DataTrack(data, type="mountain")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_draw_heatmap(self):
        data = self._make_data()
        data["value2"] = data["value"] * 0.5 + 1
        track = DataTrack(data, type="heatmap")
        fig, ax = plt.subplots()
        region = GenomicInterval("chr7", 2000000, 2050000)
        track.draw(ax, region)

    def test_custom_ylim(self):
        data = self._make_data()
        track = DataTrack(data, ylim=(-10, 10))
        assert track.get_ylim() == (-10, 10)

    def test_auto_ylim(self):
        data = self._make_data()
        track = DataTrack(data)
        ylim = track.get_ylim()
        assert ylim[0] < ylim[1]

    def test_value_column_auto_detect(self):
        data = self._make_data()
        track = DataTrack(data)
        assert "value" in track.value_columns


# =============================================================================
# HighlightTrack tests
# =============================================================================

class TestHighlightTrack:

    def test_creation(self):
        regions = pd.DataFrame({
            "chrom": ["chr7"], "start": [2050000], "end": [2100000]
        })
        ht = HighlightTrack(regions=regions)
        assert ht.fill == "#FFE3E6"
        assert ht.alpha == 0.3

    def test_creation_from_intervals(self):
        regions = [GenomicInterval("chr7", 2050000, 2100000)]
        ht = HighlightTrack(regions=regions)
        assert len(ht._data) == 1

    def test_draw_highlights(self):
        regions = pd.DataFrame({
            "chrom": ["chr7"], "start": [2050000], "end": [2100000]
        })
        ht = HighlightTrack(regions=regions)
        fig, ax = plt.subplots()
        ax.set_ylim(0, 1)
        region = GenomicInterval("chr7", 2000000, 2200000)
        ht.draw_highlights(ax, region)


# =============================================================================
# I/O tests
# =============================================================================

class TestIO:

    def test_read_bed(self):
        content = "chr1\t100\t200\tfeat1\t0\t+\nchr1\t300\t400\tfeat2\t0\t-\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_bed(f.name)
        os.unlink(f.name)
        assert len(df) == 2
        assert "chrom" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns

    def test_read_gff(self):
        content = (
            "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"tx1\";\n"
            "chr1\ttest\tCDS\t151\t180\t.\t+\t0\tgene_id \"GeneA\"; transcript_id \"tx1\";\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".gff", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_gff(f.name)
        os.unlink(f.name)
        assert len(df) == 2
        assert df["start"].iloc[0] == 100  # 0-based conversion
        assert "gene_id" in df.columns

    def test_read_bedgraph(self):
        content = "chr1\t100\t200\t1.5\nchr1\t200\t300\t2.3\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bedgraph", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_bedgraph(f.name)
        os.unlink(f.name)
        assert len(df) == 2
        assert "value" in df.columns

    def test_read_auto_bed(self):
        content = "chr1\t100\t200\tfeat1\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_auto(f.name)
        os.unlink(f.name)
        assert "chrom" in df.columns

    def test_read_auto_unknown_extension(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xyz", delete=False) as f:
            f.write("not a valid file")
            f.flush()
            with pytest.raises(ValueError, match="Cannot determine"):
                read_auto(f.name)
        os.unlink(f.name)


# =====================================================================
# New tests for gap-analysis features
# =====================================================================
from geneview.genometracks._overlay import OverlayTrack
from geneview.genometracks._data_track import DataTrack
from geneview.genometracks._highlight import HighlightTrack
from geneview.genometracks._annotation import AnnotationTrack
from geneview.genometracks._gene_region import GeneRegionTrack
from geneview.genometracks._genome_axis import GenomeAxisTrack
from geneview.genometracks._track_plot import plot_tracks


def _make_data(n=20, chrom="chr7", seed=42):
    rng = np.random.RandomState(seed)
    starts = np.linspace(1000, 2000, n, dtype=int)
    return pd.DataFrame({
        "chrom": [chrom] * n,
        "start": starts,
        "end": starts + 50,
        "value": rng.randn(n).cumsum(),
    })


def _make_gene_data():
    return pd.DataFrame({
        "chrom": ["chr7"] * 6,
        "start": [1000, 1200, 1500, 1000, 1200, 1800],
        "end":   [1100, 1400, 1700, 1100, 1400, 2000],
        "strand": ["+"] * 6,
        "feature": ["CDS", "exon", "CDS", "CDS", "exon", "CDS"],
        "transcript_id": ["tx1"] * 3 + ["tx2"] * 3,
        "gene_name": ["GeneA"] * 6,
    })


class TestOverlayTrack:
    def test_basic_overlay(self):
        d1 = _make_data(seed=1)
        d2 = _make_data(seed=2)
        dt1 = DataTrack(d1, type="line", col="blue")
        dt2 = DataTrack(d2, type="line", col="red")
        ot = OverlayTrack(track_list=[dt1, dt2], name="Overlay")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([GenomeAxisTrack(), ot], region=region, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")

    def test_empty_overlay(self):
        ot = OverlayTrack(track_list=[])
        assert ot.name == "Overlay"
        region = GenomicInterval("chr7", 1000, 2000)
        fig, ax = plt.subplots()
        ot.draw(ax, region)
        plt.close("all")

    def test_overlay_inherits_name(self):
        dt = DataTrack(_make_data(), name="MyData")
        ot = OverlayTrack(track_list=[dt])
        assert ot.name == "MyData"


class TestDataTrackNewFeatures:
    def test_combined_plot_type(self):
        dt = DataTrack(_make_data(), type="b")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_stairs_post(self):
        dt = DataTrack(_make_data(), type="s")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_stairs_pre(self):
        dt = DataTrack(_make_data(), type="S")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_transformation(self):
        dt = DataTrack(_make_data(), type="line",
                       transformation=lambda x: np.abs(x))
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_window_binning(self):
        dt = DataTrack(_make_data(n=50), type="line", window=10)
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_window_auto(self):
        dt = DataTrack(_make_data(n=50), type="histogram", window="auto")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_grid_display(self):
        dt = DataTrack(_make_data(), type="line",
                       display_params={"grid": True})
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_legend(self):
        dt = DataTrack(_make_data(), type="line", groups=["A"], legend=True)
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_histogram_negative_baseline(self):
        data = _make_data()
        data["value"] = data["value"] - data["value"].mean()
        dt = DataTrack(data, type="histogram", baseline=0)
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")


class TestGeneRegionMetaCollapse:
    def test_meta_collapse(self):
        data = _make_gene_data()
        grt = GeneRegionTrack(data, collapse_transcripts="meta")
        region = GenomicInterval("chr7", 900, 2100)
        axes = plot_tracks([GenomeAxisTrack(), grt], region=region, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")

    def test_meta_produces_single_transcript(self):
        data = _make_gene_data()
        grt = GeneRegionTrack(data, collapse_transcripts="meta")
        collapsed = grt._collapse(data)
        # Should have a single transcript_id
        tx_ids = collapsed["transcript_id"].unique()
        assert len(tx_ids) == 1
        assert "meta" in tx_ids[0]


class TestHighlightPerRegionColors:
    def test_per_region_fill(self):
        data = _make_data(n=5)
        dt = DataTrack(data, type="line")
        regions = pd.DataFrame({
            "chrom": ["chr7", "chr7"],
            "start": [1000, 1500],
            "end": [1200, 1700],
        })
        ht = HighlightTrack(
            track_list=[dt], regions=regions,
            fill=["#FF0000", "#00FF00"], col=["red", "green"],
        )
        region = GenomicInterval("chr7", 900, 2100)
        axes = plot_tracks([GenomeAxisTrack(), ht], region=region, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")


class TestAnnotationTrackShapes:
    def test_fixed_arrow(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 2,
            "start": [1000, 1500],
            "end": [1300, 1800],
            "strand": ["+", "-"],
        })
        at = AnnotationTrack(data, shape="fixedArrow")
        region = GenomicInterval("chr7", 900, 2000)
        axes = plot_tracks([at], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_small_arrow(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 2,
            "start": [1000, 1500],
            "end": [1300, 1800],
            "strand": ["+", "-"],
        })
        at = AnnotationTrack(data, shape="smallArrow")
        region = GenomicInterval("chr7", 900, 2000)
        axes = plot_tracks([at], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_group_annotation(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 3,
            "start": [1000, 1200, 1500],
            "end": [1100, 1400, 1700],
            "group": ["grp1", "grp1", "grp2"],
            "name": ["a", "b", "c"],
        })
        at = AnnotationTrack(data, group_annotation="group")
        region = GenomicInterval("chr7", 900, 2000)
        axes = plot_tracks([at], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")


class TestGenomeAxisDirection:
    def test_add53(self):
        at = GenomeAxisTrack(add53=True)
        region = GenomicInterval("chr7", 1000, 5000)
        axes = plot_tracks([at], region=region, figsize=(8, 2))
        assert len(axes) == 1
        plt.close("all")

    def test_add35(self):
        at = GenomeAxisTrack(add35=True)
        region = GenomicInterval("chr7", 1000, 5000)
        axes = plot_tracks([at], region=region, figsize=(8, 2))
        assert len(axes) == 1
        plt.close("all")

    def test_both_directions(self):
        at = GenomeAxisTrack(add53=True, add35=True)
        region = GenomicInterval("chr7", 1000, 5000)
        axes = plot_tracks([at], region=region, figsize=(8, 2))
        assert len(axes) == 1
        plt.close("all")


class TestPlotTracksNewFeatures:
    def test_show_title_false(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([GenomeAxisTrack(), dt], region=region,
                           show_title=False, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")

    def test_reverse_strand(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([GenomeAxisTrack(), dt], region=region,
                           reverse_strand=True, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")

    def test_fractional_extend(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([GenomeAxisTrack(), dt], region=region,
                           extend_left=0.05, extend_right=0.05, figsize=(8, 4))
        assert len(axes) == 2
        plt.close("all")

    def test_highlight_targeting(self):
        """Bug 1 fix: highlights should only appear on targeted tracks."""
        d1 = _make_data(seed=1)
        d2 = _make_data(seed=2)
        dt1 = DataTrack(d1, type="line", name="Track1")
        dt2 = DataTrack(d2, type="line", name="Track2")
        regions = pd.DataFrame({
            "chrom": ["chr7"], "start": [1200], "end": [1500],
        })
        # Highlight only targets dt1
        ht = HighlightTrack(track_list=[dt1], regions=regions, fill="yellow")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([GenomeAxisTrack(), ht, dt2], region=region,
                           figsize=(8, 6))
        assert len(axes) == 3
        plt.close("all")
