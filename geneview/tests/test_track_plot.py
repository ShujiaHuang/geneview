"""Tests for geneview.genometracks.plot_tracks integration.

Tests cover the main plot_tracks orchestrator function with various
combinations of track types, ensuring correct rendering and layout.
"""
import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._base import GenomicInterval
from geneview.genometracks._genome_axis import GenomeAxisTrack
from geneview.genometracks._annotation import AnnotationTrack
from geneview.genometracks._gene_region import GeneRegionTrack
from geneview.genometracks._data_track import DataTrack
from geneview.genometracks._highlight import HighlightTrack
from geneview.genometracks._track_plot import plot_tracks


def _make_annotation_data():
    return pd.DataFrame({
        "chrom": ["chr7"] * 4,
        "start": [2000000, 2070000, 2100000, 2160000],
        "end":   [2050000, 2130000, 2150000, 2170000],
        "strand": ["-", "+", "-", "-"],
        "name": ["feat1", "feat2", "feat3", "feat4"],
    })


def _make_gene_data():
    return pd.DataFrame({
        "chrom": ["chr7"] * 6,
        "start": [2001000, 2001500, 2002000, 2070000, 2070500, 2071000],
        "end":   [2001200, 2001700, 2002200, 2070200, 2070700, 2071200],
        "strand": ["+"] * 6,
        "feature": ["UTR", "CDS", "UTR", "CDS", "CDS", "UTR"],
        "transcript_id": ["tx1"] * 3 + ["tx2"] * 3,
        "gene_name": ["GeneA"] * 3 + ["GeneB"] * 3,
    })


def _make_data_track_data():
    np.random.seed(42)
    return pd.DataFrame({
        "chrom": ["chr7"] * 50,
        "start": np.arange(2000000, 2050000, 1000),
        "end":   np.arange(2001000, 2051000, 1000),
        "value": np.random.randn(50).cumsum(),
    })


class TestPlotTracks:

    def test_single_axis_track(self):
        """Should plot a single GenomeAxisTrack."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks(track, region=region)
        assert len(axes) == 1

    def test_multiple_tracks(self):
        """Should plot multiple tracks stacked."""
        ax_track = GenomeAxisTrack()
        ann_track = AnnotationTrack(_make_annotation_data())
        region = GenomicInterval("chr7", 1900000, 2200000)
        axes = plot_tracks([ax_track, ann_track], region=region)
        assert len(axes) == 2

    def test_auto_region_detection(self):
        """Should auto-detect region from track data."""
        ann_track = AnnotationTrack(_make_annotation_data())
        axes = plot_tracks([ann_track])
        assert len(axes) == 1

    def test_no_region_no_data_raises(self):
        """Should raise when no region can be determined."""
        ax_track = GenomeAxisTrack()
        with pytest.raises(ValueError, match="Cannot determine"):
            plot_tracks([ax_track])

    def test_with_title(self):
        """Should add a main title."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks(track, region=region, title="My Plot")
        assert len(axes) >= 1

    def test_custom_sizes(self):
        """Should respect custom track sizes."""
        ax_track = GenomeAxisTrack()
        ann_track = AnnotationTrack(_make_annotation_data())
        region = GenomicInterval("chr7", 1900000, 2200000)
        axes = plot_tracks([ax_track, ann_track], region=region, sizes=[1, 3])
        assert len(axes) == 2

    def test_custom_figsize(self):
        """Should use provided figure size."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks(track, region=region, figsize=(8, 4))
        fig = plt.gcf()
        assert abs(fig.get_size_inches()[0] - 8) < 0.1

    def test_extend_left_right(self):
        """Should extend the plotting region."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks(track, region=region,
                          extend_left=100000, extend_right=100000)
        # Axes should reflect extended region

    def test_plot_into_existing_ax(self):
        """Should plot into a provided axes."""
        fig, ax = plt.subplots()
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        result = plot_tracks(track, region=region, ax=ax)
        assert result[0] is ax

    def test_empty_track_list_raises(self):
        """Should raise ValueError for empty track list."""
        with pytest.raises(ValueError, match="No tracks"):
            plot_tracks([])


class TestPlotTracksIntegration:
    """Integration tests with various track combinations."""

    def test_axis_plus_annotation(self):
        """GenomeAxisTrack + AnnotationTrack."""
        ax_track = GenomeAxisTrack()
        ann_track = AnnotationTrack(_make_annotation_data())
        region = GenomicInterval("chr7", 1900000, 2200000)
        axes = plot_tracks([ax_track, ann_track], region=region)
        assert len(axes) == 2

    def test_axis_plus_gene_region(self):
        """GenomeAxisTrack + GeneRegionTrack."""
        ax_track = GenomeAxisTrack()
        gene_track = GeneRegionTrack(_make_gene_data())
        region = GenomicInterval("chr7", 1999000, 2073000)
        axes = plot_tracks([ax_track, gene_track], region=region)
        assert len(axes) == 2

    def test_axis_plus_data_track(self):
        """GenomeAxisTrack + DataTrack."""
        ax_track = GenomeAxisTrack()
        data_track = DataTrack(_make_data_track_data(), type="line")
        region = GenomicInterval("chr7", 2000000, 2050000)
        axes = plot_tracks([ax_track, data_track], region=region)
        assert len(axes) == 2

    def test_all_track_types(self):
        """All track types together."""
        ax_track = GenomeAxisTrack()
        ann_track = AnnotationTrack(_make_annotation_data())
        gene_track = GeneRegionTrack(_make_gene_data())
        data_track = DataTrack(_make_data_track_data(), type="histogram")
        region = GenomicInterval("chr7", 2000000, 2080000)
        axes = plot_tracks(
            [ax_track, ann_track, gene_track, data_track],
            region=region,
        )
        assert len(axes) == 4

    def test_with_highlight(self):
        """HighlightTrack wrapping other tracks."""
        ax_track = GenomeAxisTrack()
        ann_track = AnnotationTrack(_make_annotation_data())
        hl_regions = pd.DataFrame({
            "chrom": ["chr7"], "start": [2050000], "end": [2100000]
        })
        ht = HighlightTrack(
            track_list=[ax_track, ann_track],
            regions=hl_regions,
        )
        region = GenomicInterval("chr7", 1900000, 2200000)
        axes = plot_tracks([ht], region=region)
        # Expanded: 2 tracks from HighlightTrack
        assert len(axes) == 2

    def test_data_track_types(self):
        """Various DataTrack plot types."""
        data = _make_data_track_data()
        region = GenomicInterval("chr7", 2000000, 2050000)

        for plot_type in ("line", "histogram", "polygon", "points", "gradient", "mountain"):
            track = DataTrack(data, type=plot_type)
            axes = plot_tracks([track], region=region)
            plt.close("all")
            assert len(axes) == 1

    def test_annotation_stacking_modes(self):
        """AnnotationTrack with different stacking modes."""
        data = _make_annotation_data()
        region = GenomicInterval("chr7", 1900000, 2200000)

        for mode in ("squish", "pack", "dense", "full"):
            track = AnnotationTrack(data, stacking=mode)
            axes = plot_tracks([track], region=region)
            plt.close("all")
            assert len(axes) == 1

    def test_scale_bar_axis(self):
        """GenomeAxisTrack with scale bar."""
        track = GenomeAxisTrack(scale=0.1, label_pos="below")
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks([track], region=region)
        assert len(axes) == 1

    def test_main_title_parameter(self):
        """main parameter as alias for title."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks([track], region=region, main="Test Title")
        assert len(axes) >= 1


class TestPlotTracksEdgeCases:

    def test_single_feature_annotation(self):
        """AnnotationTrack with just one feature."""
        data = pd.DataFrame({
            "chrom": ["chr7"], "start": [2000000], "end": [2100000]
        })
        track = AnnotationTrack(data)
        region = GenomicInterval("chr7", 1900000, 2200000)
        axes = plot_tracks([track], region=region)
        assert len(axes) == 1

    def test_many_overlapping_features(self):
        """AnnotationTrack with many overlapping features."""
        np.random.seed(42)
        n = 100
        data = pd.DataFrame({
            "chrom": ["chr7"] * n,
            "start": np.random.randint(2000000, 2100000, n),
            "end": np.random.randint(2100000, 2200000, n),
        })
        data["end"] = data["start"] + np.random.randint(1000, 50000, n)
        track = AnnotationTrack(data, stacking="squish")
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks([track], region=region)
        assert len(axes) == 1

    def test_empty_data_track(self):
        """DataTrack with no data in region."""
        data = pd.DataFrame({
            "chrom": ["chr7"] * 5,
            "start": [100, 200, 300, 400, 500],
            "end": [150, 250, 350, 450, 550],
            "value": [1.0, 2.0, 3.0, 4.0, 5.0],
        })
        track = DataTrack(data)
        region = GenomicInterval("chr7", 5000000, 6000000)  # No data here
        axes = plot_tracks([track], region=region)
        assert len(axes) == 1

    def test_kwargs_passed_to_tracks(self):
        """Extra kwargs should be applied to all tracks."""
        track = GenomeAxisTrack()
        region = GenomicInterval("chr7", 2000000, 2200000)
        axes = plot_tracks([track], region=region, alpha=0.5)
        assert len(axes) >= 1
