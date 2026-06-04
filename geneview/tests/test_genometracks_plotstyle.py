"""Tests for genome tracks + plotstyle integration.

Verifies that ``plot_tracks()`` correctly accepts and applies journal
styles (nature, science, cell, geneview) to genome track figures.
"""
import numpy as np
import pandas as pd
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    plot_tracks, GenomeAxisTrack, AnnotationTrack, DataTrack,
    GeneRegionTrack, GenomicInterval,
)
from geneview.plotstyle import get_style, list_styles


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_data():
    """Numeric data for a DataTrack."""
    np.random.seed(42)
    return pd.DataFrame({
        "chrom": ["chr7"] * 50,
        "start": np.arange(2_000_000, 2_050_000, 1000),
        "end": np.arange(2_001_000, 2_051_000, 1000),
        "value": np.random.randn(50).cumsum(),
    })


@pytest.fixture
def sample_annotations():
    """Annotation data for an AnnotationTrack."""
    return pd.DataFrame({
        "chrom": ["chr7"] * 4,
        "start": [2_000_000, 2_010_000, 2_020_000, 2_035_000],
        "end": [2_005_000, 2_015_000, 2_030_000, 2_040_000],
        "strand": ["+", "-", "+", "-"],
        "name": ["GeneA", "GeneB", "GeneC", "GeneD"],
    })


@pytest.fixture
def region():
    return GenomicInterval("chr7", 2_000_000, 2_050_000)


def _make_tracks(sample_data, sample_annotations):
    """Create a standard 3-track set."""
    return [
        GenomeAxisTrack(),
        AnnotationTrack(sample_annotations, stacking="squish"),
        DataTrack(sample_data, type="histogram"),
    ]


# ---------------------------------------------------------------------------
# Test: style parameter accepted
# ---------------------------------------------------------------------------

class TestStyleParamAccepted:
    """plot_tracks() accepts style= without errors."""

    @pytest.mark.parametrize("style_name", ["geneview", "nature", "science", "cell"])
    def test_style_accepted(self, sample_data, sample_annotations, region, style_name):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style=style_name)
        assert len(axes) == 3
        plt.close("all")

    def test_style_none(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style=None)
        assert len(axes) == 3
        plt.close("all")

    def test_style_invalid_raises(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        with pytest.raises(ValueError, match="Unknown plot style"):
            plot_tracks(tracks, region=region, style="nonexistent")


# ---------------------------------------------------------------------------
# Test: figure size matches style definition
# ---------------------------------------------------------------------------

class TestFigureSize:
    """Figure dimensions should reflect the style's tracks_figsize_width."""

    def test_default_width(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style=None)
        fig = plt.gcf()
        assert fig.get_size_inches()[0] == pytest.approx(12.0, abs=0.5)
        plt.close("all")

    def test_nature_width(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style="nature")
        fig = plt.gcf()
        style = get_style("nature")
        assert fig.get_size_inches()[0] == pytest.approx(
            style.tracks_figsize_width, abs=0.5
        )
        plt.close("all")

    def test_science_width(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style="science")
        fig = plt.gcf()
        style = get_style("science")
        assert fig.get_size_inches()[0] == pytest.approx(
            style.tracks_figsize_width, abs=0.5
        )
        plt.close("all")

    def test_cell_width(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style="cell")
        fig = plt.gcf()
        style = get_style("cell")
        assert fig.get_size_inches()[0] == pytest.approx(
            style.tracks_figsize_width, abs=0.5
        )
        plt.close("all")

    def test_explicit_figsize_overrides_style(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        axes = plot_tracks(tracks, region=region, style="nature", figsize=(10, 5))
        fig = plt.gcf()
        assert fig.get_size_inches()[0] == pytest.approx(10.0, abs=0.1)
        assert fig.get_size_inches()[1] == pytest.approx(5.0, abs=0.1)
        plt.close("all")


# ---------------------------------------------------------------------------
# Test: track display params overridden by style
# ---------------------------------------------------------------------------

class TestTrackParamOverrides:
    """Style should cascade into track display parameters."""

    def test_title_panel_bg_nature(self, sample_data, sample_annotations, region):
        style = get_style("nature")
        tracks = _make_tracks(sample_data, sample_annotations)
        # Trigger the style override
        plot_tracks(tracks, region=region, style="nature")
        # After plot_tracks, tracks should have the style's title bg
        for track in tracks:
            assert track.get_param("background_title") == style.tracks_title_bg
        plt.close("all")

    def test_title_panel_bg_default(self, sample_data, sample_annotations, region):
        style = get_style("geneview")
        tracks = _make_tracks(sample_data, sample_annotations)
        plot_tracks(tracks, region=region, style="geneview")
        for track in tracks:
            assert track.get_param("background_title") == style.tracks_title_bg
        plt.close("all")

    def test_fontsize_overridden(self, sample_data, sample_annotations, region):
        style = get_style("nature")
        tracks = _make_tracks(sample_data, sample_annotations)
        plot_tracks(tracks, region=region, style="nature")
        for track in tracks:
            assert track.get_param("fontsize") == style.tracks_feature_fontsize
        plt.close("all")

    def test_axis_color_overridden(self, sample_data, sample_annotations, region):
        style = get_style("science")
        tracks = _make_tracks(sample_data, sample_annotations)
        plot_tracks(tracks, region=region, style="science")
        for track in tracks:
            assert track.get_param("col_axis") == style.tracks_axis_color
        plt.close("all")

    def test_lwd_overridden(self, sample_data, sample_annotations, region):
        style = get_style("cell")
        tracks = _make_tracks(sample_data, sample_annotations)
        plot_tracks(tracks, region=region, style="cell")
        for track in tracks:
            assert track.get_param("lwd") == style.tracks_linewidth
        plt.close("all")


# ---------------------------------------------------------------------------
# Test: user kwargs override style
# ---------------------------------------------------------------------------

class TestUserKwargsOverrideStyle:
    """User-provided kwargs in plot_tracks should take priority over style."""

    def test_kwarg_overrides_style(self, sample_data, sample_annotations, region):
        tracks = _make_tracks(sample_data, sample_annotations)
        plot_tracks(tracks, region=region, style="nature", fontsize=20)
        for track in tracks:
            assert track.get_param("fontsize") == 20
        plt.close("all")


# ---------------------------------------------------------------------------
# Test: to_track_params() method
# ---------------------------------------------------------------------------

class TestToTrackParams:

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_returns_dict(self, name):
        style = get_style(name)
        params = style.to_track_params()
        assert isinstance(params, dict)
        assert len(params) > 0

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_required_keys_present(self, name):
        style = get_style(name)
        params = style.to_track_params()
        required = {
            "background_title", "col_title", "fontsize_title",
            "col_border_title", "col_axis", "fontsize", "fontcolor", "lwd",
        }
        assert required.issubset(params.keys())


# ---------------------------------------------------------------------------
# Test: PlotStyle tracks_* fields exist
# ---------------------------------------------------------------------------

class TestPlotStyleTracksFields:

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_tracks_fields_exist(self, name):
        style = get_style(name)
        assert hasattr(style, "tracks_title_bg")
        assert hasattr(style, "tracks_title_color")
        assert hasattr(style, "tracks_title_fontsize")
        assert hasattr(style, "tracks_title_border")
        assert hasattr(style, "tracks_axis_color")
        assert hasattr(style, "tracks_axis_linewidth")
        assert hasattr(style, "tracks_tick_fontsize")
        assert hasattr(style, "tracks_feature_fontsize")
        assert hasattr(style, "tracks_linewidth")
        assert hasattr(style, "tracks_figsize_width")
        assert hasattr(style, "tracks_height_per_track")

    def test_nature_compact(self):
        nature = get_style("nature")
        default = get_style("geneview")
        # Nature should be more compact
        assert nature.tracks_figsize_width < default.tracks_figsize_width
        assert nature.tracks_feature_fontsize < default.tracks_feature_fontsize
        assert nature.tracks_height_per_track < default.tracks_height_per_track

    def test_journal_styles_white_title_bg(self):
        for name in ["nature", "science", "cell"]:
            style = get_style(name)
            assert style.tracks_title_bg == "white"

    def test_default_style_gray_title_bg(self):
        style = get_style("geneview")
        assert style.tracks_title_bg == "#D3D3D3"


# ---------------------------------------------------------------------------
# Test: GeneRegionTrack with style
# ---------------------------------------------------------------------------

class TestGeneRegionTrackStyle:

    def test_gene_region_with_nature(self, region):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 6,
            "start": [2_000_000, 2_010_000, 2_020_000,
                      2_030_000, 2_035_000, 2_040_000],
            "end": [2_005_000, 2_015_000, 2_025_000,
                    2_032_000, 2_038_000, 2_045_000],
            "strand": ["+"] * 6,
            "feature": ["UTR", "CDS", "CDS", "CDS", "CDS", "UTR"],
            "transcript_id": ["tx1"] * 6,
            "gene_name": ["TestGene"] * 6,
        })
        track = GeneRegionTrack(data)
        ax_track = GenomeAxisTrack()
        axes = plot_tracks([ax_track, track], region=region, style="nature")
        assert len(axes) == 2
        plt.close("all")


# ---------------------------------------------------------------------------
# Test: single-axes mode with style
# ---------------------------------------------------------------------------

class TestSingleAxesMode:

    def test_single_ax_with_style(self, sample_data, region):
        fig, ax = plt.subplots()
        track = DataTrack(sample_data, type="line")
        axes = plot_tracks([track], region=region, ax=ax, style="nature")
        assert len(axes) == 1
        plt.close("all")
