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
    GeneRegionTrack, AlignmentsTrack, GenomicInterval,
)
from geneview.plotstyle import get_style, list_styles, use_style


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


# ---------------------------------------------------------------------------
# Test: multi-category colour-blind-safe recolouring (journal styles only)
# ---------------------------------------------------------------------------

class TestPlotStyleCategoricalFields:
    """New multi-category colour slots on PlotStyle."""

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_fields_exist(self, name):
        style = get_style(name)
        assert hasattr(style, "tracks_categorical_palette")
        assert hasattr(style, "tracks_strand_fwd_color")
        assert hasattr(style, "tracks_strand_rev_color")

    def test_default_style_off(self):
        """Default geneview leaves the categorical slots disabled."""
        style = get_style("geneview")
        assert style.tracks_categorical_palette == []
        assert style.tracks_strand_fwd_color is None
        assert style.tracks_strand_rev_color is None

    @pytest.mark.parametrize("name", ["nature", "science", "cell"])
    def test_journal_styles_populated(self, name):
        style = get_style(name)
        assert style.tracks_categorical_palette  # non-empty
        assert style.tracks_categorical_palette == style.color_palette
        assert style.tracks_strand_fwd_color == "#0072B2"
        assert style.tracks_strand_rev_color == "#D55E00"


class TestPanelLabelStyleFields:
    """Panel-label styling slots on PlotStyle."""

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_fields_exist(self, name):
        style = get_style(name)
        assert hasattr(style, "panel_label_fontsize")
        assert hasattr(style, "panel_label_fontweight")
        assert hasattr(style, "panel_label_uppercase")

    def test_default_style_lowercase(self):
        style = get_style("geneview")
        assert style.panel_label_uppercase is False
        assert style.panel_label_fontsize == 12.0

    def test_nature_lowercase_8pt(self):
        style = get_style("nature")
        assert style.panel_label_uppercase is False
        assert style.panel_label_fontsize == 8.0

    @pytest.mark.parametrize("name", ["science", "cell"])
    def test_science_cell_uppercase(self, name):
        style = get_style(name)
        assert style.panel_label_uppercase is True
        assert style.panel_label_fontsize == 8.0


class TestAlignmentsStrandColors:
    """AlignmentsTrack forward/reverse strand colours follow the style."""

    def test_journal_strand_pair(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        with use_style("nature"):
            fwd, rev = track._strand_fill_colors()
        assert fwd == "#0072B2"
        assert rev == "#D55E00"

    def test_default_keeps_historical(self):
        track = AlignmentsTrack(filepath="dummy.bam")
        with use_style("geneview"):
            fwd, rev = track._strand_fill_colors()
        assert fwd == "#E89E9D"
        assert rev == "#8C8FCE"

    def test_explicit_colors_win(self):
        track = AlignmentsTrack(
            filepath="dummy.bam",
            fill_reads_fwd="red", fill_reads_rev="lime",
        )
        with use_style("nature"):
            fwd, rev = track._strand_fill_colors()
        assert fwd == "red"
        assert rev == "lime"


class TestDataTrackMultiSeriesColors:
    """DataTrack multi-series colours cycle the journal palette."""

    def _multi_track(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 5,
            "start": np.arange(0, 5000, 1000),
            "end": np.arange(1000, 6000, 1000),
            "a": np.arange(5, dtype=float),
            "b": np.arange(5, dtype=float),
            "c": np.arange(5, dtype=float),
        })
        return DataTrack(data, value_columns=["a", "b", "c"], type="line")

    def test_journal_cycles_palette(self):
        track = self._multi_track()
        pal = get_style("nature").tracks_categorical_palette
        with use_style("nature"):
            colors = track._multi_series_colors("#0080FF")
        assert colors == list(pal[:3])

    def test_default_repeats_single(self):
        track = self._multi_track()
        with use_style("geneview"):
            colors = track._multi_series_colors("#0080FF")
        assert colors == ["#0080FF", "#0080FF", "#0080FF"]

    def test_explicit_list_wins(self):
        track = self._multi_track()
        explicit = ["red", "green", "blue"]
        with use_style("nature"):
            colors = track._multi_series_colors(explicit)
        assert colors == explicit

    def test_single_series_not_recoloured(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 5,
            "start": np.arange(0, 5000, 1000),
            "end": np.arange(1000, 6000, 1000),
            "a": np.arange(5, dtype=float),
        })
        track = DataTrack(data, value_columns=["a"], type="line")
        with use_style("nature"):
            colors = track._multi_series_colors("#0080FF")
        assert colors == ["#0080FF"]


class TestAnnotationCategoricalColors:
    """Unknown annotation feature types use the journal palette."""

    def _unknown_feature_track(self):
        data = pd.DataFrame({
            "chrom": ["chr7"] * 3,
            "start": [2_000_000, 2_010_000, 2_020_000],
            "end": [2_005_000, 2_015_000, 2_025_000],
            "strand": ["+", "-", "+"],
            "feature": ["typeA", "typeB", "typeC"],
        })
        return AnnotationTrack(data)

    def test_nature_uses_palette(self, region):
        import matplotlib.colors as mcolors
        track = self._unknown_feature_track()
        axes = plot_tracks([GenomeAxisTrack(), track], region=region, style="nature")
        pal = get_style("nature").tracks_categorical_palette
        pal_rgba = {tuple(round(v, 3) for v in mcolors.to_rgba(c)) for c in pal}
        seen = set()
        for ax in axes:
            for patch in ax.patches:
                seen.add(tuple(round(v, 3) for v in patch.get_facecolor()))
        assert seen & pal_rgba, "expected at least one feature drawn with the journal palette"
        plt.close("all")

