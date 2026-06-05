"""Tests for LolliplotTrack & DandelionTrack (in geneview.genometracks)."""

import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks import (
    LolliplotTrack,
    lolliplot,
    DandelionTrack,
    dandelion_plot,
    GenomicInterval,
    GenomeAxisTrack,
    AnnotationTrack,
    plot_tracks,
)


# ---------------------------------------------------------------------------
# Test data helpers
# ---------------------------------------------------------------------------

def _make_snp_data():
    return pd.DataFrame({
        "chrom": ["chr1"] * 5,
        "start": [10, 100, 400, 600, 800],
        "end": [11, 101, 401, 601, 801],
        "score": [3, 7, 2, 5, 10],
        "label": ["rs10", "rs100", "rs400", "rs600", "rs800"],
        "fill": ["#D55E00", "#0080FF", "#009E73", "#E69F00", "#CC79A7"],
    })


def _make_feature_data():
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1"],
        "start": [1, 501, 1001],
        "end": [120, 1000, 1405],
        "name": ["Domain A", "Kinase", "DNA Binding"],
        "fill": ["#FF8833", "#51C6E6", "#DFA32D"],
        "height": [0.06, 0.05, 0.07],
    })


def _make_dense_snp_data():
    """Dense SNP data for dandelion clustering tests."""
    rng = np.random.RandomState(42)
    positions = np.sort(rng.randint(10, 1400, 30))
    return pd.DataFrame({
        "chrom": ["chr1"] * 30,
        "start": positions,
        "score": rng.randint(1, 8, 30),
        "label": [f"var{i}" for i in range(30)],
    })


def _make_region():
    return GenomicInterval("chr1", 0, 1500)


# ---------------------------------------------------------------------------
# LolliplotTrack tests
# ---------------------------------------------------------------------------

class TestLolliplotTrack:

    def test_basic_circle(self):
        """Basic circle lolliplot track should render without errors."""
        fig, ax = plt.subplots()
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data())
        ax.set_xlim(0, 1500)
        track.draw(ax, _make_region())
        assert ax.get_xlim()[0] < ax.get_xlim()[1]

    def test_no_features(self):
        """Track without features should still render."""
        fig, ax = plt.subplots()
        track = LolliplotTrack(_make_snp_data())
        ax.set_xlim(0, 1500)
        track.draw(ax, _make_region())

    def test_pie_type(self):
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data(), type="pie")
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        track.draw(ax, _make_region())

    def test_pin_type(self):
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data(), type="pin")
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        track.draw(ax, _make_region())

    def test_flag_type(self):
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data(), type="flag")
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        track.draw(ax, _make_region())

    def test_invalid_type_raises(self):
        with pytest.raises(ValueError, match="type must be one of"):
            LolliplotTrack(_make_snp_data(), type="invalid")

    def test_custom_region(self):
        """Custom region should override auto-detection."""
        fig, ax = plt.subplots()
        track = LolliplotTrack(_make_snp_data())
        region = GenomicInterval("chr1", 0, 500)
        ax.set_xlim(region.start, region.end)
        track.draw(ax, region)
        xlim = ax.get_xlim()
        assert xlim[0] == 0
        assert xlim[1] == 500

    def test_get_region(self):
        """get_region() should return genomic extent."""
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data())
        region = track.get_region()
        assert region is not None
        assert region.chrom == "chr1"
        assert region.start <= 10
        assert region.end >= 1405

    def test_into_existing_ax(self):
        """Should draw into a provided axes."""
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1500)
        track = LolliplotTrack(_make_snp_data())
        track.draw(ax, _make_region())

    def test_convenience_function(self):
        """lolliplot() convenience function should work."""
        ax = lolliplot(_make_snp_data(), features=_make_feature_data())
        assert ax is not None

    def test_convenience_with_ax(self):
        """lolliplot() with ax= should draw into provided axes."""
        fig, ax = plt.subplots()
        result_ax = lolliplot(_make_snp_data(), ax=ax)
        assert result_ax is ax

    def test_convenience_with_type(self):
        ax = lolliplot(_make_snp_data(), type="pie")
        assert ax is not None

    def test_missing_start_raises(self):
        bad_df = pd.DataFrame({"chrom": ["chr1"], "pos": [100]})
        with pytest.raises(ValueError, match="missing required columns"):
            LolliplotTrack(bad_df)

    def test_empty_data_in_region(self):
        """Should handle region with no SNPs gracefully."""
        fig, ax = plt.subplots()
        track = LolliplotTrack(_make_snp_data())
        ax.set_xlim(5000, 6000)
        track.draw(ax, GenomicInterval("chr1", 5000, 6000))

    def test_single_snp(self):
        df = pd.DataFrame({"chrom": ["chr1"], "start": [100], "score": [5]})
        fig, ax = plt.subplots()
        track = LolliplotTrack(df)
        ax.set_xlim(0, 200)
        track.draw(ax, GenomicInterval("chr1", 0, 200))

    def test_score_all_ones(self):
        df = _make_snp_data().copy()
        df["score"] = 1
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df).draw(ax, _make_region())

    def test_large_scores_switch_style(self):
        df = pd.DataFrame({
            "chrom": ["chr1"] * 3,
            "start": [100, 200, 300],
            "score": [50, 100, 200],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 500)
        LolliplotTrack(df, lollipop_style_switch_limit=10).draw(
            ax, GenomicInterval("chr1", 0, 500)
        )

    def test_cex_scaling(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(_make_snp_data(), cex=2.0).draw(ax, _make_region())

    def test_no_yaxis(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(_make_snp_data(), show_yaxis=False).draw(ax, _make_region())

    def test_label_on_feature(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(
            _make_snp_data(),
            features=_make_feature_data(),
            label_on_feature=True,
        ).draw(ax, _make_region())

    def test_repr(self):
        track = LolliplotTrack(_make_snp_data(), type="pin")
        r = repr(track)
        assert "LolliplotTrack" in r
        assert "pin" in r


# ---------------------------------------------------------------------------
# DandelionTrack tests
# ---------------------------------------------------------------------------

class TestDandelionTrack:

    def test_basic_fan(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(
            _make_dense_snp_data(), features=_make_feature_data(), type="fan"
        ).draw(ax, _make_region())

    def test_circle_type(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(
            _make_dense_snp_data(), features=_make_feature_data(), type="circle"
        ).draw(ax, _make_region())

    def test_pie_type(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(
            _make_dense_snp_data(), features=_make_feature_data(), type="pie"
        ).draw(ax, _make_region())

    def test_pin_type(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(
            _make_dense_snp_data(), features=_make_feature_data(), type="pin"
        ).draw(ax, _make_region())

    def test_invalid_type_raises(self):
        with pytest.raises(ValueError, match="type must be one of"):
            DandelionTrack(_make_dense_snp_data(), type="invalid")

    def test_no_features(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(_make_dense_snp_data()).draw(ax, _make_region())

    def test_get_region(self):
        track = DandelionTrack(_make_dense_snp_data(), features=_make_feature_data())
        region = track.get_region()
        assert region is not None
        assert region.chrom == "chr1"

    def test_custom_maxgaps(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(_make_dense_snp_data(), maxgaps=0.01).draw(ax, _make_region())
        plt.close("all")
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(_make_dense_snp_data(), maxgaps=0.5).draw(ax, _make_region())

    def test_custom_height_method(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(
            _make_dense_snp_data(),
            height_method=lambda scores: sum(scores),
        ).draw(ax, _make_region())

    def test_custom_region(self):
        fig, ax = plt.subplots()
        region = GenomicInterval("chr1", 0, 500)
        ax.set_xlim(region.start, region.end)
        DandelionTrack(_make_dense_snp_data()).draw(ax, region)
        xlim = ax.get_xlim()
        assert xlim[0] == 0
        assert xlim[1] == 500

    def test_convenience_function(self):
        ax = dandelion_plot(
            _make_dense_snp_data(), features=_make_feature_data()
        )
        assert ax is not None

    def test_convenience_with_ax(self):
        fig, ax = plt.subplots()
        result_ax = dandelion_plot(_make_dense_snp_data(), ax=ax)
        assert result_ax is ax

    def test_convenience_with_type(self):
        ax = dandelion_plot(_make_dense_snp_data(), type="circle")
        assert ax is not None

    def test_empty_data_in_region(self):
        fig, ax = plt.subplots(); ax.set_xlim(50000, 60000)
        DandelionTrack(_make_dense_snp_data()).draw(
            ax, GenomicInterval("chr1", 50000, 60000)
        )

    def test_single_snp(self):
        df = pd.DataFrame({"chrom": ["chr1"], "start": [100], "score": [5]})
        fig, ax = plt.subplots(); ax.set_xlim(0, 200)
        DandelionTrack(df).draw(ax, GenomicInterval("chr1", 0, 200))

    def test_two_snps_one_cluster(self):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 105],
            "score": [3, 7],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 300)
        DandelionTrack(df, maxgaps=0.1).draw(ax, GenomicInterval("chr1", 0, 300))

    def test_show_yaxis(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(_make_dense_snp_data(), show_yaxis=True).draw(ax, _make_region())

    def test_cex_scaling(self):
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(_make_dense_snp_data(), cex=1.5).draw(ax, _make_region())

    def test_repr(self):
        track = DandelionTrack(_make_dense_snp_data(), type="pin")
        r = repr(track)
        assert "DandelionTrack" in r
        assert "pin" in r


# ---------------------------------------------------------------------------
# Feature drawing tests (shared utility)
# ---------------------------------------------------------------------------

class TestFeatures:

    def test_draw_features_basic(self):
        from geneview.genometracks._mutation_features import draw_features
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1500)
        ax.set_ylim(0, 1)
        h = draw_features(ax, _make_feature_data(), 0, 1500)
        assert h > 0

    def test_draw_features_none(self):
        from geneview.genometracks._mutation_features import draw_features
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1500)
        ax.set_ylim(0, 1)
        h = draw_features(ax, None, 0, 1500)
        assert h > 0

    def test_draw_features_empty(self):
        from geneview.genometracks._mutation_features import draw_features
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1500)
        ax.set_ylim(0, 1)
        empty = pd.DataFrame({"start": [], "end": []})
        h = draw_features(ax, empty, 0, 1500)
        assert h > 0


# ---------------------------------------------------------------------------
# Shape drawing tests
# ---------------------------------------------------------------------------

class TestShapes:

    def test_draw_circle(self):
        from geneview.genometracks._mutation_shapes import draw_circle
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        draw_circle(ax, 500, 0.5, 0.02, facecolor="red")

    def test_draw_pie(self):
        from geneview.genometracks._mutation_shapes import draw_pie
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        draw_pie(ax, 500, 0.5, 0.02, [3, 7], ["red", "blue"])

    def test_draw_fan(self):
        from geneview.genometracks._mutation_shapes import draw_fan
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        draw_fan(ax, 500, 0.5, 0.02, 0.7, "green")

    def test_draw_flag(self):
        from geneview.genometracks._mutation_shapes import draw_flag
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1)
        fig.canvas.draw()
        draw_flag(ax, 500, 0.5, 0.02, "orange")

    def test_draw_stem(self):
        from geneview.genometracks._mutation_shapes import draw_stem
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000)
        ax.set_ylim(0, 1)
        draw_stem(ax, 500, 0.1, 0.8)


# ---------------------------------------------------------------------------
# Track integration tests (composed with other tracks via plot_tracks)
# ---------------------------------------------------------------------------

class TestTrackIntegration:

    def test_lollipop_with_genome_axis(self):
        """LolliplotTrack composed with GenomeAxisTrack in plot_tracks."""
        region = GenomicInterval("chr1", 0, 1500)
        tracks = [
            GenomeAxisTrack(),
            LolliplotTrack(_make_snp_data(), features=_make_feature_data()),
        ]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 2

    def test_lollipop_with_annotation_track(self):
        """LolliplotTrack composed with AnnotationTrack."""
        region = GenomicInterval("chr1", 0, 1500)
        ann_data = _make_feature_data().copy()
        tracks = [
            AnnotationTrack(ann_data, name="Domains"),
            LolliplotTrack(_make_snp_data(), features=_make_feature_data()),
        ]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 2

    def test_dandelion_with_genome_axis(self):
        """DandelionTrack composed with GenomeAxisTrack."""
        region = GenomicInterval("chr1", 0, 1500)
        tracks = [
            GenomeAxisTrack(),
            DandelionTrack(_make_dense_snp_data(), features=_make_feature_data()),
        ]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 2

    def test_lollipop_and_dandelion_together(self):
        """Both mutation tracks in the same plot_tracks call."""
        region = GenomicInterval("chr1", 0, 1500)
        tracks = [
            LolliplotTrack(_make_snp_data(), features=_make_feature_data()),
            DandelionTrack(_make_dense_snp_data(), features=_make_feature_data()),
        ]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 2

    def test_three_tracks_together(self):
        """GenomeAxisTrack + AnnotationTrack + LolliplotTrack."""
        region = GenomicInterval("chr1", 0, 1500)
        ann_data = _make_feature_data().copy()
        tracks = [
            GenomeAxisTrack(),
            AnnotationTrack(ann_data, name="Domains"),
            LolliplotTrack(_make_snp_data(), features=_make_feature_data()),
        ]
        axes = plot_tracks(tracks, region=region)
        assert len(axes) == 3

    def test_auto_region_from_track(self):
        """plot_tracks should auto-derive region from LolliplotTrack."""
        track = LolliplotTrack(_make_snp_data(), features=_make_feature_data())
        axes = plot_tracks([track])
        assert len(axes) == 1


# ---------------------------------------------------------------------------
# Per-SNP property tests
# ---------------------------------------------------------------------------

class TestPerSNPProperties:

    def test_per_snp_cex(self):
        """Per-SNP cex column should override global cex."""
        df = _make_snp_data().copy()
        df["cex"] = [0.5, 1.0, 1.5, 2.0, 0.8]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_per_snp_node_label(self):
        """node_label column should render text inside shapes."""
        df = _make_snp_data().copy()
        df["node_label"] = ["A", "B", "C", "D", "E"]
        df["node_label_color"] = "white"
        df["node_label_size"] = 5
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_per_snp_label_rotation(self):
        """label_rotation column should rotate external labels."""
        df = _make_snp_data().copy()
        df["label_rotation"] = [0, 45, 90, 135, 180]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_per_snp_label_color(self):
        """label_color column should set external label colour."""
        df = _make_snp_data().copy()
        df["label_color"] = ["red", "blue", "green", "orange", "purple"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_per_snp_dashline_col(self):
        """dashline_col column should set guide line colour per SNP."""
        df = _make_snp_data().copy()
        df["score"] = [50, 100, 30, 80, 120]  # large to force non-tanghulu
        df["dashline_col"] = ["red", "blue", "green", "orange", "purple"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_nan_cex_falls_back(self):
        """NaN cex should fall back to global default."""
        df = _make_snp_data().copy()
        df["cex"] = [1.0, np.nan, 1.5, np.nan, 0.8]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(), cex=1.2).draw(
            ax, _make_region())


# ---------------------------------------------------------------------------
# Caterpillar layout tests
# ---------------------------------------------------------------------------

class TestCaterpillarLayout:

    def test_caterpillar_basic(self):
        """side column should place SNPs above and below."""
        df = _make_snp_data().copy()
        df["side"] = ["top", "bottom", "top", "bottom", "top"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_caterpillar_all_top(self):
        """All side=top should behave normally."""
        df = _make_snp_data().copy()
        df["side"] = "top"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_caterpillar_all_bottom(self):
        """All side=bottom should draw below baseline."""
        df = _make_snp_data().copy()
        df["side"] = "bottom"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_caterpillar_with_scores(self):
        """Caterpillar with integer scores (Tanghulu mode)."""
        df = _make_snp_data().copy()
        df["score"] = [2, 3, 1, 4, 2]
        df["side"] = ["top", "bottom", "top", "bottom", "top"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())


# ---------------------------------------------------------------------------
# Tanghulu multi-shape tests
# ---------------------------------------------------------------------------

class TestTanghuluMultiShape:

    def test_list_valued_fill(self):
        """List-valued fill should alternate colours in stacked circles."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 3,
            "start": [100, 200, 300],
            "score": [3, 2, 4],
            "fill": [["#FF0000", "#0000FF"], ["#00FF00", "#FFFF00"],
                     ["#FF00FF", "#00FFFF"]],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 500)
        LolliplotTrack(df).draw(ax, GenomicInterval("chr1", 0, 500))

    def test_list_valued_shape(self):
        """List-valued shape should alternate shapes in stacked circles."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 2,
            "start": [100, 300],
            "score": [3, 2],
            "shape": [["circle", "square", "diamond"],
                      ["triangle_point_up", "circle"]],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 500)
        LolliplotTrack(df).draw(ax, GenomicInterval("chr1", 0, 500))


# ---------------------------------------------------------------------------
# Legend tests
# ---------------------------------------------------------------------------

class TestLegend:

    def test_legend_basic(self):
        """Basic legend with labels and fill."""
        df = _make_snp_data().copy()
        legend = {"labels": ["WT", "MUT"], "fill": ["#87CEFA", "#98CE31"]}
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       legend=legend).draw(ax, _make_region())

    def test_legend_with_ncol(self):
        """Legend with ncol parameter."""
        df = _make_snp_data().copy()
        legend = {"labels": ["A", "B", "C"], "fill": ["red", "blue", "green"],
                  "ncol": 3}
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       legend=legend).draw(ax, _make_region())

    def test_legend_empty(self):
        """Empty legend dict should not error."""
        df = _make_snp_data().copy()
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       legend={}).draw(ax, _make_region())

    def test_legend_via_lolliplot(self):
        """Legend parameter passed through lolliplot()."""
        legend = {"labels": ["X"], "fill": ["#FF0000"]}
        ax = lolliplot(_make_snp_data(), features=_make_feature_data(),
                       legend=legend)
        assert ax is not None


# ---------------------------------------------------------------------------
# Jitter, yaxis, ylab tests
# ---------------------------------------------------------------------------

class TestJitterYaxisYlab:

    def test_jitter_label(self):
        """jitter='label' should align labels at top with anti-overlap."""
        df = _make_snp_data().copy()
        df["score"] = [2, 3, 1, 4, 2]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       jitter="label").draw(ax, _make_region())

    def test_jitter_label_close_positions(self):
        """Close SNPs should be staggered to different rows."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 5,
            "start": [100, 101, 102, 500, 501],
            "score": [2, 3, 1, 2, 4],
            "label": ["a", "b", "c", "d", "e"],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 700)
        LolliplotTrack(df, jitter="label").draw(
            ax, GenomicInterval("chr1", 0, 700))

    def test_jitter_label_with_per_snp_rotation(self):
        """jitter='label' with custom per-SNP label rotation."""
        df = _make_snp_data().copy()
        df["score"] = [2, 3, 1, 4, 2]
        df["label_rotation"] = [0, 45, 90, 45, 0]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       jitter="label").draw(ax, _make_region())

    def test_jitter_label_with_per_snp_dashline(self):
        """jitter='label' with per-SNP dashline_col for connector lines."""
        df = _make_snp_data().copy()
        df["score"] = [2, 3, 1, 4, 2]
        df["dashline_col"] = ["red", "blue", "green", "orange", "purple"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       jitter="label").draw(ax, _make_region())

    def test_custom_yaxis(self):
        """Custom yaxis list should set tick positions."""
        df = _make_snp_data().copy()
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       yaxis=[0, 5, 10]).draw(ax, _make_region())

    def test_ylab(self):
        """ylab parameter should draw y-axis label."""
        df = _make_snp_data().copy()
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       ylab="# evidences").draw(ax, _make_region())

    def test_yaxis_and_ylab_combined(self):
        """yaxis + ylab together."""
        df = _make_snp_data().copy()
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data(),
                       yaxis=[0, 5, 10], ylab="Count").draw(ax, _make_region())


# ---------------------------------------------------------------------------
# pie.stack type tests
# ---------------------------------------------------------------------------

class TestPieStack:

    def test_pie_stack_basic(self):
        """pie.stack type should render stacked pies."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 4,
            "start": [100, 100, 300, 300],
            "score": [1, 1, 1, 1],
            "pie_values": [[3, 7], [5, 5], [2, 8], [6, 4]],
            "pie_colors": [["#0080FF", "#E69F00"]] * 4,
            "stack_factor": [1, 2, 1, 2],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 500)
        LolliplotTrack(df, type="pie.stack").draw(
            ax, GenomicInterval("chr1", 0, 500))

    def test_pie_stack_single_position(self):
        """pie.stack with all SNPs at same position."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 3,
            "start": [100, 100, 100],
            "score": [1, 1, 1],
            "pie_values": [[3, 7], [5, 5], [2, 8]],
            "pie_colors": [["red", "blue"]] * 3,
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 300)
        LolliplotTrack(df, type="pie.stack").draw(
            ax, GenomicInterval("chr1", 0, 300))


# ---------------------------------------------------------------------------
# Rescale tests
# ---------------------------------------------------------------------------

class TestRescale:

    def test_rescale_basic(self):
        """rescale parameter should remap coordinates."""
        df = pd.DataFrame({
            "chrom": ["chr1"] * 3,
            "start": [100, 500, 900],
            "score": [3, 5, 2],
        })
        rescale_map = [(0, 1000, 0, 500)]
        fig, ax = plt.subplots(); ax.set_xlim(0, 600)
        LolliplotTrack(df, rescale=rescale_map).draw(
            ax, GenomicInterval("chr1", 0, 600))

    def test_rescale_position_function(self):
        """rescale_position should remap correctly."""
        from geneview.genometracks._mutation_features import rescale_position
        rmap = [(0, 100, 0, 50), (200, 300, 100, 150)]
        assert rescale_position(50, rmap) == 25.0
        assert rescale_position(250, rmap) == 125.0
        # Outside all ranges
        assert rescale_position(150, rmap) == 150

    def test_rescale_none(self):
        """rescale=None should be identity."""
        from geneview.genometracks._mutation_features import rescale_position
        assert rescale_position(42, None) == 42


# ---------------------------------------------------------------------------
# Multi-layer features tests
# ---------------------------------------------------------------------------

class TestMultiLayerFeatures:

    def test_multi_layer_features(self):
        """feature_layer_id should stack features on separate baselines."""
        feats = pd.DataFrame({
            "chrom": ["chr1"] * 4,
            "start": [1, 501, 1, 501],
            "end": [120, 900, 120, 900],
            "name": ["A1", "A2", "B1", "B2"],
            "feature_layer_id": [1, 1, 2, 2],
        })
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500); ax.set_ylim(0, 1)
        from geneview.genometracks._mutation_features import draw_features
        h = draw_features(ax, feats, 0, 1500)
        assert h > 0

    def test_multi_layer_with_lolliplot(self):
        """Multi-layer features via LolliplotTrack."""
        feats = _make_feature_data().copy()
        feats["feature_layer_id"] = [1, 1, 2]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(_make_snp_data(), features=feats).draw(ax, _make_region())


# ---------------------------------------------------------------------------
# New shapes tests
# ---------------------------------------------------------------------------

class TestNewShapes:

    def test_square_shape(self):
        """shape='square' should render squares."""
        df = _make_snp_data().copy()
        df["shape"] = "square"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_diamond_shape(self):
        """shape='diamond' should render diamonds."""
        df = _make_snp_data().copy()
        df["shape"] = "diamond"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_triangle_up_shape(self):
        df = _make_snp_data().copy()
        df["shape"] = "triangle_point_up"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_triangle_down_shape(self):
        df = _make_snp_data().copy()
        df["shape"] = "triangle_point_down"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_per_snp_shape(self):
        """Different shapes per SNP."""
        df = _make_snp_data().copy()
        df["shape"] = ["circle", "square", "diamond",
                       "triangle_point_up", "triangle_point_down"]
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        LolliplotTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_draw_shape_dispatcher(self):
        """draw_shape dispatcher should call correct functions."""
        from geneview.genometracks._mutation_shapes import draw_shape
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000); ax.set_ylim(0, 1)
        fig.canvas.draw()
        for name in ["circle", "square", "diamond",
                     "triangle_point_up", "triangle_point_down"]:
            draw_shape(ax, name, 500, 0.5, 0.02, facecolor="red")

    def test_draw_node_label(self):
        """draw_node_label should render text."""
        from geneview.genometracks._mutation_shapes import draw_node_label
        fig, ax = plt.subplots()
        ax.set_xlim(0, 1000); ax.set_ylim(0, 1)
        draw_node_label(ax, 500, 0.5, 0.02, "X", color="white", fontsize=6)


# ---------------------------------------------------------------------------
# DandelionTrack new feature tests
# ---------------------------------------------------------------------------

class TestDandelionNewFeatures:

    def test_dandelion_rescale(self):
        """DandelionTrack with rescale parameter."""
        rmap = [(0, 1500, 0, 750)]
        fig, ax = plt.subplots(); ax.set_xlim(0, 800)
        DandelionTrack(
            _make_dense_snp_data(), features=_make_feature_data(),
            rescale=rmap,
        ).draw(ax, GenomicInterval("chr1", 0, 800))

    def test_dandelion_per_snp_cex(self):
        """DandelionTrack with per-SNP cex."""
        df = _make_dense_snp_data().copy()
        df["cex"] = np.random.uniform(0.5, 2.0, len(df))
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_dandelion_per_snp_label_rotation(self):
        """DandelionTrack with per-SNP label_rotation."""
        df = _make_dense_snp_data().copy()
        df["label_rotation"] = 45
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(df, features=_make_feature_data()).draw(ax, _make_region())

    def test_dandelion_per_snp_dashline_col(self):
        """DandelionTrack with per-SNP dashline_col."""
        df = _make_dense_snp_data().copy()
        df["dashline_col"] = "red"
        fig, ax = plt.subplots(); ax.set_xlim(0, 1500)
        DandelionTrack(df, features=_make_feature_data()).draw(ax, _make_region())
