"""Tests for extended DataTrack plot types and export_tracks."""
import os
import tempfile

import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.genometracks._data_track import DataTrack, _PLOT_TYPES
from geneview.genometracks._base import GenomicInterval
from geneview.genometracks._track_plot import plot_tracks
from geneview.genometracks._export import export_tracks
from geneview.genometracks._annotation import AnnotationTrack
from geneview.genometracks._io import read_wig


def _make_data(n=30, chrom="chr7", seed=42):
    rng = np.random.RandomState(seed)
    starts = np.linspace(1000, 2000, n, dtype=int)
    return pd.DataFrame({
        "chrom": [chrom] * n,
        "start": starts,
        "end": starts + 50,
        "value": rng.randn(n).cumsum(),
    })


def _make_multi_sample_data(n=30, seed=42):
    """Multi-sample data for confint/average tests."""
    rng = np.random.RandomState(seed)
    starts = np.linspace(1000, 2000, n, dtype=int)
    return pd.DataFrame({
        "chrom": [f"chr7"] * n,
        "start": starts,
        "end": starts + 50,
        "sample_A": rng.randn(n).cumsum(),
        "sample_B": rng.randn(n).cumsum() + 2,
        "sample_C": rng.randn(n).cumsum() - 1,
    })


class TestDataTrackNewPlotTypes:

    def test_plot_types_include_new(self):
        for t in ("a", "confint", "smooth", "horizon", "g", "r"):
            assert t in _PLOT_TYPES

    def test_average_type(self):
        data = _make_multi_sample_data()
        dt = DataTrack(data, type="a")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_confint_type(self):
        data = _make_multi_sample_data()
        dt = DataTrack(data, type="confint")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_smooth_type(self):
        data = _make_data()
        dt = DataTrack(data, type="smooth")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_smooth_custom_span(self):
        data = _make_data()
        dt = DataTrack(data, type="smooth", smooth_span=0.5)
        assert dt.smooth_span == 0.5
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_horizon_type(self):
        data = _make_data()
        data["value"] = data["value"] - data["value"].mean()
        dt = DataTrack(data, type="horizon", baseline=0)
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_grid_type(self):
        data = _make_data()
        dt = DataTrack(data, type="g")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_regression_type(self):
        data = _make_data()
        dt = DataTrack(data, type="r")
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")


class TestDataTrackCompositeTypes:

    def test_composite_list_type(self):
        data = _make_data()
        dt = DataTrack(data, type=["line", "g"])
        assert dt._composite_types == ["line", "g"]
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_composite_boxplot_average_grid(self):
        data = _make_multi_sample_data()
        dt = DataTrack(data, type=["boxplot", "a", "g"])
        region = GenomicInterval("chr7", 1000, 2050)
        axes = plot_tracks([dt], region=region, figsize=(8, 3))
        assert len(axes) == 1
        plt.close("all")

    def test_invalid_composite_raises(self):
        with pytest.raises(ValueError, match="Plot type must be"):
            DataTrack(_make_data(), type=["line", "invalid_type"])


class TestExportTracks:

    def test_export_bed(self):
        data = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "start": [100, 500],
            "end": [200, 600],
            "name": ["feat1", "feat2"],
            "strand": ["+", "-"],
        })
        track = AnnotationTrack(data)
        with tempfile.NamedTemporaryFile(suffix=".bed", delete=False) as f:
            result = export_tracks(track, f.name, fmt="bed")
        assert os.path.exists(result)
        with open(result) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        os.unlink(result)

    def test_export_gff(self):
        data = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "name": ["feat1"],
        })
        track = AnnotationTrack(data)
        with tempfile.NamedTemporaryFile(suffix=".gff", delete=False) as f:
            result = export_tracks(track, f.name, fmt="gff")
        assert os.path.exists(result)
        with open(result) as fh:
            content = fh.read()
        assert "##gff-version" in content
        os.unlink(result)

    def test_export_bedgraph(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        with tempfile.NamedTemporaryFile(suffix=".bdg", delete=False) as f:
            result = export_tracks(dt, f.name, fmt="bedgraph")
        assert os.path.exists(result)
        os.unlink(result)

    def test_export_wig(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        with tempfile.NamedTemporaryFile(suffix=".wig", delete=False) as f:
            result = export_tracks(dt, f.name, fmt="wig")
        assert os.path.exists(result)
        with open(result) as fh:
            content = fh.read()
        assert "variableStep" in content
        os.unlink(result)

    def test_export_invalid_format(self):
        data = _make_data()
        dt = DataTrack(data, type="line")
        with pytest.raises(ValueError, match="Format must be"):
            export_tracks(dt, "/tmp/test.xyz", fmt="xyz")

    def test_export_dataframe(self):
        df = pd.DataFrame({
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "value": [1.5],
        })
        with tempfile.NamedTemporaryFile(suffix=".bdg", delete=False) as f:
            result = export_tracks(df, f.name, fmt="bedgraph")
        assert os.path.exists(result)
        os.unlink(result)


class TestReadWig:

    def test_read_fixedStep(self):
        content = (
            "fixedStep chrom=chr1 start=100 step=50 span=50\n"
            "1.0\n2.0\n3.0\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".wig", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_wig(f.name)
        os.unlink(f.name)
        assert len(df) == 3
        assert df["chrom"].iloc[0] == "chr1"
        assert df["start"].iloc[0] == 99  # 0-based
        assert "value" in df.columns

    def test_read_variableStep(self):
        content = (
            "variableStep chrom=chr1 span=25\n"
            "100\t1.5\n"
            "200\t2.5\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".wig", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_wig(f.name)
        os.unlink(f.name)
        assert len(df) == 2
        assert df["start"].iloc[0] == 99  # 0-based

    def test_read_wig_empty(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".wig", delete=False) as f:
            f.write("")
            f.flush()
            df = read_wig(f.name)
        os.unlink(f.name)
        assert len(df) == 0

    def test_read_wig_nrows(self):
        content = (
            "fixedStep chrom=chr1 start=100 step=10 span=10\n"
            "1.0\n2.0\n3.0\n4.0\n5.0\n"
        )
        with tempfile.NamedTemporaryFile(mode="w", suffix=".wig", delete=False) as f:
            f.write(content)
            f.flush()
            df = read_wig(f.name, nrows=3)
        os.unlink(f.name)
        assert len(df) == 3
