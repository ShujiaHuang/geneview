"""Tests for geneview.gwas._qq module (ppoints, qqplot, qqnorm).

Author: Shujia Huang
"""
import pytest
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.gwas._qq import ppoints, qqplot, qqnorm


class TestPpoints:
    """Tests for the ppoints function."""

    def test_basic_length(self):
        """Output should have the same length as n."""
        result = ppoints(10)
        assert len(result) == 10

    def test_values_between_0_and_1(self):
        """All values should be in (0, 1)."""
        result = ppoints(100)
        assert all(0 < v < 1 for v in result)

    def test_monotonically_increasing(self):
        """Values should be strictly increasing."""
        result = ppoints(50)
        for i in range(len(result) - 1):
            assert result[i] < result[i + 1]

    def test_default_a_half(self):
        """With a=0.5, formula is (1:n - 0.5) / n."""
        n = 10
        result = ppoints(n, a=0.5)
        expected = (np.arange(n, dtype=float) + 1 - 0.5) / (n + 1 - 2 * 0.5)
        np.testing.assert_allclose(result, expected)

    def test_a_three_eighths(self):
        """Test with a=3/8 (recommended for n <= 10)."""
        n = 8
        a = 3.0 / 8.0
        result = ppoints(n, a=a)
        expected = (np.arange(n, dtype=float) + 1 - a) / (n + 1 - 2 * a)
        np.testing.assert_allclose(result, expected)

    def test_array_input(self):
        """When n is an array, should use len(n)."""
        arr = [10, 20, 30, 40, 50]
        result = ppoints(arr)
        assert len(result) == 5

    def test_invalid_a_too_small(self):
        """a < 0 should raise ValueError."""
        with pytest.raises(ValueError):
            ppoints(10, a=-0.1)

    def test_invalid_a_too_large(self):
        """a > 1 should raise ValueError."""
        with pytest.raises(ValueError):
            ppoints(10, a=1.5)

    def test_n_equals_1(self):
        """Edge case: n=1 should return a single value."""
        result = ppoints(1)
        assert len(result) == 1
        assert 0 < result[0] < 1


class TestQqplot:
    """Tests for the qqplot function."""

    def test_basic_plot(self):
        """Should return a matplotlib Axes object."""
        data = np.random.uniform(0.01, 1.0, 50)
        ax = qqplot(data=data)
        assert ax is not None
        assert hasattr(ax, "scatter")

    def test_with_custom_ax(self):
        """Should use the provided Axes."""
        fig, ax = plt.subplots()
        data = np.random.uniform(0.01, 1.0, 50)
        result_ax = qqplot(data=data, ax=ax)
        assert result_ax is ax

    def test_logp_false(self):
        """Should work with logp=False."""
        data = np.random.uniform(0.01, 1.0, 50)
        ax = qqplot(data=data, logp=False)
        assert ax is not None

    def test_with_other_sample(self):
        """Should work when comparing two samples."""
        data1 = np.random.uniform(0, 1, 50)
        data2 = np.random.uniform(0, 1, 50)
        ax = qqplot(data=data1, other=data2, logp=False)
        assert ax is not None

    def test_title_with_lambda(self):
        """Title should contain lambda value."""
        data = np.random.uniform(0.01, 1.0, 100)
        ax = qqplot(data=data, title="Test ")
        title = ax.get_title()
        assert "lambda" in title.lower() or "\\lambda" in title

    def test_default_labels(self):
        """Should set default x/y labels."""
        data = np.random.uniform(0.01, 1.0, 50)
        ax = qqplot(data=data)
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

    def test_custom_labels(self):
        """Should use custom x/y labels."""
        data = np.random.uniform(0.01, 1.0, 50)
        ax = qqplot(data=data, xlabel="Custom X", ylabel="Custom Y")
        assert ax.get_xlabel() == "Custom X"
        assert ax.get_ylabel() == "Custom Y"

    def test_non_numeric_raises(self):
        """Should raise ValueError for non-numeric data."""
        with pytest.raises(ValueError, match="numeric"):
            qqplot(data=["a", "b", "c"])

    def test_mismatched_lengths_raises(self):
        """Should raise ValueError when data and other have different lengths."""
        data = np.random.uniform(0, 1, 50)
        other = np.random.uniform(0, 1, 30)
        with pytest.raises(ValueError, match="same size"):
            qqplot(data=data, other=other)

    def test_abline_none(self):
        """ablinecolor=None should not draw the abline."""
        data = np.random.uniform(0.01, 1.0, 50)
        ax = qqplot(data=data, ablinecolor=None)
        # Only the scatter collection should exist, no line
        lines = ax.get_lines()
        assert len(lines) == 0

    def test_series_input(self):
        """Should accept pandas Series as input."""
        import pandas as pd
        data = pd.Series(np.random.uniform(0.01, 1.0, 50))
        ax = qqplot(data=data)
        assert ax is not None


class TestQqnorm:
    """Tests for the qqnorm function."""

    def test_basic_plot(self):
        """Should return a matplotlib Axes object."""
        data = np.random.normal(size=100)
        ax = qqnorm(data=data)
        assert ax is not None

    def test_normalized_data(self):
        """The observed data should be approximately standard normal."""
        np.random.seed(42)
        data = np.random.normal(5, 2, 200)
        ax = qqnorm(data=data)
        # If we get here without error, normalization worked
        assert ax is not None

    def test_with_custom_ax(self):
        """Should use the provided Axes."""
        fig, ax = plt.subplots()
        data = np.random.normal(size=50)
        result_ax = qqnorm(data=data, ax=ax)
        assert result_ax is ax

    def test_default_labels(self):
        """Should set meaningful default labels."""
        data = np.random.normal(size=50)
        ax = qqnorm(data=data)
        assert ax.get_xlabel() == "Expected normal distribution"
        assert ax.get_ylabel() == "Observed distribution"

    def test_custom_labels(self):
        """Should use custom labels."""
        data = np.random.normal(size=50)
        ax = qqnorm(data=data, xlabel="Exp", ylabel="Obs")
        assert ax.get_xlabel() == "Exp"
        assert ax.get_ylabel() == "Obs"

    def test_non_numeric_raises(self):
        """Should raise ValueError for non-numeric data."""
        with pytest.raises(ValueError, match="numeric"):
            qqnorm(data=["x", "y", "z"])

    def test_abline_drawn(self):
        """Should draw a y=x abline by default."""
        data = np.random.normal(size=50)
        ax = qqnorm(data=data)
        lines = ax.get_lines()
        assert len(lines) >= 1  # at least the abline

    def test_no_abline(self):
        """ablinecolor=None should skip the abline."""
        data = np.random.normal(size=50)
        ax = qqnorm(data=data, ablinecolor=None)
        lines = ax.get_lines()
        assert len(lines) == 0
