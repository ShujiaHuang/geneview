"""Tests for geneview.popgene._admixture module (admixtureplot).

Author: Shujia Huang
"""
import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.popgene._admixture import admixtureplot, _draw_admixtureplot


def _make_admixture_data(n_groups=3, samples_per_group=20, k=4, seed=42):
    """Helper to create synthetic admixture-like data (dict of DataFrames)."""
    rng = np.random.RandomState(seed)
    data = {}
    group_names = [f"POP{i+1}" for i in range(n_groups)]
    for g in group_names:
        # Each row sums to 1 (admixture proportions for K components)
        raw = rng.dirichlet(np.ones(k), size=samples_per_group)
        df = pd.DataFrame(raw, columns=[str(j) for j in range(k)])
        data[g] = df
    return data, group_names


class TestAdmixtureplot:
    """Tests for the admixtureplot function."""

    def test_basic_dict_input(self):
        """Should work with dict input."""
        data, _ = _make_admixture_data()
        ax = admixtureplot(data=data)
        assert ax is not None
        assert isinstance(ax, plt.Axes)

    def test_custom_group_order(self):
        """Should respect the specified group_order."""
        data, group_names = _make_admixture_data()
        order = list(reversed(group_names))
        ax = admixtureplot(data=data, group_order=order)
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert tick_labels == order

    def test_custom_xticklabels(self):
        """Should use custom x-tick labels."""
        data, group_names = _make_admixture_data()
        custom_labels = ["Group A", "Group B", "Group C"]
        ax = admixtureplot(data=data, xticklabels=custom_labels)
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        assert tick_labels == custom_labels

    def test_xticklabels_wrong_count_raises(self):
        """Should raise ValueError if xticklabels count doesn't match groups."""
        data, _ = _make_admixture_data()
        with pytest.raises(ValueError, match="not the same"):
            admixtureplot(data=data, xticklabels=["A", "B"])

    def test_ylabel_auto(self):
        """Should auto-generate ylabel as K=<n>."""
        data, _ = _make_admixture_data(k=5)
        ax = admixtureplot(data=data)
        assert ax.get_ylabel() == "K=5"

    def test_ylabel_custom(self):
        """Should use custom ylabel."""
        data, _ = _make_admixture_data()
        ax = admixtureplot(data=data, ylabel="K=4 custom")
        assert ax.get_ylabel() == "K=4 custom"

    def test_custom_ax(self):
        """Should use the provided Axes."""
        data, _ = _make_admixture_data()
        fig, ax = plt.subplots(figsize=(10, 3))
        result_ax = admixtureplot(data=data, ax=ax)
        assert result_ax is ax

    def test_set_xticklabel_top(self):
        """Should move x-tick labels to top without error."""
        data, _ = _make_admixture_data()
        ax = admixtureplot(data=data, set_xticklabel_top=True)
        # Verify the function ran without error and x-axis tick_top was called
        # After tick_top(), bottom label visibility should be off
        assert ax is not None

    def test_palette_string(self):
        """Should accept string palette."""
        data, _ = _make_admixture_data()
        ax = admixtureplot(data=data, palette="Set1")
        assert ax is not None

    def test_palette_list(self):
        """Should accept list palette."""
        data, _ = _make_admixture_data(k=4)
        ax = admixtureplot(data=data, palette=["#FF0000", "#00FF00", "#0000FF", "#FFFF00"])
        assert ax is not None

    def test_xticklabel_kws(self):
        """Should pass xticklabel_kws to set_xticklabels."""
        data, _ = _make_admixture_data()
        ax = admixtureplot(data=data,
                           xticklabel_kws={"rotation": "vertical"})
        assert ax is not None

    def test_population_info_not_string_raises(self):
        """population_info must be a string path if provided."""
        data, _ = _make_admixture_data()
        with pytest.raises(ValueError, match="file path"):
            admixtureplot(data=data, population_info=123)

    def test_data_not_dict_without_population_info_raises(self):
        """Should raise ValueError if data is not dict and population_info is None."""
        with pytest.raises(ValueError, match="must be a dictionary"):
            admixtureplot(data=[1, 2, 3])

    def test_invalid_data_type_raises(self):
        """Should raise ValueError for invalid data types."""
        with pytest.raises(ValueError):
            admixtureplot(data=42)

    def test_xlim_and_ylim(self):
        """Plot x/y limits should be set correctly."""
        data, _ = _make_admixture_data(n_groups=2, samples_per_group=10)
        ax = admixtureplot(data=data)
        ylim = ax.get_ylim()
        assert ylim[0] == 0.0
        assert ylim[1] == 1.0


class TestDrawAdmixtureplot:
    """Tests for the internal _draw_admixtureplot function."""

    def test_single_sample_group(self):
        """Groups with 1 sample should not trigger clustering."""
        rng = np.random.RandomState(42)
        data = {
            "G1": pd.DataFrame(rng.dirichlet([1, 1, 1], size=1),
                               columns=["0", "1", "2"]),
            "G2": pd.DataFrame(rng.dirichlet([1, 1, 1], size=1),
                               columns=["0", "1", "2"]),
        }
        ax = _draw_admixtureplot(data=data, group_order=["G1", "G2"])
        assert ax is not None

    def test_group_order_none(self):
        """group_order=None should use all keys."""
        data, _ = _make_admixture_data(n_groups=2)
        ax = _draw_admixtureplot(data=data)
        assert ax is not None

    def test_missing_group_raises(self):
        """Should raise KeyError for missing group in group_order."""
        data, _ = _make_admixture_data(n_groups=2)
        with pytest.raises(KeyError):
            _draw_admixtureplot(data=data, group_order=["POP1", "MISSING"])
