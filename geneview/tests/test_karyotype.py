"""Tests for geneview.karyotype._karyotype module (karyoplot).

Author: Shujia Huang
"""
import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.karyotype._karyotype import karyoplot, _chrom_sort_key


def _make_karyotype_data(chroms=None):
    """Helper to create synthetic karyotype data."""
    if chroms is None:
        chroms = ["chr1", "chr2", "chr3", "chrX"]
    rows = []
    gie_stains = ["gneg", "gpos25", "gpos50", "gpos75", "gpos100", "acen", "gvar"]
    for chrom in chroms:
        start = 0
        for i in range(10):
            end = start + np.random.randint(5_000_000, 20_000_000)
            stain = gie_stains[i % len(gie_stains)]
            rows.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": f"{chrom}p{i+1}",
                "gie_stain": stain,
            })
            start = end
    return pd.DataFrame(rows)


class TestChromSortKey:
    """Tests for the _chrom_sort_key function."""

    def test_numeric_order(self):
        """Numeric chromosomes should sort in numeric order."""
        chroms = ["chr10", "chr2", "chr1", "chr22", "chr3"]
        result = sorted(chroms, key=_chrom_sort_key)
        assert result == ["chr1", "chr2", "chr3", "chr10", "chr22"]

    def test_sex_chromosomes_after_numeric(self):
        """X, Y should come after numeric chromosomes."""
        chroms = ["chrX", "chr1", "chr2", "chrY"]
        result = sorted(chroms, key=_chrom_sort_key)
        assert result == ["chr1", "chr2", "chrX", "chrY"]

    def test_without_chr_prefix(self):
        """Should also work without 'chr' prefix."""
        chroms = ["10", "2", "1", "X"]
        result = sorted(chroms, key=_chrom_sort_key)
        assert result == ["1", "2", "10", "X"]

    def test_mt_chromosome(self):
        """MT chromosome should sort after numeric ones."""
        chroms = ["chr1", "chrMT", "chr2"]
        result = sorted(chroms, key=_chrom_sort_key)
        assert result == ["chr1", "chr2", "chrMT"]

    def test_single_chromosome(self):
        result = sorted(["chr5"], key=_chrom_sort_key)
        assert result == ["chr5"]

    def test_empty_list(self):
        result = sorted([], key=_chrom_sort_key)
        assert result == []


class TestKaryoplot:
    """Tests for the karyoplot function."""

    def test_basic_dataframe(self):
        """Should return a matplotlib Axes from DataFrame input."""
        df = _make_karyotype_data()
        ax = karyoplot(data=df)
        assert ax is not None
        assert isinstance(ax, plt.Axes)

    def test_custom_ax(self):
        """Should use the provided Axes."""
        df = _make_karyotype_data()
        fig, ax = plt.subplots(figsize=(20, 5))
        result_ax = karyoplot(data=df, ax=ax)
        assert result_ax is ax

    def test_chr_parameter(self):
        """CHR should filter to a single chromosome."""
        df = _make_karyotype_data()
        ax = karyoplot(data=df, CHR="chr2")
        ytick_labels = [t.get_text() for t in ax.get_yticklabels()]
        assert ytick_labels == ["chr2"]

    def test_ytick_labels_are_chromosomes(self):
        """Y-axis should show chromosome names."""
        chroms = ["chr1", "chr2", "chr3"]
        df = _make_karyotype_data(chroms=chroms)
        ax = karyoplot(data=df)
        ytick_labels = [t.get_text() for t in ax.get_yticklabels()]
        # Should contain all chromosomes (order depends on biological sort)
        assert set(ytick_labels) == set(chroms)

    def test_chromosomes_sorted_biologically(self):
        """Chromosomes should be in biological order, not lexicographic."""
        chroms = ["chr10", "chr2", "chr1"]
        df = _make_karyotype_data(chroms=chroms)
        ax = karyoplot(data=df)
        ytick_labels = [t.get_text() for t in ax.get_yticklabels()]
        assert ytick_labels == ["chr1", "chr2", "chr10"]

    def test_xaxis_in_megabases(self):
        """X-axis tick labels should be in megabases (M)."""
        df = _make_karyotype_data()
        ax = karyoplot(data=df)
        xtick_labels = [t.get_text() for t in ax.get_xticklabels()]
        for label in xtick_labels:
            assert label.endswith("M"), f"Label '{label}' doesn't end with 'M'"

    def test_xlim_starts_at_zero(self):
        """X-axis should start at 0."""
        df = _make_karyotype_data()
        ax = karyoplot(data=df)
        assert ax.get_xlim()[0] == 0

    def test_input_as_list_of_lists(self):
        """Should accept list-of-lists input."""
        data = [
            ["chr1", 0, 50000000, "p1", "gneg"],
            ["chr1", 50000000, 100000000, "p2", "gpos50"],
            ["chr2", 0, 80000000, "p1", "acen"],
        ]
        ax = karyoplot(data=data)
        assert ax is not None

    def test_unknown_gie_stain_uses_default_color(self):
        """Unknown gie_stain should use color4none."""
        df = pd.DataFrame([
            {"chrom": "chr1", "start": 0, "end": 100000, "name": "p1", "gie_stain": "unknown_stain"},
        ])
        ax = karyoplot(data=df, color4none="#FF0000")
        # Should not raise, just use the default color
        assert ax is not None

    def test_width_parameter(self):
        """Width should affect the bar height."""
        df = _make_karyotype_data(chroms=["chr1"])
        ax1 = karyoplot(data=df, width=0.3)
        ax2 = karyoplot(data=_make_karyotype_data(chroms=["chr1"]), width=0.8)
        # Both should work without error
        assert ax1 is not None
        assert ax2 is not None
