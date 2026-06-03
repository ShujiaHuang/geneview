"""Tests for geneview.gwas._manhattan module (manhattanplot + helpers).

Author: Shujia Huang
"""
import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.gwas._manhattan import (
    manhattanplot,
    _find_top_snp,
    _sign_snp_regions,
    _find_SNPs_which_overlap_sign_neighbour_region,
)


def _make_gwas_data(n_chroms=3, snps_per_chrom=100, seed=42):
    """Helper to create synthetic GWAS-like data."""
    rng = np.random.RandomState(seed)
    rows = []
    for chrom in range(1, n_chroms + 1):
        chrom_name = f"chr{chrom}"
        positions = np.sort(rng.randint(1, 1_000_000, size=snps_per_chrom))
        pvalues = rng.uniform(1e-10, 1.0, size=snps_per_chrom)
        for i, (pos, pv) in enumerate(zip(positions, pvalues)):
            rows.append({
                "#CHROM": chrom_name,
                "POS": pos,
                "P": pv,
                "ID": f"rs{chrom}_{i}",
            })
    return pd.DataFrame(rows)


class TestManhattanplot:
    """Tests for the manhattanplot function."""

    def test_basic_plot(self):
        """Should return a matplotlib Axes without error."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df)
        assert ax is not None
        assert hasattr(ax, "scatter")

    def test_returns_axes(self):
        """Should create and return an Axes when ax=None."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df)
        assert isinstance(ax, plt.Axes)

    def test_custom_ax(self):
        """Should use the provided Axes."""
        df = _make_gwas_data()
        fig, ax = plt.subplots()
        result_ax = manhattanplot(data=df, ax=ax)
        assert result_ax is ax

    def test_title_and_labels(self):
        """Should set title and axis labels."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, title="GWAS", xlabel="Chrom", ylabel="P")
        assert ax.get_title() == "GWAS"
        assert ax.get_xlabel() == "Chrom"
        assert ax.get_ylabel() == "P"

    def test_xtick_label_set(self):
        """Should filter ticks to the specified set."""
        df = _make_gwas_data(n_chroms=5)
        xtick = {"chr1", "chr3", "chr5"}
        ax = manhattanplot(data=df, xtick_label_set=xtick)
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        for lbl in tick_labels:
            assert lbl in xtick

    def test_chr_parameter(self):
        """CHR should plot only one chromosome."""
        df = _make_gwas_data(n_chroms=5)
        ax = manhattanplot(data=df, CHR="chr2")
        # Should have data only from chr2
        tick_labels = [t.get_text() for t in ax.get_xticklabels()]
        # When CHR is set, ticks are raw positions, not chromosome labels
        assert ax.get_xlabel() == "Chromosome"

    def test_chr_and_xtick_label_set_raises(self):
        """Should raise when both CHR and xtick_label_set are set."""
        df = _make_gwas_data()
        with pytest.raises(ValueError, match="can't be set simultaneously"):
            manhattanplot(data=df, CHR="chr1", xtick_label_set={"chr1"})

    def test_logp_false(self):
        """Should work with logp=False."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, logp=False)
        assert ax is not None

    def test_significance_lines(self):
        """Should draw suggestive and genome-wide significance lines."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, suggestiveline=1e-5, genomewideline=5e-8)
        # The axhlines are stored as Lines2D in the axes
        hlines = [line for line in ax.get_lines()]
        assert len(hlines) >= 2  # at least 2 significance lines

    def test_no_significance_lines(self):
        """Setting lines to None should not draw them."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, suggestiveline=None, genomewideline=None)
        hlines = [line for line in ax.get_lines()]
        assert len(hlines) == 0

    def test_sign_marker_p(self):
        """Should mark significant SNPs with sign_marker_color."""
        df = _make_gwas_data()
        # Force some very significant SNPs
        df.loc[0, "P"] = 1e-12
        ax = manhattanplot(data=df, sign_marker_p=1e-6)
        assert ax is not None

    def test_annotate_topsnp(self):
        """Should annotate top SNPs when is_annotate_topsnp=True."""
        df = _make_gwas_data()
        df.loc[0, "P"] = 1e-12
        df.loc[1, "P"] = 1e-10
        ax = manhattanplot(data=df, sign_marker_p=1e-6,
                           is_annotate_topsnp=True, ld_block_size=50000)
        assert ax is not None

    def test_invalid_data_type(self):
        """Should raise ValueError for non-DataFrame input."""
        with pytest.raises(ValueError, match="pandas.DataFrame"):
            manhattanplot(data=[1, 2, 3])

    def test_missing_chrom_column(self):
        """Should raise ValueError for missing chromosome column."""
        df = pd.DataFrame({"POS": [1], "P": [0.5]})
        with pytest.raises(ValueError, match="Column"):
            manhattanplot(data=df)

    def test_missing_pos_column(self):
        """Should raise ValueError for missing position column."""
        df = pd.DataFrame({"#CHROM": ["chr1"], "P": [0.5]})
        with pytest.raises(ValueError, match="Column"):
            manhattanplot(data=df)

    def test_missing_p_column(self):
        """Should raise ValueError for missing p-value column."""
        df = pd.DataFrame({"#CHROM": ["chr1"], "POS": [1]})
        with pytest.raises(ValueError, match="Column"):
            manhattanplot(data=df)

    def test_color_string_split(self):
        """Comma-separated color string should produce alternating colors."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, color="#3B5488,#53BBD5")
        assert ax is not None

    def test_hline_kws(self):
        """Should accept hline_kws for line style."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df, hline_kws={"linestyle": "--", "lw": 1.3})
        assert ax is not None

    def test_spines_hidden(self):
        """Top and right spines should be invisible."""
        df = _make_gwas_data()
        ax = manhattanplot(data=df)
        assert not ax.spines["top"].get_visible()
        assert not ax.spines["right"].get_visible()


class TestFindTopSnp:
    """Tests for _find_top_snp helper function."""

    def test_single_block(self):
        """All SNPs within one LD block should return one top SNP."""
        data = [[100, 5.0, "rs1", "chr1"],
                [120, 3.0, "rs2", "chr1"],
                [140, 8.0, "rs3", "chr1"]]
        result = _find_top_snp(data, ld_block_size=50000, is_get_biggest=True)
        assert len(result) == 1
        assert result[0][1] == 8.0  # highest y-value

    def test_single_block_smallest(self):
        """With is_get_biggest=False, should pick the smallest y-value."""
        data = [[100, 5.0, "rs1", "chr1"],
                [120, 3.0, "rs2", "chr1"],
                [140, 8.0, "rs3", "chr1"]]
        result = _find_top_snp(data, ld_block_size=50000, is_get_biggest=False)
        assert len(result) == 1
        assert result[0][1] == 3.0

    def test_two_blocks(self):
        """SNPs far apart should form two blocks."""
        data = [[100, 5.0, "rs1", "chr1"],
                [200000, 9.0, "rs2", "chr1"]]
        result = _find_top_snp(data, ld_block_size=50000, is_get_biggest=True)
        assert len(result) == 2

    def test_cross_chromosome_boundary(self):
        """SNPs on different chromosomes should NOT be grouped."""
        data = [[99990, 7.0, "rsA", "chr1"],
                [100010, 6.0, "rsB", "chr2"]]
        result = _find_top_snp(data, ld_block_size=50000, is_get_biggest=True)
        assert len(result) == 2  # one per chromosome

    def test_last_block_uses_is_get_biggest(self):
        """The last block should respect is_get_biggest parameter."""
        data = [[100, 2.0, "rs1", "chr1"],
                [120, 9.0, "rs2", "chr1"]]
        result_biggest = _find_top_snp(data, ld_block_size=50000, is_get_biggest=True)
        result_smallest = _find_top_snp(data, ld_block_size=50000, is_get_biggest=False)
        assert result_biggest[0][1] == 9.0
        assert result_smallest[0][1] == 2.0

    def test_empty_input(self):
        """Empty input should return empty list."""
        result = _find_top_snp([], ld_block_size=50000)
        assert result == []

    def test_backward_compat_3_element_items(self):
        """Should still work with 3-element items (no chrom ID)."""
        data = [[100, 5.0, "rs1"],
                [120, 3.0, "rs2"],
                [140, 8.0, "rs3"]]
        result = _find_top_snp(data, ld_block_size=50000, is_get_biggest=True)
        assert len(result) == 1


class TestSignSnpRegions:
    """Tests for _sign_snp_regions helper."""

    def test_single_snp(self):
        """Single SNP should create one region."""
        data = [[100000, 5.0, "rs1"]]
        regions = _sign_snp_regions(data, ld_block_size=50000)
        assert len(regions) == 1
        assert regions[0][0] == 50000   # 100000 - 50000
        assert regions[0][1] == 150000  # 100000 + 50000

    def test_two_distant_snps(self):
        """Two far-apart SNPs should create two regions."""
        data = [[100000, 5.0, "rs1"], [500000, 3.0, "rs2"]]
        regions = _sign_snp_regions(data, ld_block_size=50000)
        assert len(regions) == 2

    def test_two_close_snps_merge(self):
        """Two close SNPs should merge into one region."""
        data = [[100000, 5.0, "rs1"], [120000, 3.0, "rs2"]]
        regions = _sign_snp_regions(data, ld_block_size=50000)
        assert len(regions) == 1

    def test_empty_input(self):
        """Empty input should return empty list."""
        regions = _sign_snp_regions([], ld_block_size=50000)
        assert regions == []


class TestFindSnpsOverlapRegion:
    """Tests for _find_SNPs_which_overlap_sign_neighbour_region."""

    def test_basic_overlap(self):
        """Should find indices of SNPs within regions."""
        regions = [[50000, 150000]]
        x = [10000, 80000, 100000, 200000]
        result = _find_SNPs_which_overlap_sign_neighbour_region(regions, x)
        assert result == [1, 2]

    def test_no_overlap(self):
        """No SNPs in region should return empty list."""
        regions = [[50000, 100000]]
        x = [10000, 20000, 200000]
        result = _find_SNPs_which_overlap_sign_neighbour_region(regions, x)
        assert result == []

    def test_multiple_regions(self):
        """Should find overlaps with multiple regions."""
        regions = [[10000, 30000], [80000, 120000]]
        x = [20000, 50000, 100000, 150000]
        result = _find_SNPs_which_overlap_sign_neighbour_region(regions, x)
        assert result == [0, 2]
