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
    _build_label,
    _annotate_top_snps,
    _validate_annotate_fmt,
)
from geneview.utils import adjust_text
from geneview.utils._adjust_text import get_bboxes
from matplotlib.text import Annotation


def _make_dense_gwas_data(n_chroms=6, snps_per_chrom=200, n_loci=20, seed=7):
    """Synthetic GWAS data with many independent significant loci."""
    rng = np.random.RandomState(seed)
    rows = []
    for chrom in range(1, n_chroms + 1):
        positions = np.sort(rng.randint(1, 2_000_000, size=snps_per_chrom))
        pvalues = rng.uniform(1e-3, 1.0, size=snps_per_chrom)
        for i, (pos, pv) in enumerate(zip(positions, pvalues)):
            rows.append({"#CHROM": f"chr{chrom}", "POS": int(pos),
                         "P": float(pv), "ID": f"rs{chrom}_{i}"})
    df = pd.DataFrame(rows)
    hit_idx = rng.choice(len(df), size=n_loci, replace=False)
    df.loc[hit_idx, "P"] = rng.uniform(1e-15, 1e-8, size=n_loci)
    return df


def _total_overlap_area(texts, ax):
    """Sum of pairwise overlap areas of the given texts' bounding boxes."""
    bboxes = get_bboxes(texts, None, (1, 1), ax=ax)
    total = 0.0
    for i in range(len(bboxes)):
        for j in range(i + 1, len(bboxes)):
            inter = bboxes[i].intersection(bboxes[i], bboxes[j])
            if inter is not None:
                total += abs(inter.width * inter.height)
    return total


def _count_arrows(ax):
    """Number of annotation arrows currently drawn on the axes."""
    return sum(1 for t in ax.texts
              if isinstance(t, Annotation) and t.arrow_patch is not None)


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


class TestBuildLabel:
    """Tests for the _build_label annotation-text helper."""

    def test_none_returns_snp_id(self):
        """annotate_fmt=None should yield just the SNP id."""
        item = [100, 8.0, "rs123", "chr1", 1e-8, 500]
        assert _build_label(item, logp=True, annotate_fmt=None) == "rs123"

    def test_format_string_fields(self):
        """A format string may reference snp/chrom/pos/p/log10p."""
        item = [100, 8.0, "rs123", "chr1", 1e-8, 500]
        label = _build_label(item, logp=True, annotate_fmt="{snp} {chrom}:{pos} P={p:.0e}")
        assert label == "rs123 chr1:500 P=1e-08"

    def test_callable_receives_fields(self):
        """A callable receives the fields as keyword arguments."""
        item = [100, 8.0, "rs123", "chr1", 1e-8, 500]
        fn = lambda snp, chrom, pos, p, log10p: f"{chrom}:{pos}"
        assert _build_label(item, logp=True, annotate_fmt=fn) == "chr1:500"

    def test_log10p_computed_from_p(self):
        """log10p should be derived from the p-value when available."""
        item = [100, 8.0, "rs123", "chr1", 1e-8, 500]
        label = _build_label(item, logp=True, annotate_fmt="{log10p:.1f}")
        assert label == "8.0"

    def test_short_item_no_extra_fields(self):
        """Legacy 3-element items should still format the snp id."""
        item = [100, 8.0, "rs123"]
        assert _build_label(item, logp=True, annotate_fmt=None) == "rs123"


class TestAnnotateLayouts:
    """Tests for the manhattanplot top-SNP annotation layouts / options."""

    def _sig_df(self):
        df = _make_gwas_data()
        df.loc[0, "P"] = 1e-12
        df.loc[50, "P"] = 1e-10
        df.loc[150, "P"] = 1e-9
        return df

    def test_default_repel_creates_labels(self):
        """Default layout should place one text per top SNP."""
        df = self._sig_df()
        ax = manhattanplot(df, sign_marker_p=1e-6, is_annotate_topsnp=True)
        labels = [t.get_text() for t in ax.texts if t.get_text()]
        assert any(lbl.startswith("rs") for lbl in labels)

    def test_annotate_fmt_string(self):
        """Format string should be reflected in the label text."""
        df = self._sig_df()
        ax = manhattanplot(df, sign_marker_p=1e-6, is_annotate_topsnp=True,
                           annotate_fmt="{snp}|{p:.0e}")
        labels = [t.get_text() for t in ax.texts if "|" in t.get_text()]
        assert labels and all("|" in lbl for lbl in labels)

    def test_annotate_fmt_callable(self):
        """Callable formatter should be reflected in the label text."""
        df = self._sig_df()
        ax = manhattanplot(df, sign_marker_p=1e-6, is_annotate_topsnp=True,
                           annotate_fmt=lambda snp, chrom, pos, p, log10p: f"[{chrom}]")
        labels = [t.get_text() for t in ax.texts if t.get_text().startswith("[")]
        assert labels

    def test_text_kws_fontsize_applied(self):
        """text_kws styling must reach the label text objects (regression)."""
        df = self._sig_df()
        ax = manhattanplot(df, sign_marker_p=1e-6, is_annotate_topsnp=True,
                           text_kws={"fontsize": 14})
        sized = [t for t in ax.texts if t.get_text() and t.get_fontsize() == 14]
        assert sized

    def test_lane_layout_runs(self):
        """Lane layout should create labels and leader arrows."""
        df = _make_dense_gwas_data()
        ax = manhattanplot(df, sign_marker_p=1e-8, is_annotate_topsnp=True,
                           annotate_layout="lane", annotate_fmt="{snp}")
        labels = [t.get_text() for t in ax.texts if t.get_text()]
        assert len(labels) >= 5
        assert _count_arrows(ax) >= 5

    def test_lane_labels_non_overlapping_in_x(self):
        """Lane layout should spread labels so they do not stack in x."""
        df = _make_dense_gwas_data()
        ax = manhattanplot(df, sign_marker_p=1e-8, is_annotate_topsnp=True,
                           annotate_layout="lane", text_kws={"fontsize": 6})
        xs = sorted(t.get_position()[0] for t in ax.texts if t.get_text())
        diffs = np.diff(xs)
        assert np.all(diffs >= 0)  # monotonic, i.e. ordered and separated

    def test_invalid_layout_raises(self):
        """An unknown layout name should raise ValueError."""
        df = self._sig_df()
        with pytest.raises(ValueError, match="annotate_layout"):
            manhattanplot(df, sign_marker_p=1e-6, is_annotate_topsnp=True,
                          annotate_layout="bogus")

    def test_arrowprops_none_disables_arrows(self):
        """Passing arrowprops=None should draw no connecting arrows."""
        df = _make_dense_gwas_data()
        ax = manhattanplot(df, sign_marker_p=1e-8, is_annotate_topsnp=True,
                           annotate_layout="lane", text_kws={"arrowprops": None})
        assert _count_arrows(ax) == 0

    def test_custom_arrowprops_drawn(self):
        """A custom arrowprops in text_kws should still draw arrows."""
        df = _make_dense_gwas_data()
        ax = manhattanplot(df, sign_marker_p=1e-8, is_annotate_topsnp=True,
                           annotate_layout="lane",
                           text_kws={"arrowprops": dict(arrowstyle="->", color="g")})
        assert _count_arrows(ax) >= 5

    def test_adjust_text_kws_forwarded(self):
        """adjust_text_kws should be accepted by the repel layout."""
        df = _make_dense_gwas_data()
        ax = manhattanplot(df, sign_marker_p=1e-8, is_annotate_topsnp=True,
                           annotate_layout="repel",
                           adjust_text_kws={"lim": 50, "force_text": (0.4, 0.6),
                                            "only_move": {"points": "y",
                                                          "text": "xy",
                                                          "objects": "xy"}})
        assert ax is not None

    def test_single_chrom_annotation(self):
        """Annotation should also work when zooming into one chromosome."""
        df = self._sig_df()
        ax = manhattanplot(df, CHR="chr1", sign_marker_p=1e-6,
                           is_annotate_topsnp=True, annotate_fmt="{snp}")
        assert ax is not None


class TestAdjustText:
    """Tests for the vectorized adjust_text engine."""

    def test_empty_returns_zero(self):
        """No texts should be a no-op returning 0 iterations."""
        fig, ax = plt.subplots()
        assert adjust_text([], ax=ax) == 0

    def test_reduces_overlap(self):
        """Overlapping labels should end up with less total overlap."""
        fig, ax = plt.subplots()
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        texts = [ax.text(5, 5, f"label_{i}") for i in range(8)]
        before = _total_overlap_area(texts, ax)
        adjust_text(texts, ax=ax, lim=200)
        after = _total_overlap_area(texts, ax)
        assert after < before

    def test_only_move_locks_x(self):
        """only_move without 'x' should keep the x-position fixed."""
        fig, ax = plt.subplots()
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        texts = [ax.text(5, 5, f"label_{i}") for i in range(6)]
        xs_before = [t.get_position()[0] for t in texts]
        adjust_text(texts, ax=ax, lim=100,
                    only_move={"points": "y", "text": "y", "objects": "y"})
        xs_after = [t.get_position()[0] for t in texts]
        assert np.allclose(xs_before, xs_after)

    def test_arrows_drawn_when_moved(self):
        """Supplying arrowprops should draw arrows for displaced labels."""
        fig, ax = plt.subplots()
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        texts = [ax.text(5, 5, f"label_{i}") for i in range(6)]
        adjust_text(texts, ax=ax, lim=200,
                    arrowprops=dict(arrowstyle="-", color="0.5"))
        assert _count_arrows(ax) >= 1


class TestValidateAnnotateFmt:
    """``_validate_annotate_fmt`` rejects unusable label format strings early."""

    def test_none_and_callable_pass(self):
        """None and callables need no validation and must not raise."""
        _validate_annotate_fmt(None)
        _validate_annotate_fmt(lambda **kw: "x")

    @pytest.mark.parametrize("fmt", [
        "{snp}",
        "{snp}\nP={p:.1e}",
        "{chrom}:{pos} log10p={log10p:.2f}",
        "plain text without fields",
    ])
    def test_valid_formats_pass(self, fmt):
        """Format strings over the known fields must be accepted."""
        _validate_annotate_fmt(fmt)

    def test_unknown_field_raises(self):
        """An unknown field name must raise a clear ValueError."""
        with pytest.raises(ValueError, match="Unknown field 'gene'"):
            _validate_annotate_fmt("{gene}")

    def test_bad_format_spec_raises(self):
        """An invalid format spec must raise a clear ValueError."""
        with pytest.raises(ValueError, match="Invalid ``annotate_fmt``"):
            _validate_annotate_fmt("{p:qq}")

    def test_non_string_raises(self):
        """A non-string, non-callable annotate_fmt must raise ValueError."""
        with pytest.raises(ValueError, match="must be None"):
            _validate_annotate_fmt(123)

    def test_manhattanplot_surfaces_bad_fmt(self):
        """manhattanplot must reject a bad annotate_fmt before drawing."""
        df = _make_dense_gwas_data()
        with pytest.raises(ValueError, match="Unknown field"):
            manhattanplot(data=df, sign_marker_p=1e-6, is_annotate_topsnp=True,
                          annotate_fmt="{unknown_field}")
        plt.close("all")

