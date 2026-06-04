"""Tests for the geneview.plotstyle module.

Covers:
- PlotStyle dataclass and style registry
- apply_style / use_style (global + context manager)
- Each built-in style produces valid rcParams
- Plot functions accept ``style=`` without error
"""
import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.plotstyle import (
    PlotStyle,
    apply_style,
    use_style,
    get_style,
    list_styles,
    register_style,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_gwas_data(n_chroms=2, snps_per_chrom=50, seed=0):
    """Tiny synthetic GWAS data for smoke-testing plot functions."""
    rng = np.random.RandomState(seed)
    rows = []
    for chrom in range(1, n_chroms + 1):
        positions = np.sort(rng.randint(1, 500_000, size=snps_per_chrom))
        pvalues = rng.uniform(1e-8, 1.0, size=snps_per_chrom)
        for i, (pos, pv) in enumerate(zip(positions, pvalues)):
            rows.append({"#CHROM": f"chr{chrom}", "POS": pos, "P": pv,
                         "ID": f"rs{chrom}_{i}"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Registry tests
# ---------------------------------------------------------------------------

class TestRegistry:

    def test_builtin_styles_registered(self):
        """All four built-in styles must be present after import."""
        styles = list_styles()
        for name in ("geneview", "nature", "science", "cell"):
            assert name in styles, f"'{name}' not found in registry"

    def test_get_style_returns_plotstyle(self):
        for name in list_styles():
            s = get_style(name)
            assert isinstance(s, PlotStyle)

    def test_get_style_invalid_name_raises(self):
        with pytest.raises(ValueError, match="Unknown plot style"):
            get_style("nonexistent_style")

    def test_register_custom_style(self):
        custom = PlotStyle(name="_test_custom", description="tmp")
        register_style(custom)
        assert "_test_custom" in list_styles()
        assert get_style("_test_custom") is custom


# ---------------------------------------------------------------------------
# PlotStyle.to_rc_params()
# ---------------------------------------------------------------------------

class TestStyleToRcParams:

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_to_rc_params_returns_dict(self, name):
        style = get_style(name)
        params = style.to_rc_params()
        assert isinstance(params, dict)
        assert len(params) > 0

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_to_rc_params_contains_required_keys(self, name):
        style = get_style(name)
        params = style.to_rc_params()
        for key in ("font.family", "axes.titlesize", "pdf.fonttype",
                    "xtick.major.size", "figure.dpi"):
            assert key in params, f"Missing key '{key}' in {name} style"

    @pytest.mark.parametrize("name", ["geneview", "nature", "science", "cell"])
    def test_to_rc_params_values_are_valid(self, name):
        """Applying the rcParams dict must not raise."""
        style = get_style(name)
        params = style.to_rc_params()
        # Save and restore to avoid side-effects across tests
        saved = matplotlib.rcParams.copy()
        try:
            matplotlib.rcParams.update(params)
        finally:
            matplotlib.rcParams.update(saved)


# ---------------------------------------------------------------------------
# apply_style
# ---------------------------------------------------------------------------

class TestApplyStyle:

    def test_apply_style_sets_rcparams(self):
        saved = matplotlib.rcParams.copy()
        try:
            apply_style("nature")
            assert matplotlib.rcParams["axes.titlesize"] == 7.0
            assert matplotlib.rcParams["pdf.fonttype"] == 42
        finally:
            matplotlib.rcParams.update(saved)

    def test_apply_style_returns_style(self):
        result = apply_style("geneview")
        assert isinstance(result, PlotStyle)
        assert result.name == "geneview"

    def test_apply_style_invalid_raises(self):
        with pytest.raises(ValueError):
            apply_style("not_a_style")

    def test_apply_style_with_object(self):
        custom = PlotStyle(name="_apply_obj_test", font_size_title=99.0)
        saved = matplotlib.rcParams.copy()
        try:
            apply_style(custom)
            assert matplotlib.rcParams["axes.titlesize"] == 99.0
        finally:
            matplotlib.rcParams.update(saved)


# ---------------------------------------------------------------------------
# use_style (context manager)
# ---------------------------------------------------------------------------

class TestUseStyle:

    def test_use_style_restores_on_exit(self):
        saved_title_size = matplotlib.rcParams["axes.titlesize"]
        with use_style("nature"):
            assert matplotlib.rcParams["axes.titlesize"] == 7.0
        # After exiting the context, the original value must be restored.
        assert matplotlib.rcParams["axes.titlesize"] == saved_title_size

    def test_use_style_none_is_noop(self):
        """use_style(None) must not change anything."""
        before = dict(matplotlib.rcParams)
        with use_style(None) as s:
            assert s is None
        after = dict(matplotlib.rcParams)
        # rcParams should be identical (ignoring non-hashable values)
        for k, v in before.items():
            try:
                assert after[k] == v, f"rcParam '{k}' changed with use_style(None)"
            except (ValueError, TypeError):
                pass  # skip cycler objects etc.

    def test_use_style_invalid_raises(self):
        with pytest.raises(ValueError):
            with use_style("nope"):
                pass


# ---------------------------------------------------------------------------
# Plot function integration (smoke tests)
# ---------------------------------------------------------------------------

class TestPlotFunctionStyleParam:
    """Verify that each plot function accepts style= without crashing."""

    def setup_method(self):
        self.gwas = _make_gwas_data()

    def teardown_method(self):
        plt.close("all")

    def test_manhattanplot_style_nature(self):
        from geneview import manhattanplot
        ax = manhattanplot(data=self.gwas, style="nature")
        assert ax is not None

    def test_manhattanplot_style_none(self):
        from geneview import manhattanplot
        ax = manhattanplot(data=self.gwas, style=None)
        assert ax is not None

    def test_qqplot_style_science(self):
        from geneview import qqplot
        ax = qqplot(data=self.gwas["P"], style="science")
        assert ax is not None

    def test_qqnorm_style_cell(self):
        from geneview import qqnorm
        ax = qqnorm(data=self.gwas["P"].values, style="cell")
        assert ax is not None

    def test_venn_style_nature(self):
        from geneview import venn
        data = {"A": {1, 2, 3, 4}, "B": {3, 4, 5, 6}, "C": {4, 5, 6, 7}}
        ax = venn(data, style="nature")
        assert ax is not None

    def test_admixtureplot_style_geneview(self):
        from geneview import admixtureplot
        # Minimal admixture-like data
        rng = np.random.RandomState(0)
        df = pd.DataFrame(rng.dirichlet([1, 1, 1], size=20))
        data = {"pop1": df}
        ax = admixtureplot(data=data, style="geneview")
        assert ax is not None


# ---------------------------------------------------------------------------
# Style-specific assertions
# ---------------------------------------------------------------------------

class TestNatureStyle:

    def test_wong_palette_length(self):
        style = get_style("nature")
        assert len(style.color_palette) == 8

    def test_max_font_size(self):
        style = get_style("nature")
        assert style.font_size_title <= 8.0
        assert style.font_size_tick >= 5.0

    def test_no_grid(self):
        style = get_style("nature")
        assert style.axes_grid is False

    def test_spines(self):
        style = get_style("nature")
        assert style.axes_spines_top is False
        assert style.axes_spines_right is False


class TestScienceStyle:

    def test_figure_width(self):
        style = get_style("science")
        # Single-column ≈ 2.36 in
        assert 2.0 <= style.figure_figsize[0] <= 3.0

    def test_savefig_dpi(self):
        style = get_style("science")
        assert style.savefig_dpi >= 600


class TestCellStyle:

    def test_figure_width(self):
        style = get_style("cell")
        # Single-column ≈ 3.35 in
        assert 3.0 <= style.figure_figsize[0] <= 4.0

    def test_font_sizes(self):
        style = get_style("cell")
        assert style.font_size_label == 7.0
        assert style.font_size_tick >= 6.0


class TestGeneviewDefault:

    def test_legacy_font_list(self):
        style = get_style("geneview")
        assert "Arial" in style.font_sans_serif
        assert "DejaVu Sans" in style.font_sans_serif

    def test_palette_starts_with_legacy_colors(self):
        style = get_style("geneview")
        assert style.color_palette[0] == "#3B5488"
        assert style.color_palette[1] == "#53BBD5"

    def test_spines_hidden(self):
        style = get_style("geneview")
        assert style.axes_spines_top is False
        assert style.axes_spines_right is False
