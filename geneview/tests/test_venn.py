"""Tests for geneview.baseplot._venn module (venn, generate_petal_labels, etc.).

Author: Shujia Huang
"""
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geneview.baseplot._venn import (
    venn,
    vennx,
    generate_petal_labels,
    generate_colors,
    is_valid_dataset_dict,
    is_already_venn_dataset,
    _generate_logics,
    _generate_petal_sizes,
    _resolve_palette,
    _draw_venn,
)
from geneview.baseplot._venn_proportional import (
    _lens_area,
    _solve_distance,
    _layout_venn2,
    _layout_venn3,
)
from geneview.plotstyle import use_style


class TestGenerateLogics:
    """Tests for _generate_logics helper."""

    def test_2_sets(self):
        logics = list(_generate_logics(2))
        assert logics == ["01", "10", "11"]

    def test_3_sets(self):
        logics = list(_generate_logics(3))
        assert len(logics) == 7  # 2^3 - 1
        assert "111" in logics
        assert "001" in logics

    def test_n_sets_count(self):
        """Should produce 2^n - 1 logics for n sets."""
        for n in range(2, 7):
            logics = list(_generate_logics(n))
            assert len(logics) == 2**n - 1

    def test_all_logics_have_correct_length(self):
        logics = list(_generate_logics(4))
        for logic in logics:
            assert len(logic) == 4


class TestGeneratePetalLabels:
    """Tests for generate_petal_labels function."""

    def test_basic_3_sets(self):
        datasets = [{"A", "B"}, {"B", "C"}, {"C", "D"}]
        labels = generate_petal_labels(datasets)
        assert isinstance(labels, dict)
        assert len(labels) == 7  # 2^3 - 1

    def test_format_size(self):
        datasets = [{"A", "B"}, {"B", "C"}]
        labels = generate_petal_labels(datasets, fmt="{size}")
        for key, val in labels.items():
            # Each value should be a string representation of an integer
            assert val.isdigit() or val == "0"

    def test_format_percentage(self):
        datasets = [{"A", "B", "C"}, {"B", "C", "D"}]
        labels = generate_petal_labels(datasets, fmt="{percentage:.1f}%")
        for key, val in labels.items():
            assert "%" in val

    def test_format_logic(self):
        datasets = [{"A"}, {"B"}]
        labels = generate_petal_labels(datasets, fmt="{logic}")
        for key, val in labels.items():
            assert val == key

    def test_disjoint_sets(self):
        """Disjoint sets should have zero intersections."""
        datasets = [{"A"}, {"B"}, {"C"}]
        labels = generate_petal_labels(datasets, fmt="{size}")
        # All intersection petals should be 0
        assert labels["011"] == "0"
        assert labels["101"] == "0"
        assert labels["110"] == "0"
        assert labels["111"] == "0"

    def test_identical_sets(self):
        """Identical sets should show all in the intersection."""
        datasets = [{"A", "B", "C"}, {"A", "B", "C"}, {"A", "B", "C"}]
        labels = generate_petal_labels(datasets, fmt="{size}")
        assert labels["111"] == "3"  # All elements in intersection
        # Unique petals should be 0
        assert labels["100"] == "0"
        assert labels["010"] == "0"
        assert labels["001"] == "0"

    def test_invalid_n_sets_too_few(self):
        with pytest.raises(ValueError, match="between 2 and 6"):
            generate_petal_labels([{"A"}])

    def test_invalid_n_sets_too_many(self):
        with pytest.raises(ValueError, match="between 2 and 6"):
            generate_petal_labels([set() for _ in range(8)])


class TestGenerateColors:
    """Tests for generate_colors function."""

    def test_returns_list(self):
        colors = generate_colors(n_colors=3)
        assert isinstance(colors, list)
        assert len(colors) == 3

    def test_n_colors_range(self):
        """n_colors must be between 2 and 6."""
        with pytest.raises(ValueError):
            generate_colors(n_colors=1)
        with pytest.raises(ValueError):
            generate_colors(n_colors=7)

    def test_n_colors_not_int(self):
        with pytest.raises(ValueError):
            generate_colors(n_colors=3.5)


class TestIsValidDatasetDict:
    """Tests for is_valid_dataset_dict function."""

    def test_valid(self):
        data = {"A": {1, 2}, "B": {3, 4}}
        assert is_valid_dataset_dict(data) is True

    def test_not_dict(self):
        assert is_valid_dataset_dict([1, 2, 3]) is False

    def test_values_not_sets(self):
        data = {"A": [1, 2], "B": [3, 4]}
        assert is_valid_dataset_dict(data) is False

    def test_empty_dict(self):
        """An empty dict should be valid (no non-set values)."""
        assert is_valid_dataset_dict({}) is True


class TestIsAlreadyVennDataset:
    """Tests for is_already_venn_dataset function."""

    def test_valid_precomputed(self):
        data = {"01": "3", "10": "5", "11": "2"}
        assert is_already_venn_dataset(data, ["A", "B"]) is True

    def test_valid_precomputed_no_names(self):
        """Should work when names is None (infers n_sets from keys)."""
        data = {"01": "3", "10": "5", "11": "2"}
        assert is_already_venn_dataset(data, None) is True

    def test_invalid_key_length(self):
        data = {"011": "3", "10": "5"}  # inconsistent key lengths
        assert is_already_venn_dataset(data, None) is False

    def test_not_dict(self):
        assert is_already_venn_dataset([1, 2], ["A"]) is False

    def test_non_string_values(self):
        data = {"01": 3, "10": "5", "11": "2"}
        assert is_already_venn_dataset(data, ["A", "B"]) is False

    def test_names_not_list(self):
        data = {"01": "3", "10": "5", "11": "2"}
        assert is_already_venn_dataset(data, "AB") is False

    def test_empty_dict_no_names(self):
        assert is_already_venn_dataset({}, None) is False


class TestVenn:
    """Tests for the venn() main entry point."""

    def test_basic_3_sets(self):
        data = {
            "Set A": {"A", "B", "D", "E"},
            "Set B": {"C", "F", "B", "G"},
            "Set C": {"J", "C", "K"},
        }
        ax = venn(data)
        assert ax is not None

    def test_with_names(self):
        data = {"A": {1, 2, 3}, "B": {2, 3, 4}, "C": {3, 4, 5}}
        ax = venn(data, names=["X", "Y", "Z"])
        assert ax is not None

    def test_precomputed_data(self):
        """Should accept pre-computed petal labels."""
        data = {"011": "ns", "101": "48", "110": "50", "111": "ns"}
        ax = venn(data, names=["a", "b", "g"])
        assert ax is not None

    def test_precomputed_data_no_names(self):
        """Pre-computed data should work without explicit names."""
        data = {"01": "3", "10": "5", "11": "2"}
        ax = venn(data)
        assert ax is not None

    def test_invalid_data_type(self):
        with pytest.raises(TypeError, match="dictionaries of sets"):
            venn([1, 2, 3])

    def test_custom_palette(self):
        data = {"A": {1, 2}, "B": {2, 3}}
        ax = venn(data, palette="plasma")
        assert ax is not None

    def test_custom_fmt(self):
        data = {"A": {1, 2, 3}, "B": {2, 3, 4}}
        ax = venn(data, fmt="{size}\n({percentage:.0f}%)")
        assert ax is not None

    def test_legend_use_petal_color(self):
        data = {"A": {1, 2}, "B": {2, 3}}
        ax = venn(data, legend_use_petal_color=True)
        assert ax is not None

    def test_legend_loc(self):
        data = {"A": {1, 2}, "B": {2, 3}}
        ax = venn(data, legend_loc="upper left")
        assert ax is not None

    def test_custom_ax(self):
        fig, ax = plt.subplots()
        data = {"A": {1, 2}, "B": {2, 3}}
        result_ax = venn(data, ax=ax)
        assert result_ax is ax

    def test_2_to_6_sets(self):
        """Venn should work for 2 through 6 sets."""
        from numpy.random import choice
        for n in range(2, 7):
            names = [chr(65 + i) for i in range(n)]
            data = {name: set(choice(500, 100, replace=False)) for name in names}
            ax = venn(data)
            assert ax is not None
            plt.close("all")


class TestVennx:
    """Tests for the vennx() function."""

    def test_basic(self):
        data = {"01": "3", "10": "5", "11": "2"}
        ax = vennx(data, names=["A", "B"])
        assert ax is not None

    def test_names_not_list_raises(self):
        data = {"01": "3", "10": "5", "11": "2"}
        with pytest.raises(ValueError, match="should be a list"):
            vennx(data, names="AB")

    def test_invalid_data_raises(self):
        data = {"01": 3, "10": 5}  # values not strings
        with pytest.raises(TypeError, match="not a dict"):
            vennx(data, names=["A", "B"])


class TestGeneratePetalSizes:
    """Tests for the numeric _generate_petal_sizes helper."""

    def test_matches_labels(self):
        datasets = [{1, 2, 3}, {2, 3, 4}]
        sizes = _generate_petal_sizes(datasets)
        assert sizes == {"10": 1, "01": 1, "11": 2}

    def test_three_sets_counts(self):
        datasets = [set(range(0, 50)), set(range(30, 80)), set(range(20, 60))]
        sizes = _generate_petal_sizes(datasets)
        # Total per set must equal the original set size.
        t0 = sizes["100"] + sizes["110"] + sizes["101"] + sizes["111"]
        assert t0 == 50


class TestResolvePalette:
    """Tests for palette resolution and PlotStyle fallback."""

    def test_explicit_palette_wins(self):
        colors = _resolve_palette("plasma", 3, 0.4)
        assert len(colors) == 3

    def test_default_fallback_no_style(self):
        colors = _resolve_palette(None, 2, 0.4)
        assert len(colors) == 2
        # Fresh lists (mutating one must not corrupt the module constant).
        colors[0][-1] = 1.0
        again = _resolve_palette(None, 2, 0.4)
        assert again[0][-1] != 1.0

    def test_active_style_palette_used(self):
        with use_style("nature"):
            colors = _resolve_palette(None, 3, 0.4)
        # Nature's Wong palette begins with black.
        assert tuple(round(v, 3) for v in colors[0][:3]) == (0.0, 0.0, 0.0)


class TestProportionalGeometry:
    """Tests for the area-proportional geometry solver."""

    def test_lens_area_disjoint(self):
        assert _lens_area(1.0, 1.0, 5.0) == 0.0

    def test_lens_area_contained(self):
        # Small circle fully inside the large one -> area = small circle area.
        import math
        assert _lens_area(2.0, 0.5, 1.0) == pytest.approx(math.pi * 0.25)

    def test_solve_distance_reproduces_overlap(self):
        r1, r2, target = 3.0, 2.5, 4.0
        d = _solve_distance(r1, r2, target)
        assert _lens_area(r1, r2, d) == pytest.approx(target, rel=1e-3)

    def test_layout_venn2_shapes(self):
        centers, radii = _layout_venn2({"10": 40, "01": 40, "11": 20})
        assert len(centers) == 2 and len(radii) == 2
        assert radii[0] == pytest.approx(radii[1])  # equal totals -> equal radii

    def test_layout_venn3_shapes(self):
        sizes = _generate_petal_sizes(
            [set(range(0, 50)), set(range(30, 80)), set(range(20, 60))])
        centers, radii = _layout_venn3(sizes)
        assert len(centers) == 3 and len(radii) == 3


class TestVennProportional:
    """Tests for venn(..., proportional=True)."""

    def test_proportional_2_sets(self):
        data = {"A": set(range(0, 60)), "B": set(range(40, 100))}
        ax = venn(data, proportional=True)
        assert ax is not None
        assert len(ax.patches) == 2

    def test_proportional_3_sets(self):
        data = {"A": set(range(0, 50)), "B": set(range(30, 80)), "C": set(range(20, 60))}
        ax = venn(data, proportional=True, legend_use_petal_color=True)
        assert ax is not None
        assert len(ax.patches) == 3

    def test_proportional_too_many_sets_warns_and_falls_back(self):
        data = {chr(65 + i): {i, i + 1} for i in range(4)}
        with pytest.warns(UserWarning, match="2 or 3 sets"):
            ax = venn(data, proportional=True)
        assert ax is not None

    def test_proportional_precomputed_warns_and_falls_back(self):
        data = {"01": "3", "10": "5", "11": "2"}
        with pytest.warns(UserWarning, match="raw sets"):
            ax = venn(data, names=["A", "B"], proportional=True)
        assert ax is not None

    def test_proportional_legend_loc(self):
        data = {"A": {1, 2, 3}, "B": {2, 3, 4}}
        ax = venn(data, proportional=True, legend_loc="upper left")
        assert ax is not None
