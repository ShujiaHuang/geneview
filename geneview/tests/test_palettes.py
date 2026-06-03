"""Tests for geneview.palette module (color generation, xkcd_rgb, circos).

Author: Shujia Huang
"""
import pytest

from geneview.palette import generate_colors_palette, xkcd_rgb, circos


class TestGenerateColorsPalette:
    """Tests for generate_colors_palette function."""

    def test_returns_list_of_colors(self):
        """Should return a list of RGBA color tuples."""
        colors = generate_colors_palette(n_colors=5)
        assert isinstance(colors, list)
        assert len(colors) == 5

    def test_each_color_has_4_channels(self):
        """Each color should be an RGBA tuple/list with 4 values."""
        colors = generate_colors_palette(n_colors=3)
        for c in colors:
            assert len(c) == 4

    def test_alpha_parameter(self):
        """Alpha parameter should affect the alpha channel of all colors."""
        colors_opaque = generate_colors_palette(n_colors=3, alpha=1.0)
        colors_transparent = generate_colors_palette(n_colors=3, alpha=0.3)
        for c in colors_opaque:
            assert abs(c[3] - 1.0) < 1e-6
        for c in colors_transparent:
            assert abs(c[3] - 0.3) < 1e-6

    def test_different_colormaps(self):
        """Different colormaps should produce different color sets."""
        colors_viridis = generate_colors_palette(cmap="viridis", n_colors=5)
        colors_plasma = generate_colors_palette(cmap="plasma", n_colors=5)
        # They should not be identical
        assert colors_viridis != colors_plasma

    def test_list_cmap_passthrough(self):
        """When cmap is a list, colors should be derived from that list."""
        custom_colors = ["#FF0000", "#00FF00", "#0000FF"]
        result = generate_colors_palette(cmap=custom_colors, n_colors=3)
        assert len(result) == 3
        # Each color should be RGBA
        for c in result:
            assert len(c) == 4

    def test_list_cmap_alpha(self):
        """Alpha should apply when cmap is a list."""
        custom_colors = ["#FF0000", "#00FF00"]
        result = generate_colors_palette(cmap=custom_colors, n_colors=2, alpha=0.5)
        for c in result:
            assert abs(c[3] - 0.5) < 1e-6

    def test_n_colors_1(self):
        """Should work with n_colors=1."""
        colors = generate_colors_palette(n_colors=1)
        assert len(colors) == 1

    def test_n_colors_large(self):
        """Should work with a large number of colors."""
        colors = generate_colors_palette(n_colors=50)
        assert len(colors) == 50


class TestXkcdRgb:
    """Tests for the xkcd_rgb color dictionary."""

    def test_is_dict(self):
        assert isinstance(xkcd_rgb, dict)

    def test_has_entries(self):
        assert len(xkcd_rgb) > 0

    def test_values_are_hex_strings(self):
        """Each value should be a hex color string."""
        for name, color in xkcd_rgb.items():
            assert isinstance(color, str)
            assert color.startswith("#"), f"Color '{name}' value '{color}' doesn't start with #"


class TestCircos:
    """Tests for the circos color dictionary."""

    def test_is_dict(self):
        assert isinstance(circos, dict)

    def test_has_gie_stain_keys(self):
        """Should contain standard G-banding stain levels."""
        expected_keys = ["gneg", "gpos25", "gpos50", "gpos75", "gpos100",
                         "acen", "gvar", "stalk"]
        for key in expected_keys:
            assert key in circos, f"Missing expected key: {key}"

    def test_values_are_hex_strings(self):
        """Each value should be a hex color string."""
        for name, color in circos.items():
            assert isinstance(color, str)
            assert color.startswith("#"), f"Circos color '{name}' value '{color}' doesn't start with #"
