"""Tests for geneview.utils._misc module (is_numeric, etc.).

Author: Shujia Huang
"""
import pytest
import numpy as np

from geneview.utils import is_numeric


class TestIsNumeric:
    """Tests for is_numeric utility function."""

    # Positive cases
    @pytest.mark.parametrize("value", [
        0, 1, -1, 100,
        0.0, 1.5, -3.14, 1e10, 1e-10,
        "0", "1", "-2", "3.14", "-0.5", "1e5", "1E-3",
        ".5", "-.5",
        np.int32(5), np.float64(3.14), np.int64(-10),
    ])
    def test_numeric_values(self, value):
        """All these values should be recognized as numeric."""
        assert is_numeric(value) is True

    # Negative cases
    @pytest.mark.parametrize("value", [
        "a", "abc", "", " ", "1.2.3", "--5", "12a",
    ])
    def test_non_numeric_values(self, value):
        """All these values should NOT be recognized as numeric."""
        assert is_numeric(value) is False

    def test_none_raises_type_error(self):
        """None is not convertible to float, raises TypeError."""
        with pytest.raises(TypeError):
            is_numeric(None)

    def test_boolean_is_numeric(self):
        """Booleans can be cast to float, so is_numeric returns True."""
        assert is_numeric(True) is True
        assert is_numeric(False) is True

    def test_special_float_strings(self):
        """Test special float string representations."""
        assert is_numeric("inf") is True
        assert is_numeric("-inf") is True
        assert is_numeric("nan") is True

    def test_complex_types_raise_type_error(self):
        """Complex types raise TypeError (not caught by is_numeric)."""
        with pytest.raises(TypeError):
            is_numeric([1, 2])
        with pytest.raises(TypeError):
            is_numeric({"a": 1})
        with pytest.raises(TypeError):
            is_numeric((1, 2))
