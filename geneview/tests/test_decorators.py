"""
Author: Shujia Huag
Date: 2021-05-11 11:49:59
"""
import pytest

from ..utils import deprecate_positional_args


# This test was adapted from scikit-learn
# github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/tests/test_validation.py
def test_deprecate_positional_args_warns_for_function():
    @deprecate_positional_args
    def f1(a, b, *, c=1, d=1):
        return a, b, c, d

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variable as a keyword arg: c\."
    ):
        assert f1(1, 2, 3) == (1, 2, 3, 1)

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variables as keyword args: c, d\."
    ):
        assert f1(1, 2, 3, 4) == (1, 2, 3, 4)

    @deprecate_positional_args
    def f2(a=1, *, b=1, c=1, d=1):
        return a, b, c, d

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variable as a keyword arg: b\.",
    ):
        assert f2(1, 2) == (1, 2, 1, 1)

    # The * is placed before a keyword only argument without a default value
    @deprecate_positional_args
    def f3(a, *, b, c=1, d=1):
        return a, b, c, d

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variable as a keyword arg: b\.",
    ):
        assert f3(1, 2) == (1, 2, 1, 1)


def test_deprecate_positional_args_warns_for_class():
    class A1:
        @deprecate_positional_args
        def __init__(self, a, b, *, c=1, d=1):
            self.a = a, b, c, d

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variable as a keyword arg: c\."
    ):
        assert A1(1, 2, 3).a == (1, 2, 3, 1)

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variables as keyword args: c, d\."
    ):
        assert A1(1, 2, 3, 4).a == (1, 2, 3, 4)

    class A2:
        @deprecate_positional_args
        def __init__(self, a=1, b=1, *, c=1, d=1):
            self.a = a, b, c, d

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variable as a keyword arg: c\.",
    ):
        assert A2(1, 2, 3).a == (1, 2, 3, 1)

    with pytest.warns(
            FutureWarning,
            match=r"Pass the following variables as keyword args: c, d\.",
    ):
        assert A2(1, 2, 3, 4).a == (1, 2, 3, 4)
