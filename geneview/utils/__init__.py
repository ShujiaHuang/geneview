"""
Utility functionality (:mod:`geneview.util`)
============================================

This package provides general exception/warning definitions used throughout
geneview, as well as various utility functionality, including plotting and
unit-testing convenience functions.

"""
from ._misc import is_numeric, is_integer, iqr, _kde_support, \
    categorical_order, freedman_diaconis_bins
from ._plott import General
from ._palette import desaturate, set_hls_values, get_color_cycle
from ._dataset import get_dataset_names, load_dataset
from ._adjust_text import adjust_text
from ._decorators import deprecate_positional_args

__all__ = ["is_numeric",
           "is_integer",
           "iqr",
           "_kde_support",
           "General",
           "desaturate",
           "get_dataset_names",
           "load_dataset",
           "set_hls_values",
           "categorical_order",
           "get_color_cycle",
           "freedman_diaconis_bins",
           "adjust_text",
           "deprecate_positional_args"]
