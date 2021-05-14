"""
Utility functionality (:mod:`geneview.utils`)
============================================

This package provides general exception/warning definitions used throughout
geneview, as well as various utility functionality, including plotting and
unit-testing convenience functions.

"""
from ._misc import is_numeric
from ._dataset import get_dataset_names, load_dataset
from ._adjust_text import adjust_text
from ._decorators import deprecate_positional_args

__all__ = ["is_numeric",
           "get_dataset_names",
           "load_dataset",
           "adjust_text",
           "deprecate_positional_args"]
