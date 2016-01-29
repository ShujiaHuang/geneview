"""
Utility functionality (:mod:`geneview.util`)
============================================

This package provides general exception/warning definitions used throughout
geneview, as well as various utility functionality, including plotting and
unit-testing convenience functions.

"""
from ._misc import is_numeric
from ._plott import axlabel, despine, offset_spines

__all__ = ['is_numeric', 'axlabel', 'despine', 'offset_spines']
