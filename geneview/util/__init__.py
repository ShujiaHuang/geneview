"""
Utility functionality (:mod:`geneview.util`)
============================================

This package provides general exception/warning definitions used throughout
geneview, as well as various utility functionality, including plotting and
unit-testing convenience functions.

"""
from ._misc import chr_id_cmp, is_numeric, is_integer, iqr, _kde_support 
from ._palette import desaturate, set_hls_values, get_color_cycle
from ._dataset import get_dataset_names, load_dataset
from ._plott import axlabel, despine, offset_spines, set_hls_values 

__all__ = ['chr_id_cmp', 'is_numeric', 'is_integer', 'iqr', '_kde_support', 
           'axlabel', 'despine', 'offset_spines', 'desaturate', 'set_hls_values', 
           'get_color_cycle', 'get_dataset_names', 'load_dataset', 
           'set_hls_values']
