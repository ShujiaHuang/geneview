"""
Read more about seaborn color palette
http://stanford.edu/~mwaskom/software/seaborn-dev/tutorial/color_palettes.html
"""
from ._xkcd_rgb import xkcd_rgb
from ._crayons import crayons
from ._circos import circos

from ._rcmod import setup, reset_default, reset_orig, axes_style, set_style, \
    plotting_context, set_context, set_palette  

from ._palettes import color_palette, hls_palette, husl_palette, mpl_palette, \
    dark_palette, light_palette, diverging_palette, blend_palette, \
    xkcd_palette, crayon_palette, circos_palette, cubehelix_palette, \
    set_color_codes

