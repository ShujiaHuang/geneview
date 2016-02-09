from . import gwas
from . import util
from .palette import set, reset_default, reset_orig, axes_style, set_style, \
    plotting_context, set_context, set_palette, color_palette, hls_palette, \
    husl_palette, mpl_palette, dark_palette, light_palette, diverging_palette, \
    blend_palette, xkcd_palette, crayon_palette, cubehelix_palette, \
    set_color_codes

from .ext.miscplot import puppyplot, palplot

# set the default aesthetics plot style and make the font size as 1.2 times of
# the default setting in ``plotting_context()``
set(font_scale=1.2)

__version__ = '0.0.1.dev7'
