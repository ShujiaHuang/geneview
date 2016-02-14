from . import gwas
from . import util
from .palette import *

from .ext.miscplot import puppyplot, palplot

# set the default aesthetics plot style and make the font size as 1.2 times of
# the default setting in ``plotting_context()`` in ``.palette._rcmod``
set(font_scale=1.2) # this function is in the ``.palette`` package

__version__ = '0.0.1.dev7'
