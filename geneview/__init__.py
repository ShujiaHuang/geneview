from . import gwas
from . import util
from .palette import *

from .ext.miscplot import puppyplot, palplot

# Set aesthetic parameters in one step and make the font size as 1.2 times
# in the default setting in ``plotting_context()`` in ``.palette._rcmod``
set_all(font_scale=1.2) # this function is in the ``.palette`` package

__version__ = '0.0.1.dev7'
