from .palette import *
from .ext import *
from .io import *
from .baseplot import *
from .karyotype import *

from ._fastqplot import fqqualplot, fastqreport
from ._manhattan import manhattanplot
from ._qq import qqplot, qqnorm

# Set aesthetic parameters in one step and make the font size as 1.2 times
# in the default setting in ``plotting_context()`` in ``.palette._rcmod``
setup(font_scale=1.2)  # this function is in the ``.palette`` package

__version__ = '0.0.2'
