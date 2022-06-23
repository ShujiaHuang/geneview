import matplotlib

from .palette import *
from .karyotype import karyoplot
from .baseplot import venn, generate_petal_labels
from .gwas import manhattanplot, qqplot, qqnorm
from .popgene import admixtureplot


matplotlib.rcParams['ps.fonttype']     = 42
matplotlib.rcParams['pdf.fonttype']    = 42
matplotlib.rcParams['font.sans-serif'] = ["Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']     = 'sans-serif'
