__version__ = "0.7.0"

from .palette import *
from .utils import load_dataset, get_dataset_names
from .karyotype import karyoplot
from .baseplot import venn, generate_petal_labels
from .gwas import manhattanplot, qqplot, qqnorm
from .popgene import admixtureplot
from .genometracks import (
    plot_tracks, GenomeAxisTrack, AnnotationTrack,
    GeneRegionTrack, DataTrack, HighlightTrack, GenomicInterval,
)
from .plotstyle import apply_style, use_style, list_styles, PlotStyle

# Apply the default geneview style at import time.
# This replaces the previous hardcoded rcParams block and sets up fonts,
# colours, export settings, etc.  Users can switch to a journal-specific
# style afterwards, e.g.  ``apply_style("nature")``.
apply_style("geneview")
