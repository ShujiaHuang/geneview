__version__ = "0.8.2"

from .palette import *
from .utils import load_dataset, get_dataset_names
from .karyotype import karyoplot
from .baseplot import venn, generate_petal_labels
from .gwas import manhattanplot, qqplot, qqnorm
from .popgene import admixtureplot
from .genometracks import (
    plot_tracks, GenomeAxisTrack, AnnotationTrack,
    GeneRegionTrack, DataTrack, HighlightTrack, GenomicInterval,
    LolliplotTrack, DandelionTrack, lolliplot, dandelion_plot,
    plot_tracks_grid, plot_tracks_multi, add_panel_labels, save_figure,
)
from .mtdna import (
    mito_genome_map, heteroplasmy_scatter, heteroplasmy_heatmap,
    mito_coverage_plot, mito_copynumber_plot,
    read_mito_vcf, read_mito_copynumber, read_mito_coverage,
    get_mt_genes, gene_at, MT_LENGTH,
)
from .plotstyle import (
    apply_style, use_style, list_styles, PlotStyle, set_style, get_active_style,
)

# Apply the default geneview style at import time.
# This replaces the previous hardcoded rcParams block and sets up fonts,
# colours, export settings, etc.  Users can switch to a journal-specific
# style afterwards, e.g.  ``apply_style("nature")``.
apply_style("geneview")
