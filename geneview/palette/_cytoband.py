"""
Canonical cytoband (chromosome ideogram band) colors.

Single source of truth for the ``gieStain`` -> color mapping shared by
:func:`geneview.karyotype.karyoplot` and
:class:`geneview.genometracks.IdeogramTrack`.

The values are pulled from the Circos karyotype palette (see ``_circos.py``)
so the two rendering paths always agree instead of maintaining two copies.
"""
from ._circos import circos

# gieStain codes used in UCSC/Gviz cytoband tables, ordered light -> dark.
_CYTOBAND_STAINS = (
    "gneg", "gpos25", "gpos33", "gpos50", "gpos66",
    "gpos75", "gpos100", "gpos", "acen", "stalk", "gvar",
)

# Canonical gieStain -> hex color map (derived from the Circos palette).
CYTOBAND_COLORS = {stain: circos[stain] for stain in _CYTOBAND_STAINS}

# Fallback color for undefined / unknown stains (Gviz uses light grey).
CYTOBAND_DEFAULT_COLOR = "#C8C8C8"


def get_cytoband_color(gie_stain, default=CYTOBAND_DEFAULT_COLOR):
    """Return the hex color for a cytoband ``gieStain`` code.

    Parameters
    ----------
    gie_stain : str
        A cytoband stain code such as ``"gneg"``, ``"gpos50"``, ``"acen"``.
    default : matplotlib color, optional
        Color to return for unknown / undefined stains. Defaults to a light
        grey (``CYTOBAND_DEFAULT_COLOR``).

    Returns
    -------
    str
        A hex color string.
    """
    return CYTOBAND_COLORS.get(str(gie_stain).strip(), default)
