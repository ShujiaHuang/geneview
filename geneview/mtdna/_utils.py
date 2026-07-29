"""Internal helpers shared across the :mod:`geneview.mtdna` plots.

Keeps colour resolution and variant-table coercion in one place so the
individual plot modules stay focused on drawing.

Author: Shujia Huang
"""
from typing import Dict, Optional

import numpy as np
import pandas as pd

from ..plotstyle import get_active_style
from ._reference import MT_FEATURE_COLORS, gene_at

# Order in which feature types are mapped onto a style palette.
_FEATURE_ORDER = ("protein_coding", "rRNA", "tRNA", "control_region")


def resolve_feature_colors(feature_colors: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Resolve the per-feature-type fill colours.

    Priority: an explicit ``feature_colors`` dict (merged over the defaults) >
    the active :class:`~geneview.plotstyle.PlotStyle` palette (its first four
    colours mapped to protein-coding / rRNA / tRNA / control-region, so the
    map follows ``nature``/``science``/``cell`` themes) > the built-in
    :data:`~geneview.mtdna._reference.MT_FEATURE_COLORS`.

    Parameters
    ----------
    feature_colors : dict, optional
        Partial or full ``{feature_type: color}`` overrides.

    Returns
    -------
    dict
        A complete mapping for all four feature types.
    """
    colors = dict(MT_FEATURE_COLORS)
    active = get_active_style()
    if active is not None and len(active.color_palette) >= len(_FEATURE_ORDER):
        for ftype, col in zip(_FEATURE_ORDER, active.color_palette):
            colors[ftype] = col
    if feature_colors:
        colors.update(feature_colors)
    return colors


def coerce_variant_frame(variants, vaf_col="vaf", pos_col="pos") -> pd.DataFrame:
    """Return a tidy variant frame with guaranteed ``pos``/``vaf`` columns.

    Accepts either the long frame produced by
    :func:`geneview.mtdna.read_mito_vcf` or any DataFrame / mapping carrying a
    position column.  Adds a ``feature_type`` column (looked up from the rCRS
    gene map) when it is absent, so colouring by region always works.

    Parameters
    ----------
    variants : pandas.DataFrame or mapping
        Variant records.  Must contain a position column (``pos`` by default).
    vaf_col, pos_col : str
        Column names for the heteroplasmy fraction and position.

    Returns
    -------
    pandas.DataFrame
        A copy with normalised ``pos`` (int), ``vaf`` (float, NaN when
        missing) and ``feature_type`` columns.
    """
    df = pd.DataFrame(variants).copy()
    if pos_col not in df.columns:
        raise ValueError(
            "variants must contain a '%s' column (got: %s)."
            % (pos_col, ", ".join(map(str, df.columns)))
        )
    df["pos"] = df[pos_col].astype(int)
    if vaf_col in df.columns:
        df["vaf"] = pd.to_numeric(df[vaf_col], errors="coerce")
    else:
        df["vaf"] = np.nan
    if "feature_type" not in df.columns:
        df["feature_type"] = df["pos"].apply(
            lambda p: (gene_at(p) or {}).get("feature_type", "control_region")
        )
    return df
