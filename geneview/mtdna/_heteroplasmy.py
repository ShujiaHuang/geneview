"""Heteroplasmy landscape plots for mtDNA variants.

Two complementary views of the single most important mtDNA measurement, the
heteroplasmy fraction (VAF):

* :func:`heteroplasmy_scatter` — a linear position-vs-VAF scatter across the
  16,569 bp genome, with an optional gene strip and heteroplasmy threshold
  line.  Ideal for one sample or a modest cohort.
* :func:`heteroplasmy_heatmap` — a samples x variant-sites matrix of VAFs.
  Ideal for spotting shared vs private / low-level heteroplasmy in a cohort.

Author: Shujia Huang
"""
import numpy as np

from .._core.decorators import styled_plot
from ..plotstyle import get_active_style
from ._reference import MT_LENGTH, MT_FEATURE_LABELS, get_mt_genes
from ._utils import resolve_feature_colors, coerce_variant_frame

__all__ = ["heteroplasmy_scatter", "heteroplasmy_heatmap"]

# Colours for the categorical hue encodings shared by the scatter.
_STATUS_COLORS = {"HET": "#C44E52", "HOM": "#4C72B0"}
_VARTYPE_COLORS = {"SNV": "#4C72B0", "INS": "#55A868", "DEL": "#DD8452",
                   "MNV": "#8172B3", "REF": "#999999"}


def _sample_palette(labels):
    """Return a ``{label: color}`` map, following the active style palette."""
    labels = list(dict.fromkeys(labels))
    active = get_active_style()
    if active is not None and active.color_palette:
        base = list(active.color_palette)
    else:
        import matplotlib.pyplot as plt
        base = list(plt.rcParams["axes.prop_cycle"].by_key().get("color", []))
    if not base:
        base = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3"]
    return {lab: base[i % len(base)] for i, lab in enumerate(labels)}


def _draw_gene_band(ax, colors, y0, y1):
    """Draw the rCRS gene strip as coloured rectangles in ``[y0, y1]``."""
    import matplotlib.patches as mpatches
    for g in get_mt_genes():
        ax.add_patch(mpatches.Rectangle(
            (g["start"], y0), g["end"] - g["start"], y1 - y0,
            facecolor=colors.get(g["feature_type"], "#999999"),
            edgecolor="none", zorder=1))


@styled_plot(figsize=(12, 3.8), apply_spines=True)
def heteroplasmy_scatter(
    data,
    *,
    hue="feature",
    feature_colors=None,
    het_threshold=0.03,
    show_gene_band=True,
    marker_size=26,
    alpha=0.85,
    legend=True,
    vaf_col="vaf",
    pos_col="pos",
    ax=None,
    style=None,
):
    """Plot heteroplasmy fraction (VAF) against mtDNA position.

    Parameters
    ----------
    data : pandas.DataFrame
        Variant records, typically from :func:`geneview.mtdna.read_mito_vcf`.
        Must contain position and VAF columns.
    hue : {"feature", "status", "var_type", "sample", None}, optional
        How to colour the points.  ``"feature"`` (default) uses the gene
        feature type; ``"sample"`` colours by the ``sample`` column (cohort
        view); ``None`` uses a single colour.
    feature_colors : dict, optional
        Override feature-type colours (see
        :func:`geneview.mtdna._utils.resolve_feature_colors`).
    het_threshold : float, optional
        Draw a dashed horizontal line at this VAF to mark the
        heteroplasmy-calling threshold (``mitoquest caller -j``).  ``None``
        hides it.  Default 0.03.
    show_gene_band : bool, optional
        Draw the rCRS gene strip just below the axis.  Default True.
    marker_size : float, optional
        Marker area (pt^2).  Default 26.
    alpha : float, optional
        Marker alpha.  Default 0.85.
    legend : bool, optional
        Draw a legend for the hue encoding.  Default True.
    vaf_col, pos_col : str, optional
        Column names for the VAF and position.
    ax : matplotlib Axes, optional
        Target axes; created when omitted.
    style : str or PlotStyle, optional
        Journal style applied for the call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    df = coerce_variant_frame(data, vaf_col=vaf_col, pos_col=pos_col)
    colors = resolve_feature_colors(feature_colors)

    # Resolve the per-point colours and legend entries for the chosen hue.
    color_map, group_col, legend_title = _resolve_hue(df, hue, colors)

    band_lo, band_hi = -0.06, -0.015
    if show_gene_band:
        _draw_gene_band(ax, colors, band_lo, band_hi)
        y_bottom = band_lo - 0.01
    else:
        y_bottom = 0.0

    if group_col is None:
        ax.scatter(df["pos"], df["vaf"], s=marker_size, alpha=alpha,
                   color=next(iter(color_map.values())) if color_map else "#4C72B0",
                   edgecolor="white", linewidth=0.3, zorder=3)
    else:
        for key, sub in df.groupby(group_col):
            ax.scatter(sub["pos"], sub["vaf"], s=marker_size, alpha=alpha,
                       color=color_map.get(key, "#444444"), label=str(key),
                       edgecolor="white", linewidth=0.3, zorder=3)

    if het_threshold is not None:
        ax.axhline(het_threshold, color="#888888", lw=0.8, ls="--", zorder=2)

    ax.set_xlim(0, MT_LENGTH)
    ax.set_ylim(y_bottom, 1.02)
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Position on rCRS (bp)")
    ax.set_ylabel("Heteroplasmy fraction (VAF)")

    if legend and group_col is not None:
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=7, title=legend_title)
    return ax


def _resolve_hue(df, hue, feature_colors):
    """Return (``{key: color}``, group_column, legend_title) for *hue*."""
    if hue is None:
        return {"_": "#4C72B0"}, None, ""
    if hue == "status" and "status" in df.columns:
        return _STATUS_COLORS, "status", "Zygosity"
    if hue == "var_type" and "var_type" in df.columns:
        keys = list(dict.fromkeys(df["var_type"]))
        return {k: _VARTYPE_COLORS.get(k, "#444444") for k in keys}, "var_type", "Variant type"
    if hue == "sample" and "sample" in df.columns:
        return _sample_palette(df["sample"]), "sample", "Sample"
    # Default: feature type.
    labels = {"protein_coding", "rRNA", "tRNA", "control_region"}
    # Re-key by the display label so the legend reads nicely.
    df["_feat_label"] = df["feature_type"].map(MT_FEATURE_LABELS)
    color_by_label = {MT_FEATURE_LABELS[k]: feature_colors[k] for k in labels}
    return color_by_label, "_feat_label", "Feature"


@styled_plot(figsize=(12, 5), apply_spines=False)
def heteroplasmy_heatmap(
    data,
    *,
    cmap="rocket_r",
    site_label="pos",
    vmin=0.0,
    vmax=1.0,
    colorbar=True,
    grid=True,
    vaf_col="vaf",
    pos_col="pos",
    sample_col="sample",
    ax=None,
    style=None,
):
    """Draw a samples x variant-sites heatmap of heteroplasmy fractions.

    Parameters
    ----------
    data : pandas.DataFrame
        Long-format variant records (one row per sample x site x allele),
        typically from :func:`geneview.mtdna.read_mito_vcf`.
    cmap : str, optional
        Matplotlib/geneview colormap name for the VAF scale.  Default
        ``"rocket_r"`` (falls back to ``"Reds"`` if unavailable).
    site_label : {"pos", "variant"}, optional
        Column label style: bare position (default) or ``pos ref>alt``.
    vmin, vmax : float, optional
        Colour scale bounds.  Default 0..1.
    colorbar : bool, optional
        Draw a colorbar.  Default True.
    grid : bool, optional
        Draw thin cell grid lines.  Default True.
    vaf_col, pos_col, sample_col : str, optional
        Column names for VAF, position and sample.
    ax : matplotlib Axes, optional
        Target axes; created when omitted.
    style : str or PlotStyle, optional
        Journal style applied for the call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    df = coerce_variant_frame(data, vaf_col=vaf_col, pos_col=pos_col)
    if sample_col not in df.columns:
        raise ValueError("data must contain a '%s' column for the heatmap." % sample_col)

    # Build a site label column.
    if site_label == "variant" and {"ref", "alt"}.issubset(df.columns):
        df["_site"] = (df["pos"].astype(str) + " "
                       + df["ref"].astype(str) + ">" + df["alt"].astype(str))
    else:
        df["_site"] = df["pos"].astype(int)

    # Order sites by genomic position; samples by first appearance.
    site_order = (df.drop_duplicates("_site")
                    .sort_values("pos")["_site"].tolist())
    sample_order = list(dict.fromkeys(df[sample_col]))
    matrix = (df.pivot_table(index=sample_col, columns="_site",
                             values="vaf", aggfunc="max")
                .reindex(index=sample_order, columns=site_order))

    values = matrix.to_numpy(dtype=float)
    cmap_obj = _resolve_cmap(cmap)
    cmap_obj.set_bad(color="#f0f0f0")
    masked = np.ma.masked_invalid(values)

    mesh = ax.pcolormesh(masked, cmap=cmap_obj, vmin=vmin, vmax=vmax,
                         edgecolors="white" if grid else "none",
                         linewidth=0.4 if grid else 0)

    ax.set_xticks(np.arange(len(site_order)) + 0.5)
    ax.set_xticklabels([str(s) for s in site_order], rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(sample_order)) + 0.5)
    ax.set_yticklabels([str(s) for s in sample_order], fontsize=7)
    ax.set_xlabel("Variant site")
    ax.set_ylabel("Sample")
    ax.invert_yaxis()
    ax.set_xlim(0, len(site_order))
    ax.set_ylim(0, len(sample_order))

    if colorbar:
        cbar = ax.figure.colorbar(mesh, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("Heteroplasmy fraction (VAF)", fontsize=7)
    return ax


def _resolve_cmap(name):
    """Return a Colormap for *name*, tolerating geneview/seaborn names."""
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    try:
        return mpl.colormaps[name].copy()
    except (KeyError, AttributeError):
        pass
    # ``rocket_r`` etc. may not exist in bare matplotlib; degrade gracefully.
    fallback = "Reds" if name.endswith("_r") else "viridis"
    return plt.get_cmap(fallback).copy()
