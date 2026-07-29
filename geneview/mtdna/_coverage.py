"""Coverage-depth and copy-number plots for mtDNA analysis.

* :func:`mito_coverage_plot` — per-base / binned sequencing depth across the
  mitochondrial contig for one or many samples.  Depth is the QC backbone of
  any mtDNA call set and exposes uneven coverage / NUMT artefacts.
* :func:`mito_copynumber_plot` — per-sample mtDNA copy number with its 95%
  confidence interval, straight from ``mitoquest copynum``.

Author: Shujia Huang
"""
import numpy as np

from .._core.decorators import styled_plot
from ._reference import MT_LENGTH
from ._utils import resolve_feature_colors
from ._heteroplasmy import _sample_palette, _draw_gene_band

__all__ = ["mito_coverage_plot", "mito_copynumber_plot"]


@styled_plot(figsize=(12, 3.8), apply_spines=True)
def mito_coverage_plot(
    coverage,
    *,
    hue="sample",
    log=False,
    fill=True,
    show_gene_band=True,
    feature_colors=None,
    linewidth=1.0,
    legend=True,
    pos_col="pos",
    depth_col="depth",
    sample_col="sample",
    ax=None,
    style=None,
):
    """Plot mtDNA sequencing depth across the genome.

    Parameters
    ----------
    coverage : pandas.DataFrame
        Long-format depth records, typically from
        :func:`geneview.mtdna.read_mito_coverage`.  Must contain position and
        depth columns; a ``sample`` column enables multi-sample overlay.
    hue : {"sample", None}, optional
        Colour lines by sample (default) or use a single colour (``None``).
    log : bool, optional
        Use a logarithmic depth axis.  Default False.
    fill : bool, optional
        Fill the area under each single-sample curve.  Ignored (lines only)
        when more than one sample is present.  Default True.
    show_gene_band : bool, optional
        Draw the rCRS gene strip below the axis.  Default True.
    feature_colors : dict, optional
        Override feature-type colours for the gene band.
    linewidth : float, optional
        Depth line width.  Default 1.0.
    legend : bool, optional
        Draw a per-sample legend when colouring by sample.  Default True.
    pos_col, depth_col, sample_col : str, optional
        Column names for position, depth and sample.
    ax : matplotlib Axes, optional
        Target axes; created when omitted.
    style : str or PlotStyle, optional
        Journal style applied for the call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import pandas as pd

    df = pd.DataFrame(coverage).copy()
    for col in (pos_col, depth_col):
        if col not in df.columns:
            raise ValueError("coverage must contain a '%s' column." % col)

    has_samples = sample_col in df.columns and df[sample_col].nunique() > 0
    samples = list(dict.fromkeys(df[sample_col])) if has_samples else ["_all"]
    multi = len(samples) > 1
    palette = _sample_palette(samples) if (hue == "sample" and has_samples) else None

    colors = resolve_feature_colors(feature_colors)

    max_depth = float(np.nanmax(df[depth_col])) if len(df) else 1.0
    # Reserve a band below zero for the gene strip, scaled to the data.
    band_hi = 0.0
    band_lo = -0.06 * (max_depth if not log else 1.0)
    if log:
        band_lo, band_hi = None, None  # gene band disabled on log axes

    for s in samples:
        sub = (df[df[sample_col] == s] if has_samples else df).sort_values(pos_col)
        col = palette[s] if palette else "#4C72B0"
        ax.plot(sub[pos_col], sub[depth_col], color=col, lw=linewidth,
                label=str(s) if has_samples else None, zorder=3)
        if fill and not multi:
            ax.fill_between(sub[pos_col], sub[depth_col], color=col,
                            alpha=0.25, zorder=2)

    if log:
        ax.set_yscale("log")
    elif show_gene_band and band_lo is not None:
        _draw_gene_band(ax, colors, band_lo, band_hi)

    ax.set_xlim(0, MT_LENGTH)
    if not log:
        ax.set_ylim(band_lo if (show_gene_band and band_lo is not None) else 0,
                    max_depth * 1.05 if max_depth > 0 else 1.0)
    ax.set_xlabel("Position on rCRS (bp)")
    ax.set_ylabel("Depth")

    if legend and hue == "sample" and multi:
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=7, title="Sample")
    return ax


@styled_plot(figsize=(7, 4.5), apply_spines=True)
def mito_copynumber_plot(
    data,
    *,
    orient="v",
    color=None,
    sort=True,
    show_ci=True,
    baseline=None,
    sample_col="sample",
    cn_col="copy_number",
    ci_low_col="ci_low",
    ci_high_col="ci_high",
    ax=None,
    style=None,
):
    """Plot per-sample mtDNA copy number with 95% confidence intervals.

    Parameters
    ----------
    data : pandas.DataFrame
        Per-sample copy-number records, typically from
        :func:`geneview.mtdna.read_mito_copynumber`.
    orient : {"v", "h"}, optional
        Bar orientation: vertical (default) or horizontal (better for many
        samples).
    color : str, optional
        Bar colour.  Defaults to the first colour of the active style palette.
    sort : bool, optional
        Sort samples by copy number (descending).  Default True.
    show_ci : bool, optional
        Draw 95% CI error bars when the CI columns are present.  Default True.
    baseline : float, optional
        Draw a reference line at this copy number (e.g. a cohort median).
        ``None`` hides it.
    sample_col, cn_col, ci_low_col, ci_high_col : str, optional
        Column names.
    ax : matplotlib Axes, optional
        Target axes; created when omitted.
    style : str or PlotStyle, optional
        Journal style applied for the call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import pandas as pd

    df = pd.DataFrame(data).copy()
    for col in (sample_col, cn_col):
        if col not in df.columns:
            raise ValueError("data must contain a '%s' column." % col)
    if sort:
        df = df.sort_values(cn_col, ascending=(orient == "h")).reset_index(drop=True)

    if color is None:
        active = _first_palette_color()
        color = active or "#4C72B0"

    n = len(df)
    idx = np.arange(n)
    cn = df[cn_col].to_numpy(dtype=float)

    err = None
    if show_ci and {ci_low_col, ci_high_col}.issubset(df.columns):
        lo = df[ci_low_col].to_numpy(dtype=float)
        hi = df[ci_high_col].to_numpy(dtype=float)
        if np.isfinite(lo).any() and np.isfinite(hi).any():
            err = np.vstack([np.clip(cn - lo, 0, None), np.clip(hi - cn, 0, None)])

    if orient == "h":
        ax.barh(idx, cn, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        if err is not None:
            ax.errorbar(cn, idx, xerr=err, fmt="none", ecolor="#333333",
                        elinewidth=0.8, capsize=2, zorder=4)
        ax.set_yticks(idx)
        ax.set_yticklabels(df[sample_col].astype(str), fontsize=7)
        ax.set_xlabel("mtDNA copy number")
        ax.set_ylabel("Sample")
        if baseline is not None:
            ax.axvline(baseline, color="#888888", ls="--", lw=0.8, zorder=2)
    else:
        ax.bar(idx, cn, color=color, edgecolor="white", linewidth=0.5, zorder=3)
        if err is not None:
            ax.errorbar(idx, cn, yerr=err, fmt="none", ecolor="#333333",
                        elinewidth=0.8, capsize=2, zorder=4)
        ax.set_xticks(idx)
        ax.set_xticklabels(df[sample_col].astype(str), rotation=90, fontsize=7)
        ax.set_ylabel("mtDNA copy number")
        ax.set_xlabel("Sample")
        if baseline is not None:
            ax.axhline(baseline, color="#888888", ls="--", lw=0.8, zorder=2)

    return ax


def _first_palette_color():
    """Return the first colour of the active style palette, if any."""
    from ..plotstyle import get_active_style
    active = get_active_style()
    if active is not None and active.color_palette:
        return active.color_palette[0]
    return None
