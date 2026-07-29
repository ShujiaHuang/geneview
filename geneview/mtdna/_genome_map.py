"""Circular map of the human mitochondrial genome (the flagship mtDNA plot).

Draws the 16,569 bp rCRS as a ring of 37 genes (plus the D-loop control
region) coloured by feature type, with optional variant *lollipops* pointing
inward whose length encodes the heteroplasmy fraction (VAF).  This is the
canonical way to communicate *where* on the mitochondrial genome a study's
variants fall and *how heteroplasmic* they are.

Author: Shujia Huang
"""
from typing import Dict, Optional

import numpy as np

from .._core.decorators import styled_plot
from ._reference import (
    MT_LENGTH, MT_FEATURE_LABELS, get_mt_genes,
)
from ._utils import resolve_feature_colors, coerce_variant_frame

__all__ = ["mito_genome_map"]


def _theta(pos):
    """Map a 1-based mtDNA position to a polar angle (radians)."""
    return 2.0 * np.pi * (np.asarray(pos, dtype=float) - 1.0) / MT_LENGTH


def _screen_angle_deg(theta):
    """Screen angle (deg, CCW from east) for a polar axes with N-zero, CW."""
    return np.degrees(np.pi / 2.0 - theta)


@styled_plot(figsize=(8, 8), apply_spines=False,
             subplot_kws={"subplot_kw": {"projection": "polar"}})
def mito_genome_map(
    variants=None,
    *,
    gene_label="large",
    color_by="feature",
    feature_colors: Optional[Dict[str, str]] = None,
    gene_ring=(0.86, 0.96),
    lollipop_base=0.84,
    lollipop_min=0.42,
    marker_size=28,
    tick_interval=2000,
    title=None,
    legend=True,
    vaf_col="vaf",
    pos_col="pos",
    ax=None,
    style=None,
):
    """Draw a circular rCRS map of the mitochondrial genome.

    Parameters
    ----------
    variants : pandas.DataFrame, optional
        Variant records to overlay as inward lollipops.  Accepts the tidy
        frame from :func:`geneview.mtdna.read_mito_vcf` (or any frame with a
        position column).  Lollipop length encodes ``vaf`` when present; when
        ``vaf`` is missing every stem uses the full length.  ``None`` draws
        the bare gene ring.
    gene_label : {"large", "all", "none"}, optional
        Which gene arcs to label.  ``"large"`` (default) labels protein-coding
        genes, rRNAs and the control region; ``"all"`` also labels tRNAs;
        ``"none"`` labels nothing.
    color_by : {"feature", "status", "var_type"}, optional
        How to colour lollipop markers.  ``"feature"`` (default) uses the gene
        feature type at each position; ``"status"`` uses HET/HOM; ``"var_type"``
        uses SNV/INS/DEL.
    feature_colors : dict, optional
        Override ``{feature_type: color}`` for the gene arcs (and, when
        ``color_by="feature"``, the markers).  See
        :func:`geneview.mtdna._utils.resolve_feature_colors` for resolution.
    gene_ring : tuple of float, optional
        ``(inner_radius, outer_radius)`` of the gene ring in axes units
        (0..1).  Default ``(0.86, 0.96)``.
    lollipop_base : float, optional
        Radius at which lollipop stems start (just inside the ring).
    lollipop_min : float, optional
        Radius reached by a VAF of 1.0 (stems point inward, so smaller =
        longer).  A VAF of 0 stays at ``lollipop_base``.
    marker_size : float, optional
        Base area (pt^2) of the lollipop head markers.
    tick_interval : int, optional
        Spacing (bp) between position tick labels around the ring. ``0`` hides
        them.  Default 2000.
    title : str, optional
        Text drawn in the centre of the ring.
    legend : bool, optional
        Draw a legend for the colour encoding.  Default True.
    vaf_col, pos_col : str, optional
        Column names for the VAF and position in *variants*.
    ax : matplotlib polar Axes, optional
        Target axes (must be polar).  Created automatically when omitted.
    style : str or PlotStyle, optional
        Journal style applied for the duration of the call.

    Returns
    -------
    matplotlib.axes.Axes
        The polar axes with the map drawn on it.
    """
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    colors = resolve_feature_colors(feature_colors)
    r_in, r_out = gene_ring

    # --- Gene arcs -------------------------------------------------------
    genes = get_mt_genes()
    label_types = {
        "large": {"protein_coding", "rRNA", "control_region"},
        "all": {"protein_coding", "rRNA", "control_region", "tRNA"},
        "none": set(),
    }.get(gene_label, {"protein_coding", "rRNA", "control_region"})

    for g in genes:
        t0, t1 = _theta(g["start"]), _theta(g["end"])
        width = (t1 - t0) if t1 >= t0 else (t1 - t0 + 2 * np.pi)
        center = t0 + width / 2.0
        ax.bar(center, r_out - r_in, width=width, bottom=r_in,
               color=colors.get(g["feature_type"], "#999999"),
               edgecolor="white", linewidth=0.4, align="center", zorder=2)

        if g["feature_type"] in label_types:
            rot = _screen_angle_deg(center) - 90.0
            # Flip labels on the lower half so none are upside down.
            if 90 < (rot % 360) < 270:
                rot += 180
            ax.text(center, r_out + 0.035, g["name"],
                    rotation=rot, rotation_mode="anchor",
                    ha="center", va="center", fontsize=6.5, zorder=5)

    # --- Position ticks --------------------------------------------------
    if tick_interval and tick_interval > 0:
        for bp in range(0, MT_LENGTH, int(tick_interval)):
            t = _theta(bp if bp > 0 else 1)
            ax.plot([t, t], [r_out, r_out + 0.012], color="#555555",
                    lw=0.6, zorder=3)
            rot = _screen_angle_deg(t) - 90.0
            if 90 < (rot % 360) < 270:
                rot += 180
            label = "0" if bp == 0 else "%.1fk" % (bp / 1000.0)
            ax.text(t, r_out + 0.075, label, rotation=rot,
                    rotation_mode="anchor", ha="center", va="center",
                    fontsize=5.5, color="#555555", zorder=3)

    # --- Variant lollipops ----------------------------------------------
    legend_handles = []
    if variants is not None and len(variants) > 0:
        df = coerce_variant_frame(variants, vaf_col=vaf_col, pos_col=pos_col)
        color_map, group_col = _lollipop_color_map(df, color_by, colors)

        for _, row in df.iterrows():
            t = _theta(row["pos"])
            vaf = row["vaf"]
            frac = float(vaf) if np.isfinite(vaf) else 1.0
            frac = min(max(frac, 0.0), 1.0)
            r_tip = lollipop_base - frac * (lollipop_base - lollipop_min)
            key = row.get(group_col, None) if group_col else None
            col = color_map.get(key, "#444444")
            ax.plot([t, t], [lollipop_base, r_tip], color=col, lw=0.8,
                    zorder=4, alpha=0.9)
            ax.scatter([t], [r_tip], s=marker_size, color=col,
                       edgecolor="white", linewidth=0.4, zorder=6)

        # A dedicated marker legend is only needed when the lollipops are
        # coloured by something other than feature type (otherwise the gene
        # feature legend below already explains the colours).
        if legend and color_by != "feature":
            import matplotlib.lines as mlines
            for key, col in color_map.items():
                legend_handles.append(
                    mlines.Line2D([], [], color=col, marker="o", linestyle="",
                                  markersize=6, label=str(key)))

    # --- Feature-type legend (always meaningful) ------------------------
    if legend:
        import matplotlib.patches as mpatches
        feat_handles = [
            mpatches.Patch(color=colors[f], label=MT_FEATURE_LABELS[f])
            for f in ("protein_coding", "rRNA", "tRNA", "control_region")
        ]
        feat_legend = ax.legend(handles=feat_handles, loc="upper left",
                                bbox_to_anchor=(1.02, 1.0), frameon=False,
                                fontsize=7, title="Feature")
        if legend_handles:
            # Keep the feature legend and add a second one for the marker
            # encoding (status / variant type).
            ax.add_artist(feat_legend)
            ax.legend(handles=legend_handles, loc="lower left",
                      bbox_to_anchor=(1.02, 0.0), frameon=False,
                      fontsize=7, title=_legend_title(color_by))

    if title:
        ax.text(0, 0, title, ha="center", va="center", fontsize=10,
                fontweight="bold", transform=ax.transData)

    return ax


def _lollipop_color_map(df, color_by, feature_colors):
    """Return (``{key: color}``, group_column) for the requested encoding."""
    if color_by == "status" and "status" in df.columns:
        palette = {"HET": "#C44E52", "HOM": "#4C72B0"}
        keys = [k for k in ("HET", "HOM") if k in set(df["status"])]
        return {k: palette.get(k, "#444444") for k in keys}, "status"
    if color_by == "var_type" and "var_type" in df.columns:
        palette = {"SNV": "#4C72B0", "INS": "#55A868", "DEL": "#DD8452",
                   "MNV": "#8172B3", "REF": "#999999"}
        keys = list(dict.fromkeys(df["var_type"]))
        return {k: palette.get(k, "#444444") for k in keys}, "var_type"
    # Default: colour by feature type.
    return dict(feature_colors), "feature_type"


def _legend_title(color_by):
    return {"status": "Zygosity", "var_type": "Variant type",
            "feature": "Feature"}.get(color_by, "")
