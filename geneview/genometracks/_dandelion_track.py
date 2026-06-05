"""
DandelionTrack — clustered variant track for genometracks.

Ported from trackViewer's ``dandelion.plot()`` (R/Bioconductor).  Nearby
variants are clustered into groups; each group is drawn as a "dandelion"
with stems fanning out from a central base point.  Works both as a
composable track inside :func:`plot_tracks` and as a standalone plot via
:func:`dandelion_plot`.
"""

from typing import Optional, Dict, Any, List, Callable, Tuple

import numpy as np
import pandas as pd

from ._base import Track, GenomicInterval
from ._mutation_features import draw_features, rescale_position
from ._mutation_shapes import (
    draw_circle, draw_pie, draw_fan, draw_pin, draw_stem, draw_stem_bent,
    _blended, _data_radius,
)
from ._lollipop_track import _validate_snp_data, _validate_features, _snp_val


# ---------------------------------------------------------------------------
# DandelionTrack
# ---------------------------------------------------------------------------

class DandelionTrack(Track):
    """Dandelion-style clustered variant track.

    Nearby variants are grouped into clusters.  Each cluster is drawn as a
    dandelion: a central stem rises from the feature baseline to a height
    proportional to the cluster size (or a custom ``height_method``), and
    individual variants fan out from the stem top.

    Can be used inside :func:`plot_tracks` alongside other tracks, or
    standalone via the :func:`dandelion_plot` convenience function.

    Parameters
    ----------
    snp_data : pd.DataFrame
        Variant data with columns: ``chrom``, ``start``.
        Optional: ``score``, ``color``, ``fill``, ``border``, ``label``.
    features : pd.DataFrame or None
        Gene-feature data (same format as :class:`LolliplotTrack`).
    type : str
        Shape type: ``'fan'``, ``'circle'``, ``'pie'``, ``'pin'``.
    maxgaps : float
        Maximum gap between clustered variants as a fraction of the total
        region width.  Default ``1/50`` (2%).
    height_method : callable or None
        ``f(scores) -> float`` used to compute stem height from a cluster's
        score array.  Default is ``len`` (count of variants).
    cex : float
        Shape-size scaling factor.
    show_yaxis : bool
        Show a y-axis.
    label_on_feature : bool
        Place feature labels on rectangles.
    rescale : list or None
        Coordinate remapping: list of ``(from_start, from_end, to_start,
        to_end)`` tuples.
    name : str
        Track name.
    height : float
        Relative track height.
    """

    TYPES = ("fan", "circle", "pie", "pin")

    def __init__(
        self,
        snp_data: pd.DataFrame,
        features: Optional[pd.DataFrame] = None,
        type: str = "fan",
        maxgaps: float = 1 / 50,
        height_method: Optional[Callable] = None,
        cex: float = 1.0,
        show_yaxis: bool = False,
        label_on_feature: bool = False,
        rescale: Optional[list] = None,
        name: str = "Dandelion",
        height: float = 3.0,
        display_params: Optional[Dict[str, Any]] = None,
    ):
        if type not in self.TYPES:
            raise ValueError(f"type must be one of {self.TYPES}, got '{type}'")
        super().__init__(name=name, height=height, display_params=display_params)
        self.snp_data = _validate_snp_data(snp_data)
        self.features = _validate_features(features)
        self.type = type
        self.maxgaps = maxgaps
        self.height_method = height_method or len
        self.cex = cex
        self.show_yaxis = show_yaxis
        self.label_on_feature = label_on_feature
        self.rescale = rescale

    # ------------------------------------------------------------------
    # Clustering
    # ------------------------------------------------------------------

    def _cluster_snps(self, positions: np.ndarray, region_width: int):
        """Group SNP positions into clusters based on ``maxgaps``."""
        if len(positions) == 0:
            return []
        sorted_idx = np.argsort(positions)
        sorted_pos = positions[sorted_idx]
        max_gap = region_width * self.maxgaps

        clusters = []
        current_cluster = [sorted_idx[0]]
        for i in range(1, len(sorted_pos)):
            if sorted_pos[i] - sorted_pos[i - 1] > max_gap:
                clusters.append(np.array(current_cluster))
                current_cluster = [sorted_idx[i]]
            else:
                current_cluster.append(sorted_idx[i])
        clusters.append(np.array(current_cluster))
        return clusters

    # ------------------------------------------------------------------
    # Track interface
    # ------------------------------------------------------------------

    def get_region(self) -> Optional[GenomicInterval]:
        """Return the full genomic extent of the data + features."""
        all_starts = []
        all_ends = []
        if len(self.snp_data) > 0:
            all_starts.append(int(self.snp_data["start"].min()))
            all_ends.append(int(self.snp_data["end"].max()))
        if self.features is not None and len(self.features) > 0:
            all_starts.append(int(self.features["start"].min()))
            all_ends.append(int(self.features["end"].max()))
        if not all_starts:
            return None
        chrom = str(self.snp_data["chrom"].iloc[0])
        return GenomicInterval(
            chrom=chrom,
            start=min(all_starts),
            end=max(all_ends),
        )

    def draw(self, ax, region: GenomicInterval) -> None:
        """Render the dandelion track onto the given Axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        region : GenomicInterval
        """
        r_start = region.start
        r_end = region.end
        region_width = r_end - r_start

        ax.set_ylim(0, 1)
        ax.get_yaxis().set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        # ---- Draw features ----
        baseline_y = 0.08
        feature_height = draw_features(
            ax, self.features, r_start, r_end,
            baseline_y=baseline_y,
            label_on_feature=self.label_on_feature,
            rescale_map=self.rescale,
        )

        # ---- Filter SNPs to region ----
        snps = self.snp_data[
            (self.snp_data["start"] >= r_start)
            & (self.snp_data["start"] <= r_end)
        ].copy()

        if len(snps) == 0:
            return

        snps = snps.sort_values("start").reset_index(drop=True)

        # ---- Apply rescale ----
        if self.rescale is not None:
            snps["_draw_x"] = snps["start"].apply(
                lambda p: rescale_position(float(p), self.rescale))
        else:
            snps["_draw_x"] = snps["start"].astype(float)

        # ---- Cluster SNPs ----
        positions = snps["start"].values
        clusters = self._cluster_snps(positions, region_width)

        # ---- Compute layout parameters ----
        feature_top = baseline_y + feature_height
        stem_base_y = feature_top + 0.005
        radius_axes = 0.04 * self.cex

        cluster_heights = []
        for cluster_idx in clusters:
            scores = snps.iloc[cluster_idx]["score"].values
            cluster_heights.append(self.height_method(scores))

        max_height = max(cluster_heights) if cluster_heights else 1
        max_stem_top = 0.80
        label_space = 0.12
        available_height = max_stem_top - stem_base_y - label_space
        yaxis_max = max_height

        for ci, cluster_idx in enumerate(clusters):
            cluster_snps = snps.iloc[cluster_idx]
            n = len(cluster_snps)
            raw_h = cluster_heights[ci]

            if max_height > 0:
                norm_h = raw_h / max_height
            else:
                norm_h = 0.0
            stem_top_y = stem_base_y + norm_h * available_height

            center_x = float(cluster_snps["_draw_x"].mean())

            # Draw central stem — from feature top for visual connection
            draw_stem(
                ax, center_x, feature_top, stem_top_y,
                color="black", linewidth=1.0, linestyle="-",
            )

            # Fan out individual stems
            if n == 1:
                snp = cluster_snps.iloc[0]
                snp_cex = float(_snp_val(snp, "cex", self.cex))
                snp_r = radius_axes * (snp_cex / self.cex) if self.cex else radius_axes
                self._draw_shape(ax, center_x, stem_top_y, snp, snp_r)
                label = _snp_val(snp, "label", _snp_val(snp, "name", ""))
                if label and str(label).strip():
                    label_rot = float(_snp_val(snp, "label_rotation", 90))
                    label_col = _snp_val(snp, "label_color", "black")
                    _, ry = _data_radius(ax, snp_r)
                    ax.text(
                        center_x, stem_top_y + ry + 0.005,
                        str(label),
                        ha="center", va="bottom", rotation=label_rot,
                        fontsize=6 * snp_cex, color=label_col,
                        zorder=6, clip_on=False,
                    )
            else:
                fan_r = radius_axes * 2.5 * min(n / 3.0, 3.0)
                _, ry_fan = _data_radius(ax, fan_r)

                for si in range(n):
                    snp = cluster_snps.iloc[si]
                    snp_cex = float(_snp_val(snp, "cex", self.cex))
                    snp_r = radius_axes * 0.8 * (snp_cex / self.cex) if self.cex else radius_axes * 0.8
                    frac = si / (n - 1) if n > 1 else 0.5
                    angle = (frac - 0.5) * np.pi * 0.75

                    rx_fan, ry_fan = _data_radius(ax, fan_r)
                    x2 = center_x + rx_fan * np.sin(angle)
                    y2 = stem_top_y + ry_fan * np.cos(angle)

                    draw_stem_bent(
                        ax,
                        x0=center_x, y0=stem_top_y,
                        x1=center_x, y1=stem_top_y,
                        x2=x2, y2=y2,
                        color="black", linewidth=0.8,
                    )

                    self._draw_shape(ax, x2, y2, snp, snp_r)

                    label = _snp_val(snp, "label", _snp_val(snp, "name", ""))
                    if label and str(label).strip():
                        label_rot = float(_snp_val(snp, "label_rotation", 90))
                        label_col = _snp_val(snp, "label_color", "black")
                        _, ry_s = _data_radius(ax, snp_r)
                        ax.text(
                            x2, y2 + ry_s + 0.005,
                            str(label),
                            ha="center", va="bottom", rotation=label_rot,
                            fontsize=5 * snp_cex, color=label_col,
                            zorder=6, clip_on=False,
                        )

            # Dashed guide line from stem top to top of plot area
            dash_col = "#CCCCCC"
            # Use first SNP's dashline_col if available
            first_snp = cluster_snps.iloc[0]
            dash_col = _snp_val(first_snp, "dashline_col", "#CCCCCC")
            _, ry_guide = _data_radius(ax, radius_axes)
            ax.plot(
                [center_x, center_x],
                [stem_top_y + ry_guide, max_stem_top + 0.05],
                color=dash_col, linewidth=0.5, linestyle=":",
                zorder=1, clip_on=False,
            )

        # ---- Y-axis ----
        if self.show_yaxis and max_height > 1:
            ax.text(
                -0.02, stem_base_y, "0",
                transform=ax.transAxes,
                ha="right", va="center", fontsize=6, color="gray",
            )
            ax.text(
                -0.02, max_stem_top, str(int(yaxis_max)),
                transform=ax.transAxes,
                ha="right", va="center", fontsize=6, color="gray",
            )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _draw_shape(self, ax, x, y, snp, radius_axes):
        """Draw a shape at (x, y) based on ``self.type``."""
        fill = snp.get("fill", "#0080FF") if "fill" in snp.index else "#0080FF"
        border = snp.get("border", "black") if "border" in snp.index else "black"
        score = float(snp.get("score", 1)) if "score" in snp.index else 1.0

        if pd.isna(fill):
            fill = "#0080FF"
        if pd.isna(border):
            border = "black"

        if self.type == "fan":
            fan_score = min(score, 1.0) if score <= 1.0 else score / max(score, 1.0)
            draw_fan(ax, x, y, radius_axes, fan_score, fill, border)
        elif self.type == "circle":
            draw_circle(ax, x, y, radius_axes, fill, border)
        elif self.type == "pie":
            pie_vals = snp.get("pie_values", [1, 1]) if "pie_values" in snp.index else [1, 1]
            pie_cols = snp.get("pie_colors", ["#0080FF", "#E69F00"]) if "pie_colors" in snp.index else ["#0080FF", "#E69F00"]
            if isinstance(pie_vals, str):
                pie_vals = [float(v) for v in pie_vals.split(",")]
            if isinstance(pie_cols, str):
                pie_cols = pie_cols.split(",")
            draw_pie(ax, x, y, radius_axes, pie_vals, pie_cols)
        elif self.type == "pin":
            draw_pin(
                ax, x, y, radius_axes,
                facecolor=fill if fill else "#D55E00",
                edgecolor=border,
            )

    def __repr__(self):
        return (
            f"DandelionTrack(name='{self.name}', height={self.height}, "
            f"type='{self.type}', n_variants={len(self.snp_data)})"
        )


# ---------------------------------------------------------------------------
# Standalone convenience function
# ---------------------------------------------------------------------------

def dandelion_plot(
    snp_data: pd.DataFrame,
    features: Optional[pd.DataFrame] = None,
    type: str = "fan",
    region: Optional[GenomicInterval] = None,
    figsize: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    show_title: bool = False,
    ax=None,
    **kwargs,
):
    """Convenience function: create a :class:`DandelionTrack` and render it.

    When called without ``ax``, creates a new figure via
    :func:`plot_tracks`.  When ``ax`` is provided, draws directly
    onto it (standalone mode).

    Parameters
    ----------
    snp_data : pd.DataFrame
        Variant data (see :class:`DandelionTrack`).
    features : pd.DataFrame or None
        Gene-feature data.
    type : str
        Shape type: ``'fan'``, ``'circle'``, ``'pie'``, ``'pin'``.
    region : GenomicInterval, optional
        Genomic region.  Auto-determined from data if *None*.
    figsize : tuple, optional
        Figure size in inches.  Default ``(12, 4)``.
    title : str, optional
        Plot title.
    show_title : bool
        Whether to show the track title panel.  Default ``False``.
    ax : matplotlib Axes, optional
        Draw into this axes instead of creating a new figure.
    **kwargs
        Additional keyword arguments passed to :class:`DandelionTrack`.

    Returns
    -------
    matplotlib.axes.Axes or list of Axes
    """
    track = DandelionTrack(
        snp_data, features=features, type=type, **kwargs,
    )

    if ax is not None:
        if region is None:
            region = track.get_region()
            if region is None:
                raise ValueError("Cannot determine region from data.")
            span = region.end - region.start
            region = GenomicInterval(
                region.chrom,
                max(0, region.start - int(span * 0.02)),
                region.end + int(span * 0.02),
            )
        ax.set_xlim(region.start, region.end)
        track.draw(ax, region)
        return ax

    from ._track_plot import plot_tracks

    if region is None:
        region = track.get_region()
        if region is not None:
            span = region.end - region.start
            region = GenomicInterval(
                region.chrom,
                max(0, region.start - int(span * 0.02)),
                region.end + int(span * 0.02),
            )

    axes = plot_tracks(
        [track],
        region=region,
        figsize=figsize or (12, 4),
        title=title,
        show_title=show_title,
    )
    return axes[0] if isinstance(axes, list) else axes
