"""
Color scheme presets for genome tracks.

Provides named color schemes that can be applied to GeneRegionTrack and
AnnotationTrack objects, similar to Gviz's ``setScheme()`` / display
parameter presets.
"""

from typing import Dict, Any, Optional

# Categorical palettes for gene/transcript coloring
_GENE_PALETTE = [
    "#3C5488", "#E64B35", "#00A087", "#F39B7F",
    "#8491B4", "#91D1C2", "#B09C85", "#7E6148",
    "#DC0000", "#0080FF", "#00BA38", "#FDB462",
]

_TRANSCRIPT_PALETTE = [
    "#F0E442", "#56B4E9", "#009E73", "#CC79A7",
    "#E69F00", "#D55E00", "#0072B2", "#999999",
]


_SCHEMES: Dict[str, Dict[str, Any]] = {
    "default": {
        "fill": "lightgray",
        "col": "darkgray",
    },
    "genes": {
        "palette": _GENE_PALETTE,
        "mode": "gene",
    },
    "transcripts": {
        "palette": _TRANSCRIPT_PALETTE,
        "mode": "transcript",
    },
}


def get_scheme_names():
    """Return available scheme names."""
    return list(_SCHEMES.keys())


def apply_scheme(track, scheme_name: str) -> None:
    """Apply a named color scheme to a track.

    Parameters
    ----------
    track : Track
        A GeneRegionTrack or AnnotationTrack instance.
    scheme_name : str
        One of ``"default"``, ``"genes"``, ``"transcripts"``.

    Raises
    ------
    ValueError
        If the scheme name is not recognized.
    """
    if scheme_name not in _SCHEMES:
        raise ValueError(
            f"Unknown scheme '{scheme_name}'. "
            f"Available: {get_scheme_names()}"
        )

    scheme = _SCHEMES[scheme_name]

    if scheme_name == "default":
        track.set_param("fill", scheme["fill"])
        track.set_param("col", scheme["col"])
        return

    # For genes/transcripts schemes, assign colors per group
    palette = scheme["palette"]
    mode = scheme["mode"]

    data = getattr(track, "data", None)
    if data is None:
        return

    # Determine grouping column
    group_col = None
    if mode == "gene":
        for col in ["gene_name", "gene_id", "gene"]:
            if col in data.columns:
                group_col = col
                break
    elif mode == "transcript":
        for col in ["transcript_id", "transcript"]:
            if col in data.columns:
                group_col = col
                break

    if group_col is None:
        # No grouping column, just apply the first palette color
        track.set_param("fill", palette[0])
        return

    # Build a color map for each unique group
    unique_groups = data[group_col].dropna().unique()
    color_map = {}
    for i, grp in enumerate(unique_groups):
        color_map[grp] = palette[i % len(palette)]

    # Store the color map on the track for use during drawing
    track._scheme_color_map = color_map
    track._scheme_group_col = group_col
