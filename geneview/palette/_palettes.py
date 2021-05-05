from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable


def generate_colors_palette(cmap="viridis", n_colors=10, alpha=1.0):
    """Generate colors from matplotlib colormap; pass list to use exact colors"""
    if isinstance(cmap, list):
        colors = [list(to_rgba(color, alpha=alpha)) for color in cmap]
    else:
        scalar_mappable = ScalarMappable(cmap=cmap)
        colors = scalar_mappable.to_rgba(range(n_colors), alpha=alpha).tolist()
    return colors
