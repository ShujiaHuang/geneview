from matplotlib.colors import to_rgba
from matplotlib.cm import ScalarMappable


def generate_colors_palette(cmap="viridis", n_colors=10):
    """Generate colors from matplotlib colormap; pass list to use exact colors"""
    if isinstance(cmap, list):
        colors = [list(to_rgba(color)) for color in cmap]
    else:
        scalar_mappable = ScalarMappable(cmap=cmap)
        colors = scalar_mappable.to_rgba(range(n_colors)).tolist()
    return colors
