"""
AtlasXomics-style spatial ATAC grid plotting.

Draws each ATAC spot as a non-overlapping square tile at its true physical pitch,
rather than a scatter dot of arbitrary size — tiles abut edge-to-edge with no gaps
or overlaps only once coordinates are true microns (see mosaic_run2-2.ipynb's
corrected-calibration section). Useful both as a visualization and as a visual QC
check of that calibration: gaps/overlaps between tiles indicate the coordinates are
not in true microns.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


def plot_atac_grid(ax, coords_um, values=None, pitch_um=10.0, cmap="viridis",
                    color="lightgrey", edgecolor="none", vmin=None, vmax=None,
                    alpha=1.0):
    """
    Draw each ATAC spot as a non-overlapping pitch_um x pitch_um square tile,
    centered at coords_um, instead of a scatter dot.

    coords_um : (N, 2) array, true microns (post calibration fix).
    values    : optional (N,) array to color-map tiles (else uniform `color`).

    Returns the PatchCollection (for colorbar attachment).
    """
    half = pitch_um / 2.0
    rects = [Rectangle((x - half, y - half), pitch_um, pitch_um) for x, y in coords_um]
    pc = PatchCollection(rects, edgecolor=edgecolor, linewidth=0.0)
    if values is not None:
        pc.set_array(np.asarray(values))
        pc.set_cmap(cmap)
        if vmin is not None or vmax is not None:
            pc.set_clim(vmin, vmax)
    else:
        pc.set_facecolor(color)
    pc.set_alpha(alpha)
    ax.add_collection(pc)
    ax.set_xlim(coords_um[:, 0].min() - pitch_um, coords_um[:, 0].max() + pitch_um)
    ax.set_ylim(coords_um[:, 1].min() - pitch_um, coords_um[:, 1].max() + pitch_um)
    ax.set_aspect("equal", adjustable="box")
    return pc


def plot_atac_grid_with_xenium(ax, atac_coords_um, xenium_coords_um, atac_values=None,
                                xenium_point_size=1.5, xenium_color="black",
                                xenium_alpha=0.4, pitch_um=10.0, highlight_mask=None,
                                highlight_color="red", **grid_kwargs):
    """
    Overlay: ATAC squares (plot_atac_grid) + Xenium cells as small scatter points
    (Xenium cells are sub-pitch scale, so points rather than tiles).

    highlight_mask : optional bool array, same length as atac_coords_um — recolors
    a subset of tiles, e.g. to compare NN- ("Method A") vs soft-threshold-
    ("Method B"/"Voronoi") assigned spots.
    """
    pc = plot_atac_grid(ax, atac_coords_um, values=atac_values, pitch_um=pitch_um, **grid_kwargs)
    if highlight_mask is not None:
        plot_atac_grid(ax, atac_coords_um[highlight_mask], color=highlight_color,
                        pitch_um=pitch_um, alpha=0.6)
    ax.scatter(xenium_coords_um[:, 0], xenium_coords_um[:, 1], s=xenium_point_size,
               c=xenium_color, alpha=xenium_alpha, linewidths=0, zorder=5)
    return pc


def zoom_to_region(ax, center_um, half_width_um):
    """Set axis limits for a zoomed-in view (e.g. half_width_um=150 for ~15x15 spots)."""
    cx, cy = center_um
    ax.set_xlim(cx - half_width_um, cx + half_width_um)
    ax.set_ylim(cy - half_width_um, cy + half_width_um)
    ax.set_aspect("equal", adjustable="box")
