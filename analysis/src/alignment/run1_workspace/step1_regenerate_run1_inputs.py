#!/usr/bin/env python3
"""
Regenerate the TRUE original run1 inputs (atac_affine_aligned.h5ad / xenium_affine_aligned.h5ad)
using the exact original MOSAICField.ipynb (archived) logic: eyeballed scalefac=1.35 applied to
Xenium, manually-cropped Giotto Xenium export, ATAC rightmost-50% trim, free (unconstrained-scale)
affine_align. This reproduces cells 8-25 of analysis/src/alignment/archive/MOSAICField.ipynb
verbatim (read-only reference, not modified).

Cropping happens exactly ONCE, before affine_align -- no second crop/refit pass afterward. That
two-pass pattern belongs to a separate, newer rigid-transform effort and would not faithfully
reproduce this original methodology.

All outputs are written into run1_workspace/mosaicfield_outputs/ -- an isolated directory, never
touching the shared analysis/src/alignment/mosaicfield_outputs/ used by other work.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

MOSAICFIELD_PKG = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/MOSAICField"
sys.path.append(os.path.abspath(f"{MOSAICFIELD_PKG}/src/MOSAICField"))
sys.path.append(os.path.abspath(MOSAICFIELD_PKG))
from src.MOSAICField.affine_alignment import affine_align

outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mosaicfield_outputs")
os.makedirs(outdir, exist_ok=True)
print(f"Isolated output directory: {outdir}")


def rotate_coords(coords, angle_deg, center=None):
    theta = np.deg2rad(angle_deg)
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    if center is None:
        center = np.array([0, 0])
    return (coords - center) @ R.T + center


def plot_slices_overlap(coords1, coords2, title="Overlap"):
    plt.figure(figsize=(12, 6))
    plt.scatter(coords1[:, 0], coords1[:, 1], s=1, alpha=0.6, label="ATAC", c="red")
    plt.scatter(coords2[:, 0], coords2[:, 1], s=1, alpha=0.6, label="Xenium", c="blue")
    plt.gca().set_aspect("equal")
    plt.legend()
    plt.title(title)
    plt.savefig(os.path.join(outdir, f"run1_{title.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '')}.png"),
                dpi=100, bbox_inches="tight")
    plt.close()


# ============================================================================
# Cell 8: load raw ATAC (unchanged from archive/MOSAICField.ipynb)
# ============================================================================
print("Loading raw ATAC...")
atac = sc.read_h5ad(
    "/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/"
    "ultima1942optimize/set1_vf25000-cr1-0-vi1-nc30/combined.h5ad"
)

# ============================================================================
# Cells 12-13: Xenium h5ad already exists on disk (built previously by this
# exact original pipeline) -- reuse it as-is, do NOT rebuild (rebuilding would
# needlessly overwrite the shared xenium_cells_features_coords.h5ad used
# elsewhere in the project).
# ============================================================================
xenium_path = ("/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/"
               "multiomic/Xenium_488B/giotto_output/xenium_cells_features_coords.h5ad")
print(f"Loading existing Xenium h5ad: {xenium_path}")
xenium = sc.read_h5ad(xenium_path)
print(f"Raw Xenium: {xenium.shape}")

# ============================================================================
# Cell 18: ATAC trim to rightmost 50% (tissue coverage) -- ONE crop, before affine
# ============================================================================
atac_coords = np.asarray(atac.obsm["spatial"], dtype=float)
x_min, x_max = atac_coords[:, 0].min(), atac_coords[:, 0].max()
x_range = x_max - x_min
x_threshold = x_max - 0.5 * x_range
mask_rightmost = atac_coords[:, 0] > x_threshold
atac_trimmed = atac[mask_rightmost].copy()
print(f"Original ATAC cells: {atac.shape[0]}")
print(f"Trimmed ATAC cells (rightmost tissue): {atac_trimmed.shape[0]}")
atac_coords_trimmed = atac_coords[mask_rightmost]
atac_trimmed.obsm["spatial"] = atac_coords_trimmed

# ============================================================================
# Cell 19: rotate ATAC 270 degrees
# ============================================================================
atac_coords_trimmed = np.asarray(atac_trimmed.obsm["spatial"], dtype=float)
atac_coords_rotated = rotate_coords(atac_coords_trimmed, 270)
atac_trimmed.obsm["spatial"] = atac_coords_rotated
print("ATAC rotated by 270 degrees")

# ============================================================================
# Cell 20: mirror Xenium y-axis (Xenium crop already baked into the h5ad --
# this is the ONLY Xenium crop, applied once, before affine_align)
# ============================================================================
xenium_coords = xenium.obs[["sdimx", "sdimy"]].to_numpy()
xenium_coords_flipped = xenium_coords.copy()
xenium_coords_flipped[:, 1] = -xenium_coords_flipped[:, 1]
xenium.obsm["spatial"] = xenium_coords_flipped
print("Xenium mirrored (y-axis flipped)")

# ============================================================================
# Cell 22: eyeballed scalefac=1.35 applied to Xenium (the original approach,
# reproduced verbatim -- NOT the corrected calibration used elsewhere), mirror
# ATAC y-axis, center both
# ============================================================================
xenium_coords = xenium.obs[["sdimx", "sdimy"]].to_numpy()
atac_coords = atac_trimmed.obsm["spatial"].copy()
atac_coords[:, 1] *= -1

scalefac = 1.35
xenium_coords_scaled = xenium_coords * scalefac

xenium_aligned = xenium.copy()
atac_aligned = atac_trimmed.copy()
xenium_aligned.obsm["spatial"] = xenium_coords_scaled
atac_aligned.obsm["spatial"] = atac_coords

atac_center = atac_aligned.obsm["spatial"].mean(axis=0)
xenium_center = xenium_aligned.obsm["spatial"].mean(axis=0)
atac_aligned.obsm["spatial"] -= atac_center
xenium_aligned.obsm["spatial"] -= xenium_center

plot_slices_overlap(atac_aligned.obsm["spatial"], xenium_aligned.obsm["spatial"],
                     "Before Affine Alignment mirrored scaled centered")

# ============================================================================
# Cell 23: PCA on both (normalize_total + log1p + HVG + PCA)
# ============================================================================
sc.pp.normalize_total(atac_aligned, inplace=True)
sc.pp.log1p(atac_aligned)
sc.pp.highly_variable_genes(atac_aligned, flavor="seurat", n_top_genes=2000, inplace=True, subset=True)
sc.pp.pca(atac_aligned, n_comps=50)

sc.pp.normalize_total(xenium_aligned, inplace=True)
sc.pp.log1p(xenium_aligned)
sc.pp.highly_variable_genes(xenium_aligned, flavor="seurat", n_top_genes=2000, inplace=True, subset=True)
sc.pp.pca(xenium_aligned, n_comps=50)
print("PCA completed for both ATAC and Xenium datasets")

# ============================================================================
# Cell 24: free (unconstrained-scale) affine_align -- the ORIGINAL alignment
# method, single pass, no crop/refit afterward
# ============================================================================
np.random.seed(0)
n_subsample = 2000
atac_n = min(n_subsample, atac_aligned.shape[0])
xenium_n = min(n_subsample, xenium_aligned.shape[0])

atac_indices = np.random.choice(atac_aligned.shape[0], atac_n, replace=False)
xenium_indices = np.random.choice(xenium_aligned.shape[0], xenium_n, replace=False)

atac_sub = atac_aligned[atac_indices].copy()
xenium_sub = xenium_aligned[xenium_indices].copy()

print(f"  ATAC: {atac_aligned.shape[0]} -> {atac_sub.shape[0]} cells")
print(f"  Xenium: {xenium_aligned.shape[0]} -> {xenium_sub.shape[0]} cells")

print("Computing affine transformation (free/unconstrained scale, original method)...")
atac_aligned_sub, xenium_aligned_sub, T, P = affine_align(
    atac_sub, xenium_sub, obsm_name="X_pca", max_iter=10, alpha=0.9
)


def affine_transform(coords, T):
    homogeneous = np.vstack([coords.T, np.ones((1, coords.shape[0]))])
    return (T @ homogeneous)[:2, :].T


atac_coords_full = np.asarray(atac_aligned.obsm["spatial"], dtype=float)
atac_coords_transformed = affine_transform(atac_coords_full, T)

atac_affine_aligned = atac_aligned.copy()
atac_affine_aligned.obsm["spatial"] = atac_coords_transformed
print("Affine alignment complete (computed on 2000-pt subsample, applied to full ATAC dataset)")

plot_slices_overlap(atac_affine_aligned.obsm["spatial"], xenium_aligned.obsm["spatial"],
                     "After Affine Alignment Manual 2000pt Subsample")

# Sanity check (informational only -- expect the ORIGINAL, uncorrected-scale
# numbers here, ~5-9 cells/spot, NOT the rigid-transform fix's ~1-3)
from scipy.spatial import cKDTree
import pandas as pd
tree = cKDTree(atac_affine_aligned.obsm["spatial"])
_, nn_idx = tree.query(xenium_aligned.obsm["spatial"], k=1)
cells_per_spot = pd.Series(nn_idx).value_counts()
print(f"\n[Sanity check] median cells/spot = {cells_per_spot.median():.2f}, "
      f"mean = {cells_per_spot.mean():.2f} "
      f"(expect original ~5-9, NOT the rigid-fix's ~1-3, to confirm this is genuinely run1)")

# ============================================================================
# Cell 25: save
# ============================================================================
atac_affine_aligned.write_h5ad(os.path.join(outdir, "atac_affine_aligned.h5ad"))
xenium_aligned.write_h5ad(os.path.join(outdir, "xenium_affine_aligned.h5ad"))
print(f"\nSaved run1's original atac_affine_aligned.h5ad / xenium_affine_aligned.h5ad to {outdir}/")
