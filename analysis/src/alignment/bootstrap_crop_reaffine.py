#!/usr/bin/env python3
"""
Bootstrap crop + re-affine step.

Problem: Xenium imaged a larger tissue region than the ATAC chip covers. The
first-pass affine_align (already run, see mosaicfield_outputs/atac_affine_aligned.h5ad)
had no way to know this and stretched ATAC ~2.2x to span roughly the same
footprint as the full (uncropped) Xenium point cloud, inflating spot spacing
from the true 10um pitch to ~22um.

Fix: use the first-pass affine transform's own output (atac_affine_aligned,
already expressed in Xenium's coordinate frame) to find where in Xenium's
frame the true ATAC-corresponding tissue sits. Crop Xenium to that region
(with margin), then rerun affine_align FROM ATAC'S TRUE CALIBRATED (pre-affine)
COORDINATES against only the cropped Xenium subset -- this time both point
clouds should have comparable physical extents, so there's no need to stretch
ATAC to fill excess space.

This script recomputes ATAC's true pre-affine coordinates fresh (cheap: just
load + trim + rotate + calibrate, no PCA/affine yet needed for that part) so
pass 2 fits a clean transform directly from true microns, not compounded on
top of pass 1's distortion.

Stops BEFORE rasterization/nonlinear alignment (per instruction) -- only
reports the resulting spot-spacing sanity check and plots, for review before
committing to the expensive downstream steps.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

sys.path.append(os.path.abspath("./MOSAICField/src/MOSAICField"))
sys.path.append(os.path.abspath("./MOSAICField"))
from src.MOSAICField.affine_alignment import affine_align

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"
outdir = "./mosaicfield_outputs"

# ============================================================================
# Step A: recompute ATAC's true pre-affine (calibrated) coordinates fresh
# (mirrors mosaic_run2-2.ipynb Step -1a/-1b/-1c exactly)
# ============================================================================
print("Recomputing ATAC's true pre-affine calibrated coordinates...")
atac = sc.read_h5ad(
    "/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/"
    "ultima1942optimize/set1_vf25000-cr1-0-vi1-nc30/combined.h5ad"
)
atac_coords_px = np.asarray(atac.obsm["spatial"], dtype=float)
x_min, x_max = atac_coords_px[:, 0].min(), atac_coords_px[:, 0].max()
x_threshold = x_max - 0.5 * (x_max - x_min)
mask_rightmost = atac_coords_px[:, 0] > x_threshold
atac_trimmed = atac[mask_rightmost].copy()
atac_trimmed.obsm["spatial"] = atac_coords_px[mask_rightmost]

def rotate_coords(coords, angle_deg, center=None):
    theta = np.deg2rad(angle_deg)
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    if center is None:
        center = np.array([0, 0])
    return (coords - center) @ R.T + center

atac_trimmed.obsm["spatial"] = rotate_coords(
    np.asarray(atac_trimmed.obsm["spatial"], dtype=float), 270
)

tissue_positions = pd.read_csv(
    os.path.join(PROJECT_ROOT, "Data/01_inputs/spatial/tissue_positions_list.csv")
)
import statsmodels.api as sm
X = sm.add_constant(tissue_positions[["array_row", "array_col"]])
fit_x = sm.OLS(tissue_positions["x_spatial"], X).fit()
fit_y = sm.OLS(tissue_positions["y_spatial"], X).fit()
PIXELS_PER_ARRAY_STEP = float(np.mean([abs(fit_x.params["array_row"]), abs(fit_y.params["array_col"])]))
UM_PER_PIXEL_ATAC = 10.0 / PIXELS_PER_ARRAY_STEP

atac_coords_true = atac_trimmed.obsm["spatial"].copy().astype(float) * UM_PER_PIXEL_ATAC
atac_coords_true[:, 1] *= -1

atac_aligned = atac_trimmed.copy()
atac_aligned.obsm["spatial"] = atac_coords_true
atac_center = atac_aligned.obsm["spatial"].mean(axis=0)
atac_aligned.obsm["spatial"] -= atac_center

print(f"ATAC true calibrated footprint: {atac_aligned.shape[0]:,} spots, "
      f"x=[{atac_aligned.obsm['spatial'][:,0].min():.0f}, {atac_aligned.obsm['spatial'][:,0].max():.0f}], "
      f"y=[{atac_aligned.obsm['spatial'][:,1].min():.0f}, {atac_aligned.obsm['spatial'][:,1].max():.0f}]")

sc.pp.normalize_total(atac_aligned, inplace=True)
sc.pp.log1p(atac_aligned)
sc.pp.highly_variable_genes(atac_aligned, flavor="seurat", n_top_genes=2000, inplace=True, subset=True)
sc.pp.pca(atac_aligned, n_comps=50)

# ============================================================================
# Step B: use pass-1 output to find the ATAC-corresponding Xenium region, crop
# ============================================================================
print("\nLoading first-pass affine-aligned objects (already on disk)...")
atac_pass1 = ad.read_h5ad(os.path.join(outdir, "atac_affine_aligned.h5ad"))
xenium_pass1 = ad.read_h5ad(os.path.join(outdir, "xenium_affine_aligned.h5ad"))

atac_coords_pass1 = atac_pass1.obsm["spatial"]  # already in Xenium's frame (post first-pass T)
xenium_coords = xenium_pass1.obsm["spatial"]    # Xenium's own frame, unchanged

print(f"ATAC (pass 1, in Xenium frame): x=[{atac_coords_pass1[:,0].min():.0f}, {atac_coords_pass1[:,0].max():.0f}], "
      f"y=[{atac_coords_pass1[:,1].min():.0f}, {atac_coords_pass1[:,1].max():.0f}]")
print(f"Xenium (full, uncropped): {xenium_coords.shape[0]:,} cells, "
      f"x=[{xenium_coords[:,0].min():.0f}, {xenium_coords[:,0].max():.0f}], "
      f"y=[{xenium_coords[:,1].min():.0f}, {xenium_coords[:,1].max():.0f}]")

MARGIN_FRAC = 0.25
x_min, x_max = atac_coords_pass1[:, 0].min(), atac_coords_pass1[:, 0].max()
y_min, y_max = atac_coords_pass1[:, 1].min(), atac_coords_pass1[:, 1].max()
x_pad = (x_max - x_min) * MARGIN_FRAC
y_pad = (y_max - y_min) * MARGIN_FRAC
crop_x = (x_min - x_pad, x_max + x_pad)
crop_y = (y_min - y_pad, y_max + y_pad)

crop_mask = (
    (xenium_coords[:, 0] >= crop_x[0]) & (xenium_coords[:, 0] <= crop_x[1]) &
    (xenium_coords[:, 1] >= crop_y[0]) & (xenium_coords[:, 1] <= crop_y[1])
)
print(f"\nCrop region (ATAC pass-1 bbox + {int(MARGIN_FRAC*100)}% margin): x={crop_x}, y={crop_y}")
print(f"Xenium cells kept after crop: {crop_mask.sum():,} / {len(crop_mask):,} ({100*crop_mask.mean():.1f}%)")

xenium_cropped = xenium_pass1[crop_mask].copy()

fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(xenium_coords[~crop_mask, 0], xenium_coords[~crop_mask, 1],
           s=0.5, alpha=0.15, c="lightgrey", label="Xenium (excluded, excess tissue)")
ax.scatter(xenium_coords[crop_mask, 0], xenium_coords[crop_mask, 1],
           s=0.5, alpha=0.3, c="royalblue", label="Xenium (kept, cropped region)")
ax.scatter(atac_coords_pass1[:, 0], atac_coords_pass1[:, 1],
           s=1, alpha=0.6, c="red", label="ATAC (pass-1 affine footprint)")
ax.set_aspect("equal")
ax.legend(markerscale=10)
ax.set_title("Bootstrap crop: Xenium region kept for pass-2 affine_align")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "bootstrap_crop_region.png"), dpi=150)
print("Saved bootstrap_crop_region.png")

# ============================================================================
# Step C: pass-2 affine_align -- true calibrated ATAC vs cropped Xenium
# ============================================================================
print("\nRecomputing PCA on cropped Xenium subset...")
sc.pp.pca(xenium_cropped, n_comps=50)

np.random.seed(1)
n_subsample = 2000
atac_n = min(n_subsample, atac_aligned.shape[0])
xenium_n = min(n_subsample, xenium_cropped.shape[0])

atac_indices = np.random.choice(atac_aligned.shape[0], atac_n, replace=False)
xenium_indices = np.random.choice(xenium_cropped.shape[0], xenium_n, replace=False)

atac_sub = atac_aligned[atac_indices].copy()
xenium_sub = xenium_cropped[xenium_indices].copy()

print(f"\nPass 2: affine_align, true-calibrated ATAC ({atac_aligned.shape[0]:,} spots) "
      f"vs cropped Xenium ({xenium_cropped.shape[0]:,} cells)...")
atac_aligned_sub2, xenium_aligned_sub2, T2, P2 = affine_align(
    atac_sub, xenium_sub, obsm_name="X_pca", max_iter=10, alpha=0.9
)

def affine_transform(coords, T):
    homogeneous = np.vstack([coords.T, np.ones((1, coords.shape[0]))])
    return (T @ homogeneous)[:2, :].T

atac_coords_pass2_input = np.asarray(atac_aligned.obsm["spatial"], dtype=float)
atac_coords_pass2 = affine_transform(atac_coords_pass2_input, T2)

print(f"Pass 2 ATAC footprint: x=[{atac_coords_pass2[:,0].min():.0f}, {atac_coords_pass2[:,0].max():.0f}], "
      f"y=[{atac_coords_pass2[:,1].min():.0f}, {atac_coords_pass2[:,1].max():.0f}]")

# ============================================================================
# Step D: sanity check -- inter-spot spacing + cells/spot after pass 2
# ============================================================================
tree2 = cKDTree(atac_coords_pass2)
inter_spot_dist2, _ = tree2.query(atac_coords_pass2, k=2)
median_inter_spot2 = np.median(inter_spot_dist2[:, 1])
print(f"\nPass 2: median inter-spot distance = {median_inter_spot2:.2f} um "
      f"(true pitch = 10 um; pass-1 gave 22.1 um)")

xenium_coords_cropped = xenium_cropped.obsm["spatial"]
nn_dist2, nn_idx2 = tree2.query(xenium_coords_cropped, k=1)
cells_per_spot2 = pd.Series(nn_idx2).value_counts()
print(f"Pass 2: median cells/spot = {cells_per_spot2.median():.2f}, "
      f"mean = {cells_per_spot2.mean():.2f} "
      f"(true expectation ~1-3; pass-1 full-pipeline result was median=5, mean=5.2)")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].scatter(xenium_coords[crop_mask, 0], xenium_coords[crop_mask, 1],
                s=0.5, alpha=0.2, c="royalblue", label="Xenium (cropped)")
axes[0].scatter(atac_coords_pass1[:, 0], atac_coords_pass1[:, 1],
                s=1, alpha=0.6, c="red", label="ATAC (pass 1)")
axes[0].set_aspect("equal"); axes[0].legend(markerscale=10)
axes[0].set_title("Pass 1 (uncropped Xenium target)\nmedian inter-spot = 22.1 um")

axes[1].scatter(xenium_coords[crop_mask, 0], xenium_coords[crop_mask, 1],
                s=0.5, alpha=0.2, c="royalblue", label="Xenium (cropped)")
axes[1].scatter(atac_coords_pass2[:, 0], atac_coords_pass2[:, 1],
                s=1, alpha=0.6, c="red", label="ATAC (pass 2)")
axes[1].set_aspect("equal"); axes[1].legend(markerscale=10)
axes[1].set_title(f"Pass 2 (cropped Xenium target)\nmedian inter-spot = {median_inter_spot2:.1f} um")

plt.tight_layout()
plt.savefig(os.path.join(outdir, "bootstrap_pass1_vs_pass2.png"), dpi=150)
print("Saved bootstrap_pass1_vs_pass2.png")

print("\nDone. Review bootstrap_crop_region.png and bootstrap_pass1_vs_pass2.png "
      "before proceeding to rasterization/nonlinear alignment.")
