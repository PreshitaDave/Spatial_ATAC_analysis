#!/usr/bin/env python3
"""
Rigid-transform (rotation + translation, scale locked to 1) alignment test.

Both ATAC (after the pixel->micron calibration fix) and Xenium (already true
microns) are now in the same physical units. If ATAC only covers a small
sub-region of the larger Xenium-imaged tissue, the correct registration
should be close to scale=1 (just rotate + translate the small ATAC patch to
its correct location) -- there's no legitimate reason for a global stretch.

affine_align (MOSAICField's default) fits a free, unconstrained 2x2 matrix
with no scale constraint, so it's free to stretch ATAC to better match
Xenium's PCA-feature landscape via the FGW optimal-transport correspondence,
which is what caused the ~2.2x inflation seen in earlier tests.

This script reuses the exact same FGW optimal-transport correspondence-finding
as MOSAICField's affine_align, but replaces the final least-squares affine fit
with a weighted orthogonal Procrustes (rigid: rotation + translation only,
scale fixed at 1) fit.

No cropping of Xenium is needed here -- a rigid transform can't over-stretch
ATAC to fill excess space, so this is tested directly against the FULL
(uncropped) Xenium point cloud.

Stops BEFORE rasterization/nonlinear alignment -- only reports the resulting
spot-spacing sanity check and plots, for review before committing to the
expensive downstream steps.
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
import ot
from scipy.spatial import distance, cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"
outdir = "./mosaicfield_outputs"

# ============================================================================
# Step A: recompute ATAC's true pre-affine calibrated coordinates fresh
# (identical to bootstrap_crop_reaffine.py Step A / mosaic_run2-2.ipynb Step -1)
# ============================================================================
print("Recomputing ATAC's true pre-affine calibrated coordinates...")
atac = sc.read_h5ad(
    "/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/"
    "ultima1942optimize/set1_vf25000-cr1-0-vi1-nc30/combined.h5ad"
)
atac_coords_px = np.asarray(atac.obsm["spatial"], dtype=float)
# Keep only the rightmost 50% of the ATAC array: there is a second, separate
# tissue section present in this combined.h5ad (a different physical piece,
# not the one imaged by this Xenium slide) that must be excluded -- both
# regions pass the 'on_off' tissue-QC flag, so that flag alone isn't enough
# to distinguish them. This positional trim is required, not optional.
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
X_design = sm.add_constant(tissue_positions[["array_row", "array_col"]])
fit_x = sm.OLS(tissue_positions["x_spatial"], X_design).fit()
fit_y = sm.OLS(tissue_positions["y_spatial"], X_design).fit()
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
# Step B: load Xenium (full, uncropped) -- reuse pass-1 processed object
# (already normalized/log1p/HVG/PCA'd, coords are native true-um, unchanged)
# ============================================================================
print("\nLoading Xenium (full, uncropped, already processed)...")
xenium_pass1 = ad.read_h5ad(os.path.join(outdir, "xenium_affine_aligned.h5ad"))
# xenium_pass1.obsm['spatial'] currently holds pass-1's SCALED/shifted coords
# (scalefac artifact from the OLD pipeline) -- reload native sdimx/sdimy instead
xenium_native = sc.read_h5ad(
    os.path.join(PROJECT_ROOT, "Data/04_analysis/multiomic/Xenium_488B/giotto_output",
                 "xenium_cells_features_coords.h5ad")
)
xenium_coords_native = xenium_native.obs[["sdimx", "sdimy"]].to_numpy().astype(float)
xenium_aligned = xenium_pass1.copy()  # reuse already-computed X_pca (HVG/PCA on same cells)
xenium_aligned.obsm["spatial"] = xenium_coords_native
xenium_center = xenium_aligned.obsm["spatial"].mean(axis=0)
xenium_aligned.obsm["spatial"] -= xenium_center

print(f"Xenium true native footprint: {xenium_aligned.shape[0]:,} cells, "
      f"x=[{xenium_aligned.obsm['spatial'][:,0].min():.0f}, {xenium_aligned.obsm['spatial'][:,0].max():.0f}], "
      f"y=[{xenium_aligned.obsm['spatial'][:,1].min():.0f}, {xenium_aligned.obsm['spatial'][:,1].max():.0f}]")

# ============================================================================
# Step C: FGW correspondence + RIGID (rotation+translation, scale=1) fit
# ============================================================================
def rigid_transformation(X, Y, P):
    """Weighted orthogonal Procrustes: rotation + translation only, scale=1.
    Minimizes sum_ij P_ij |R x_i + t - y_j|^2."""
    w = P.sum(axis=1)                      # (n,) row weights
    denom = np.maximum(w, 1e-12)
    Y_matched = (P @ Y) / denom[:, None]   # weighted-average target per source point
    w_sum = w.sum()
    x_bar = (w[:, None] * X).sum(axis=0) / w_sum
    y_bar = (w[:, None] * Y_matched).sum(axis=0) / w_sum
    Xc = X - x_bar
    Yc = Y_matched - y_bar
    H = (Xc * w[:, None]).T @ Yc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, d])
    R = Vt.T @ D @ U.T
    t = y_bar - R @ x_bar
    T = np.eye(3)
    T[:2, :2] = R
    T[:2, 2] = t
    transformed = (R @ X.T).T + t
    return transformed, T

def rigid_transform_apply(coords, T):
    homogeneous = np.vstack([coords.T, np.ones((1, coords.shape[0]))])
    return (T @ homogeneous)[:2, :].T

def FGW_rigid(X, Y, FX, FY, max_iter=10, alpha=0.9, device="cpu"):
    X = X.astype(np.float32)
    Y = Y.astype(np.float32)
    D_X = distance.cdist(FX, FX)
    D_Y = distance.cdist(FY, FY)
    D_X /= D_X.max()
    D_Y /= D_Y.max()
    D_X = torch.tensor(D_X).to(device)
    D_Y = torch.tensor(D_Y).to(device)

    T_net = np.eye(3)
    for it in range(max_iter):
        print("Iter:", it)
        C = distance.cdist(X, Y)
        C /= C[C > 0].max()
        C = torch.tensor(C).to(device)
        P = ot.gromov.fused_gromov_wasserstein(C, D_X, D_Y, alpha=alpha).to(torch.float32)
        X, T = rigid_transformation(X=X, Y=Y, P=P.cpu().numpy())
        T_net = T @ T_net

    C = distance.cdist(X, Y)
    C /= C[C > 0].max()
    C = torch.tensor(C).to(device)
    P = ot.gromov.fused_gromov_wasserstein(C, D_X, D_Y, alpha=alpha).to(torch.float32)
    return X, Y, T_net, P

np.random.seed(2)
n_subsample = 2000
atac_n = min(n_subsample, atac_aligned.shape[0])
xenium_n = min(n_subsample, xenium_aligned.shape[0])
atac_indices = np.random.choice(atac_aligned.shape[0], atac_n, replace=False)
xenium_indices = np.random.choice(xenium_aligned.shape[0], xenium_n, replace=False)
atac_sub = atac_aligned[atac_indices].copy()
xenium_sub = xenium_aligned[xenium_indices].copy()

print(f"\nRigid FGW alignment: true-calibrated ATAC ({atac_aligned.shape[0]:,} spots, "
      f"subsampled {atac_n:,}) vs FULL Xenium ({xenium_aligned.shape[0]:,} cells, "
      f"subsampled {xenium_n:,})...")

X0 = atac_sub.obsm["spatial"].astype(np.float32)
Y0 = xenium_sub.obsm["spatial"].astype(np.float32)
FX = atac_sub.obsm["X_pca"].astype(np.float32)
FY = xenium_sub.obsm["X_pca"].astype(np.float32)

_, _, T_rigid, _ = FGW_rigid(X0, Y0, FX, FY, max_iter=10, alpha=0.9)

atac_coords_full = np.asarray(atac_aligned.obsm["spatial"], dtype=float)
atac_coords_rigid = rigid_transform_apply(atac_coords_full, T_rigid)

print(f"\nRigid-fit ATAC footprint: x=[{atac_coords_rigid[:,0].min():.0f}, {atac_coords_rigid[:,0].max():.0f}], "
      f"y=[{atac_coords_rigid[:,1].min():.0f}, {atac_coords_rigid[:,1].max():.0f}]")

# ============================================================================
# Step D: sanity check -- inter-spot spacing + cells/spot after rigid fit
# ============================================================================
tree = cKDTree(atac_coords_rigid)
inter_spot_dist, _ = tree.query(atac_coords_rigid, k=2)
median_inter_spot = np.median(inter_spot_dist[:, 1])
print(f"\nRigid fit: median inter-spot distance = {median_inter_spot:.2f} um "
      f"(true pitch = 10 um; free-affine pass-1 gave 22.1 um)")

xenium_coords_full = xenium_aligned.obsm["spatial"]
nn_dist, nn_idx = tree.query(xenium_coords_full, k=1)
cells_per_spot = pd.Series(nn_idx).value_counts()
print(f"Rigid fit: median cells/spot = {cells_per_spot.median():.2f}, "
      f"mean = {cells_per_spot.mean():.2f} "
      f"(true expectation ~1-3; free-affine result was median=5, mean=5.2)")

# ============================================================================
# Plots
# ============================================================================
fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(xenium_coords_full[:, 0], xenium_coords_full[:, 1],
           s=0.5, alpha=0.15, c="royalblue", label="Xenium (full, native um)")
ax.scatter(atac_coords_rigid[:, 0], atac_coords_rigid[:, 1],
           s=1, alpha=0.6, c="red", label="ATAC (rigid fit, scale=1)")
ax.set_aspect("equal")
ax.legend(markerscale=10)
ax.set_title(f"Rigid (rotation+translation only) fit\nmedian inter-spot = {median_inter_spot:.1f} um")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "rigid_fit_overview.png"), dpi=150)
print("Saved rigid_fit_overview.png")

# zoomed view near ATAC centroid
cx, cy = atac_coords_rigid.mean(axis=0)
fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
ax.scatter(xenium_coords_full[:, 0], xenium_coords_full[:, 1],
           s=1, alpha=0.3, c="royalblue", label="Xenium")
ax.scatter(atac_coords_rigid[:, 0], atac_coords_rigid[:, 1],
           s=2, alpha=0.7, c="red", label="ATAC (rigid)")
ax.set_xlim(cx - 150, cx + 150)
ax.set_ylim(cy - 150, cy + 150)
ax.set_aspect("equal")
ax.legend(markerscale=5)
ax.set_title(f"Zoomed: rigid fit near ATAC centroid\nmedian inter-spot = {median_inter_spot:.1f} um")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "rigid_fit_zoom.png"), dpi=150)
print("Saved rigid_fit_zoom.png")

print("\nDone. Review rigid_fit_overview.png and rigid_fit_zoom.png before "
      "proceeding to rasterization/nonlinear alignment.")
