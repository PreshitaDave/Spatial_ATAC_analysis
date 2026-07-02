#!/usr/bin/env python3
"""
Two-pass rigid alignment + principled auto-crop, validated against the TRUE
raw (uncropped) Xenium data -- not the already-manually-cropped Giotto export.

Background: `xenium_cells_features_coords.h5ad` (used everywhere else in this
pipeline) is NOT the true raw Xenium capture -- it is already the output of a
manual, eyeballed crop in R/Giotto (`5_Xenium_giotto_analysis.Rmd` line 74:
`subsetGiottoLocs(xenium488B, x_min=5000, x_max=10600)`), done because Xenium
imaged more tissue than the ATAC chip covers. The earlier rigid_align_test.py
validation was accidentally circular (it used this already-cropped set as its
"full, uncropped" baseline). This script builds a genuine baseline by loading
10x's raw cells.csv.gz + cell_feature_matrix.h5 directly, bypassing Giotto/R
entirely, then derives an automatic crop from the rigid fit's own (correctly
sized, not stretched) result instead of the manual x_min/x_max eyeball.

Steps:
  0. Load true raw Xenium (cells.csv.gz + cell_feature_matrix.h5) -- 110,809
     cells, no crop.
  A. Pass-1 rigid FGW fit: true-calibrated, correctly-trimmed ATAC vs the full
     raw Xenium (110,809 cells) -- the real, honest baseline.
  B. Derive an automatic crop: NN-distance mask (not a bounding box) from
     pass-1's ATAC footprint, with a cutoff justified by measuring seed-to-seed
     registration jitter.
  C. Pass-2 rigid FGW fit: true-calibrated ATAC vs the auto-cropped Xenium,
     with PCA recomputed fresh on the cropped subset.

Stops BEFORE rasterization/nonlinear alignment -- only reports sanity checks
and plots, for review before committing to expensive downstream steps or
porting into mosaic_run2-2.ipynb.
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

sys.path.append(os.path.abspath("./MOSAICField/src/MOSAICField"))
sys.path.append(os.path.abspath("./MOSAICField"))

PROJECT_ROOT = "/projectnb/paxlab/presh/projects/spatial_atac"
XENIUM_RAW_DIR = ("/projectnb/paxlab/DATA/DriesSpatial/Xenium/"
                   "output-XETG00253__0048833__488B__20241217__182301")
outdir = "./mosaicfield_outputs/crop_fix_review"
os.makedirs(outdir, exist_ok=True)


# ============================================================================
# Shared functions (identical to rigid_align_test.py -- reused verbatim)
# ============================================================================
def rotate_coords(coords, angle_deg, center=None):
    theta = np.deg2rad(angle_deg)
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    if center is None:
        center = np.array([0, 0])
    return (coords - center) @ R.T + center


def rigid_transformation(X, Y, P):
    """Weighted orthogonal Procrustes: rotation + translation only, scale=1.
    Minimizes sum_ij P_ij |R x_i + t - y_j|^2."""
    w = P.sum(axis=1)
    denom = np.maximum(w, 1e-12)
    Y_matched = (P @ Y) / denom[:, None]
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


def rotation_angle_deg(T):
    return np.degrees(np.arctan2(T[1, 0], T[0, 0]))


def sanity_check(atac_coords, xenium_coords, label):
    tree = cKDTree(atac_coords)
    inter_spot_dist, _ = tree.query(atac_coords, k=2)
    median_inter_spot = np.median(inter_spot_dist[:, 1])
    nn_dist, nn_idx = tree.query(xenium_coords, k=1)
    cells_per_spot = pd.Series(nn_idx).value_counts()
    print(f"[{label}] median inter-spot distance = {median_inter_spot:.2f} um "
          f"(true pitch = 10 um)")
    print(f"[{label}] median cells/spot = {cells_per_spot.median():.2f}, "
          f"mean = {cells_per_spot.mean():.2f} (target ~1-3)")
    return median_inter_spot, cells_per_spot


# ============================================================================
# Step 0: load TRUE raw, uncropped Xenium (bypass Giotto/R entirely)
# ============================================================================
print("=" * 70)
print("STEP 0: Loading TRUE raw Xenium (cells.csv.gz + cell_feature_matrix.h5)")
print("=" * 70)

xenium_raw_expr = sc.read_10x_h5(os.path.join(XENIUM_RAW_DIR, "cell_feature_matrix.h5"))
xenium_raw_expr.var_names_make_unique()

cells_meta = pd.read_csv(os.path.join(XENIUM_RAW_DIR, "cells.csv.gz"))
cells_meta = cells_meta.set_index("cell_id")
cells_meta = cells_meta.loc[xenium_raw_expr.obs_names]

xenium_raw = xenium_raw_expr
xenium_raw.obsm["spatial"] = cells_meta[["x_centroid", "y_centroid"]].to_numpy()

print(f"Raw Xenium: {xenium_raw.shape[0]:,} cells, {xenium_raw.shape[1]} genes")
print(f"Raw Xenium spatial range: x=[{xenium_raw.obsm['spatial'][:,0].min():.0f}, "
      f"{xenium_raw.obsm['spatial'][:,0].max():.0f}], "
      f"y=[{xenium_raw.obsm['spatial'][:,1].min():.0f}, {xenium_raw.obsm['spatial'][:,1].max():.0f}]")
print(f"(Compare: existing Giotto-cropped set used elsewhere = 61,151 cells, "
      f"x=[5000,10600] -- should be a strict subset of the above)")

sc.pp.normalize_total(xenium_raw, inplace=True)
sc.pp.log1p(xenium_raw)
sc.pp.highly_variable_genes(xenium_raw, flavor="seurat", n_top_genes=2000, inplace=True, subset=True)
sc.pp.pca(xenium_raw, n_comps=50)

# ============================================================================
# Recompute ATAC's true pre-affine calibrated coordinates (unchanged logic)
# ============================================================================
print("\n" + "=" * 70)
print("Recomputing ATAC's true pre-affine calibrated coordinates")
print("=" * 70)
atac = sc.read_h5ad(
    "/projectnb/paxlab/DATA/DriesSpatial/Atlasxomics/D1942/D1942_deepseq/"
    "ultima1942optimize/set1_vf25000-cr1-0-vi1-nc30/combined.h5ad"
)
atac_coords_px = np.asarray(atac.obsm["spatial"], dtype=float)
# Keep only the rightmost 50%: there is a second, separate physical tissue
# section in this combined.h5ad that must be excluded from this alignment
# (confirmed by user) -- this positional trim is required, not optional.
x_min, x_max = atac_coords_px[:, 0].min(), atac_coords_px[:, 0].max()
x_threshold = x_max - 0.5 * (x_max - x_min)
mask_rightmost = atac_coords_px[:, 0] > x_threshold
atac_trimmed = atac[mask_rightmost].copy()
atac_trimmed.obsm["spatial"] = atac_coords_px[mask_rightmost]
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

# Center Xenium too (so both are in a comparable, origin-centered frame like
# the rest of the pipeline)
xenium_center = xenium_raw.obsm["spatial"].mean(axis=0)
xenium_raw.obsm["spatial"] = xenium_raw.obsm["spatial"] - xenium_center

# ============================================================================
# Step A: pass-1 rigid fit -- true ATAC vs TRUE FULL raw Xenium (110,809 cells)
# ============================================================================
print("\n" + "=" * 70)
print("STEP A: Pass-1 rigid fit -- ATAC vs TRUE FULL raw Xenium (no crop)")
print("=" * 70)

def run_rigid_fit(atac_obj, xenium_obj, seed, n_subsample=2000, max_iter=10, alpha=0.9):
    rng = np.random.default_rng(seed)
    atac_n = min(n_subsample, atac_obj.shape[0])
    xenium_n = min(n_subsample, xenium_obj.shape[0])
    atac_indices = rng.choice(atac_obj.shape[0], atac_n, replace=False)
    xenium_indices = rng.choice(xenium_obj.shape[0], xenium_n, replace=False)
    atac_sub = atac_obj[atac_indices].copy()
    xenium_sub = xenium_obj[xenium_indices].copy()
    X0 = atac_sub.obsm["spatial"].astype(np.float32)
    Y0 = xenium_sub.obsm["spatial"].astype(np.float32)
    FX = atac_sub.obsm["X_pca"].astype(np.float32)
    FY = xenium_sub.obsm["X_pca"].astype(np.float32)
    _, _, T, _ = FGW_rigid(X0, Y0, FX, FY, max_iter=max_iter, alpha=alpha)
    return T

T1 = run_rigid_fit(atac_aligned, xenium_raw, seed=2)
atac_coords_full = np.asarray(atac_aligned.obsm["spatial"], dtype=float)
atac_coords_pass1 = rigid_transform_apply(atac_coords_full, T1)
print(f"Pass-1 ATAC footprint: x=[{atac_coords_pass1[:,0].min():.0f}, {atac_coords_pass1[:,0].max():.0f}], "
      f"y=[{atac_coords_pass1[:,1].min():.0f}, {atac_coords_pass1[:,1].max():.0f}]")
sanity_check(atac_coords_pass1, xenium_raw.obsm["spatial"], "Pass 1 (true full Xenium, no crop)")

# ============================================================================
# Step A(jitter): estimate registration jitter across seeds, to justify the
# crop cutoff in Step B
# ============================================================================
print("\n" + "=" * 70)
print("Estimating seed-to-seed registration jitter (to justify crop cutoff)")
print("=" * 70)
centroids = []
angles = []
for seed in [10, 11, 12]:
    T_s = run_rigid_fit(atac_aligned, xenium_raw, seed=seed)
    coords_s = rigid_transform_apply(atac_coords_full, T_s)
    centroids.append(coords_s.mean(axis=0))
    angles.append(rotation_angle_deg(T_s))
    print(f"  seed {seed}: centroid={coords_s.mean(axis=0)}, rotation={angles[-1]:.2f} deg")
centroids = np.array(centroids)
centroid_jitter = np.linalg.norm(centroids - centroids.mean(axis=0), axis=1).max()
angle_jitter = np.ptp(angles)
print(f"Max centroid displacement across seeds: {centroid_jitter:.1f} um")
print(f"Rotation angle spread across seeds: {angle_jitter:.2f} deg")

CROP_CUTOFF_UM = max(100.0, 8.0 * centroid_jitter)
print(f"\nChosen crop cutoff (>= 8x centroid jitter, floor 100um): {CROP_CUTOFF_UM:.1f} um")

# ============================================================================
# Step B: NN-distance-based auto-crop from pass-1's footprint
# ============================================================================
print("\n" + "=" * 70)
print("STEP B: Deriving auto-crop mask from pass-1 footprint")
print("=" * 70)
tree_pass1 = cKDTree(atac_coords_pass1)
xenium_coords_full = xenium_raw.obsm["spatial"]
dist_to_atac, _ = tree_pass1.query(xenium_coords_full, k=1)
crop_mask = dist_to_atac <= CROP_CUTOFF_UM
print(f"Xenium cells kept after auto-crop (<= {CROP_CUTOFF_UM:.0f}um from pass-1 ATAC): "
      f"{crop_mask.sum():,} / {len(crop_mask):,} ({100*crop_mask.mean():.1f}%)")

# Compare against the old manual crop [5000, 10600] in native (uncentered) sdimx frame
xenium_native_x = xenium_coords_full[:, 0] + xenium_center[0]
old_manual_mask = (xenium_native_x >= 5000) & (xenium_native_x <= 10600)
overlap = (crop_mask & old_manual_mask).sum()
print(f"Overlap with old manual crop [5000,10600]: {overlap:,} cells "
      f"({100*overlap/max(crop_mask.sum(),1):.1f}% of auto-crop, "
      f"{100*overlap/max(old_manual_mask.sum(),1):.1f}% of manual crop)")
print(f"Auto-crop native-x range: [{xenium_native_x[crop_mask].min():.0f}, "
      f"{xenium_native_x[crop_mask].max():.0f}]")

fig, ax = plt.subplots(figsize=(12, 5))
ax.scatter(xenium_coords_full[~crop_mask, 0], xenium_coords_full[~crop_mask, 1],
           s=0.5, alpha=0.15, c="lightgrey", label="Xenium (excluded, excess tissue)")
ax.scatter(xenium_coords_full[crop_mask, 0], xenium_coords_full[crop_mask, 1],
           s=0.5, alpha=0.3, c="royalblue", label="Xenium (auto-crop, kept)")
ax.scatter(atac_coords_pass1[:, 0], atac_coords_pass1[:, 1],
           s=1, alpha=0.6, c="red", label="ATAC (pass-1 rigid footprint)")
ax.set_aspect("equal")
ax.legend(markerscale=10)
ax.set_title(f"Auto-crop from pass-1 rigid fit (NN cutoff={CROP_CUTOFF_UM:.0f}um)\n"
             f"vs. TRUE full raw Xenium ({xenium_raw.shape[0]:,} cells)")
plt.tight_layout()
plt.savefig(os.path.join(outdir, "rigid_crop_refine_step_b.png"), dpi=150)
print("Saved rigid_crop_refine_step_b.png")

# ============================================================================
# Step C: pass-2 rigid fit on auto-cropped Xenium, PCA recomputed fresh
# ============================================================================
print("\n" + "=" * 70)
print("STEP C: Pass-2 rigid fit on auto-cropped Xenium (fresh PCA)")
print("=" * 70)

xenium_cropped_raw = xenium_raw_expr[crop_mask].copy()
xenium_cropped_raw.obsm["spatial"] = xenium_coords_full[crop_mask]
sc.pp.normalize_total(xenium_cropped_raw, inplace=True)
sc.pp.log1p(xenium_cropped_raw)
sc.pp.highly_variable_genes(xenium_cropped_raw, flavor="seurat", n_top_genes=2000, inplace=True, subset=True)
sc.pp.pca(xenium_cropped_raw, n_comps=50)

T2 = run_rigid_fit(atac_aligned, xenium_cropped_raw, seed=2)
atac_coords_pass2 = rigid_transform_apply(atac_coords_full, T2)
print(f"Pass-2 ATAC footprint: x=[{atac_coords_pass2[:,0].min():.0f}, {atac_coords_pass2[:,0].max():.0f}], "
      f"y=[{atac_coords_pass2[:,1].min():.0f}, {atac_coords_pass2[:,1].max():.0f}]")
sanity_check(atac_coords_pass2, xenium_cropped_raw.obsm["spatial"], "Pass 2 (auto-cropped Xenium)")

centroid_shift = np.linalg.norm(atac_coords_pass2.mean(axis=0) - atac_coords_pass1.mean(axis=0))
angle_shift = abs(rotation_angle_deg(T2) - rotation_angle_deg(T1))
print(f"\nPass 1 -> Pass 2 centroid shift: {centroid_shift:.1f} um "
      f"(flag if >> seed jitter of {centroid_jitter:.1f} um)")
print(f"Pass 1 -> Pass 2 rotation shift: {angle_shift:.2f} deg "
      f"(flag if >> seed jitter of {angle_jitter:.2f} deg)")

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
axes[0].scatter(xenium_coords_full[:, 0], xenium_coords_full[:, 1],
                s=0.3, alpha=0.15, c="royalblue", label="Xenium (full raw)")
axes[0].scatter(atac_coords_pass1[:, 0], atac_coords_pass1[:, 1],
                s=1, alpha=0.6, c="red", label="ATAC (pass 1)")
axes[0].set_aspect("equal"); axes[0].legend(markerscale=10)
axes[0].set_title("Pass 1: rigid fit vs full raw Xenium")

axes[1].scatter(xenium_cropped_raw.obsm["spatial"][:, 0], xenium_cropped_raw.obsm["spatial"][:, 1],
                s=0.5, alpha=0.2, c="royalblue", label="Xenium (auto-cropped)")
axes[1].scatter(atac_coords_pass2[:, 0], atac_coords_pass2[:, 1],
                s=1, alpha=0.6, c="red", label="ATAC (pass 2)")
axes[1].set_aspect("equal"); axes[1].legend(markerscale=10)
axes[1].set_title("Pass 2: rigid fit vs auto-cropped Xenium")

plt.tight_layout()
plt.savefig(os.path.join(outdir, "rigid_crop_refine_pass1_vs_pass2.png"), dpi=150)
print("Saved rigid_crop_refine_pass1_vs_pass2.png")

print("\nDone. Review rigid_crop_refine_step_b.png and "
      "rigid_crop_refine_pass1_vs_pass2.png before proceeding.")
