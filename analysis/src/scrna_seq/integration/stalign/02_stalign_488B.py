#!/usr/bin/env python3
"""
02_stalign_488B.py
STAlign / PASTE2: RNA-guided spatial re-alignment of ATAC pseudo-RNA → Xenium.
Each ATAC spot has imputed RNA (GeneIntegrationMatrix from Option 1) + spatial coords.
STAlign finds an optimal transport plan that minimises expression + spatial discrepancy,
producing a refined coordinate estimate for each ATAC spot.

Compare against MOSAICField nonlinear alignment to quantify the improvement.

Prerequisites:
  - 00_export_for_stalign.R (writes atac_pseudorna + spatial coords)
  - Xenium h5ad at mosaicfield_outputs/xenium_nonlinear_aligned.h5ad

Install (if needed): pip install paste-bio   (PASTE2, closely related to STAlign)
STAlign paper: Shi et al. 2023  https://doi.org/10.1101/2023.08.13.553130
PASTE2 paper:  Chen et al. 2023 https://doi.org/10.1101/2023.06.07.544059
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix

TISSUE = "488B"
BASE_OUT = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration"
STALIGN_DIR = os.path.join(BASE_OUT, "stalign", TISSUE)
XENIUM_H5AD = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/mosaicfield_outputs/xenium_nonlinear_aligned.h5ad"

os.makedirs(STALIGN_DIR, exist_ok=True)

print(f"=== STAlign: tissue {TISSUE} ===")

# ── Load ATAC pseudo-RNA + spatial coords ─────────────────────────────────────
print("Loading ATAC pseudo-RNA (GeneIntegrationMatrix)...")
atac_mat = csr_matrix(mmread(os.path.join(STALIGN_DIR, f"atac_pseudorna_{TISSUE}.mtx")).T)
atac_genes = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_gene_names_{TISSUE}.csv"))["gene"].tolist()
atac_cells = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_cell_names_{TISSUE}.csv"))["cell"].tolist()
atac_coords = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_spatial_coords_{TISSUE}.csv"), index_col=0)
atac_meta   = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_metadata_{TISSUE}.csv"), index_col=0)

atac = ad.AnnData(X=atac_mat)
atac.obs_names = atac_cells
atac.var_names = atac_genes

# Attach spatial coordinates
coord_index = atac_coords.index if atac_coords.index.dtype == object else pd.Index(atac_cells)
atac.obsm["spatial"] = atac_coords.loc[atac_cells][["x_um","y_um"]].values
atac.obs["predicted_type"] = atac_meta.loc[atac_cells, "predicted_type"].values if "predicted_type" in atac_meta.columns else "Unknown"

print(f"ATAC pseudo-RNA: {atac.shape[0]} spots × {atac.shape[1]} genes")

# ── Load Xenium spatial reference ─────────────────────────────────────────────
print("Loading Xenium spatial h5ad...")
xenium = sc.read_h5ad(XENIUM_H5AD)
print(f"Xenium: {xenium.shape[0]} cells × {xenium.shape[1]} genes")

# ── Preprocessing ─────────────────────────────────────────────────────────────
shared_genes = list(set(atac.var_names) & set(xenium.var_names))
print(f"Shared genes (ATAC pseudo-RNA ∩ Xenium): {len(shared_genes)}")

# Subset to shared
atac_s = atac[:, shared_genes].copy()
xen_s  = xenium[:, shared_genes].copy()

# Normalize both to same scale
sc.pp.normalize_total(atac_s, target_sum=1e4)
sc.pp.log1p(atac_s)
sc.pp.normalize_total(xen_s, target_sum=1e4)
sc.pp.log1p(xen_s)

# Subsample Xenium if very large (PASTE OT is O(n²))
MAX_XENIUM = 5000
if xen_s.n_obs > MAX_XENIUM:
    print(f"Subsampling Xenium to {MAX_XENIUM} cells for OT computation...")
    sc.pp.subsample(xen_s, n_obs=MAX_XENIUM, random_state=42)

print(f"Using {atac_s.n_obs} ATAC spots × {xen_s.n_obs} Xenium cells")

# ── Try PASTE2 / STAlign OT alignment ─────────────────────────────────────────
try:
    import paste as pst

    print("Running PASTE2 pairwise alignment (ATAC pseudo-RNA ↔ Xenium)...")
    # alpha: weight of spatial term (0=expression only, 1=spatial only); 0.1 is typical
    pi, log_dict = pst.pairwise_align(
        atac_s,
        xen_s,
        alpha          = 0.1,
        dissimilarity  = "kl",
        use_rep        = None,
        return_obj     = True,
        numItermax     = 1000,
        verbose        = True
    )

    print(f"Transport plan shape: {pi.shape}  (ATAC spots × Xenium cells)")

    # RNA-guided coordinates: weighted average of Xenium positions by transport plan
    xen_coords = xen_s.obsm["spatial"]   # Xenium spatial coords (µm)
    row_sums = pi.sum(axis=1, keepdims=True) + 1e-9
    new_coords = (pi @ xen_coords) / row_sums
    print(f"New ATAC coordinates computed (RNA-guided): {new_coords.shape}")

    # ── Compare MOSAICField vs STAlign ────────────────────────────────────────
    orig_coords = atac_s.obsm["spatial"]
    displacement = new_coords - orig_coords
    residual_mag = np.sqrt((displacement ** 2).sum(axis=1))
    print(f"Displacement from MOSAICField → STAlign:")
    print(f"  Median: {np.median(residual_mag):.1f} µm")
    print(f"  Mean:   {np.mean(residual_mag):.1f} µm")
    print(f"  90th:   {np.percentile(residual_mag, 90):.1f} µm")

    # Save outputs
    new_coord_df = pd.DataFrame(new_coords, columns=["x_um_stalign","y_um_stalign"],
                                 index=atac_s.obs_names)
    new_coord_df["x_um_orig"] = orig_coords[:, 0]
    new_coord_df["y_um_orig"] = orig_coords[:, 1]
    new_coord_df["displacement_um"] = residual_mag
    new_coord_df.to_csv(os.path.join(STALIGN_DIR, f"stalign_atac_new_coords_{TISSUE}.csv"))
    print(f"Saved new coordinates → stalign_atac_new_coords_{TISSUE}.csv")

    np.savez_compressed(os.path.join(STALIGN_DIR, f"stalign_transport_plan_{TISSUE}.npz"),
                        pi=pi, atac_obs=np.array(atac_s.obs_names),
                        xenium_obs=np.array(xen_s.obs_names))
    print("Saved transport plan.")

    # Quick summary stats for comparison script
    summary = {
        "tissue": TISSUE,
        "n_atac_spots": atac_s.n_obs,
        "n_xenium_cells": xen_s.n_obs,
        "n_shared_genes": len(shared_genes),
        "alpha": 0.1,
        "median_displacement_um": float(np.median(residual_mag)),
        "mean_displacement_um": float(np.mean(residual_mag)),
        "p90_displacement_um": float(np.percentile(residual_mag, 90))
    }
    pd.DataFrame([summary]).to_csv(
        os.path.join(STALIGN_DIR, f"stalign_summary_{TISSUE}.csv"), index=False)

except ImportError:
    print("WARNING: paste-bio not installed. Install with: pip install paste-bio")
    print("Alternatively: pip install stalign   (if available in your env)")
    print("\nWriting stub outputs for pipeline continuity...")
    stub = pd.DataFrame({
        "x_um_orig": atac_s.obsm["spatial"][:, 0],
        "y_um_orig": atac_s.obsm["spatial"][:, 1],
        "x_um_stalign": atac_s.obsm["spatial"][:, 0],
        "y_um_stalign": atac_s.obsm["spatial"][:, 1],
        "displacement_um": np.zeros(atac_s.n_obs)
    }, index=atac_s.obs_names)
    stub.to_csv(os.path.join(STALIGN_DIR, f"stalign_atac_new_coords_{TISSUE}.csv"))
    pd.DataFrame([{"tissue": TISSUE, "status": "stub_paste_not_installed"}]).to_csv(
        os.path.join(STALIGN_DIR, f"stalign_summary_{TISSUE}.csv"), index=False)

print(f"\nSTAlign ({TISSUE}) complete.")
