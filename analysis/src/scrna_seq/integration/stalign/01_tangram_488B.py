#!/usr/bin/env python3
"""
01_tangram_488B.py
Tangram: map non-spatial scRNA-seq cells onto Xenium spatial reference (tissue 488B).
Optimal transport maps each scRNA cell to a Xenium location, giving spatially resolved
cell type composition at Xenium resolution, aggregated per ATAC spot.

Prerequisites:
  - 00_export_for_stalign.R must have been run (provides scrna and atac_pseudorna exports)
  - Xenium h5ad at mosaicfield_outputs/xenium_nonlinear_aligned.h5ad

Install (if not in env): pip install tangram-sc
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix

TISSUE = "488B"
BASE_OUT = "/projectnb/paxlab/presh/projects/spatial_atac/Data/04_analysis/multiomic/scrna_integration"
STALIGN_DIR = os.path.join(BASE_OUT, "stalign", TISSUE)
XENIUM_H5AD = "/projectnb/paxlab/presh/projects/spatial_atac/analysis/src/alignment/mosaicfield_outputs/xenium_nonlinear_aligned.h5ad"

os.makedirs(STALIGN_DIR, exist_ok=True)

print(f"=== Tangram: tissue {TISSUE} ===")

# ── Load scRNA-seq ────────────────────────────────────────────────────────────
print("Loading scRNA-seq...")
scrna_mat = csr_matrix(mmread(os.path.join(STALIGN_DIR, f"scrna_normexpr_{TISSUE}.mtx")).T)
scrna_genes = pd.read_csv(os.path.join(STALIGN_DIR, f"scrna_gene_names_{TISSUE}.csv"))["gene"].tolist()
scrna_cells = pd.read_csv(os.path.join(STALIGN_DIR, f"scrna_cell_names_{TISSUE}.csv"))["cell"].tolist()
scrna_meta = pd.read_csv(os.path.join(STALIGN_DIR, f"scrna_metadata_{TISSUE}.csv"), index_col=0)

scrna = ad.AnnData(X=scrna_mat, obs=scrna_meta.loc[scrna_cells] if len(scrna_meta) == len(scrna_cells) else scrna_meta)
scrna.obs_names = scrna_cells
scrna.var_names = scrna_genes
scrna.obs["cell_type"] = scrna.obs.get("cell_type", "Unknown").fillna("Unknown")
print(f"scRNA: {scrna.shape[0]} cells × {scrna.shape[1]} genes")
print("Cell types:", dict(scrna.obs["cell_type"].value_counts()))

# ── Load Xenium spatial reference ─────────────────────────────────────────────
print("Loading Xenium spatial h5ad...")
xenium = sc.read_h5ad(XENIUM_H5AD)
print(f"Xenium: {xenium.shape[0]} cells × {xenium.shape[1]} genes")

# ── Find marker genes for mapping ────────────────────────────────────────────
# Use HVGs shared between scRNA and Xenium
sc.pp.normalize_total(scrna)
sc.pp.log1p(scrna)

sc.pp.normalize_total(xenium)
sc.pp.log1p(xenium)

shared_genes = list(set(scrna.var_names) & set(xenium.var_names))
print(f"Shared genes: {len(shared_genes)}")

# Use cell type markers as training genes for Tangram (more reliable than all genes)
marker_genes = [
    "EPCAM","KRT8","KRT18","KRT14","ESR1","PGR","ERBB2",     # Tumor
    "CD3E","CD3D","CD8A","CD4","GZMB","PRF1",                 # T cell
    "GNLY","NKG7","KLRD1",                                    # NK
    "MS4A1","CD79A","CD79B",                                  # B cell
    "CD14","LYZ","FCGR3A","CST3","CD68","CD163",              # Myeloid
    "COL1A1","COL1A2","DCN","FAP","VIM",                      # Fibroblast
    "PECAM1","VWF","CDH5","CLDN5"                             # Endothelial
]
training_genes = [g for g in marker_genes if g in shared_genes]
print(f"Training genes for Tangram: {len(training_genes)}")

try:
    import tangram as tg

    tg.pp_adatas(scrna, xenium, genes=training_genes)

    print("Running Tangram mapping (cells → Xenium space)...")
    ad_map = tg.map_cells_to_space(
        adata_sc     = scrna,
        adata_sp     = xenium,
        mode         = "cells",
        density_prior = "rna_count_based",
        num_epochs   = 500,
        device       = "cpu",
        verbose      = True
    )

    # Transfer cell type annotations onto Xenium cells
    tg.project_cell_annotations(ad_map, xenium, annotation="cell_type")
    print("Cell type probabilities transferred to Xenium.")

    # Save Xenium with deconvolved cell types
    xenium_out = os.path.join(STALIGN_DIR, f"tangram_xenium_celltypes_{TISSUE}.h5ad")
    xenium.write_h5ad(xenium_out)
    print(f"Saved: {xenium_out}")

    # ── Aggregate cell type composition per ATAC spot ──────────────────────
    # Load ATAC spatial coords
    atac_coords = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_spatial_coords_{TISSUE}.csv"))
    atac_meta   = pd.read_csv(os.path.join(STALIGN_DIR, f"atac_metadata_{TISSUE}.csv"))

    xenium_coords = pd.DataFrame(xenium.obsm["spatial"], columns=["x_um","y_um"],
                                 index=xenium.obs_names)

    # For each Xenium cell, find nearest ATAC spot (within 50 µm)
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=1, algorithm="ball_tree")
    nn.fit(atac_coords[["x_um","y_um"]].values)
    dists, idxs = nn.kneighbors(xenium_coords.values)

    xenium.obs["nearest_atac_spot"] = atac_coords["cell"].values[idxs[:, 0]]
    xenium.obs["nearest_atac_dist"]  = dists[:, 0]
    xenium_close = xenium[xenium.obs["nearest_atac_dist"] < 50].copy()

    # Cell type probability columns
    ct_cols = [c for c in xenium.obs.columns if c.startswith("cell_type_")]
    if not ct_cols:
        ct_cols = [c for c in xenium.obs.columns if c in scrna.obs["cell_type"].unique()]

    comp = xenium_close.obs.groupby("nearest_atac_spot")[ct_cols].mean()
    comp.to_csv(os.path.join(STALIGN_DIR, f"tangram_atac_spot_composition_{TISSUE}.csv"))
    print(f"Saved ATAC spot cell type composition ({comp.shape[0]} spots).")

    # Save transport matrix
    np.save(os.path.join(STALIGN_DIR, f"tangram_transport_plan_{TISSUE}.npy"),
            ad_map.X.toarray() if hasattr(ad_map.X, "toarray") else ad_map.X)
    print("Saved transport plan.")

except ImportError:
    print("WARNING: tangram-sc not installed. Install with: pip install tangram-sc")
    print("Writing stub outputs for pipeline continuity...")
    pd.DataFrame(columns=["nearest_atac_spot","Tumor","T_cell","B_cell","Myeloid","NK_cell","Fibroblast","Endothelial"]).to_csv(
        os.path.join(STALIGN_DIR, f"tangram_atac_spot_composition_{TISSUE}.csv"), index=False)

print(f"\nTangram ({TISSUE}) complete.")
