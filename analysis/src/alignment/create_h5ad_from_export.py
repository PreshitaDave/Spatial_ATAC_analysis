#!/usr/bin/env python3
"""
Build h5ad object from exported ArchR data (CSVs + MTX).
Usage: python3 create_h5ad_from_export.py <export_dir> <output_path>
"""

import sys
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

def create_atac_h5ad(export_dir, output_path):
    """
    Create AnnData h5ad from ArchR export directory.

    Expected files:
      - atac_peak_matrix.mtx (peaks × cells in MTX format)
      - atac_peak_names.csv
      - atac_cell_names.csv
      - atac_coords.csv
      - atac_obs.csv
      - atac_var.csv
      - atac_peak_coords.csv (optional, for uns)
    """

    export_dir = Path(export_dir)
    print(f"Reading from: {export_dir}")

    # --- Load matrix (peaks × cells in MTX, transpose to cells × peaks) ---
    print("Loading matrix...")
    mat = mmread(export_dir / "atac_peak_matrix.mtx")
    print(f"  MTX shape (peaks × cells): {mat.shape}")

    # Convert to CSR and transpose
    if not isinstance(mat, csr_matrix):
        mat = csr_matrix(mat)
    mat = mat.T  # now cells × peaks
    print(f"  After transpose (cells × peaks): {mat.shape}")

    # --- Load names ---
    print("Loading names...")
    peak_names = pd.read_csv(export_dir / "atac_peak_names.csv")['peak'].values
    cell_ids = pd.read_csv(export_dir / "atac_cell_names.csv")['cell_id'].values

    print(f"  Peaks: {len(peak_names)}")
    print(f"  Cells: {len(cell_ids)}")

    # Verify dimensions
    assert mat.shape[0] == len(cell_ids), f"Cell count mismatch: {mat.shape[0]} vs {len(cell_ids)}"
    assert mat.shape[1] == len(peak_names), f"Peak count mismatch: {mat.shape[1]} vs {len(peak_names)}"

    # --- Create AnnData object ---
    print("Creating AnnData object...")
    adata = ad.AnnData(X=mat)
    adata.obs_names = cell_ids
    adata.var_names = peak_names

    # --- Load obs (metadata per cell/spot) ---
    print("Loading obs...")
    obs_df = pd.read_csv(export_dir / "atac_obs.csv", index_col=0)
    # Reindex to match AnnData order
    obs_df = obs_df.loc[adata.obs_names]
    adata.obs = obs_df
    print(f"  Obs columns: {list(adata.obs.columns)}")

    # --- Load var (statistics per peak) ---
    print("Loading var...")
    var_df = pd.read_csv(export_dir / "atac_var.csv", index_col=0)
    # Reindex to match AnnData order
    var_df = var_df.loc[adata.var_names]
    adata.var = var_df
    print(f"  Var columns: {list(adata.var.columns)}")

    # --- Load spatial coordinates ---
    print("Loading spatial coordinates...")
    coords_df = pd.read_csv(export_dir / "atac_coords.csv", index_col=0)
    coords_df = coords_df.loc[adata.obs_names]

    # Extract x, y columns (or handle different names)
    if 'x' in coords_df.columns and 'y' in coords_df.columns:
        spatial = coords_df[['x', 'y']].values
    elif 'xcor' in coords_df.columns and 'ycor' in coords_df.columns:
        spatial = coords_df[['xcor', 'ycor']].values
    else:
        # Try to find numeric columns that look like coordinates
        numeric_cols = coords_df.select_dtypes(include=[np.number]).columns.tolist()
        if len(numeric_cols) >= 2:
            spatial = coords_df[numeric_cols[:2]].values
        else:
            raise ValueError("Could not find coordinate columns in coords CSV")

    adata.obsm['spatial'] = spatial
    print(f"  Spatial coords: {spatial.shape}")
    print(f"    x=[{spatial[:,0].min():.1f}, {spatial[:,0].max():.1f}]")
    print(f"    y=[{spatial[:,1].min():.1f}, {spatial[:,1].max():.1f}]")

    # --- Add peak genomic coordinates to uns ---
    if (export_dir / "atac_peak_coords.csv").exists():
        print("Loading peak genomic coordinates...")
        peak_coords = pd.read_csv(export_dir / "atac_peak_coords.csv", index_col=0)
        peak_coords = peak_coords.loc[adata.var_names]
        adata.uns['peak_coords'] = peak_coords.to_dict('list')
        print(f"  Peaks with genomic coords: {len(peak_coords)}")

    # --- Basic QC ---
    print("\nData summary:")
    print(f"  Shape: {adata.shape}")
    print(f"  X sparsity: {1 - adata.X.nnz / np.prod(adata.shape):.3f}")
    print(f"  X range: [{adata.X.data.min():.3f}, {adata.X.data.max():.3f}]")
    print(f"  Obs columns: {list(adata.obs.columns)}")
    print(f"  Var columns: {list(adata.var.columns)}")

    # --- Save ---
    print(f"\nWriting to: {output_path}")
    adata.write_h5ad(output_path)
    print("Done!")

    return adata

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 create_h5ad_from_export.py <export_dir> [output_path]")
        print("  export_dir: directory with CSVs and MTX from create_h5ad_from_ArchR.R")
        print("  output_path: path to save h5ad (default: atac_from_archR.h5ad)")
        sys.exit(1)

    export_dir = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) > 2 else "atac_from_archR.h5ad"

    adata = create_atac_h5ad(export_dir, output_path)
