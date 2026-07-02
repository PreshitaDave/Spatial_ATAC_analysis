#!/usr/bin/env python3
"""
Export all data needed for R loss-function comparison.
Computes the affine-only mapping (baseline, unique to this script -- no equivalent
exists elsewhere) and reuses the nonlinear-stage mapping already computed and saved
by mosaic_run2-2.ipynb's own STEP 4b (xenium_to_atac_mapping_voronoi.csv), instead of
recomputing a duplicate cKDTree query for the nonlinear stage.
Run from the mosaic_field directory.
"""
import os, numpy as np, pandas as pd, anndata as ad, re
from scipy.spatial import cKDTree
from scipy.sparse import issparse, csr_matrix
from scipy.io import mmwrite

outdir = "./mosaicfield_outputs"
rdir = os.path.join(outdir, "r_comparison")
os.makedirs(rdir, exist_ok=True)

# --- Load aligned objects ---
atac_nl = ad.read_h5ad(os.path.join(outdir, 'atac_nonlinear_aligned.h5ad'))
xenium = ad.read_h5ad(os.path.join(outdir, 'xenium_nonlinear_aligned.h5ad'))

warped_um = atac_nl.obsm['spatial']
affine_um = atac_nl.obsm['spatial_affine']
xenium_um = xenium.obsm['spatial']

atac_names = np.array(atac_nl.obs_names)
xenium_names = np.array(xenium.obs_names)

print(f"ATAC spots: {len(atac_names):,}")
print(f"Xenium cells: {len(xenium_names):,}")

# --- Compute mappings for a given set of ATAC coordinates ---
def compute_mappings(atac_coords, xenium_coords, atac_names, xenium_names,
                     cutoff=20.0, soft_cutoff=30.0):
    tree = cKDTree(atac_coords)
    nn_dist, nn_idx = tree.query(xenium_coords, k=1)

    nn_map = pd.DataFrame({
        'xenium_cell': xenium_names,
        'xenium_x_um': xenium_coords[:, 0],
        'xenium_y_um': xenium_coords[:, 1],
        'atac_spot': atac_names[nn_idx],
        'atac_spot_idx': nn_idx,
        'distance_um': nn_dist,
        'assigned': nn_dist <= cutoff,
    })

    vor_map = pd.DataFrame({
        'xenium_cell': xenium_names,
        'xenium_x_um': xenium_coords[:, 0],
        'xenium_y_um': xenium_coords[:, 1],
        'atac_spot': atac_names[nn_idx],
        'atac_spot_idx': nn_idx,
        'distance_um': nn_dist,
        'distant': nn_dist > soft_cutoff,
    })

    n_nn = nn_map['assigned'].sum()
    n_conf = (~vor_map['distant']).sum()
    print(f"  NN ({cutoff}um): {n_nn:,} assigned ({100*n_nn/len(nn_map):.1f}%)")
    print(f"  Voronoi (+{soft_cutoff}um): {len(vor_map):,} total, {n_conf:,} confident")
    return nn_map, vor_map

print("\n--- Affine-only alignment (baseline; no equivalent computed elsewhere) ---")
aff_nn, aff_vor = compute_mappings(affine_um, xenium_um, atac_names, xenium_names)

print("\n--- Nonlinear alignment (reusing mosaic_run2-2.ipynb STEP 4b output) ---")
_voronoi_csv = os.path.join(outdir, "xenium_to_atac_mapping_voronoi.csv")
if not os.path.exists(_voronoi_csv):
    raise FileNotFoundError(
        f"{_voronoi_csv} not found -- run mosaic_run2-2.ipynb's STEP 4b first "
        f"(this script no longer recomputes the nonlinear-stage mapping itself)."
    )
_nl_map = pd.read_csv(_voronoi_csv)
CUTOFF, SOFT_CUTOFF = 20.0, 30.0
nl_nn = _nl_map.copy()
nl_nn["assigned"] = nl_nn["distance_um"] <= CUTOFF
nl_nn = nl_nn[["xenium_cell", "xenium_x_um", "xenium_y_um", "atac_spot",
               "atac_spot_idx", "distance_um", "assigned"]]
nl_vor = _nl_map[["xenium_cell", "xenium_x_um", "xenium_y_um", "atac_spot",
                  "atac_spot_idx", "distance_um", "distant"]].copy()
n_nn = nl_nn["assigned"].sum()
n_conf = (~nl_vor["distant"]).sum()
print(f"  NN ({CUTOFF}um): {n_nn:,} assigned ({100*n_nn/len(nl_nn):.1f}%)")
print(f"  Voronoi (+{SOFT_CUTOFF}um): {len(nl_vor):,} total, {n_conf:,} confident")

# --- Save 4 mapping CSVs ---
aff_nn.to_csv(os.path.join(rdir, 'mapping_affine_nn.csv'), index=False)
aff_vor.to_csv(os.path.join(rdir, 'mapping_affine_voronoi.csv'), index=False)
nl_nn.to_csv(os.path.join(rdir, 'mapping_nonlinear_nn.csv'), index=False)
nl_vor.to_csv(os.path.join(rdir, 'mapping_nonlinear_voronoi.csv'), index=False)
print("\nSaved 4 mapping CSVs")

# --- ATAC coordinates (both spaces) ---
atac_coords_df = pd.DataFrame({
    'spot_id': atac_names,
    'spot_idx': range(len(atac_names)),
    'affine_x': affine_um[:, 0],
    'affine_y': affine_um[:, 1],
    'nonlinear_x': warped_um[:, 0],
    'nonlinear_y': warped_um[:, 1],
})
for col in ['row', 'col', 'cluster', 'n_fragment', 'tsse']:
    if col in atac_nl.obs:
        atac_coords_df[col] = atac_nl.obs[col].values
atac_coords_df.to_csv(os.path.join(rdir, 'atac_coords.csv'), index=False)
print("Saved atac_coords.csv")

# --- Xenium coordinates ---
xenium_coords_df = pd.DataFrame({
    'cell_id': xenium_names,
    'cell_idx': range(len(xenium_names)),
    'x_um': xenium_um[:, 0],
    'y_um': xenium_um[:, 1],
})
xenium_coords_df.to_csv(os.path.join(rdir, 'xenium_coords.csv'), index=False)
print("Saved xenium_coords.csv")

# --- ATAC peak matrix (sparse, Matrix Market) ---
X_atac = atac_nl.X
if not issparse(X_atac):
    X_atac = csr_matrix(X_atac)
mmwrite(os.path.join(rdir, 'atac_peaks.mtx'), X_atac)
pd.DataFrame({'peak': atac_nl.var_names}).to_csv(
    os.path.join(rdir, 'atac_peak_names.csv'), index=False)
pd.DataFrame({'spot_id': atac_nl.obs_names}).to_csv(
    os.path.join(rdir, 'atac_spot_names.csv'), index=False)
print(f"Saved atac_peaks.mtx ({X_atac.shape})")

# --- Parse peak genomic coordinates ---
peak_info = []
for p in atac_nl.var_names:
    m = re.match(r'(chr\w+):(\d+)-(\d+)', p)
    if m:
        peak_info.append({
            'peak': p, 'chr': m.group(1),
            'start': int(m.group(2)), 'end': int(m.group(3))
        })
peak_df = pd.DataFrame(peak_info)
peak_df.to_csv(os.path.join(rdir, 'atac_peak_coords.csv'), index=False)
print(f"Saved atac_peak_coords.csv ({len(peak_df)} peaks)")

# --- Xenium expression matrix (sparse, Matrix Market) ---
X_xen = xenium.X
if not issparse(X_xen):
    X_xen = csr_matrix(X_xen)
mmwrite(os.path.join(rdir, 'xenium_expression.mtx'), X_xen)
pd.DataFrame({'gene': xenium.var_names}).to_csv(
    os.path.join(rdir, 'xenium_gene_names.csv'), index=False)
pd.DataFrame({'cell_id': xenium.obs_names}).to_csv(
    os.path.join(rdir, 'xenium_cell_names.csv'), index=False)
print(f"Saved xenium_expression.mtx ({X_xen.shape})")

# --- Summary ---
print(f"\n{'='*70}")
print(f"All exports in: {rdir}/")
for f in sorted(os.listdir(rdir)):
    sz = os.path.getsize(os.path.join(rdir, f))
    print(f"  {f}  ({sz/1024:.0f} KB)")
print(f"\nTransfer r_comparison/ to HPC for R analysis.")
