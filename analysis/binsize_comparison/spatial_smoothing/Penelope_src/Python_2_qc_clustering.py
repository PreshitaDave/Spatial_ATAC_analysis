"""
Python_2_qc_clustering.py — D02475 PanNET spatial ATAC-seq
===========================================================

Step 2: QC, spectral embedding, spatial embedding smoothing, Leiden clustering.

Clustering strategy (embedding-space smoothing):
  The spectral embedding is averaged over each spot's SPATIAL_K nearest spatial
  neighbors for SMOOTH_ITERS rounds before KNN/Leiden runs.  This makes the KNN
  graph spatially aware from the start, producing contiguous spatial domains by
  construction rather than correcting labels after the fact.

  Two Leiden resolutions are swept (0.5 and 1.0); leiden_sm_1.0 is the primary
  output used for peak calling and all downstream analysis.

Input:  D02475_out/step1_imported.h5ad      (from Python_1)
Output: D02475_out/step2_clustered.h5ad     (checkpoint for Python_3)
        D02475_out/cell_metadata.csv         (input for Signac_1)
        D02475_out/cluster_barcodes/         (per-cluster TXTs for Signac)
        D02475_out/tsse.png, umap_*.png, spatial_*.png
"""

import snapatac2 as snap
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import scipy.sparse as sp
import os
from resource_logger import ResourceLogger
rl = ResourceLogger("Python_2_qc_clustering")

sc.settings.verbosity = 3

# ============================================================
# Configuration
# ============================================================

OUTPUT_DIR = "/projectnb/rd-spat/HOME/pvarela/New_Pipeline/SnapATAC2_Signac_2/D02475_out"
SPATIAL_CSV = (
    "/projectnb/rd-spat/DATA/Lab/NET_project/Spatial_ATACseq/"
    "spatials/D2475/spatial/tissue_positions_list.csv"
)

# QC filter
MIN_TSSE = 4
# Composite artifact score: z(log n_fragment) - z(tsse).
# High scores flag spots with many fragments relative to TSS quality
# (contamination, tissue folds, barcode collisions).
# Set to None to skip this filter.
ARTIFACT_CUTOFF = 2.5

# Spectral embedding
# 5 kb bins: 10× larger than the 500 bp default, giving far less sparsity per
# feature for FFPE spatial ATAC where each tixel has ~500 usable nuclear frags.
N_FEATURES = 50000
N_COMPS    = 30
# Component 0 correlates with sequencing depth (r ≈ -0.99 on D02474, same chemistry).
# Drop it so KNN/UMAP/Leiden cluster on biology, not coverage.
USE_DIMS   = list(range(1, N_COMPS))   # dims 1..29, exclude depth component 0
# ENCODE hg38 blacklist: removes NUMT pseudo-gene bins and other artifact regions
# that inflate coverage and drive spurious clusters.
BLACKLIST  = "/projectnb/rd-spat/PROJECT/Collabs/panNET/Results/Jeff/20260608_Ruben_Grant/refs/hg38-blacklist.v2.bed"

# Leiden resolutions to sweep.  Both outputs are stored; leiden_sm_1.0 is the
# primary label used for peak calling and downstream analysis.
LEIDEN_RESOLUTIONS = [0.5, 1.0]

# Embedding-space spatial smoothing.
# Spectral vectors are averaged over SPATIAL_K nearest spatial neighbors
# for SMOOTH_ITERS rounds before KNN/Leiden runs.  This makes the Leiden
# graph itself spatially aware, producing contiguous domains by construction
# rather than correcting labels after the fact.
SPATIAL_K    = 8   # spatial neighbors for averaging
SMOOTH_ITERS = 2   # rounds of neighbor-averaging

SPOT_SIZE = 3   # marker size for spatial scatter plots

CLUSTER_COLORS = [
    "#C0392B", "#8E44AD", "#1A3A6B", "#E67E22", "#1ABC9C",
    "#2ECC71", "#D35400", "#7F8C8D", "#2980B9", "#F39C12",
    "#27AE60", "#E74C3C",
]

# ============================================================
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/cluster_barcodes", exist_ok=True)

# ============================================================
# Load checkpoint from Python_1
# ============================================================

print("Loading step1_imported.h5ad...")
data = snap.read(f"{OUTPUT_DIR}/step1_imported.h5ad", backed=None)
print(f"  Cells: {data.n_obs}")

# Re-read spatial CSV for coordinate lookup
spatial = pd.read_csv(
    SPATIAL_CSV,
    header=None,
    names=["barcode", "in_tissue", "array_row", "array_col", "pxl_col", "pxl_row"],
)
spatial = spatial[spatial["in_tissue"] == 1].copy()
spatial.index = spatial["barcode"]

# ============================================================
# STEP 2a — TSS enrichment score + filter
# ============================================================

print("Step 2a: TSS enrichment...")
snap.metrics.tsse(data, snap.genome.hg38)

fig = snap.pl.tsse(data, interactive=False, show=False)
fig.update_layout(
    width=600, height=500,
    plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)",
    xaxis=dict(
        showline=True, showgrid=False, linewidth=1, linecolor="black",
        ticks="outside", tickwidth=2, tickcolor="black", ticklen=10, mirror=True,
        title_font=dict(family="Arial", size=18, color="black"),
        tickfont=dict(family="Arial", size=16, color="black"),
    ),
    yaxis=dict(
        showline=True, showgrid=False, linewidth=1, linecolor="black",
        ticks="outside", tickwidth=2, tickcolor="black", ticklen=10, mirror=True,
        title_font=dict(family="Arial", size=18, color="black"),
        tickfont=dict(family="Arial", size=16, color="black"),
    ),
)
fig.write_image(f"{OUTPUT_DIR}/tsse.png", scale=4)

snap.pp.filter_cells(data, min_tsse=MIN_TSSE)
print(f"  Cells after TSSE >= {MIN_TSSE} filter: {data.n_obs}")
rl.step("TSS enrichment + filter")

# ============================================================
# STEP 2a2 — Composite artifact score filter
#
# artifact_score = z(log1p n_fragment) - z(tsse)
# High = many fragments relative to TSS quality.
# This catches contamination, tissue folds, and barcode collisions
# without discarding clean high-coverage spots that have good TSS.
# ============================================================

if ARTIFACT_CUTOFF is not None:
    from scipy.stats import zscore as _zscore

    log_frag = np.log1p(data.obs["n_fragment"].values.astype(float))
    z_frag   = _zscore(log_frag)
    z_tsse   = _zscore(data.obs["tsse"].values.astype(float))
    artifact = z_frag - z_tsse
    data.obs["artifact_score"] = artifact

    # Spatial plot of artifact score before filtering so the cutoff can be tuned.
    _coords = spatial.reindex(data.obs_names)
    _art_fig = go.Figure(go.Scatter(
        x=_coords["pxl_col"].tolist(),
        y=_coords["pxl_row"].tolist(),
        mode="markers",
        marker=dict(
            color=artifact.tolist(),
            colorscale="RdYlBu_r",
            size=SPOT_SIZE,
            opacity=0.9,
            showscale=True,
            colorbar=dict(title="artifact_score"),
        ),
    ))
    _art_fig.update_layout(
        title="Artifact score (z_frag − z_tsse) — pre-filter",
        yaxis=dict(scaleanchor="x"),
        plot_bgcolor="white", paper_bgcolor="white",
    )
    _art_fig.write_image(f"{OUTPUT_DIR}/artifact_score_spatial.png", scale=4)
    print(f"  Saved: artifact_score_spatial.png")

    n_before = data.n_obs
    keep     = data.obs["artifact_score"] < ARTIFACT_CUTOFF
    data     = data[keep].copy()
    print(f"  Artifact filter (score < {ARTIFACT_CUTOFF}): "
          f"removed {n_before - data.n_obs} spots, {data.n_obs} remain")
    rl.step("artifact score filter")

# ============================================================
# STEP 2b — Tile matrix, feature selection, spectral embedding, UMAP
# No spatial signal added here — pure chromatin embedding.
# ============================================================

print("Step 2b: tile matrix (5 kb bins) + spectral embedding...")
snap.pp.add_tile_matrix(data, bin_size=5000)
snap.pp.select_features(data, n_features=N_FEATURES, blacklist=BLACKLIST, inplace=True)

snap.tl.spectral(data, n_comps=N_COMPS)

# Check whether component 0 is a depth artifact before dropping it.
import scipy.stats as _stats
_comp0 = data.obsm["X_spectral"][:, 0]
_r, _p = _stats.pearsonr(_comp0, data.obs["n_fragment"])
print(f"Spectral component 0 vs n_fragment: r = {_r:.3f}, p = {_p:.2e}")
print(f"  {'Dropping component 0 (|r| > 0.5).' if abs(_r) > 0.5 else 'Component 0 does NOT appear to be a depth artifact (|r| <= 0.5) — consider keeping it.'}")
rl.step("tile matrix + spectral embedding")

# QC UMAP on raw spectral embedding (stored as X_umap).
# Used only for n_fragment / tsse diagnostic plots.
# The clustering UMAP uses the spatially-smoothed embedding built in Step 2c.
snap.tl.umap(data, use_dims=USE_DIMS, min_dist=0.05)

sc.settings.set_figure_params(dpi=300, facecolor="white", fontsize=12)

fig = snap.pl.umap(data, color="n_fragment", show=False)
fig.write_image(f"{OUTPUT_DIR}/umap_nfrag.png", scale=4)
fig = snap.pl.umap(data, color="tsse", show=False)
fig.write_image(f"{OUTPUT_DIR}/umap_tsse.png", scale=4)

# ============================================================
# STEP 2c — Embedding-space spatial smoothing + clustering
#
# Rather than correcting labels after clustering, we smooth the
# spectral embedding itself over each spot's SPATIAL_K nearest
# spatial neighbors for SMOOTH_ITERS rounds.  The resulting
# X_spectral_smooth is then used for KNN / UMAP / Leiden so that
# the graph is spatially informed from the start.
# ============================================================

print(f"Step 2c: spatial embedding smoothing (k={SPATIAL_K}, iters={SMOOTH_ITERS})...")

aligned = spatial.reindex(data.obs_names)
coords_xy = np.column_stack([
    aligned["pxl_col"].values.astype(float),
    aligned["pxl_row"].values.astype(float),
])

# Build row-normalised spatial averaging operator (n × n sparse)
tree = cKDTree(coords_xy)
_, idx = tree.query(coords_xy, k=SPATIAL_K + 1)   # col 0 = self
n = coords_xy.shape[0]
row_idx = np.repeat(np.arange(n), SPATIAL_K + 1)
col_idx = idx.ravel()
A = sp.csr_matrix((np.ones(row_idx.size), (row_idx, col_idx)), shape=(n, n))
A = A.multiply(1.0 / A.sum(1))   # row-normalised
A = sp.csr_matrix(A)

# Extract depth-free spectral embedding (drop comp 0) and smooth iteratively
Z  = np.asarray(data.obsm["X_spectral"])[:, USE_DIMS]
Zs = Z.copy()
for _ in range(SMOOTH_ITERS):
    Zs = A @ Zs
data.obsm["X_spectral_smooth"] = Zs
print(f"  Smoothed embedding stored as X_spectral_smooth {Zs.shape}")
rl.step("spatial embedding smoothing")

print("  KNN + UMAP + Leiden on smoothed embedding...")
snap.pp.knn(data, use_rep="X_spectral_smooth")
snap.tl.umap(data, use_rep="X_spectral_smooth", key_added="umap_smooth", min_dist=0.05)

cluster_id_sets = {}
for r in LEIDEN_RESOLUTIONS:
    key = f"leiden_sm_{r}"
    snap.tl.leiden(data, resolution=r, key_added=key)
    labels = np.array([str(x) for x in list(data.obs[key])])
    cluster_id_sets[key] = (labels, sorted(set(labels.tolist()), key=lambda x: int(x)))
    print(f"  {key}: {len(cluster_id_sets[key][1])} clusters  {cluster_id_sets[key][1]}")
rl.step("KNN + UMAP + Leiden clustering")

# Primary output: leiden_sm_1.0 (finer resolution, spatially smoothed).
# leiden_sm_0.5 kept as the coarse-resolution companion.
primary_key = f"leiden_sm_{LEIDEN_RESOLUTIONS[-1]}"
primary_labels, cluster_ids_smooth = cluster_id_sets[primary_key]
coarse_labels,  cluster_ids_raw    = cluster_id_sets[f"leiden_sm_{LEIDEN_RESOLUTIONS[0]}"]

# Build obs DataFrame for plots + export.
# "leiden" and "leiden_smooth" are compat aliases read by Python_3/Python_4/Signac.
obs = pd.DataFrame(
    {"leiden_sm_0.5": coarse_labels,
     "leiden_sm_1.0": primary_labels,
     "leiden":        coarse_labels,    # compat alias: res 0.5
     "leiden_smooth": primary_labels,   # compat alias: res 1.0 (primary)
     "tsse":          list(data.obs["tsse"]),
     "n_fragment":    list(data.obs["n_fragment"])},
    index=data.obs_names,
)
obs["array_row"]  = aligned["array_row"].values
obs["array_col"]  = aligned["array_col"].values
obs["spatial_x"]  = aligned["pxl_col"].values
obs["spatial_y"]  = aligned["pxl_row"].values
obs_filt = obs.copy()   # all on-tissue spots retained; no per-cluster outlier removal

# UMAP of primary clusters (smoothed embedding)
umap_sm = np.asarray(data.obsm["X_umap_smooth"])
fig_um, ax_um = plt.subplots(figsize=(8, 6))
for i, cid in enumerate(cluster_ids_smooth):
    mask = primary_labels == cid
    ax_um.scatter(umap_sm[mask, 0], umap_sm[mask, 1],
                  c=[CLUSTER_COLORS[i % len(CLUSTER_COLORS)]], s=3, label=cid)
ax_um.legend(title="cluster", bbox_to_anchor=(1.05, 1), loc="upper left", markerscale=3)
ax_um.set_xlabel("UMAP1"); ax_um.set_ylabel("UMAP2")
ax_um.set_title(f"Smoothed embedding — leiden_sm_1.0 ({len(cluster_ids_smooth)} clusters)")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/umap_leiden_smooth.png", dpi=200, bbox_inches="tight")
plt.close()
print(f"  Saved: umap_leiden_smooth.png")

# ============================================================
# STEP 2f — Spatial scatter plots
# ============================================================

print("Step 2f: spatial plots...")


def _categorical_spatial(obs_df, label_col, cluster_list, title, out_file):
    fig = go.Figure()
    for i, cid in enumerate(cluster_list):
        mask = obs_df[label_col] == cid
        fig.add_trace(go.Scatter(
            x=obs_df.loc[mask, "spatial_x"].tolist(),
            y=obs_df.loc[mask, "spatial_y"].tolist(),
            mode="markers", name=f"Cluster {cid}",
            marker=dict(color=CLUSTER_COLORS[i % len(CLUSTER_COLORS)],
                        size=SPOT_SIZE, opacity=0.8),
        ))
    fig.update_layout(title=title, yaxis=dict(scaleanchor="x"),
                      xaxis_title="pixel_x", yaxis_title="pixel_y",
                      plot_bgcolor="white", paper_bgcolor="white")
    fig.write_image(out_file, scale=4)


def _continuous_spatial(obs_df, value_col, title, out_file, colorscale="Viridis", vmax=None):
    fig = go.Figure(go.Scatter(
        x=obs_df["spatial_x"].tolist(), y=obs_df["spatial_y"].tolist(),
        mode="markers",
        marker=dict(color=obs_df[value_col].tolist(), colorscale=colorscale,
                    size=SPOT_SIZE, opacity=0.8, cmax=vmax, showscale=True,
                    colorbar=dict(title=value_col)),
    ))
    fig.update_layout(title=title, yaxis=dict(scaleanchor="x"),
                      plot_bgcolor="white", paper_bgcolor="white")
    fig.write_image(out_file, scale=4)


_categorical_spatial(obs_filt, "leiden_sm_0.5", cluster_ids_raw,
                     "Spatially smoothed clusters (res=0.5)",
                     f"{OUTPUT_DIR}/spatial_leiden_sm0.5.png")
_categorical_spatial(obs_filt, "leiden_sm_1.0", cluster_ids_smooth,
                     "Spatially smoothed clusters (res=1.0, primary output)",
                     f"{OUTPUT_DIR}/spatial_leiden_smooth.png")
_continuous_spatial(obs_filt, "n_fragment", "n_fragment",
                    f"{OUTPUT_DIR}/spatial_nfrag.png", vmax=10000)
_continuous_spatial(obs_filt, "tsse", "TSSE",
                    f"{OUTPUT_DIR}/spatial_tsse.png")

for i, cid in enumerate(cluster_ids_smooth):
    colors = [CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
              if c == cid else "gainsboro"
              for c in obs_filt["leiden_smooth"]]
    fig = go.Figure(go.Scatter(
        x=obs_filt["spatial_x"].tolist(), y=obs_filt["spatial_y"].tolist(),
        mode="markers", marker=dict(color=colors, size=SPOT_SIZE, opacity=0.8),
    ))
    fig.update_layout(title=f"Cluster #{cid} (smoothed)",
                      yaxis=dict(scaleanchor="x"),
                      plot_bgcolor="white", paper_bgcolor="white")
    fig.write_image(f"{OUTPUT_DIR}/spatial_cluster_{cid}.png", scale=2)

# ============================================================
# Export metadata CSV + per-cluster barcode lists for Signac
# ============================================================

print("Exporting metadata CSV and per-cluster barcode lists...")

meta_out = obs.copy()
meta_out.index.name = "cellid"
meta_out.to_csv(f"{OUTPUT_DIR}/cell_metadata.csv")

barcode_dir = f"{OUTPUT_DIR}/cluster_barcodes"
for cid in cluster_ids_smooth:
    mask     = obs_filt["leiden_smooth"] == cid
    barcodes = obs_filt.index[mask].tolist()
    out_path = f"{barcode_dir}/cluster_{cid}_barcodes.txt"
    with open(out_path, "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    print(f"  Cluster {cid}: {len(barcodes)} barcodes → cluster_{cid}_barcodes.txt")

# ============================================================
# Save checkpoint for Python_3_gene_activity.py
# ============================================================

data.write(f"{OUTPUT_DIR}/step2_clustered.h5ad")
print(f"\nCheckpoint saved: {OUTPUT_DIR}/step2_clustered.h5ad")
rl.step("checkpoint write")
rl.done()
print(f"Metadata CSV:      {OUTPUT_DIR}/cell_metadata.csv")
print("Run Signac_1_processing.R and Python_3_gene_activity.py next.")
