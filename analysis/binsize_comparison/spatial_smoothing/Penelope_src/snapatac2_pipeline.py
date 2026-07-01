#!/usr/bin/env python

"""
SnapATAC2 spatial ATAC-seq pipeline — PanNET sample D02475  [FIXED v3]
=======================================================================

Architecture change from v2
----------------------------
v2 appended spatial coordinates to the spectral embedding before KNN/Leiden.
This caused Leiden to tile the tissue geographically rather than find chromatin
communities — producing spatially coherent but biologically overclustered maps.

v3 uses a TWO-STAGE approach:

  Stage 1 — Chromatin-only clustering (Steps 1–4, unchanged from v2)
    KNN and Leiden run on the pure spectral embedding with NO spatial signal.
    This produces clusters that reflect chromatin state, not geography.
    Resolution is deliberately conservative (LEIDEN_RES = 0.3) to get broad,
    biologically meaningful groupings for a PanNET (~5–8 clusters expected).

  Stage 2 — Spatial majority-vote smoothing (Step 4b, new)
    For each spot, look at its K nearest spatial neighbors (array grid neighbors,
    not chromatin neighbors). Assign the spot the majority cluster label among
    its spatial neighborhood. This smooths away isolated spots that are
    chromatinically assigned to a distant cluster while surrounded by a different
    tissue region — giving spatial coherence WITHOUT forcing Leiden to tile.

    SMOOTH_K controls the neighborhood radius:
      - SMOOTH_K = 0  → no smoothing (pure chromatin clusters, may be noisy)
      - SMOOTH_K = 15 → mild smoothing, removes isolated spots
      - SMOOTH_K = 30 → moderate smoothing, recommended starting point
      - SMOOTH_K = 50 → strong smoothing, may over-merge boundaries
    The smoothed label is stored as obs["leiden_smooth"]; the original chromatin
    label is kept as obs["leiden"] for comparison.

Other fixes retained from v2:
  - Doublet filter (3× median fragment count)
  - Per-cluster barcode export for Signac pseudo-bulk peak calling
"""

import snapatac2 as snap
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import cKDTree
from collections import Counter
import subprocess
import os

sc.settings.verbosity = 3

# =============================================================================
# Configuration
# =============================================================================

FRAGMENT_FILE = "/projectnb/rd-spat/DATA/Lab/NET_project/Spatial_ATACseq/chromap_outputs/D02475_NG08213/chromap_output/fragments.tsv.gz"
SPATIAL_CSV   = "/projectnb/rd-spat/DATA/Lab/NET_project/Spatial_ATACseq/spatials/D2475/spatial/tissue_positions_list.csv"
SAMPLE_ID     = "D02475"
OUTPUT_DIR    = "./snapatac2_out"
CACHE_DIR     = os.path.join(OUTPUT_DIR, "cache")

MIN_FRAGS             = 1000
MIN_TSSE              = 4
N_FEATURES            = 50000
N_COMPS               = 30

# Stage 1: chromatin-only Leiden.
# For a PanNET slide you likely have 5–8 biologically meaningful compartments.
# Start at 0.5; raise toward 0.8 if you want more sub-structure.
LEIDEN_RES = 1.0

# Per-cluster fragment outlier filter (Step 4c)
# Spots with n_fragment > N × their cluster median are removed.
# 3.0 is conservative; lower to 2.5 if cluster 3 still dominates DA peaks.
PER_CLUSTER_MAX_MULTIPLIER = 3.0

# Stage 2: spatial smoothing.
# Number of nearest spatial neighbors used for majority-vote label smoothing.
# Uses pixel coordinates (same scale as pxl_col / pxl_row in tissue_positions).
# 30 is a good starting point for Visium (55 µm spot pitch).
SMOOTH_K = 15
SPOT_SIZE = 3      # marker size for spatial scatter plots — decrease for denser slides

# Marker genes for the dot-plot in Step 6b.
# These are PanNET-relevant genes; edit to match your biology.
# Genes absent from the gene activity matrix are silently skipped.
MARKER_GENES = [
    # neuroendocrine
    "CHGA", "CHGB", "SYP", "ENO2",
    # islet cell lineage
    "INS", "GCG", "SST", "PPY",
    # transcription factors
    "ARX", "PDX1", "NKX2-2", "NEUROD1", "PAX6",
    # PanNET tumor suppressors
    "MEN1", "DAXX", "ATRX",
    # receptor / signalling
    "SSTR2",
    # proliferation
    "MKI67",
]

CLUSTER_COLORS = [
        "#C0392B", "#8E44AD", "#1A3A6B",
        "#E67E22", "#1ABC9C", "#2ECC71",
        "#D35400", "#7F8C8D", "#2980B9",
        "#F39C12", "#27AE60", "#E74C3C"
]

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR,  exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "cluster_barcodes"), exist_ok=True)

CP = {
    1: os.path.join(CACHE_DIR, "step1_imported.h5ad"),
    2: os.path.join(CACHE_DIR, "step2_qc.h5ad"),
    3: os.path.join(CACHE_DIR, "step3_embedded.h5ad"),
    4: os.path.join(CACHE_DIR, "step4_clustered.h5ad"),
}

# Cache invalidation: steps 3 and 4 depend on N_COMPS and LEIDEN_RES.
# If you change either parameter, the old checkpoints are stale and must be
# deleted so the pipeline recomputes the embedding and clustering.
# This block writes a small params file next to the cache and clears steps
# 3 and 4 automatically if the params have changed since last run.
_params_file = os.path.join(CACHE_DIR, "embedding_params.txt")
_current_params = f"N_COMPS={N_COMPS},LEIDEN_RES={LEIDEN_RES}"
_cached_params  = open(_params_file).read().strip() if os.path.exists(_params_file) else ""

if _cached_params and _cached_params != _current_params:
    print(f"  Parameter change detected: {_cached_params} → {_current_params}")
    print("  Deleting stale embedding/clustering checkpoints (steps 3 and 4)...")
    for _step in [3, 4]:
        if os.path.exists(CP[_step]):
            os.remove(CP[_step])
            print(f"    Deleted: {os.path.basename(CP[_step])}")

with open(_params_file, "w") as _f:
    _f.write(_current_params)


def _load(step):
    path = CP[step]
    print(f"  Resuming from checkpoint: {os.path.basename(path)}")
    return snap.read(path, backed=None)


def _save(data, step):
    data.write(CP[step])
    print(f"  Checkpoint saved: {os.path.basename(CP[step])}")


def _obs_to_df(data, cols):
    return pd.DataFrame(
        {c: list(data.obs[c]) for c in cols},
        index=data.obs_names,
    )


# =============================================================================
# Read spatial positions (always needed)
# =============================================================================

spatial = pd.read_csv(
    SPATIAL_CSV,
    header=None,
    names=["barcode", "in_tissue", "array_row", "array_col", "pxl_col", "pxl_row"],
)
spatial = spatial[spatial["in_tissue"] == 1].copy()
spatial.index = spatial["barcode"] + "-1"
in_tissue_barcodes = spatial.index.tolist()

# =============================================================================
# STEP 1 — Import fragments
# =============================================================================

if os.path.exists(CP[1]):
    data = _load(1)
else:
    print("Step 1: importing fragments...")
    data = snap.pp.import_fragments(
        FRAGMENT_FILE,
        chrom_sizes=snap.genome.hg38,
        sorted_by_barcode=False,
        whitelist=in_tissue_barcodes,
        min_num_fragments=MIN_FRAGS,
    )
    print(f"  Cells after import: {data.n_obs}")
    _save(data, 1)

# =============================================================================
# STEP 2 — QC: TSS enrichment filter 
# =============================================================================

if os.path.exists(CP[2]):
    data = _load(2)
else:
    print("Step 2: QC metrics + filters...")
    fig = snap.pl.frag_size_distr(data, show=False)
    fig.write_image(f"{OUTPUT_DIR}/fragment_size_dist.png", scale=4)

    snap.metrics.tsse(data, snap.genome.hg38)

    fig = snap.pl.tsse(data, interactive=False, show=False)
    fig.write_image(f"{OUTPUT_DIR}/tsse.png", scale=4)

    snap.pp.filter_cells(data, min_tsse=MIN_TSSE)
    print(f"  Cells after TSSE >= {MIN_TSSE} filter: {data.n_obs}")

   

    _save(data, 2)

# =============================================================================
# STEP 3 — Tile matrix, PURE CHROMATIN spectral embedding, UMAP
# No spatial signal injected here (v3 change).
# =============================================================================

if os.path.exists(CP[3]):
    data = _load(3)
else:
    print("Step 3: tile matrix, spectral embedding (chromatin only), UMAP...")
    snap.pp.add_tile_matrix(data)
    snap.pp.select_features(data, n_features=N_FEATURES, inplace=True)

    # Pure chromatin spectral embedding — spatial coords NOT appended here
    snap.tl.spectral(data, n_comps=N_COMPS)

    snap.tl.umap(data, min_dist=0.05)

    fig = snap.pl.umap(data, color="n_fragment", show=False)
    fig.write_image(f"{OUTPUT_DIR}/umap_nfrag.png", scale=4)
    fig = snap.pl.umap(data, color="tsse", show=False)
    fig.write_image(f"{OUTPUT_DIR}/umap_tsse.png", scale=4)

    _save(data, 3)

# =============================================================================
# STEP 4 — KNN + Leiden on chromatin embedding only (Stage 1)
# =============================================================================

if os.path.exists(CP[4]):
    data = _load(4)
else:
    print("Step 4: KNN + Leiden (chromatin only, Stage 1)...")
    snap.pp.knn(data)
    snap.tl.leiden(data, resolution=LEIDEN_RES, key_added="leiden")

    fig = snap.pl.umap(data, color="leiden", show=False)
    fig.write_image(f"{OUTPUT_DIR}/umap_leiden_chromatin.png", scale=4)

    _save(data, 4)

# =============================================================================
# STEP 4b — Spatial majority-vote smoothing (Stage 2, new in v3)
#
# Goal: remove isolated spots that are chromatinically assigned to a distant
# cluster while physically surrounded by a different tissue region.
#
# Method:
#   1. Build a KD-tree on pixel coordinates of all retained barcodes.
#   2. For each spot, find its SMOOTH_K nearest spatial neighbors.
#   3. Assign the majority cluster label among those neighbors.
#      Ties are broken by keeping the original label.
#   4. Store as obs["leiden_smooth"]; keep obs["leiden"] for comparison.
#
# This step does NOT use the spectral embedding at all — it only uses the
# tissue_positions pixel coordinates.  It runs fast (seconds) and does not
# require re-running the checkpoint steps above.
# =============================================================================

print("Step 4b: spatial majority-vote label smoothing...")

obs = _obs_to_df(data, ["leiden", "tsse", "n_fragment"])
aligned = spatial.reindex(data.obs_names)
obs["spatial_x"] = aligned["pxl_col"].values
obs["spatial_y"] = aligned["pxl_row"].values

# Drop spots with missing spatial coords (shouldn't happen after whitelist filter)
valid_mask = ~(np.isnan(obs["spatial_x"]) | np.isnan(obs["spatial_y"]))
if (~valid_mask).sum() > 0:
    print(f"  Warning: {(~valid_mask).sum()} spots with missing spatial coords — kept with original label")

coords = obs[["spatial_x", "spatial_y"]].values
tree   = cKDTree(coords)

leiden_labels = obs["leiden"].values.copy()
smoothed      = leiden_labels.copy()

# Query SMOOTH_K+1 neighbors (first result is the spot itself)
_, neighbor_idx = tree.query(coords, k=SMOOTH_K + 1)

for i in range(len(leiden_labels)):
    if not valid_mask.iloc[i]:
        continue
    neighbor_labels = leiden_labels[neighbor_idx[i, 1:]]   # exclude self
    vote, count     = Counter(neighbor_labels).most_common(1)[0]
    smoothed[i]     = vote

obs["leiden_smooth"] = smoothed

# How many spots changed label?
n_changed = (smoothed != leiden_labels).sum()
print(f"  Spots relabeled by smoothing: {n_changed} / {len(leiden_labels)} "
      f"({100 * n_changed / len(leiden_labels):.1f}%)")

# Write smoothed labels back to data.obs for downstream steps
# SnapATAC2 obs is polars-backed; safest to keep obs separately and use it
# for all downstream steps rather than writing back.
# (We pass obs["leiden_smooth"] to Signac export below.)

cluster_ids_raw    = sorted(set(leiden_labels.tolist()),  key=lambda x: int(x))
cluster_ids_smooth = sorted(set(smoothed.tolist()),       key=lambda x: int(x))
print(f"  Clusters (chromatin): {len(cluster_ids_raw)}  |  "
      f"Clusters (smoothed): {len(cluster_ids_smooth)}")

# =============================================================================
# STEP 4c — Per-cluster fragment outlier removal
#
# Why per-cluster rather than global?
#   A global upper fragment cutoff (e.g. 3× overall median) removes spots from
#   cellularly dense tissue regions (e.g. tumor core) that genuinely have more
#   fragments than stromal spots — those aren't doublets, they're biology.
#   Per-cluster filtering is more principled: it removes spots that are outliers
#   *within their own cluster*, which are the true technical artifacts (doublets,
#   ambient DNA) regardless of which tissue compartment they come from.
#
# Effect on Signac:
#   Cluster 3 in the previous run had median fragments ~5–10× higher than other
#   clusters, with spots reaching 125k fragments vs ~5–10k elsewhere. This caused
#   DESeq2 to call ~7k spurious DA peaks for cluster 3 and saturated all motif
#   p-values at machine epsilon. Per-cluster filtering removes those extreme spots
#   while keeping the legitimate high-coverage spots that define the cluster.
#
# PER_CLUSTER_MAX_MULTIPLIER: spots with n_fragment > N × cluster_median are
# removed. 3.0 is conservative; lower toward 2.5 if cluster 3 still dominates
# after this step, or raise toward 4.0 to be more permissive.
# =============================================================================
 
print("Step 4c: per-cluster fragment outlier removal...")
 
n_frags       = np.array(list(data.obs["n_fragment"]))
# Align fragment counts to the obs dataframe index (obs was built from data.obs_names)
frag_series   = pd.Series(n_frags, index=data.obs_names)
frag_aligned  = frag_series.reindex(obs.index).values
 
keep_mask = np.ones(len(obs), dtype=bool)
 
for cid in cluster_ids_smooth:
    cl_mask   = obs["leiden_smooth"].values == cid
    cl_frags  = frag_aligned[cl_mask]
    cl_median = np.median(cl_frags)
    cl_cutoff = PER_CLUSTER_MAX_MULTIPLIER * cl_median
    outliers  = cl_mask & (frag_aligned > cl_cutoff)
    keep_mask[outliers] = False
    n_out = outliers.sum()
    print(f"  Cluster {cid}: median={cl_median:.0f}  cutoff={cl_cutoff:.0f}  "
          f"outliers removed={n_out} ({100*n_out/cl_mask.sum():.1f}%)")
 
n_before = len(obs)
obs      = obs[keep_mask].copy()
print(f"  Total spots retained: {len(obs)} / {n_before} "
      f"({100*len(obs)/n_before:.1f}%)")
 
# Recompute cluster list after filtering (a cluster could theoretically vanish)
cluster_ids_smooth = sorted(set(obs["leiden_smooth"].tolist()), key=lambda x: int(x))
print(f"  Clusters after outlier removal: {len(cluster_ids_smooth)}")

# =============================================================================
# STEP 5 — Spatial plots: both raw chromatin and smoothed labels
# =============================================================================

print("Step 5: spatial plots...")


def _spatial_scatter(x, y, color_vals, color_label, title="", out_file=None,
                     scale=4, colorscale="Viridis", vmax=None):
    fig = go.Figure(go.Scatter(
        x=x, y=y, mode="markers",
        marker=dict(color=color_vals, colorscale=colorscale, size=SPOT_SIZE, opacity=0.8,
                    cmax=vmax, showscale=True, colorbar=dict(title=color_label)),
    ))
    fig.update_layout(title=title, yaxis=dict(scaleanchor="x"),
                      xaxis_title="spatial_x", yaxis_title="spatial_y",
                      plot_bgcolor="white", paper_bgcolor="white")
    if out_file:
        fig.write_image(out_file, scale=scale)
    return fig


def _categorical_spatial(obs_df, label_col, cluster_list, title, out_file, scale=4):
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
                      xaxis_title="spatial_x", yaxis_title="spatial_y",
                      plot_bgcolor="white", paper_bgcolor="white")
    fig.write_image(out_file, scale=scale)


# Chromatin-only clusters (for comparison)
_categorical_spatial(obs, "leiden", cluster_ids_raw,
                     "Spatial clusters — chromatin only (pre-smoothing)",
                     f"{OUTPUT_DIR}/spatial_leiden_chromatin.png")

# Smoothed clusters (primary output)
_categorical_spatial(obs, "leiden_smooth", cluster_ids_smooth,
                     f"Spatial clusters — smoothed (K={SMOOTH_K})",
                     f"{OUTPUT_DIR}/spatial_leiden_smooth.png")

_spatial_scatter(obs["spatial_x"], obs["spatial_y"], obs["n_fragment"],
                 "n_fragment", vmax=10000, out_file=f"{OUTPUT_DIR}/spatial_nfrag.png")
_spatial_scatter(obs["spatial_x"], obs["spatial_y"], obs["tsse"],
                 "TSSE", out_file=f"{OUTPUT_DIR}/spatial_tsse.png")

# Per-cluster highlight plots using smoothed labels
for i, cid in enumerate(cluster_ids_smooth):
    colors = [CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
              if c == cid else "gainsboro"
              for c in obs["leiden_smooth"]]
    fig = go.Figure(go.Scatter(
        x=obs["spatial_x"].tolist(), y=obs["spatial_y"].tolist(),
        mode="markers", marker=dict(color=colors, size=SPOT_SIZE, opacity=0.8),
    ))
    fig.update_layout(title=f"Cluster #{cid} (smoothed)",
                      yaxis=dict(scaleanchor="x"),
                      plot_bgcolor="white", paper_bgcolor="white")
    fig.write_image(f"{OUTPUT_DIR}/spatial_cluster_{cid}.png", scale=2)

# =============================================================================
# STEP 6 — Gene activity matrix and differential analysis
# Uses smoothed cluster labels.
# =============================================================================

GENE_MATRIX_CACHE = os.path.join(CACHE_DIR, "step6_gene_matrix.h5ad")

if os.path.exists(GENE_MATRIX_CACHE):
    print("Step 6: loading gene activity from cache...")
    gene_matrix = sc.read(GENE_MATRIX_CACHE)
else:
    print("Step 6: gene activity matrix...")
    gene_matrix = snap.pp.make_gene_matrix(data, snap.genome.hg38)

    sc.pp.filter_genes(gene_matrix, min_cells=20)

    # Drop MT genes before normalize_total. Without removal, normalize_total
    # sets the scale using mostly MT counts (D02475 is ~60-72% mitochondrial),
    # compressing all nuclear marker scores toward zero and biasing MAGIC
    # toward MT-driven cell similarity.
    is_mt = gene_matrix.var_names.str.startswith(("MT-", "MT.", "mt-"))
    n_mt = is_mt.sum()
    gene_matrix = gene_matrix[:, ~is_mt].copy()
    print(f"  Dropped {n_mt} mitochondrial genes before normalization")

    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)

    # MAGIC imputation — reduces dropout in sparse spatial ATAC gene activity
    # scores before differential testing.  Requires: pip install magic-impute
    try:
        sc.external.pp.magic(gene_matrix, solver="approximate")
        print("  MAGIC imputation complete")
    except Exception as e:
        print(f"  MAGIC skipped ({e}); install with: pip install magic-impute")

    sc.pp.highly_variable_genes(gene_matrix, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # Use smoothed labels; drop cells removed during per-cluster outlier filter
    gene_matrix.obs["leiden_smooth"] = obs.reindex(gene_matrix.obs_names)["leiden_smooth"].values
    gene_matrix = gene_matrix[gene_matrix.obs["leiden_smooth"].notna()].copy()

    sc.tl.rank_genes_groups(
        gene_matrix,
        groupby="leiden_smooth",
        method="wilcoxon",
        key_added="rank_genes_groups_leiden",
    )

    rgg = gene_matrix.uns["rank_genes_groups_leiden"]
    for key in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        pd.DataFrame(rgg[key]).to_csv(
            f"{OUTPUT_DIR}/rank_genes_groups_leiden.{key}.csv"
        )

    gene_matrix.write(GENE_MATRIX_CACHE)
    print(f"  Gene activity AnnData cached: {os.path.basename(GENE_MATRIX_CACHE)}")

# =============================================================================
# STEP 6b — Dot plot for marker gene activity
# =============================================================================

print("Step 6b: dot plot for marker gene activity...")
available_markers = [g for g in MARKER_GENES if g in gene_matrix.var_names]

if not available_markers:
    print("  No MARKER_GENES found in gene matrix — skipping dot plot.")
    print("  Edit MARKER_GENES in the configuration section to match your genes.")
else:
    if len(available_markers) < len(MARKER_GENES):
        missing = [g for g in MARKER_GENES if g not in gene_matrix.var_names]
        print(f"  {len(available_markers)}/{len(MARKER_GENES)} markers found; "
              f"missing from matrix: {missing}")

    cmap_dot = LinearSegmentedColormap.from_list(
        "dot_cmap", ["w", "lightblue", "orangered"], N=256
    )
    sc.settings.set_figure_params(dpi=200, facecolor="white", fontsize=12)
    sc.pl.dotplot(
        gene_matrix, var_names=available_markers, groupby="leiden_smooth",
        standard_scale="var", expression_cutoff=0.85, swap_axes=True,
        cmap=cmap_dot, show=False,
    )
    plt.savefig(f"{OUTPUT_DIR}/dotplot_marker_genes.png", dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Dot plot saved: dotplot_marker_genes.png")

# =============================================================================
# STEP 6c — Cluster correlation matrix on gene activity scores
# =============================================================================

print("Step 6c: cluster correlation matrix...")
sc.settings.set_figure_params(dpi=200, facecolor=(0, 0, 0, 0), fontsize=12)
sc.pl.correlation_matrix(
    gene_matrix,
    groupby="leiden_smooth",
    dendrogram=True,
    show_correlation_numbers=True,
    show=False,
)
plt.savefig(f"{OUTPUT_DIR}/cluster_correlation_matrix.png", dpi=200, bbox_inches="tight")
plt.close("all")
print("  Correlation matrix saved: cluster_correlation_matrix.png")

# =============================================================================
# STEP 7 — Export metadata + per-cluster barcodes (smoothed labels) for Signac
# =============================================================================

print("Step 7: exporting metadata and per-cluster barcodes for Signac...")

meta_out = _obs_to_df(data, ["leiden", "tsse", "n_fragment"])
meta_out["leiden_smooth"] = obs.reindex(meta_out.index)["leiden_smooth"].values
aligned = spatial.reindex(data.obs_names)
meta_out["array_row"] = aligned["array_row"].values
meta_out["array_col"] = aligned["array_col"].values
meta_out["spatial_x"] = aligned["pxl_col"].values
meta_out["spatial_y"] = aligned["pxl_row"].values
meta_out.index.name = "cellid"
meta_out.to_csv(f"{OUTPUT_DIR}/cell_metadata.csv")

# Write per-cluster barcode files using SMOOTHED labels
barcode_dir = os.path.join(OUTPUT_DIR, "cluster_barcodes")
for cid in cluster_ids_smooth:
    mask     = meta_out["leiden_smooth"] == cid
    barcodes = meta_out.index[mask].tolist()
    out_path = os.path.join(barcode_dir, f"cluster_{cid}_barcodes.txt")
    with open(out_path, "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    print(f"  Cluster {cid}: {len(barcodes)} barcodes → {os.path.basename(out_path)}")

# =============================================================================
# STEP 8 — Cluster-specific fragment BED files for genome browser visualization
#
# Splits the bulk fragment file into one sorted BED per cluster, then bgzips
# and tabix-indexes each for direct loading in IGV or UCSC.  Uses the
# per-cluster outlier-filtered obs barcodes (same set as the barcode TXT files
# written above) so the browser tracks are consistent with the DA analysis.
#
# Requires bgzip and tabix on $PATH (part of htslib / samtools packages).
# If they are absent the uncompressed BED files are still written and the step
# exits cleanly.
# =============================================================================

print("Step 8: splitting fragments into per-cluster BED files for browser visualization...")
BED_DIR = os.path.join(OUTPUT_DIR, "cluster_fragments")
os.makedirs(BED_DIR, exist_ok=True)

try:
    print(f"  Reading fragment file: {FRAGMENT_FILE}")
    frag_df = pd.read_csv(
        FRAGMENT_FILE, sep="\t", header=None, comment="#",
        names=["chr", "start", "end", "barcode", "count"],
        dtype={"chr": str, "start": int, "end": int, "barcode": str, "count": int},
    )
    print(f"  Total fragments: {len(frag_df):,}")

    for cid in cluster_ids_smooth:
        cluster_barcodes = set(obs.index[obs["leiden_smooth"] == cid].tolist())
        subset = frag_df[frag_df["barcode"].isin(cluster_barcodes)]

        bed_path = os.path.join(BED_DIR, f"cluster_{cid}.bed")
        subset.to_csv(bed_path, sep="\t", header=False, index=False)
        print(f"  Cluster {cid}: {len(subset):,} fragments → {os.path.basename(bed_path)}")

        gz_path = bed_path + ".gz"
        try:
            subprocess.run(["bgzip", "-f", bed_path], check=True, capture_output=True)
            subprocess.run(["tabix", "-p", "bed", gz_path], check=True, capture_output=True)
            print(f"    bgzip + tabix: {os.path.basename(gz_path)}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            print(f"    bgzip/tabix not found — BED left uncompressed")

except Exception as e:
    print(f"  Fragment BED splitting failed: {e}")

print(f"\nPipeline complete. Outputs written to {OUTPUT_DIR}/")
print(f"  Total cells: {data.n_obs}")
print(f"  Chromatin clusters: {len(cluster_ids_raw)}  |  "
      f"Smoothed clusters: {len(cluster_ids_smooth)}")
print(f"  Spots relabeled by smoothing: {n_changed} ({100*n_changed/len(leiden_labels):.1f}%)")
print(f"  Per-cluster barcode lists (smoothed): {OUTPUT_DIR}/cluster_barcodes/")
print(f"  Per-cluster fragment BEDs (browser): {BED_DIR}/")